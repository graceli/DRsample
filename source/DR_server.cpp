/*
 *  This file is part of Distributed Replica.
 *  Copyright May 9 2009
 *
 *  Distributed Replica manages a series of simulations that separately sample phase space
 *  and coordinates their efforts under the Distributed Replica Potential Energy Function.
 *  See, for example T. Rodinger, P.L. Howell, and R. Pomès, "Distributed Replica Sampling"    
 *  J. Chem. Theory Comput., 2:725 (2006).
 *
 *  Distributed Replica is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Distributed Replica is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Distributed Replica.  If not, see <http://www.gnu.org/licenses/>.
 */

/*******************************************************************************************************************

 DR_server

 Creator: Tomas Rodinger (tom.rodinger@utoronto.ca)
 Program Suite Maintainer: Chris Neale (chris.neale@utoronto.ca, candrewn@gmail.com)
 Programmers: Tomas Rodinger, Ching-Hsing Yu, Chris Neale

 Description: Runs a distributed replica simulation in a spatial or temperature coordinate
	Performs either Monte Carlo moves, Boltzmann jumping in discrete space, or Boltzmann jumping in continuous space
 
 Usage: nohup DR_server scriptfile &
 scriptfile should have a format like title.script where title is two characters
 ex. t1.script

 These are the different type of simulations that can be run:

 coordinate_type   move_type          move data                           sizeof(move_data)
==============================================================================================
 Spatial           MonteCarlo         new_corr, system dE                 2 floats 
 Spatial           BoltjmannJumping   discrete Es                         Nreplicas floats
 Spatial           Continuous         Error: not supported                n/a
 Spatial           vRE                Error: not supported                n/a
 Spatial           NoMoves            None                                0
 Temperature       MonteCarlo         system E 	                          1 float
 Temperature       BoltjmannJumping   system E 	                          1 float
 Temperature       Continuous         system E 	                          1 float
 Temperature       vRE                system E                            1 float
 Temperature       NoMoves            None                                0
 Umbrella          MonteCarlo         position of what umbrella acts on   1 float
 Umbrella          BoltjmannJumping   position of what umbrella acts on   1 float
 Umbrella          Continuous         position of what umbrella acts on   1 float
 Umbrella          vRE                position of what umbrella acts on   1 float
 Umbrella          NoMoves            None                                0

*******************************************************************************************************************/

#define SNAPSHOT_VERSION 2.0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <sys/statvfs.h>
#include <signal.h>

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include "DR_protocol.h"
#include "force_database_class.h"
#include "read_input_script_file.h"
#include "vre.h"

#include <netinet/in.h>
#if defined(__ICC)
// icc doesn't like GNU __extension__ functions
// this has to happen AFTER !!
#undef htons
#undef ntohs
#undef htonl
#undef ntohl
#define htons(x) __bswap_constant_16(x)
#define ntohs(x) __bswap_constant_16(x)
#define htonl(x) __bswap_constant_32(x)
#define ntohl(x) __bswap_constant_32(x)
#endif // defined(_ICC)

//submit a job every 1 hour, no matter what
#define QUEUE_INTERVAL  3600 
#define DISK_ALMOST_FULL_CHECK_SECONDS 600 
#define FINISH_ON_AVERAGE_CHECK_SECONDS 600 
//MIN_DISK_SPACE_TO_RUN is in GB
#define MIN_DISK_SPACE_TO_RUN 1.0 
#define NODE_DISPLAY_SECONDS 600 
#define MOBILITY_CHECK_SECONDS 600
#define MAX_FAILURES_FOR_SUBMISSION 1000

//BUFFER_SIZE should be at least MAX_FILENAME_SIZE+1
#define BUFFER_SIZE 4096  
//REPLICA_MICRODIVISIONS is for continuous boltzmann jumping; this should be an odd nummber
#define REPLICA_MICRODIVISIONS 51 

#define USER_RESPONSIBILITY_STRING "I_TAKE_RESPONSIBILITY"
#define currentProgrammerName "Chris Neale"
#define currentProgrammerEmail "chris.neale@utoronto.ca"

const char *DEFINED_START_POS_FILE = "./switchStart.txt";

enum simulation_status_enum {Running, DiskAlmostFull, Finished, AllottedTimeOver};
enum energy_cancellation_status_enum {Disabled, Pending, Active, Active_and_Printed};

struct server_option_struct{
	bool loadSnapshot;
	char snapshotName[100];
	bool userTakesResponsibility;
	char userResponsibilityString[100];
	char title[3];
	int mobile_timeClientStarted;
	char logdir[500];
	int verbose;
};
#define DEFAULT_SERVER_OPTION_STRUCT {false,"",false,"","  ",0,"",0}

struct server_variable_struct{
	char working_directory[200];
	double beta;
	unsigned int Ncrashed_jobs;
	bool save_snapshot_now;
	enum energy_cancellation_status_enum energy_cancellation_status;
	int min_running_replica;
	int max_running_replica;
	unsigned int Nreserved_queue_slots;
	enum simulation_status_enum simulation_status;
	int Nconnected_clients;
	int Natoms;
	int nfailedsubinarow;
};
#define DEFAULT_SERVER_VARIABLE_STRUCT {"",1.0,0,false,Disabled,0,0,0,Running,0,0,0}

struct node_struct{
	bool active;
	char ip[50];
	unsigned int start_time;
	bool awaitingDump;
	char serverMessageForClient[500];
	bool messageWaitingIndicator;
};

struct client_bundle{
	//used to pass variables to wait_for_clients(), then client_interaction() via pthread_create()
	struct client_struct* client;
	struct server_option_struct *opt;
	struct server_variable_struct *var;
	struct script_struct *script;
	struct node_struct *node;
};

class force_database_class *force_database;

struct client_struct{
	int fd;
	struct timeval time;
	char *ptr;
	char log[10000];
	char ip[50];
};

#define MESSAGE_GLOBALVAR_LENGTH 10000               //reduce with caution. There is no overflow test
char logFile_globalVar[10];
// The replica_mutex is also used to control access to the node structure
pthread_mutex_t replica_mutex;
pthread_mutex_t log_mutex;
pthread_mutex_t queue_mutex;
pthread_mutex_t database_mutex;
pthread_mutex_t vre_mutex;

char months[12][4]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

float circularCancelLow, circularCancelHigh;
float newCircularCancelLow, newCircularCancelHigh;


//======================================================================================================================================

void debugflush(void){
	fflush(stdout);
	fflush(stderr);
}

double sqr(double x){
	return(x*x);
}

// Appends the given text to the log file
// code here is used to specify the file descriptor of the connected client (specifying -1 prints no code)
void append_log_entry(int code, const char *entry){
	int fd;
	time_t t;
	struct tm *tt;
	char time_string[40];
	int writecheck;

	pthread_mutex_lock(&log_mutex);

	t=time(NULL);
	tt=localtime(&t);
	if(code<0) sprintf(time_string,"[%s/%.2d/%.4d %.2d:%.2d:%.2d] ",months[tt->tm_mon],tt->tm_mday,tt->tm_year+1900,tt->tm_hour,tt->tm_min,tt->tm_sec);
	else       sprintf(time_string,"[%s/%.2d/%.4d %.2d:%.2d:%.2d] <%d> ",months[tt->tm_mon],tt->tm_mday,tt->tm_year+1900,tt->tm_hour,tt->tm_min,tt->tm_sec,code);

	if( (fd=open(logFile_globalVar,O_WRONLY|O_CREAT|O_APPEND, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))==-1 ){
		perror("log file");
		fprintf(stderr,"Error: cannot open log file for appending: %s\n",logFile_globalVar);
		fflush(stderr);
		exit(1);
	}
	fchmod(fd,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);

	// we catch the variable to avoid compiler warnings, but its just a log file so who cares...
	writecheck=write(fd,time_string,strlen(time_string));
	writecheck=write(fd,entry,strlen(entry));

	close(fd);
	
	pthread_mutex_unlock(&log_mutex);
}

// prints an error message to the log file and exits. If a log file has not been defined yet, the message goes to stdout
void error_warning(const char *error_message){
	char message[MESSAGE_GLOBALVAR_LENGTH];

	sprintf(message,"Error ERROR error: %s\n",error_message);
	if(logFile_globalVar[0]!=0)
		append_log_entry(-1,message);
	perror(message);
	fflush(stderr);
}

void error_quit(const char *error_message){
	error_warning(error_message);
	exit(1);
}

//YUPYUPYUP

int findRepByNode(const struct script_struct *script, int nodeSlot){
	//not sure how well this will work for NNI!=1
	int i;
	
	for(i=0;i<script->Nreplicas; i+=script->Nsamesystem_uncoupled){
		//script->replica[i].nodeSlot>=0 required for consistency with other code as a first check... not really important here though
		if(script->replica[i].nodeSlot>=0 && script->replica[i].nodeSlot==nodeSlot){
			//no check to see if multiple replicas think they are on the same node
			return i;
		}
	}
	//not an error, but replica is expected to ne not running...
	return -1;
}

void displayNodes(const struct script_struct *script, const struct node_struct *node){
	char message[MESSAGE_GLOBALVAR_LENGTH];
	int i,r;

	for(i=0;i<div(script->Nreplicas,script->Nsamesystem_uncoupled).quot; i++){
		if(node[i].active){
			r=findRepByNode(script,i);
			if(r>=0){
				sprintf(message,"NODE[%d] %s started %d has replica %d\n",i,node[i].ip,node[i].start_time,r);
			}else{
				sprintf(message,"NODE[%d] %s started %d has replica NONE\n",i,node[i].ip,node[i].start_time);
			}
		}else{
			sprintf(message,"NODE[%d] index is not currently being used\n",i);
		}
		append_log_entry(-1,message);
	}
}


int findInactiveNodeSlot(const struct script_struct *script, const struct node_struct *node){
	int i;

	for(i=0;i<div(script->Nreplicas,script->Nsamesystem_uncoupled).quot; i++){
		if(!node[i].active){
			return i;
		}
	}
	//massive error if I get to here
	return -1;
}

int findNodeByIP(const struct script_struct *script, const struct node_struct *node, const char *ip){
	int i;
	
	for(i=0;i<div(script->Nreplicas,script->Nsamesystem_uncoupled).quot; i++){
		if(node[i].active && strcmp(node[i].ip,ip)==0){
			return i;
		}
	}
	//it's ok if I get to here
	return -1;
}

void connectNodeToReplica(int *repSaveNodeSlot, int myNewNodeSlot){
	*(repSaveNodeSlot)=myNewNodeSlot;
}

void obtainNode(struct node_struct *node, int *repSaveNodeSlot, const char *ip, int myNewNodeSlot, int start_time){
	node->active=true;
	strcpy(node->ip,ip);
	if(start_time>0){
		//this came through from the client
		node->start_time=start_time;
	}else{
		node->start_time=time(NULL);
	}
	node->awaitingDump=false;
	connectNodeToReplica(repSaveNodeSlot,myNewNodeSlot);
}

void disconnectNodeFromReplica(int *repSaveNodeSlot){
	*(repSaveNodeSlot)=-1;
}

void releaseNode(struct node_struct *node, int *repSaveNodeSlot){
	node->active=false;
	//node->ip[0]='\0'; // this is not required and it complicates reporting
	disconnectNodeFromReplica(repSaveNodeSlot);
}

/* Usage:
 * 
 * obtainNode(&(B->node[myNewNodeSlot]),&(B->script->replica[replicaN[0]].nodeSlot),B->client->ip,myNewNodeSlot,(int)*((float*)(tcs->data)));
 * connectNodeToReplica(&(B->script->replica[replicaN[0]].nodeSlot),myNewNodeSlot);
 * 
 * releaseNode(&(node[script->replica[replicaN].nodeSlot]),&(script->replica[replicaN].nodeSlot));
 * disconnectNodeFromReplica(&(script->replica[replicaN].nodeSlot));
 */


int saveVREtoFile(const char *fnam,int Nreplicas){
	FILE *f;
	f=fopen(fnam,"w");
	if(f==NULL){
		perror("Unable to write the VRE file -- unable to open file\n");
		return 1;
	}
	showVRE(Nreplicas,f);
	fclose(f);
	return 0;
}

// Saves a snapshot of the state of the distributed replica simulation to a file
// The snapshot contains the coordinate positions of all replicas, their current sequence number,
// a restart file, and the averaged coordinates at each discrete replica position.
// Also save vRE structure.
char * save_snapshot(const struct script_struct *script, const struct server_variable_struct *var, const struct server_option_struct *opt){
	//CN moved the replica mutex lock to something that must be done outside of this routine -- necessary for mobile server
	static char filename[30];
	int fd;
	unsigned int size;
	float version=SNAPSHOT_VERSION;
	
	append_log_entry(-1,"Saving a state snapshot\n");

	sprintf(filename,"%s.%d.snapshot",opt->title,(int)time(NULL));
	if( (fd=open(filename, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))==-1 ) error_quit("cannot open file for writing");

	size=script->Nreplicas*sizeof(struct replica_struct);

	if( write(fd,&version,sizeof(version))!=sizeof(version) ) error_quit("cannot write to file");
	if( write(fd,&(script->Nreplicas),sizeof(script->Nreplicas))!=sizeof(script->Nreplicas) ) error_quit("cannot write to file");
	if( write(fd,&(var->Natoms),sizeof(var->Natoms))!=sizeof(var->Natoms) ) error_quit("cannot write to file");

	if( write(fd,script->replica,size)!=size ) error_quit("cannot write to file");
	for(int i=0;i<script->Nreplicas;i++){
		size=script->replica[i].restart.data_size;
		if( write(fd,script->replica[i].restart.data,size)!=size ) error_quit("cannot write to file");

		size=var->Natoms*sizeof(struct atom_struct);
		if( write(fd,script->replica[i].atom,size)!=size ) error_quit("cannot write to file");

		size=N_PRESENCE_BITS/8;
		if( write(fd,script->replica[i].presence,size)!=size ) error_quit("cannot write to file");
	}

	if(script->replica_move_type==vRE){
		long int nallocated,nlastused,nrecyclepush;
		float *val;
		struct vre_item_struct *item;

		for(int i=0;i<script->Nreplicas;i++){
			save_vre_pointers(i,&nallocated,&nlastused,&item);
			if( write(fd,&nallocated,sizeof(nallocated))!=sizeof(nallocated) ) error_quit("cannot write to file -- vre nallocated");
			if( write(fd,&nlastused,sizeof(nlastused))!=sizeof(nlastused) ) error_quit("cannot write to file -- vre nlastused");
			if( write(fd,item,sizeof(struct vre_item_struct)*(nlastused+1))!=sizeof(struct vre_item_struct)*(nlastused+1) ) error_quit("cannot write to file -- vre item");
		}
		for(int i=0;i<script->Nreplicas;i++){
			save_secvre_pointers(i,&nallocated,&nlastused,&nrecyclepush,&val);
			if( write(fd,&nallocated,sizeof(nallocated))!=sizeof(nallocated) ) error_quit("cannot write to file -- vre secondary nallocated");
			if( write(fd,&nlastused,sizeof(nlastused))!=sizeof(nlastused) ) error_quit("cannot write to file -- vre secondary nlastused");
			if( write(fd,&nrecyclepush,sizeof(nrecyclepush))!=sizeof(nrecyclepush) ) error_quit("cannot write to file -- vre secondary nrecyclepush");
			if( write(fd,val,sizeof(float)*(nlastused+1))!=sizeof(float)*(nlastused+1) ) error_quit("cannot write to file -- vre secondary val");
		}
		//FOR DEBUGGING PURPOSES ONLY
		//char saveName[30];
		//sprintf(saveName,"VRE_saved.txt");
		//saveVREtoFile(saveName,script->Nreplicas);
	}

	close(fd);
	return (filename);
}

//Loads a snapshot of the state of the distributed replica simulation from a file
// The snapshot contains the coordinate positions of all replicas, their current sequence number,
// a restart file, and the averaged coordinates at each discrete replica position.
// Also load vRE structure.
void load_snapshot(char *filename, struct script_struct *script, struct server_variable_struct *var){
	int fd;
	unsigned int size;
	int Nreplicas_in_snapshot;
	int i;
	float version;
	unsigned int sampling_runs;
	char message[MESSAGE_GLOBALVAR_LENGTH];
	
	sprintf(message,"Loading state snapshot: %s\n", filename);
	append_log_entry(-1,message);

	if( (fd=open(filename,O_RDONLY))==-1 ) error_quit("cannot open file for reading");
	if( read(fd,&version,sizeof(version))!=sizeof(version) ) error_quit("cannot read from file");
	if(version!=SNAPSHOT_VERSION){
		if(!(script->replica_move_type!=vRE && version==1.0 && SNAPSHOT_VERSION==2.0)){
			error_quit("this program cannot read this version of the snapshot"); 
		}
		append_log_entry(-1,"ERROR error Error: The current SNAPSHOT_VERSION is 2.0, and your snapshot is version 1.0. However, you are not using vRE so this is allowed. Note: use at your own risk!!! (talk to Chris Neale if you want some assistance here).\n");
	}
	if( read(fd,&Nreplicas_in_snapshot,sizeof(Nreplicas_in_snapshot))!=sizeof(Nreplicas_in_snapshot) ) error_quit("cannot read from file");
	if(Nreplicas_in_snapshot!=script->Nreplicas) error_quit("number of replicas in the snapshot and in script file don't match"); 
	if( read(fd,&(var->Natoms),sizeof(var->Natoms))!=sizeof(var->Natoms) ) error_quit("cannot read from file");
	size=sizeof(struct replica_struct);
	if(script->replica==NULL){ 
		//CN wonders why this does not invoke an error_quit
		script->replica=new replica_struct[script->Nreplicas];
	}
	for(i=0;i<script->Nreplicas;i++){
		sampling_runs=script->replica[i].sampling_runs;  
		// we want to keep the quantity from the script file, not the one in the snapshot
		if( read(fd,&(script->replica[i]),size)!=size ) error_quit("cannot read from file");
		script->replica[i].sampling_runs=sampling_runs;
	}

	printf("Read in initial snapshot data, Nreplicas: %u   Natoms: %u\n",script->Nreplicas,var->Natoms); //##DEBUG

	for(i=0;i<script->Nreplicas;i++){
		printf("Reading replica %d, status: %c , sample_count: %u ,  sampling_runs: %u\n",i,script->replica[i].status, script->replica[i].sample_count, script->replica[i].sampling_runs); //##DEBUG

		if( (script->replica[i].status=='S') || (script->replica[i].status=='R') ) script->replica[i].status='N';
		
		size=script->replica[i].restart.data_size;
		printf("allocating memory with size %u\n",size); //##DEBUG
		script->replica[i].restart.data=new unsigned char[size];
		script->replica[i].restart.allocated_memory=script->replica[i].restart.data_size;
		if( read(fd,script->replica[i].restart.data,size)!=size ) error_quit("cannot read from file");

		size=var->Natoms*sizeof(struct atom_struct);
		printf("allocating memory with size %u\n",size); //##DEBUG
		script->replica[i].atom=new atom_struct[var->Natoms];
		if( read(fd,script->replica[i].atom,size)!=size ) error_quit("cannot read from file");

		size=N_PRESENCE_BITS/8;
		printf("allocating memory with size %u\n",size); //##DEBUG
		script->replica[i].presence=new unsigned int[N_PRESENCE_BITS/32];
		if( read(fd,script->replica[i].presence,size)!=size ) error_quit("cannot read from file");
	}

	if(script->replica_move_type==vRE){
		long int *nallocated,*nlastused, *nrecyclepush;
		float **val;
		struct vre_item_struct **item;
		long int vreAllocated;

		vreAllocated=getVREallocationLevel(script->Nreplicas);
		for(int i=0;i<script->Nreplicas;i++){
			load_vre_pointers(i,&nallocated,&nlastused,&item);
			if(read(fd,nallocated,sizeof(*nallocated))!=sizeof(*nallocated) ) error_quit("cannot read from file");
			if(read(fd,nlastused,sizeof(*nlastused))!=sizeof(*nlastused) ) error_quit("cannot read from file");
			if(*nlastused>=vreAllocated) error_quit("Loading in a vRE structure for which memory is not available.");
			if(read(fd,*item,sizeof(struct vre_item_struct)*(*nlastused+1))!=sizeof(struct vre_item_struct)*(*nlastused+1) ) error_quit("cannot read from file");
		}
		vreAllocated=getSECVREallocationLevel(script->Nreplicas);
		for(int i=0;i<script->Nreplicas;i++){
			load_secvre_pointers(i,&nallocated,&nlastused,&nrecyclepush,&val);
			if(read(fd,nallocated,sizeof(*nallocated))!=sizeof(*nallocated) ) error_quit("cannot read from file");
			if(read(fd,nlastused,sizeof(*nlastused))!=sizeof(*nlastused) ) error_quit("cannot read from file");
			if(*nlastused>=vreAllocated) error_quit("Loading in a vRE structure for which memory is not available.");
			if(read(fd,nrecyclepush,sizeof(*nrecyclepush))!=sizeof(*nrecyclepush) ) error_quit("cannot read from file");
			if(read(fd,*val,sizeof(float)*(*nlastused+1))!=sizeof(float)*(*nlastused+1) ) error_quit("cannot read from file");
		}
		//FOR DEBUGGING PURPOSES ONLY
		//char saveName[30];
		//sprintf(saveName,"VRE_loaded.txt");
		//saveVREtoFile(saveName,script->Nreplicas);
	}
	close(fd);
}

//frees all memory used to store replica information
void free_all_replicas(struct script_struct *script){
	for(int i=0;i<script->Nreplicas;i++){
		delete[] script->replica[i].restart.data;
		delete[] script->replica[i].atom;
		delete[] script->replica[i].presence;
	}
	delete[] script->replica;
}

// For now basically checks that the received restart file has a size greater than zero
unsigned char check_restart_file_integrity(struct client_struct *client, struct buffer_struct restart){
	if(restart.data_size>0){
		return(1);
	}else{
		client->ptr+=sprintf(client->ptr,"Restart data fails integrity check; size is %u\n",restart.data_size);
		return(0);
	}
}
	
// Accepts the given restart file as the latest restart file for the specified replica
void commit_restart_file(int replicaN, struct buffer_struct *restart, struct script_struct *script){
	printf("commit_restart_file: freeing memory at: %p\n",script->replica[replicaN].restart.data); //##DEBUG
	delete[] script->replica[replicaN].restart.data;
	script->replica[replicaN].restart=*restart;
	printf("comitted restart file now at: %p\n",script->replica[replicaN].restart.data); //##DEBUG
	restart->data=NULL;
	restart->data_size=0;
	restart->allocated_memory=0;
}

// Averages the given coordinate file into the master coordinate file for the given replica
void commit_coordinate_data(int replicaN, int bin_number, struct buffer_struct *coordinate_v, const struct script_struct *script, struct server_variable_struct *var){
	unsigned int address;
	unsigned char bit;
	unsigned int presence;
	float *coordinate=(float*)(coordinate_v->data);
	int i;

	if(var->Natoms==0){
		printf("Natoms is 0, therefore we need to allocate memory for the first time\n");   //##DEBUG
		var->Natoms=coordinate_v->data_size/sizeof(float)/3; // 3 coordinates per atom
		printf("Natoms is now %d\n",var->Natoms);                                                //##DEBUG
		
		for(i=0;i<script->Nreplicas;i++){
			script->replica[i].atom=new atom_struct[var->Natoms];
			printf("Replica %d: memory successfully allocated\n",i);                    //##DEBUG
			memset(script->replica[i].atom,0,var->Natoms*sizeof(struct atom_struct));
			printf("%d bytes filled with zeros\n",(int)(var->Natoms*sizeof(struct atom_struct)));   //##DEBUG
		}	
	}	

	printf("--------> committing coordinate data for replica %d, bin: %d\n",replicaN, bin_number); //##DEBUG
	
	bit=script->replica[replicaN].sequence_number&31;
	address=script->replica[replicaN].sequence_number>>5;
	
	printf("calculating presence location: sequence_number: %u   address: %u   bit: %hhu\n",script->replica[replicaN].sequence_number, address, bit); //##DEBUG
	presence=script->replica[replicaN].presence[address];
	
	printf("current presence: %u\n",presence); //##DEBUG
	
	if( ((presence>>bit)&1)==0 ){
		printf("the presence bit was 0\n"); //##DEBUG

		if(coordinate_v->data_size/sizeof(float)/3==var->Natoms) // this condition will be false only on one very rare circumstance: when one or more corrupt coordinate files are received before the Natoms variable has been updated with the correct number
			for(i=0;i<var->Natoms;i++){
				script->replica[bin_number].atom[i].x+=coordinate[i*3+0];
				script->replica[bin_number].atom[i].y+=coordinate[i*3+1];
				script->replica[bin_number].atom[i].z+=coordinate[i*3+2];
				script->replica[bin_number].atom[i].weight++;
			}

		printf("new atoms 0 data: %lf %lf %lf  %u\n",script->replica[bin_number].atom[0].x,script->replica[bin_number].atom[0].y,script->replica[bin_number].atom[0].z,script->replica[bin_number].atom[0].weight); //##DEBUG

		presence|=(1<<bit);
		script->replica[replicaN].presence[address]=presence;
		printf("new presence: %u\n",presence); //##DEBUG
	}
}

// Makes a duplicate copy of the given restart data
void copy_restart_data(struct buffer_struct *destination, struct buffer_struct *source){
	printf("copy_restart_data: destination pointer is: %p\n",destination->data);  //##DEBUG
	unsigned int size=source->data_size;
	destination->data_size=size;
	destination->allocated_memory=size;
	destination->data=new unsigned char[size];
	memcpy(destination->data,source->data,size);
}

// choses the best replica to run next such as to minimize network trafic
// this algorithm chooses the replica with the lowest sequence number, however, if more than one replica
// have the same lowest sequence number then preference will be given to the replica that just ran on that Node,
// this is so that a restart file need not be sent to the client over the nework.
void find_replica_to_run(int *replicaN, const struct script_struct *script){
	unsigned int min_sequence_number=2000000000;
	int i;

	if(*replicaN>=0){
		if(script->replica[*replicaN].status=='R'){ 
			return;
		}else if (script->replica[*replicaN].status=='N'){
			min_sequence_number=script->replica[*replicaN].sequence_number;
		}
	}
	for(i=0;i<script->Nreplicas;i++){
		if(script->replica[i].status=='N' && script->replica[i].sequence_number<min_sequence_number){
			min_sequence_number=script->replica[i].sequence_number;
			*replicaN=i;
		}
	}
	if(script->replica[*replicaN].status!='N') *replicaN=-1;
}

// Finds the oldest client and drops it if it is older than 50% of its runtime.
// This is useful if a new job has just started, but there is nothing left to run -- don't waste the new walltime
// HOWEVER: can only allow one new job to be waiting on each old job.
int drop_one_old_node(struct client_struct *client, const struct script_struct *script, struct node_struct *node){
	int i, onode, otime, age;

	//The replica_mutex must be on already!!!

	// find the oldest Node
	otime=-1;
        onode=-1;
	for(i=0;i<script->Nreplicas;i+=script->Nsamesystem_uncoupled){
		if(script->replica[i].status=='R' && script->replica[i].nodeSlot>=0){
			// if replica.status=='R' && replica.nodeSlot<0, it represents a huge error in any event
			if(!node[script->replica[i].nodeSlot].awaitingDump && (otime==-1 || node[script->replica[i].nodeSlot].start_time<otime)){
				onode=script->replica[i].nodeSlot;
				otime=node[onode].start_time;
			}
		}
	}
	if(onode<0){
		//can happen if all nodes are currently awaiting a dump
		return (onode);
	}
	age=time(NULL)-otime;
	// decide if it is old enough (currently, only allow stoppage if at least script->cycleClients of the allocated time has been used)
	client->ptr+=sprintf(client->ptr,"DUMP: Considering dump of age %d based on Nodetime %d and eval %d \n", age,script->node_time, (int)ceil((double)(script->node_time)*script->cycleClients));
	if(age>(int)ceil((double)(script->node_time)*script->cycleClients)){
		// it's old enough, so change its clock to cause it to termiante early.
		node[onode].awaitingDump=true;
		node[onode].start_time-=script->node_time;
	}else{
		// not old enough... return -1
		onode=-1;
	}
	if(onode!=-1){
		// wait for the dumped replica to finish its communication by other methods
		pthread_mutex_unlock(&replica_mutex);
		while(node[onode].active && node[onode].awaitingDump){
			//node[onode].awaitingDump is necessary because a different node could pick this one up (that would deadlock this one)
			sleep(1);
			fprintf(stderr,"drop_one_old_node (WAITING): Waiting for node[%d] (ip=%s) to become inactive.\n",onode, node[onode].ip);fflush(stderr);
		}
		pthread_mutex_lock(&replica_mutex);
		//no need to node[onode].awaitingDump=false because that is done in connection to a new node
	}
	return (onode);
}

// checks that the received energy data has the correct size
// CN thinks that a const pointer should be sent for energy, but will leave it alone for now
unsigned char check_energy_file_integrity(struct client_struct *client, struct buffer_struct energy, const struct script_struct *script){
	unsigned int expected_file_size;
	
	// see the table at the top of this file to see how these conditions were derived
	if(script->replica_move_type==NoMoves) expected_file_size=0;
	else if(script->coordinate_type==Temperature || script->coordinate_type==Umbrella ) expected_file_size=sizeof(float);
	else if(script->replica_move_type==MonteCarlo) expected_file_size=2*sizeof(float);
	else expected_file_size=script->Nreplicas*sizeof(float);
	printf("checking energy file integrity: correct size is %u ; actual size is: %u\n",expected_file_size,energy.data_size); //##DEBUG
	if(energy.data_size==expected_file_size){
		return(1);
	}else{
		client->ptr+=sprintf(client->ptr,"Energy data fails integrity check; expected size is %u; acutal size is %u\n",expected_file_size,energy.data_size);
		return(0);
	}
}

// CN thinks that a const pointer should be sent for sample, but will leave it alone for now
unsigned char check_sample_file_integrity(struct client_struct *client, struct buffer_struct sample, unsigned int sampleoradditional, const struct script_struct *script){
/* check_sample_file_integrity checks that the sample (forces) and additional (arbitrary) data has the correct size
 * set sampleoradditional=0 for sample data (size based on Nsamples_per_run*Nligands)
 * set sampleoradditional>0 for additional data (size based on Nsamples_per_run)
 */
	unsigned int expected_file_size;
	unsigned int multiplier;

	if(sampleoradditional==0){
		multiplier=script->Nligands;
	}else{
		multiplier=1;
	}
	expected_file_size=script->Nsamples_per_run*multiplier*sizeof(float);

	printf("checking sample file integrity: correct size is %u ; actual size is: %u\n",expected_file_size,sample.data_size); //##DEBUG
	
	if(sample.data_size==expected_file_size){
		return(1);
	}else{
		client->ptr+=sprintf(client->ptr,"Sample data fails integrity check; expected size is %u; acutal size is %u\n",expected_file_size,sample.data_size);
		return(0);
	}
}

// checks that the received coordinate data has information for the right number of atoms
// CN thinks that a const pointer should be sent for coordinate, but will leave it alone for now
unsigned char check_coordinate_file_integrity(struct client_struct *client, struct buffer_struct coordinate, const struct server_variable_struct *var){
	unsigned int expected_file_size;
	
	expected_file_size=var->Natoms*sizeof(float)*3;

	printf("checking coordinate file integrity: correct size is %u ; actual size is: %u\n",expected_file_size,coordinate.data_size); //##DEBUG
	
	for(int i=0;i<2;i++){ //##DEBUG
		printf("%f %f %f\n",((float*)coordinate.data)[i*3+0],((float*)coordinate.data)[i*3+1], //##DEBUG
			((float*)coordinate.data)[i*3+2]); //##DEBUG
	} //##DEBUG
	printf("last 1:\n"); //##DEBUG
	for(int i=var->Natoms-2;i<var->Natoms;i++){ //##DEBUG
		printf("%f %f %f\n",((float*)coordinate.data)[i*3+0],((float*)coordinate.data)[i*3+1], //##DEBUG
			((float*)coordinate.data)[i*3+2]); //##DEBUG
	} //##DEBUG
	
	if( (coordinate.data_size==0) || ( (coordinate.data_size!=expected_file_size) && (var->Natoms!=0) ) ){
		client->ptr+=sprintf(client->ptr,"Coordinate data fails integrity check; expected size is %u; acutal size is %u\n",expected_file_size,coordinate.data_size);
		return(0);
	}
	
	return(1);
}

// transforms the given coordinate value ('w') into a coordinate space with uniform spacing between the nominal positions of replicas
float replica_linearizing_function(float w, const struct script_struct *script){
	int i;
	float fraction;
	
	for(i=0;i<script->Nreplicas;i++){
		if(w<script->replica[i].w_nominal) break;
	}
		
	if(i==0) i++;
	if(i==script->Nreplicas) i--;
	fraction=(w-script->replica[i-1].w_nominal)/(script->replica[i].w_nominal-script->replica[i-1].w_nominal);
	
	return((float)(i-1)+fraction);
}

// calculates the potential of replica districution. ie. this comupted the DRPE (see the paper)
double replica_potential(const struct script_struct *script){
	int i,j;
	double E_total;
	float swap_temp;

	for(i=0;i<script->Nreplicas;i++){
		script->replica[i].w_sorted=script->replica[i].w;
	}

	for(i=0;i<script->Nreplicas-1;i++){
		for(j=i+1;j<script->Nreplicas;j++){
			if(script->replica[i].w_sorted>script->replica[j].w_sorted){
				swap_temp=script->replica[i].w_sorted;
				script->replica[i].w_sorted=script->replica[j].w_sorted;
				script->replica[j].w_sorted=swap_temp;
			}
		}
	}

	for(i=0;i<script->Nreplicas;i++){
		script->replica[i].w_sorted=replica_linearizing_function(script->replica[i].w_sorted,script);
	}

	E_total=0.0;
	for(i=0;i<script->Nreplicas;i++){
		for(j=0;j<script->Nreplicas;j++){
			E_total+=pow((double)(script->replica[i].w_sorted-script->replica[j].w_sorted)-(double)(i-j),2.0);
		}
	}

	E_total*=script->replica_potential_scalar1;

	double w_sum=0.0;
	double w_nominal_sum=((double)script->Nreplicas-1.0)/2*script->Nreplicas; 
	//this is really the sum over i from 0 to Nreplicas-1 of replica_linearizing_function(replica[i].w_nominal)
	//it's just faster to calculate this way
	for(i=0;i<script->Nreplicas;i++){
		w_sum+=script->replica[i].w_sorted;
	}

	E_total+=pow(w_sum-w_nominal_sum,2.0)*script->replica_potential_scalar2;
	
	return(E_total);
}

// determine where we will attempt to move to next given that we currently are at 'w'
float calculate_monte_carlo_move(float w, const struct script_struct *script){
	int i;
	int integer;
	double fraction;
	double circReturn;

	for(i=0;i<script->Nreplicas;i++){
		if(w<script->replica[i].w_nominal) break;
	}

	if(i==0) i++;
	else if(i==script->Nreplicas) i--;

	fraction=(w-script->replica[i-1].w_nominal)/(script->replica[i].w_nominal-script->replica[i-1].w_nominal);

	printf("Calculating change in coordinate, w: %lf   i: %d   fraction: %lf\n",w,i,fraction); //##DEBUG

	if(drand48()>=0.5){
		fraction+=script->replica_step_fraction;
		integer=(int)fraction;
		if(i+integer>script->Nreplicas-1) integer=script->Nreplicas-1-i;
		i+=integer;
		fraction-=(double)integer;
	}else{
		fraction-=script->replica_step_fraction;
		integer=-(int)(fraction-1);
		if(i-integer<1) integer=i-1;
		i-=integer;
		fraction+=(double)integer;
	}

	if(script->circular_replica_coordinate){
		circReturn=fraction*(script->replica[i].w_nominal-script->replica[i-1].w_nominal)+script->replica[i-1].w_nominal;
		if(circReturn>script->circular_greater_equality){
			return(script->replica_move_type==MonteCarlo)?script->replica[0].w_nominal:circReturn-script->circular_equality_distance;
		}
		if(circReturn<script->circular_lesser_equality){
			return(script->replica_move_type==MonteCarlo)?script->replica[script->Nreplicas-1].w_nominal:circReturn+script->circular_equality_distance;
		}
		return(circReturn);
	}

	return( fraction*(script->replica[i].w_nominal-script->replica[i-1].w_nominal)+script->replica[i-1].w_nominal );

	//printf("new_w = %lf   new_w2 = %lf\n",*new_w,*new_w2);  //##DEBUG
}

// determine which discrete nominal coordinate position 'w' is closes to
int find_bin_from_w(double w, const struct script_struct *script){
	double left_space, right_space;
	double left_limit, right_limit;
	int i;

	if(script->Nreplicas==1) return(0);

	for(i=0;i<script->Nreplicas;i++){
		if(i>0){
			left_space=script->replica[i].w_nominal-script->replica[i-1].w_nominal;
		}
		if(i<script->Nreplicas-1){
			right_space=script->replica[i+1].w_nominal-script->replica[i].w_nominal;
		}
		if(i==0){
			left_space=right_space;
		}else if(i==script->Nreplicas-1){
			right_space=left_space;
		}
		
		left_limit=script->replica[i].w_nominal-left_space*0.5;
		right_limit=script->replica[i].w_nominal+right_space*0.5;
		
		if( (w>=left_limit) && (w<=right_limit) ) return(i);
	}
	return(-1);
}

// determine where the given replica will move to next using the Monte Carlo move scheme
// energy data contains the required data to do the move attempt, see comments below
float determine_new_replica_position_monte_carlo_or_vre(struct client_struct *client, int replicaN, const float *energy_data, const struct script_struct *script, const struct server_variable_struct *var){
	double new_coor, old_coor;
	int new_bin, old_bin;
	double new_cancellation, old_cancellation;
	double old_DRPE, new_DRPE, system_energy_change, DRPE_change, total_energy_change, probability;
	double random_number;
	double newdist,olddist;
	float vrepop;
	int vrecheck,vresource;

	old_coor=script->replica[replicaN].w;

	if(script->coordinate_type==Spatial){
		// for spatial simulations, the first item in the array is the new w coordinate
		new_coor=energy_data[0];
	}else{
		new_coor=calculate_monte_carlo_move(old_coor,script);
	}
	
	new_bin=find_bin_from_w(new_coor,script);
	old_bin=find_bin_from_w(old_coor,script);

	if(script->replica_move_type==vRE){
		if(script->replica[replicaN].sequence_number>=script->vRE_initial_noSave){
			pthread_mutex_lock(&vre_mutex);
			pushVRE(old_bin,replicaN,energy_data[0]);
			pthread_mutex_unlock(&vre_mutex);
		}else{
			client->ptr+=sprintf(client->ptr,"Not storing value for later use as a virtual move for the first %ld steps in vRE\n",script->vRE_initial_noSave);
		}
		if(script->replica[replicaN].sequence_number<script->vRE_initial_noMoves){
			client->ptr+=sprintf(client->ptr,"Not allowing movement for the first %ld steps in vRE\n",script->vRE_initial_noMoves);
			return(old_coor);
		}
	}

	if(new_bin<var->min_running_replica || new_bin>var->max_running_replica){
		client->ptr+=sprintf(client->ptr,"Move rejected because it crosses a suspend boundary, keeping coordinate: %lf\n",old_coor);
		return(old_coor);
	}
		
	old_DRPE=replica_potential(script);
	script->replica[replicaN].w=new_coor;
	new_DRPE=replica_potential(script);
	DRPE_change=new_DRPE-old_DRPE;
		
	new_cancellation=script->replica[new_bin].cancellation_energy;
	old_cancellation=script->replica[old_bin].cancellation_energy;

	if(script->replica_move_type==vRE){
		pthread_mutex_lock(&vre_mutex);
		vrecheck=popVRE(new_bin,replicaN,&vrepop,&vresource);
		//fprintf(stderr,"RECEIVED POP: %f\n",vrepop);
		pthread_mutex_unlock(&vre_mutex);
		if(vrecheck!=0){
			client->ptr+=sprintf(client->ptr,"VRE_INFO:(noMove) the vre structure has no values. Move not allowed. Keeping coordinate: %lf\n",old_coor);
			script->replica[replicaN].w=old_coor;
			return(old_coor);
		}
		client->ptr+=sprintf(client->ptr,"VRE_INFO:(replica,nominal,value)_real_then_virtual %d %lf %f %d %lf %f\n",replicaN,old_coor,energy_data[0],vresource,new_coor,vrepop);
	}


	if(script->coordinate_type==Spatial){
		//Note that circular coordinate energy correction is the responsibility 
		//of the client program for Spatial coordinates
		system_energy_change=energy_data[1];   
		//for spatial simulations, the second item in the array is 
		//the system energy change to go to the new w
		system_energy_change+=new_cancellation-old_cancellation;
		system_energy_change*=var->beta;
		DRPE_change*=var->beta;
		total_energy_change=DRPE_change+system_energy_change;      // energies here are dimensionless
		client->ptr+=sprintf(client->ptr,"[system change, cancellation_change, DRPE change, total change]: %lf %lf %lf %lf\n",system_energy_change,new_cancellation-old_cancellation,DRPE_change,total_energy_change);
	}else if(script->coordinate_type==Temperature){
		system_energy_change=(new_coor-old_coor)*energy_data[0];
		if(script->replica_move_type==vRE){
			system_energy_change+=(old_coor-new_coor)*(double)vrepop;
		}
		// for temperature simulations, the first item in the array is the system energy; 
		// the system_energy_change calculated here is dimensionless
		system_energy_change+=old_coor*old_cancellation-new_coor*new_cancellation;
		total_energy_change=DRPE_change+system_energy_change; // calculate total dimensionless energy change for the move
	}else{ // "Umbrella" simulation
		double exact_coordinate_position=energy_data[0];  
		//for umbrealla simulations, the first item in the array is the exact 
		//coordinate value of the entity that the umbrella acts on
		newdist=exact_coordinate_position-new_coor;
		olddist=exact_coordinate_position-old_coor;
		if(script->circular_replica_coordinate){
			//Ensure that we use the closest periodic distance. 
			if(newdist>=script->circular_equality_distance/2.0)newdist-=script->circular_equality_distance;
			if(newdist<=-script->circular_equality_distance/2.0)newdist+=script->circular_equality_distance;
			if(olddist>=script->circular_equality_distance/2.0)olddist-=script->circular_equality_distance;
			if(olddist<=-script->circular_equality_distance/2.0)olddist+=script->circular_equality_distance;
		}
		// CN says INCORRECT: system_energy_change=0.5*script->replica[replicaN].force*(sqr(newdist)-sqr(olddist));
		system_energy_change=0.5*(script->replica[new_bin].force*sqr(newdist) - script->replica[old_bin].force*sqr(olddist));

client->ptr+=sprintf(client->ptr,"MEOW: real pos %f energy %f ",exact_coordinate_position,system_energy_change);

		if(script->replica_move_type==vRE){
			newdist=(double)vrepop-old_coor;
			olddist=(double)vrepop-new_coor;
			if(script->circular_replica_coordinate){
				if(newdist>=script->circular_equality_distance/2.0)newdist-=script->circular_equality_distance;
				if(newdist<=-script->circular_equality_distance/2.0)newdist+=script->circular_equality_distance;
				if(olddist>=script->circular_equality_distance/2.0)olddist-=script->circular_equality_distance;
				if(olddist<=-script->circular_equality_distance/2.0)olddist+=script->circular_equality_distance;
			}
			system_energy_change+=0.5*(script->replica[old_bin].force*sqr(newdist) - script->replica[new_bin].force*sqr(olddist));

client->ptr+=sprintf(client->ptr," -- virtual pos %f energy %f\n",vrepop,0.5*(script->replica[old_bin].force*sqr(newdist) - script->replica[new_bin].force*sqr(olddist)));

		}

		system_energy_change+=new_cancellation-old_cancellation;
		system_energy_change*=var->beta;
		DRPE_change*=var->beta;
		total_energy_change=DRPE_change+system_energy_change;      // energies here are dimensionless
		client->ptr+=sprintf(client->ptr,"[system change, cancellation_change, DRPE change, total change]: %lf %lf %lf %lf\n",system_energy_change,new_cancellation-old_cancellation,DRPE_change,total_energy_change);
	}

	probability=exp(-total_energy_change);
	client->ptr+=sprintf(client->ptr,"Attempting monte carlo move ( %f to %f ) using these dimensionless quantities: [system change, DRPE change, total change, probability]: %lf %lf %lf %lf\n",old_coor,new_coor,system_energy_change,DRPE_change,total_energy_change,probability);
	
	random_number=drand48();
	if(probability>random_number){
		client->ptr+=sprintf(client->ptr,"Move accepted, new coordinate: %lf\n",new_coor);
		return(new_coor);
	}else{
		client->ptr+=sprintf(client->ptr,"Move rejected, keeping coordinate: %lf\n",old_coor);
		script->replica[replicaN].w=old_coor;
		return(old_coor);
	}
}

// determine where the given replica will move to next using the Boltzmann jumping scheme
// energy data contains the required data to do the move attempt, see comments below
float determine_new_replica_position_boltzmann_jump(struct client_struct *client, int replicaN, const float *energy_data, const struct script_struct *script, const struct server_variable_struct *var){
	int i;
	double new_coor, old_coor;
	double DRPE, system_energy, probability;
	float total_energy[script->Nreplicas];
	float energy_min;
	double random_number;
	double probability_sum;

	old_coor=script->replica[replicaN].w;

	energy_min=1e20;		
	for(i=0;i<script->Nreplicas;i++){
		if(i<var->min_running_replica || i>var->max_running_replica){
			total_energy[i]=1e20;
		}else{
			new_coor=script->replica[i].w_nominal;
			script->replica[replicaN].w=new_coor;
			DRPE=replica_potential(script);
			
			if(script->coordinate_type==Spatial){
				system_energy=energy_data[i];
				system_energy+=script->replica[i].cancellation_energy;
				total_energy[i]=var->beta*(DRPE+system_energy);
			}else if(script->coordinate_type==Temperature){
				system_energy=energy_data[0];   
				// for temperature space simulations, the first item in the array is the system energy
				system_energy+=script->replica[i].cancellation_energy;
				total_energy[i]=DRPE+new_coor*system_energy;
			}else{ // "Umbrella" simulation
				double exact_coordinate_position=energy_data[0];  
				//for umbrealla simulations, the first item in the array is the exact 
				//coordinate value of the entity that the umbrella acts on
				//CN says INCORRECT: system_energy=0.5*script->replica[replicaN].force*sqr(exact_coordinate_position-new_coor);
				system_energy=0.5*script->replica[i].force*sqr(exact_coordinate_position-new_coor);
				system_energy+=script->replica[i].cancellation_energy;
				total_energy[i]=var->beta*(DRPE+system_energy);
			}
		}
		if(total_energy[i]<energy_min) energy_min=total_energy[i];
	}

	for(i=0;i<script->Nreplicas;i++){
		total_energy[i]-=energy_min;
		fprintf(stderr,"energy after shifting: %lf\n",total_energy[i]); //##DEBUG
	}

	probability_sum=0.0;			
	for(i=0;i<script->Nreplicas;i++){
		probability=exp(-total_energy[i]);
		probability_sum+=probability;
		total_energy[i]=probability;   // total_energy now contains the probabilities
		fprintf(stderr,"probability[%d]=%lf\n",i,probability); //##DEBUG
	}

	for(i=0;i<script->Nreplicas;i++){
		total_energy[i]/=probability_sum;  // normalize the probability distribution
	}

	client->ptr+=sprintf(client->ptr,"Boltzmann jump possibilities:");
	for(i=0;i<script->Nreplicas;i++){
		if(total_energy[i]>0.001){
			client->ptr+=sprintf(client->ptr," (%d:%0.3f)",i,total_energy[i]);
		}
	}
	*(client->ptr++)='\n';

	probability_sum=0.0;
	random_number=drand48();
	for(i=0;i<script->Nreplicas;i++){
		probability=total_energy[i];
		probability_sum+=probability;
		if(probability_sum>=random_number) break;
	}

	fprintf(stderr,"random_number=%lf\n",random_number); //##DEBUG
	fprintf(stderr,"Probability sum: %0.20lf\n",probability_sum); //##DEBUG


	new_coor=script->replica[i].w_nominal;
	script->replica[replicaN].w=new_coor;

	double jump_distance=fabs(new_coor-old_coor);
	if(jump_distance<0.0001){
		client->ptr+=sprintf(client->ptr,"Boltzmann jump was unproductive\n");
	}else{
		client->ptr+=sprintf(client->ptr,"Boltzmann jump to new coordinate: %lf (distance of %lf)\n",new_coor,jump_distance);
	}

	return(new_coor);
}

// determine where the given replica will move to next using the Boltzmann jumping scheme in a continuous space
// energy data contains the required data to do the move attempt, see comments below
float determine_new_replica_position_continuous(struct client_struct *client, int replicaN, const float *energy_data, const struct script_struct *script, const struct server_variable_struct *var){
	int i,ii,j;
	int bin;
	int index;
	float left_w;
	float right_w;
	float division;
	double right_micro_w,left_micro_w;
	float old_coor, new_coor;
	double DRPE, system_energy, cancellation_energy, total_dimensionless_energy, probability;
	double random_number=drand48();
	double total_energy[script->Nreplicas*REPLICA_MICRODIVISIONS];

	old_coor=script->replica[replicaN].w;
	
	double probability_sum=0.0;
	float energy_min=1.0e20;

	for(i=0;i<script->Nreplicas;i++){
		ii=i;
		if(ii>=script->Nreplicas-1) ii--;
		left_w=script->replica[ii].w_nominal;
		right_w=script->replica[ii+1].w_nominal;
		division=(right_w-left_w)/REPLICA_MICRODIVISIONS;

		float left_cancellation_energy=script->replica[ii].cancellation_energy;
		float right_cancellation_energy=script->replica[ii+1].cancellation_energy;
		float cancellation_energy_division=(right_cancellation_energy-left_cancellation_energy)/REPLICA_MICRODIVISIONS;

		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;			

			fprintf(stderr,"new_coor=%f\n",new_coor); //##DEBUG

			if( (i<var->min_running_replica) || (i>var->max_running_replica) ){
				total_energy[index]=1e20;
			}else{
				new_coor=script->replica[i].w_nominal+(j-(REPLICA_MICRODIVISIONS-1)/2)*division;
				script->replica[replicaN].w=new_coor;
				DRPE=replica_potential(script);
				fprintf(stderr,"DRPE: %lf\n",DRPE); //##DEBUG			

				cancellation_energy=script->replica[i].cancellation_energy+(j-(REPLICA_MICRODIVISIONS-1)/2)*cancellation_energy_division;
				
				if(script->coordinate_type==Umbrella){
					double exact_coordinate_position=energy_data[0];  
					//for umbrealla simulations, the first item in the array is the exact 
					//coordinate value of the entity that the umbrella acts on
					//INCORRECT: system_energy=0.5*script->replica[replicaN].force*sqr(exact_coordinate_position-new_coor);
					//CN chose to correct the above line by simply applying the force to which the REPLICA_MICRODIVISION belongs, although interpolation is also possible
					system_energy=0.5*script->replica[i].force*sqr(exact_coordinate_position-new_coor);
					fprintf(stderr,"system_energy: %lf\tcancellation: %f\texact_coordinate_position %lf\tnew_coor %f\tscript->replica[i].force %f\n",system_energy, cancellation_energy,exact_coordinate_position,new_coor,script->replica[i].force); //##DEBUG
					total_dimensionless_energy=var->beta*(DRPE+system_energy+cancellation_energy);
				}else{
					system_energy=energy_data[0];  
					//for temperature space simulations, the first item in the array 
					//is the system energy
					fprintf(stderr,"system_energy: %lf   beta: %lf   cancellation: %f\n",system_energy, new_coor, cancellation_energy); //##DEBUG
					total_dimensionless_energy=DRPE+new_coor*(system_energy+cancellation_energy);
				}
				if(total_dimensionless_energy<energy_min) energy_min=total_dimensionless_energy;
				total_energy[index]=total_dimensionless_energy;
				fprintf(stderr,"total_energy: %lf\n",total_energy[index]); //##DEBUG
			}
		}
	}

	for(i=0;i<script->Nreplicas;i++){
		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;			
			total_energy[index]-=energy_min;
			fprintf(stderr,"total_energy after shifting: %lf\n",total_energy[index]); //##DEBUG
		}
	}

	for(i=0;i<script->Nreplicas;i++){
		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;			
			probability=exp(-total_energy[index]);
			total_energy[index]=probability;   // total_energy now contains the probabilities
		}
	}

	unsigned char first=1;
	for(i=0;i<script->Nreplicas;i++){
		ii=i;
		if(ii>=script->Nreplicas-1) ii--;
		left_w=script->replica[ii].w_nominal;
		right_w=script->replica[ii+1].w_nominal;
		division=(right_w-left_w)/REPLICA_MICRODIVISIONS;

		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;
			right_micro_w=script->replica[i].w_nominal+(j-(REPLICA_MICRODIVISIONS-1)/2)*division;
			if(!first){
				probability=(right_micro_w-left_micro_w)*(total_energy[index-1]+total_energy[index])/2;
				probability_sum+=probability;
				fprintf(stderr,"left_micro_w: %lf   right_micro_w: %lf   probability: %lf\n",left_micro_w,right_micro_w,probability); //##DEBUG
			}
			left_micro_w=right_micro_w;
			first=0;
		}
	}

	for(i=0;i<script->Nreplicas;i++){
		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;			
			total_energy[index]/=probability_sum;  // normalize the probability distribution
		}
	}

	fprintf(stderr,"random_number=%lf\n",random_number); //##DEBUG
	fprintf(stderr,"Probability sum: %0.20lf\n",probability_sum); //##DEBUG

	double probability_sum2=0.0;
	first=1;
	for(i=0;i<script->Nreplicas;i++){
		ii=i;
		if(ii>=script->Nreplicas-1) ii--;
		left_w=script->replica[ii].w_nominal;
		right_w=script->replica[ii+1].w_nominal;
		division=(right_w-left_w)/REPLICA_MICRODIVISIONS;

		for(j=0;j<REPLICA_MICRODIVISIONS;j++){
			index=i*REPLICA_MICRODIVISIONS+j;
			right_micro_w=script->replica[i].w_nominal+(j-(REPLICA_MICRODIVISIONS-1)/2)*division;
			if(!first){
				probability=(right_micro_w-left_micro_w)*(total_energy[index-1]+total_energy[index])/2;
				probability_sum2+=probability;
				fprintf(stderr,"left_micro_w: %lf   right_micro_w: %lf   probability: %lf  probability_sum2: %lf   p1: %lf   p2: %lf\n",left_micro_w,right_micro_w,probability,probability_sum2,total_energy[index-1],total_energy[index]); //##DEBUG
				if(probability_sum2>=random_number) goto done_loop;
			}
			left_micro_w=right_micro_w;
			first=0;
		}
	}
	done_loop:
	double fraction=random_number-(probability_sum2-probability);
	double m=((double)total_energy[index]-(double)total_energy[index-1])/(right_micro_w-left_micro_w);
	double y_intercept=(double)total_energy[index-1]-m*left_micro_w;
	//double x=11;
	//double R=m*0.5*(x*x-left_micro_w*left_micro_w)+y_intercept*(x-left_micro_w);
	//fprintf(stderr,"slope: %lf    y-intercept: %lf   R: %lf\n",m,y_intercept,R);
	
	double a=m*0.5;
	double b=y_intercept;
	double c=-m*0.5*left_micro_w*left_micro_w-y_intercept*left_micro_w-fraction;
	new_coor=(-b+sqrt(b*b-4*a*c))/(2*a);
	fprintf(stderr,"new_coor 1: %lf   new_coor 2: %lf\n",new_coor,(double)(-b-sqrt(b*b-4*a*c))/(2*a));  //##DEBUG

	fprintf(stderr,"random_number: %lf   probability: %lf   fraction: %f   new_coor: %f\n",random_number,probability,fraction,new_coor); //##DEBUG

	/* 
	client->ptr+=sprintf(client->ptr,"Boltzmann jump possibilities:");

	for(i=0;i<Nreplicas;i++)
		if(total_energy[i]>0.001)
			client->ptr+=sprintf(client->ptr," (%d:%0.3f)",i,total_energy[i]);
	*(client->ptr++)='\n';
	*/

	script->replica[replicaN].w=new_coor;

	double jump_distance=fabs(new_coor-old_coor);
	//if(jump_distance<0.0001)
	//	client->ptr+=sprintf(client->ptr,"Boltzmann jump was unproductive\n");
	//else
	client->ptr+=sprintf(client->ptr,"Boltzmann jump to new beta: %lf (distance of %lf)\n",new_coor,jump_distance);

	//fprintf(stderr,"i=%d\n",i); //##DEBUG
	//fprintf(stderr,"Updating replica %u: to w: %f\n",script->replicaN,script->replica[script->replicaN].w); //##DEBUG

	return(new_coor);
}

// runs the specified system command and appends a record to the log file
int execute(char *command, char quiet){
	int val;
	char message[MESSAGE_GLOBALVAR_LENGTH];
	if(quiet!='q'){
		sprintf(message,"Executing command: %s\n",command);
		append_log_entry(-1,message);
	}

	val=system(command);
	if(val!=0){
		sprintf(message,"ERROR: nonzero return value (%d) from system command: %s\n",val,command);
		append_log_entry(-1,message);
	}
	return val;
}

// submits a job to the queue system on the cluster
void start_queue_shell(bool conditional, const struct script_struct *script, struct server_variable_struct *var){
	char command[1000];
	int tries;
	int val;
	char message[MESSAGE_GLOBALVAR_LENGTH];

	pthread_mutex_lock(&queue_mutex);
	if( (conditional==false) || (var->Nreserved_queue_slots<ldiv(script->Nreplicas,script->Nsamesystem_uncoupled).quot) ){
		sprintf(message,"The current reserved queue slots count is %u ; starting another queue shell\n",var->Nreserved_queue_slots);
		append_log_entry(-1,message);
		sprintf(command,"%s/drsub %s",var->working_directory,var->working_directory); 
		tries=0;
		val=-1;
		while(val!=0 && tries++<2){
			//try twice for each submission in case there is an error
			val=execute(command,'v');
		}
		if(val==0){
			if(conditional==true){
				var->Nreserved_queue_slots++;
				if(var->Ncrashed_jobs>0) var->Ncrashed_jobs--;
			}
			var->nfailedsubinarow=0;
		}else{
			var->nfailedsubinarow++;
		}
	}
	pthread_mutex_unlock(&queue_mutex);
}

/*
 * This is no longer used since the drsub script is now responsible for all of this
// determines if the script for submitting jobs to the queue exists on disk
bool check_if_queue_shell_exists(void){
	int fd;
	int size;
	
	if( (fd=open(QUEUE_SHELL_FILENAME,O_RDONLY))==-1 ) return(false);
	size=lseek(fd,0,SEEK_END);
	close(fd);
	if(size>0) return(true);
	else return(false);
}
*/

// decrement our counter of the number of jobs currently in the queue
void decrement_Nreserved_queue_slots(struct server_variable_struct *var){
	pthread_mutex_lock(&queue_mutex);
	if(var->Nreserved_queue_slots>0)var->Nreserved_queue_slots--;
	pthread_mutex_unlock(&queue_mutex);
}

// determines the amount of disk space left on the disk and if there isn't enough, temporarily stop running simulations
//  when more disk space becomes available, simulations are allowed to run again
void check_if_disk_almost_full(struct server_variable_struct *var){
	struct statvfs data;
	double available_disk_GB;
	char message[MESSAGE_GLOBALVAR_LENGTH];

	if (statvfs("./", &data) == -1) error_quit("determination of disk usage failed");
	available_disk_GB=(double)data.f_bavail*(double)data.f_frsize/(1024*1024*1024);
	if(available_disk_GB<MIN_DISK_SPACE_TO_RUN){
		sprintf(message,"Not enough disk space (GB): %0.3lf, jobs will temporarily not run\n",available_disk_GB);
		if(var->simulation_status==Running) var->simulation_status=DiskAlmostFull;
	}else{
		sprintf(message,"Disk space available (GB): %0.3lf\n",available_disk_GB);
		if(var->simulation_status==DiskAlmostFull) var->simulation_status=Running;
	}
	append_log_entry(-1,message);
}

void check_if_finish_on_average(const struct script_struct *script, struct server_variable_struct *var){
	int i, s=0, d=0;

	for(i=0;i<script->Nreplicas;i++){
		//Note that the sequence numbers DO get properly incremented when NNI!=1 so this should work
		s+=script->replica[i].sequence_number;
		d+=script->replica[i].sampling_runs;
	}
	if(s>d){
		var->simulation_status=Finished;
		append_log_entry(-1,"Simulation completing because the AVERAGE simulation time was met.\n");
	}
}

// this is for two ligand simulations or simulations using the second ligand as the "funnel"
// this function allows the coordinate positions of the two ligands to be coupled.
// ie. when 'w' changes then this function will return the new value of 'w2'
double calculate_w2_from_w(double w, const struct script_struct *script){
	unsigned int i;
	double fraction;
	double w2;
	
	if(script->Nreplicas==1) return(script->replica[0].w2_nominal);
	if(script->Nligands<2) return(NAN);
	
	for(i=0;i<script->Nreplicas;i++){
		if(w < script->replica[i].w_nominal) break;
	}
	
	if(i<1) i=1;
	if(i>=script->Nreplicas) i=script->Nreplicas-1;

	fraction=(w-script->replica[i-1].w_nominal)/(script->replica[i].w_nominal-script->replica[i-1].w_nominal);

	w2=fraction*(script->replica[i].w2_nominal-script->replica[i-1].w2_nominal)+script->replica[i-1].w2_nominal;
	
	return(w2);
}

// suspend any replicas which have reached the desired sequence number
// however, only suspend replicas from the ends of the coordinate in and only do so if there is only one replica at the end discrete coordinate value
void check_termination_conditions(struct client_struct *client, const struct script_struct *script, struct server_variable_struct *var, const struct server_option_struct *opt){
	int i, bin;
	unsigned char test[script->Nreplicas];
	bool running=false;

	for(i=0;i<script->Nreplicas;i++) test[i]=0;
	for(i=0;i<script->Nreplicas;i++) test[find_bin_from_w(script->replica[i].w,script)]++;
	//for(i=0;i<script->Nreplicas;i++) printf("test[%d]=%hhu\n",i,test[i]); //##DEBUG
	for(i=0;i<script->Nreplicas;i++) if(script->replica[i].sample_count<script->replica[i].sampling_runs) { test[i]=0; running=true; }
	if(running==false) var->simulation_status=Finished;
	//printf("Set simulation_status=Finished in check_termination_conditions\n"); //##DEBUG
	//for(i=0;i<script->Nreplicas;i++) printf("test[%d]=%hhu\n",i,test[i]); //##DEBUG
	for(i=script->min_unsuspended_replica;i<=script->max_unsuspended_replica;i++) test[i]=0;

	//if(client!=NULL) for(i=0;i<script->Nreplicas;i++) client->ptr+=sprintf(client->ptr,"test[%d]=%hhu\n",i,test[i]); //##DEBUG
	
	for(i=0;i<script->Nreplicas;i++) if(test[i]!=1) break;
	var->min_running_replica=i;
	for(i=script->Nreplicas-1;i>=0;i--) if(test[i]!=1) break;
	var->max_running_replica=i;

	for(i=0;i<script->Nreplicas;i++){
		bin=find_bin_from_w(script->replica[i].w,script);
		if( (bin<var->min_running_replica) || (bin>var->max_running_replica) ){
			if(script->replica[i].status!='S') client->ptr+=sprintf(client->ptr,"Suspending replica %d\n", i);
			script->replica[i].status='S';
		}
	}

	if(opt->verbose!=0){
		client->ptr+=sprintf(client->ptr,"Checking suspend conditions: current min/max running replicas are: %d/%d ; counts at these boundaries are: %u/%u\n", var->min_running_replica, var->max_running_replica, script->replica[var->min_running_replica].sample_count, script->replica[var->max_running_replica].sample_count);
	}
}

// prints the number of clients that are currently connected
void print_number_of_connected_clients(const struct server_variable_struct *var){
	char message[MESSAGE_GLOBALVAR_LENGTH];
	pthread_mutex_lock(&replica_mutex);
	sprintf(message,"Number of currently connected clients is: %d\n",var->Nconnected_clients);
	pthread_mutex_unlock(&replica_mutex);
	append_log_entry(-1,message);
}

// change the count of the number of clients that are currently connected
void change_number_of_connected_clients(signed char x,struct server_variable_struct *var){
	pthread_mutex_lock(&replica_mutex);
	var->Nconnected_clients+=x;
	pthread_mutex_unlock(&replica_mutex);
}

// check if any jobs have timed out, if so, reset the flag so that they are resubmitted
void check_for_crash(int replicaNinput, const struct script_struct *script, struct server_variable_struct *var, struct node_struct *node){
	char message[MESSAGE_GLOBALVAR_LENGTH];
	int repstart, repend, t, mytime, replicaN;
	
	if(replicaNinput<0){
		//search through all replicas
		repstart=0;	
		repend=script->Nreplicas;
	}else{
		repstart=repend=replicaNinput;
	}
	pthread_mutex_lock(&replica_mutex);

	//this function canonly handle NNI-leader replicaN values. Other usage will cause an error
	mytime=time(NULL);
	for(replicaN=repstart;replicaN<repend;replicaN+=script->Nsamesystem_uncoupled){	
		if(script->replica[replicaN].status=='R'){
			t=mytime-script->replica[replicaN].last_activity_time;
			// if the replica has been gone too long:
			if(t>=script->job_timeout){
				// turn off the replica
				script->replica[replicaN].status='N';
				decrement_Nreserved_queue_slots(var);
				var->Ncrashed_jobs++;
				// update the node manager that this node no longer exists -- must be careful if it wakes up later and tries to get used
				if(script->replica[replicaN].nodeSlot<0){
					error_quit("Massive error in node management (A): attempt to release an unused node index in check_for_crash()\n");
				}
				releaseNode(&(node[script->replica[replicaN].nodeSlot]),&(script->replica[replicaN].nodeSlot));
				sprintf(message,"Replica %u: no client transaction activity for %d seconds, %u second timeout exceeded, restarting replica\n",replicaN,t,script->job_timeout);
				append_log_entry(-1,message);
			}
			
		}
	}

	pthread_mutex_unlock(&replica_mutex);
}

// if there are enough samples in each bin then compute the cancellation energy and activate the cancellation function
void conditionally_activate_energy_cancellation(struct script_struct *script, struct server_variable_struct *var){
	int i;
	double cancellation_energy[script->Nreplicas];
	char message[MESSAGE_GLOBALVAR_LENGTH];

	for(i=0;i<script->Nreplicas;i++){
		if(script->replica[i].cancellation_count<script->cancellation_threshold) break;
	}
	if(i<script->Nreplicas) return;

	for(i=0;i<script->Nreplicas;i++){
		script->replica[i].cancellation_energy=0.0;
	}

	if(script->coordinate_type==Spatial||script->coordinate_type==Umbrella){
		float dw;
		double F, old_F;
		double dE;
				
		for(unsigned char l=0;l<script->Nligands;l++){
			sprintf(message,"Forces have been averaged:\n");
			append_log_entry(-1,message);
			for(i=0;i<script->Nreplicas;i++){
				cancellation_energy[i]=script->replica[i].cancellation_accumulator[l]/script->replica[i].cancellation_count;
				sprintf(message,"N: %d\tFavg: %f\t(ligand): %d\n",i,cancellation_energy[i],l);
				append_log_entry(-1,message);
			}

			old_F=cancellation_energy[0];
			cancellation_energy[0]=0.0;
			for(i=1;i<script->Nreplicas;i++){
				if(l==0) dw=script->replica[i].w_nominal-script->replica[i-1].w_nominal;
				else     dw=script->replica[i].w2_nominal-script->replica[i-1].w2_nominal;
				F=cancellation_energy[i];
				dE=-(F+old_F)/2*dw;
				old_F=F;
				cancellation_energy[i]=cancellation_energy[i-1]+dE;
				script->replica[i].cancellation_energy-=cancellation_energy[i];
			}
		}
		//Note: only works with one ligand
		//Note: uses the final dw, which may not be correct (but it is for my alanine dipeptide run)
		if(script->circular_replica_coordinate){
			float temp;
			temp=-((script->replica[script->Nreplicas-1].cancellation_accumulator[0]/script->replica[script->Nreplicas-1].cancellation_count)+(script->replica[0].cancellation_accumulator[0]/script->replica[0].cancellation_count))/2*dw;
			circularCancelLow=-(0-temp);
			circularCancelHigh=-(-script->replica[script->Nreplicas-1].cancellation_energy+temp);
			//Now extend the error over each frame
			for(i=1;i<script->Nreplicas;i++){
				sprintf(message,"TESTER: %d was %f ",i,script->replica[i].cancellation_energy);
				append_log_entry(-1,message);
				script->replica[i].cancellation_energy-=(circularCancelHigh-0.0)*((float)i/(float)script->Nreplicas);
				sprintf(message,"and now is %f\n",script->replica[i].cancellation_energy);
				append_log_entry(-1,message);
			}
			//This is just to output it for interest
			temp=-((script->replica[script->Nreplicas-1].cancellation_accumulator[0]/script->replica[script->Nreplicas-1].cancellation_count)+(script->replica[0].cancellation_accumulator[0]/script->replica[0].cancellation_count))/2*dw;
			newCircularCancelLow=-(0-temp);
			newCircularCancelHigh=-(-script->replica[script->Nreplicas-1].cancellation_energy+temp);
		}
	}else{
		for(unsigned char l=0;l<script->Nligands;l++){
			for(i=0;i<script->Nreplicas;i++){
				cancellation_energy[i]=script->replica[i].cancellation_accumulator[l]/script->replica[i].cancellation_count;
				script->replica[i].cancellation_energy-=cancellation_energy[i];
			}
		}
	}

	var->energy_cancellation_status=Active;
	script->replica_potential_scalar1=script->replica_potential_scalar1_after_threshold;
	script->replica_potential_scalar2=script->replica_potential_scalar2_after_threshold;
}

// prints out a summary of the cancellation energies that are being used
void print_energy_cancellation_summary(const struct script_struct *script, const struct server_variable_struct *var){
	char message[MESSAGE_GLOBALVAR_LENGTH];
	char *summary;
	char *ptr;
	int i;

	summary=new char[script->Nreplicas*100+3*100];
	ptr=summary;
	ptr+=sprintf(ptr,"The following cancellation energies are now being used:\n");
	ptr+=sprintf(ptr,"-------------------------------------------------------\n");
	for(i=0;i<script->Nreplicas;i++){
		ptr+=sprintf(ptr,"coord: %f   energy: %f\n",script->replica[i].w_nominal,script->replica[i].cancellation_energy);
	}
	ptr+=sprintf(ptr,"LowCircularCoord: energy: %f\n",circularCancelLow);
	ptr+=sprintf(ptr,"HighCircularCoord: energy: %f\n",circularCancelHigh);
	ptr+=sprintf(ptr,"NEW LowCircularCoord: energy: %f\n",newCircularCancelLow);
	ptr+=sprintf(ptr,"NEW HighCircularCoord: energy: %f\n",newCircularCancelHigh);
	ptr+=sprintf(ptr,"-------------------------------------------------------\n");
	append_log_entry(-1,summary);
	delete[] summary;
	sprintf(message,"Using beta=%f\n",var->beta);
	append_log_entry(-1,message);
}

// read the specified number of byte from the given clients TCP/IP connection and put them into buffer
// if this fails then write an appropriate error message into the log file
unsigned char read_bytes_from_socket(struct client_struct *client, const char *failure_description, void *buff, unsigned int number_to_read){
	int Nread;
	int Nread_total=0;
	unsigned int n=number_to_read;  //##DEBUG
	
	//printf("number to read is %u\n",number_to_read);  //##DEBUG
	
	do{
		//printf("attempting to read %d bytes\n",number_to_read-Nread_total); //##DEBUG
		Nread=read(client->fd, (char*)buff+Nread_total, number_to_read-Nread_total);
		//printf("read returned %d bytes\n",Nread); //##DEBUG
		Nread_total+=Nread;
	}while( (Nread_total<number_to_read) && (Nread>0) );

	if(Nread<=0){
		client->ptr+=sprintf(client->ptr,"%s: failure reading %u bytes from the socket; read returned %d, errno is %d\n",failure_description, number_to_read, Nread, errno);
		return(0);
	}
	//printf("read a total of %d ************\n",Nread_total); //##DEBUG
	
	
	//for(int j=0;j<n;j++)                                 //##DEBUG
	//	printf("%hhu ",((unsigned char*)buff)[j]);   //##DEBUG
	//printf("\n");                                        //##DEBUG
	
	return(1);
}

// check to see if the server and client are suing the same protocol for communicating
unsigned char check_protocol_version(struct client_struct *client){
	unsigned int protocol_version;
	
	if(!read_bytes_from_socket(client, "Warning: cannot read protocol version", &protocol_version, PROTOCOL_VERSION_SIZE)) return(0);
	if(protocol_version != PROTOCOL_VERSION){
		client->ptr+=sprintf(client->ptr,"Warning: client and server protocol mismatch. Expecting: %u   Received: %u\n",PROTOCOL_VERSION,protocol_version);
		return(0);
	}
	return(1);
}

// receive a key and a command from the specified client
enum command_enum receive_key_and_command(struct client_struct *client){
	char buff[KEY_SIZE+COMMAND_SIZE];
	enum command_enum command;
	
	if(!read_bytes_from_socket(client, "Warning: reading key and command", buff, KEY_SIZE+COMMAND_SIZE)) return(InvalidCommand);
	printf("KEY: [%c%c%c%c%c%c%c%c%c%c%c] \n",buff[0],buff[1],buff[2],buff[3],buff[4],buff[5],buff[6],buff[7],buff[8],buff[9],buff[10]);  //##DEBUG
	printf("KEY: [%s]  command: %hhu\n",buff+KEY_LOCATION,buff[COMMAND_LOCATION]);  //##DEBUG

	command=(enum command_enum)buff[COMMAND_LOCATION];
	if(command<Exit){
		if( strncmp(buff+KEY_LOCATION,COMMAND_KEY,KEY_SIZE)!=0 ){
			client->ptr+=sprintf(client->ptr,"Warning: the key did not match\n");
			return(InvalidCommand);
		}
	}else if(command<InvalidCommand){
		if( strncmp(buff+KEY_LOCATION,COMMAND_KEY2,KEY_SIZE)!=0 ){
			client->ptr+=sprintf(client->ptr,"Warning: attempted a restricted function without the correct key!\n");
			return(InvalidCommand);
		}
	}else{
		return(InvalidCommand);
	}

	return(command);
}

// receive a file of the type specified by 'command' and place its contents into a buffer
unsigned char receive_file(struct client_struct *client, enum command_enum command, struct buffer_struct *data_buffer){
	char buffer[BUFFER_SIZE];
	int bytes_to_read;
	int bytes_left_to_read;
	int filename_size=-1;
	int file_fd=-1;
	
	printf("Receiving file...\n");  //##DEBUG
	
	if(!read_bytes_from_socket(client, "Getting size of file", &bytes_left_to_read, sizeof(bytes_left_to_read))) return(0);

	if( (command!=TakeThisFile) && (data_buffer->allocated_memory==0) )
	{
		data_buffer->data=new unsigned char[bytes_left_to_read];
		data_buffer->allocated_memory=bytes_left_to_read;
		printf("Allocating memory for file, size is: %d\n",bytes_left_to_read); //##DEBUG
	}

	//if(command==TakeRestartFile) printf("************* PRINTFING RESTART FILE FROM SOCKET ********************\n");	//##DEBUG

	//int fd;	                                                            //##DEBUG
	//if(command==TakeRestartFile){  //##DEBUG
	//	if( (fd=open("RESTART_IN2",O_WRONLY|O_CREAT|O_TRUNC,0644))==-1 ) exit(1);   //##DEBUG
	//}  //##DEBUG
	//char *c=(char *)data_buffer->data;
	
	while(bytes_left_to_read>0){
		//printf("%d bytes left to read\n",bytes_left_to_read);  //##DEBUG
		bytes_to_read=(bytes_left_to_read<BUFFER_SIZE) ? bytes_left_to_read : BUFFER_SIZE;
		if(!read_bytes_from_socket(client, "Getting file contents", buffer, bytes_to_read)) return(0);
		
		if( (command==TakeThisFile) && (file_fd==-1) ){
			int i;
			for(i=0;i<MAX_FILENAME_SIZE+1;i++)
				if(buffer[i]==0) break;   // find null terminating character
			if(i==MAX_FILENAME_SIZE+1) return(0);
			filename_size=i;
			//this for loop checks if there are any funny characters in the filename, 
			//of course to prevent malicious attacks by hackers
			for(i=0;i<filename_size;i++){ 
				if( ((buffer[i]<'a') || (buffer[i]>'z')) && ((buffer[i]<'A') || (buffer[i]>'Z')) && 
						(buffer[i]!='_') && (buffer[i]!='.') )
				{
					 return(0);
				}
			}
			printf("writing to file: [%s]\n",buffer);  //##DEBUG
			if( (file_fd=open(buffer,O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))==-1 ) return(0);
			fchmod(file_fd,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
			filename_size++;  // the null character is considered part of the filename
			bytes_left_to_read-=filename_size;
		}
		
		//if(command==TakeRestartFile)     //##DEBUG
		//	for(int j=0;j<bytes_to_read;j++) printf("%c",buffer[j]);   //##DEBUG
		//write(fd,ptr-bytes_to_read,bytes_to_read);   //##DEBUG
		//usleep(1000);
		
		if(command==TakeThisFile){
			if( write(file_fd,buffer+filename_size,bytes_to_read)!=bytes_to_read){
				client->ptr+=sprintf(client->ptr,"Error: cannot writing to file\n");
				return(0);
			}
			filename_size=0; 
		}else{
			data_buffer->data_size+=bytes_to_read;
			if( data_buffer->data_size > data_buffer->allocated_memory ){
				client->ptr+=sprintf(client->ptr,"Invalid data size\n");
				return(0);
			}
			//printf("Read in %d bytes ***************\n",bytes_to_read); //##DEBUG
			memcpy( data_buffer->data + data_buffer->data_size - bytes_to_read, buffer, bytes_to_read );
		}
		bytes_left_to_read-=bytes_to_read;
	}
	if(file_fd!=-1) close(file_fd);

	//if(command==TakeRestartFile) close(fd);  //##DEBUG

	//if(command==TakeRestartFile){  //##DEBUG
	//	int fd;	                                                            //##DEBUG
	//	if(command==TakeRestartFile) //##DEBUG
	//	{//##DEBUG
	//		if( (fd=open("RESTART_OUT",O_WRONLY|O_CREAT|O_TRUNC,0644))==-1 ) exit(1);	//##DEBUG
	//	}//##DEBUG
	//	write(fd, data_buffer->data, data_buffer->data_size );  //##DEBUG
	//	close(fd);  //##DEBUG
	//}  //##DEBUG
	
	return(1);
}

// receive the replica ID from the specified client
void receive_replica_ID(struct client_struct *client, int *replicaN, unsigned int *sequence_number, const struct server_option_struct *opt, const struct script_struct *script){
	struct ID_struct ID;
	
	printf("Attempting to receive replica ID\n");  //##DEBUG

	*replicaN=-1;
	if(!read_bytes_from_socket(client, "Getting repica ID", &ID, sizeof(ID))) return;

	if(strncmp(ID.title,opt->title,2)!=0){
		if(strncmp(ID.title,"**",2)==0){
			*replicaN=-2;
		}else{
			*replicaN=-1;
		}
		return;
	}

	if( (ID.replica_number>=0) && (ID.replica_number<script->Nreplicas) ) *replicaN=ID.replica_number;
	*sequence_number=ID.sequence_number;
}

// send the replica ID to the given client using the specified socket
unsigned char send_replica_ID(int socket_fd, int replicaN, unsigned int sequence_number, const struct server_option_struct *opt){
	struct ID_struct *ID;
	char buffer[KEY_SIZE+COMMAND_SIZE+sizeof(*ID)];
	
	ID=(struct ID_struct *)(buffer+KEY_SIZE+COMMAND_SIZE);

	memcpy(buffer,COMMAND_KEY,sizeof(COMMAND_KEY));
	buffer[KEY_SIZE]=ReplicaID;

	memcpy(ID->title,opt->title,4);
	ID->replica_number=replicaN;
	ID->sequence_number=sequence_number;
	
	return(write(socket_fd, buffer, sizeof(buffer))==sizeof(buffer));
}

// send the simulations parameters for the replica that is about to run to the given client
void send_simulation_parameters(struct client_struct *client, const struct replica_struct *rep, const struct script_struct *script, struct node_struct *node, int replicaN){
	//CN notes that previously rep was an actual copy so that it could be modified here and that would not change
	//the values Maintained by the other routines. Since CN now has to send an array of replica_struct, this will be
	//done by way of extra variables and no changes will be made to the replica_struct
	//replicaN neesd to be sent so that I can access the correct message structure...
	char buffer[KEY_SIZE+COMMAND_SIZE+500];
	unsigned int data_size;
	int random_seed=time(NULL);
	float *w2, *new_w;
	float w;
	int nni;
	int writecheck;

	w2=(float *)malloc(script->Nsamesystem_uncoupled*sizeof(float));
	if(w2==NULL){
		fprintf(stderr,"Error: unable to allocate memory for w2 in send_simulation_parameters().\n");
		exit(1);
	}
	new_w=(float *)malloc(script->Nsamesystem_uncoupled*sizeof(float));
	if(new_w==NULL){
		fprintf(stderr,"Error: unable to allocate memory for new_w in send_simulation_parameters().\n");
		exit(1);
	}

	if(script->coordinate_type==Temperature){
		//script->Nsamesystem_uncoupled is not compatable with temperature
		w=(float)1.0/(rep[0].w*BOLTZMANN_CONSTANT);
		w2[0]=NAN;
		client->ptr+=sprintf(client->ptr,"Sending simulation parameters: temperature=%0.1f; steps=%u; seed=%u\n", w, rep[0].sampling_steps, random_seed);
	}else{
		for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
			w2[nni]=calculate_w2_from_w(rep[nni].w,script);
			if(isnan(w2[nni])){
				client->ptr+=sprintf(client->ptr,"Sending simulation parameters: coord: %f; steps: %u; seed: %u\n", rep[nni].w, rep[nni].sampling_steps, random_seed);
			}else{
				client->ptr+=sprintf(client->ptr,"Sending simulation parameters: coord: %f; coord2: %f; steps: %u; seed: %u\n", rep[nni].w, w2[nni], rep[nni].sampling_steps, random_seed);
			}
		}
	}
	
	memcpy(buffer,COMMAND_KEY,sizeof(COMMAND_KEY));
	buffer[KEY_SIZE]=TakeSimulationParameters;
	data_size=KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int); // leave space for the file size

	if(!isnan(rep[0].force)){
		//CN doesn't think that a test for !isnan on all nni's is necessary
		data_size+=sprintf(buffer+data_size,"force");
		for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
			data_size+=sprintf(buffer+data_size," %f",rep[nni].force);
		}
		data_size+=sprintf(buffer+data_size,"\n");
	}

	data_size+=sprintf(buffer+data_size,"wref");
	if(script->coordinate_type==Temperature){
		data_size+=sprintf(buffer+data_size," %f\n",w);
	}else{
		for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
			data_size+=sprintf(buffer+data_size," %f",rep[nni].w);
		}
		data_size+=sprintf(buffer+data_size,"\n");
	}

	if(!isnan(w2[0])){
		data_size+=sprintf(buffer+data_size,"wref2");
		for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
			data_size+=sprintf(buffer+data_size," %f",w2[nni]);
		}
		data_size+=sprintf(buffer+data_size,"\n");
	}
	
	if(script->coordinate_type==Spatial && script->replica_move_type==MonteCarlo){
		data_size+=sprintf(buffer+data_size,"wrefchange");
		for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
			new_w[nni]=calculate_monte_carlo_move(rep[nni].w,script);
			data_size+=sprintf(buffer+data_size," %f",new_w[nni]);
		}
		data_size+=sprintf(buffer+data_size,"\n");

		w2[nni]=calculate_w2_from_w(new_w[nni],script);
		if(!isnan(w2[0])){
			data_size+=sprintf(buffer+data_size,"wrefchange2");
			for(nni=0;nni<script->Nsamesystem_uncoupled;nni++){
				data_size+=sprintf(buffer+data_size," %f",w2[nni]);
			}
			data_size+=sprintf(buffer+data_size,"\n");
		}
	}

	//sampling_steps can not be different for nni's in the same system
	data_size+=sprintf(buffer+data_size,"sampNsteps %u\n",rep[0].sampling_steps);
	data_size+=sprintf(buffer+data_size,"rnd %d\n",random_seed);

	//Send the special message -- this will usually mean that all of the other information is useless, but it doesn't hurt for now
	if(script->replica[replicaN].nodeSlot<0){
		error_quit("Massive error in node management (B): script->replica[replicaN].nodeSlot < 0 in send_simulation_parameters() before usage\n");
	}
		
	if(node[script->replica[replicaN].nodeSlot].messageWaitingIndicator){
		data_size+=sprintf(buffer+data_size,"MESSAGE %s\n",node[script->replica[replicaN].nodeSlot].serverMessageForClient);
		pthread_mutex_lock(&replica_mutex);
		node[script->replica[replicaN].nodeSlot].messageWaitingIndicator=false;
		pthread_mutex_unlock(&replica_mutex);
	}

	*(unsigned int*)(buffer+KEY_SIZE+COMMAND_SIZE)=data_size-(KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int));

	//writecheck=write(client->fd, buffer, data_size);
	// should catch like above, but I'll leave the warning as it should really be fixed
	write(client->fd, buffer, data_size);
	free(w2);free(new_w);
}

// send the restart file for the replica that is about to run to the given client
void send_restart_file(int socket_fd, struct buffer_struct restart){
	char buffer[KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int)];
	int block_size;
	int ptr;
	char b[BUFFER_SIZE];
	int writecheck;

	printf("Sending restart data, the size of the file is: %u\n",restart.data_size);  //##DEBUG

	memcpy(buffer,COMMAND_KEY,sizeof(COMMAND_KEY));
	buffer[KEY_SIZE]=TakeRestartFile;
	*(unsigned int*)(buffer+KEY_SIZE+COMMAND_SIZE)=restart.data_size;

	//writecheck=write(socket_fd, buffer, sizeof(buffer));
	// should catch like above, but I'll elave the warning as it should really be fixed
	write(socket_fd, buffer, sizeof(buffer));
	
	printf("Sending TakeRestartFileHeader: %d bytes\n",(int)sizeof(buffer)); 		//##DEBUG
	//int fd;	                                                        	    	//##DEBUG
	//if( (fd=open("RESTART_IN",O_RDONLY))==-1 ) exit(1);					//##DEBUG
	//printf("************** SENDING RESTART FILE FROM DATABASE *****************\n"); 	//##DEBUG

	block_size=BUFFER_SIZE;
	for(ptr=0;ptr<restart.data_size;ptr+=block_size)
	{
		if(ptr+block_size>restart.data_size) block_size=restart.data_size-ptr;
		//printf("Sending data chunk: %d bytes\n",block_size);             	//##DEBUG
		//for(int j=0;j<block_size;j++)                                    	//##DEBUG
		//	printf("%hhu ",((unsigned char*)restart.data+ptr)[j]);   	//##DEBUG
		//printf("\n");                                                    	//##DEBUG
		//read(fd,b,block_size);
		//write(socket_fd, b, block_size);

		//writecheck=write(socket_fd, (char*)restart.data+ptr, block_size);
		// should catch like above, but I'll elave the warning as it should really be fixed
		write(socket_fd, (char*)restart.data+ptr, block_size);
		
		//for(int j=0;j<block_size;j++) printf("%c",*((char*)restart.data+ptr+j));  //##DEBUG
		
		//write(fd,(char*)restart.data+ptr, block_size);	//##DEBUG
	}
	//close(fd);
}

// This function runs as a thread for each client that connects
// Handles all interaction with the client:
// receives replica ID, move energy data, sample data, coordinate data, restart file
// checks the integrity of all the files that are received
// calls the appropriate routines to make a move in the coordinate for the replica
// sends new replica ID, simulations parameters as well as possible a restart file to the client
// updates current restart file and averaged coordainte data
// saves sample data to a database

void *client_interaction(struct client_bundle *B){
	enum client_status_enum {Communicating, NewNode, ReplicaFinished, Error};
	//char buff[KEY_SIZE+COMMAND_SIZE];
	enum command_enum command;		
	int *replicaN;
	int *old_replicaN;
	int replica_time_running=0;
	int node_time_running=0;
	int *bin;
	struct replica_struct *current_replica;
	struct buffer_struct *energy;
	struct buffer_struct *sample;
	struct buffer_struct **additional_data; 
	struct buffer_struct *coordinate;
	struct buffer_struct *tcs; //time client started
	struct buffer_struct *jid;
	enum client_status_enum client_status;
	struct buffer_struct *sampleOrAdditional;
	int NsampleOrAdditional=-1; // -1 case is for sample, rest are for additional_data
	int *save_sample_data_replicaN;
	unsigned int *save_sample_data_sequence_number;
	float *save_sample_data_w;
	unsigned int ai,i,nni;
	double average;
	unsigned char l;
	int node_just_reanimated;
	int dumpnode;
	bool allow_write_record=true;
	//client_status=NewNode is overloaded so it can't be used to inhibit force database access for a node that is really new 
        // For a neally new node, there is nothing to write anyway, but also it won't have a message and so might lead to a segfault
	bool newConnection=false;
	bool unexpectedClient=false;

	//CN wonders if there is a way to avoid allocating this memory every time.
	energy=(buffer_struct *)malloc(B->script->Nsamesystem_uncoupled*sizeof(buffer_struct));
	sample=(buffer_struct *)malloc(B->script->Nsamesystem_uncoupled*sizeof(buffer_struct));
	additional_data=(buffer_struct **)malloc(B->script->Nsamesystem_uncoupled*sizeof(buffer_struct *));
	coordinate=(buffer_struct *)malloc(B->script->Nsamesystem_uncoupled*sizeof(buffer_struct));
	tcs=(buffer_struct *)malloc(sizeof(buffer_struct));
	jid=(buffer_struct *)malloc(sizeof(buffer_struct));
	current_replica=(replica_struct *)malloc(B->script->Nsamesystem_uncoupled*sizeof(replica_struct));
	replicaN=(int *)malloc(B->script->Nsamesystem_uncoupled*sizeof(int));
	old_replicaN=(int *)malloc(B->script->Nsamesystem_uncoupled*sizeof(int));
	save_sample_data_replicaN=(int *)malloc(B->script->Nsamesystem_uncoupled*sizeof(int));
	save_sample_data_sequence_number=(unsigned int *)malloc(B->script->Nsamesystem_uncoupled*sizeof(unsigned int));
	save_sample_data_w=(float *)malloc(B->script->Nsamesystem_uncoupled*sizeof(float));
	bin=(int *)malloc(B->script->Nsamesystem_uncoupled*sizeof(int));
	if(energy==NULL||sample==NULL||additional_data==NULL||coordinate==NULL||tcs==NULL||jid==NULL||current_replica==NULL||replicaN==NULL||old_replicaN==NULL||save_sample_data_replicaN==NULL||save_sample_data_sequence_number==NULL||save_sample_data_w==NULL||bin==NULL){
		fprintf(stderr,"Error: unable to allocate memory for energy, sample, additional_data, coordinate, tcs, jid, current_replica, replicaN, old_replicaN, save_sample_data_replicaN, save_sample_data_sequence_number, save_sample_data_w, or bin during client_interaction()\n");
		exit(1);
	}
	for(i=0;i<B->script->Nsamesystem_uncoupled;i++){
		energy[i].data=sample[i].data=coordinate[i].data=NULL;
		energy[i].data_size=sample[i].data_size=coordinate[i].data_size=0;
		energy[i].allocated_memory=sample[i].allocated_memory=coordinate[i].allocated_memory=0;
		additional_data[i]=(buffer_struct *)malloc(B->script->Nadditional_data*sizeof(buffer_struct));
		if(additional_data[i]==NULL){
			fprintf(stderr,"Error: unable to allocate memory for additional_data[i] in client_interaction()\n");
			exit(1);
		}
		for(ai=0; ai<B->script->Nadditional_data; ai++){
			additional_data[i][ai].data=NULL;
			additional_data[i][ai].data_size=additional_data[i][ai].allocated_memory=0;
		}
		current_replica[i].restart.data=NULL;
		current_replica[i].restart.data_size=current_replica[i].restart.allocated_memory=0;
		replicaN[i]=old_replicaN[i]=save_sample_data_replicaN[i]=-1;
	}
	tcs->data=NULL;
	tcs->data_size=tcs->allocated_memory=0;
	jid->data=NULL;
	jid->data_size=jid->allocated_memory=0;

	printf("thread successfully started\n"); //##DEBUG

	client_status=Communicating;
	if(!check_protocol_version(B->client)) client_status=Error;

	nni=0;
	while(client_status==Communicating){
		client_status=Error;

		if( (command=receive_key_and_command(B->client))==InvalidCommand ) break;

		switch(command)
		{
		case ReplicaID:
			printf("ReplicaID command received\n");  //##DEBUG
			if(nni!=0){
				B->client->ptr+=sprintf(B->client->ptr,"Programming error: replica ID received on other than the first NNI\n");
				break;
			}
			receive_replica_ID(B->client, &replicaN[0], &current_replica[0].sequence_number,B->opt,B->script);
			old_replicaN[0]=replicaN[0];
			if(replicaN[0]>=0){
				client_status=Communicating;
				B->client->ptr+=sprintf(B->client->ptr,"Replica ID received: %2sw%d.%u\n",B->opt->title,replicaN[nni],current_replica[nni].sequence_number);
				//ReplicaID is only sent for the first nni so we need to fill in some values that we know to be true
				for(i=1;i<B->script->Nsamesystem_uncoupled;i++){
					old_replicaN[i]=replicaN[i]=replicaN[0]+i; // we enforce that they are consecutively numbered (rep# not W)
					current_replica[i].sequence_number=current_replica[0].sequence_number;
				}
			}
			else if(replicaN[nni]==-2){
				B->client->ptr+=sprintf(B->client->ptr,"A Node has just been occupied\n");
				decrement_Nreserved_queue_slots(B->var);
				client_status=NewNode;
				newConnection=true;
			}else{
				B->client->ptr+=sprintf(B->client->ptr,"Warning: invalid replica ID received\n");
			}
			break;
		case TakeTCS:
			printf("TakeTCS command received\n");  //##DEBUG
			if(receive_file(B->client, command, tcs)){
				client_status=Communicating;
				// tcs gets sent around as a float since there are already routines for sending floats
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"This client started running (in seconds) : %d\n",(int)*((float*)(tcs->data)));  
				}
			}
			break;
		case TakeJID:
			printf("TakeJID command received\n");  //##DEBUG
			if(receive_file(B->client, command, jid)){
				client_status=Communicating;
				// jid gets sent around as a float since there are already routines for sending floats
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"This client is running as JOB: %d\n",(int)*((float*)(jid->data)));
				}
			}
			break;
		case TakeThisFile:
			printf("TakeThisFile command received\n");  //##DEBUG                      
			if(receive_file(B->client, command, NULL)){
				client_status=Communicating;
				if(B->opt->verbose){	
					B->client->ptr+=sprintf(B->client->ptr,"A file was successfully received and written to disk\n");
				}
			}
			break;
		case TakeMoveEnergyData:
			printf("TakeMoveEnergyData command received\n");  //##DEBUG
			if(B->script->replica_move_type==NoMoves){
				B->client->ptr+=sprintf(B->client->ptr,"Client attempting to send energy data but we aren't attempting moves along the coordinate\n");
				break;
			}
			if(receive_file(B->client, command, &energy[nni])){
				client_status=Communicating;
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"Energy data was successfully received\n");
				}
			}
			break;
		case TakeSampleData:
			printf("TakeSampleData command received\n");  //##DEBUG
			if(!B->script->need_sample_data){
				B->client->ptr+=sprintf(B->client->ptr,"Client attempting to send sample data but we don't want it\n");
				break;
			}
			//need typecast to compare to NsampleOrAdditional (can be = -1)
			if(NsampleOrAdditional>=(int)B->script->Nadditional_data){
				B->client->ptr+=sprintf(B->client->ptr,"Client attempting to send additional data but we don't want it\n");
				fprintf(stderr,"NsampleOrAdditional=%d and Nadditional_data=%d\n",NsampleOrAdditional,B->script->Nadditional_data);
				break;
			}
			if(NsampleOrAdditional==-1){
				//sample data must come first
				sampleOrAdditional=&sample[nni];
			}else{
				sampleOrAdditional=&additional_data[nni][NsampleOrAdditional];
			}
			if( receive_file(B->client, command, sampleOrAdditional) && 
				(sampleOrAdditional->data_size==sampleOrAdditional->allocated_memory) )
			{
				client_status=Communicating;
				if(B->opt->verbose){
					if(NsampleOrAdditional==-1){
						B->client->ptr+=sprintf(B->client->ptr,"A sample data file was successfully received\n");
					}else{
						B->client->ptr+=sprintf(B->client->ptr,"An additional data file was successfully received\n");
					}
				}
				//for(i=0;i<B->script->Nsamples_per_run;i++){
				//	B->client->ptr+=sprintf(B->client->ptr,"CN LOOKEY %f\n",((float *)sampleOrAdditional->data)[i]);
				//}
	
			}
			++NsampleOrAdditional;
			break;
		case TakeCoordinateData:
			printf("TakeCoordinateData command received\n");  //##DEBUG
			if(!B->script->need_coordinate_data){
				B->client->ptr+=sprintf(B->client->ptr,"Client attempting to send coordinate data but we don't want it\n");
				break;
			}
			if( receive_file(B->client, command, &coordinate[nni]) && 
				(coordinate[nni].data_size==coordinate[nni].allocated_memory) )
			{
				client_status=Communicating;
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"A coordinate data file was successfully received\n");
				}
			}
			break;
		case TakeRestartFile:
		case NextNonInteracting:
			printf("TakeRestartFile or NextNonInteracting command received\n");  //##DEBUG
			if(command==TakeRestartFile && !receive_file(B->client, command, &current_replica[nni].restart)){
				break;
			}
			if(B->opt->verbose){
				B->client->ptr+=sprintf(B->client->ptr,"A restart file or indication of NextNonInteracting was successfully received\n");
			}
			if(replicaN[nni]<0) break;
			if(command==TakeRestartFile && !check_restart_file_integrity(B->client, current_replica[nni].restart)) break;
			if( B->script->need_sample_data && !check_sample_file_integrity(B->client, sample[nni], 0, B->script) ) break;
			if( B->script->need_sample_data && B->script->Nadditional_data>0){
				if((int)B->script->Nadditional_data!=NsampleOrAdditional){
					B->client->ptr+=sprintf(B->client->ptr,"Expected %d additional data files, but only found %d\n",B->script->Nadditional_data,NsampleOrAdditional);
					break;
				}
				for(ai=0; ai<B->script->Nadditional_data; ai++){
					if(!check_sample_file_integrity(B->client, additional_data[nni][ai], 1, B->script)) break;
				}
			}
			if(B->script->need_coordinate_data && !check_coordinate_file_integrity(B->client,coordinate[nni],B->var)) break;
			if(B->script->replica_move_type!=NoMoves && !check_energy_file_integrity(B->client,energy[nni],B->script)) break;
			if(nni+1==B->script->Nsamesystem_uncoupled){
				client_status=ReplicaFinished;
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"All integrity checks passed\n");
				}
			}else{
				client_status=Communicating;
				nni++;
				NsampleOrAdditional=-1;
			}
			NsampleOrAdditional=-1;
			break;
		case Exit:
			printf("Exit command received\n");  //##DEBUG
			B->var->simulation_status=Finished;
			B->client->ptr+=sprintf(B->client->ptr,"Exit command was received; DR_server will exit\n");
			break;
		case Snapshot:
			B->var->save_snapshot_now=true;
			B->client->ptr+=sprintf(B->client->ptr,"Snapshot command was received; DR_server will write out a snapshot\n");
			break;
		default:
			B->client->ptr+=sprintf(B->client->ptr,"Warning: an unknown command was received from the client\n");
		}
	};

	if(replicaN[nni]!=-2 && nni+1!=B->script->Nsamesystem_uncoupled){
		// When replicaN[nni]==-2, it is a new replica and only sends the replicaID
		B->client->ptr+=sprintf(B->client->ptr,"Only received information for %d noninteracting copies from the client but there should be %d.\n",nni+1,B->script->Nsamesystem_uncoupled);
		client_status=Error;
	}
	//nni will now become a general purpose index of B->script->Nsamesystem_uncoupled in for loops

	//fprintf(stderr,"Trying to get a lock after first comm round\n");fflush(stderr);    //CN FIND PROBLEM 
	pthread_mutex_lock(&replica_mutex);
	node_just_reanimated=0;
	switch(client_status)
	{
	case ReplicaFinished:
		for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
			if(B->script->replica[replicaN[nni]].status!='R'){
				B->client->ptr+=sprintf(B->client->ptr,"Warning: this replica is apparently not running, what's going on?\n");
				client_status=Error;
				unexpectedClient=true;
				break; //only breaks for loop
			}
			if(current_replica[nni].sequence_number!=B->script->replica[replicaN[nni]].sequence_number){
				B->client->ptr+=sprintf(B->client->ptr,"Warning: the sequence number for this job is invalid\n");
				client_status=Error;
				unexpectedClient=true;
				break; //only breaks for loop
			}
		}
		if(client_status==Error){
			if(!B->script->allow_requeue){
				break; //breaks switch
			}
			//leave unexpectedClient=true because I don't want to write any of this to the database
			client_status=NewNode;
			B->client->ptr+=sprintf(B->client->ptr,"This run ALLOWS REQUEUEING of resurected Nodes, this will occur now...\n");
			// must set replicaN[0]=-1 or else it will double submit the job that it had previously
			replicaN[0]=-1;
			node_just_reanimated=1;
			//Allow flow through to case NewNode if allow_requeue is set. However, since there was some type
			//of error, we don't want to do the rest of the checking and updating.
		}else{
			for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
				save_sample_data_replicaN[nni]=replicaN[nni];
				save_sample_data_sequence_number[nni]=B->script->replica[replicaN[nni]].sequence_number;
				save_sample_data_w[nni]=B->script->replica[replicaN[nni]].w;
				bin[nni]=find_bin_from_w(B->script->replica[replicaN[nni]].w,B->script);
				if(B->var->energy_cancellation_status==Pending){
					if(B->script->coordinate_type==Spatial||B->script->coordinate_type==Umbrella){
						for(l=0;l<B->script->Nligands;l++){
							average=0.0;
							//CN worries about B->script->Nsamples_per_run being only a single value
							for(i=0;i<B->script->Nsamples_per_run;i++){
								average+=((float*)sample[nni].data)[i*B->script->Nligands+l];
							}
							average/=B->script->Nsamples_per_run;
							B->script->replica[bin[nni]].cancellation_accumulator[l]+=average;
						}
					}else{
						B->script->replica[bin[nni]].cancellation_accumulator[0]+=*((float*)(energy[nni].data));
					}
					B->script->replica[bin[nni]].cancellation_count++;
				}
			}
			//Intentionally moved this outside of the above for-loop
			if(B->var->energy_cancellation_status==Pending){
				conditionally_activate_energy_cancellation(B->script,B->var);
				if(B->var->energy_cancellation_status==Active){
					B->client->ptr+=sprintf(B->client->ptr,"The energy cancellation feature has been activated; please stand by for a summary\n");
				}
			}
	
			//only commit restart file for first nni
			commit_restart_file(replicaN[0],&current_replica[0].restart,B->script);
			printf("Just out of commit pointer is: %p\n",current_replica[0].restart.data);  //##DEBUG
		
			for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){	
				if(coordinate[nni].data!=NULL) commit_coordinate_data(replicaN[nni],bin[nni],&coordinate[nni],B->script,B->var);
				if(B->opt->verbose){
					if(B->script->coordinate_type==Temperature){
						B->client->ptr+=sprintf(B->client->ptr,"Bin %d at temperature %0.1f will be incremented\n", bin[nni], (float)1.0/(B->script->replica[replicaN[nni]].w*BOLTZMANN_CONSTANT));
					}else{
						B->client->ptr+=sprintf(B->client->ptr,"Bin %d at w %f will be incremented\n", bin[nni], B->script->replica[replicaN[nni]].w);
					}
				}
				B->script->replica[bin[nni]].sample_count++;
				B->script->replica[replicaN[nni]].sequence_number++;
				if(B->opt->verbose){
					B->client->ptr+=sprintf(B->client->ptr,"Incrementing sequence number of replica %d to %hu\n", replicaN[nni], B->script->replica[replicaN[nni]].sequence_number);
				}
				if(B->script->replica_move_type==MonteCarlo||B->script->replica_move_type==vRE){
					determine_new_replica_position_monte_carlo_or_vre(B->client, replicaN[nni], (float*)energy[nni].data,B->script,B->var);
				}else if(B->script->replica_move_type==BoltzmannJumping){
					determine_new_replica_position_boltzmann_jump(B->client, replicaN[nni], (float*)energy[nni].data,B->script,B->var);
				}else if(B->script->replica_move_type==Continuous){
					determine_new_replica_position_continuous(B->client, replicaN[nni], (float*)energy[nni].data,B->script,B->var);
				}
			}

			if(B->script->replica[replicaN[0]].nodeSlot<0){
				error_quit("Massive error in node management (C): B->script->replica[replicaN[0]].nodeSlot < 0 in client_interaction() before usage\n");
			}
			
			node_time_running=time(NULL)-B->node[B->script->replica[replicaN[0]].nodeSlot].start_time;
			if(B->opt->verbose){
				B->client->ptr+=sprintf(B->client->ptr,"Node %s has been running for (seconds): %d\n",B->node[B->script->replica[replicaN[0]].nodeSlot].ip, node_time_running);
			}
		
			check_termination_conditions(B->client,B->script,B->var,B->opt);
			for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){                                     //##DEBUG
				printf("the simulation status is %c\n",B->script->replica[replicaN[nni]].status);  //##DEBUG
			}                                                                                          //##DEBUG

			if(node_time_running<0){
				B->client->ptr+=sprintf(B->client->ptr,"WARNING: negative runtime. There is some lack of synchronization in your clocks!\n");
			}else if(node_time_running>=B->script->node_time){
				//Release the node -- cleanup will be done based on client_status=Error a bit later
				if(B->script->replica[replicaN[0]].status=='R'){
					B->script->replica[replicaN[0]].status='N';
				}
				B->client->ptr+=sprintf(B->client->ptr,"Time limit on Node exceeded; Node will be freed\n");
				client_status=Error;
				break;
			}
			replica_time_running=time(NULL)-B->script->replica[replicaN[0]].start_time_on_current_node;
			if(replica_time_running>=B->script->replica_change_time && B->script->replica[replicaN[0]].status=='R'){
				//just allow the replicas to exchange... it could be turned on again
				disconnectNodeFromReplica(&(B->script->replica[replicaN[0]].nodeSlot));
				B->script->replica[replicaN[0]].status='N';
			}
		} // ends long else -- there should not be anything between this and case NewNode
	case NewNode:
		//fprintf(stderr,"finding replica to run\n");fflush(stderr);             //CN FIND PROBLEM 
		find_replica_to_run(&replicaN[0],B->script);       // Determine the system to run based on the first nni
		//fprintf(stderr,"find_replica returned %d\n",replicaN[0]);fflush(stderr); //CN FIND PROBLEM 

		if(replicaN[0]<0 && B->var->simulation_status==Running && B->script->cycleClients>-0.001){
			//There are no new jobs to run, but it is possible that this Node should be picked up and another one dropped for efficiency
			// drop_one_old_node() must not return until that one has actually come back and been turned off (comm done).
			// NOTE: there is currently a mutex on the replicas...
			//fprintf(stderr,"trying to drop\n");fflush(stderr);    //CN FIND PROBLEM 
			dumpnode=drop_one_old_node(B->client,B->script,B->node);
			//fprintf(stderr,"finished trying to drop\n");fflush(stderr); //CN FIND PROBLEM 
			if(dumpnode>=0){
				B->client->ptr+=sprintf(B->client->ptr,"DUMP (part 1): Shut down node [%d] (IP=%s) to make room for a new job\n",dumpnode,B->node[dumpnode].ip);
				find_replica_to_run(&replicaN[0],B->script);
				if(replicaN[0]>=0){
					B->client->ptr+=sprintf(B->client->ptr,"DUMP (part 2): A new job should start using replica %d\n",replicaN[0]);
				}else{
					B->client->ptr+=sprintf(B->client->ptr,"DUMP (part 2): Error: could not find a new job to run after dumping! (no action taken on error)\n");
				}
			}else{
				//if dumpnode<0 then the decision was to not dump any Nodes
				B->client->ptr+=sprintf(B->client->ptr,"DUMP (part 0): A Node dump was considered, but was not enacted as it would not be efficient at this time.\n");
			}
		}

		if(replicaN[0]<0 || B->var->simulation_status!=Running){
			B->client->ptr+=sprintf(B->client->ptr,"No jobs to run at this time\n");
			if(B->var->simulation_status!=Running){
				B->client->ptr+=sprintf(B->client->ptr,"  because simulation_status!=Running\n");
				if(B->var->simulation_status==DiskAlmostFull){
					B->client->ptr+=sprintf(B->client->ptr,"    because simulation_status==DiskAlmostFull\n");
				}
			}
			client_status=Error;
			break;
		}

		// Getting to here means that a node was just occupied and will be sent a job
		// there is already a mutex lock on replicas
		// YUPYUPYUP
		int myNewNodeSlot;

		/*
		B->client->ptr+=sprintf(B->client->ptr,"NODE: what node has ip %s ...\n",B->client->ip);
		displayNodes(B->script,B->node);
		*/

		myNewNodeSlot=findNodeByIP(B->script, B->node, B->client->ip);

		if(myNewNodeSlot>=0){
			B->client->ptr+=sprintf(B->client->ptr,"NODE: found [%d] already active\n",myNewNodeSlot);
			connectNodeToReplica(&(B->script->replica[replicaN[0]].nodeSlot),myNewNodeSlot);
		}else{
			B->client->ptr+=sprintf(B->client->ptr,"NODE: didn't find one already active\n");
			// it's a new node
			myNewNodeSlot = findInactiveNodeSlot(B->script, B->node);
			if(myNewNodeSlot<0){
				B->client->ptr+=sprintf(B->client->ptr,"NODE: can't find an inactive, dropping this client\n");
				client_status=Error;
				break;
			}
			B->client->ptr+=sprintf(B->client->ptr,"NODE: picked up inactive node %d\n",myNewNodeSlot);
			obtainNode(&(B->node[myNewNodeSlot]),&(B->script->replica[replicaN[0]].nodeSlot),B->client->ip,myNewNodeSlot,(int)*((float*)(tcs->data)));
		}

		for(nni=1;nni<B->script->Nsamesystem_uncoupled;nni++){
			replicaN[nni]=replicaN[0]+nni;
		}

		for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
			B->script->replica[replicaN[nni]].status='R';
			current_replica[nni]=B->script->replica[replicaN[nni]];
			current_replica[nni].restart.data=NULL;
			current_replica[nni].restart.data_size = current_replica[nni].restart.allocated_memory = 0;
		}

		//Again, this next check is only based on the first nni and we only use the restart for this first nni
		if(node_just_reanimated || (replicaN[0]!=old_replicaN[0] && current_replica[0].sequence_number>0)){
			copy_restart_data(&current_replica[0].restart,&(B->script->replica[replicaN[0]].restart));
		}
		break;
	} // end swtich

	if(client_status==Error && !unexpectedClient){
		//Do this before the replica mutex ever comes off
		//unexpected clients should simply exit
		if(replicaN[0]>=0 && B->script->replica[replicaN[0]].nodeSlot>=0){
			// but don't release a node that was never taken
			if(B->node[B->script->replica[replicaN[0]].nodeSlot].messageWaitingIndicator){
				// There are only messages left to run the mobile server. Therefore the forcedatabase is already closed
				// must turn off the messageWaitingIndicator since this client will never go through send_simulation_parameters()
				allow_write_record=false;
				B->node[B->script->replica[replicaN[0]].nodeSlot].messageWaitingIndicator=false;
			}
			releaseNode(&(B->node[B->script->replica[replicaN[0]].nodeSlot]),&(B->script->replica[replicaN[0]].nodeSlot));
                        if(B->script->replica[replicaN[0]].status=='R'){
                                B->script->replica[replicaN[0]].status='N';
                        }
		}
	}
	
	pthread_mutex_unlock(&replica_mutex);

	for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
		printf("freeing energy: pointer before is: %p\n",energy[nni].data);  //##DEBUG
		delete[] energy[nni].data;
		printf("freeing energy: pointer after is: %p\n",energy[nni].data);  //##DEBUG

		printf("freeing coordinate: pointer before is: %p\n",coordinate[nni].data);  //##DEBUG
		delete[] coordinate[nni].data;
		printf("freeing coordinate: pointer after is: %p\n",coordinate[nni].data);  //##DEBUG
	}

	if(client_status!=Error){
		//only send the replica ID of the first nni
		send_replica_ID(B->client->fd, replicaN[0], current_replica[0].sequence_number,B->opt);
		B->client->ptr+=sprintf(B->client->ptr,"Replica ID sent: %2sw%d.%u",B->opt->title,replicaN[0],current_replica[0].sequence_number);
		//only send the restart of the first nni
		if(current_replica[0].restart.data!=NULL){
			send_restart_file(B->client->fd, current_replica[0].restart);
			B->client->ptr+=sprintf(B->client->ptr,", restart file sent");
		}
		
		B->client->ptr[0]='\n'; B->client->ptr++;
		for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
			B->script->replica[replicaN[nni]].start_time_on_current_node=B->script->replica[replicaN[nni]].last_activity_time=time(NULL);
		}
		if(B->node[B->script->replica[replicaN[0]].nodeSlot].messageWaitingIndicator){
			// There are only messages left to run the mobile server. Therefore the forcedatabase is already closed
			allow_write_record=false;
		}
		send_simulation_parameters(B->client,current_replica,B->script,B->node,replicaN[0]);
	}

	close(B->client->fd);

	if(allow_write_record && !newConnection && !node_just_reanimated && !unexpectedClient){
		// can't allow this if we're in the process of exiting and have closed the force_database
		// if node_just_reanimated then I don't want to save this stuff
		pthread_mutex_lock(&database_mutex);
		for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){
			if(save_sample_data_replicaN[nni]>=0){
				force_database->set_record_parameters(save_sample_data_replicaN[nni], save_sample_data_sequence_number[nni], save_sample_data_w[nni]);
				if(sample[nni].data!=NULL) force_database->add_all_forces_at_once(sample[nni].data);
				for(ai=0; ai<B->script->Nadditional_data; ai++){
					if(additional_data[nni][ai].data!=NULL){
						force_database->add_all_additionals_at_once(additional_data[nni][ai].data,ai);
						/*
						 * This section (with variables above), has proved to CN
						 * that additional_data[ai].data is correct
						 * int z;
						 * float *p;
						 * p=(float *)additional_data[ai].data;
						 * B->client->ptr+=sprintf(B->client->ptr,"CN INFO\n");
						 * for(z=0; z<41; z++){
						 * 	B->client->ptr+=sprintf(B->client->ptr,"%d\t%d\t%f\n",ai,z,p[z]);
						 * }
						 */
					}else{
						B->client->ptr+=sprintf(B->client->ptr,"Error: additional_data[%d][%d].data==NULL\n",nni,ai);
					}
				}
				force_database->write_record();
			}
		}
		pthread_mutex_unlock(&database_mutex);
		if(B->opt->verbose){
			B->client->ptr+=sprintf(B->client->ptr,"Sample data saved to database\n");
		}
	}else{
		if(B->opt->verbose){
			if(newConnection){
				B->client->ptr+=sprintf(B->client->ptr,"Sample data NOT saved to database -- it's a new connection so there is no data\n");
			}else if (node_just_reanimated){
				B->client->ptr+=sprintf(B->client->ptr,"Sample data NOT saved to database -- unexpected client (being resurected)\n");
			}else if (unexpectedClient){
				B->client->ptr+=sprintf(B->client->ptr,"Sample data NOT saved to database -- unexpected client (being released)\n");
			}else{
				B->client->ptr+=sprintf(B->client->ptr,"Sample data NOT saved to database -- message passing for Mobile server instead\n");
			}
		}
	}


	for(nni=0;nni<B->script->Nsamesystem_uncoupled;nni++){	
		printf("freeing restart: pointer before is: %p\n",current_replica[nni].restart.data);  //##DEBUG
		delete[] current_replica[nni].restart.data;
		printf("freeing restart: pointer after is: %p\n",current_replica[nni].restart.data);  //##DEBUG

		printf("freeing sample: pointer before is: %p\n",sample[nni].data);  //##DEBUG
		delete[] sample[nni].data;
		printf("freeing sample: pointer after is: %p\n",sample[nni].data);   //##DEBUG

		for(ai=0; ai<B->script->Nadditional_data; ai++){
			printf("freeing additional[%d][%d]: pointer before is: %p\n",nni,ai,additional_data[nni][ai].data);  //##DEBUG
			delete[] additional_data[nni][ai].data;
			printf("freeing additional[%d][%d]: pointer after is: %p\n",nni,ai,additional_data[nni][ai].data);   //##DEBUG
		}
	}
	free(energy);
	free(sample);
	for(nni=0; nni<B->script->Nsamesystem_uncoupled; nni++){
		free(additional_data[nni]);
	}
	free(additional_data);
	free(coordinate);
	free(current_replica);
	free(replicaN);
	free(old_replicaN);
	free(save_sample_data_replicaN);
	free(save_sample_data_sequence_number);
	free(save_sample_data_w);
	free(bin);

	struct timeval t;
	gettimeofday(&t,NULL);
	B->client->ptr+=sprintf(B->client->ptr,"-  - --- Client interaction ends; time elapsed is %ld ms -------------------------------- -  -\n",(t.tv_sec-B->client->time.tv_sec)*1000+(t.tv_usec-B->client->time.tv_usec)/1000);
	
	append_log_entry(B->client->fd, B->client->log);
	
	delete B->client;
	change_number_of_connected_clients(-1,B->var);
	// B was created in the calling function
	delete B;
	
	return(NULL);
}

// This runs as a thread and simply listens on the port for incoming client connections
// when a client connects, a thread is created for the new client and then this begins listening again on the port for further clients
void *wait_for_clients(struct client_bundle *B){
	struct sockaddr_in serv_addr;
	struct sockaddr_in client_addr;
	int server_sockfd;
	pthread_t server_handle;
	int pthread_failure=0;

	signal(SIGPIPE, SIG_IGN);

	if ((server_sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) error_quit("cannot open socket");
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	serv_addr.sin_port = htons(B->script->port);
	if(bind(server_sockfd, (sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) error_quit("cannot bind socket");
	if(listen(server_sockfd,100)<0) error_quit("cannot listen to socket");

	append_log_entry(-1,"Waiting for clients to connect...\n");
	while(B->var->simulation_status!=Finished){
		int len;
		int client_sockfd;
		unsigned int Nmoves=0;

		len = sizeof(client_addr);
		client_sockfd = accept(server_sockfd,(sockaddr *)&client_addr,(socklen_t *)&len);

		printf("simulation_status=%d\n",(int)B->var->simulation_status); //##DEBUG
		if(B->var->simulation_status!=Finished){
			printf("About to allocate memory\n"); //##DEBUG
			struct client_struct* client_data=new struct client_struct;

			printf("Allocated memory for client\n"); //##DEBUG

			client_data->fd=client_sockfd;
			gettimeofday(&client_data->time,NULL);
			client_data->ptr=client_data->log;

			printf("Added time of day to the clock\n"); //##DEBUG

			change_number_of_connected_clients(+1,B->var);

			printf("Number of connected clients increased\n"); //##DEBUG

			sprintf(client_data->ip,"%s",inet_ntoa(client_addr.sin_addr));
			client_data->ptr+=sprintf(client_data->ptr,"-  - --- Client has connected from IP address %s --- -  -\n",client_data->ip);

			printf("Client has connected\n"); //##DEBUG

			struct client_bundle* Blocal=new struct client_bundle;
			Blocal->client=client_data;
			Blocal->opt=B->opt;
			Blocal->var=B->var;
			Blocal->script=B->script;
			Blocal->node=B->node;
			//client_interaction() will delete Blocal
			pthread_failure=0;
			if(pthread_create(&server_handle,NULL, (void* (*)(void*))client_interaction, Blocal)!=0){
				error_warning("pthread_create failed in DR_server");
				pthread_failure=1;
			}else{
				if(pthread_detach(server_handle)!=0){
					error_warning("pthread_detach failed");
					pthread_failure=1;
				}
			}
			if(pthread_failure){
				client_data->ptr+=sprintf(client_data->ptr,"*** Client interaction from IP address %s aborted... no thread for client_interaction().\n", inet_ntoa(client_addr.sin_addr));
			}
		}
	}
	printf("server_sockfd is %d\n",server_sockfd);  //##DEBUG
	close(server_sockfd);

	return(NULL);
}

// prints out the details of the distributed replica simulation that is about to run
void print_simulation_details(const struct script_struct *script, const struct server_option_struct *opt){
	//Why the heck is this is DR_server and not read_input_script_file.h wonders CN
	char message[MESSAGE_GLOBALVAR_LENGTH];
	int i,j;

	append_log_entry(-1,                "*  * *********************************************************************** *  *\n");
	append_log_entry(-1,        "          D I S T R I B U T E D   R E P L I C A   S I M U L A T I O N\n");
	if(script->coordinate_type==Spatial)
		append_log_entry(-1,        "                           along a spatial coordinate\n");
	else if(script->coordinate_type==Temperature)
		append_log_entry(-1,        "                         along a temperature coordinate\n");
	else
		append_log_entry(-1,        "                     along an umbrella position coordinate\n");
	if(script->replica_move_type==MonteCarlo){
		sprintf(message,  "     using Monte Carlo moves with a step size of %0.3f linearized unit(s)\n",script->replica_step_fraction);
		append_log_entry(-1,message);
	}else if(script->replica_move_type==vRE){
		sprintf(message,   "        using vRE moves with a step size of %0.3f linearized unit(s)\n",script->replica_step_fraction);
		append_log_entry(-1,message);
	}else if(script->replica_move_type==BoltzmannJumping){
		sprintf(message,  "               using discrete Boltzmann jumps (%u possibilities)\n",script->Nreplicas);
		append_log_entry(-1,message);
	}else if(script->replica_move_type==Continuous)
		append_log_entry(-1,        "                           using a continuous space\n");
	else if(script->replica_move_type==NoMoves)
		append_log_entry(-1,        "                            with no move attempts\n");
	if(script->circular_replica_coordinate)
		append_log_entry(-1,        "          moves will be attempted directly between first and last replicas\n");

	append_log_entry(-1,                "*  * *********************************************************************** *  *\n");

	append_log_entry(-1,"Input parameters:\n");
	sprintf(message,"Title of this run: %s\n",opt->title);
	append_log_entry(-1,message);
	if(script->coordinate_type!=Temperature){
		//CN changed 052908l used to be if ==Spatial
		sprintf(message,"Temperature: %f\n",script->temperature);
		append_log_entry(-1,message);
	}
	sprintf(message,"TCP/IP communications server port: %hu\n",script->port);
	append_log_entry(-1,message);
	sprintf(message,"Name of log file: %s\n",logFile_globalVar);
	append_log_entry(-1,message);

	sprintf(message,"Distributed replica potential scalars: %f %f\n",script->replica_potential_scalar1, script->replica_potential_scalar2);
	append_log_entry(-1,message);
	if(script->cancellation_threshold>0){
		sprintf(message,"Energy cancellation will be done when sequence number %d is reached for all replicas\n", script->cancellation_threshold);
		append_log_entry(-1,message);
		sprintf(message," at that time, the following replica potential scalars will be used: %f %f\n",script->replica_potential_scalar1_after_threshold, script->replica_potential_scalar2_after_threshold);
		append_log_entry(-1,message);
	}
	if(script->replica_move_type==vRE){
		sprintf(message,"Virtual Replica Exchange (vRE) moves will not be attempted for the first %ld steps, during which a pool of virtual reverse moves will be accumulated.\n",script->vRE_initial_noMoves);
		append_log_entry(-1,message);
		sprintf(message,"Virtual Replica Exchange (vRE) moves will not be recorded for the first %ld steps, we assume because you think that the starting structures are not properly equilibrated.\n",script->vRE_initial_noSave);
		append_log_entry(-1,message);
		sprintf(message,"Virtual Replica Exchange (vRE) secondary list has a length of %ld values at each nominal position.\n",script->vRE_secvre_size);
		append_log_entry(-1,message);
	}
	sprintf(message,"A Node will be occupied for a maximum of %d seconds\n",script->node_time);
	append_log_entry(-1,message);
	sprintf(message,"A change of replica running on a Node can occur as often as %d seconds\n",script->replica_change_time);
	append_log_entry(-1,message);
	sprintf(message,"A snapshot of the state of all replicas will be saved every %d seconds\n",script->snapshot_save_interval);
	append_log_entry(-1,message);
	sprintf(message,"Maximum time allowed for one sequence number to finish is %d seconds\n",script->job_timeout);
	append_log_entry(-1,message);
	if(script->allow_requeue){
		sprintf(message,"\tNOTE: when a client does return to the server after this timeout, it may be assigned a new job\n");
		append_log_entry(-1,message);
		sprintf(message,"\t      This may be inefficient when a Node is truely down, but it is good in the case of a preemptive test queue that may preempt a job making it appear to have crashed when it really has not.\n");
		append_log_entry(-1,message);
	}else{
		append_log_entry(-1,"NOTE: if a node returns after the timeout period has passed then it will be released.\n");
	}
	sprintf(message,"Replicas between %d and %d will be left unsuspended until sampling is complete for all replicas\n",script->min_unsuspended_replica,script->max_unsuspended_replica);
	append_log_entry(-1,message);
	if(script->stopOnAverageTimeExceeded){
		sprintf(message,"  ** Actually, a testing feature has been enabled in your xx.script file so that the run will terminate when the average number of sampling steps has been completed, although not necessarily what you asked for for each nominal position.\n");
		append_log_entry(-1,message);
	}
	if(script->allotted_time_for_server==0){
		sprintf(message,"The server will be allowed to run indefinately. If you are submitting the server to a machine with a wallclock limit, then you would be better to predefine a period of time after which the server will be shut down via ALLOTTED_TIME_FOR_SERVER\n");
		append_log_entry(-1,message);
	}else{
		sprintf(message,"The server will save a snapshot and exit after %u seconds. This is a good idea if you are submitting the server to a machine with a wallclock limit\n", script->allotted_time_for_server);
		append_log_entry(-1,message);
	}
	if(script->defineStartPos){
		sprintf(message,"The nominal positions are defined in your script, however, replicas will start at positions as defined in %s\n",DEFINED_START_POS_FILE);
		append_log_entry(-1,message);
	}else{
		sprintf(message,"The nominal positions are defined in your script and one replica will start at each nominal position\n");
		append_log_entry(-1,message);
	}
	if(script->cycleClients>-0.001){
		sprintf(message,"When a complete set of replicas is running and a new client connects, it may kick out an older client if that older client has already finished %f proportion of its allocated runtime.\n",script->cycleClients);
		append_log_entry(-1,message);
	}else{
		sprintf(message,"When a complete set of replicas is running and a new client connects, it will be dropped.\n");
		append_log_entry(-1,message);
	}

	sprintf(message,"%u Replicas will run with the following parameters:\n",script->Nreplicas);
	append_log_entry(-1,message);
	append_log_entry(-1,"-----------------------------------------------------------------------------------------\n");

	message[0]=0;
	sprintf(message+strlen(message),"       coord");
	if(!isnan(script->replica[0].w2_nominal)) sprintf(message+strlen(message),"      coord2");
	if(!isnan(script->replica[0].force)) sprintf(message+strlen(message),"       force");
	sprintf(message+strlen(message),"       moves");
	sprintf(message+strlen(message),"       steps");
	sprintf(message+strlen(message),"  cancel_energy");
	sprintf(message+strlen(message),"  nominal_coord");
	sprintf(message+strlen(message),"  starting_coord");
	sprintf(message+strlen(message),"\n");
	append_log_entry(-1,message);
	for(i=0;i<script->Nreplicas;i++){
		message[0]=0;
		sprintf(message+strlen(message),"%12.3f",script->replica[i].w_nominal);
		if(!isnan(script->replica[i].w2_nominal)) sprintf(message+strlen(message),"%12.3f",script->replica[i].w2_nominal);
		if(!isnan(script->replica[i].force)) sprintf(message+strlen(message),"%12.3f",script->replica[i].force);
		sprintf(message+strlen(message),"%12u",script->replica[i].sampling_runs);
		sprintf(message+strlen(message),"%12u",script->replica[i].sampling_steps);
		sprintf(message+strlen(message),"%15.3f",script->replica[i].cancellation_energy);
		sprintf(message+strlen(message),"%15.3f",script->replica[i].w_start);
		sprintf(message+strlen(message),"%16.3f",script->replica[i].w);
		sprintf(message+strlen(message),"\n");
		append_log_entry(-1,message);
	}
	append_log_entry(-1,"-----------------------------------------------------------------------------------------\n");
	if(script->Nsamesystem_uncoupled!=1){
		sprintf(message,"%u non-interacting replicas will be simulated as part of the same system.\n",script->Nsamesystem_uncoupled);
		for(i=0;i<ldiv(script->Nreplicas,script->Nsamesystem_uncoupled).quot;i++){
			sprintf(message+strlen(message),"\tSystem %d\t%12.3f",i,script->replica[i*script->Nsamesystem_uncoupled].w_nominal);
			for(j=1;j<script->Nsamesystem_uncoupled;j++){
				sprintf(message+strlen(message),"\t&%12.3f",script->replica[i*script->Nsamesystem_uncoupled+j].w_nominal);
			}
			sprintf(message+strlen(message),"\n");
		}
		append_log_entry(-1,message);
	}
	if(script->submit_jobs){
		append_log_entry(-1,"Queue shells will automatically be submitted using drsub\n");
	}else{
		append_log_entry(-1,"Queue shells will NOT be submitted\n");
	}
}

void showUserResponsibility(const char *c){
	fprintf(stderr,"WARNING: the following options are no longer supported:\n");
	fprintf(stderr,"         a) 2 ligands\n");
	fprintf(stderr,"In order to use one of these options, please use an earlier version of %s\n",c);
	fprintf(stderr,"The current program Maintainer, %s, has either found or introduced ", currentProgrammerName);
	fprintf(stderr,"problems with these methods and has chosen to leave them unsupported for lack of use.\n");
	fprintf(stderr,"If you would like to utilize one of these options, please contact %s at ",currentProgrammerName);
	fprintf(stderr,"%s as s/he may be interested in upgrading the code for your use.\n",currentProgrammerEmail);
	fprintf(stderr,"Nevertheless, you may continue at your own peril in spite of this warning by adding ");
	fprintf(stderr,"'-r %s' to the command line of %s\n",USER_RESPONSIBILITY_STRING,c);
}

void showUsage(const char *c){
	fprintf(stderr,"Usage: %s tt.script [-stdv]\n",c);
	fprintf(stderr,"       -s [string] to restart from a snapshot (e.g. tt.283429.snapshot)\n");
	fprintf(stderr,"       -t [integer] time server node started (for the mobile server)\n");
	fprintf(stderr,"       -d [string] directory in which to put the log file (e.g. /dev/shm)\n");
	fprintf(stderr,"       -v [integer] non-zero to have a more verbose log file\n");
	fflush(stderr);
}

int parseCommandLine(int argc, char * const argv[], struct server_option_struct *opt){
	int i;
	int gots=0;
	int gotr=0;
	int gott=0;
	int gotd=0;
	int gotv=0;

	if( (argc<2) ){
		fprintf(stderr,"Error: the script filename was not provided\n");
		return 1;
	}

	if(strlen(argv[1])<2){
		fprintf(stderr,"Error: script filename must have at least 2 characters (the first 2 are used for the title)\n");
		return 1;
	}
	opt->title[0]=argv[1][0];
	opt->title[1]=argv[1][1];
	opt->title[2]=0;

	if(isupper(opt->title[0])||isupper(opt->title[1])){
		//DR_client_comm pushed these vals to lowercase and they will not match unless the user uses lowercase --- so force this
		fprintf(stderr,"Error: The first 2 characters of the script filename must be lowercase\n");
		return 1;
	}

	if(argc==3){
		fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[2]);
		return 1;
	} 
	for(i=3; i<argc; i+=2){
		if(argv[i-1][0]!='-'){
			// Should only be dashed arguments from this point on.
			fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[i-1]);
			return 1;
		}
		if(argv[i-1][1]=='s'){
			if(gots){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->loadSnapshot=true;
			sscanf(argv[i],"%s",opt->snapshotName);
			gots=1;
		}else if(argv[i-1][1]=='r'){
			if(gotr){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			sscanf(argv[i],"%s",opt->userResponsibilityString);
			if(strcmp(USER_RESPONSIBILITY_STRING,opt->userResponsibilityString)!=0){
				fprintf(stderr,"Error: that is not the correct string.\n");
				showUserResponsibility(argv[0]);
				return 1;
			}
			opt->userTakesResponsibility=true;
			gotr=1;
		}else if(argv[i-1][1]=='t'){
			if(gott){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			sscanf(argv[i],"%d",&(opt->mobile_timeClientStarted));
			gott=1;
                }else if(argv[i-1][1]=='d'){
                        if(gotd){
                                fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
                                return 1;
                        }
                        sscanf(argv[i],"%s",opt->logdir);
			if(opt->logdir[strlen(opt->logdir)-1]!='/'){
				strcat(opt->logdir,"/");
			}
                        gotd=1;
                }else if(argv[i-1][1]=='v'){
                        if(gotv){
                                fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
                                return 1;
                        }
                        sscanf(argv[i],"%d",&(opt->verbose));
                        gotv=1;
		}else{
			fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[i-1]);
			return 1;
		}
	}
	return 0;
}

int defineStartingPositions(struct script_struct *script){
	FILE *ssf;
	char linein[1001];
	int *ssray;
	int ssi;
	char message[MESSAGE_GLOBALVAR_LENGTH];

	ssray=(int *)malloc(script->Nreplicas*sizeof(int));
	if(ssray==NULL){
		fprintf(stderr,"Error: unable to allocate memory for ssray\n");
		return 1;
	}
	ssf=fopen(DEFINED_START_POS_FILE,"r");
	if(ssf==NULL){
		sprintf(message,"Error: unable to open file %s in your working directory. You spedified DEFINE_STARTING_POSITIONS in the script file, so this file must exist.\n", DEFINED_START_POS_FILE);
		append_log_entry(-1,message);
		return 1;
	}
	ssi=0;
	while(fgets(linein,1000,ssf)!=NULL){
		if(ssi>=script->Nreplicas||sscanf(linein,"%d",&(ssray[ssi]))!=1){
			if(ssi>=script->Nreplicas){
				sprintf(message,"Error: Too many values in %s\n",DEFINED_START_POS_FILE);
				append_log_entry(-1,message);
			}else{
				sprintf(message,"Error: could not find value in %s -- could there be an empty or non-numeric line?\n",DEFINED_START_POS_FILE);
				append_log_entry(-1,message);
			}
			fclose(ssf);
			free(ssray);
			return 1;
		}
		ssi++;
	}
	fclose(ssf);
	if(ssi!=script->Nreplicas){
		sprintf(message,"Error: Not enough values in %s\n",DEFINED_START_POS_FILE);
		 append_log_entry(-1,message);
	}
	for(ssi=0;ssi<script->Nreplicas;ssi++){
		script->replica[ssi].w=script->replica[ssray[ssi]].w_nominal;
		//sprintf(message,"RE-DEFINE_STARTING_POSITIONS: assigining replica %d a starting value of %f\n",ssi, script->replica[ssi].w);
		//append_log_entry(-1,message);
	}
	free(ssray);
	return 0;
}

int findClientForServer(struct script_struct *script, struct server_variable_struct *var, const struct server_option_struct *opt, node_struct *node, int server_start_time){
	int recentnode;
	int recentrep;
	int timegain;
	char *snapshotname;
	char message[MESSAGE_GLOBALVAR_LENGTH];

	pthread_mutex_lock(&replica_mutex);
	append_log_entry(-1,"MOBILITY: FIND A NEW CLIENT FOR THE SERVER\n");
	recentrep=-1;
	for(int i=0;i<script->Nreplicas;i+=script->Nsamesystem_uncoupled){
		if(script->replica[i].nodeSlot>=0 && node[script->replica[i].nodeSlot].active){
			if (recentrep==-1 || node[script->replica[i].nodeSlot].start_time>node[script->replica[recentrep].nodeSlot].start_time){
				recentrep=i;
			}
		}
	}
	if(recentrep==-1){
		//no active nodes found at all
		pthread_mutex_unlock(&replica_mutex);
		return 1;
	}
	timegain=node[script->replica[recentrep].nodeSlot].start_time - server_start_time;
	if(timegain < script->mobility_requiredTimeGain){
		//a move would not be worth it
		pthread_mutex_unlock(&replica_mutex);
		return 1;
	}
	recentnode=script->replica[recentrep].nodeSlot;
	sprintf(message,"MOBILITY: The server will move to client node[%d] with ip %s that is currently running replica %d\n",recentnode, node[recentnode].ip, recentrep);
	append_log_entry(-1,message);

	// save a snapshot and call the forcedatabase destructor -- this will close it so that the new server can safely open it
	snapshotname=save_snapshot(script,var,opt);
	delete force_database;

	// send one node a message to become a new server and also sent it the name of the snapshot that it should use
	node[recentnode].messageWaitingIndicator=true;
	sprintf(node[recentnode].serverMessageForClient,"BECOME_NEW_SERVER %s", snapshotname);

	// send a message to all other ACTIVE nodes to start contacting the new server, and let them know at what ip they can do that
	for(int i=0;i<script->Nreplicas;i+=script->Nsamesystem_uncoupled){
		if(i!=recentnode && node[i].active){
			node[i].messageWaitingIndicator=true;
			sprintf(node[i].serverMessageForClient,"HOLD_AND_CONTACT %s", node[recentnode].ip);
		}
	}

	pthread_mutex_unlock(&replica_mutex);
	//the serverMessageForClient will be sent to the client during the send_simulation_parameters() called from client_interaction(),
	//  and then the messageWaitingIndicator will be turned off.

	//wait for all message indicators to be turned off. Must check for active nodes only as a node may be released by a check_for_crash call

	for(int i=0;i<script->Nreplicas;i+=script->Nsamesystem_uncoupled){
		while(node[i].active && node[i].messageWaitingIndicator){
			sleep(1);
			check_for_crash(i,script,var,node);
			sprintf(message,"MOBILITY (WAITING): mobility trying to have node[%d] get the message.\n",i);
			append_log_entry(-1,message);
			displayNodes(script,node);
		}
		sprintf(message,"MOBILITY (one stopped WAITING): node %d has checked in and got the message.\n",i);
		append_log_entry(-1,message);
	}
	append_log_entry(-1,"MOBILITY (all stopped WAITING): all clients have reported in regarding the mobility\n");
	return 0;
}


int main(int argc, char *argv[]){
	int i,j;
	unsigned int N;
	int last_crash_check;
	int last_unconditional_start_queue_shell;
	int last_start_queue_shell_check;
	int last_snapshot_save;
	int last_disk_almost_full_check;
	int last_finish_on_average_check;
	int last_node_display;
	int last_mobility_check;
	int start_time;
	int this_server_start_time;
	pthread_t server_handle;
	int skipFinalSnapshot=0;
	int crash_checking_interval=0;
	char message[MESSAGE_GLOBALVAR_LENGTH];

	struct script_struct script;
	struct server_option_struct opt=DEFAULT_SERVER_OPTION_STRUCT;
	struct server_variable_struct var=DEFAULT_SERVER_VARIABLE_STRUCT;
	struct database_struct db=EMPTY_DATABASE_STRUCT;
	struct node_struct *node;

	int check;
	
	check=parseCommandLine(argc,argv,&opt);
	if(check!=0){
		fprintf(stderr,"Error: parseCommandLine() returned non-zero\n");
		showUsage(argv[0]);
		exit(check);
	}

	if(getcwd(var.working_directory,198)==NULL){
		fprintf(stderr,"Error: getcwd() returned non-zero\n");
		exit(check);
	}

	read_input_script_file_class *input_script = new read_input_script_file_class;
	input_script->read_input_script_file(argv[1], &script);
	delete input_script;

	printf("Input script read in\n"); //##DEBUG
	sprintf(logFile_globalVar,"%s%s.log",opt.logdir,opt.title);

	if(opt.loadSnapshot && script.defineStartPos){
		error_quit("You have both loaded a snapshot (-s option to DR_server) and defined the starting positions (in the script file). If you are restarting a run, you probably want to remove the DEFINE_STARTING_POSITIONS flag from the script file. Exiting to try to help you.\n");
	}

	char forcedbnam[64];
	sprintf(forcedbnam,"%s.forcedatabase",opt.title);

	FILE *ffd = fopen(forcedbnam,"r");
	if(ffd) {
		// exists
		fclose(ffd);
		if(!opt.loadSnapshot){
			//ensure that the forcedatabase file doesn't already exist if starting a new run
			//#define DATABASE_FILENAME "%s.forcedatabase"
			error_quit("You are not loading a snaphot, therefore this appears to be a new run, so why does a .forcedatabase file exist? Exiting to avoid problems.\n");
		}
	}

	var.beta=(double)1.0/(script.temperature*BOLTZMANN_CONSTANT);
	var.min_running_replica=0;

	var.max_running_replica=script.Nreplicas-1;

	printf("Running=%d - Finished=%d - DiskAlmostFull=%d\n",Running,Finished,DiskAlmostFull); //##DEBUG

	int s=time(NULL);
	srand48(s);
	sprintf(message,"Seeding random number generator with %d\n",s);
	append_log_entry(-1,message);

	pthread_mutex_init(&replica_mutex,NULL);
	pthread_mutex_init(&log_mutex,NULL);
	pthread_mutex_init(&queue_mutex,NULL);
	pthread_mutex_init(&database_mutex,NULL);
	pthread_mutex_init(&vre_mutex,NULL);

	if(script.replica_move_type==vRE){
		set_secvre_size(script.vRE_secvre_size); //must be called before allocateVRE()
		if(allocateVRE(script.Nreplicas,-1)!=0){
			error_quit("unable to allocate memory for vRE structure\n");
		}
		//CN intentionally loads in these files prior to loading a snapshot
		for(int zz=0;zz<script.Nreplicas;zz++){
			if(loadFileIntoVREforStartup(zz, script.replica[zz].vREfile)!=0){
				error_quit("unable to load vRE initialization data\n");
			}
		}
	}

	//Messy routine to allow loading new cancellation values when using a snapshot (continued later by loading temp back)
	float *temp;
	temp=(float *)malloc(script.Nreplicas*sizeof(float));
	if(temp==NULL)exit(-999);
	int tempi;
	for(tempi=0;tempi<script.Nreplicas; tempi++){
		temp[tempi]=script.replica[tempi].cancellation_energy;
	}

	if(opt.loadSnapshot){
		load_snapshot(opt.snapshotName,&script,&var);
		printf("Snapshot loaded\n"); //##DEBUG
	}

	for(tempi=0;tempi<script.Nreplicas; tempi++){
		script.replica[tempi].cancellation_energy=temp[tempi];
	}

	struct client_struct *client=new client_struct;
	client->ptr=client->log;
	client->ptr+=sprintf(client->ptr,"Doing an initial check to see which replicas should be suspended...\n");
	check_termination_conditions(client,&script,&var,&opt);
	append_log_entry(-1,client->log);
	delete client;

	force_database=new force_database_class(opt.title,0);
	db.Nligands=script.Nligands;
	if(script.need_sample_data){
		db.Nforces=script.Nsamples_per_run;
	}else{
		db.Nforces=0;
	}
	//C. Neale wonders why Nenergies is always zero -- would be useful in a temperature run
	db.Nenergies=0;
	db.Nadditional_data=script.Nadditional_data;

	if(!force_database->set_header_information(&db)) error_quit("the force database header has inapropriate parameters");

	sprintf(message,"The database file has been successfully opened or created: number of records: %u; number of ligands: %hhu; number of samples per run: %u; number of energy data per run: %u; number of additional data types: %u\n",force_database->get_number_of_records(),db.Nligands,db.Nforces,db.Nenergies,db.Nadditional_data);
	append_log_entry(-1,message);

	if(script.cancellation_threshold>0){
		var.energy_cancellation_status=Pending;
		conditionally_activate_energy_cancellation(&script,&var);
		if(var.energy_cancellation_status==Active){
			append_log_entry(-1,"The energy cancellation feature has been activated; please stand by for a summary\n");
			if(var.simulation_status==Finished){
				print_energy_cancellation_summary(&script,&var);
			}
		}
	}

	// Need to load the saved cancellations again
	for(tempi=0;tempi<script.Nreplicas; tempi++){
		script.replica[tempi].cancellation_energy=temp[tempi];
	}

	if(script.loadedCancel){
		print_energy_cancellation_summary(&script,&var);
		script.loadedCancel=false;
	}

	if(script.defineStartPos){
		if(defineStartingPositions(&script)!=0){
			//CN notes that this will cause problems for the Mobile server when the user selects this option -- b/c error_quit in routine
			error_quit("defineStartingPositions failed in DR_server Main.\n");
		}
	}

	print_simulation_details(&script,&opt);

	//struct *node_struct node;
	if((node=(struct node_struct *)malloc(div(script.Nreplicas,script.Nsamesystem_uncoupled).quot*sizeof(struct node_struct)))==NULL){
		error_quit("Unable to allocate memory for struct node_struct in DR_server Main. This is a top level error, try restarting your server.\n");
	}
	for(tempi=0;tempi<div(script.Nreplicas,script.Nsamesystem_uncoupled).quot; tempi++){
		node[tempi].active=false;
		node[tempi].ip[0]='\0';
		node[tempi].start_time=0;
		node[tempi].awaitingDump=false;
		node[tempi].serverMessageForClient[0]='\0';
		node[tempi].messageWaitingIndicator=false;
	}

	last_crash_check=
	last_unconditional_start_queue_shell=
	last_start_queue_shell_check=
	last_snapshot_save=
	last_disk_almost_full_check=
	last_finish_on_average_check=
	last_node_display=
	last_mobility_check=
	this_server_start_time= time(NULL);

	if(opt.mobile_timeClientStarted>0){
		//this server has already run for some period of time
		start_time=opt.mobile_timeClientStarted;
	}else{
		start_time=this_server_start_time;
	}

	struct client_bundle* B=new struct client_bundle;
	B->client=(struct client_struct *)NULL;
	B->opt=&opt;
	B->var=&var;
	B->script=&script;
	B->node=node;
	if(pthread_create(&server_handle,NULL,(void* (*)(void*))wait_for_clients,B)!=0){
		error_quit("pthread_create failed in DR_server Main. This is a top level error, simply try restarting your server.");
	}
	sleep(10); //don't want the clients to start before we are ready for them

	if(script.submit_jobs){
		for(i=0;i<ldiv(script.Nreplicas,script.Nsamesystem_uncoupled).quot;i++){
			for(j=0;j<script.Nsamesystem_uncoupled; j++){
				//can't base this test on =='N' because replicas might start before I get through submitting all
				if(script.replica[i*script.Nsamesystem_uncoupled+j].status!='S'){
					start_queue_shell(true,&script,&var);
					break;
				}
			}
		}
	}

	crash_checking_interval=div(script.job_timeout,2).quot;

	while(var.simulation_status!=Finished && var.simulation_status!=AllottedTimeOver){
		sleep(1);
		
		if(time(NULL)-last_crash_check>=crash_checking_interval){
			last_crash_check=time(NULL);
			print_number_of_connected_clients(&var);
			append_log_entry(-1,"Checking if any replicas crashed...\n");
			check_for_crash(-1,&script,&var,node);
		}
		if(script.submit_jobs && 
		  ((time(NULL)-last_start_queue_shell_check>=script.node_time/(ldiv(script.Nreplicas,script.Nsamesystem_uncoupled).quot)) || 
		  (var.Ncrashed_jobs>0))){ 
			last_start_queue_shell_check=time(NULL);
			start_queue_shell(true,&script,&var); // conditionally start a shell
		}
		if(script.submit_jobs && (time(NULL)-last_unconditional_start_queue_shell>=QUEUE_INTERVAL) ){
			last_unconditional_start_queue_shell=time(NULL);
			start_queue_shell(false,&script,&var); // unconditionally start a shell
		}
		if(time(NULL)-last_disk_almost_full_check>=DISK_ALMOST_FULL_CHECK_SECONDS){
			last_disk_almost_full_check=time(NULL);
			check_if_disk_almost_full(&var);
		}
		if(script.stopOnAverageTimeExceeded && time(NULL)-last_finish_on_average_check>=FINISH_ON_AVERAGE_CHECK_SECONDS){
			last_finish_on_average_check=time(NULL);
			check_if_finish_on_average(&script,&var);
		}
		if(time(NULL)-last_snapshot_save>=script.snapshot_save_interval){
			last_snapshot_save=time(NULL);
			var.save_snapshot_now=true;
		}
		if(var.save_snapshot_now){
			pthread_mutex_lock(&replica_mutex);
			save_snapshot(&script,&var,&opt);
			pthread_mutex_unlock(&replica_mutex);
			var.save_snapshot_now=false;
		}
		if(var.energy_cancellation_status==Active){
			print_energy_cancellation_summary(&script,&var);
			var.energy_cancellation_status=Active_and_Printed;
		}
		if(script.allotted_time_for_server>0 && (time(NULL)-start_time>script.allotted_time_for_server)){
			var.simulation_status=AllottedTimeOver;
			sprintf(message,"Allotted server simulation time of %u seconds has been consumed. The run will now save a snapshot and exit\n",script.allotted_time_for_server);
			append_log_entry(-1,message);
		}
		if(script.mobility_time>0 && (script.allotted_time_for_server - (time(NULL)-start_time) < script.mobility_time) && (time(NULL)-this_server_start_time > script.job_timeout*2) && (time(NULL)-last_mobility_check > MOBILITY_CHECK_SECONDS)){
			last_mobility_check=time(NULL);
			// the condition based on this_server_start_time is to allow new clients to connect before mobility is possible
			if((findClientForServer(&script,&var,&opt,node,start_time))==0){
				var.simulation_status=Finished;
				skipFinalSnapshot=1;
				break;
			}
		}

		if(script.submit_jobs==true && var.nfailedsubinarow>MAX_FAILURES_FOR_SUBMISSION){
			script.submit_jobs=false;
			sprintf(message,"ERROR error Error: failed to submit a new client %d times in a row. Turning off submission.\n",MAX_FAILURES_FOR_SUBMISSION);
			append_log_entry(-1,message);
		}

		if(time(NULL)-last_node_display>=NODE_DISPLAY_SECONDS){
			last_node_display=time(NULL);
			pthread_mutex_lock(&replica_mutex);
			displayNodes(&script,node);
			pthread_mutex_unlock(&replica_mutex);
		}

	}

//	pthread_cancel(server_handle);

	if(!skipFinalSnapshot){
		pthread_mutex_lock(&replica_mutex);
		save_snapshot(&script,&var,&opt);
		pthread_mutex_unlock(&replica_mutex);
		delete force_database;
	}

	delete B;
	free_all_replicas(&script);

	append_log_entry(-1,"======================- Session End -======================\n");

	pthread_mutex_destroy(&replica_mutex);
	pthread_mutex_destroy(&log_mutex);
	pthread_mutex_destroy(&queue_mutex);
	pthread_mutex_destroy(&database_mutex);
	pthread_mutex_destroy(&vre_mutex);

	return(0);
}

