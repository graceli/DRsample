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

 DR_tester

 By Tomas Rodinger
 tom.rodinger@utoronto.ca

 Contributions by
 Chris Neale (distributed replica in circular space)

 v1.1 (Tom Rodinger June 22 2007)
	 -code to handle the "Umbrella" simulation method has been added to Monte Carlo and Boltzmann jumping functions

 v1.1.1 (Chris Neale July 24 2007)
	 - Modified to use read_input_script_file.h because that is now the most up to date version
		 The changes in that file do not affect this program and were intended for analyse_force_database.cpp

 C. Neale Oct 15 2007. 
	 - Tester Umbrella code instituted in v1.1 is incorrect. In order to fix this, change particle_x[] definition.
		 Was: particle_x[0]=xref_of_particle0=x_of_particle0
					particle_x[1]=x_of_particle1
		 New: particle_x[0]=xref_of_particle0
					particle_x[1]=x_of_particle1
					particle_x[2]=x_of_particle0
		 For Spatial particle_x[0]==particle_x[2]

 C. Neale May 30 2008.
	 - Modified the code to allow script->Nsamesystem_uncoupled to work

 Chris Neale June 6 2008
	 - Added free() statements for all mallocs

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <arpa/inet.h>
#include <fstream>
#include <string>
#include <zlib.h>
#include <netinet/in.h>
#include <errno.h>
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

#include "DR_protocol.h"
#include "read_input_script_file.h"

unsigned int protocol_version=PROTOCOL_VERSION;
// #define RST_FILE_SIZE (2*sizeof(float))
// But now we must send 3 pieces of information in the restart file since 
// x[2] is the actual x coord for Umbrella simulation (also Temperature)
#define RST_FILE_SIZE (3*sizeof(float))
// This will be further modified by script->Nsamesystem_uncoupled

float T;
float kT;
float B;
#define INITIAL_T 298.0

using namespace std;

// Start of block of variables locked by the replica_mutex
pthread_mutex_t replica_mutex;
unsigned int Nfinished_replicas=0;
// End of block of variables locked by the replica_mutex

struct sockaddr_in server_address; // address information of server we are connecting to

struct tester_option_struct{
	int includeNoise;
	int sleepTime;
	char title[3];
	int numtosubmit;
	char exactInputFile[500]; //There is no error checking for overflow
};
#define DEFAULT_TESTER_OPTION_STRUCT {1,100000,"  ",-1,""}
static int verbose_globalVar; //not part of struct since it is a debugging feature only

struct client_bundle{
	struct script_struct *script;
	struct tester_option_struct *opt;
};

//------------------------------------------------------------------------------------------------------------------

#define sqr(x) ((x)*(x))

void error_quit(const char *message){
	fprintf(stderr,"%s\n",message); 
	exit(1);
}

#define INTERPARTICLE_FORCE_CONSTANT 0.5
#define COEF6 0.0025694
#define COEF5 0.0562981
#define COEF4 0.4240118
#define COEF3 1.2804560
#define COEF2 1.6494678
#define COEF1 2.4635373
#define COEF0 4.956993

//One could make a more simple test by: /E = x*x; F = -2*x;

double compute_energy(const double particle_x[], double fc, const struct script_struct *script, const struct tester_option_struct *opt){
	double E;
	double x;

/*
 *        Was: particle_x[0]=xref_of_particle0=x_of_particle0
 *             particle_x[1]=x_of_particle1
 *        New: particle_x[0]=xref_of_particle0 (or temperature)
 *             particle_x[1]=x_of_particle1
 *             particle_x[2]=x_of_particle0
 *        For Spatial particle_x[0]==particle_x[2]
 */       

	//Interaction of primary particle with free energy surface
	x=particle_x[2];
	E=COEF6*(x*x*x*x*x*x) - COEF5*(x*x*x*x*x) + COEF4*(x*x*x*x) - COEF3*(x*x*x) + COEF2*(x*x) - COEF1*x + COEF0;

	//Interaction of primary particle with its own umbrella
	if(script->coordinate_type==Umbrella){
		// can not allow this to occur for temperature since particle_x[0] is the temperature
		// For spatial, there is no additional energy term here
		E+=0.5*fc*sqr(particle_x[2]-particle_x[0]);
	}

	//Interaction of secondary particle with primary particle (noise)
	if(opt->includeNoise){
		E += (double)0.5*INTERPARTICLE_FORCE_CONSTANT*sqr(particle_x[2]-particle_x[1]);
	}

	return(E);
}

double compute_force(const double particle_x[],double fc, const struct script_struct *script, const struct tester_option_struct *opt){
	double F;
	double x;

	if(script->coordinate_type==Spatial || script->coordinate_type==Temperature){
		//We send the entire system force for Spatial and Temperature simulations
		//Interaction of primary particle with free energy surface	
		x=particle_x[2];
		F = -( COEF6*6*(x*x*x*x*x) - COEF5*5*(x*x*x*x) + COEF4*4*(x*x*x) - COEF3*3*(x*x) + COEF2*2*(x) - COEF1*1 + COEF0*0 );

		//Interaction of secondary particle with primary particle (noise)
		//Not sure why T. Rodinger had F+= and not F-= below 
		//CN changed it to F-= Oct 18 2007 based on assumption that the two particles are attractive
		if(opt->includeNoise){
			//F+= (double)INTERPARTICLE_FORCE_CONSTANT*(particle_x[2]-particle_x[1]);
			F-= (double)INTERPARTICLE_FORCE_CONSTANT*(particle_x[2]-particle_x[1]);
		}
	}else if(script->coordinate_type==Umbrella){
		//We send only the negative umbrella force for Umbrella simulations
		//Interaction of primary particle with its own umbrella
		//We send the negative umbrella force (averages to the system force)
		//Since umbrella force (F=-dE/dx) is neg sign, we send the positive sign.
		F=fc*(particle_x[2]-particle_x[0]);
	}

	return(F);
}

void run_simulation(unsigned int n_steps, double x[], float *sample_data, float *additional_data, double fc, const struct script_struct *script, const struct tester_option_struct *opt){
	unsigned int i,j;
	double old_E, new_E;
	double old_x1,old_x0;
	double probability;

/*
 *        Was: particle_x[0]=xref_of_particle0=x_of_particle0
 *             particle_x[1]=x_of_particle1
 *        New: particle_x[0]=xref_of_particle0
 *             particle_x[1]=x_of_particle1
 *             particle_x[2]=x_of_particle0
 * For Spatial particle_x[0]==particle_x[2]
 */

	//printf("START %f %f %f\n",x[0],x[1],x[2]);

	j=0;
	for(i=0;i<n_steps;i++){
		fprintf(stderr,"CN: i=%d\tx=%f\t%f\t%f\n",i,x[0],x[1],x[2]); //##DEBUG
		old_E=compute_energy(x,fc,script,opt);
		//move secondary particle
		old_x1=x[1];
		old_x0=x[2];
		x[1]+=(drand48()-0.5)*0.1;
		if(script->coordinate_type==Umbrella || script->coordinate_type==Temperature){
			//move primary particle
			x[2]+=(drand48()-0.5)*0.1;
		}
		new_E=compute_energy(x,fc,script,opt);
		probability=exp(-(new_E-old_E)*B);
		if(probability<drand48()){
			//change not accepted
			x[1]=old_x1;
			x[2]=old_x0;
		}
		if(sample_data!=NULL){
			sample_data[j++]=compute_force(x,fc,script,opt);
			//fprintf(stderr,"sample force of %f\n",compute_force(x,fc,script,opt));
			//CN Note that fc below is clearly the incorrect one. However, I am not supporting
			//multiple ligands, so leave it be.
			if(script->Nligands==2) sample_data[j++]=compute_force(x,fc,script,opt);
			if(additional_data!=NULL){
				additional_data[j-script->Nligands]=x[2];
			}
		}else{
			fprintf(stderr,"sample_data was NULL\n");
		}
	}

	//printf("END %f %f %f\n",x[0],x[1],x[2]);
}

void read4K(int sockfd, void *buff, int nbytes){
	int count,fourK=4096,offset=0;
	int iread,total=nbytes;

	while(nbytes>0){
		count=nbytes<fourK?nbytes:fourK;
		//fprintf(stderr,"attempting to read %d bytes\n",count); //##DEBUG
		iread=read(sockfd,(char *)buff+offset,count);
		//fprintf(stderr,"bytes actually read: %d\n",iread); //##DEBUG
		if(iread==0){
			fprintf(stderr,"read4K: read return 0, expecting %d more bytes\n",nbytes);
			exit(1);
		}else if(iread<0){
			fprintf(stderr,"read4K: read return -1\n");
			perror("read");
			exit(1);
		}
		nbytes-=iread;
		offset+=iread;
	}
}

void send_replica_ID(int sockfd, struct ID_struct *ID){
	unsigned int send_size;
	unsigned char buff[KEY_SIZE+COMMAND_SIZE+sizeof(*ID)];
	
	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=ReplicaID;
	memcpy(buff+KEY_SIZE+COMMAND_SIZE,ID,sizeof(*ID));
	send_size=KEY_SIZE+COMMAND_SIZE+sizeof(*ID);
	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the replica ID\n");
		exit(1);
	}
	if(verbose_globalVar)fprintf(stderr,"%u bytes sent\n",send_size);
}

void send_NextNonInteracting(int sockfd){
	unsigned int send_size;
	unsigned char buff[KEY_SIZE+COMMAND_SIZE];

	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=NextNonInteracting;
	send_size=KEY_SIZE+COMMAND_SIZE;
	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the identifier for NextNonInteracting\n");
		exit(1);
	}
	if(verbose_globalVar)fprintf(stderr,"NextNonInteracting command sent\n",send_size);
}


void send_energy_file(int sockfd, double particle_x[], float new_coord, const struct script_struct *script, const struct tester_option_struct *opt){
	unsigned int send_size,data_size;
	double saved_x0;
	float *energy;
	int i;
	
	if(script->replica_move_type==NoMoves) return;
	else if(script->coordinate_type==Temperature || script->coordinate_type==Umbrella ) data_size=sizeof(float);
	else if(script->replica_move_type==MonteCarlo) data_size=2*sizeof(float);
	else data_size=script->Nreplicas*sizeof(float);

	send_size=KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int)+data_size;

	unsigned char buff[send_size];
	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=TakeMoveEnergyData;
	
	*(unsigned int *)(buff+COMMAND_LOCATION+COMMAND_SIZE)=data_size;
	energy=(float *)(buff+COMMAND_LOCATION+COMMAND_SIZE+sizeof(unsigned int));

	saved_x0=particle_x[2];
	switch(script->replica_move_type)
	{
	case MonteCarlo:
	case vRE:
		if(script->coordinate_type==Umbrella){
			//Send the current position
			*energy=particle_x[2];
		}else if(script->coordinate_type==Spatial){
			*energy=new_coord;
			*(energy+1)=-compute_energy(particle_x,0.0,script,opt); //fc=0 for spatial
			particle_x[2]=(double)new_coord;
			*(energy+1)+=compute_energy(particle_x,0.0,script,opt); //fc=0 for spatial
		}else if (script->coordinate_type==Temperature){
			*energy=compute_energy(particle_x,0.0,script,opt); //fc=0 for temperature
		}
		break;
	case BoltzmannJumping:
	case Continuous:
		if(script->coordinate_type==Umbrella){
			//Send the current position
			*energy=particle_x[2];
		}else if(script->coordinate_type==Spatial){
			//Continuous Spatial can not get here. If changes are made to allow that then this will need to be updated
			for(i=0;i<script->Nreplicas;i++){
				particle_x[0]=script->replica[i].w_nominal;
				*(energy++) = compute_energy(particle_x,0.0,script,opt); //fc=0 for spatial
			}
		}else if (script->coordinate_type==Temperature){
			*energy=compute_energy(particle_x,0.0,script,opt); //fc=0 for temperature
		}
		break;	
	}
	particle_x[2]=saved_x0;
	
	if(verbose_globalVar)fprintf(stderr,"sending %u bytes\n",send_size);

	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the energy file\n");
		exit(1);
	}
}

void send_sample_data(int sockfd, float force_data[], const struct script_struct *script){
	unsigned int send_size,data_size;
	int i;
	data_size=script->replica[0].sampling_steps*script->Nligands*sizeof(float);

	send_size=KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int)+data_size;

	unsigned char buff[send_size];
	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=TakeSampleData;
	
	*(unsigned int *)(buff+COMMAND_LOCATION+COMMAND_SIZE)=data_size;
	memcpy(buff+COMMAND_LOCATION+COMMAND_SIZE+sizeof(unsigned int),force_data,data_size);

	if(verbose_globalVar)fprintf(stderr,"sending %u bytes\n",send_size);

	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the sample data\n");
		exit(1);
	}
}

#define COORDINATE_DATA_SIZE 3
void send_coordinate_file(int sockfd){
	unsigned int send_size,data_size;
	int i;
	data_size=COORDINATE_DATA_SIZE*sizeof(float);

	send_size=KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int)+data_size;

	unsigned char buff[send_size];
	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=TakeCoordinateData;
	
	*(unsigned int *)(buff+COMMAND_LOCATION+COMMAND_SIZE)=data_size;
	for(i=0;i<COORDINATE_DATA_SIZE;i++) ((float *)(buff+COMMAND_LOCATION+COMMAND_SIZE+sizeof(unsigned int)))[i] = 0.0;

	if(verbose_globalVar)fprintf(stderr,"sending %u bytes\n",send_size);

	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the coordinate data\n");
		exit(1);
	}
}

void send_restart_file(int sockfd, double *restart_data, struct script_struct *script){
	unsigned int send_size;
	unsigned char Nfloats;
	int i;
	
	send_size=KEY_SIZE+COMMAND_SIZE+sizeof(unsigned int)+RST_FILE_SIZE*script->Nsamesystem_uncoupled;

	unsigned char buff[send_size];
	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=TakeRestartFile;
	
	*(unsigned int *)(buff+COMMAND_LOCATION+COMMAND_SIZE)=RST_FILE_SIZE*script->Nsamesystem_uncoupled;
	memcpy( buff+COMMAND_LOCATION+COMMAND_SIZE+sizeof(unsigned int), restart_data, RST_FILE_SIZE*script->Nsamesystem_uncoupled);

	fprintf(stderr,"SENDING RESTART:\t");                    //##DEBUG LENGTHY
	for(i=0;i<3*script->Nsamesystem_uncoupled; i++){         //##DEBUG LENGTHY
		fprintf(stderr,"%f\t",((double *)restart_data)[i]);    //##DEBUG LENGTHY
	}                                                        //##DEBUG LENGTHY
	fprintf(stderr,"\n");                                    //##DEBUG LENGTHY

	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the restart data\n");
		exit(1);
	}
}

enum command_enum receive_command(int sockfd){
	char buff[KEY_SIZE+COMMAND_SIZE];
	
	if( read(sockfd, buff, KEY_SIZE+COMMAND_SIZE)!=KEY_SIZE+COMMAND_SIZE )
		return(InvalidCommand);

	if( strncmp(buff+KEY_LOCATION, COMMAND_KEY, KEY_SIZE)!=0 ){
		fprintf(stderr,"Error: the key did not match\n");
		exit(1);
	}
	
	return (enum command_enum)buff[COMMAND_LOCATION];
} 

void receive_replica_ID(int sockfd, struct ID_struct *ID, const struct script_struct *script){
	read4K(sockfd,ID,sizeof(*ID));
	if(ID->replica_number>=script->Nreplicas){
		fprintf(stderr,"Error: received an invalid replica number: %d\n",ID->replica_number);
		exit(1);
	}
}

void receive_restart_file(int sockfd, double *restart_data, const struct script_struct *script){
	int file_size;

	read4K(sockfd,&file_size,sizeof(file_size));
	if(file_size!=RST_FILE_SIZE*script->Nsamesystem_uncoupled){
		fprintf(stderr,"Error: the restart file that was received (%d bytes) is not the right size (%u bytes)\n",file_size,(int)RST_FILE_SIZE*script->Nsamesystem_uncoupled);
		exit(1);
	}
	read4K(sockfd,restart_data,file_size);	
	int i;
	fprintf(stderr,"RECEIVING RESTART:\t");                  //##DEBUG LENGTHY
	for(i=0;i<3*script->Nsamesystem_uncoupled; i++){         //##DEBUG LENGTHY
		fprintf(stderr,"%f\t",((double *)restart_data)[i]);    //##DEBUG LENGTHY
	}                                                        //##DEBUG LENGTHY
	fprintf(stderr,"\n");                                    //##DEBUG LENGTHY

}

#define MAX_PARAMETERS_SIZE 500
void receive_simulation_parameters(int sockfd, double *x, double *xchange, double *fc, const struct script_struct *script){
	int file_size;
	char buff[MAX_PARAMETERS_SIZE+1];
	char *start;
	char *p;
	int i,j;
	unsigned char variable_mode;
	float *x2,*xchange2;          
	// CN notes that x2 and xchange2 don't appear to be used at all (except for verbose output) and that could be
	//   problematic since it is outputting something separately from where it really determines it.
	int sampling_steps=0;
	int rnd_seed;
	int nni;

	x2=(float *)malloc(script->Nsamesystem_uncoupled*sizeof(float));
	xchange2=(float *)malloc(script->Nsamesystem_uncoupled*sizeof(float));
	if(x2==NULL||xchange2==NULL){
		fprintf(stderr,"Error: unable to allocate memory for x2 or xchange2 in receive_simulation_parameters()\n");
		exit(1);
	}

	read4K(sockfd,&file_size,sizeof(file_size));
	if(file_size>MAX_PARAMETERS_SIZE){
		fprintf(stderr,"Error: the size of the simulation parameters is too large: %d\n",file_size);
		free(x2);free(xchange2);
		exit(1);
	}

	read4K(sockfd,buff,file_size);
	
	//buff[file_size]=0;
	fprintf(stderr,"CN: received a file containing the simulation parameters: [%s]\n",buff); //##DEBUG
	
	start=p=buff;
	variable_mode=0;

	for(i=0;i<file_size;i++){
		if(buff[i]=='\n' || buff[i]==' '){
			buff[i]=0;
			if(buff+i>start){
				//Note that it is ESSENTIAL that the longer ones are checked first
				if(strncmp(start,"wrefchange2",11)==0){
					variable_mode=4;nni=0;
				}else if(strncmp(start,"wrefchange",10)==0){
					variable_mode=3;nni=0;
				}else if(strncmp(start,"wref2",5)==0){
					variable_mode=2;nni=0;
				}else if(strncmp(start,"wref",4)==0){
					variable_mode=1;nni=0;
				}else if(strncmp(start,"force",5)==0){
					variable_mode=5;nni=0;
				}else if(strncmp(start,"sampNsteps",10)==0){
					variable_mode=6;nni=0;
				}else if(strncmp(start,"rnd",3)==0){
					variable_mode=7;nni=0;
				}else{
					switch(variable_mode)
					{
					case 1:	
						*(x+3*nni)=atof(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: wref(%d): %f\n",nni,*(x+3*nni));
						nni++;
						break;
					case 2:	
						*(x2+nni)=atof(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: wref2(%d): %f\n",nni,*(x2+nni)); 
						nni++;
						break;
					case 3:	
						*(xchange+nni)=atof(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: wrefchange(%d): %f\n",nni,*(xchange+nni)); 
						nni++;
						break;
					case 4:	
						*(xchange2+nni)=atof(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: wrefchange2(%d): %f\n",nni,*(xchange2+nni)); 
						nni++;
						break;
					case 5:	
						*(fc+nni)=atof(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: force(%d): %f\n",nni,*(fc+nni)); 
						nni++;
						break;
					case 6:	
						sampling_steps=atoi(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: sampNsteps: %i\n",sampling_steps); 
						break;
					case 7:	
						rnd_seed=atoi(start); 
						if(verbose_globalVar)fprintf(stderr,"PAR: rnd: %d\n",rnd_seed); 
						break;
					}
				}
			}
			start=buff+i+1;
		}
	}
	
	if(sampling_steps!=script->replica[0].sampling_steps){
		fprintf(stderr,"Error: the number of sampling steps requested from the server (%d) is not the expected number (%d)\n",sampling_steps,script->replica[0].sampling_steps);
		free(x2);free(xchange2);
		exit(1);
	}
	//CN thinks that there should be some test here for x2 based on x similar to the one above
	free(x2);free(xchange2);
}

void *client(struct client_bundle *arg){
	int sockfd;
	struct ID_struct ID;
	enum command_enum command;
	bool replica_done;
	bool client_done;
	double *particle_x;
	//particle_x contains 0,1,2 for Nsamesystem_uncoupled=1 then 0,1,2 for Nsamesystem_uncoupled=2 ...
	double *new_coordinate;
	float **sample_data=(float **)NULL;
	float **additional_data=(float **)NULL;
	double *fc;
	bool FirstTimeThrough=true;
	int i;

	struct script_struct *script=arg->script;
	struct tester_option_struct *opt=arg->opt;

	ID.title[0]=ID.title[1]='*';
	// valgrind complains about the first call to send_replica_ID(sockfd,&ID); around line 640 of this file
	// CN expected the line below to solve that problem, but it did not
	ID.replica_number=ID.sequence_number=ID.title[2]=ID.title[3]=0; 

	if(script->need_sample_data){
		sample_data=(float **)malloc(script->Nsamesystem_uncoupled*sizeof(float *));
		if(sample_data==NULL){
			fprintf(stderr,"Error: unable to allocate memory for sample_data in client().\n");fflush(stderr);
			exit(1);
		}
		for(i=0;i<script->Nsamesystem_uncoupled;i++){
			sample_data[i]=(float *)malloc(script->replica[0].sampling_steps*script->Nligands*sizeof(float));
			if(sample_data[i]==NULL){
				fprintf(stderr,"Error: unable to allocate memory for sample_data[i] in client().\n");fflush(stderr);
				exit(1);
			}
		}
	}

	if(script->Nadditional_data>0){
		//There is only a single additional data for the test and therefore the extra 
		//dimension to the array that would be present in the client is not even allocated here
		additional_data=(float **)malloc(script->Nsamesystem_uncoupled*sizeof(float *));
		if(additional_data==NULL){
			fprintf(stderr,"Error: unable to allocate memory for additional_data in client().\n");fflush(stderr);
			exit(1);
		}
		for(i=0;i<script->Nsamesystem_uncoupled;i++){
			additional_data[i]=(float *)malloc(script->replica[0].sampling_steps*sizeof(float));
			if(additional_data[i]==NULL){
				fprintf(stderr,"Error: unable to allocate memory for additional_data[i] in client().\n");fflush(stderr);
				exit(1);
			}
		}
	}

	if((particle_x=(double *)malloc(script->Nsamesystem_uncoupled*3*sizeof(double)))==NULL){
		fprintf(stderr,"Error: unable to allocate memory for particle_x\n");
		exit(1);
	}
	if((new_coordinate=(double *)malloc(script->Nsamesystem_uncoupled*sizeof(double)))==NULL){
		fprintf(stderr,"Error: unable to allocate memory for new_coordinate\n");
		exit(1);
	}
	if((fc=(double *)malloc(script->Nsamesystem_uncoupled*sizeof(double)))==NULL){
		fprintf(stderr,"Error: unable to allocate memory for fc\n");
		exit(1);
	}
	for(i=0;i<script->Nsamesystem_uncoupled*3;i++){
		particle_x[i]=0.0;
	}

	replica_done=false;
	while(!replica_done){
		if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0){
			fprintf(stderr,"Error: cannot open socket\n");
			if(script->need_sample_data){
				for(i=0;i<script->Nsamesystem_uncoupled;i++){
					free(sample_data[i]);
					if(script->Nadditional_data>0) free(additional_data[i]);
				}
				free(sample_data);
				if(script->Nadditional_data>0) free(additional_data);
			}
			free(particle_x);
			free(new_coordinate);
			free(fc);
			exit(1);
		}

		errno=0;
		if(connect(sockfd, (struct sockaddr *)&server_address, sizeof(struct sockaddr)) != -1){
			// first thing we do is send the protocol version we are using
			write(sockfd,&protocol_version,PROTOCOL_VERSION_SIZE); 
			fprintf(stderr,"sending replica ID\n");  //##DEBUG
			send_replica_ID(sockfd,&ID);
			
			if(ID.title[0]!='*'){
				for(i=0;i<script->Nsamesystem_uncoupled;i++){
					fprintf(stderr,"sending energy file\n"); //##DEBUG
					send_energy_file(sockfd,&(particle_x[i*3]),new_coordinate[i],script,opt);
					if(script->need_sample_data){
						fprintf(stderr,"sending force file\n"); //##DEBUG
						send_sample_data(sockfd,sample_data[i],script);
						fprintf(stderr,"FORCE DATA FOR %s%d.%d:\n",ID.title,ID.replica_number,ID.sequence_number);  //##DEBUG LENGTHY
						for(int zzz=0;zzz<script->replica[0].sampling_steps;zzz++){                                 //##DEBUG LENGTHY
							fprintf(stderr,"%f\n",sample_data[i][zzz]);                                               //##DEBUG LENGTHY
						}                                                                                           //##DEBUG LENGTHY
					}
					if(script->Nadditional_data>0){
						//This will cause an error with Nligands>1 probably
						fprintf(stderr,"sending additional data file\n"); //##DEBUG
						fprintf(stderr,"ADDITIONAL DATA FOR %s%d.%d:\n",ID.title,ID.replica_number,ID.sequence_number); //##DEBUG LENGTHY
						for(int zzz=0;zzz<script->replica[0].sampling_steps;zzz++){                                     //##DEBUG LENGTHY
							fprintf(stderr,"%f\n",additional_data[i][zzz]);                                               //##DEBUG LENGTHY
						}                                                                                               //##DEBUG LENGTHY
						//Only one additional data file is allowed with the tester so there is no loop required here
						send_sample_data(sockfd,additional_data[i],script);
					}
					if(script->need_coordinate_data){
						fprintf(stderr,"sending crd file\n"); //##DEBUG
						//this is just a file-passing place holder. The coords are not really sent
						send_coordinate_file(sockfd);
					}
					fprintf(stderr,"sending rst file\n"); //##DEBUG
					//script->Nsamesystem_uncoupled get sent together
					if(i==0){
						send_restart_file(sockfd, particle_x, script);
						fprintf(stderr,"CN (send restart): particle_x=[%f,%f,%f]\n",particle_x[0],particle_x[1],particle_x[2]); //##DEBUG
					}else{
						send_NextNonInteracting(sockfd);
					}
				}
			}

			client_done=false;
			while(!client_done){
				command=receive_command(sockfd);
				fprintf(stderr,"received a command: %hhu\n",command); //##DEBUG
				switch(command)
				{
				case ReplicaID:
					fprintf(stderr,"received RepicaID command\n");fflush(stderr); //##DEBUG
					receive_replica_ID(sockfd, &ID,script);
					fprintf(stderr,"The replica ID is %sw%d.%u\n",ID.title,ID.replica_number,ID.sequence_number);fflush(stderr); //##DEBUG
					break;
				case TakeRestartFile:
					fprintf(stderr,"received TakeRestartFile command\n");fflush(stderr); //##DEBUG
					receive_restart_file(sockfd, particle_x, script);
					fprintf(stderr,"CN (receive restart): particle_x=[%f,%f,%f]\n",particle_x[0],particle_x[1],particle_x[2]); //##DEBUG
					break;
				case TakeSimulationParameters:
					fprintf(stderr,"received TakeSimulationParameters command\n");fflush(stderr); //##DEBUG
					receive_simulation_parameters(sockfd,particle_x,new_coordinate,fc,script);
					fprintf(stderr,"CN (receive sim param): particle_x=[%f,%f,%f]\n",particle_x[0],particle_x[1],particle_x[2]); //##DEBUG
					client_done=true;
					break;
				case InvalidCommand:
					client_done=true;
					replica_done=true;
					break;
				default:
					fprintf(stderr,"received an unexpected command\n");
					if(script->need_sample_data){
						for(i=0;i<script->Nsamesystem_uncoupled;i++){
							free(sample_data[i]);
							if(script->Nadditional_data>0) free(additional_data[i]);
						}
						free(sample_data);
						if(script->Nadditional_data>0) free(additional_data);
					}
					free(particle_x);
					free(new_coordinate);
					free(fc);
					exit(1);
				}
			}
		}else{
			//fprintf(stderr,"Error: cannot connect to the server (errno=%d)\n",errno);
			perror("Error: cannot connect to the server");
			replica_done=true;
		}
		close(sockfd);

		if(FirstTimeThrough){
			if(script->coordinate_type==Spatial || script->coordinate_type==Umbrella){
				for(i=0;i<script->Nsamesystem_uncoupled;i++){
					particle_x[i*3+1]=particle_x[i*3+2]=particle_x[i*3+0];
				}
			}else if (script->coordinate_type==Temperature){
				//particle_x[0] is temperature
				//Nsamesystem_uncoupled must be 1 for temperature
				particle_x[1]=particle_x[2]=0.0;
			}
			FirstTimeThrough=false;
		}
		if(script->coordinate_type==Spatial){
			for(i=0;i<script->Nsamesystem_uncoupled;i++){
				particle_x[i*3+2]=particle_x[i*3+0];
			}
		}
		if(script->coordinate_type==Temperature){
			T=particle_x[0];
			kT=(8.31451*T/4184.0);
			B=(1/kT);
			fprintf(stderr,"MOO: T=%f\tB=1/kT=%f\n",T,B); //##DEBUG LENGTHY
		}
		// This sleep statment may be required in order to avoid the run stopping 
		// with dr.log indicating "Error: pthread_create failed"
		usleep((useconds_t)opt->sleepTime); 
		
		if(!replica_done){
			for(i=0;i<script->Nsamesystem_uncoupled;i++){
				//fprintf(stderr,"ACTUALLY RUNNING SIMULATION %s%d.%d nni = %d:\n",ID.title,ID.replica_number,ID.sequence_number,i);
				run_simulation(script->replica[ID.replica_number+i].sampling_steps,&(particle_x[i*3]),sample_data[i],additional_data[i],fc[i],script,opt);
			}
		}
	}
	
	pthread_mutex_lock(&replica_mutex);
	Nfinished_replicas+=script->Nsamesystem_uncoupled;
	pthread_mutex_unlock(&replica_mutex);
	if(script->need_sample_data){
		for(i=0;i<script->Nsamesystem_uncoupled;i++){
			free(sample_data[i]);
			if(script->Nadditional_data>0) free(additional_data[i]);
		}
		free(sample_data);
		if(script->Nadditional_data>0) free(additional_data);
	}
	free(particle_x);
	free(new_coordinate);
	free(fc);
	return((void *)NULL);
}

void outputExactFile(const struct script_struct *script, const struct tester_option_struct *opt){
  double exact[3], din;
  FILE *e, *ein;
  char exactfilename[100], linein[1001];
	int i;

  sprintf(exactfilename,"%s.exact",opt->title);
  e=fopen(exactfilename,"w");
  if(e!=NULL){
    if(opt->exactInputFile[0]!='\0'){
      ein=fopen(opt->exactInputFile,"r");
      if(ein!=NULL){
        while(fgets(linein,1000,ein)!=NULL){
          if(sscanf(linein,"%lf",&din)!=1)continue;
          exact[0]=exact[1]=exact[2]=din;
          fprintf(e,"%lf\t%lf\n",din,compute_energy(exact,0.0,script,opt));
        }
        fclose(ein);
      }else{
        fprintf(stderr,"Unable to open %s to read the exact input file (special mode).\n",opt->exactInputFile);
      }
    }else{
      if(script->coordinate_type==Spatial || script->coordinate_type==Umbrella){
        for(i=0;i<script->Nreplicas;i++){
          exact[0]=exact[1]=exact[2]=script->replica[i].w_nominal;
          fprintf(e,"%lf\t%lf\n",script->replica[i].w_nominal,compute_energy(exact,0.0,script,opt));
        }
      }else if(script->coordinate_type==Temperature){
        //No way to know what range is sampled (apart from knowing the energy function)
        //Therefore sample some points and the analysis program will later call with -e option anyway
        for(exact[0]=0.0; exact[0]<=8.7; exact[0]+=0.3){
          exact[1]=exact[2]=exact[0];  //exact[0] is not used
          fprintf(e,"%lf\t%lf\n",exact[0],compute_energy(exact,0.0,script,opt));
        }
      }
    }
    fclose(e);
  }else{
    fprintf(stderr,"Unable to open %s to write the exact file.\n",exactfilename);
  }
}

void showUsage(const char *c, const struct tester_option_struct *opt){
	fprintf(stderr,"Usage: %s IP-address script-file [-nsv]\n",c);
	fprintf(stderr,"OR:    %s localhost  script-file [-nsv]\n",c);
	fprintf(stderr,"       -n [int] include noise (default = %d)\n",opt->includeNoise);
	fprintf(stderr,"          ( =0) no noise\n");
	fprintf(stderr,"          (!=0) with noise\n");
	fprintf(stderr,"       -s [int] microseconds to sleep (default = %d)\n",opt->sleepTime);
	fprintf(stderr,"       -v [int] verbose (default = %d)\n",verbose_globalVar);
	fprintf(stderr,"          ( =0) not verbose\n");
	fprintf(stderr,"          (!=0) verbose (a debugging feature)\n");
	fprintf(stderr,"       -r [int] number of replicas to submit (default = %d)\n",opt->numtosubmit);
	fprintf(stderr,"                Negative is a flag for Nreplicas/Nsamesystem_uncoupled\n");
	fprintf(stderr,"       -e [string] SPECIAL USAGE (no actual test) specify the filename containing\n");
	fprintf(stderr,"                   positions for which the exact solution is desired. A file te.exact\n");
	fprintf(stderr,"                   will be written containing these values.\n");
	exit(1);
}

int parseCommandLine(int argc, char * const argv[], struct tester_option_struct *opt){
	int i;
	int gotn=0;
	int gots=0;
	int gotv=0;
	int gotr=0;
	int gote=0;

	opt->exactInputFile[0]='\0';

	if(argc<3){
		fprintf(stderr,"Error: the script filename and port were not both provided\n");
		return 1;
	}
	if(strlen(argv[2])<2){
		fprintf(stderr,"Error: script filename must have at least 2 characters (the first 2 are used for the title)\n");
		return 1;
	}
	opt->title[0]=argv[2][0];
	opt->title[1]=argv[2][1];
	opt->title[2]=0;

	for(i=4; i<argc; i+=2){
		if(argv[i-1][0]!='-'){
			// Should only be dashed arguments from this point on.
			fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[i-1]);
			return 1;
		}
		if(argv[i-1][1]=='n'){
			if(gotn){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->includeNoise=atoi(argv[i]);
			gotn=1;
		}else if(argv[i-1][1]=='s'){
			if(gots){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->sleepTime=atoi(argv[i]);
			gots=1;
		}else if(argv[i-1][1]=='v'){
			if(gotv){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			verbose_globalVar=atoi(argv[i]);
			gotv=1;
		}else if(argv[i-1][1]=='r'){
			if(gotr){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->numtosubmit=atoi(argv[i]);
			gotr=1;
		}else if(argv[i-1][1]=='e'){
			if(gote){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			sscanf(argv[i],"%s",opt->exactInputFile);
			gote=1;
		}else{
			fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[i-1]);
			return 1;
		}
	}
	if(i-1<argc){
		fprintf(stderr,"Error: incorrect command line format. Extraneous arguments begin at %s\n",argv[i-1]);
		return 1;
	}
	return 0;
}
	
int main(int argc, char *argv[]){
	//unsigned short port;
	//enum coordinate_type_enum coordinate_type;
	//unsigned int Nmoves=10000;
	//unsigned int protocol_version=PROTOCOL_VERSION;
	pthread_t server_handle;
	int i;

	struct script_struct script;
	struct tester_option_struct opt=DEFAULT_TESTER_OPTION_STRUCT;

	T=INITIAL_T;
	kT=(8.31451*T/4184.0);
	B=(1/kT);

	if(argc==1){
		showUsage(argv[0],&opt);
		exit(1);
	}

	parseCommandLine(argc,argv,&opt);

	if(opt.exactInputFile[0]=='\0'){
		fprintf(stderr,"*** STARTING TEST OF THE DISTRIBUTED REPLICA SYSTEM ***\n");
	}else{
		//Not a real run, the analysis tool is just asking for exact solution values
		fprintf(stderr,"*** obtaining exact solution values from the tester.\n");
	}

	read_input_script_file_class *input_script = new read_input_script_file_class;
	input_script->read_input_script_file(argv[2],&script);
	delete input_script;
	if(opt.numtosubmit<0){
		opt.numtosubmit=(int)ldiv(script.Nreplicas,script.Nsamesystem_uncoupled).quot;
	}

	if(script.coordinate_type!=Spatial && script.coordinate_type!=Umbrella && script.coordinate_type!=Temperature){
		fprintf(stderr,"Error: coordinate_type must be one of Spatial, Umbrella, or Temperature\n");
		exit(1);
	}

	if(script.Nsamples_per_run>1&&!opt.includeNoise){
		fprintf(stderr,"\nWARNING\nWARNING\nWARNING: It is a waste of time and disk space to "); 
		fprintf(stderr,"have more than 1 sample per run while not including noise. ");
		fprintf(stderr,"Recommend that you reduce your number of steps from %d to 1 ",script.Nsamples_per_run);
		fprintf(stderr,"in your script file or include noise via the command line flag.\n");
		fprintf(stderr,"Nevertheless, the run will not be terminated.\n");
		sleep(5);
	}

	if(script.Nsamples_per_run==1&&opt.includeNoise){
		fprintf(stderr,"ERROR: including noise requires more than one sample per run. ");
		fprintf(stderr,"Without more than one sample, noise is impossible. Turn off noise or ");
		fprintf(stderr,"increase your sample number per run in your script file\n");
		exit(1);
	}

	if(script.Nadditional_data>1){
		fprintf(stderr,"ERROR: the test suite can not currently accept more than 1 piece of ADDITIONAL_DATA\n");
		exit(1);
	}

	server_address.sin_family = AF_INET;    // host byte order 
	if(strcmp(argv[1],"localhost")==0) server_address.sin_addr.s_addr = inet_addr("127.0.0.1");
	else server_address.sin_addr.s_addr = inet_addr(argv[1]);
	server_address.sin_port = htons(script.port);  // short, network byte order 
	memset(&(server_address.sin_zero), 0, 8);  // zero the rest of the struct 

	pthread_mutex_init(&replica_mutex,NULL);
	srand48(3454545);

	outputExactFile(&script,&opt);

	if(opt.exactInputFile[0]!='\0'){
		pthread_mutex_destroy(&replica_mutex);
		exit(0);
	}

	struct client_bundle c_bundle;
	c_bundle.script=&script;
	c_bundle.opt=&opt;
	for(i=0;i<opt.numtosubmit;i++){
		if(pthread_create(&server_handle,NULL, (void* (*)(void*))client, &c_bundle)!=0){
			error_quit("Error: pthread_create failed in DR_tester");
		}
		if(pthread_detach(server_handle)!=0) error_quit("Error: pthread_detach failed");
	}

	while(1){
		sleep(1);
		
		pthread_mutex_lock(&replica_mutex);
		i=Nfinished_replicas;
		pthread_mutex_unlock(&replica_mutex);
		
		fprintf(stderr,"FINISHED %d REPLICAS\n",i);
		
		if(i==script.Nreplicas) break;
	}

	pthread_mutex_destroy(&replica_mutex);

	fprintf(stderr,"*** THE TEST IS COMPLETED ***\n");
}

