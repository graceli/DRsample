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

#ifndef _READ_INPUT_SCRIPT_FILE_H
#define _READ_INPUT_SCRIPT_FILE_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


#define BOLTZMANN_CONSTANT (8.31451/4184.0)

#define MAX_COLUMNS 9
#define MAX_PARAMETER_CHARACTERS 50
#define N_PRESENCE_BITS 100000

enum coordinate_type_enum {CoordinateTypeUndefined,Spatial,Temperature,Umbrella};
enum replica_move_type_enum {MoveTypeUndefined,MonteCarlo,BoltzmannJumping,Continuous,NoMoves,vRE};

struct buffer_struct{
	unsigned char *data;
	unsigned int data_size;
	unsigned int allocated_memory;
};

struct atom_struct{
	double x,y,z;
	unsigned int weight;
};

struct replica_struct{
	//w_nominal is fixed and stores the start position of each replica
	//w refers to the coordinate. w can be a spatial coordinate (which can be the fourth dimension), 
	//or w can be lambda, or can be beta if a temperature coordinate is used.
	char status;
	float w;
	float w_nominal;
	float w_start;
	float w2_nominal;
	float w_sorted;
	float force;
	unsigned int sequence_number;
	unsigned int sample_count;
	unsigned int sampling_runs;
	unsigned int sampling_steps;
	double cancellation_accumulator[2];
	unsigned short cancellation_count;
	float cancellation_energy;
	unsigned int last_activity_time;
	unsigned int start_time_on_current_node;
	struct buffer_struct restart;
	struct atom_struct *atom;
	unsigned int *presence;
	char vREfile[500];
	int nodeSlot;
};

struct script_struct{
	struct replica_struct *replica;
	unsigned int Nreplicas;
	coordinate_type_enum coordinate_type;
	enum replica_move_type_enum replica_move_type;
	float temperature;
	unsigned int Nligands;
	unsigned int node_time;
	unsigned int replica_change_time;
	unsigned int snapshot_save_interval;
	unsigned int job_timeout;
	unsigned int port;
	float replica_potential_scalar1;
	float replica_potential_scalar2;
	float replica_potential_scalar1_after_threshold;
	float replica_potential_scalar2_after_threshold;
	unsigned int cancellation_threshold;
	float replica_step_fraction;
	bool need_sample_data;
	bool need_coordinate_data;
	bool submit_jobs;
	unsigned int min_unsuspended_replica;
	unsigned int max_unsuspended_replica;
	bool circular_replica_coordinate;
	float circular_lesser_equality;
	float circular_greater_equality;
	float circular_equality_distance; //greater < lesser sets up failure condition on purpose
	unsigned int Nsamples_per_run;
	unsigned int Nadditional_data;
	unsigned int Nsamesystem_uncoupled;
	bool loadedCancel;
	bool stopOnAverageTimeExceeded;
	long int vRE_initial_noMoves;
	long int vRE_initial_noSave;
	long int vRE_secvre_size;
	bool allow_requeue;
	unsigned int allotted_time_for_server;
	bool defineStartPos;
	float cycleClients;
	int mobility_time;
	int mobility_requiredTimeGain;
};

class read_input_script_file_class{
private:
	//C. Neale changed it so that replica is no longer internal
	unsigned int Nreplicas;

	void error_quit(const char *error_message){
		fprintf(stderr,"Error: %s\n",error_message);
		fflush(stderr);
		exit(1);
	}

	void new_replica(replica_struct **replica, float w, float w2, float force, unsigned int sampling_runs, unsigned int sampling_steps, float cancel_energy, float startingNominal, const char *vreFile){
		struct replica_struct *temp;

		if(Nreplicas>0){
			temp=*replica;
		}
		*replica=new replica_struct[Nreplicas+1];
		if(Nreplicas>0){
			for(int i=0;i<Nreplicas;i++){
				(*replica)[i]=temp[i];
			}
			delete[] temp;
		}
		
		(*replica)[Nreplicas].status='N';
		(*replica)[Nreplicas].w=w;
		(*replica)[Nreplicas].w_nominal=w;
		if(isnan(startingNominal)){
//fprintf(stderr,"NAN setting start to %f\n",(*replica)[Nreplicas].w);
			(*replica)[Nreplicas].w_start=(*replica)[Nreplicas].w;
		}else{
//fprintf(stderr,"NON NAN Replacing %f with %f\n",(*replica)[Nreplicas].w,(*replica)[Nreplicas].w_start);
			(*replica)[Nreplicas].w_start=startingNominal;
			(*replica)[Nreplicas].w=(*replica)[Nreplicas].w_start;
		}
		(*replica)[Nreplicas].w2_nominal=w2;
		(*replica)[Nreplicas].w_sorted=0;   //set to zero to avoid writing undefined variable
		(*replica)[Nreplicas].force=force;
		(*replica)[Nreplicas].sequence_number=0;
		(*replica)[Nreplicas].sample_count=0;
		(*replica)[Nreplicas].sampling_runs=sampling_runs;
		(*replica)[Nreplicas].sampling_steps=sampling_steps;
		(*replica)[Nreplicas].cancellation_accumulator[0]=0.0;
		(*replica)[Nreplicas].cancellation_accumulator[1]=0.0;
		(*replica)[Nreplicas].cancellation_count=0;
		(*replica)[Nreplicas].cancellation_energy=cancel_energy;
		(*replica)[Nreplicas].last_activity_time=time(NULL);
		(*replica)[Nreplicas].start_time_on_current_node=time(NULL);
		(*replica)[Nreplicas].restart.data=NULL;
		(*replica)[Nreplicas].restart.data_size=0;
		(*replica)[Nreplicas].restart.allocated_memory=0;
		sprintf((*replica)[Nreplicas].vREfile,"%s",vreFile);
		(*replica)[Nreplicas].nodeSlot=-1;

		(*replica)[Nreplicas].atom=NULL;

		//This is pointed to by script->replica.atom when coordinate data is not saved.
		//Otherwise, valgrind memory checker complained and it took C. Neale many hours
		//to determine that the error was inconsequential
		//C. Neale now thinks that the problem was .status somehow, therefore changed back
		//(*replica)[Nreplicas].atom=new atom_struct[1];
		//(*replica)[Nreplicas].atom[0].x=0.0;
		//(*replica)[Nreplicas].atom[0].y=0.0;
		//(*replica)[Nreplicas].atom[0].z=0.0;
		//(*replica)[Nreplicas].atom[0].weight=0;

		//There is a memory leak related to .presence. However, it is relatively small so leave it alone	
		(*replica)[Nreplicas].presence=new unsigned int[N_PRESENCE_BITS/32];
		memset((*replica)[Nreplicas].presence,0,N_PRESENCE_BITS/8); //C. Neale changed divisor from 32 to 8
	
		Nreplicas++;
	}


	unsigned int count_parameters(char *line){
		unsigned int i;
		bool b, old_b;
		unsigned int count;
		
		old_b=false;
		i=0;
		count=0;
		while(line[i]!=0){
			b=(line[i]!=' ') && (line[i]!='\t');
			if(b && !old_b) count++;
			old_b=b;
			i++;
		}
		return(count);
	}
	
	void parse_line(char *line, unsigned char paramN, char *param){
		unsigned int i,j,k;

		//printf("Getting parameter number %u from line: [%s]\n",paramN,line);
	
		j=0;
		for(i=0;i<=paramN;i++){
			k=0;
			while( ( (line[j]==' ') || (line[j]=='\t') ) && (line[j]!=0) ) j++;
			//printf("START AT: [%s]\n",&line[j]);
			while( (line[j]!=' ') && (line[j]!='\t') && (line[j]!=0) && (k<MAX_PARAMETER_CHARACTERS) ){
				param[k]=line[j];
				j++;
				k++;
			}
			param[k]=0;
			//printf("PARMAETER: [%s]\n",param);
		}
	}

public:
	read_input_script_file_class(void){
		//C. Neale changed it so that replica is no longer internal
		Nreplicas=0;
	}

	void read_input_script_file(char *filename, struct script_struct *script){
		FILE *fd;
		char buffer[500];
		char command[MAX_PARAMETER_CHARACTERS+1];
		char param[MAX_PARAMETER_CHARACTERS+1];
		enum column_enum {ColumnW1,ColumnW2,ColumnForce,ColumnMoves,ColumnSteps,ColumnCancelEnergy,ColumnStartW1,vREfile} column[MAX_COLUMNS];
		signed char n_columns=-1;
		bool spec_node_time=false;
		bool spec_replica_change_time=false;
		bool spec_snapshot_save_interval=false;
		bool spec_job_timeout=false;
		bool spec_unsuspended_replica=false;

		script->replica=(replica_struct *)NULL;	//may cause problems if script is read multiple times
		script->coordinate_type=CoordinateTypeUndefined;
		script->replica_move_type=MoveTypeUndefined;
		script->temperature=-1;
		script->Nligands=1;
		script->node_time=0;
		script->replica_change_time=0;
		script->snapshot_save_interval=0;
		script->job_timeout=0;
		script->port=0;
		script->replica_potential_scalar1=-1.0;
		script->replica_potential_scalar2=-1.0;
		script->replica_potential_scalar1_after_threshold=0.0;
		script->replica_potential_scalar2_after_threshold=0.0;
		script->cancellation_threshold=0;
		script->replica_step_fraction=-1.0;
		script->need_sample_data=false;
		script->need_coordinate_data=false;
		script->submit_jobs=false;
		script->min_unsuspended_replica=0;
		script->max_unsuspended_replica=0;
		script->circular_replica_coordinate=false;
		script->circular_lesser_equality=1.0;
		script->circular_greater_equality=-1.0;
		script->circular_equality_distance=-1.0; //greater < lesser sets up failure condition on purpose
		script->Nsamples_per_run=0;
		script->Nadditional_data=0;
		script->Nsamesystem_uncoupled=1;
		script->loadedCancel=false;
		script->stopOnAverageTimeExceeded=false;
		script->vRE_initial_noMoves=0;
		script->vRE_initial_noSave=0;
		script->vRE_secvre_size=-1;
		script->allow_requeue=false;
		script->allotted_time_for_server=0;
		script->defineStartPos=false;
		script->cycleClients=-1.0;
		script->mobility_time=0;
		script->mobility_requiredTimeGain=0;
		
		if((fd=fopen(filename,"r"))==NULL){
			fprintf(stderr,"Error: cannot open input script file %s\n",filename);
			fflush(stderr);
			exit(1);
		}
	
		while (!feof(fd)){
			if(fgets(buffer, 499, fd)==NULL) continue;
			buffer[strlen(buffer)-1]=0;
	
			parse_line(buffer, 0, command);
			if(command[0]==0) continue;
			if( (command[0]=='/') && (command[1]=='/') ) continue;
	
			if(strcasecmp(command,"SIMULATION")==0){
				parse_line(buffer, 1, param);
				if(strcasecmp(param,"spatial")==0) script->coordinate_type=Spatial;
				else if(strcasecmp(param,"temperature")==0) script->coordinate_type=Temperature;
				else if(strcasecmp(param,"umbrella")==0) script->coordinate_type=Umbrella;
	
				parse_line(buffer, 2, param);
				if(strcasecmp(param,"montecarlo")==0) script->replica_move_type=MonteCarlo;
				else if(strcasecmp(param,"boltzmann")==0) script->replica_move_type=BoltzmannJumping;
				else if(strcasecmp(param,"continuous")==0) script->replica_move_type=Continuous;
				else if(strcasecmp(param,"nomoves")==0) script->replica_move_type=NoMoves;
				else if(strcasecmp(param,"vre")==0) script->replica_move_type=vRE;
			}else if(strcasecmp(command,"PORT")==0){
				sscanf(buffer,"%*s %u",&(script->port));
			}else if(strcasecmp(command,"TEMPERATURE")==0){
				sscanf(buffer,"%*s %f",&(script->temperature));
			}else if(strcmp(command,"REPLICASTEP")==0){
				sscanf(buffer,"%*s %f",&(script->replica_step_fraction));
			}else if(strcasecmp(command,"POTENTIALSCALAR")==0){
				sscanf(buffer,"%*s %f %f",&(script->replica_potential_scalar1), &(script->replica_potential_scalar2));
			}else if(strcasecmp(command,"CANCELLATION")==0){
				sscanf(buffer,"%*s %f %f %u",&(script->replica_potential_scalar1_after_threshold), &(script->replica_potential_scalar2_after_threshold), &(script->cancellation_threshold));
			}else if(strcasecmp(command,"NODETIME")==0){
				sscanf(buffer,"%*s %d",&(script->node_time));
				spec_node_time=true;
			}else if(strcasecmp(command,"REPLICACHANGETIME")==0){
				sscanf(buffer,"%*s %d",&(script->replica_change_time));
				spec_replica_change_time=true;
			}else if(strcasecmp(command,"SNAPSHOTTIME")==0){
				sscanf(buffer,"%*s %d",&(script->snapshot_save_interval));
				spec_snapshot_save_interval=true;
			}else if(strcasecmp(command,"TIMEOUT")==0){
				sscanf(buffer,"%*s %d",&(script->job_timeout));
				spec_job_timeout=true;
			}else if(strcasecmp(command,"RUNNINGREPLICAS")==0){
				sscanf(buffer,"%*s %d %d",&(script->min_unsuspended_replica), &(script->max_unsuspended_replica));
				spec_unsuspended_replica=true;
			}else if(strcasecmp(command,"NEEDSAMPLEDATA")==0){
				script->need_sample_data=true;
			}else if(strcasecmp(command,"NEEDCOORDINATEDATA")==0){
				script->need_coordinate_data=true;
			}else if(strcasecmp(command,"SUBMITJOBS")==0){
				script->submit_jobs=true;
			}else if(strcasecmp(command,"CIRCULAR")==0){
				if(sscanf(buffer,"%*s %f %f",&(script->circular_lesser_equality),&(script->circular_greater_equality))==2){
					script->circular_replica_coordinate=true;
					script->circular_equality_distance=script->circular_greater_equality-script->circular_lesser_equality;
				}
			}else if(strcasecmp(command,"ADDITIONALDATA")==0){
				sscanf(buffer,"%*s %u",&(script->Nadditional_data));
			}else if(strcasecmp(command,"N_SAMESYSTEM_UNCOUPLED")==0){
				sscanf(buffer,"%*s %u",&(script->Nsamesystem_uncoupled));
			}else if(strcasecmp(command,"STOP_ON_AVERAGE_TIME_EXCEEDED")==0){
				script->stopOnAverageTimeExceeded=true;
			}else if(strcasecmp(command,"VRE_INITIAL_NOMOVES")==0){
				sscanf(buffer,"%*s %ld",&(script->vRE_initial_noMoves));
                        }else if(strcasecmp(command,"VRE_INITIAL_NOSAVE")==0){
                                sscanf(buffer,"%*s %ld",&(script->vRE_initial_noSave));
			}else if(strcasecmp(command,"VRE_SECONDARY_LIST_LENGTH")==0){
                                sscanf(buffer,"%*s %ld",&(script->vRE_secvre_size));
			}else if(strcasecmp(command,"ALLOW_REQUEUE")==0){
				script->allow_requeue=true;
			}else if(strcasecmp(command,"ALLOTTED_TIME_FOR_SERVER")==0){
				sscanf(buffer,"%*s %u",&(script->allotted_time_for_server));
			}else if(strcasecmp(command,"DEFINE_STARTING_POSITIONS")==0){
				script->defineStartPos=true;
			}else if(strcmp(command,"CYCLE_CLIENTS")==0){
				sscanf(buffer,"%*s %f",&(script->cycleClients));
			}else if(strcmp(command,"SERVER_TIMELEFT_ENTER_MOBILE_STATE")==0){
				sscanf(buffer,"%*s %d",&(script->mobility_time));
			}else if(strcmp(command,"SERVER_TIMEGAIN_ENTER_MOBILE_STATE")==0){
				sscanf(buffer,"%*s %d",&(script->mobility_requiredTimeGain));
			}else if(strcasecmp(command,"COLUMNS")==0){
				bool W1_defined=false;
				if(n_columns!=-1){
					fprintf(stderr,"Error: ""COLUMN"" was specified more than once\n");
					exit(1);
				}
				n_columns=count_parameters(buffer)-1;
				if(n_columns>MAX_COLUMNS){
					fprintf(stderr,"Error: there are too many parameters specified for ""COLUMN""\n");
					exit(1);
				}
				for(int i=0;i<n_columns;i++){
					parse_line(buffer, i+1, param);
					if( (strcasecmp(param,"LIGAND1")==0) || (strcasecmp(param,"TEMPERATURE")==0) ){
						column[i]=ColumnW1;
						W1_defined=true;
					}
					else if(strcasecmp(param,"LIGAND2")==0) column[i]=ColumnW2;
					else if(strcasecmp(param,"FUNNEL")==0) column[i]=ColumnW2;
					else if(strcasecmp(param,"FORCE")==0) column[i]=ColumnForce;
					else if(strcasecmp(param,"MOVES")==0) column[i]=ColumnMoves;
					else if(strcasecmp(param,"STEPS")==0) column[i]=ColumnSteps;
					else if(strcasecmp(param,"CANCEL")==0){
						column[i]=ColumnCancelEnergy;
						script->loadedCancel=true;
					}else if(strcasecmp(param,"STARTL1")==0) column[i]=ColumnStartW1;
					else if (strcasecmp(param,"VREFILE")==0) column[i]=vREfile;
					else{
						fprintf(stderr,"Error: unexpected parameter in the following line: %s\n",command);
						exit(1);
					}
				}
				if(!W1_defined) error_quit("one of the columns must give the coordinate position of the replica");
			}
			else if(strcasecmp(command,"JOB")==0){
				float w;
				float w2=NAN;
				float force=NAN;
				float cancel_energy=0.0;
				unsigned int sampling_runs=1;
				unsigned int sampling_steps=1;
				float startw1=NAN;
				char vreFile[500];
				vreFile[0]='\0';
	
				if(n_columns==-1){
					fprintf(stderr,"Error: ""COLUMNS"" needs to be specified before ""JOB""\n");
					exit(1);
				}
				
				if(n_columns!=count_parameters(buffer)-1){
					fprintf(stderr,"Error: the number of parameters in the ""JOB"" definition does not match that defined by ""COLUMNS""\n");
					exit(1);
				}
				
				for(int i=0;i<n_columns;i++){
					parse_line(buffer, i+1, param);
					switch(column[i])
					{
					case ColumnW1:
						w=atof(param);
						break;
					case ColumnW2:
						w2=atof(param);
						script->Nligands=2;
						break;
					case ColumnForce:
						force=atof(param);
						break;
					case ColumnMoves:
						sampling_runs=atoi(param);
						break;
					case ColumnSteps:
						sampling_steps=atoi(param);
						break;
					case ColumnCancelEnergy:
						cancel_energy=atof(param);
						break;
					case ColumnStartW1:
						startw1=atof(param);
						break;	
					case vREfile:
						sprintf(vreFile,"%s",param);
						break;
					}
				}
				new_replica(&(script->replica), w, w2, force, sampling_runs, sampling_steps, cancel_energy, startw1,vreFile);
			}else{
				fprintf(stderr,"Error: Extraneous command found in input script: [%s]\n",buffer);
				exit(1);
			}
		}
	
		fclose(fd);
		
		if(script->coordinate_type==CoordinateTypeUndefined) error_quit("need to specify the type of simulation using the SIMULATION key word");
		else if( (script->coordinate_type==Spatial) || (script->coordinate_type==Umbrella) ){
			if(script->temperature<0) error_quit("Need to specify a temperature for a spatial or umbrella simulation");	
		}else{
			if(script->temperature!=-1.0) error_quit("It is nonsensical to specify a TEMPERATURE if the coordinate is temperature");	
		}
		if(script->replica_move_type==MoveTypeUndefined) error_quit("need to specify the move type using the SIMULATION key word");
		else if(script->replica_move_type==MonteCarlo||script->replica_move_type==vRE){
			if(script->replica_step_fraction<=0.0) error_quit("the REPLICASTEP should be greater than zero to do Monte Carlo or vRE moves");
		}else{
			if(script->replica_step_fraction!=-1.0) error_quit("the REPLICASTEP should only be specified for Monte Carlo or vRE simulations");
		}
		if( (Nreplicas==1) && (script->replica_move_type!=NoMoves) ) error_quit("Moves cannot be performed if the simulation has only one replica");
		if( (script->coordinate_type==Spatial) && (script->replica_move_type==Continuous) ) error_quit("a ""Spatial"" simulation is not currently compatible with the continuous boltzmann jumping method");
		if( (script->coordinate_type==Temperature) && (script->Nligands!=1) ) error_quit("a ""Temperature"" simulation is compatible with one ligand only");
		if( (script->coordinate_type==Umbrella) && (script->Nligands!=1) ) error_quit("an ""Umbrella"" simulation is compatible with one ligand only");
		if(script->vRE_initial_noMoves!=0&&script->replica_move_type!=vRE) error_quit("the ""vRE_initial_nomoves"" option is only compatible with a ""vRE"" simulation.\n");
                if(script->vRE_initial_noSave!=0&&script->replica_move_type!=vRE) error_quit("the ""vRE_initial_save"" option is only compatible with a ""vRE"" simulation.\n");
                if(script->vRE_secvre_size!=-1&&script->replica_move_type!=vRE) error_quit("the ""vRE_secvre_size"" option is only compatible with a ""vRE"" simulation.\n");
		if(script->replica_move_type==vRE && script->defineStartPos){
			fprintf(stderr,"Warning: using vRE and defineStartPos together will require you to load in an initial VRE list or else a replica will never be able to move occupy a previously unoccupied nominal position!\n");
		}
		if(script->port<1) error_quit("Must specify a valid port between 1 and 65535");
		if(Nreplicas==0) error_quit("No replicas specified");
		if( (script->replica_potential_scalar1<0) || (script->replica_potential_scalar2<0) ) error_quit("must specify two potential scalars both greater than or equat to 0");
		if( (script->replica_potential_scalar1_after_threshold<0) || (script->replica_potential_scalar2_after_threshold<0) ) error_quit("invalid 'after threshold' potential scalars");
		if( (script->cancellation_threshold>0) && (script->coordinate_type==Spatial) && !(script->need_sample_data) ) error_quit("cannot do energy cancellation without sample data; please specify NEEDSAMPLEDATA in the script file");
		if(!spec_node_time || !spec_replica_change_time || !spec_snapshot_save_interval || !spec_job_timeout) error_quit("must specify node time, replica change time, snapshot save interval, and job timeout");
		if(script->circular_replica_coordinate){
			if(script->coordinate_type==Temperature) error_quit("Circular replica coordinate is nonsensical with temperature replicas");
			if(Nreplicas==1) error_quit("Circular replica coordinate is nonsensical with only one replica");
			if(Nreplicas==2) error_quit("Circular replica coordinate is unnecessary with only two replicas");
			if(script->replica_move_type==NoMoves) error_quit("Circular replica coordinate is nonsensical without exchanges (REPLICASTEP was 0.0)");
			if( (script->replica_potential_scalar2!=0.0) || (script->replica_potential_scalar2_after_threshold!=0.0) ) error_quit("for CIRCULAR simulations, the second parameters of POTENTIALSCALAR and CANCELLATION should be 0.0");
			if(script->Nligands!=1) error_quit("Circular replica coordinate does not currently support more than one ligand");
			if(script->replica_step_fraction>=1000.0) error_quit("Circular replica coordinate does not currently support Boltzmann jumping. T. Rodinger thinks it should work. C. Neale is worried about the image convention not being currently implemented there");
			//TR says: Chris, I think that circular should work with boltzmann jumping no problem
			if(script->circular_equality_distance<0.0) error_quit("Circular replica requires that circular_greater_equality>circular_lesser_equality. These numbers must be attainable by your simulation and represent the same point in coordinate space. For example, dihedral sampling with w_nominal positions at 0,10,20,...330,340,350 would require circular_greater_equality=355 and circular_lesser_equality=-5. More specifically, it is essential that replica[0].w_nominal-(replica[1].w_nominal-replica[0].w_nominal)*0.5<=circular_lesser_equality) && replica[Nreplicas-1].w_nominal+(replica[Nreplicas-1].w_nominal-replica[Nreplicas-2].w_nominal)*0.5>=circular_greater_equality");
			if(script->replica_move_type != MonteCarlo && 
			   (script->replica[1].w_nominal-script->replica[0].w_nominal != script->replica[Nreplicas-1].w_nominal-script->replica[Nreplicas-2].w_nominal ||
			    script->replica[1].w_nominal-script->replica[0].w_nominal != script->replica[0].w_nominal-script->replica[Nreplicas-1].w_nominal+script->circular_equality_distance
				 )
			  ){
					error_quit("In order to use a circular coordinate with anything other than MonteCarlo moves, the difference between the first two nominal positions must equal the difference between the last two nominal positions. Further, these must both equal the circularized difference between the first and last nominal positions.\n");
			}
			if(script->replica[0].w_nominal-(script->replica[1].w_nominal-script->replica[0].w_nominal)>=script->circular_lesser_equality || script->replica[Nreplicas-1].w_nominal+(script->replica[Nreplicas-1].w_nominal-script->replica[Nreplicas-2].w_nominal)<=script->circular_greater_equality) error_quit("Circular replica settings will not generate any first-to-last moves. It is essential that replica[0].w_nominal-(replica[1].w_nominal-replica[0].w_nominal)<circular_lesser_equality) && replica[Nreplicas-1].w_nominal+(replica[Nreplicas-1].w_nominal-replica[Nreplicas-2].w_nominal)>circular_greater_equality");
		}

		if(script->allotted_time_for_server>0 && script->mobility_time>script->allotted_time_for_server){
			// default mobility_time=0
			// default allotted_time_for_server=0
			error_quit("In order to use the mobile server, you must also specify an allotted_time_for_server. In addition, the mobility_time must be less than the allotted_time_for_server.\n");
		}
		if(script->allotted_time_for_server>0 && script->mobility_requiredTimeGain>script->allotted_time_for_server){
			error_quit("You selected a mobile server, but it can never go mobile since mobility_requiredTimeGain > allotted_time_for_server.\n");
		}
	
		if(script->coordinate_type==Temperature){
			for(int i=0;i<Nreplicas-1;i++){
				if(script->replica[i].w_nominal-script->replica[i+1].w_nominal<0.011) error_quit("replica temperature must be unique and go in descending order");
			}
		}else{
			for(int i=0;i<Nreplicas-1;i++){
				if(script->replica[i+1].w_nominal-script->replica[i].w_nominal<0.011) error_quit("replica w coordinates must be unique and go in ascending order");
			}
		}
		for(int i=0;i<Nreplicas-1;i++){
			if(script->replica[i].w_nominal!=script->replica[i].w_start){
				fprintf(stderr,"NOTE: replica %d defines nominal %f and yet it will start as %f (as requested in your script file).\n",i,script->replica[i].w_nominal,script->replica[i].w_start);
			}
		}
		script->Nsamples_per_run=script->replica[0].sampling_steps;
		for(int i=1;i<Nreplicas;i++) if(script->replica[i].sampling_steps!=script->Nsamples_per_run) error_quit("replicas must all have the same number of sample steps per run");
		if(script->coordinate_type==Temperature){
			for(int i=0;i<Nreplicas;i++){
				script->replica[i].w_nominal=(float)1.0/(script->replica[i].w_nominal*BOLTZMANN_CONSTANT);
				script->replica[i].w=script->replica[i].w_nominal;
			}
		}
		if(!spec_unsuspended_replica){
			script->min_unsuspended_replica=0;
			script->max_unsuspended_replica=Nreplicas-1;
		}
		if(script->max_unsuspended_replica>Nreplicas-1){
			script->max_unsuspended_replica=Nreplicas-1;
		}

		if(script->Nsamesystem_uncoupled!=1){
			if(script->Nsamesystem_uncoupled==0) error_quit("It is nonsensical to use DR with N_SAMESYSTEM_UNCOUPLED equal to zero. This would mean that you don't have a reaction coordinate at all.\n");
			if(script->coordinate_type==Temperature) error_quit("It is nonsensical to use DR with N_SAMESYSTEM_UNCOUPLED greater than one while using TEMPERATURE as a reaction coordinate.\n");
			if(ldiv(Nreplicas,script->Nsamesystem_uncoupled).rem!=0) error_quit("Currently, N_SAMESYSTEM_UNCOUPLED is only supported when the number of non-interacting samples is an exact integer factor of the number of nominal positions.\n");
			if(script->min_unsuspended_replica!=0||script->max_unsuspended_replica!=Nreplicas-1) error_quit("It is not possible to suspend any replicas while using N_SAMESYSTEM_UNCOUPLED != 1.\n");
			for(int i=1;i<Nreplicas;i++){
				if(script->replica[i].sampling_steps!=script->replica[0].sampling_steps) error_quit("The number of sampling steps at each nominal position must be equal when using N_SAMESYSTEM_UNCOUPLED != 1.\n");
				if(script->replica[i].sampling_runs!=script->replica[0].sampling_runs) error_quit("The number of sampling runs at each nominal position must be equal when using N_SAMESYSTEM_UNCOUPLED != 1.\n");
				//CN notes that technically the number of sampling runs need not be the same, but this allows us to
				//only deal with the sampling_runs of the first nni
				//I may have already set everything up so that it would work, but for now: disallow it
			}
		}
		
		script->Nreplicas=Nreplicas;
	}
};

#endif /* read_input_script_file.h */
