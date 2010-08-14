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

/* CN is not sure if this still stands???
 * C. Neale is currently editing this file
 * He has reached the read_in_database section and it is a bit of a quagmire
 * Currently he is moving the records to a section of the db->records
 * but that has not been updated yet in the struct_database
 */

/*
 * This program analyses a force database file and produces a postscript file as output
 */

#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "force_database_class.h"
#include "read_input_script_file.h"

#define N_FORCE_POINTS 9  // this should be an odd number
#define N_ENERGY_POINTS 101
#define AVERAGING_WINDOW 10
#define EQUILIBRATION_FRACTION 0.0 // this fraction of data is considered equilibration and is not taken into account in the average force calculation on page 1 and 2
#define MAX_REPLICAS 1000
#define MAX_FILENAME_LENGTH 30
#define MAX_SEQUENCE_NUMBER 50000  //was 30000
#define NUMVALUES 10  //number of values output to define the y-axis
#define MAX_DATABASE_SIZE 2000000000

#define BOLTZMANN_CONSTANT (8.31451/4184.0)

#define DASHON "[1 2] 0 setdash\n"
#define DASHONSTRONG "[1 3] 0 setdash\n"
#define DASHONWEAK "[1 1] 0 setdash\n"
#define DASHOFF "[] 0 setdash\n"

struct plotinfo_struct{
	float scale_x;
	float scale_y;
	float translate_x;
	float translate_y;
	float yscalemult;
	char xtitle[500];
	char ytitle[500];
	char title[500];
	int xa,xb,xc; //variables for looping over xgrid
	int ya,yb,yc; //variables for looping over ygrid
	unsigned char ligand_number;
	float *xlabels; //set to null to label x according to nominal positions or send labels
};
#define EMPTY_PLOTINFO {0.0,0.0,0.0,0.0,0.0,"","","",0,0,0,0,0,0,0,(float *)NULL}

struct pageinfo_struct{
	float width;
	float height;
	unsigned int font_size;
	unsigned int title_font_size;
	unsigned int numvalues;
	bool showGrid;
};
#define FULLSIZE_PAGEINFO {10.5,8.0,8,24,NUMVALUES,true}
#define SMALLSIZE_PAGEINFO {5.0,3.5,8,8,NUMVALUES,true}

void setPageinfoFull(struct pageinfo_struct *a){
	a->width=10.5;
	a->height=8.0;
	a->font_size=8;
	a->title_font_size=24;
	a->numvalues=NUMVALUES;
	a->showGrid=true;
}

void setPageinfoSmall(struct pageinfo_struct *a){
	a->width=5.0;
	a->height=3.5;
	a->font_size=8;
	a->title_font_size=8;
	a->numvalues=NUMVALUES;
	a->showGrid=true;
}

struct sampleDensity_struct{
	unsigned int **b;
	float min;
	float max;
	unsigned int Ncol;      //number of histogram bins
	float binWidth;
};
#define EMPTY_SAMPLEDENSITY {(unsigned int **)NULL,0.0,0.0,0,0.0}

struct pmf_struct{
	float *f;
	unsigned int *n;
	float min;
	float max;
	unsigned int Ncol;      //number of histogram bins
	float binWidth;
};
#define EMPTY_PMF {(float *)NULL,(unsigned int *)NULL,0.0,0.0,0,0.0}

struct exact_struct{
	float *e;		//the exact pmf value
	float *b;		//the bin
	unsigned int n;		//num bins
};
#define EMPTY_EXACT {(float *)NULL,(float *)NULL,0}

struct detailedBalance_struct{
	float w;
	double moveUp;		//probability
	double moveDown;	//probability
	double moveSame;	//probability
	long numSteps;
};
#define EMPTY_DETAILEDBALANCE {0.0,0.0,0.0,0.0,0}

struct analysis_option_struct{
	char title[3];
	int sequence_number_limit;
	int writeDatabase;
	int writeDatabaseFull;
	char writeDatabaseName[100];
	int userAskedForFitting;
	int verbose;
	int useCancellation;
	int additionalDataWithSampling;
	int histNcol;
	int discardOutside;
	int useExact;
	float equilFraction;
	int justwriteDatabase;
};
#define DEFAULT_ANALYSIS_OPTION_STRUCT {"",-1,0,1,"",0,0,1,1,0,1,1,EQUILIBRATION_FRACTION,0}

struct stats_struct{
	unsigned int max_time;
	unsigned int max_sequence_number;
	float max_samples;
	float max_force[2],min_force[2];
};
#define EMPTY_STATS_STRUCT {0,0,0,{0,0},{0,0}}

struct nominal_struct{
	float w[2];
};

#if(N_FORCE_POINTS>(N_ENERGY_POINTS+AVERAGING_WINDOW))
	#define GRAPH_ARRAY N_FORCE_POINTS
#else
	#define GRAPH_ARRAY N_ENERGY_POINTS+AVERAGING_WINDOW
#endif

struct graph_struct{
	float w[2];
	float replica_position[MAX_SEQUENCE_NUMBER];
	float rc_position[MAX_SEQUENCE_NUMBER];
	double point[2][GRAPH_ARRAY];
	double weight_sum[GRAPH_ARRAY];
	double average[2];
	double avg_weight_sum[2];
	double samples;
};


void showUsage(const char *c, const struct analysis_option_struct *opt);
int parseCommandLine(int argc,char * const argv[], struct analysis_option_struct *opt);
int checkOptions(struct analysis_option_struct *opt, const struct script_struct *script);
//Main

void print_postscript_forceAverage(unsigned char ligand_number, const struct graph_struct *graph, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct stats_struct *stats, const struct analysis_option_struct *opt);
void print_postscript_trajectory(int printSelection,int circularFlag, int ignorerounds, const struct pageinfo_struct *pageinfo, int plotRCinstead, struct detailedBalance_struct *detbal, const struct database_struct *db, const struct graph_struct *graph, const struct analysis_option_struct *opt, const struct script_struct *script, const struct stats_struct *stats, int storecode);
void print_postscript_sequenceDensity(const struct pageinfo_struct *pageinfo, unsigned int * const *sequenceDensity, const struct script_struct *script, const struct graph_struct *graph);
void print_postscript_sampleDensity(const struct pageinfo_struct *pageinfo, const struct sampleDensity_struct *sd, const struct script_struct *script, const struct graph_struct *graph);
void print_postscript_pmf(unsigned char ligand_number, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const exact_struct *exact, const struct pmf_struct *pmf, const struct graph_struct *graph);
void print_postscript_dGandA(unsigned char ligand_number, const float *cancel, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct graph_struct *graph);
void print_postscript_dGminusA(unsigned char ligand_number, const float *cancel, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct graph_struct *graph);
void print_postscript_samplingOverlap(int whichData, const struct pageinfo_struct *pageinfo, const struct script_struct *script, struct database_struct *db, const struct nominal_struct *nominal, const struct graph_struct *graph, const struct analysis_option_struct *opt, const struct stats_struct *stats);
void print_postscript_Pexchange(const struct pageinfo_struct *pageinfo, struct detailedBalance_struct *detbal, const struct script_struct *script, const struct graph_struct *graph);
void print_postscript_Pupdown(const struct pageinfo_struct *pageinfo,struct detailedBalance_struct *detbal, const struct script_struct *script, const struct graph_struct *graph);

void showPlot(const struct plotinfo_struct *plot, const struct pageinfo_struct *page, const struct script_struct *script,const struct graph_struct *graph);
int fileExists(const char *filename);
signed char compare_records(struct record_struct *r1, struct record_struct *r2);
void quicksort_database(struct database_struct *db, int left, int right);
struct record_struct* rec(unsigned int record_number,const struct database_struct *db);
int read_in_database(const struct analysis_option_struct *opt, const struct script_struct *script, struct database_struct *db, const struct nominal_struct *nominal);
char get_next_force(float w[2], float force[2], float *rc, int whichRC, unsigned int *time, int *replica_number, unsigned short *sequence_number, const struct database_struct *db);
void get_data_statistics(struct stats_struct *stats, const struct database_struct *db);
void condense_forces(unsigned int N_values_to_average, struct nominal_struct *nominal, struct graph_struct *graph, const struct script_struct *script, const analysis_option_struct *opt, struct stats_struct *stats, const struct database_struct *db);
exact_struct * getExactFromFile(const char *title, exact_struct *exact, const struct script_struct *script);
int getCancellationFromLog(const char *title, float *cancel, const struct script_struct *script);
int find_bin_from_w(double w, int l, const struct nominal_struct *nominal, const struct script_struct *script);
int round_up(float x);
int round_down(float x);
void rot_trans_regular(void);
char *get_RGB_colour(int colour_num);
void setup_regular(unsigned char *current_page);
void calcSequenceDensity(unsigned int **sequenceDensity, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const analysis_option_struct *opt, const struct stats_struct *stats);
void getMinMaxSampleDensity(float *min, float *max, unsigned int whichData, unsigned int whichReplica, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script);
int calcSampleDensity(struct sampleDensity_struct *sd, unsigned int whichData, const struct graph_struct *graph, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const struct analysis_option_struct *opt, const struct stats_struct *stats);
int calcPMF(struct pmf_struct *pmf, unsigned int whichData, const struct graph_struct *graph, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const struct analysis_option_struct *opt);
void setFont(unsigned int size);
int floatEqual(float i,float j);
void showFirstComboPageText(int px, int py);
void showSampleDensityNotAvailable(int px,int py,const struct script_struct *script,const struct analysis_option_struct *opt);
void showSecondComboPageText(int px,int py);


void showUsage(const char *c, const struct analysis_option_struct *opt){
	printf("This program creates .ps graphs based on forcedatabase.\n");
	//verbose option hidden from [list]
	printf("Usage: %s tt.script [-lcdtfahme] > analysis.ps\n",c);
	printf("       -l [int] sequence-number-limit; negative indicates no limit (default = %d)\n",opt->sequence_number_limit);
	printf("       -c [int] plot cancellation data from tt.log file (default = %d)\n",opt->useCancellation);
	printf("          ( =0) do not attempt\n");
	printf("          (!=0) plot cancellation (values must exist in tt.log file)\n");
	printf("       -d [string] text-database-name (default = not created)\n");
	printf("       -t [int] text-database-type (default = %d)\n",opt->writeDatabaseFull);
	printf("          ( =0) only output replica#, sequence#, w\n");
	printf("          (!=0) also output all forces and additional_data\n");
	printf("       -j [int] set nonzero to only write the text database but not produce a ps file (must use with -d).\n");
	printf("       -f [int] do-fitting-convergence-analysis (default = %d)\n",opt->userAskedForFitting);
	printf("          ( =0) do not attempt\n");
	printf("          (!=0) do fitting (user beware)\n"); 
	printf("       -a [int0 or >=1] which additional_data holds your sampling value; zero indicates not to use (default = %d)\n",opt->additionalDataWithSampling); 
	printf("       -h [int 0 or >=1] number of histogram bins for customizable plots; use zero for = Nreplicas (default = %d)\n",opt->histNcol);
	printf("       -m [int] discard data outside JOB range for for customizable plots (default = %d)\n",opt->discardOutside);
	printf("          ( =0) do not discard\n");
	printf("          (!=0) discard (useful for comparisons using DR_tester)\n");
	printf("       -e [real] initial fraction of data to discard (default = %f)\n",opt->equilFraction);
}

int parseCommandLine(int argc,char * const argv[], struct analysis_option_struct *opt){
	int i;
	int gotl=0;
	int gotd=0;
	int gott=0;
	int gotf=0;
	int gotv=0;
	int gotc=0;
	int gota=0;
	int goth=0;
	int gotm=0;
	int gote=0;
	int gotj=0;

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

	for(i=3; i<argc; i+=2){
		if(argv[i-1][0]!='-'){
			fprintf(stderr,"Error: incorrect command line format. Command %s not understood.\n",argv[i-1]);
			return 1;
		}
		if(argv[i-1][1]=='l'){
			if(gotl){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->sequence_number_limit=atoi(argv[i]);
			gotl=1;
		}else if(argv[i-1][1]=='d'){
			if(gotd){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->writeDatabase=1;
			sscanf(argv[i],"%s",opt->writeDatabaseName);
			gotd=1;
		}else if(argv[i-1][1]=='j'){
			if(gotj){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->justwriteDatabase=atoi(argv[i]);
			gotj=1;
		}else if(argv[i-1][1]=='t'){
			if(gott){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->writeDatabaseFull=atoi(argv[i]);
			gott=1;
		}else if(argv[i-1][1]=='f'){
			if(gotf){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->userAskedForFitting=atoi(argv[i]);
			gotf=1;
		}else if(argv[i-1][1]=='v'){
			if(gotv){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->verbose=atoi(argv[i]);
			gotv=1;
		}else if(argv[i-1][1]=='c'){
			if(gotc){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->useCancellation=atoi(argv[i]);
			gotc=1;
		}else if(argv[i-1][1]=='a'){
			if(gota){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->additionalDataWithSampling=atoi(argv[i]);
			gota=1;
		}else if(argv[i-1][1]=='h'){
			if(goth){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->histNcol=atoi(argv[i]);
			goth=1;
		}else if(argv[i-1][1]=='m'){
			if(gotm){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->discardOutside=atoi(argv[i]);
			gotm=1;
		}else if(argv[i-1][1]=='e'){
			if(gote){
				fprintf(stderr,"Error: argument %s given multiple times.\n",argv[i-1]);
				return 1;
			}
			opt->equilFraction=(float)atof(argv[i]);
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

int checkOptions(struct analysis_option_struct *opt, const struct script_struct *script){
	if(opt->additionalDataWithSampling>0&&opt->additionalDataWithSampling>script->Nadditional_data){
		fprintf(stderr,"Error: Your run ony has %d additional data columns, and you asked to use column %d with the -a option on the command line\n",script->Nadditional_data,opt->additionalDataWithSampling);
		return 1;
	}
	if(opt->histNcol<0){
		fprintf(stderr,"Error: number of customizable histogram bins must be >= 0 (0 is flag for using Nreplicas)\n");
		return 1;
	}
	if(opt->histNcol==0){
		opt->histNcol=script->Nreplicas;
	}
	if(opt->equilFraction<0.0||opt->equilFraction>1.0){
		fprintf(stderr,"Error: equilFraction must be on [0,1]\n");
		return 1;
	}
	if(opt->justwriteDatabase!=0 && opt->writeDatabase==0){
		fprintf(stderr,"Error: Can not *only* write the database when not writing the database at all.\n");
		return 1;
	}
	return 0;
}

int main(int argc, char *argv[]){
	unsigned int i,l;
	unsigned char current_page=1;
	int check;

	struct analysis_option_struct opt=DEFAULT_ANALYSIS_OPTION_STRUCT;
	struct script_struct script;
	struct pageinfo_struct pageinfo=FULLSIZE_PAGEINFO;
	struct stats_struct stats=EMPTY_STATS_STRUCT;
	struct database_struct db=EMPTY_DATABASE_STRUCT;

	float *cancel=(float *)NULL;
	struct exact_struct *exact=(exact_struct *)NULL;
	unsigned int **sequenceDensity=(unsigned int **)NULL;
	struct sampleDensity_struct sampleDensity=EMPTY_SAMPLEDENSITY;
	struct pmf_struct pmf=EMPTY_PMF;
	struct detailedBalance_struct *detailedBalance=(detailedBalance_struct *)NULL;
	struct nominal_struct *nominal=(struct nominal_struct *)NULL;
	struct graph_struct *graph=(struct graph_struct *)NULL;

	check=parseCommandLine(argc,argv,&opt);
	if(check!=0){
		showUsage(argv[0],&opt);
		exit(check);
	}

	read_input_script_file_class *input_script = new read_input_script_file_class;
	input_script->read_input_script_file(argv[1], &script);
	delete input_script;

	nominal=(struct nominal_struct *)malloc(MAX_REPLICAS*sizeof(nominal_struct));
	if(nominal==NULL){
		fprintf(stderr,"Error: unable to allocate memory for nominal\n");
		exit(1);
	}
	for(i=0;i<script.Nreplicas;i++){
		nominal[i].w[0]=script.replica[i].w_nominal;
		nominal[i].w[1]=script.replica[i].w2_nominal;
	}
	delete[] script.replica;

	check=checkOptions(&opt,&script);
	if(check!=0){
		showUsage(argv[0],&opt);
		exit(check);
	}
	sampleDensity.Ncol=pmf.Ncol=(unsigned int)opt.histNcol;
	fprintf(stderr,"Configurable histograms will have %d columns\n",opt.histNcol);

	check=read_in_database(&opt,&script,&db,nominal);
	if(check!=0){
		fprintf(stderr,"Error: read_in_database() returned non-zero\n");
		exit(check);
	}

	if(opt.justwriteDatabase!=0){
		//EARLY EXIT BECAUSE THE USER DOESN'T WANT THE POSTSCRIPT FILE
		exit(0);
	}

	if(db.record[db.Nrecords-1]->replica_number-db.replica_offset>script.Nreplicas){
		fprintf(stderr,"Error: there are more replicas in the database than there are in the script file\n");
		exit(1);
	}

	fprintf(stderr,"Reading in statistics\n");
	get_data_statistics(&stats,&db);
	fprintf(stderr,"Maximum force points found to be %u\n",stats.max_time);
	fprintf(stderr,"Maximum sequence number found to be %u\n",stats.max_sequence_number);
	if(stats.max_sequence_number>MAX_SEQUENCE_NUMBER){
		fprintf(stderr,"Error: too many sequence numbers, %hu is the limit\n", (unsigned short)MAX_SEQUENCE_NUMBER);
		exit(1);
	}
	
	printf("%%!PS-Adobe-2.0\n%%%%Created by program analyse_force_database\n\n");

	fprintf(stderr,"Reading data for force plots\n");
	graph=(struct graph_struct *)malloc(MAX_REPLICAS*sizeof(graph_struct));
	if(graph==NULL){
		fprintf(stderr,"Error: unable to allocate memory for graph\n");
		exit(1);
	}
	condense_forces(stats.max_time/N_FORCE_POINTS+1, nominal, graph, &script, &opt, &stats, &db);

	sequenceDensity=(unsigned int **)malloc((script.Nreplicas+1)*sizeof(unsigned int *));
	if(sequenceDensity==NULL){
		fprintf(stderr,"Error: Unable to allocate memory for sequenceDensity\n");
		exit(1);
	}
	for(i=0; i<=script.Nreplicas; i++){
		sequenceDensity[i]=(unsigned int *)malloc((script.Nreplicas+1)*sizeof(unsigned int));
		if(sequenceDensity[i]==NULL){
			fprintf(stderr,"Error: Unable to allocate memory for sequenceDensity[i]\n");
			exit(1);
		}
	}	
	calcSequenceDensity(sequenceDensity,&db,nominal,&script,&opt,&stats);

	if(opt.additionalDataWithSampling>0){
		sampleDensity.b=(unsigned int **)malloc((script.Nreplicas+1)*sizeof(unsigned int *));
		if(sampleDensity.b==NULL){
			fprintf(stderr,"Error: Unable to allocate memory for sampleDensity.b\n");
			exit(1);
		}
		for(i=0; i<=script.Nreplicas; i++){
			sampleDensity.b[i]=(unsigned int *)malloc((sampleDensity.Ncol+1)*sizeof(unsigned int));
			if(sampleDensity.b[i]==NULL){
				fprintf(stderr,"Error: Unable to allocate memory for sampleDensity.b[i]\n");
				exit(1);
			}
		}
		check=calcSampleDensity(&sampleDensity,opt.additionalDataWithSampling,graph,&db,nominal,&script,&opt,&stats);
		if(check!=0){
			fprintf(stderr,"Error: calcSampleDensity() returned non-zero\n");
			exit(check);
		}

		pmf.f=(float *)malloc((pmf.Ncol+1)*sizeof(float));
		if(pmf.f==NULL){
			fprintf(stderr,"Error: Unable to allocate memory for pmf.f\n");
			exit(1);
		}
		pmf.n=(unsigned int *)malloc((pmf.Ncol+1)*sizeof(unsigned int));
		if(pmf.n==NULL){
			fprintf(stderr,"Error: Unable to allocate memory for pmf.n\n");
			exit(1);
		}
		check=calcPMF(&pmf,opt.additionalDataWithSampling,graph,&db,nominal,&script,&opt);
		if(check!=0){
			fprintf(stderr,"Error: calcPMF() returned non-zero\n");
			exit(check);
		}
	}

	if(opt.useCancellation){
		if((cancel=(float *)malloc((script.Nreplicas+1)*sizeof(float)))==NULL){
			fprintf(stderr,"Error: Unable to allocate memory for cancellation terms, skipping\n");
			opt.useCancellation=0;
		}
		if(opt.useCancellation && getCancellationFromLog(opt.title,cancel,&script)!=0){
			fprintf(stderr,"Unable to determine cancellation from the log file, skipping\n");
			opt.useCancellation=0;
			free(cancel);
		}
	}

	exact=(exact_struct *)malloc(sizeof(exact_struct));
	if(exact==NULL){
		fprintf(stderr,"Error: Unable to allocate memory for extract\n");
		exit(1);
	}

	//Malloc for exact_struct in the getExactFromFile() function
	if(opt.useExact){
		if((exact=getExactFromFile(opt.title,exact,&script))==NULL){
			//No error message since this will only exist if this is from DR_tester run
			opt.useExact=0;
		}
	}

	detailedBalance=(struct detailedBalance_struct *)malloc(script.Nreplicas*sizeof(struct detailedBalance_struct));
	if(detailedBalance==NULL){
		fprintf(stderr,"Error: Unable to allocate memory for detailedBalance\n");
		exit(1);
	}

	int px,py;
	if(db.Nforces>0 && script.coordinate_type!=Temperature){
		// Combination page
		setPageinfoSmall(&pageinfo);
		setup_regular(&current_page);
		showFirstComboPageText(400,190);

		if((script.coordinate_type==Umbrella || script.coordinate_type==Temperature) && 
			script.Nadditional_data>0 && opt.additionalDataWithSampling>0)
		{
			print_postscript_sampleDensity(&pageinfo,&sampleDensity,&script,graph);
		}else{
			showSampleDensityNotAvailable(80,140,&script,&opt);
		}
		printf("0.0 180.0 translate\n");
		print_postscript_sequenceDensity(&pageinfo,sequenceDensity,&script,graph);
		printf("0.0 180.0 translate\n");
		print_postscript_pmf(0,&pageinfo,&script,exact,&pmf,graph);

		printf("400.0 0.0 translate\n");
		if(opt.useCancellation){
			print_postscript_dGandA(0,cancel,&pageinfo,&script,graph);
			printf("0.0 -180.0 translate\n");
			print_postscript_dGminusA(0,cancel,&pageinfo,&script,graph);
		}else{
			setFont(12);
			px=80;py=140;
			printf("%d %d moveto\n",px,py);
			printf("(Cancellation not available in log file) show\n");
			printf("0.0 -180.0 translate\n");
			printf("%d %d moveto\n",px,py);
			printf("(Cancellation not available in log file) show\n");
		}
		printf("showpage\n");
		setPageinfoFull(&pageinfo);

		if(script.coordinate_type!=Temperature){
			for(l=0;l<script.Nligands;l++){
				setup_regular(&current_page);
				print_postscript_pmf(l,&pageinfo,&script,exact,&pmf,graph);
				printf("showpage\n");
				setup_regular(&current_page);
				print_postscript_forceAverage(l,graph,&pageinfo,&script,&stats,&opt);
				printf("showpage\n");
			}
		}
	}
	setPageinfoFull(&pageinfo);

	if(db.Nforces>0&&script.Nadditional_data>0&&script.coordinate_type==Umbrella){
		setup_regular(&current_page);
		bool tempBool=pageinfo.showGrid;
		pageinfo.showGrid=false;
		print_postscript_samplingOverlap(opt.additionalDataWithSampling,&pageinfo,&script,&db,nominal,graph,&opt,&stats);
		pageinfo.showGrid=tempBool;
		printf("showpage\n");
	}

	setup_regular(&current_page);
	print_postscript_sequenceDensity(&pageinfo,sequenceDensity,&script,graph);
	printf("showpage\n");

	if((script.coordinate_type==Umbrella || script.coordinate_type==Temperature) &&
			script.Nadditional_data>0 && opt.additionalDataWithSampling>0)
	{
		setup_regular(&current_page);
		print_postscript_sampleDensity(&pageinfo,&sampleDensity,&script,graph);
		printf("showpage\n");
	}

	if(script.circular_replica_coordinate){
		setup_regular(&current_page);
		print_postscript_trajectory(0,-1,-1,&pageinfo,0,detailedBalance,&db,graph,&opt,&script,&stats,0);
		printf("showpage\n");
		setup_regular(&current_page);
		print_postscript_trajectory(0,1,0,&pageinfo,0,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,0);
		printf("showpage\n");
		if(script.Nadditional_data>0){
			setup_regular(&current_page);
			print_postscript_trajectory(0,-1,-1,&pageinfo,opt.additionalDataWithSampling,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,0);
			printf("showpage\n");
			setup_regular(&current_page);
			print_postscript_trajectory(0,1,0,&pageinfo,opt.additionalDataWithSampling,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,0);
			printf("showpage\n");
		}
	}else{
		setup_regular(&current_page);
		print_postscript_trajectory(0,0,-1,&pageinfo,0,detailedBalance,&db,graph,&opt,&script,&stats,0);
		printf("showpage\n");
		if(script.Nadditional_data>0){
			setup_regular(&current_page);
			print_postscript_trajectory(0,0,-1,&pageinfo,opt.additionalDataWithSampling,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,1);
			printf("showpage\n");

			for(int xq=0;xq<script.Nreplicas;xq++){
				setup_regular(&current_page);
				print_postscript_trajectory(xq*-1 -1,0,-1,&pageinfo,opt.additionalDataWithSampling,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,2);
				printf("showpage\n");
			}

			if(script.coordinate_type==Temperature){
				setup_regular(&current_page);
				print_postscript_trajectory(0,0,-1,&pageinfo,-opt.additionalDataWithSampling,(struct detailedBalance_struct *)NULL,&db,graph,&opt,&script,&stats,0);
				printf("showpage\n");
			}
		}
	}

	setPageinfoSmall(&pageinfo);
	setup_regular(&current_page);
	showSecondComboPageText(400,190);
	print_postscript_Pexchange(&pageinfo,detailedBalance,&script,graph);
	printf("0.0 180.0 translate\n");
	print_postscript_Pupdown(&pageinfo,detailedBalance,&script,graph);
	printf("showpage\n");

	if(db.Nforces>0){
		// Want this procedure always to be last
		if(opt.userAskedForFitting){
			fprintf(stderr,"Error: The fitting procedure has been removed.\n");
		}
	}

	free(db.unsorted_record);
	free(db.record);
	free(cancel);
}

/***********************************************************************************************/

void print_postscript_forceAverage(unsigned char ligand_number, const struct graph_struct *graph, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct stats_struct *stats, const struct analysis_option_struct *opt){
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;


	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(stats->max_force[ligand_number]-stats->min_force[ligand_number]);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(stats->min_force[ligand_number]*plotinfo.scale_y);
	plotinfo.yscalemult=(stats->max_force[ligand_number]-stats->min_force[ligand_number])/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
	plotinfo.ya=round_up(stats->min_force[ligand_number]/plotinfo.yscalemult);
	plotinfo.yb=round_down(stats->max_force[ligand_number]/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"Force Average");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{  
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"Average System Force\\(kcal/mol/unit\\)");
	plotinfo.ligand_number=ligand_number;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	printf("0.8 setgray\n");
	first=1;
	for(i=plotinfo.xa;i<plotinfo.xb;i++){
		w=graph[i].w[ligand_number];
		force=graph[i].average[ligand_number];
		
		if(graph[i].avg_weight_sum[0]>0.0001){
			printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,force*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
		}
	}
	printf("4 setlinewidth\n");
	printf("stroke\n");

	printf("0.0 setgray\n");
	for(i=plotinfo.xa;i<plotinfo.xb;i++){
		w=graph[i].w[ligand_number];
		printf("%% Start of new data at w = %f\n",w);
		printf("newpath\n");
		first=1;
		unsigned int max_force_points=0;
		for(j=0;j<N_FORCE_POINTS;j++)
			if(graph[i].weight_sum[j]>1.0e-5) max_force_points=j;

		for(j=(unsigned int)((float)max_force_points*opt->equilFraction);j<N_FORCE_POINTS;j++){
		//for(j=0;j<N_FORCE_POINTS;j++){
			force=graph[i].point[ligand_number][j];
			//printf("%% WEIGHT SUM: [i,j]=[%d %u]  %lf\n",i,j,graph[i].weight_sum[j]);
			if(graph[i].weight_sum[j]>=0.1){
				//fprintf(stderr,"FORCE: %u: %f\n",j,force);
		
				printf("%f %f ",((float)i*(N_FORCE_POINTS+1)+1+j)*plotinfo.scale_x+plotinfo.translate_x,force*plotinfo.scale_y+plotinfo.translate_y);
				if(first) printf("moveto\n");
				else printf("lineto\n");
				first=0;
			}else{
				first=1;
			}
		}
		printf("0.25 setlinewidth\n");
		printf("stroke\n");
	}
}

void print_postscript_trajectory(int printSelection,int circularFlag, int ignorerounds, const struct pageinfo_struct *pageinfo, int plotRCinstead, struct detailedBalance_struct *detbal, const struct database_struct *db, const struct graph_struct *graph, const struct analysis_option_struct *opt, const struct script_struct *script, const struct stats_struct *stats, int storecode){
	// Usage:
	// printSelection = 0 prints all
	// printSelection > 0 prints replica numbers divisible by printSelection (good for showing only a few)
	// printSelection < 0 prints only a single replica ( - printSelection - 1 ); the minus 1 is required to allow access to replica 0
	// circularFlag = 0 prints exactly at replicas
	// circularFlag > 0 uses periodic boundary conditions to expand the plot and avoid "jumps"
	// circularFlag < 0 uses knowledge of periodic boundary conditions to avoid the lines that span the whole graph during "jumps"
	// ignorerounds = 0 will display all rounds
	// ignorerounds > 0 will ignore the rounds from prior to the specified number 
	// ignorerounds < 0 will determine the number of rounds to ignore at runtime
	//
	// If whichData is nonzero, then the actual sampling will be plotted as opposed to the wref position

	static int dataflag=1; //ensures that this data is only printed the first time this routine is run
        static float storemax;
        static float storemin;
	int i,j,k,oldk,kmod,maxkmod,minkmod,didwrap;
	float oldkfrac,thiskfrac;
	unsigned char l;
	unsigned char first;
	float w;
	float space,fraction;
	char text[10];
	float max_f=-1000000.0,min_f=1000000.0;
	long *moveup,*movedown,*movesame;
	long *moveupdel,*movedowndel,*movesamedel;
	int autodetectignore,autodetectnotyet,didmove;
	float old_position;
	float *binlabels,rcmin,rcmax;
	int b,highEnd;
	float heatColor,heatR,heatG,heatB;
	float thex,they,old_thex,old_they;

	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;

	if(ignorerounds<0){
		ignorerounds=0;
		autodetectignore=1;
	}else{
		autodetectignore=0;
	}
	
	for(l=0;l<script->Nligands;l++){
		if(stats->max_force[l]>max_f) max_f=stats->max_force[l];
		if(stats->min_force[l]<min_f) min_f=stats->min_force[l];
	}

	//Determine how far outside of the boundaries any one replica went
	//Allow determination to be by w even if I am graphing additonal data
	maxkmod=minkmod=0;
	if(circularFlag>0){
		if(plotRCinstead){
			rcmin=rcmax=graph[0].rc_position[0];
		}

		for(i=0;i<script->Nreplicas;i++){
			if(printSelection<0&&i!=-printSelection-1){
				continue;
			}else if (printSelection>0 && div(i,printSelection).rem!=0){
				continue;
			}
			first=1;kmod=0;
			for(j=(unsigned int)((float)stats->max_sequence_number*opt->equilFraction);j<=stats->max_sequence_number;j++){
			//for(j=0;j<=stats->max_sequence_number;j++){
				if(plotRCinstead){
					w=graph[i].rc_position[j];
				}else{
					w=graph[i].replica_position[j];
				}
				if(w>1.0e10) continue;
				if(plotRCinstead){
					if(w<rcmin)rcmin=w;
					if(w>rcmax)rcmax=w;
				}

				if(script->coordinate_type==Temperature && plotRCinstead){
					for(k=0; k<opt->histNcol; ++k){
						if(w<binlabels[b])break;
					}
					if(k==0) k++;
					if(k==opt->histNcol) k--;
					space=binlabels[k]-binlabels[k-1];
					fraction=(w-binlabels[k-1])/space;
					thiskfrac=(float)k+fraction;
					highEnd=opt->histNcol;
				}else{
					for(k=0;k<script->Nreplicas;k++){
						if(w<graph[k].w[0]) break;
					}
					if(k==0) k++;
					if(k==script->Nreplicas) k--;
					space=graph[k].w[0]-graph[k-1].w[0];
					fraction=(w-graph[k-1].w[0])/space;
					thiskfrac=(float)k+fraction;
					highEnd=script->Nreplicas;
				}

				if(!first){
					if(thiskfrac-oldkfrac>(float)highEnd/2.0){
						if(j>=ignorerounds)kmod--;
					}else if(oldkfrac-thiskfrac>(float)highEnd/2.0){
						if(j>=ignorerounds)kmod++;
					}
				}
				oldkfrac=(float)k+fraction;
				oldk=(int)round(oldkfrac);
				k+=kmod;
				first=0;
				old_position=(float)k+fraction+kmod*(highEnd-1);
			}
			if(kmod>maxkmod)maxkmod=kmod;
			if(kmod<minkmod)minkmod=kmod;
		}
		//Make it symmetric will assist printout
		//if(-minkmod>maxkmod)maxkmod=-minkmod;
		//if(-maxkmod<minkmod)minkmod=-maxkmod;
	}else if(plotRCinstead){
		rcmin=rcmax=graph[0].rc_position[0];
		for(i=0;i<script->Nreplicas;i++){
			if(printSelection<0&&i!=-printSelection-1){
				continue;
			}else if(printSelection>0 && div(i,printSelection).rem!=0){
				continue;
			}
			for(j=0;j<=stats->max_sequence_number;j++){
				w=graph[i].rc_position[j];
				if(w>1.0e10) continue;
				if(w<rcmin)rcmin=w;
				if(w>rcmax)rcmax=w;
			}
		}
	}

        if(storecode==2){
          //revert to stored values
          rcmax=storemax;
          rcmin=storemin;
        }
        if(storecode==1){
          //store these values
          storemax=rcmax;
          storemin=rcmin;
        }

	//if(script->coordinate_type==Temperature && plotRCinstead){
	if(plotRCinstead){
		plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(opt->histNcol*(N_FORCE_POINTS+1));
		plotinfo.translate_x=(1.0*72);
		plotinfo.xa=0;
		plotinfo.xb=opt->histNcol;
		plotinfo.xc=(int)(ceil)((float)(opt->histNcol)/(float)script->Nreplicas);
		binlabels=(float *)calloc((opt->histNcol),sizeof(float));
		if(binlabels==NULL){
			fprintf(stderr,"Error: Memory allocation error for binlabels\n");
			exit(1);
		}
		for(b=0; b<opt->histNcol; ++b){
			binlabels[b]=rcmin+((float)b+0.5)*((rcmax-rcmin)/(float)opt->histNcol);
		}
		plotinfo.xlabels=binlabels;
	}else{
		plotinfo.scale_x=((float)pageinfo->width-1.0-0.5)*72/(float)(script->Nreplicas*(maxkmod-minkmod+1)*(N_FORCE_POINTS+1));
		if(circularFlag>0&&maxkmod-minkmod>0){
			plotinfo.translate_x=(1.0*72)*((float)(maxkmod-minkmod)+0.5);
		}else{
			plotinfo.translate_x=(1.0*72)*(float)(-minkmod+1);
		}
		plotinfo.xa=0;
		plotinfo.xb=script->Nreplicas;
		plotinfo.xc=1;
		plotinfo.xlabels=(float *)NULL;
	}
	plotinfo.scale_y=(((float)pageinfo->height-1.0-0.5-0.1-0.1)*72)/(float)(stats->max_sequence_number-ignorerounds);
	plotinfo.translate_y=(1.0*72)-((float)ignorerounds*plotinfo.scale_y);
	plotinfo.yscalemult=(float)(stats->max_sequence_number-ignorerounds)/(float)NUMVALUES;
	plotinfo.ya=round_up((0+ignorerounds)/plotinfo.yscalemult);
	plotinfo.yb=round_down(stats->max_sequence_number/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"Trajectory movement");
	if(plotRCinstead){
		if(script->coordinate_type==Temperature && plotRCinstead<0){
			sprintf(plotinfo.xtitle,"Reaction Coordinate (sampled) -- color by temperature");
		}else{
			sprintf(plotinfo.xtitle,"Reaction Coordinate (sampled)");
		}
	}else{
		if(script->coordinate_type==Temperature){
			sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
		}else{
			sprintf(plotinfo.xtitle,"Reaction Coordinate (restraint)");
		}
	}
	sprintf(plotinfo.ytitle,"Sequence Number");
	plotinfo.ligand_number=0;

	showPlot(&plotinfo,pageinfo,script,graph);
			 
	if((moveup=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for moveup in print_postscript_page4()\n");
	}
	if((movedown=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for movedown in print_postscript_page4()\n");
	}
	if((movesame=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for movesame in print_postscript_page4()\n");
	}
	if((moveupdel=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for moveupdel in print_postscript_page4()\n");
	}
	if((movedowndel=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for movedowndel in print_postscript_page4()\n");
	}
	if((movesamedel=(long *)calloc(script->Nreplicas+1,sizeof(long)))==NULL){
		fprintf(stderr,"Memory allocation error for movesamedel in print_postscript_page4()\n");
	}

	printf("%% Start of replica positions plots\n");
	double productivity_ratio=0.0;
	unsigned int productivity_ratio_count=0;
	if(script->coordinate_type==Temperature && plotRCinstead<0){
		printf(DASHONSTRONG);
	}

	FILE *g;
	int mycolors[1000];
	int myi=-1;
	g=fopen("colors.txt","r");
	if(g!=NULL){
		int myj;
		for(myj=0;myj<1000;myj++){
			mycolors[myj]=-1;
		}
		char linein[1001],firststring[1001];
		while(fgets(linein,1000,g)!=NULL && sscanf(linein,"%s",firststring)!=0){
			++myi;
			sscanf(linein,"%d",&(mycolors[myi]));
		}
		fclose(g);
	}

	for(i=0;i<script->Nreplicas;i++){
		if(autodetectignore==1){
			autodetectnotyet=1;
		}else{
			autodetectnotyet=0;
		}

		if(printSelection<0&&i!=-printSelection-1)continue;
		else if(printSelection>0 && div(i,printSelection).rem!=0)continue;
		if(myi>=0 && mycolors[i]>=0){
			printf("%s setrgbcolor\n",get_RGB_colour(mycolors[i]));
		}else{
			printf("%s setrgbcolor\n",get_RGB_colour(i));
		}
		printf("newpath\n");

		first=1;kmod=0;
		for(j=(unsigned int)((float)stats->max_sequence_number*opt->equilFraction);j<=stats->max_sequence_number;j++){
		//for(j=0;j<=stats->max_sequence_number;j++){
			if(plotRCinstead){
				w=graph[i].rc_position[j];
			}else{
				w=graph[i].replica_position[j];
			}
			//fprintf(stderr,"REPLICA W: %f\n",w);
			if(w>1.0e10) continue;
			//if(script->coordinate_type==Temperature && plotRCinstead){
			if(plotRCinstead){
				for(k=0; k<opt->histNcol; ++k){
					if(w<binlabels[b])break;
				}
				if(k==0) k++;
				if(k==opt->histNcol) k--;
				space=binlabels[k]-binlabels[k-1];
				fraction=(w-binlabels[k-1])/space;
				thiskfrac=(float)k+fraction;
				highEnd=opt->histNcol;
			}else{
				for(k=0;k<script->Nreplicas;k++){
					if(w<graph[k].w[0]) break;
				}
				if(k==0) k++;
				if(k==script->Nreplicas) k--;
				space=graph[k].w[0]-graph[k-1].w[0];
				fraction=(w-graph[k-1].w[0])/space;
				thiskfrac=(float)k+fraction;
				highEnd=script->Nreplicas;
			}

			didwrap=0;
			didmove=0;
			if(!first){
				if(thiskfrac-oldkfrac>(float)highEnd/2.0){
					//Just wrapped from low end to high end therefore moved down
					kmod-=highEnd-1;
					if(j>=ignorerounds) movedown[oldk]++;
					didwrap=-1;
					if(autodetectignore&&autodetectnotyet){
						movedowndel[oldk]++;
						didmove=1;
					}
				}else if(oldkfrac-thiskfrac>(float)highEnd/2.0){
					//Just wrapped from high end to low end therefore moved up
					kmod+=highEnd-1;
					if(j>=ignorerounds) moveup[oldk]++;
					didwrap=1;
					if(autodetectignore&&autodetectnotyet){
						moveupdel[oldk]++;
						didmove=1;
					}
				}else if(thiskfrac>oldkfrac){
					if(j>=ignorerounds) moveup[oldk]++;
					if(autodetectignore&&autodetectnotyet){
						moveupdel[oldk]++;
						didmove=1;
					}
				}else if(thiskfrac<oldkfrac){
					if(j>=ignorerounds) movedown[oldk]++;
					if(autodetectignore&&autodetectnotyet){
						movedowndel[oldk]++;
						didmove=1;
					}
				}else{
					if(j>=ignorerounds) movesame[oldk]++;
					if(autodetectignore&&autodetectnotyet){
						movesamedel[oldk]++;
					}
				}
			}
			oldkfrac=(float)k+fraction;
			oldk=(int)round(oldkfrac);
			if(circularFlag>0)k+=kmod;

			if(autodetectignore&&autodetectnotyet&&didmove){
				autodetectnotyet=0;
			}

			if(j>=ignorerounds){
				if(script->coordinate_type==Temperature && plotRCinstead<0 && !first){
					printf("%f %f moveto\n",old_thex,old_they);
				}
				thex=old_thex=((float)k-0.5+fraction)*(N_FORCE_POINTS+1.0)*plotinfo.scale_x+plotinfo.translate_x;
				they=old_they=(float)j*plotinfo.scale_y+plotinfo.translate_y;
				printf("%f %f ",thex,they);
				if(first){
					 printf("moveto\n");
				}else{
					if((didwrap!=0 && circularFlag<0)||j==ignorerounds) printf("moveto\n");
					else printf("lineto\n");
					if(autodetectignore==0||autodetectnotyet==0){
						if(circularFlag<0){
						//CN is not sure if (Nreplicas-1) is sufficient for boltzmann jumping or continuous
							productivity_ratio+=fabs((float)k+fraction-old_position+(float)(didwrap*(script->Nreplicas-1)));
						}else{
							//non-circuar runs and expanded plots (already consider pbc jumps)
							productivity_ratio+=fabs((float)k+fraction-old_position);	
						}
						productivity_ratio_count++;
					}
				}
				if(script->coordinate_type==Temperature && plotRCinstead<0 && !first){
					//plotRCinstead<0 is a flag to color by temperature
					//increasing [i] is decreasing temperature
					//increasing heatColor is increasing Temperature
					//increasing heatColor is increasing Temperature
					heatColor=(graph[i].replica_position[j]-graph[script->Nreplicas-1].w[0])/(graph[0].w[0]-graph[script->Nreplicas-1].w[0]);
					heatR=heatColor-0.5;
					if(heatR<0.0)heatR=0.0;
					heatR*=2.0;
					heatB=1.0-heatColor-0.5;
					if(heatB<0.0)heatB=0.0;
					heatB*=2.0;
					heatG=1.0-heatR-heatB;
					printf("%f %f %f setrgbcolor\n",heatR,heatG,heatB);
					//printf("%% heatcolor %f (%f,%f,%f) for a temperature of %f (%f)\n",heatColor,heatR,heatG,heatB,graph[i].replica_position[j],4184.0*(1/graph[i].replica_position[j])/8.31451);
					printf("stroke\n");
					printf("newpath\n");
				}
			}
			first=0;
			old_position=(float)k+fraction;
		}
		printf("0.5 setlinewidth\n");
		printf("stroke\n");
	}
	if(script->coordinate_type==Temperature && plotRCinstead<0){
		printf(DASHOFF);
	}

	if(productivity_ratio_count>0) productivity_ratio/=productivity_ratio_count;
	if(dataflag)printf("%% PRODUCTIVITY RATIO: %lf   N_MOVES: %u\n",productivity_ratio,productivity_ratio_count);
	double numtotattempt;
	for(i=0;i<script->Nreplicas;i++){
		if(printSelection<0&&i!=-printSelection-1)continue;
		else if(printSelection>0 && div(i,printSelection).rem!=0)continue;
		numtotattempt=(double)(moveup[i+1]+movedown[i+1]+movesame[i+1]-moveupdel[i+1]-movedowndel[i+1]-movesamedel[i+1]);
		if(dataflag){
			printf("%% REPLICA W: %f \tMOVEUP: %f \tMOVEDOWN: %f \tMOVESAME: %f \tNUMSTEPS: %ld\n",graph[i].w[0],(double)(moveup[i+1]-moveupdel[i+1])/numtotattempt,(double)(movedown[i+1]-movedowndel[i+1])/numtotattempt,(double)(movesame[i+1]-movesamedel[i+1])/numtotattempt,(long)numtotattempt);
			if(detbal!=NULL){
				detbal[i].w=graph[i].w[0];
				detbal[i].moveUp=(double)(moveup[i+1]-moveupdel[i+1])/numtotattempt;
				detbal[i].moveDown=(double)(movedown[i+1]-movedowndel[i+1])/numtotattempt;
				detbal[i].moveSame=(double)(movesame[i+1]-movesamedel[i+1])/numtotattempt;
				detbal[i].numSteps=(long)numtotattempt;
			}
		}
	}
	dataflag=0;
}

void print_postscript_sequenceDensity(const struct pageinfo_struct *pageinfo, unsigned int * const *sequenceDensity, const struct script_struct *script, const struct graph_struct *graph){
	int i,j;
	unsigned char first;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/1.0;
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(0.0*plotinfo.scale_y);
	plotinfo.yscalemult=1.0/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
	plotinfo.ya=0;
	plotinfo.yb=round_down(1.0/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"Relative Intended Sample Density");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"Relative Intended Sample Density");
	plotinfo.ligand_number=0;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	for(j=0;j<=script->Nreplicas;j++){	
		first=1;
		for(i=0;i<script->Nreplicas;i++){
			printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(float)sequenceDensity[j][i]/(float)sequenceDensity[j][script->Nreplicas]*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
		}
		if(j==script->Nreplicas){
			printf("4 setlinewidth\n");
			printf("0 setgray\n");
			printf("stroke\n");

		}else{
			printf("0.5 setlinewidth\n");
			printf("%s setrgbcolor\n",get_RGB_colour(j));
			printf("stroke\n");
		}
	}
}


void print_postscript_sampleDensity(const struct pageinfo_struct *pageinfo, const struct sampleDensity_struct *sd, const struct script_struct *script, const struct graph_struct *graph){
	int i,j;
	unsigned char first;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;
	float absmax,absmin;
	float *binlabels;

	absmax=1.0;
	absmin=0.0;

	binlabels=(float *)calloc((sd->Ncol),sizeof(float));
	if(binlabels==NULL){
		fprintf(stderr,"Error: Memory allocation error for binlabels\n");
		exit(1);
	}
	for(j=0; j<sd->Ncol; ++j){
		binlabels[j]=sd->min+((float)j+0.5)*sd->binWidth;
	}

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(sd->Ncol*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(absmax-absmin);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(absmin*plotinfo.scale_y);
	plotinfo.yscalemult=(absmax-absmin)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=(int)sd->Ncol;
	plotinfo.xc=(int)(ceil)((float)(sd->Ncol)/(float)script->Nreplicas);
	plotinfo.ya=round_up(absmin/plotinfo.yscalemult);
	plotinfo.yb=round_down(absmax/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"Relative Actual Sample Density");
	//This is plotted vs. reaction coordinate even for temperature simulations
	sprintf(plotinfo.xtitle,"Reaction Coordinate");
	sprintf(plotinfo.ytitle,"Relative Actual Sample Density");
	plotinfo.ligand_number=0;
	plotinfo.xlabels=binlabels;

	showPlot(&plotinfo,pageinfo,script,graph);

	for(i=0;i<=script->Nreplicas;i++){	
		first=1;
		for(j=0;j<(int)sd->Ncol;j++){
			printf("%f %f ",((((float)j)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(float)sd->b[i][j]/(float)sd->b[i][sd->Ncol]*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
		}
		if(i==script->Nreplicas){
			printf("4 setlinewidth\n");
			printf("0 setgray\n");
			printf("stroke\n");

		}else{
			printf("0.5 setlinewidth\n");
			printf("%s setrgbcolor\n",get_RGB_colour(i));
			printf("stroke\n");
		}
	}
}

void print_postscript_pmf(unsigned char ligand_number, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const exact_struct *exact, const struct pmf_struct *pmf, const struct graph_struct *graph){
	static int alreadyran=0;
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force,sumf,minsumf,maxsumf,offset;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;
	float exactMin,exactMax;
	float nsumf,nminsumf,nmaxsumf,noffset;
	float *binlabels;

	minsumf=maxsumf=offset=0.0;

	if(script->coordinate_type==Spatial||script->coordinate_type==Umbrella){
		sumf=0.0;
		for(i=1;i<script->Nreplicas;i++){
			sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
			if(sumf>maxsumf)maxsumf=sumf;
			if(sumf<minsumf)minsumf=sumf;
		}
		offset=minsumf;
		//fprintf(stderr,"PMF graph yrange = [%f to %f] with offset %f\n",minsumf,maxsumf,offset);
	}

	if(exact!=NULL){
		exactMin=exactMax=exact->e[0];
		for(i=1;i<exact->n;i++){
			if(exact->e[i]<exactMin)exactMin=exact->e[i];
			if(exact->e[i]>exactMax)exactMax=exact->e[i];
		}
		if(exactMax>maxsumf)maxsumf=exactMax;
		//fprintf(stderr,"PMF graph yrange = [%f to %f] with offset %f\n",minsumf,maxsumf,offset);
		//if(exactMax+offset>maxsumf)maxsumf=exactMax+offset;
	}

	//min of both is zero, so no need to correct minsumf	
	if(pmf->f!=NULL){
		nsumf=nminsumf=nmaxsumf=noffset=0.0;
		for(i=1;i<pmf->Ncol;i++){
			nsumf-=(pmf->f[i]+pmf->f[i-1])/2*pmf->binWidth;
			if(nsumf>nmaxsumf)nmaxsumf=nsumf;
			if(nsumf<nminsumf)nminsumf=nsumf;
			//fprintf(stderr,"Temporary PMF graph yrange = [%f to %f]\n",nminsumf,nmaxsumf); 
		}
		noffset=nminsumf;
		if(script->coordinate_type==Temperature){
			//There is no original graph for Temperature
			offset=noffset;
			minsumf=nminsumf;
		}
	}

	fprintf(stderr,"PMF graph yrange = [%f to %f] with offset %f\n",minsumf,maxsumf,offset);
	if(script->coordinate_type==Temperature){
		plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(pmf->Ncol*(N_FORCE_POINTS+1));
		plotinfo.translate_x=(1.0*72);
		plotinfo.xa=0;
		plotinfo.xb=pmf->Ncol;
		plotinfo.xc=(int)(ceil)((float)(pmf->Ncol)/(float)script->Nreplicas);
		binlabels=(float *)calloc((pmf->Ncol),sizeof(float));
		if(binlabels==NULL){
			fprintf(stderr,"Error: Memory allocation error for binlabels\n");
			exit(1);
		}
		for(int b=0; b<pmf->Ncol; ++b){
			binlabels[b]=pmf->min+((float)b+0.5)*((pmf->max-pmf->min)/(float)pmf->Ncol);
		}
		plotinfo.xlabels=binlabels;
	}else{
		plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
		plotinfo.translate_x=(1.0*72);
		plotinfo.xa=0;
		plotinfo.xb=script->Nreplicas;
		plotinfo.xc=1;
		plotinfo.xlabels=(float *)NULL;
	}
	//CN notes that "(maxsumf-(minsumf-offset)" can equal zero and then there is a div by zero nan error
	//This can be fixed with the line below, although that causes lots of other problems:
	//offset=0;
	// note that the above line is not an optimal solution
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(maxsumf-(minsumf-offset));
	plotinfo.translate_y=(1.1*72)-(0.0*plotinfo.scale_y);
	plotinfo.yscalemult=(maxsumf-(minsumf-offset))/(float)pageinfo->numvalues;
	plotinfo.ya=round_up((minsumf-offset)/plotinfo.yscalemult);
	plotinfo.yb=round_down((maxsumf)/plotinfo.yscalemult); //CN: used to be (maxsumf-offset)
	plotinfo.yc=1;

	//TRY THIS:
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(maxsumf-minsumf);
	//The following line was resulting in a pmf that was too high up
	//plotinfo.translate_y=(1.1*72)-(offset*plotinfo.scale_y);
	plotinfo.translate_y=(1.1*72)-(0.0*plotinfo.scale_y);

	plotinfo.yscalemult=(maxsumf-minsumf)/(float)pageinfo->numvalues;
	plotinfo.ya=round_up((minsumf-offset)/plotinfo.yscalemult);
	plotinfo.yb=round_down((maxsumf-offset)/plotinfo.yscalemult);
	plotinfo.yc=1;

	sprintf(plotinfo.title,"Force Integral");
	//CN doesn't need this any more: sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	sprintf(plotinfo.xtitle,"Reaction Coordinate");
	sprintf(plotinfo.ytitle,"Force Integral \\(kcal/mol\\)");
	plotinfo.ligand_number=ligand_number;

	//fprintf(stderr,"maxsumf=%f\tminsumf=%f\toffset=%f\tscale_y=%f\ttranslate_y=%f\tyscalemult=%f\tya=%f\tyb=%f\tyc=%f\n",maxsumf,minsumf,offset,plotinfo.scale_y,plotinfo.translate_y,plotinfo.yscalemult,plotinfo.ya,plotinfo.yb,plotinfo.yc); 

	showPlot(&plotinfo,pageinfo,script,graph);

	if(script->coordinate_type==Spatial||script->coordinate_type==Umbrella){
		first=1;
		for(i=0;i<script->Nreplicas;i++){
			w=graph[i].w[ligand_number];
			force=graph[i].average[ligand_number];
			if(i==0){
				sumf=minsumf=maxsumf=0.0;
			}else{
				sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
			}
			if(graph[i].avg_weight_sum[0]>0.0001)
			{
				printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(sumf-offset)*plotinfo.scale_y+plotinfo.translate_y);
				if(first) printf("moveto\n");
				else printf("lineto\n");
				if(!alreadyran)printf("%% PMF: %f %f\n",graph[i].w[ligand_number],sumf-offset);
				first=0;
			}
		}
		printf("1 setlinewidth\n");
		printf("stroke\n");
	}
	//if(script->coordinate_type==Spatial||script->coordinate_type==Umbrella){
		if(exact!=NULL){
			printf(DASHON);
			first=1;
			for(i=0;i<exact->n;i++){
				printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(exact->e[i]-exactMin)*plotinfo.scale_y+plotinfo.translate_y);
				if(first) printf("moveto\n");
				else printf("lineto\n");
				first=0;
			}
			printf("1 setlinewidth\n");
			printf("1 0 0 setrgbcolor\n");
			printf("stroke\n");
			printf(DASHOFF);
		}
	//}

	if(pmf->f!=NULL){
		first=1;
		for(i=0;i<(int)pmf->Ncol;i++){
			if(i==0){
				nsumf=nminsumf=nmaxsumf=0.0;
			}else{
				nsumf-=(pmf->f[i]+pmf->f[i-1])/2*pmf->binWidth;
			}
			printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x*script->Nreplicas/(float)pmf->Ncol+plotinfo.translate_x,(nsumf-noffset)*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
			if(!alreadyran)printf("%% PMF_rehisto: %f %f %d\n",pmf->min+pmf->binWidth*(float)i,nsumf-noffset,pmf->n[i]);
		}
		printf("1 setlinewidth\n");
		printf("0 0 1 setrgbcolor\n");
		printf("stroke\n");
	}

	printf("0 0 0 setrgbcolor\n");
	alreadyran=1;
}

void print_postscript_dGandA(unsigned char ligand_number, const float *cancel, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct graph_struct *graph){
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force,sumf,minsumf,maxsumf,offset;
	float maxcancel,mincancel,absmin,absmax;
	int pointAtMin;

	sumf=minsumf=maxsumf=0.0;
	mincancel=maxcancel=cancel[0];
	pointAtMin=0;
	for(i=1;i<script->Nreplicas;i++){
		sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
		if(sumf>maxsumf)maxsumf=sumf;
		if(sumf<minsumf){
			minsumf=sumf;
			pointAtMin=i;
		}
		if(cancel[i]>maxcancel)maxcancel=cancel[i];
		if(cancel[i]<mincancel)mincancel=cancel[i];
	}
	offset=minsumf;
	absmax=(maxcancel>maxsumf-offset)?maxcancel:maxsumf-offset;
	absmin=(mincancel<minsumf-offset)?mincancel:minsumf-offset;

	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(absmax-absmin);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(absmin*plotinfo.scale_y);
	plotinfo.yscalemult=(absmax-absmin)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
	plotinfo.ya=round_up(absmin/plotinfo.yscalemult);
	plotinfo.yb=round_down(absmax/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"Force Integral");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"Force Integral \\(kcal/mol\\)");
	plotinfo.ligand_number=ligand_number;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	first=1;
	for(i=0;i<script->Nreplicas;i++){
		w=graph[i].w[ligand_number];
		force=graph[i].average[ligand_number];
		if(i==0){
			sumf=minsumf=maxsumf=0.0;
		}else{
			sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
		}
		if(graph[i].avg_weight_sum[0]>0.0001){
			printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(sumf-offset)*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
		}
	}
	printf("4 setlinewidth\n");
	printf("stroke\n");

	first=1;
	for(i=0;i<script->Nreplicas;i++){
		printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,cancel[i]*plotinfo.scale_y+plotinfo.translate_y);
		if(first) printf("moveto\n");
		else printf("lineto\n");
		first=0;
	}
	printf("1 setlinewidth\n");
	printf("stroke\n");

	first=1;
	for(i=0;i<script->Nreplicas;i++){
		printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(cancel[i]*-1.0+cancel[pointAtMin])*plotinfo.scale_y+plotinfo.translate_y);
		if(first) printf("moveto\n");
		else printf("lineto\n");
		first=0;
	}
	printf("1 setlinewidth\n");
	printf("1 0 0 setrgbcolor\n");
	printf("stroke\n");

}

void print_postscript_dGminusA(unsigned char ligand_number, const float *cancel, const struct pageinfo_struct *pageinfo, const struct script_struct *script, const struct graph_struct *graph){
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force,sumf,offset;
	float absmin,absmax;
	float val;
	float diffs,diffc,oldsumf,diffspcmin,diffspcmax;
	int r;

	sumf=oldsumf=0.0;
	absmax=absmin=sumf+cancel[0];
	diffspcmin=diffspcmax=0.0;
	for(i=1;i<script->Nreplicas;i++){
		sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
		diffs=sumf-oldsumf;
		diffc=cancel[i]-cancel[i-1];
		if(diffs+diffc<diffspcmin)diffspcmin=diffs+diffc;
		if(diffs+diffc>diffspcmax)diffspcmax=diffs+diffc;
		val=sumf+cancel[i];
		if(val>absmax)absmax=val;
		if(val<absmin)absmin=val;
		oldsumf=sumf;
	}
	offset=absmin;
	absmax-=offset;
	absmin-=offset;
	if(diffspcmin<absmin)absmin=diffspcmin;
	if(diffspcmax>absmax)absmax=diffspcmax;

	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(absmax-absmin);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(absmin*plotinfo.scale_y);
	plotinfo.yscalemult=(absmax-absmin)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
		plotinfo.ya=round_up(absmin/plotinfo.yscalemult);
		plotinfo.yb=round_down(absmax/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"dG-A");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"dG-A");
	plotinfo.ligand_number=ligand_number;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	for(r=1; r<=2; r++){
		first=1;
		for(i=0;i<script->Nreplicas;i++){
			w=graph[i].w[ligand_number];
			force=graph[i].average[ligand_number];
			if(i==0){
				sumf=oldsumf=0.0;
			}else{
				sumf-=(graph[i].average[ligand_number]+graph[i-1].average[ligand_number])/2*(graph[i].w[ligand_number]-graph[i-1].w[ligand_number]);
			}
			diffs=sumf-oldsumf;
			diffc=cancel[i]-cancel[i-1];
			if(graph[i].avg_weight_sum[0]>0.0001){
				if(r==1){
					printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,((sumf-offset)+cancel[i])*plotinfo.scale_y+plotinfo.translate_y);
				}else{
					printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(diffs+diffc)*plotinfo.scale_y+plotinfo.translate_y);
				}
				if(first) printf("moveto\n");
				else printf("lineto\n");
				first=0;
			}
			oldsumf=sumf;
		}
		if(r==1){
			printf("4 setlinewidth\n");
		}else{
			printf("1 setlinewidth\n");
		}
		printf("stroke\n");
	}
}

#define numBins 100
void print_postscript_samplingOverlap(int whichData, const struct pageinfo_struct *pageinfo, const struct script_struct *script, struct database_struct *db, const struct nominal_struct *nominal, const struct graph_struct *graph, const struct analysis_option_struct *opt, const struct stats_struct *stats){
	//use whichData==1 for the first additionalData
	int i,j;
	int **histo;
	bool first=true;
	double min,max;
	int *tot;
	float scale_x,scale_y;
	float translate_x,translate_y;
	float w;
	float force,sumf,minsumf,maxsumf,offset;
	char text[10];
	float yscalemult;
	float absmin,absmax;
	float binWidth;
	float *binlabels;

	//Format:
	//force[0]=record[recordN]->generic_data[sampleN*Nligands]

	first=true;
	for(i=0;i<db->Nrecords;i++){
		for(j=script->Nsamples_per_run*(whichData); j<script->Nsamples_per_run*(whichData+1); j++){
			if(first){
				min=max=db->record[i]->generic_data[j];
				first=false;
			}		
			if(db->record[i]->generic_data[j]>max)max=db->record[i]->generic_data[j];
			if(db->record[i]->generic_data[j]<min)min=db->record[i]->generic_data[j];
		}
	}
	binWidth=(max-min)/(float)numBins;
	//fprintf(stderr,"DATA_INFO: max:%f min:%f binWidth:%f\n",max,min,binWidth);

	histo=(int **)malloc((script->Nreplicas)*sizeof(int *));
	for(i=0; i<script->Nreplicas; i++){
		histo[i]=(int *)calloc((numBins+1),sizeof(int));
		//require +1 to handle max value, it will later be dumped down
		if(histo[i]==NULL){
			fprintf(stderr,"Error: Memory allocation error for histo\n");
			exit(1);
		}
	}
	tot=(int *)calloc((script->Nreplicas),sizeof(int));
	if(tot==NULL){
		fprintf(stderr,"Error: Memory allocation error for tot\n");
		exit(1);
	}
	binlabels=(float *)calloc((numBins),sizeof(float));
	if(binlabels==NULL){
		fprintf(stderr,"Error: Memory allocation error for binlabels\n");
		exit(1);
	}
	for(j=0; j<numBins; ++j){
		binlabels[j]=min+((float)j+0.5)*binWidth;
		//printf("%% Label[i]=%f\n",binlabels[i]);
	}

	minsumf=0.0;
	for(i=0;i<db->Nrecords;i++){
		for(j=script->Nsamples_per_run*(whichData); j<script->Nsamples_per_run*(whichData+1); j++){
			if(db->record[i]->sequence_number<(unsigned int)(opt->equilFraction*stats->max_sequence_number)){
				continue;
			}
			//fprintf(stderr,"Putting w:%f val:%f in histo[%d][%d]\n",record[i]->w,record[i]->generic_data[j],find_bin_from_w(record[i]->w),(int)floor((record[i]->generic_data[j]-min)/binWidth));
			++histo[find_bin_from_w(db->record[i]->w,0,nominal,script)][(int)floor((db->record[i]->generic_data[j]-min)/binWidth)];
			++tot[find_bin_from_w(db->record[i]->w,0,nominal,script)];
		}
	}
	for(i=0;i<script->Nreplicas; i++){
		histo[i][numBins-1]+=histo[i][numBins];
		histo[i][numBins]=0;
	}

	absmin=absmax=0.0;
	for(i=0; i<script->Nreplicas; i++){
		for(j=0; j<numBins; ++j){
			if((float)histo[i][j]/(float)tot[i]>absmax)absmax=(float)histo[i][j]/(float)tot[i];
		}
	}
	offset=0.0;

	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(numBins*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(absmax-absmin);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(absmin*plotinfo.scale_y);
	plotinfo.yscalemult=(absmax-absmin)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=(int)numBins;
	plotinfo.xc=(int)(ceil)((float)(numBins)/(float)script->Nreplicas);
	plotinfo.ya=round_up(absmin/plotinfo.yscalemult);
	plotinfo.yb=round_down(absmax/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"P\\(rc\\) one line per umbrella center");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"P\\(rc\\) one line per umbrella center");
	plotinfo.ligand_number=0;
	plotinfo.xlabels=binlabels;

	showPlot(&plotinfo,pageinfo,script,graph);

	for(i=0; i<script->Nreplicas; i++){
		first=true;
		for(j=0; j<numBins; ++j){
			printf("%f %f ",((((float)j)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(float)histo[i][j]/(float)tot[i]*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=false;
		}
		printf("1 setlinewidth\n");
		printf("%s setrgbcolor\n",get_RGB_colour(i));
		printf("stroke\n");
	}
}

void print_postscript_Pexchange(const struct pageinfo_struct *pageinfo, struct detailedBalance_struct *detbal, const struct script_struct *script, const struct graph_struct *graph){
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force,minsumf,maxsumf,offset;
	double val;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;
	int which;

	minsumf=maxsumf=detbal[0].moveSame;
	for(i=0;i<script->Nreplicas;i++){
		if(detbal[i].moveSame>maxsumf)maxsumf=detbal[i].moveSame;
		if(detbal[i].moveUp>maxsumf)maxsumf=detbal[i].moveUp;
		if(detbal[i].moveDown>maxsumf)maxsumf=detbal[i].moveDown;
		if(detbal[i].moveSame<minsumf)minsumf=detbal[i].moveSame;
		if(detbal[i].moveUp<minsumf)minsumf=detbal[i].moveUp;
		if(detbal[i].moveDown<minsumf)minsumf=detbal[i].moveDown;
	}
	offset=minsumf;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(maxsumf-minsumf);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-(0.0*plotinfo.scale_y);
	plotinfo.yscalemult=(maxsumf-minsumf)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
	plotinfo.ya=round_up((minsumf-offset)/plotinfo.yscalemult);
	plotinfo.yb=round_down((maxsumf-offset)/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"P(exchange)");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"P\\(exchange\\)");
	plotinfo.ligand_number=0;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	for(which=0; which<=2;which++){
		first=1;
		for(i=0;i<script->Nreplicas;i++){
			if(which==0){
				val=detbal[i].moveUp;
			}else if(which==1){
				val=detbal[i].moveDown;
			}else{
				val=detbal[i].moveSame;
			}
			printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(val-offset)*plotinfo.scale_y+plotinfo.translate_y);
			if(first) printf("moveto\n");
			else printf("lineto\n");
			first=0;
		}
		if(which==0){
			printf("4 setlinewidth\n");
			printf("1 0 0 setrgbcolor\n");
		}else if(which==1){
			printf("1 setlinewidth\n");
			printf("0 1 0 setrgbcolor\n");
		}else{
			printf("1 setlinewidth\n");
			printf("0 0 1 setrgbcolor\n");
			printf(DASHON);
		}
		//up=red; down=green; same=blue
		printf("stroke\n");
	}
	printf(DASHOFF);
}

void print_postscript_Pupdown(const struct pageinfo_struct *pageinfo,struct detailedBalance_struct *detbal, const struct script_struct *script, const struct graph_struct *graph){
	int i;
	unsigned int j;
	unsigned char first;
	float w;
	float force,minsumf,maxsumf,offset;
	double val;
	struct plotinfo_struct plotinfo=EMPTY_PLOTINFO;
	long totalSamples=0;
	float pi,pj,pij,pji,ioverj;

	for(i=0;i<script->Nreplicas;i++){
		totalSamples+=detbal[i].numSteps;
	}

	//Doesn't currently account for circular reaction coordinate
	first=1;
	for(i=0;i<script->Nreplicas-1;i++){
		pi=(float)detbal[i].numSteps/(float)totalSamples;
		pij=detbal[i].moveUp;
		pj=(float)detbal[i+1].numSteps/(float)totalSamples;
		pji=detbal[i+1].moveDown;
		ioverj=(pi*pij)/(pj*pji);
		if(first){
			minsumf=maxsumf=ioverj;
			first=0;
		}
		if(ioverj<minsumf)minsumf=ioverj;
		if(ioverj>maxsumf)maxsumf=ioverj;
		printf("%% [p(%d)*p(%d->%d)]/[p(%d)*p(%d->%d)] = %f\n",i,i,i+1,i+1,i+1,i,ioverj);
	}
	//offset=minsumf;
	// The y-scale values are entirely incorrect here
	offset=0.0;

	plotinfo.scale_x=(pageinfo->width-1.0-0.5)*72/(script->Nreplicas*(N_FORCE_POINTS+1));
	plotinfo.scale_y=((pageinfo->height-1.0-0.5-0.1-0.1)*72)/(maxsumf-minsumf);
	plotinfo.translate_x=(1.0*72);
	plotinfo.translate_y=(1.1*72)-((minsumf-offset)*plotinfo.scale_y);
	plotinfo.yscalemult=(maxsumf-minsumf)/(float)pageinfo->numvalues;
	plotinfo.xa=0;
	plotinfo.xb=script->Nreplicas;
	plotinfo.xc=1;
	plotinfo.ya=round_up((minsumf-offset)/plotinfo.yscalemult);
	plotinfo.yb=round_down((maxsumf-offset)/plotinfo.yscalemult);
	plotinfo.yc=1;
	sprintf(plotinfo.title,"\\[P\\(i\\)*P\\(i->i+1\\)\\]\\/\\[P\\(i+1\\)*P\\(i+1->i\\)\\]");
	if(script->coordinate_type==Temperature){
		sprintf(plotinfo.xtitle,"Temperature \\(K\\)");
	}else{
		sprintf(plotinfo.xtitle,"Reaction Coordinate");
	}
	sprintf(plotinfo.ytitle,"\\[P\\(i\\)*P\\(i->i+1\\)\\]\\/\\[P\\(i+1\\)*P\\(i+1->i\\)\\]");
	plotinfo.ligand_number=0;
	plotinfo.xlabels=(float *)NULL;

	showPlot(&plotinfo,pageinfo,script,graph);

	first=1;
	for(i=0;i<script->Nreplicas-1;i++){
		pi=(float)detbal[i].numSteps/(float)totalSamples;
		pij=detbal[i].moveUp;
		pj=(float)detbal[i+1].numSteps/(float)totalSamples;
		pji=detbal[i+1].moveDown;
		ioverj=(pi*pij)/(pj*pji);
		printf("%f %f ",((((float)i)+0.5)*(N_FORCE_POINTS+1))*plotinfo.scale_x+plotinfo.translate_x,(ioverj-offset)*plotinfo.scale_y+plotinfo.translate_y);
		if(first) printf("moveto\n");
		else printf("lineto\n");
		first=0;
	}
	printf("4 setlinewidth\n");
	printf("stroke\n");
}


void showPlot(const struct plotinfo_struct *plot, const struct pageinfo_struct *page, const struct script_struct *script, const struct graph_struct *graph){
	int i,j;
	char text[10];
	float val,range;

	if(page->showGrid){
		printf(DASHON);
		printf("%% Print X grid\n");
		for(i=plot->xa;i<plot->xb;i+=plot->xc){
			printf("0.5 setgray\n");
			printf("newpath\n");
			printf("%f %f moveto\n",((float)i+1)*(N_FORCE_POINTS+1)*plot->scale_x+plot->translate_x,(float)(1.0*72));
			printf("%f %f lineto\n",((float)i+1)*(N_FORCE_POINTS+1)*plot->scale_x+plot->translate_x,(((float)page->height-0.5)*72));
			printf("0.25 setlinewidth\n");
			printf("stroke\n");
		}

		printf("%% Print Y grid\n");
		range=(plot->yb-plot->ya)*plot->yscalemult;
		//CN was not able to round this well... he gave up
		for(i=plot->ya;i<=plot->yb;i+=plot->yc){
			printf("0.5 setgray\n");
			printf("newpath\n");
			val=(float)i*plot->yscalemult*plot->scale_y+plot->translate_y;
			printf("%f %f moveto\n",(1.0*72),val);
			printf("%f %f lineto\n",((page->width-0.5)*72),val);
			printf("0.25 setlinewidth\n");
			printf("stroke\n");
		}
		printf(DASHOFF);

		printf("%% Print Y=0 line\n");
		printf("0.5 setgray\n");
		printf("newpath\n");
		printf("%f %f moveto\n",(float)(1.0*72),plot->translate_y);
		printf("%f %f lineto\n",((page->width-0.5)*72),plot->translate_y);
		printf("1 setlinewidth\n");
		printf("stroke\n");
	}

	printf("%% Print graph box\n");
	printf("0.0 setgray\n");
	printf("newpath\n");
	printf("%f %f moveto\n",(float)(1.0*72),(float)(1.0*72));
	printf("%f %f lineto\n",((page->width-0.5)*72),(float)(1.0*72));
	printf("%f %f lineto\n",((page->width-0.5)*72),((page->height-0.5)*72));
	printf("%f %f lineto\n",(float)(1.0*72),((page->height-0.5)*72));
	printf("closepath\n");
	printf("1 setlinewidth\n");
	printf("stroke\n");

	printf("%% Print X scale\n");
	for(i=plot->xa;i<plot->xb;i+=plot->xc){
		printf("0.0 setgray\n");
		setFont(page->font_size);
		printf("%f %f moveto\n",(float)i*(N_FORCE_POINTS+1)*plot->scale_x+plot->translate_x,(float)(1.0*72)-page->font_size);
		if(i/2!=(i+1)/2) printf("0.0 -%f rmoveto\n",(float)page->font_size);
		if(plot->xlabels==NULL){
			/*
			 * T=particle_x[0];
			 * kT=(8.31451*T/4184.0);
			 * B=(1/kT);
			 * kT=(1/B);
			 * T=4184.0*kT/8.31451;
			 */
			if(script->coordinate_type==Temperature && plot->ligand_number==0){
				sprintf(text,"%0.2f",4184.0*(1/graph[i].w[plot->ligand_number])/8.31451);
			}else{
				sprintf(text,"%0.2f",graph[i].w[plot->ligand_number]);
			}
		}else{
			sprintf(text,"%0.2f",plot->xlabels[i]);
		}
		if(text[strlen(text)-1]=='0')text[strlen(text)-1]=0;
		if(text[strlen(text)-1]=='0')text[strlen(text)-1]=0;
		if(text[strlen(text)-1]=='.')text[strlen(text)-1]=0;
		printf("(%s) stringwidth pop\n",text);
		printf("%f exch sub 2 div\n",((float)N_FORCE_POINTS+1)*plot->scale_x);
		printf("0 rmoveto\n");
		printf("(%s) show\n",text);
	}
	printf("%% X title\n");
	setFont(page->title_font_size);
	printf("%f %f moveto\n",(float)(plot->xb-plot->xa)/2.0*(N_FORCE_POINTS+1)*plot->scale_x+plot->translate_x,(float)((0.7*72)-page->title_font_size));
	
	printf("(%s)dup stringwidth pop 2 div neg 0 rmoveto show\n",plot->xtitle);

	printf("%% Print Y scale\n");
	for(i=plot->ya;i<=plot->yb;i+=plot->yc){
		printf("0.0 setgray\n");
		setFont(page->font_size);
		val=(float)i*plot->yscalemult*plot->scale_y+plot->translate_y;
		printf("%f %f moveto\n",(1.0*72)-5,val-(float)page->font_size/2+1);
		if((plot->yb-plot->ya)*plot->yscalemult<1.0){
			printf("(%4.3f) stringwidth pop\n",(float)i*plot->yscalemult);
		}else if((plot->yb-plot->ya)*plot->yscalemult<10.0){
			printf("(%4.2f) stringwidth pop\n",(float)i*plot->yscalemult);
		}else if((plot->yb-plot->ya)*plot->yscalemult<100.0){
			printf("(%4.1f) stringwidth pop\n",(float)i*plot->yscalemult);
		}else{
			printf("(%4.0f) stringwidth pop\n",(float)i*plot->yscalemult);
		}
		printf("0.0 exch sub\n");
		printf("0 rmoveto\n");
		if((plot->yb-plot->ya)*plot->yscalemult<1.0){
			printf("(%4.3f) show\n",(float)i*plot->yscalemult);
		}else if((plot->yb-plot->ya)*plot->yscalemult<10.0){
			printf("(%4.2f) show\n",(float)i*plot->yscalemult);
		}else if((plot->yb-plot->ya)*plot->yscalemult<100.0){
			printf("(%4.1f) show\n",(float)i*plot->yscalemult);
		}else{
			printf("(%4.0f) show\n",(float)i*plot->yscalemult);
		}
	}
	printf("%% Y title\n");
	setFont(page->title_font_size);
	printf("%f %f moveto\n",(float)((0.7*72)-page->title_font_size),(float)(plot->yb+plot->ya)*plot->yscalemult/2.0*plot->scale_y+plot->translate_y-(float)page->title_font_size/2+1);
	printf("currentpoint gsave 90 rotate\n");
	printf("(%s)dup stringwidth pop 2 div neg 0 rmoveto show\n",plot->ytitle);
	printf("grestore\n");

	printf("%% Start of %s\n",plot->title);
	printf("0.0 setgray\n");
	printf("newpath\n");
}

/***********************************************************************************************/

int fileExists(const char *filename){
	return access(filename,F_OK);
}

signed char compare_records(struct record_struct *r1, struct record_struct *r2)
{
	if     (r1->replica_number  > r2->replica_number) return(1);
	else if(r1->replica_number  < r2->replica_number) return(-1);
	else if(r1->sequence_number > r2->sequence_number) return(1);
	else if(r1->sequence_number < r2->sequence_number) return(-1);
	else return(0);
}

void quicksort_database(struct database_struct *db, int left, int right){
	struct record_struct *pivot;
	int l_hold, r_hold;

	l_hold = left;
	r_hold = right;
	pivot = db->record[left];
	while (left < right){
		while ((compare_records(db->record[right],pivot) >= 0) && (left < right)) right--;
		if (left != right){
			db->record[left] = db->record[right];
			left++;
		}
		while ((compare_records(db->record[left],pivot) <= 0) && (left < right)) left++;
		if (left != right){
			db->record[right] = db->record[left];
			right--;
		}
	}
	db->record[left] = pivot;
	right = r_hold;
	r_hold = left;
	left = l_hold;
	if (left < r_hold)  quicksort_database(db, left, r_hold-1);
	if (right > r_hold) quicksort_database(db, r_hold+1, right);
}

struct record_struct* rec(unsigned int record_number,const struct database_struct *db){
	return((struct record_struct*)(((char*)db->unsorted_record)+db->record_size*record_number));
}

int read_in_database(const struct analysis_option_struct *opt, const struct script_struct *script, struct database_struct *db, const struct nominal_struct *nominal){
	class force_database_class *database;
	int i,j;
	unsigned long size;
	

	database=new force_database_class(opt->title,1);
	database->get_header_information(db);

	if(db->Nrecords==0){
		fprintf(stderr,"Error: the force database does not contain any records\n");
		return 1;
	}
	if(db->Nadditional_data!=script->Nadditional_data){
		fprintf(stderr,"Error: The forcedatabase has %u additional data types, while the input script specifies %u\n",db->Nadditional_data,script->Nadditional_data);
		return 1;
	}
	if(db->Nligands!=script->Nligands)
	{
		fprintf(stderr,"Error: the number of ligands in the database (%hhu) is not the number of ligands specified in the script file (%hhd)\n",db->Nligands,script->Nligands);
		return 1;
	}
	
	db->record_size=sizeof(record_struct)+(db->Nforces*db->Nligands+db->Nenergies+db->Nforces*db->Nadditional_data)*sizeof(float);
	size=db->Nrecords*db->record_size;
	fprintf(stderr,"The database size should be %lu excluding the header\n",size);
	/*
	if(size>MAX_DATABASE_SIZE){
		fprintf(stderr,"The database is too big; only the first %.2fGB will be considered\n",MAX_DATABASE_SIZE/1000000);	
		db->Nrecords=MAX_DATABASE_SIZE/db->record_size;
		size=db->Nrecords*db->record_size;
		fprintf(stderr,"Nrecords adjusted to %u\n",db->Nrecords);	
		fprintf(stderr,"%lu bytes of the database will be considered\n",size);
	}
	*/
	
	if( (db->unsorted_record=(struct record_struct*)malloc(db->Nrecords*db->record_size))==NULL ){
		fprintf(stderr,"Error: cannot allocated enough memory for the database\n");
		return 1;
	}

	if( (db->record=(struct record_struct**)malloc((db->Nrecords+1)*sizeof(struct record_struct *)))==NULL ){
		fprintf(stderr,"Error: cannot allocated enough memory for the database\n");
		return 1;
	}

	fprintf(stderr,"Reading records");
	for(i=0;i<db->Nrecords;i++){
		database->read_record_to_given_memory(i,rec(i,db));
		db->record[i]=rec(i,db);
		if( (i%1000)==999 ) fprintf(stderr,"."); fflush(stderr);
	}
	fprintf(stderr,"\n");
	
	fprintf(stderr,"Quicksorting database... ");
	quicksort_database(db, 0, db->Nrecords-1);
	fprintf(stderr,"done.\n");

	fprintf(stderr,"Removing duplicate records... ");
	for(i=0;i<db->Nrecords-1;i++){
		if(compare_records(db->record[i],db->record[i+1])==0){
			//void *memmove(void *s1, const void *s2, size_t n);
			//copies n bytes from the object pointed to by s2 into the object pointed to by s1
			//Therefore this method copies over the second instance
			//There is an assumption that it is the first instance that was used... is this true?
			fprintf(stderr, "\nEliminating duplicate record at position %d   ",i);
			memmove(db->record+i+1, db->record+i+2, ((int)db->Nrecords-i-2)*sizeof(struct record_struct *));
			i--;
			db->Nrecords--;
		}
	}
	fprintf(stderr,"done.\n");

	if(opt->sequence_number_limit>=0){
		unsigned int count;
		fprintf(stderr,"Removing records with sequence numbers beyond the limit... ");
		for(i=0;i<db->Nrecords;i++)
			if(db->record[i]->sequence_number>opt->sequence_number_limit){
				memmove(db->record+i, db->record+i+1, ((int)db->Nrecords-i-1)*sizeof(struct record_struct *));
				i--;
				db->Nrecords--;
				count++;
			}
		fprintf(stderr,"removed %u records.\n",count);
	}

	fprintf(stderr,"Confirming new order... ");
	for(i=0;i<db->Nrecords-1;i++){
		if(compare_records(db->record[i],db->record[i+1])!=-1){
			fprintf(stderr,"There is a problem with the new order when comparing the following records:\n");
			fprintf(stderr,"Record 1: replica number: %d   sequence_number: %u\n",db->record[i]->replica_number,db->record[i]->sequence_number);
			fprintf(stderr,"Record 2: replica number: %d   sequence_number: %u\n",db->record[i+1]->replica_number,db->record[i+1]->sequence_number);
			return 1;
		}
	}
	fprintf(stderr,"done.\n");

	db->replica_offset=db->record[0]->replica_number;

	if(opt->writeDatabase){
		// record_size=sizeof(record_struct)+(Nforces*Nligands+Nenergies+Nforces*Nadditional_data)*sizeof(float);
		FILE *f;
		int k;
		f=fopen(opt->writeDatabaseName,"w");
		if(f==NULL){
			fprintf(stderr,"Error: unable to open text database file for output, skipping\n");
			return 1;
		}
		for(i=0;i<db->Nrecords;i++){
			fprintf(f,"record: replica#: %d   sequence#: %u   w: %f   w_nominal: %d\n",db->record[i]->replica_number,db->record[i]->sequence_number,db->record[i]->w,find_bin_from_w(db->record[i]->w,0,nominal,script));
			if(opt->writeDatabaseFull){
				for(j=0;j<db->Nforces*db->Nligands;j+=db->Nligands){
					db->record[i]->generic_data[j];
					if(db->Nligands==1)      fprintf(f,"Force: %f\n",db->record[i]->generic_data[j]);
					else if(db->Nligands==2) fprintf(f,"Force: %f %f\n",db->record[i]->generic_data[j],db->record[i]->generic_data[j+1]);
				}
				if(db->Nadditional_data>0){
					for(j=db->Nforces*db->Nligands;j<db->Nforces*db->Nligands+db->Nforces;j++){
						fprintf(f,"Additional data: ");
						for(k=0;k<db->Nadditional_data;k++){
						        fprintf(f,"%f ",db->record[i]->generic_data[j+k*db->Nforces]);
						}
						fprintf(f,"\n");
					}
				}
			}
		}

		fclose(f);
	}
	return 0;
}

char get_next_force(float w[2], float force[2], float *rc, int whichRC, unsigned int *time, int *replica_number, unsigned short *sequence_number, const struct database_struct *db){
	static unsigned int recordN=0;
	static unsigned int sampleN=0;

	if(recordN>=db->Nrecords){
		recordN=0;
		return('f');
	}
	
	w[0]=db->record[recordN]->w;
	force[0]=0.0;
	force[1]=0.0;

	//BUG ALERT:
	//C. Neale doesn't think that this handles two ligands correctly
	if(db->Nforces>0){
		force[0]=db->record[recordN]->generic_data[sampleN*db->Nligands];
		if(whichRC>0){
			*rc=db->record[recordN]->generic_data[sampleN+db->Nforces*(db->Nligands+whichRC-1)];
		}
		if(db->Nligands==2){
			force[1]=db->record[recordN]->generic_data[sampleN*db->Nligands+1];
		}
	}
	*replica_number=db->record[recordN]->replica_number-db->replica_offset;
	*sequence_number=db->record[recordN]->sequence_number;
	*time=db->Nforces * (*sequence_number) + sampleN;
		
	sampleN++;

//BUG ALERT:
//C. Neale is sure that this won't work with 2 ligands now that all of
//the data from both ligands (not just 1st half) is stored in the forcedatabase
//as per the changes that C. Neale introduced in what he believes was a 'correction'
	if(sampleN*db->Nligands+1>=db->Nforces){
		sampleN=0;
		recordN++;
	}
	return('p');
}

void get_data_statistics(struct stats_struct *stats, const struct database_struct *db){
	float w[2];
	float force[2];
	unsigned int time;
	unsigned short sequence_number;
	int replica_number;
	int dummy;

	stats->max_time=0;
	
	while(1){
		if(get_next_force(w,force,(float *)NULL,0,&time,&replica_number,&sequence_number,db)=='f') break;
		if(time>stats->max_time) stats->max_time=time;
		if(sequence_number>stats->max_sequence_number) stats->max_sequence_number=sequence_number;
	}
}

void condense_forces(unsigned int N_values_to_average, struct nominal_struct *nominal, struct graph_struct *graph, const struct script_struct *script, const analysis_option_struct *opt, struct stats_struct *stats, const struct database_struct *db){
	unsigned int i,j;
	unsigned char l;
	float w[2];
	int replicaN;
	float force[2];
	float rc;
	float wa,wb;
	float fraction_a,fraction_b;
	unsigned int time;
	double weight_sum;
	unsigned short sequence_number;

	for(i=0;i<script->Nreplicas;i++){
		graph[i].w[0]=nominal[i].w[0];
		graph[i].w[1]=nominal[i].w[1];
		memset(graph[i].weight_sum,0,sizeof(graph[i].weight_sum));
		if(opt->verbose)fprintf(stderr,"Adding: Nreplicas: %u   w: %lf   w2: %lf\n",i,graph[i].w[0],graph[i].w[1]);
	}

	for(i=0;i<script->Nreplicas;i++){
		for(j=0;j<=stats->max_sequence_number;j++){	
			graph[i].replica_position[j]=graph[i].rc_position[j]=1e20;
		}
		graph[i].samples=0.0;
	}
	
	fprintf(stderr,"Number of discrete w positions: %u\n",script->Nreplicas);
	fprintf(stderr,"max time is: %u\n",stats->max_time);

	while(1){
		if(get_next_force(w,force,&rc,opt->additionalDataWithSampling,&time,&replicaN,&sequence_number,db)=='f') break;
		graph[replicaN].replica_position[sequence_number]=w[0];
		graph[replicaN].rc_position[sequence_number]=rc;

		time/=N_values_to_average;

		for(i=0;i<script->Nreplicas;i++){
			if(graph[i].w[0]>w[0]) break;
		}

		if(i==0)            wa=graph[0].w[0]-(graph[1].w[0]-graph[0].w[0]);
		else                wa=graph[i-1].w[0];
		if(i==script->Nreplicas) wb=graph[script->Nreplicas-1].w[0]+(graph[script->Nreplicas-1].w[0]-graph[script->Nreplicas-2].w[0]);
		else                wb=graph[i].w[0];
		fraction_b=(w[0]-wa)/(wb-wa);
		fraction_a=(double)1.0-fraction_b;
		if( (fraction_a>=-1e-5) && (fraction_a<=1.0+1e-5) ){
			if(i>0){
				graph[i-1].weight_sum[time]+=fraction_a;
				graph[i-1].samples+=fraction_a;
			}
			if(i<script->Nreplicas){
				graph[i].weight_sum[time]+=fraction_b;
				graph[i].samples+=fraction_b;
			}

			for(l=0;l<script->Nligands;l++){
				if(i>0){
					weight_sum=graph[i-1].weight_sum[time];
					if(weight_sum>1e-5) graph[i-1].point[l][time]=graph[i-1].point[l][time]*(1.0-fraction_a/weight_sum)+force[l]*(fraction_a/weight_sum);
				}

				if(i<script->Nreplicas){
					weight_sum=graph[i].weight_sum[time];
					if(weight_sum>1e-5) graph[i].point[l][time]=graph[i].point[l][time]*(1.0-fraction_b/weight_sum)+force[l]*(fraction_b/weight_sum);
				}
			}
		}
	}

	for(l=0;l<script->Nligands;l++){
		stats->max_force[l]=-1e10;
		stats->min_force[l]=1e10;
	}
	stats->max_samples=0.0;
	for(i=0;i<script->Nreplicas;i++){
		if(graph[i].samples>stats->max_samples) stats->max_samples=graph[i].samples;
		for(l=0;l<script->Nligands;l++){
			graph[i].average[l]=0.0;
			graph[i].avg_weight_sum[l]=0.0;
			weight_sum=0.0;

			//fprintf(stderr,"Weight sum of graph points at w = %f: %u\n",graph[i].w[0],graph[i].weight_sum);
			unsigned int max_force_points=0;
			for(j=0;j<N_FORCE_POINTS;j++) 
				if(graph[i].weight_sum[j]>1.0e-5) max_force_points=j;

			for(j=(unsigned int)((float)max_force_points*opt->equilFraction);j<N_FORCE_POINTS;j++){
				force[l]=graph[i].point[l][j];
				if(graph[i].weight_sum[j]>1.0e-5){
					//fprintf(stderr,"force is %f\n",force[l]);
					graph[i].average[l]+=force[l]*graph[i].weight_sum[j];
					graph[i].avg_weight_sum[l]+=graph[i].weight_sum[j];
					weight_sum+=1.0;
					//fprintf(stderr,"Total weight sum: %lf\n",avg_weight_sum[l]);
					if(force[l]<stats->min_force[l]) stats->min_force[l]=force[l];
					if(force[l]>stats->max_force[l]) stats->max_force[l]=force[l];
				}
			}
			//graph[i].avg_weight_sum[l]/=weight_sum;
			graph[i].average[l]/=graph[i].avg_weight_sum[l];
		}
	}
	
	for(l=0;l<script->Nligands;l++){
		fprintf(stderr,"max_force[%hhu]: %f   min_force[%hhu]: %f\n",l,stats->max_force[l],l,stats->min_force[l]);
	}
}


exact_struct * getExactFromFile(const char *title, exact_struct *exact, const struct script_struct *script){
	char command[500];
	char linein[1000];
	char test[50];
	int i;
	FILE *f;
	float val,eval;

	sprintf(test,"./%s.exact",title);
	if(access(test,F_OK)!=0)return (exact_struct *)NULL;   //return without error message

	f=fopen(test,"r");
	if(f==NULL){
		fprintf(stderr,"Error: unable to open the .exact file\n");
		return (exact_struct *)NULL;
	}
	i=0;
	while(fgets(linein,1000,f)!=NULL){
		if(linein[0]=='\0'||sscanf(linein,"%f %f",&val,&eval)!=2)continue;
		i++;
	}
	fclose(f);
	exact->n=i;

	/*
	if(i!=script->Nreplicas){
		if(script->coordinate_type!=Temperature){
			//file is expected to be of unknown size for temperature simulations
			fprintf(stderr,"Error: unexpected size of .exact file\n");
			return (exact_struct *)NULL;
		}
	}
	*/

	exact->e=(float *)malloc(i*sizeof(float));
	if(exact->e==NULL){
		fprintf(stderr,"Error: Unable to allocate memory for exact.e\n");
		return (exact_struct *)NULL;
	}
	exact->b=(float *)malloc(i*sizeof(float));
	if(exact->b==NULL){
		free(exact->e);
		fprintf(stderr,"Error: Unable to allocate memory for exact.b\n");
		return (exact_struct *)NULL;
	}

	f=fopen(test,"r");
	if(f==NULL){
		fprintf(stderr,"Error: unable to open the .exact file\n");
		return (exact_struct *)NULL;
	}
	i=0;
	while(fgets(linein,1000,f)!=NULL){
		if(linein[0]=='\0'||sscanf(linein,"%f %f",&(exact->b[i]),&(exact->e[i]))!=2)continue;
		i++;
	}
	fclose(f);

	return exact;
}

int getCancellationFromLog(const char *title, float *cancel, const struct script_struct *script){
	char command[500];
	char linein[1000];
	char test[50];
	int i;
	FILE *f;
	int multipleOccurrence=0;
	float val;

	sprintf(test,"%s.log",title);
	if(access(test,F_OK)!=0)return 1;   //return without error message
	sprintf(test,"./%s.cancel",title);

	sprintf(command,"grep ^coord: ./%s.log | awk '{print $4}' > %s",title,test);	
	//expected output format = "coord: 0.000000   energy: 0.000000"
	system(command);

	f=fopen(test,"r");
	if(f==NULL){
		fprintf(stderr,"Error: could not open the .cancel file\n");
		return 1;
	}
	i=0;
	while(fgets(linein,1000,f)!=NULL){
		if(linein[0]=='\0'||sscanf(linein,"%f",&val)!=1)continue;
		i++;
	}
	fclose(f);

	if(i==0)return 1; //there were no cancellation values

	if(i!=script->Nreplicas){
		if(i<script->Nreplicas||div(i,script->Nreplicas).rem!=0){
			fprintf(stderr,"Error: unexpected size of .cancel file\n");
			return 1;
		}
		//will need to test that all copies are identical
		multipleOccurrence=1;
	}

	f=fopen(test,"r");
	if(f==NULL){
		fprintf(stderr,"Error: could not open the .cancel file\n");
		return 1;
	}
	i=0;
	while(fgets(linein,1000,f)!=NULL){
		if(linein[0]=='\0'||sscanf(linein,"%f",&val)!=1)continue;
		if(i<script->Nreplicas){
			cancel[i]=val;
		}else{
			if(!floatEqual(cancel[div(i,script->Nreplicas).rem],val)){
				fprintf(stderr,"Error: repeated values in .cancel file do not match\n");
				return 1;
			}
		}
		i++;
	}
	fclose(f);
	return 0;
}

//Copied directly from DR_server.cpp (C. Neale Oct 17 2007)
int find_bin_from_w(double w, int l, const struct nominal_struct *nominal, const struct script_struct *script){
	double left_space, right_space;
	double left_limit, right_limit;
	int i;

	if(script->Nreplicas==1) return(0);

	for(i=0;i<script->Nreplicas;i++)
	{
		if(i>0)                   left_space=nominal[i].w[l]-nominal[i-1].w[l];
		if(i<script->Nreplicas-1) right_space=nominal[i+1].w[l]-nominal[i].w[l];

		if(i==0)                        left_space=right_space;
		else if(i==script->Nreplicas-1) right_space=left_space;

		left_limit=nominal[i].w[l]-left_space*0.5;
		right_limit=nominal[i].w[l]+right_space*0.5;

		if( (w>=left_limit) && (w<=right_limit) ) return(i);
	}
	return(-1);
}

int round_up(float x)
{
	if(x<0.0) return( (int)x );
	else      return( (int)(x+0.999999) );
}

int round_down(float x)
{
	if(x<0.0) return( (int)(x-0.999999) );
	else      return( (int)x );
}

void rot_trans_regular(void){
	printf("90 rotate\n");
	printf("0.0 -612.0 translate\n");
}

#define N_COLOURS 1000
char *get_RGB_colour(int colour_num)
{
	static float R_color_table[N_COLOURS];
	static float G_color_table[N_COLOURS];
	static float B_color_table[N_COLOURS];
	float R,G,B;
	static char colour_str[20];
	static unsigned char first=1;
	int i;

	colour_num=colour_num % N_COLOURS;

	if(first)
	{
		srand48(324531);
		first=0;
		for(i=0;i<N_COLOURS;i++)
		{
			do
			{
				R=drand48();
				G=drand48();
				B=drand48();
			} while ((R>0.5) && (G>0.5) && (B>0.5));
			
			R_color_table[i]=R;
			G_color_table[i]=G;
			B_color_table[i]=B;
		}
	}
	
	sprintf(colour_str,"%3.1f %3.1f %3.1f ",R_color_table[colour_num],G_color_table[colour_num],B_color_table[colour_num]);

	return(colour_str);
}

void setup_regular(unsigned char *current_page){
	printf("%%%%Page: %hhu %hhu\n\n",*current_page,*current_page);
	(*current_page)++;
	rot_trans_regular();
}


//*********************************************************************


void calcSequenceDensity(unsigned int **sequenceDensity, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const analysis_option_struct *opt, const struct stats_struct *stats){
// sequenceDensity[replica][sampling position]
// sequenceDensity[][0 to Nreplicas-1] stores values, [][Nreplicas] stores max value
// sequenceDensity[0 to Nreplicas-1][] stores individual replicas, [Nreplicas][] stores summation
	int i,j,bin;

	for(i=0;i<=script->Nreplicas; i++){
		for(j=0;j<=script->Nreplicas; j++){
			sequenceDensity[i][j]=0;
		}
	}

	for(i=0;i<db->Nrecords;i++){
		if(db->record[i]->sequence_number<(unsigned int)(opt->equilFraction*stats->max_sequence_number)){
			continue;
		}
		bin=find_bin_from_w(db->record[i]->w,0,nominal,script);
		++sequenceDensity[db->record[i]->replica_number][bin];
	}
	for(i=0;i<script->Nreplicas; i++){
		for(j=0;j<script->Nreplicas; j++){
			if(sequenceDensity[i][j]>sequenceDensity[i][script->Nreplicas])sequenceDensity[i][script->Nreplicas]=sequenceDensity[i][j];
			sequenceDensity[script->Nreplicas][j]+=sequenceDensity[i][j];
			if(sequenceDensity[script->Nreplicas][j]>sequenceDensity[script->Nreplicas][script->Nreplicas])sequenceDensity[script->Nreplicas][script->Nreplicas]=sequenceDensity[script->Nreplicas][j];
		}
	}
}

void getMinMaxSampleDensity(float *min, float *max, unsigned int whichData, unsigned int whichReplica, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script){
	//whichReplica == Nreplicas is a flag to get global min and max
	int i,j,bin;
	bool first=true;

	//Might be able to speed this up significantly since they are ordered
	//But I only use whichReplica=Nreplica flag so don't worry
	for(i=0;i<db->Nrecords;i++){
		bin=find_bin_from_w(db->record[i]->w,0,nominal,script);
		if(whichReplica!=script->Nreplicas && bin!=whichReplica)continue;
		for(j=script->Nsamples_per_run*(whichData); j<script->Nsamples_per_run*(whichData+1); j++){
			if(db->record[i]->generic_data[j]<*min||first)*min=db->record[i]->generic_data[j];
			if(db->record[i]->generic_data[j]>*max||first)*max=db->record[i]->generic_data[j];
			first=false;
		}
	}
}

int calcSampleDensity(struct sampleDensity_struct *sd, unsigned int whichData, const struct graph_struct *graph, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const struct analysis_option_struct *opt, const struct stats_struct *stats){
// sd->b[replica][histogram bin]
// sd->b[][0 to Nhisto-1] stores values, [][Nhisto] stores max value
// sd->b[0 to Nreplicas-1][] stores individual replicas, 
//      [Nreplicas][] stores summation
//      [Nreplicas][Nhisto] stores global max
// Use whichData == 1 for first additional data
// Data must be allocated for sd->b on [Nreplicas+1][sd->Ncol+1]

	int i,j;
	float val;
	float range;

	//Initialize sd->b and determine min and max 'sample' values
	for(i=0;i<=script->Nreplicas; i++){
		for(j=0;j<=sd->Ncol; j++){
			sd->b[i][j]=0;
		}
	}
	if(opt->discardOutside&&(script->coordinate_type==Spatial||script->coordinate_type==Umbrella)){
		range=graph[script->Nreplicas-1].w[0]-graph[0].w[0];
		sd->min=graph[0].w[0]-(range/sd->Ncol/2.0);
		sd->max=graph[script->Nreplicas-1].w[0]+(range/sd->Ncol/2.0);
	}else{
		getMinMaxSampleDensity(&(sd->min),&(sd->max),whichData,script->Nreplicas,db,nominal,script);
	}
	if(sd->max==sd->min){
		fprintf(stderr,"Error: sample max (%f) == sample min (%f). Exiting to avoid div by 0\n",sd->max,sd->min);
		return 1;
	}
	sd->binWidth=(sd->max-sd->min)/sd->Ncol;
	
	for(i=0;i<db->Nrecords;i++){
		for(j=script->Nsamples_per_run*(whichData); j<script->Nsamples_per_run*(whichData+1); j++){
			if(db->record[i]->sequence_number<(unsigned int)(opt->equilFraction*stats->max_sequence_number)){
				continue;
			}
			val=db->record[i]->generic_data[j];
			if(val<sd->min||val>sd->max)continue; //will only ocur if values are discarded
			++sd->b[db->record[i]->replica_number][(int)floor((val-sd->min)/sd->binWidth)];
		}
	}
	for(i=0;i<script->Nreplicas; i++){
		//Dump max value down to bin below
		sd->b[i][sd->Ncol-1]+=sd->b[i][sd->Ncol];
		sd->b[i][sd->Ncol]=0;
		//sum over all replicas
		for(j=0;j<sd->Ncol;j++){
			sd->b[script->Nreplicas][j]+=sd->b[i][j];
		}
	}
	//get the max values
	for(i=0;i<=script->Nreplicas; i++){
		sd->b[i][sd->Ncol]=sd->b[i][0];
		for(j=1;j<sd->Ncol;j++){
			if(sd->b[i][j]>sd->b[i][sd->Ncol])sd->b[i][sd->Ncol]=sd->b[i][j];
		}
	}
	return 0;
}

int calcPMF(struct pmf_struct *pmf, unsigned int whichData, const struct graph_struct *graph, const struct database_struct *db, const struct nominal_struct *nominal, const struct script_struct *script, const struct analysis_option_struct *opt){
// pmf->f[histogram bin]
// pmf->f[0 to Nhisto-1] stores values
// Use whichData == 1 for first additional data
// Data must be allocated for pmf->f and pmf->n on [pmf->Ncol+1]
//Important note: this is not going to be perfect for temperature since it includes all data without kB factor

	int i,j,k;
	float val;
	float range;

	//Initialize pmf->f and determine min and max 'sample' values
	for(i=0;i<=pmf->Ncol; i++){
		pmf->f[i]=0.0;
		pmf->n[i]=0;
	}
	if(opt->discardOutside&&(script->coordinate_type==Spatial||script->coordinate_type==Umbrella)){
		range=graph[script->Nreplicas-1].w[0]-graph[0].w[0];
		pmf->min=graph[0].w[0]-(range/(pmf->Ncol-1)/2.0);
		pmf->max=graph[script->Nreplicas-1].w[0]+(range/(pmf->Ncol-1)/2.0);
		fprintf(stderr,"PMF RANGE = %f (%f to %f)\n",range,pmf->min,pmf->max);
	}else{
		getMinMaxSampleDensity(&(pmf->min),&(pmf->max),whichData,script->Nreplicas,db,nominal,script);
		fprintf(stderr,"PMF RANGE = (%f to %f)\n",pmf->min,pmf->max);
	}
	if(floatEqual(pmf->max,pmf->min)){
		fprintf(stderr,"Error: sample max (%f) == sample min (%f). Exiting to avoid div by 0\n",pmf->max,pmf->min);
		return 1;
	}
	pmf->binWidth=(pmf->max-pmf->min)/pmf->Ncol;
	
	for(i=0;i<db->Nrecords;i++){
		for(j=script->Nsamples_per_run*(whichData),k=0; j<script->Nsamples_per_run*(whichData+1); j++,k++){
			val=db->record[i]->generic_data[j];
			if(val<pmf->min||val>pmf->max)continue; //will only ocur if values are discarded
			pmf->f[(int)floor((val-pmf->min)/pmf->binWidth)]+=db->record[i]->generic_data[k];
			++pmf->n[(int)floor((val-pmf->min)/pmf->binWidth)];
		}
	}
	//Dump max value down to bin below
	pmf->f[pmf->Ncol-1]+=pmf->f[pmf->Ncol];
	pmf->n[pmf->Ncol-1]+=pmf->n[pmf->Ncol];
	pmf->f[pmf->Ncol]=0.0;

	//Now create the average. Leave the n-value so can test for div zero later
	for(i=0;i<pmf->Ncol; i++){
		if(pmf->n[i]!=0){
			pmf->f[i]/=(float)pmf->n[i];
		}
	}
	return 0;
}

void setFont(unsigned int size){
	printf("/Helvetica findfont\n");
	printf("%d scalefont\n",size);
	printf("setfont\n");
}

#define TOL 0.000001
int floatEqual(float i,float j){
	if(i+TOL<j||i-TOL>j){
		return 0;
	}else{
		return 1;       //they are equal with the tolerance
	}
}

void showFirstComboPageText(int px, int py){
	px=400;py=190;
	setFont(12);
	printf("%d %d moveto\n",px,py);
	printf("(Left Top: )show\n");
	printf("currentpoint gsave\n");
	printf("(Potential of Mean Force \\(PMF\\)) show\n");
	printf("grestore\n");
	setFont(10);
	printf("0 -15 rmoveto\n");py-=15;
	printf("currentpoint gsave\n");
	printf("(simple integration \\(black\\) and full data re-integration \\(blue\\)) show\n");
	printf("grestore\n");
	printf("0 -15 rmoveto\n");py-=15;
	printf("(the dotted red line is the exact solution from the DR_tester if this is a test) show\n");
	setFont(12);
	py-=20;
	printf("%d %d moveto\n",px,py);
	printf("(Left Middle: relative INTENDED sample density \\(from wref\\)) show\n");
	py-=20;
	printf("%d %d moveto\n",px,py);
	printf("(Left Bottom: relative ACTUAL sample density \\(from additional data\\)) show\n");
	py-=20;
	printf("%d %d moveto\n",px,py);
	printf("(Right Top: ) show\n");
	printf("currentpoint gsave\n");
	printf("(PMF \\(thick line\\) and cancellation \\(thin black line\\)) show\n");
	printf("grestore\n");
	setFont(10);
	printf("0 -15 rmoveto\n");py-=15;
	printf("(the thin red line is a reflection/translation of the cancellation) show\n");
	setFont(12);
	py-=20;
	printf("%d %d moveto\n",px,py);
	printf("(Right Middle: ) show\n");
	printf("currentpoint gsave\n");
	printf("(\\(thick line\\) Sum of black lines from top right) show\n");
	printf("grestore\n");
	printf("0 -20 rmoveto\n");py-=20;
	printf("(\\(thin line\\) Derivative of thick line) show\n");
}

void showSampleDensityNotAvailable(int px,int py,const struct script_struct *script,const struct analysis_option_struct *opt){
	setFont(12);
	px=80;py=140;	
	printf("%d %d moveto\n",px,py);
	printf("(Actual sampling not available) show\n");
	setFont(10);
	py-=20;printf("%d %d moveto\n",px,py);
	if(!(script->coordinate_type==Umbrella || script->coordinate_type==Temperature)){
		printf("(Coordinate type must be umbrella or temperature) show\n");
	}
	py-=20;printf("%d %d moveto\n",px,py);
	if(!script->Nadditional_data>0){
		printf("(Nadditional_data must be > 0) show\n");
	}
	py-=20;printf("%d %d moveto\n",px,py);
	if(!opt->additionalDataWithSampling>0){
		printf("(command line option selecting additional data must be > 0) show\n");
	}
}

void showSecondComboPageText(int px,int py){
	setFont(12);
	printf("%d %d moveto\n",px,py);
	printf("(Left Middle: Detailed balance satisfaction)show\n");
	py-=20;
	printf("%d %d moveto\n",px,py);
	printf("(Left Bottom: ) show\n");
	printf("currentpoint gsave\n");
	printf("(Exchange probability for) show\n");
	printf("grestore\n");
	printf("0 -15 rmoveto\n");
	printf("currentpoint gsave\n");
	printf("(i -> i + 1 \\(thick red\\),) show\n");
	printf("grestore\n");
	printf("0 -15 rmoveto\n");
	printf("currentpoint gsave\n");
	printf("(i -> i - 1 \\(thin green\\),) show\n");
	printf("grestore\n");
	printf("0 -15 rmoveto\n");
	printf("(i -> i \\(broken blue\\)) show\n");
}
