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

#include <stdlib.h>
#include <stdio.h>

#define MAXJOBS 250

void showUsage(const char *c){
	fprintf(stderr,"Usage: %s <script file> <database text file> <sequence numbers to skip> <max sequence number to use> <data to use> > wham.input\n",c);
	fprintf(stderr,"       - The database text file can be obtained from analyze_force_database -d database.txt -t 1\n");
	fprintf(stderr,"       - the redirected output (wham.input above) can be used for Alan Grossfield's WHAM\n");
	fprintf(stderr,"       - to use the full data, set skip=0 and max<0\n");
	fprintf(stderr,"       - <data to use> : 1 is first additional data, 2 is second additional data...\n");
}

int main(int argn, char *args[]){
	char command[1001];
	char linein[1001];
	FILE *f;
	float w,data;
	int s,n,i,r;
	FILE **g;
	char name[MAXJOBS][100];
	float jobs[MAXJOBS];
	double fc[MAXJOBS];
	int njobs;
	int nskip=0;
	int nmax=-1;
	int dtu=1;

	if(argn!=6){
		showUsage(args[0]);
		exit(1);
	}
	sscanf(args[3],"%d",&nskip);
	sscanf(args[4],"%d",&nmax);
	sscanf(args[5],"%d",&dtu);
	

	f=fopen(args[1],"r");
	if(f==NULL){
		printf("Error: unable to open file to read\n");
		exit(1);
	}
	njobs=0;
	while(fgets(linein,1000,f)!=NULL){
		//JOB   1.30    10000    11      119.50286806883365200764
		if(linein[0]!='J' || linein[1]!='O' || linein[2]!='B')continue;
		if(sscanf(linein,"%*s %f %*d %*d %lf",&jobs[njobs],&fc[njobs])!=2){
			printf("Error: improperly formatted line in script <<%s>>\n",linein);
			printf("  expect it to look like this: JOB   1.30    10000    11      119.50286\n");
			fclose(f);
			exit(1);
		}
		njobs++;
	}
	fclose(f);

	f=fopen(args[2],"r");
	if(f==NULL){
		printf("Error: unable to open file to read\n");
		exit(1);
	}

	g=(FILE **)malloc(njobs*sizeof(FILE *));
	if(g==NULL){
		printf("Error: unable to malloc\n");
		exit(1);
	}

	for(i=0;i<njobs;i++){
		sprintf(name[i],"%.6f.wham.data",jobs[i]);	
		g[i]=fopen(name[i],"w");
		if(g[i]==NULL){
			printf("Error: unable to open file to write\n");
			exit(1);
		}
	}
	while(fgets(linein,1000,f)!=NULL){
		if(linein[0]=='r'&&linein[1]=='e'&&linein[2]=='c'){
			//record: replica#: 0   sequence#: 0   w: 0.200000   w_nominal: 0
			if(sscanf(linein,"%*s %*s %d %*s %d %*s %f %*s %d\n",&r,&s,&w,&n)!=4){
				printf("Error: unable to find r,s,w and n\n");
				exit(1);
			}
			continue;
		}
		if(linein[0]=='A'&&linein[1]=='d'&&linein[2]=='d'){
			//Additional data: 0.299169 0.298097 0.204570 0.332778 0.336268 0.370613
			switch(dtu)
			{
			case 1:
				if(sscanf(&(linein[17]),"%f",&data)!=1){
					printf("Error, did not find datapoint\n");
					exit(1);
				}
				break;
			case 2:
				if(sscanf(&(linein[17]),"%*f %f",&data)!=1){ 
					printf("Error, did not find datapoint\n"); 
					exit(1); 
				} 
				break;
			case 3:
				if(sscanf(&(linein[17]),"%*f %*f %f",&data)!=1){
					printf("Error, did not find datapoint\n"); 
					exit(1); 
				} 
				break;
			case 4:
				if(sscanf(&(linein[17]),"%*f %*f %*f %f",&data)!=1){
					printf("Error, did not find datapoint\n"); 
					exit(1); 
				} 
				break;
			case 5:
				if(sscanf(&(linein[17]),"%*f %*f %*f %*f %f",&data)!=1){
					printf("Error, did not find datapoint\n"); 
					exit(1); 
				} 
				break;
			case 6:
				if(sscanf(&(linein[17]),"%*f %*f %*f %*f %*f %f",&data)!=1){
					printf("Error, did not find datapoint\n"); 
					exit(1); 
				} 
				break;
			default:
				printf("Error: program can not handle more than 6 additional data types. Must edit the code and recompile.\n");
				exit(1);
			}
			if(s>nskip && (nmax<0 || s<nmax)){
				fprintf(g[n],"1 %f %d %d\n",data,r,s);
			}
			continue;
		}
	}
	fclose(f);
	for(i=0;i<njobs;i++){
		fclose(g[i]);
	}

	//Now sort the files by the 4th column
	for(i=0;i<njobs;i++){
		//sprintf(command,"sort -n -k 4 ./%s | awk '{print $1,$2}' > sorting.tmp",name[i]);
		// the above line gave me strange sorting on the valued column...
		sprintf(command,"cat -n ./%s|sort -n -k 5 | awk '{print $2,$3}' > sorting.tmp",name[i]);
		system(command);
		sprintf(command,"mv sorting.tmp ./%s",name[i]);
		system(command);

	}

	for(i=0;i<njobs;i++){
		printf("%s %f %.20lf\n",name[i],jobs[i],fc[i]);
	}

	fprintf(stderr,"Try this: /hpf/projects1/pomes/cneale/exe/wham_AlanGrossfield/wham/wham %f %f 100 0.00001 298.0 0 wham.input wham.output\n",jobs[0],jobs[njobs-1]);
}
