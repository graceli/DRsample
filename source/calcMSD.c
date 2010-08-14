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

#include <stdio.h>
#include <stdlib.h>

#define LINESIZE 1000

int main(int argn, char *args[]){
	FILE *f;
	char linein[LINESIZE];
	int maxs,mins,i,j,maxr,n,s,oldr,oldn,minn,r,waitingN,olds;
	float w;
	float **rec;
	float *sd;
	int *count;

	f=fopen("database.txt","r");
	if(f==NULL){
		printf("Error, unable to open database.txt\n");
		exit(1);
	}
	oldr=0;
	oldn=0;
	olds=0;
	maxs=0;
	mins=9999999;
	minn=0;
	waitingN=1;
	while(fgets(linein,LINESIZE,f)!=NULL){
		//record: replica#: 0   sequence#: 0   w: 0.672750   w_nominal: 0
		if(sscanf(linein,"%*s %*s %d %*s %d %*s %f %*s %d",&r,&s,&w,&n)!=4)continue;
    if(waitingN&&n!=oldn){
      waitingN=0;
      if(s>minn)minn=s;
    }
		if(r!=oldr){
			if(olds>maxs)maxs=olds;
			if(olds<mins)mins=olds;
			waitingN=1;
		}
		oldr=r;
		oldn=n;
		olds=s;
	}
	fclose(f);
	maxr=r;
	//minn=203        mins=0  maxs0   maxr=30

	fprintf(stderr,"minn=%d\tmins=%d\tmaxs=%d\tmaxr=%d\n",minn,mins,maxs,maxr);

	rec=(float **)malloc((maxr+1)*sizeof(float *));
	if(rec==NULL){
		printf("Error: unable to allocate memory for rec\n");
		exit(1);
	}
	for(r=0;r<=maxr;r++){
		rec[r]=(float *)malloc((maxs+1)*sizeof(float));
		if(rec[r]==NULL){
			printf("Error: unable to allocate memory for rec[r]\n");
  	  exit(1);
	  }
	}
	count=(int *)calloc(mins-minn+1,sizeof(int));
	sd=(float *)calloc(mins-minn+1,sizeof(float));
	if(count==NULL||sd==NULL){
		printf("Error: unable to allocate memory for count or sd (asked for %d-%d slots)\n",mins,minn);
    exit(1);
  }

  f=fopen("database.txt","r");
  if(f==NULL){
    printf("Error, unable to open database.txt\n");
    exit(1);
  }
	while(fgets(linein,LINESIZE,f)!=NULL){
    //record: replica#: 0   sequence#: 0   w: 0.672750   w_nominal: 0
    if(sscanf(linein,"%*s %*s %d %*s %d %*s %f %*s %d",&r,&s,&w,&n)!=4)continue;
		rec[r][s]=w;
	}
	fclose(f);

	for(r=0;r<=maxr;r++){
		for(i=minn;i<mins;i++){
			for(j=i+1;j<=mins;j++){
				sd[j-i]+=(rec[r][i]-rec[r][j])*(rec[r][i]-rec[r][j]);
				count[j-i]++;
			}
		}
	}

	printf("N\tMSD\tpoints\n");
	for(i=1;i<=mins-minn;i++){
		if(count[i]!=0){
			printf("%d\t%f\t%d\n",i,sd[i]/(float)count[i],count[i]);
		}else{
			printf("%d\t%f\t%d\n",i,-1.0,count[i]);
		}
	}
}
