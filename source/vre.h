#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#define DEFAULT_NUMSAVES_PRIMARY 100000
#define DEFAULT_NUMSAVES_SECONDARY 1000

long int actual_numsaves_secondary=DEFAULT_NUMSAVES_SECONDARY;

// Begin vre Primary data storage ---------------------------

struct vre_item_struct{
	float val;
	int source; 
};

struct vre_struct{
	struct vre_item_struct *n;
	long int nallocated;
	long int nlastused;
};

static struct vre_struct *vre;

// End vre Primary data storage -----------------------------

// Begin vre Secondary data storage -------------------------

struct secondary_vre_struct{
	float *val;
	long int nallocated;
	long int nlastused;
	long int nrecyclepush;
};

static struct secondary_vre_struct *secv;

// End vre Secondary data storage ---------------------------

void save_vre_pointers(int numnominal, long int *nallocated, long int *nlastused, struct vre_item_struct **item);
void load_vre_pointers(int numnominal, long int **nallocated, long int **nlastused, struct vre_item_struct ***item);
void save_secvre_pointers(int numnominal, long int *nallocated, long int *nlastused, long int *nrecyclepush, float **val);
void load_secvre_pointers(int numnominal, long int **nallocated, long int **nlastused, long int **nrecyclepush, float ***val);
int allocateVRE_primary(int numnominal, int numsaves);
int allocateVRE_secondary(int numnominal, int numsaves);
int allocateVRE(int numnominal, int numsaves);
int popVRE(int moveto, int thisrep, float *popped, int *source);
void pushVRE(int thisnominal, int thisrep, float pushed);
long int getVREallocationLevel(int numReplicas);
long int getSECVREallocationLevel(int numReplicas);
void showVRE(int numReplicas, FILE *f);
int loadFileIntoVREforStartup(int replicaN, const char *fnam);
void set_secvre_size(long int val);

void save_vre_pointers(int numnominal, long int *nallocated, long int *nlastused, struct vre_item_struct **item){
	*nallocated=vre[numnominal].nallocated;
	*nlastused=vre[numnominal].nlastused;
	*item=vre[numnominal].n;
}

void load_vre_pointers(int numnominal, long int **nallocated, long int **nlastused, struct vre_item_struct ***item){
	*nallocated=&(vre[numnominal].nallocated);
	*nlastused=&(vre[numnominal].nlastused);
	*item=&(vre[numnominal].n);
}

void save_secvre_pointers(int numnominal, long int *nallocated, long int *nlastused, long int *nrecyclepush, float **val){
	*nallocated=secv[numnominal].nallocated;
	*nlastused=secv[numnominal].nlastused;
	*nrecyclepush=secv[numnominal].nrecyclepush;
	*val=secv[numnominal].val;
}

void load_secvre_pointers(int numnominal, long int **nallocated, long int **nlastused, long int **nrecyclepush, float ***val){
	*nallocated=&(secv[numnominal].nallocated);
	*nlastused=&(secv[numnominal].nlastused);
	*nrecyclepush=&(secv[numnominal].nrecyclepush);
	*val=&(secv[numnominal].val);
}

int allocateVRE_primary(int numnominal, int numsaves){
	int i;
	if(numsaves<0)numsaves=DEFAULT_NUMSAVES_PRIMARY;

	vre=(struct vre_struct *)malloc(numnominal*sizeof(struct vre_struct));
	if(vre==NULL){
		return -1;
	}
	for(i=0;i<numnominal;i++){
		vre[i].nallocated=numsaves;
		vre[i].nlastused=-1;
		vre[i].n=(struct vre_item_struct *)malloc(numsaves*sizeof(struct vre_item_struct));
		if(vre[i].n==NULL){
			return -1;
		}
	}
	return 0;
}

int allocateVRE_secondary(int numnominal, int numsaves){
	int i;
	if(numsaves<0)numsaves=DEFAULT_NUMSAVES_SECONDARY;

	secv=(struct secondary_vre_struct *)malloc(numnominal*sizeof(struct secondary_vre_struct));
	if(secv==NULL){
		return -1;
	}
	for(i=0;i<numnominal;i++){
		secv[i].nallocated=numsaves;
		secv[i].nlastused=-1;
		secv[i].nrecyclepush=-1;
		secv[i].val=(float *)malloc(numsaves*sizeof(float));
		if(secv[i].val==NULL){
			return -1;
		}
	}
	return 0;
}

int allocateVRE(int numnominal, int numsaves){
  int checka,checkb;

	checka=allocateVRE_primary(numnominal,numsaves);
  checkb=allocateVRE_secondary(numnominal,numsaves);

  return (checka||checkb);
}

int popVRE(int moveto, int thisrep, float *popped, int *source){
/* popVRE finds a cancelation value and removes it from the list 
 *  - moveto is the nominal index to which a move may be made
 *  - thisrep is the replica index of the current sampling
 *  - popped is used to return the cancellation value
 *  - source is used to return the replica from which this value was derived.
 *    This is not necessary for execution, but allows more information to be output.
 *    source will be set to -1 when the popped value was from the secondary list.
 */
	long int vp,svp;

	for(vp=vre[moveto].nlastused;vp>=0 && vre[moveto].n[vp].source==thisrep; vp--);
	if(vp>=0){
		/* Success:
		 * Return the popped value via a pointer
		 * Put the popped value into the secondary data structure,
		 * Replace the popped value by the one at the end of the list
		 * Reduce the number of entries by one.
		 */ 
		*popped=vre[moveto].n[vp].val;
		*source=vre[moveto].n[vp].source;
		if(++(secv[moveto].nlastused)<secv[moveto].nallocated){
			secv[moveto].val[secv[moveto].nlastused]=*popped;
		}else{
			//secv structure is full
			secv[moveto].nlastused--;
			//New code for v2.1.2 to cycle about used values
			if(secv[moveto].nrecyclepush==-1||secv[moveto].nrecyclepush>secv[moveto].nlastused){
				secv[moveto].nrecyclepush=0;
			}
			secv[moveto].val[secv[moveto].nrecyclepush]=*popped;
			secv[moveto].nrecyclepush++;
		}
		vre[moveto].n[vp].val=vre[moveto].n[vre[moveto].nlastused].val;
		vre[moveto].n[vp].source=vre[moveto].n[vre[moveto].nlastused].source;
		vre[moveto].nlastused--;
		return 0;
	}
	//Did not find an unused entry. Now randomly pick one from the secondary array.
	if(secv[moveto].nlastused<0){
		//There are absolutely no values available
		return -1;
	}
	//CN CHECK LATER THAT THIS drand48 IS RANDOMIZED BY MAIN srand48 CALL
	svp=(long int)ceil(drand48()*(double)secv[moveto].nlastused);
	*popped=secv[moveto].val[svp];
	*source=-1;
	return 0; 
}

void pushVRE(int thisnominal, int thisrep, float pushed){
/* pushVRE adds a new cancellation value for later use
 *  - thisnominal is the current nominal index of the replica
 *  - thisrep is the replica index of the current sampling
 *  - pushed carries the cancellation value
 */

	if(++(vre[thisnominal].nlastused)<vre[thisnominal].nallocated){
		vre[thisnominal].n[vre[thisnominal].nlastused].source=thisrep;
		vre[thisnominal].n[vre[thisnominal].nlastused].val=pushed;
	}else{
		//there is no space left in the array. Not really an error. There may be more space later
		vre[thisnominal].nlastused--;
	}
}

long int getVREallocationLevel(int numReplicas){
	int i;
	long int minAllocated;

	minAllocated=vre[0].nallocated;
	for(i=0;i<numReplicas;i++){
		if(vre[i].nallocated<minAllocated) minAllocated=vre[i].nallocated;
	}

	return minAllocated;
}

long int getSECVREallocationLevel(int numReplicas){
        int i;
        long int minAllocated;

        minAllocated=secv[0].nallocated;
        for(i=0;i<numReplicas;i++){
                if(secv[i].nallocated<minAllocated) minAllocated=secv[i].nallocated;
        }

        return minAllocated;
}


void showVRE(int numReplicas, FILE *f){
	int i,j;

	fprintf(f,"\n\nTHE VRE PRIMARY STRUCTURE\n");
	for(i=0;i<numReplicas;i++){
		fprintf(f,"Nominal position %d\n",i);
		fprintf(f,"Nallocated %ld\n",vre[i].nallocated);
		fprintf(f,"NlastUsed  %ld\n",vre[i].nlastused);
		for(j=0;j<=vre[i].nlastused;j++){
			fprintf(f,"%f %d\n",vre[i].n[j].val,vre[i].n[j].source);
		}
	}
	fprintf(f,"\n\nTHE VRE SECONDARY STRUCTURE\n");
	for(i=0;i<numReplicas;i++){
		fprintf(f,"Nominal position %d\n",i);
		fprintf(f,"Nallocated %ld\n",secv[i].nallocated);
		fprintf(f,"NlastUsed  %ld\n",secv[i].nlastused);
		fprintf(f,"NrecyclePush  %ld\n",secv[i].nrecyclepush);
		for(j=0;j<=secv[i].nlastused;j++){
			fprintf(f,"%f\n",secv[i].val[j]);
		}
	}
}

int loadFileIntoVREforStartup(int replicaN, const char *fnam){
	FILE *f;
	char linein[100];
	float val;

	if(fnam[0]=='\0')return 0;

	f=fopen(fnam,"r");
	if(f==NULL)return 1;
	while(fgets(linein,100,f)!=NULL){
		if(sscanf(linein,"%f",&val)!=1)continue;
		if(++(vre[replicaN].nlastused)>=vre[replicaN].nallocated){
			--vre[replicaN].nlastused;
			fclose(f);
			return 1;
		}
		vre[replicaN].n[vre[replicaN].nlastused].val=val;
		vre[replicaN].n[vre[replicaN].nlastused].source=-9;
	}

	fclose(f);
	return 0;
}

void set_secvre_size(long int val){
        if(val<0){
                actual_numsaves_secondary=DEFAULT_NUMSAVES_SECONDARY;
        }else{
                actual_numsaves_secondary=val;
        }
}

