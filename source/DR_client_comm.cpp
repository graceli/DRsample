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

/*
 Chris Neale Oct 2 2007
   - Added ability to move additional files like the job.force file -- this currently does not impact the forcedatabase
   - SendBinFile now returns an integer to flag if the file was not found for use with additional_data files

 Chris Neale June 6 2008
   - Added ability to deal with NNI!=1

 Chris Neale June 6 2008
   - Added free() statements for all mallocs

*/

#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
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
#include "read_input_script_file.h" //this means that it requires math.h (that header has NAN)
#include "string_double.h"


#define IDSIZE sizeof(struct ID_struct)
#define INTSIZE sizeof(int)

#define COMPRESS_RESTART true

unsigned int do_compress2(int ifd, char *ocp);
void do_uncompress2(int total, char *icp, int ofd);
void sendFile(int sockfd, char *filename, enum command_enum command, bool compress);
int sendBinFile(int sockfd, char *filename, enum command_enum command);
void sendCrdFile(int sockfd, char *filename, enum command_enum command);
void sendJID(int sockfd, float jid);
void sendTCS(int sockfd, float tcs);
void sendReplicaID(int sockfd, struct ID_struct ID);
void receiveReplicaID(int sockfd, struct ID_struct *ID);
void send_NextNonInteracting(int sockfd);
enum command_enum readCommand(int socketfd);
void receiveFile(int sockfd, char *fileName, int fileSize);
void receiveFileUncompressed(int sockfd, char *fileName, int fileSize);
void receiveParCHARMM(int sockfd, struct ID_struct *ID, int fileSize);
void read4K(int sockfd, void *buff, int nbytes);

using namespace std;

void showUsage(const char *c){
	fprintf(stderr,"Usage: %s  IP-address  port  replicaIDcode  time-client-started  job-id\n",c);
	fprintf(stderr,"OR:    %s  IP-address  port  '**'           time-client-started  job-id\n",c);
	fprintf(stderr,"             * time-client-started is manditory, but is only used with a mobile server (send 0 if you don't care)\n");
	fprintf(stderr,"                                   the time should be given in seconds since January 1, 1970.\n");
	fprintf(stderr,"                                   if you send <=0 then the server will track times internally (sub-optimal for mobile server).\n");
	fprintf(stderr,"             * job-id is manditory, but is only used for tracking (send 0 if you don't know the job ID on the client)\n");
}

int main(int argc, char *argv[]){
	int sockfd; 
	int jid,tcs;
	struct sockaddr_in their_addr; // connector's address information 

	if (argc != 6) {
		showUsage(argv[0]);
		exit(1);
	}
	sscanf(argv[4],"%d",&tcs);
        sscanf(argv[5],"%d",&jid);
	if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
		perror("socket");
		exit(1);
	}

	their_addr.sin_family = AF_INET;    // host byte order 
	their_addr.sin_addr.s_addr = inet_addr(argv[1]);
	their_addr.sin_port = htons(atoi(argv[2]));  // short, network byte order 
	memset(&(their_addr.sin_zero), '\0', 8);  // zero the rest of the struct 

	if(connect(sockfd, (struct sockaddr *)&their_addr, sizeof(struct sockaddr)) == -1) {
		perror("connect");
		exit(1);
	}

	/************************************************************************/

	struct ID_struct ID;
	unsigned int baseNameLength=strlen(argv[3]);

	char rstFileName[64];
	char energyFileName[64];
	char forceFileName[64];
	char addFileName[64];
	char crdFileName[64];
	char tmpName[64];

	char IDbaseName[60];
	char IDbaseName2[60];        //decrease sequence_number by 1 for restart file
	char IDrstFileName[80];
	char IDenergyFileName[80];
	char IDforceFileName[80];
	//char IDsetupFileName[80];
	unsigned int protocol_version=PROTOCOL_VERSION;
	int old_replica_number=-1;
	int new_replica_number=-1;
	int addfile,sbf;
	char saveArgv3[50];
	int nni;

	write(sockfd,&protocol_version,PROTOCOL_VERSION_SIZE); // first thing we do is send the protocol version we are using

	if(argv[3][0]=='*' && argv[3][1]=='*'){
		ID.title[0]=ID.title[1]='*';
		sendTCS(sockfd,(float)tcs);
		sendJID(sockfd,(float)jid);
		sendReplicaID(sockfd,ID);
	}else{
		//sendID

		if(argv[3][2]!='w' && argv[3][2]!='W' ){ fprintf(stderr,"wrong id\n"); exit(1); }

		for(int i=0; i<3; ++i)
			if(isupper(argv[3][i])) (argv[3][i])+=('a'-'A');

		strcpy(rstFileName,argv[3]);
		strcpy(energyFileName,argv[3]);
		strcpy(forceFileName,argv[3]);
		strcpy(crdFileName,argv[3]);
		strcat(rstFileName,".rst");
		strcat(energyFileName,".energy");
		strcat(forceFileName,".force");
		strcat(crdFileName,".crd");
		strcpy(saveArgv3,argv[3]);   //keep it for iterations with addFileName
		argv[3][2]=' ';
		int len=strlen(argv[3]);

		for(int i=3; i<len; ++i){
			if(argv[3][i]=='.') argv[3][i]=' ';
		}

		sscanf(argv[3],"%s%d%u",ID.title,&ID.replica_number,&ID.sequence_number);
		old_replica_number=new_replica_number=ID.replica_number;
/*
//CN removed this as it is unexpected
		if(ID.sequence_number>0){
			char old_RST_filename[60];
			sprintf(old_RST_filename,"%sw%d.%u.rst",ID.title,ID.replica_number,ID.sequence_number-1);
			fprintf(stderr,"deleting old RST file [%s]\n",old_RST_filename); //##DEBUG
			int e=unlink(old_RST_filename);  // if we are receiving a new RST file then delete old one
			if(e!=0) fprintf(stderr,"delete failed with errno %d\n",errno); //##DEBUG
		}
*/

		fprintf(stderr,"sending ID\n"); //##DEBUG
		sendReplicaID(sockfd,ID);                              //send ReplicaID
		fprintf(stderr,"ID sent\n"); //##DEBUG

		sendTCS(sockfd,(float)tcs);
		sendJID(sockfd,(float)jid);
		for(nni=1;;nni++){
			sprintf(tmpName,"%s.nni%d",energyFileName,nni);
			fprintf(stderr,"sending energy file\n"); //##DEBUG
			sbf=sendBinFile(sockfd,tmpName,TakeMoveEnergyData); //send move_energy_data
			if(sbf!=0){
				//If the file did not exist then assume that all is well.
				fprintf(stderr,"  * Above file does not exist, assuming end of non-interacting systems, this is not an error message.\n"); //##DEBUG
				break;
			}
			fprintf(stderr,"energy file sent\n"); //##DEBUG
			sprintf(tmpName,"%s.nni%d",forceFileName,nni);
			fprintf(stderr,"sending force file\n"); //##DEBUG
			sendBinFile(sockfd,tmpName,TakeSampleData);      //send SampleData, force file
			fprintf(stderr,"force file sent\n"); //##DEBUG
			//It is absolutely essential that the sample data is sent before additional data
			//This loop will break when a numbered file des not exist
			for(addfile=1;;addfile++){
				sprintf(tmpName,"%s.add%d.nni%d",saveArgv3,addfile,nni);
				fprintf(stderr,"sending additional file %s\n",tmpName); //##DEBUG
				sbf=sendBinFile(sockfd,tmpName,TakeSampleData);      //send additional data file(s) 
				if(sbf!=0){
					//If the file did not exist then assume that all is well.
					fprintf(stderr,"  * Above file does not exist, assuming end of additional data, this is not an error message.\n"); //##DEBUG
					break;
				}
				fprintf(stderr,"additional file %s sent\n",tmpName); //##DEBUG
			}

			fprintf(stderr,"sending crd file\n"); //##DEBUG
			sprintf(tmpName,"%s.nni%d",crdFileName,nni);
			sendCrdFile(sockfd,tmpName,TakeCoordinateData);    //send CoordinateData, force crd file
			fprintf(stderr,"crd file sent\n"); //##DEBUG
			fprintf(stderr,"sending rst file\n"); //##DEBUG
			if(nni==1){
				sendFile(sockfd,rstFileName,TakeRestartFile,COMPRESS_RESTART);     //send restart file
				fprintf(stderr,"rst file sent\n"); //##DEBUG
			}else{
				//send indication of next NNI
				send_NextNonInteracting(sockfd);  //line added by CN August 6 2008
			}
			//send TakeThisFile
			//sendFile(sockfd,"test.123",TakeThisFile,false);
		}
	}
	enum command_enum command;
	bool done=false;
	char oneKbuff[1024];
	int  readFileSize,oneK=1024;
	//command_enum {ReplicaID, TakeThisFile, TakeRestartFile, TakeSampleData, TakeMoveEnergyData, TakeSimulationParameters, TakeCoordinateData, TakeTCS, TakeJID, NextNonInteracting, Exit, Snapshot, InvalidCommand};
	while(!done){
		fprintf(stderr,"trying to get a command\n"); //##DEBUG
		command=readCommand(sockfd);
		fprintf(stderr,"received a command: %hhu\n",command); fflush(stderr);//##DEBUG
		switch(command){
			case ReplicaID:
				fprintf(stderr,"received RepicaID command\n"); //##DEBUG
				receiveReplicaID(sockfd,&ID);
				new_replica_number=ID.replica_number;
				sprintf(IDbaseName,"%sw%d.%u",ID.title,ID.replica_number,ID.sequence_number);
				sprintf(IDenergyFileName,"%s.energy.nni%d",IDbaseName,nni);
				sprintf(IDforceFileName,"%s.force.nni%d",IDbaseName,nni);
				//sprintf(IDsetupFileName,"%s.setup",IDbaseName);
				sprintf(IDbaseName2,"%sw%d.%u",ID.title,ID.replica_number,ID.sequence_number-1);
				sprintf(IDrstFileName,"%s.rst",IDbaseName2); //CN removed a tailing .nni%d -- don't need .nni there for .rst he thinks -- also, even if we did, the , nni); was not there so it was a random int
				printf("%s\n",IDbaseName);
				//fprintf(stderr,"IDbaseName=%s\n",IDbaseName); //##DEBUG
				//fprintf(stderr,"IDbaseName2=%s\n",IDbaseName2); //##DEBUG
				//fprintf(stderr,"IDrstFileName=%s\n",IDrstFileName); //##DEBUG
				break;
			case TakeThisFile: //uncompressed, untested
				int total,count,ofd,fnSize;
				fprintf(stderr,"received TakeThisFile command\n"); //##DEBUG
				read4K(sockfd,&total,INTSIZE);
				count=total<oneK?total:oneK;
				read4K(sockfd,oneKbuff,count);
				fnSize=strlen(oneKbuff)+1;
				ofd=open(oneKbuff,O_WRONLY|O_CREAT|O_TRUNC,0644);
				if(ofd==-1){
					fprintf(stderr,"cannot open %s\n",oneKbuff); 
					exit(1);
				}
				write(ofd,oneKbuff+fnSize,count-fnSize);
				if(total-=count>0){
					char tbuf[total];
					read4K(sockfd,tbuf,total);
					write(ofd,tbuf,total);
				}
				close(ofd);
				break;
			case TakeRestartFile: //compressed
				fprintf(stderr,"received TakeRestartFile command\n"); //##DEBUG
				new_replica_number=-1;
				read4K(sockfd,&readFileSize,INTSIZE);
				if(COMPRESS_RESTART){
					receiveFile(sockfd,IDrstFileName,readFileSize);
				}else{
					receiveFileUncompressed(sockfd,IDrstFileName,readFileSize);
				}
				fprintf(stderr,"wrote restart file to file [%s]\n",IDrstFileName); //##DEBUG
				break;
			//case TakeSampleData:
				//break;
			case TakeMoveEnergyData: //uncompressed
				read4K(sockfd,&readFileSize,INTSIZE);
				receiveFileUncompressed(sockfd,IDenergyFileName,readFileSize);
				break;
			case TakeSimulationParameters:
				read4K(sockfd,&readFileSize,INTSIZE);
				//receiveFileUncompressed(sockfd,IDsetupFileName,readFileSize);
				receiveParCHARMM(sockfd,&ID,readFileSize);
				done=true;
				break;
			//case TakeCoordinateData:
				//break;
			default:
				fprintf(stderr,"received an unexpected command\n"); //##DEBUG
				done=true;
				break;
		}
	}
	fprintf(stderr,"done is true, exiting loop\n"); //##DEBUG
/*
//CN removed this as it is unexpected
	if( (new_replica_number!=old_replica_number) && (old_replica_number!=-1) ){
		fprintf(stderr,"deleting old RST file [%s]\n",rstFileName); //##DEBUG
		int e=unlink(rstFileName);   // if we are receiving a new RST file then delete old one
		if(e!=0) fprintf(stderr,"delete failed with errno %d\n",errno); //##DEBUG
	}     
*/
	close(sockfd);
	fprintf(stderr,"socket closed; helper finished\n\n\n"); //##DEBUG
	return 0;
}

void sendJID(int sockfd, float jid){
	//send jid as a float since there are already routines to handle floats in the server
	int sz1,sz2;
	unsigned int jidsize;
	sz1=KEY_SIZE+COMMAND_SIZE;
	jidsize=sizeof(jid);
	sz2=sz1+INTSIZE+jidsize;
	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=TakeJID;
	memcpy(cmd+sz1,&jidsize,INTSIZE);
	memcpy(cmd+sz1+INTSIZE,&jid,jidsize);
	int n=write(sockfd,cmd,sz2);
	fprintf(stderr,"sz2=%d wrote JID %d byte to socket %d\n",sz2,n,sockfd); //##DEBUG
}

void sendTCS(int sockfd, float tcs){
	//send tcs as a float since there are already routines to handle floats in the server
	int sz1,sz2;
	unsigned int tcssize;
	sz1=KEY_SIZE+COMMAND_SIZE;
	tcssize=sizeof(tcs);
	sz2=sz1+INTSIZE+tcssize;
	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=TakeTCS;
	memcpy(cmd+sz1,&tcssize,INTSIZE);
	memcpy(cmd+sz1+INTSIZE,&tcs,tcssize);
	int n=write(sockfd,cmd,sz2);
	fprintf(stderr,"sz2=%d wrote TCS %d byte to socket %d\n",sz2,n,sockfd); //##DEBUG
}


void sendReplicaID(int sockfd, struct ID_struct ID){
	int sz1,sz2;
	sz1=KEY_SIZE+COMMAND_SIZE;
	sz2=sz1+IDSIZE;
	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=ReplicaID;
	memcpy(cmd+sz1,&ID,IDSIZE);
	int n=write(sockfd,cmd,sz2);
	fprintf(stderr,"sz2=%d wrote ID %d byte to socket %d\n",sz2,n,sockfd); //##DEBUG
}


void receiveReplicaID(int sockfd, struct ID_struct *ID){
	read4K(sockfd,ID,IDSIZE);
}

void send_NextNonInteracting(int sockfd){
	unsigned int send_size;
	unsigned char buff[KEY_SIZE+COMMAND_SIZE];

	memcpy(buff,COMMAND_KEY,KEY_SIZE);
	buff[COMMAND_LOCATION]=NextNonInteracting;
	send_size=KEY_SIZE+COMMAND_SIZE;
	if(write(sockfd,buff,send_size)!=send_size){
		fprintf(stderr,"Error: cannot send the identifier for NextNonInteracting\n");
		//Do not exit as the server will take care of erroring out
	}
}


void sendCrdFile(int sockfd, char *filename, enum command_enum command){
	ifstream inp(filename);
	if(!inp.is_open()) return;  // line added by T. Rodinger

	string s;
	while(inp>>s){
		if(s[0]=='*')
			getline(inp,s,'\n');
		else
			break;
	}
	int natom((int)string_double(s));
	float xyz[natom*3];
	int n=0; 
	while(inp>>s>>s>>s>>s){
		inp>>xyz[n++]>>xyz[n++]>>xyz[n++];
		inp>>s>>s>>s;
	}

	unsigned int fileSize=0;
	fileSize=natom*3*sizeof(float);

	int sz1,sz1a,sz2;
	sz1=KEY_SIZE+COMMAND_SIZE;
	sz1a=INTSIZE; //fileSize
	sz2=sz1+sz1a;

	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=command;
	memcpy(cmd+sz1,&fileSize,sz1a);

	//  if(fileSize==0) while(1) sleep(10); 

	write(sockfd,cmd,sz2);
	write(sockfd,xyz,fileSize);

	inp.close();
}


int sendBinFile(int sockfd, char *filename, enum command_enum command){
	unsigned int fileSize=0;
	FILE *fp;
	float value;
	int line=0;
	char buff[100]; // added by T. Rodinger
	int i;

	//  if((fp=fopen(filename,"r"))==NULL){ fprintf(stderr,"cannot open %s\n",filename); exit(1); }
	if((fp=fopen(filename,"r"))==NULL) return 1; // the original line above was modified by T. Rodinger

	char *substring;
	i=0;
	// while loop and contents added by T. Rodinger
	while(!feof(fp)){                      
		if(fgets(buff, 99, fp)==NULL) continue;
		buff[strlen(buff)-1]=0;
		substring=strtok(buff," ");
		while(substring!=NULL){
			substring=strtok(NULL, " ");
			i++;
		}
	}
	//  fprintf(stderr,"%d data items found in the force file\n",i); //##DEBUG

	fileSize=i*sizeof(float);
	rewind(fp);
	float allV[i];
	char *string_end;
	i=0;
	// while loop and contents added by T. Rodinger
	while(!feof(fp)){
		if(fgets(buff, 99, fp)==NULL) continue;
		buff[strlen(buff)-1]=0;
		//     fprintf(stderr,"Got line: [%s]\n",buff); //##DEBUG
		substring=strtok(buff," ");
		while(substring!=NULL){
			//        fprintf(stderr,"Got substring: [%s]\n",substring); //##DEBUG
			allV[i]=strtod(substring,&string_end);
			if(string_end==substring) allV[i]=1.0e20;
			substring=strtok(NULL, " ");
			i++;
		}
	}

	int sz1,sz1a,sz2;
	sz1=KEY_SIZE+COMMAND_SIZE;
	sz1a=INTSIZE; //fileSize
	sz2=sz1+sz1a;

	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=command;
	memcpy(cmd+sz1,&fileSize,sz1a);

	//  if(fileSize==0) while(1) sleep(10); 

	write(sockfd,cmd,sz2);
	write(sockfd,allV,fileSize);

	fclose(fp);
	return 0;
}


void sendFile(int sockfd, char *filename, enum command_enum command, bool compress){
	int fd;
	unsigned int fileSize=0;

	fd=open(filename,O_RDONLY);
	if(fd==-1){
		fprintf(stderr,"cannot open %s\n",filename);
		exit(1);
	}

	int fullSize=lseek(fd,0,SEEK_END);
	lseek(fd,0,SEEK_SET);
	char buffer[fullSize];
	if(compress){
		fileSize=do_compress2(fd,buffer);
	}else{
		fileSize=fullSize;
		read4K(fd,buffer,fileSize);
	}

	int sz1,sz1a,sz2,sz2a=0;
	sz1=KEY_SIZE+COMMAND_SIZE;
	sz1a=sizeof(unsigned int); //fileSize
	sz2=sz1+sz1a;

	if(command==TakeThisFile){
		sz2a=strlen(filename);
		sz2+=(sz2a+1);
	}

	char cmd[sz2];
	memcpy(cmd,COMMAND_KEY,KEY_SIZE);
	cmd[COMMAND_LOCATION]=command;
	//memcpy(cmd+sz1,&fileSize,sz1a);

	if(command==TakeThisFile){
		unsigned int tmp=fileSize;
		tmp+=(sz2a+1);
		memcpy(cmd+sz1,&tmp,sz1a);
		memcpy(cmd+sz1+sz1a,filename,sz2a);
		char c=0;
		memcpy(cmd+sz2-1,&c,1);
		write(sockfd,cmd,sz2);
	}else{
		memcpy(cmd+sz1,&fileSize,sz1a);
		write(sockfd,cmd,sz2);
	}

	//  if(fileSize==0) while(1) sleep(10); 

	write(sockfd,buffer,fileSize);
	int fd1; //##DEBUG
	fd1=open("RESTART",O_WRONLY|O_CREAT|O_TRUNC,0644); //##DEBUG
	write(fd1,buffer,fileSize); //##DEBUG
	close(fd1); //##DEBUG
	fprintf(stderr,"********************** WROTE RESTART FILE **********************\n"); //##DEBUG

	close(fd);
}


enum command_enum readCommand(int sockfd){
	char buff[KEY_SIZE+COMMAND_SIZE];
	read4K(sockfd,buff,KEY_SIZE+COMMAND_SIZE);
	if( strncmp(buff+KEY_LOCATION,COMMAND_KEY,KEY_SIZE)!=0 ){
		fprintf(stderr,"the key did not match\n"); //##DEBUG
		exit(1);
	}

	return (enum command_enum)buff[COMMAND_LOCATION];
} 


void receiveFile(int sockfd, char *fileName, int fileSize){
	char buffer[fileSize];
	int fd;
	fd=open(fileName,O_WRONLY|O_CREAT|O_TRUNC,0644);
	if(fd==-1){
		fprintf(stderr,"cannot open %s\n",fileName);
		exit(1);
	}
	read4K(sockfd,buffer,fileSize);
	fprintf(stderr,"received file, size is: %d\n",fileSize); //##DEBUG
	do_uncompress2(fileSize,buffer,fd);
	close(fd);
}


void receiveFileUncompressed(int sockfd, char *fileName, int fileSize){
	char buffer[fileSize];
	int fd;
	fd=open(fileName,O_WRONLY|O_CREAT|O_TRUNC,0644);
	if(fd==-1){
		fprintf(stderr,"cannot open %s\n",fileName);
		exit(1);
	}
	read4K(sockfd,buffer,fileSize);
	fprintf(stderr,"received file, size is: %d\n",fileSize); //##DEBUG
	write(fd,buffer,fileSize);
	fprintf(stderr,"file written to disk\n"); //##DEBUG
	close(fd);
}


char *stringCat(char *d, const char *s){
	int len=strlen(s);
	for(int i=0; i<len; ++i)
	*(d++)=*(s++);
	return d;
}


void receiveParCHARMM(int sockfd, struct ID_struct *ID, int fileSize){
	char base1[20];
	char base2[60];
	sprintf(base1,"%sw%d.",ID->title,ID->replica_number);
	sprintf(base2,"%s%u",base1,ID->sequence_number);
	char fileName[]="setup";
	char buffer[fileSize];
	int fd;
	fd=open(fileName,O_WRONLY|O_CREAT|O_TRUNC,0644);
	if(fd==-1){
		fprintf(stderr,"cannot open %s\n",fileName);
		exit(1);
	}
	read4K(sockfd,buffer,fileSize);

	//fprintf(stderr,"start*****************************************\n"); //##DEBUG
	//buffer[fileSize]=0; //##DEBUG
	//fprintf(stderr,"%s\n",buffer); //##DEBUG
	//fprintf(stderr,"end*******************************************\n"); //##DEBUG

	char buff2[2048];
	char *c;
	c=buff2;
	c=stringCat(c,"set ");
	for(int i=0; i<fileSize; ++i){
		*(c++)=buffer[i];
		if(buffer[i]=='\n'){
			c=stringCat(c,"set ");
		}
	}

	if(ID->sequence_number==0){
		c=stringCat(c,"iob -1\n");
	}else{
		char num[40];
		sprintf(num,"%d",ID->sequence_number-1);
		c=stringCat(c,"iob \"");
		c=stringCat(c,base1);
		c=stringCat(c,num);
		c=stringCat(c,"\"\n");
	}
	c=stringCat(c,"set job \"");
	c=stringCat(c,base2);
	c=stringCat(c,"\"\n");
	*c=0;

	int len=strlen(buff2);
	if(len>2040) fprintf(stderr,"Warning: setup file may be too big.\n"); //##DEBUG
	write(fd,buff2,strlen(buff2));

	fprintf(stderr,"PARAMETER FILE:\n %s\n",buff2); //##DEBUG

	close(fd);
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

#define CHUNK 4096

unsigned int do_compress2(int ifd, char *ocp){
	z_stream z;
	char input_buffer[CHUNK];
	char output_buffer[CHUNK];
	z.avail_in = 0;
	z.next_out = (Bytef*)output_buffer;
	z.avail_out = CHUNK;
	z.zalloc = 0;
	z.zfree = 0;
	z.opaque = 0;
	deflateInit(&z,5);
	int i,count;
	unsigned int total=0;
	for(;;){
		if(z.avail_in == 0){
			z.next_in = (Bytef*)input_buffer;
			z.avail_in = read(ifd,input_buffer,sizeof(input_buffer));
		}
		if(z.avail_in == 0) break;
		deflate( &z, Z_SYNC_FLUSH );
		count = CHUNK - z.avail_out;
		total+=count;
		i=0;
		while(count--){
			*(ocp++)=output_buffer[i++];
		}
		z.next_out = (Bytef*)output_buffer;
		z.avail_out = CHUNK;
	}
	deflateEnd(&z);
	return total;
}


void do_uncompress2(int total, char *icp, int ofd){
	z_stream z;
	char input_buffer[CHUNK];
	char output_buffer[CHUNK];
	z.avail_in = 0;
	z.next_out = (Bytef*)output_buffer;
	z.avail_out = CHUNK;
	z.zalloc = 0;
	z.zfree = 0;
	z.opaque = 0;
	inflateInit(&z);
	int i,count,isz;
	isz=sizeof(input_buffer);
	for(;;){
		if(z.avail_in == 0){
			z.next_in = (Bytef*)input_buffer;
			if(total>isz){
				z.avail_in=isz;
				count=isz;
				i=0;
				while(count--) input_buffer[i++]=*(icp++);
				total-=isz;
			}else{
				z.avail_in=total;
				i=0;
				while(total--) input_buffer[i++]=*(icp++);
			}
		}
		inflate( &z, Z_SYNC_FLUSH );
		count = CHUNK - z.avail_out;
		if(count) write(ofd,output_buffer,count);
		if(total<0 && z.avail_in==0) break;
		z.next_out = (Bytef*)output_buffer;
		z.avail_out = CHUNK;
	}
	inflateEnd(&z);
}

