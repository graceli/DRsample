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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <fcntl.h>

int main(int argc, char *argv[])
{
	char dirname[40];
	char buff[8192];
	int in_fd;
	int out_fd;
	
	if(argc!=2)
	{
		fprintf(stderr,"Usage: %s package-location\n",argv[0]);
		exit(1);
	}

	if( (in_fd=open(argv[1],O_RDONLY))==-1 )
	{
		perror("Error opening source package file");
		exit(1);
	}

	sprintf(dirname,"%.10d.XXXXXX",(int)time(NULL));
	mkdtemp(dirname);
	if(dirname[0]==0)
	{
		fprintf(stderr,"Error create a unique temporary directory\n");
		exit(1);
	}
	
	sprintf(buff,"%s/package.tar.gz",dirname);
	
	if( (out_fd=open(buff, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))==-1 )
	{
		perror("Error creating destination package file");
		exit(1);
	}
	
	while(1)
	{
		int Nread=read(in_fd,buff,8192);
		if(Nread==0) break;
		if(Nread==-1)
		{
			perror("Error reading source package file");
			exit(1);
		}
		if(write(out_fd,buff,Nread)!=Nread)
		{
			perror("Error writing destination package file");
			exit(1);
		}
	}
	close(in_fd);
	close(out_fd);
	
	printf("%s\n",dirname);

	return(0);
}
