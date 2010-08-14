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


using namespace std;

int main(int argc, char *argv[])
{

  int sockfd; 
  struct sockaddr_in their_addr; // connector's address information 

  if (argc != 4) {
    fprintf(stderr,"Usage: %s IP-address port EXIT | SNAPSHOT\n",argv[0]);
    exit(1);
  }

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

  unsigned int protocol_version=PROTOCOL_VERSION;
  write(sockfd,&protocol_version,PROTOCOL_VERSION_SIZE); // first thing we do is send the protocol version we are using

  if(strcmp(argv[3],"EXIT")==0){
    int sz=KEY_SIZE+COMMAND_SIZE;
    char cmd[sz];
    memcpy(cmd,COMMAND_KEY2,KEY_SIZE);
    cmd[COMMAND_LOCATION]=Exit;
    write(sockfd,cmd,sz);
  }
  else if(strcmp(argv[3],"SNAPSHOT")==0){
    int sz=KEY_SIZE+COMMAND_SIZE;
    char cmd[sz];
    memcpy(cmd,COMMAND_KEY2,KEY_SIZE);
    cmd[COMMAND_LOCATION]=Snapshot;
    write(sockfd,cmd,sz);
    exit(1);
  }
  exit(1);
  

/************************************************************************/

  close(sockfd);
  return 0;
}
