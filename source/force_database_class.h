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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

struct header_struct
{
	unsigned int Nrecords;
	unsigned int Nligands;
	unsigned int Nforces_per_record;
	unsigned int Nenergies_per_record;
	unsigned int NadditionalColumns_per_record;
// NadditionalColumns_per_record represents what is refered to as Nadditional_data in other files and the input file.
// Here, the actual Nadditional_per_record is Nforces_per_record*NadditionalColumns_per_record
};

struct record_struct
{	
	int replica_number;
	unsigned int sequence_number;
	float w;
	float generic_data[0];
};

struct database_struct{
	unsigned int Nrecords;
	unsigned int Nforces;
	unsigned int Nenergies;
	unsigned int Nligands;
	unsigned int Nadditional_data;
	unsigned int record_size;
	int replica_offset;
	struct record_struct *unsorted_record;
	struct record_struct **record;
};
#define EMPTY_DATABASE_STRUCT {0,0,0,0,0,0,0,(struct record_struct *)NULL,(struct record_struct **)NULL}

class force_database_class
{
private:
	int fd;
	struct header_struct header;
	unsigned char verbose;

	struct record_struct *record;
	float *forces;
	float *energies;
	float *additionals;
	unsigned int Nforces;
	unsigned int Nenergies;
	unsigned int Nadditionals;
	
	unsigned int record_size;
	unsigned int Nreplicas;
	unsigned int sequence_number;
	int i;

	void header_existence_check(void)
	{
		if(!header_exists)
		{
			fprintf(stderr,"Error: cannot add data to the record or write the record into the database until a header has been created\n");
			exit(1);
		}
	}

	void initialise_record(void)
	{
		if(record==NULL)
		{
			record_size=sizeof(record_struct)+(header.Nforces_per_record*header.Nligands+header.Nenergies_per_record+header.Nforces_per_record*header.NadditionalColumns_per_record)*sizeof(float);
//			printf("Record size calculated to be %u\n",record_size);
			record=(struct record_struct*)malloc(record_size);
			if(record==NULL)
			{
				fprintf(stderr,"Error: cannot allocate memory for the record\n");
				exit(1);
			}
			forces=record->generic_data;
			energies=record->generic_data+header.Nforces_per_record*header.Nligands;
												additionals=energies+header.Nenergies_per_record;
		}
		record->replica_number=-1;
		Nforces=0;
		Nenergies=0;
								Nadditionals=0;
	}
	
	void write_header(void)
	{
		lseek(fd,0,SEEK_SET);
	
		if(write(fd,&header,sizeof(header))!=sizeof(header))
		{
			perror("Error: write of the header failed\n");
			exit(1);
		}
	}
	
	void increment_number_of_records(void)
	{
		header_existence_check();
	
		header.Nrecords++;
	
		lseek(fd,0,SEEK_SET);
	
		if(write(fd,&header.Nrecords,sizeof(header.Nrecords))!=sizeof(header.Nrecords))
		{
			fprintf(stderr,"Error: cannot increment number of records in database (write failed)\n");
			exit(1);
		}
	}
	
	void print_header_information(void)
	{
		fprintf(stderr,"Header information:\n");
		fprintf(stderr,"   Number of records: %u\n",header.Nrecords);
		fprintf(stderr,"   Number of ligands: %hhu\n",header.Nligands);
		fprintf(stderr,"   Number of force samples per record (per ligand): %u\n",header.Nforces_per_record);
		fprintf(stderr,"   Number of move data elements per record: %u\n",header.Nenergies_per_record);
		if(header.NadditionalColumns_per_record!=0){
			fprintf(stderr,"   Number of additional data columns per record: %u\n",header.NadditionalColumns_per_record);
			fprintf(stderr,"     Each having %u data elements\n",header.Nforces_per_record);
		}
	}

public:
	unsigned char header_exists;

	force_database_class(const char *title, unsigned char verbose_p)
	{
		#define DATABASE_FILENAME "%s.forcedatabase"
		char temp[sizeof(DATABASE_FILENAME)];
		header_exists=0;
		record=NULL;
		Nforces=0;
		Nenergies=0;
		Nadditionals=0;
		verbose=verbose_p;

		if(strlen(title)==2)
			sprintf(temp,DATABASE_FILENAME,title);
		else
			strcpy(temp,title);
			
		if( (fd=open(temp,O_RDWR|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))==-1 )
		{
			fprintf(stderr,"Error: cannot open status file for reading/writing: %s",temp);
			exit(1);
		}
		fchmod(fd,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);

		long filesize=lseek(fd,0,SEEK_END);

		if(verbose) fprintf(stderr,"The size of the database file is: %ld\n",filesize);

		
		if(filesize>=(long)sizeof(header))
		{
			lseek(fd,0,SEEK_SET);
			if(read(fd,&header,sizeof(header))!=sizeof(header))
			{
				fprintf(stderr,"Error: read of the header failed\n");
				exit(1);
			}
			if(verbose) print_header_information();
			header_exists=1;
			initialise_record();
			
		}
		else
			header_exists=0;			
	}

	~force_database_class()
	{
		close(fd);
		if(record!=NULL) free(record);
	}

	unsigned char set_header_information(const struct database_struct *db)
	{
		if(header_exists){
			if( (db->Nligands!=header.Nligands) || (db->Nforces!=header.Nforces_per_record) || 
			    (db->Nenergies!=header.Nenergies_per_record) || (db->Nadditional_data!=header.NadditionalColumns_per_record)){
				return(0);
			}
			return(1);
		}

		header.Nrecords=0;
		header.Nligands=db->Nligands;
		header.Nforces_per_record=db->Nforces;
		header.Nenergies_per_record=db->Nenergies;
		header.NadditionalColumns_per_record=db->Nadditional_data;
		write_header();
		fprintf(stderr,"A new header has been created\n");
		if(verbose) print_header_information();
		header_exists=1;
		initialise_record();
		return(1);
	}
	
	void set_number_of_records(unsigned int Nrecords)
	{
		if(!header_exists)
		{
			fprintf(stderr,"Error: cannot set the number of records because a header doesn't exist\n");
			exit(1);
		}
		header.Nrecords=Nrecords;
		write_header();
	}

	void get_header_information(struct database_struct *db)
	{
		header_existence_check();
		db->Nrecords=header.Nrecords;
		db->Nligands=header.Nligands;
		db->Nforces=header.Nforces_per_record;
		db->Nenergies=header.Nenergies_per_record;	
		db->Nadditional_data=header.NadditionalColumns_per_record;
	}
	
	unsigned int get_number_of_records(void)
	{
		return(header.Nrecords);
	}
	
	void set_record_parameters(unsigned int replica_number, unsigned int sequence_number, float w)
	{
		header_existence_check();
		record->replica_number=replica_number;
		record->sequence_number=sequence_number;
		record->w=w;
	}

/* These functions don't appear to be used. C. Neale commented them out Oct 13 2007 in order
 * that he doesn't need to be sure they are accurate in relation to the definition of Nforces
 * If you want to use these functions, please convince yourself that they are correct
 *
	void add_force(float force)
	{
		header_existence_check();
		if(Nforces==header.Nforces_per_record)
		{
			fprintf(stderr,"Error: too many force measurements exist in this file\n");
			exit(1);		
		}
		forces[Nforces++]=force;	
	}
	
	void add_move_data(float energy)
	{
		header_existence_check();
		if(Nenergies==header.Nenergies_per_record)
		{
			fprintf(stderr,"Error: too many move data elements added to this record\n");
			exit(1);		
		}
		energies[Nenergies++]=energy;
	}

				void add_additional(float additional)
	{
		header_existence_check();
		if(Nadditionals==header.Nforces_per_record*header.NadditionalColumns_per_record){
												fprintf(stderr,"Error: too many additional data elements added to this record\n");
												exit(1);
								}
		additionals[Nadditionals++]=additional;
	}
 *
 * End of functions that C. Neale commented out
 */
	
	unsigned char write_record(void)
	{
		unsigned char r=1;
		header_existence_check();
		if(record->replica_number<0)
		{
			fprintf(stderr,"Error: attempt to write a record with no record header information\n");
			exit(1);
		}
		if(Nforces!=header.Nforces_per_record)
		{
			fprintf(stderr,"Warning: attempt to write a record with incomplete force information\n");
			r=0;
		}
		if(Nenergies!=header.Nenergies_per_record)
		{
			fprintf(stderr,"Warning: attempt to write a record with incomplete move data information\n");
			r=0;
		}
		if(Nadditionals!=header.Nforces_per_record*header.NadditionalColumns_per_record){
			fprintf(stderr,"Warning: attempt to write a record with incomplete additional data information\n");
			fprintf(stderr,"         Nadditionals=%d and header.Nforces_per_record*header.NadditionalColumns_per_record=(%d*%d)=%d\n",Nadditionals,header.Nforces_per_record,header.NadditionalColumns_per_record,header.Nforces_per_record*header.NadditionalColumns_per_record);
			r=0;
		}

		if(r)
		{
			unsigned long int record_position=(unsigned long int)sizeof(header)+(unsigned long int)header.Nrecords*record_size;
		
			lseek(fd,record_position,SEEK_SET);
		
			if(write(fd,record,record_size)!=record_size)
			{
				fprintf(stderr,"Error: write of record failed\n");
				exit(1);
			}

			increment_number_of_records();
		}
		initialise_record();
		return(r);
	}

	void add_all_forces_at_once(void *buffer)
	{
		Nforces=header.Nforces_per_record;
		memcpy(forces,buffer,Nforces*header.Nligands*sizeof(float));
	}

	void add_all_additionals_at_once(void *buffer,int n)
	{
		Nadditionals=header.Nforces_per_record*n;
		memcpy(additionals+Nadditionals,buffer,header.Nforces_per_record*sizeof(float));
		Nadditionals+=header.Nforces_per_record; //required for test in write_record()
	}

	unsigned char read_record(unsigned int record_number)
	{
		header_existence_check();
		
		if(record_number>=header.Nrecords) return(0);

		unsigned long int record_position=(unsigned long int)sizeof(header)+(unsigned long int)record_number*record_size;
		
		lseek(fd,record_position,SEEK_SET);
		
		if(read(fd,record,record_size)!=record_size)
		{
			fprintf(stderr,"Error: read of record failed\n");
			exit(1);
		}
		return(1);
	}
	
	unsigned char read_record_to_given_memory(unsigned int record_number, struct record_struct *rec)
	{
		header_existence_check();
		
		if(record_number>=header.Nrecords) return(0);

		unsigned long int record_position=(unsigned long int)sizeof(header)+(unsigned long int)record_number*record_size;
		
		lseek(fd,record_position,SEEK_SET);
		
		if(read(fd,rec,record_size)!=record_size)
		{
			fprintf(stderr,"Error: read of record failed\n");
			exit(1);
		}
		return(1);
	}
	
	void print_record(FILE *f)
	{
		unsigned int i,j;

		header_existence_check();

		fprintf(f,"Replica: %d   Sequence number: %d   w: %f\n", record->replica_number, record->sequence_number, record->w);
		for(i=0;i<header.Nforces_per_record*header.Nligands;i+=header.Nligands)
			if(header.Nligands==1)      fprintf(f,"Force: %f\n",forces[i]);
			else if(header.Nligands==2) fprintf(f,"Force: %f %f\n",forces[i],forces[i+1]);
								for(i=0;i<header.Nforces_per_record;i++){
			fprintf(f,"Additional data: ");
									for(j=0;j<header.NadditionalColumns_per_record;i++){
				fprintf(f,"%f ",additionals[i*header.Nforces_per_record+j]);
			}
			fprintf(f,"\n");
		}

		if(header.Nenergies_per_record==2) fprintf(f,"Move data: New w: %f   Energy change: %f\n",energies[0],energies[1]);
		if(header.Nenergies_per_record==4) fprintf(f,"Move data: Inc w: %f   Dec w: %f   Inc dE: %f   Dec dE: %f\n",energies[0],energies[1],energies[2],energies[3]);
		else for(i=0;i<header.Nenergies_per_record;i++) fprintf(f,"Move energy: %f\n",energies[i]);
	}

};
