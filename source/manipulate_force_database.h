/*
 * This header file contains routines that are used by analyze_force_database and joinFDBs
 */

#include "force_database_class.h"
#include "read_input_script_file.h"


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

