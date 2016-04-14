#include "MeiSeq.h"
#include "seqUtility.h"
#include "morphError.h"
#include <fstream>


void SetMeiSeqAndNames( string & filename, vector<string> & me_names, vector<string> & me_seqs )
{
	// count # records
	std::ifstream file;
	file.open( filename.c_str() );
	if ( !file.is_open() )
		morphErrorFile(filename);
	string line;
	int n_record = 0;
	while( getline(file, line) ) {
		if ( line.size()==0 )
			continue;
		if ( line[0] == '>' )
			n_record++;
	}
	file.close();
	me_names.resize( n_record );
	me_seqs.resize( n_record );

	// add
	file.open( filename.c_str() );
	if ( !file.is_open() )
		morphErrorFile(filename);
	string seq;
	string name;
	int index = -1;
	while( getline(file, line) ) {
		if ( line.size()==0 )
			continue;
		if ( line[0] == '>' ) { // name line
			if (index>=0) {
				me_names[index] = name;
				me_seqs[index] = seq;
			}
			index++;
			seq.clear();
			int space_index = line.length();
			for(int i=1; i<line.length(); i++) {
				if ( line[i] == ' ' || line[i] == '\t' ) {
					space_index = i;
					break;
				}
			}
			name = line.substr( 1, space_index - 1 );
		}
		else // seq line
			seq += line;
	}
	// add last record
	me_names[index] = name;
	me_seqs[index] = seq;
	file.close();
}

string GetMeiNameFromIndex( int m )
{
	string name;
	if (m==0)
		name = "ALU";
	else if (m==1)
		name = "L1";
	else if (m==2)
		name = "SVA";
	else
		name = "UNKNOWN";
	return name;
}

int GetMeiIndexFromName( string & name )
{
	int m = -1;
	for( int i=0; i<name.size(); i++ )
		name[i] = std::toupper( name[i] );
	if ( name.compare("ALU") == 0 )
		m = 0;
	else if ( name.compare("L1") == 0 )
		m = 1;
	else if ( name.compare("SVA") == 0 )
		m = 2;
	else {
		string str = "Unrecognized mei name: " + name;
		morphError(str, 3);
	}
	return m;
}
