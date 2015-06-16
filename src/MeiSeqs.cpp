#include "MeiSeqs.h"
#include "Utilities.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <sstream>

using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;
using std::stringstream;

void loadMEsequence( string & me_list_name, vector< vector<string> > & MEseqs, vector< vector<string> > & MEnames )
{
// read me list
	ifstream me_list;
	me_list.open( me_list_name.c_str() );
	CheckInputFileStatus( me_list, me_list_name.c_str() );
	string line;
	vector<string> me_files;
	me_files.resize(3);
	while( getline( me_list, line ) ) {
		stringstream ss;
		ss << line;
		string me_name;
		getline( ss, me_name, '\t' );
		if ( me_name.compare("ALU") == 0 )
			getline(ss, me_files[0], '\t');
		else if ( me_name.compare("L1") == 0 )
			getline(ss, me_files[1], '\t');
		else if ( me_name.compare("SVA") == 0 )
			getline(ss, me_files[2], '\t');
		else {
			cerr << "ERROR: invalid ME name in MElist: " << me_name << endl;
			exit(1);
		}
	}
	me_list.close();
	if ( me_files[0].empty() || me_files[1].empty() || me_files[2].empty() ) {
		cerr << "ERROR: one or more ME fasta missing!" << endl;
		exit(1);
	}
// read me seqs & names
	MEseqs.resize(3);
	MEnames.resize(3);	
// set size first
	for( int i=0; i<3; i++ ) {
		ifstream infa;
		infa.open( me_files[i].c_str() );
		CheckInputFileStatus( infa, me_files[i].c_str() );
		string line;
		int rec_count = 0;
		while( getline( infa, line ) ) {
			if ( line[0] == '>' )
				rec_count++;
		}
		infa.close();
		if ( rec_count == 0 ) {
			cerr << "ERROR: fasta " << me_files[i] << " does not contain any fasta record!" << endl;
			exit(1);
		}
		MEseqs[i].resize( rec_count );
		MEnames[i].resize( rec_count);
	}
// load name & seq
	for( int mindex = 0; mindex<3; mindex++ ) {
		ifstream infa;
		infa.open( me_files[mindex].c_str() );
		CheckInputFileStatus( infa, me_files[mindex].c_str() );
		string line;
		int ri = 0;
		string seq; // sequence
		string name; // name
		while( getline( infa, line ) ) {
			if ( line.size() == 0 ) // skip empty lines
				continue;
			if ( line[0] == '>' ) { // name line
				if ( !seq.empty() ) {
					MEseqs[mindex][ri] = seq;
					MEnames[mindex][ri] = name;
					ri++;
				}
				seq.clear();
				int space_index = line.length();
				for( int i=1; i<(int)line.length(); i++ ) {
					if ( line[i] == ' ' || line[i] == '\t' ) {
						space_index = i;
						break;
					}
				}
				name = line.substr( 1, space_index - 1 );
			}
			else { // seq line
				seq += line;
			}
		}
		MEseqs[mindex][ri] = seq;
		MEnames[mindex][ri] = name;
		infa.close();
	}	
}

string ReverseCompSeq( string & seq ) {
	string rv;
	rv.resize( seq.size() );
	for( int i=0; i<(int)seq.size(); i++ )
		rv[i] = GetCompNt( seq[i] );
	return rv;
}

char GetCompNt( char nt)
{
	char c;
	switch( nt ) {
		case('A'):
			c = 'T'; break;
        case ('T'):
			c = 'A'; break;
		case ('C'):
			c = 'G'; break;
		case ('G'):
			c = 'C'; break;
		case ('a'):
			c = 'T'; break;
		case ('t'):
			c = 'A'; break;
		case ('c'):
			 c = 'G'; break;
		case ('g'):
			c = 'C'; break;
		default:
			c = 'N'; break; 
	}
	return c;
}

int GetMEtypeFromAlt( string & field)
{
	int index = field.length() - 2;
	int mtype;
	if ( field[index] == 'U' )
		mtype = 0;
	else if ( field[index] == '1' )
		mtype = 1;
	else if ( field[index] == 'A' )
		mtype =2;
	else {
		cerr << "ERROR: [GetMEtypeFromAlt] " << field << " is not a regular ALT from Morphling Genotype!" << endl;
		exit(1);
	}
	return mtype;
}