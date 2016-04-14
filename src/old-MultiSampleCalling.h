#ifndef MULTISAMPLECALLING_H
#define MULTISAMPLECALLING_H

#include "Options.h"
#include "Sites.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

struct MultiSampleInfo {
	int n;
	vector<string> sample_names;
	vector<string> bam_list;
	vector<string> discovery_dirs;
	vector<int> depth;
};


void MultiSampleCalling( int argc, char * argv[] );

void SetGlobals( Options* ptrMainOptions );

void LoadDiscoverSampleListInfo(string & sample_list_name, MultiSampleInfo & msinfo);

void setMainOptionStr( string & arg_str, string & dummy_str );

void setChrList( Options* ptrMainOptions, vector<string> & chr_list );

void setMeiSingleFileNames( string & filename, vector<string> & mei_seq_file_names );

string GetMorphBasePath();

bool str_is_float( string & str );

#endif