#ifndef MULTISAMPLECALLING_H
#define MULTISAMPLECALLING_H

#include "Options.h"
#include "Sites.h"
#include "Likelihood.h"
#include "RefStat.h"
#include "singleSample.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

void MultiSampleCalling( int argc, char * argv[] );

void SetGlobals( Options* ptrMainOptions );

void initializeLogFile( Options* ptrMainOptions );

void LoadSampleList(string & sample_list_name, vector<SingleSampleInfo> & msinfo );

// assembly part
string generateSiteVcf( Options * ptrMainOptions, vector<SingleSampleInfo> & msinfo, 
	vector< vector<string> > & MEnames, vector< vector<string> > & MEseqs );

void setMainOptionStr( string & arg_str, string & dummy_str );

void setChrList( Options* ptrMainOptions, vector<string> & chr_list );

void setMeiFileNames( vector<string> & ms_name );

string GetMorphBasePath();

bool LoadSitesFromVcf( map<string, vector<GtRec> > & can_sites, string & vcf_name, vector< vector<string> > & MEnames );

void loadRefFileMap( string & filename, map<string, string> & fv );

void printGtVcf( vector<SingleSampleInfo> & msinfo, string & hard_filter_vcf_name, 
	map<string, vector<GtRec> > & can_sites, string & out_name );

void estimateHomPower();

#endif