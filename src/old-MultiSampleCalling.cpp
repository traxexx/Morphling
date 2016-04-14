#include <iostream>
#include <fstream>
#include <math.h> // round
#include <sstream>
#include <string>
#include <unistd.h> // getpid

#include "MultiSampleCalling.h"
#include "Globals.h"
#include "MeiSeq.h"
#include "AsbSites.h"
#include "RgtSites.h"

using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::getline;

/* multi-sample calling steps:
 	1. read through all vcfs, get candidate site list with q>=q-thred
 	2. read through bam files, get assembly info. output filtering stat & vcf
 	3. for pass site, read through each vcf, get GL from bam on those un-genotyped sites. ---> ReGenotype()
 	4. make final calling based on posterior-likelihood
*/

// main function for doing multi-calling
void MultiSampleCalling( int argc, char * argv[] )
{
// set options
	string arg_str, dummy_str;
	setMainOptionStr( arg_str, dummy_str );
	Options * ptrMainOptions = new Options( argc, argv, arg_str, dummy_str );

// set globals & path
	SetGlobals(ptrMainOptions);
	vector<string> mei_seq_file_names;
	setMeiSingleFileNames(  ptrMainOptions->ArgMap["MElist"], mei_seq_file_names );
	vector< vector<string> > MEnames;
	MEnames.resize(NMEI);
	vector< vector<string> > MEseqs;
	MEseqs.resize(NMEI);
	for( int m=0; m<NMEI; m++ )
		SetMeiSeqAndNames( mei_seq_file_names[m], MEnames[m], MEseqs[m]);

// generate sample list (4 cols)
	MultiSampleInfo msinfo;
	LoadDiscoverSampleListInfo( ptrMainOptions->ArgMap["SampleList"], msinfo );

// generate site list
	map<int, map<string, map<int, int> > > candidate_sites; // mtype-> chr -> round(breakp/WIN) : breakp
	vector<string> chr_list;
	setChrList( ptrMainOptions, chr_list );
	for( int m=0; m<NMEI; m++ ) {
		Sites single_type_sites( m, chr_list );
		single_type_sites.SetSiteList( msinfo.discovery_dirs );
		single_type_sites.ExportPositions( candidate_sites[m] );
	}

// assembly
	AsbSites assembled_sites( MEnames, MEseqs );
	assembled_sites.LoadCandidateSitesFromSites( candidate_sites );
	for( int s=0; s<msinfo.bam_list.size(); s++ ) {
		assembled_sites.AddEvidenceFromSingleBam( msinfo.bam_list[s], msinfo.depth[s] );
	}
	assembled_sites.Assembly();

	string hard_filter_vcf_name = ptrMainOptions->ArgMap["Out"] + ".hardfilter.site.vcf";
	assembled_sites.PrintToVcf( hard_filter_vcf_name );

/****** assembly part finish ***/

	RgtSites regt_sites( MEnames, MEseqs );
	regt_sites.LoadSitesFromVcf( hard_filter_vcf_name );
	for( int s=0; s<msinfo.sample_names.size(); s++ ) {
		regt_sites.ReGenotypeFromSingleBam( msinfo.sample_names[s], msinfo.bam_list[s] );
	}
	regt_sites.PrintVcf( ptrMainOptions->ArgMap["Out"] );

/*
	SetGenotypeParameters(ptrMainOptions);
	std::vector<std::vector<string> > re_gt_list;
	LoadRegenotypeListFromScratch(ptrMainOptions->ArgMap["WorkDir"], re_gt_list);
	ReGenotypeOnebyOne(ptrMainOptions, all_sites, re_gt_list);
	ConsensusVcf FinalVcf( all_sites, re_gt_list );
	FinalVcf.Print(ptrMainOptions->ArgMap["Out"]);

// finish	
	cout << "Morphling finished with no error reported. Check final output at: " << ptrMainOptions->ArgMap["Out"].c_str() << endl;
*/
}

void SetGlobals( Options* ptrMainOptions )
{
	WIN = 600;
	MIN_CLIP = 20;
	NMEI = 3;
}

void LoadDiscoverSampleListInfo(string & sample_list_name, MultiSampleInfo & msinfo)
{
	string line;
	ifstream sample_list_file;
	sample_list_file.open( sample_list_name.c_str() );
	if ( !sample_list_file.is_open() ) {
		cerr << "ERROR: Unable to open sample list: " << sample_list_name << endl;
		exit(1);
	}
	while( getline( sample_list_file, line ) ) {
		std::stringstream ss;
		ss << line;
		vector<string> current_info;
		string field;
		while( getline( ss, field, '\t') )
			current_info.push_back(field);
		if ( current_info.size() != 4 ) {
			cerr << "ERROR: [LoadSampleList()] current line:\n " << line << "is not a regular sample list line!\n sample list should be 3 fields: sample-name bam discover-out-directory." << endl;
			exit(1);
		}
		msinfo.sample_names.push_back( current_info[0] );
		msinfo.bam_list.push_back( current_info[1] );
		msinfo.discovery_dirs.push_back( current_info[2] );
		int dp;
		if ( str_is_float( current_info[3] ) )
			dp = stof( current_info[3] );
		else {
			cerr << "ERROR: sample list line doesn't have a valid depth value:\n " << line << endl;
			exit(2);
		}
		msinfo.depth.push_back( dp );
	}
	sample_list_file.close();
// add slash to discover directory
	for( vector<string>::iterator t=msinfo.discovery_dirs.begin(); t!= msinfo.discovery_dirs.end(); t++ ) {
		int last = t->length() - 1;
		if ( (*t)[last] != '/' )
			(*t) += '/';
	}
}

void setMainOptionStr( string & arg_str, string & dummy_str )
{
	string ref_path = GetMorphBasePath() + "refs/";

	arg_str = "-Win=600;-Step=100;-Chr=-1;";
	arg_str += "-GenomeFasta=/net/wonderland/home/saichen/reference/archive/hs37d5.fa;";
	arg_str += "-MElist=" + ref_path + "MobileElement.list;";
	arg_str += "-SampleList= ;-Out= ";
	dummy_str = "";
}


/*
	decide from option: single chr or all chr?
	read from .fai
	return error if needed
*/
void setChrList( Options* ptrMainOptions, vector<string> & chr_list )
{
	std::map<string, bool> chr_name_tmp_map;
	for( int i=0; i<22; i++ ) {
		chr_name_tmp_map[ std::to_string(i) ] = 1;
		string s = "chr" + std::to_string(i);
		chr_name_tmp_map[ s ] = 1;
	}	
	chr_name_tmp_map["X"] = 1;
	chr_name_tmp_map["Y"] = 1;
	chr_name_tmp_map["chrX"] = 1;
	chr_name_tmp_map["chrY"] = 1;

	string fai_name = ptrMainOptions->ArgMap["GenomeFasta"] + ".fai";
	std::ifstream fai;
	fai.open( fai_name.c_str() );
	if ( !fai.is_open() ) {
		cerr << "ERROR: [setChrList] Unable to open fasta index: " << fai_name << endl;
		exit(10);
	}
	string line;
	if ( ptrMainOptions->ArgMap["Chr"].compare("-1") == 0 ) { // all chr
		while( getline( fai, line ) ) {
			stringstream ss;
			ss << line;
			string name;
			getline( ss, name, '\t' );
			if ( chr_name_tmp_map.find( name ) != chr_name_tmp_map.end() )
				chr_list.push_back( name );
		}
		if ( chr_list.empty() ) {
			cerr << "ERROR: [MultiSampleCalling->setChrList] cannot get correct chr name from fasta index: " << fai_name << endl;
			exit(12);
		}
	}
	else { // specific chr
		string ref_name = ptrMainOptions->ArgMap["Chr"];
		while( getline( fai, line ) ) {
			stringstream ss;
			ss << line;
			string name;
			getline( ss, name, '\t' );
			if ( ref_name.compare(name) == 0 )
				chr_list.push_back( name );
		}
		if ( chr_list.empty() ) {
			cerr << "ERROR: Cannot match chromosome " << ref_name << " in -Chr option to fasta index!\n" << endl;
			exit(11);
		}
	}
	fai.close();
}

void setMeiSingleFileNames( string & filename, vector<string> & mei_seq_file_names )
{
	mei_seq_file_names.resize(NMEI);
	ifstream file;
	file.open( filename.c_str() );
	if (!file.is_open()) {
		cerr << "ERROR: cannot open " << filename << ". Check if -MElist option is correcly set!\n";
		exit(1);
	}
	string line;
	while(getline(file, line)) {
		stringstream ss;
		ss << line;
		string name, fname;
		getline(ss, name, '\t');
		getline(ss, fname, '\t');
		int m = GetMeiIndexFromName( name );
		mei_seq_file_names[m] = fname;
		if (m==-1) {
			cerr << "ERROR: " << name << " is not recognized in -MElist!\n" << endl;
			exit(1);
		}
	}
	file.close();
	for( int m=0; m<NMEI; m++ ) {
		if ( mei_seq_file_names[m].length()==0 ) {
			string name = GetMeiNameFromIndex(m);
			cerr << "ERROR: Missing " << name << " seq file in -MElist!\n";
			exit(1);
		}
	}
}


// get bin path first
// then remove "bin/""
string GetMorphBasePath()
{
	string path;
	pid_t pid = getpid();
    char buf[20] = {0};
    sprintf(buf,"%d",pid);
    std::string _link = "/proc/";
    _link.append( buf );
    _link.append( "/exe");
    char proc[512];
    int ch = readlink(_link.c_str(),proc,512);
    if (ch != -1) {
        proc[ch] = 0;
        path = proc;
        std::string::size_type t = path.find_last_of("/");
        path = path.substr(0,t);
    }
    if ( path.length() < 3 ) {
    	cerr << "ERROR: Unable to get Morphling bin path!\n" << endl;
    	exit(2);
    }
    if ( path.length() == 3 )
    	path = "";
    else {
    	path = path.substr( 0, path.length() - 3 ) + '/';
    }
    return path;
}


// int is also counted as valid float number
bool str_is_float( string & str )
{
	if ( str[str.size()-1] == '.' )
		return 0;
	bool success = 1;
	bool dot = 0;
	for( int i=0; i<str.size(); i++ ) {
		if ( str[i] == '.' ) {
			if ( !dot ) {
				dot = 1;
				continue;
			}
			else {
				success = 0;
				break;
			}
		}
		if ( !isdigit( str[i] ) ) {
			success = 0;
			break;
		}
	}
	return success;
}



































