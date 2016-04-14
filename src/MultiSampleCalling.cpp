#include <iostream>
#include <fstream>
#include <math.h> // round
#include <sstream>
#include <string>
#include <unistd.h> // getpid

#include "MultiSampleCalling.h"
#include "morphError.h"
#include "Globals.h"
#include "Likelihood.h"
#include "MeiSeq.h"
#include "AsbSites.h"
#include "glCalc.h"
#include "qcCheck.h"
#include "VcfRecord.h"

using std::ifstream;
using std::ofstream;
using std::stringstream;
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

// initialize log file
	initializeLogFile(ptrMainOptions);
	morphMessage("\nAnalysis started");

// set globals & path
	PATH = GetMorphBasePath(); // PATH is global
	SetGlobals(ptrMainOptions);
	vector<string> mei_seq_file_names;
	setMeiFileNames( mei_seq_file_names );
//	setMeiSingleFileNames(  ptrMainOptions->ArgMap["MElist"], mei_seq_file_names );
	vector< vector<string> > MEnames;
	MEnames.resize(NMEI);
	vector< vector<string> > MEseqs;
	MEseqs.resize(NMEI);
	for( int m=0; m<NMEI; m++ )
		SetMeiSeqAndNames( mei_seq_file_names[m], MEnames[m], MEseqs[m]);

// generate sample list (5 cols)
	vector<SingleSampleInfo> msinfo;
	LoadSampleList( ptrMainOptions->ArgMap["SampleList"], msinfo );
	string msg = "Loaded sample list. Size = " + std::to_string(msinfo.size());
	morphMessageNoTime(msg);

	// generate site list and assembly sites across samples
	string pass_site_vcf_name;
	if (ptrMainOptions->ArgMap["SiteVcf"].compare(".") == 0) { // need to assembly from scratch
		morphMessage("Generating site vcf...");
		pass_site_vcf_name = generateSiteVcf( ptrMainOptions, msinfo, MEnames, MEseqs );
		morphMessage( "Finished generating site vcf" );
	}
	else // read directly
		pass_site_vcf_name = ptrMainOptions->ArgMap["SiteVcf"];

	if (ptrMainOptions->OptMap["siteOnly"]) // only generate site vcf
		return;

	/******** genotype part ****/
	// load site first
	map<string, vector<GtRec> > can_sites;
	bool load_status = LoadSitesFromVcf( can_sites, pass_site_vcf_name, MEnames );
	if (!load_status) {
		morphWarning("Cannot load from site vcf. Do not perform genotyping");
		if (!ptrMainOptions->OptMap["statOnly"]) {
			string cmd = "cp -f " + pass_site_vcf_name + " " + ptrMainOptions->ArgMap["OutVcf"];
			ExecuteCmd(cmd);
			return;
		}
	}

	// generate ref stats
	map<string, string> refstat_map;
	if ( ptrMainOptions->ArgMap["rSingle"].compare(".") == 0 && ptrMainOptions->ArgMap["rStat"].compare(".") == 0  ) {
		// generate from scratch
		bool ref_exclusion;
		if (ptrMainOptions->OptMap["noRefExclusion"])
			ref_exclusion = 0;
		else
			ref_exclusion = 1;
		morphMessageNoTime("Generating ref stats from scratch...");
		for( int s=0; s<msinfo.size(); s++ ) {
			string out_prefix = ptrMainOptions->ArgMap["OutVcf"] + ".sample." + std::to_string(s);
			GenerateRefStats( msinfo[s], MEseqs, MEnames, pass_site_vcf_name, out_prefix, ref_exclusion );
			refstat_map[msinfo[s].sample_name] = out_prefix;
		}
		morphMessageNoTime("Generated all ref stats");
	}
	else { // load from arguments
		if ( ptrMainOptions->ArgMap["rSingle"].compare(".") != 0 ) { // a single ref stat
			for(int s=0; s<msinfo.size(); s++)
				refstat_map[msinfo[s].sample_name] = ptrMainOptions->ArgMap["rSingle"];
		}
		else { // load from list
			loadRefFileMap( ptrMainOptions->ArgMap["rStat"], refstat_map );
			if (refstat_map.empty())
				morphError("No ref stats loaded. Check if -rStat option is properlly set?");
		}
	}

	// estimate homozygotes detection power from ctrl
	if (!ptrMainOptions->OptMap["noPower"])
		estimateHomPower();

	if (ptrMainOptions->OptMap["statOnly"]) // only generate reference stat
		return;

	// gt by sample
	for(int s=0; s<msinfo.size(); s++) {
		string sample_name = msinfo[s].sample_name;
		map<string, string>::iterator t = refstat_map.find(sample_name);
		if (t==refstat_map.end()) {
			string str = "Cannot find refstat of sample " + sample_name +  ". Check your -rStat option!";
			morphError(str, 30);
		}
		string ref_prefix = t->second;
		Likelihood lh( MEseqs, msinfo[s], ref_prefix );
		string msg = "Genotyping sample #" + std::to_string(s) + "...";
		morphMessage(msg);
		for( map<string, vector<GtRec> >::iterator c = can_sites.begin(); c != can_sites.end(); c++ ) {
			for( int i=0; i<c->second.size(); i++ ) {
				c->second[i].info.resize( msinfo.size() );
				lh.SetLhRec( c->second[i].info[s], c->first, c->second[i].position, c->second[i].mtype, c->second[i].subtype );
			}
		}
		msg = "Finish genotyping sample #" + std::to_string(s);
		morphMessage(msg);
	}

	// print final genotyped vcf
	morphMessageNoTime("Printing final vcf...");
	printGtVcf( msinfo, pass_site_vcf_name, can_sites, ptrMainOptions->ArgMap["OutVcf"] );

	string str = "Completed. Check final vcf at: " + ptrMainOptions->ArgMap["OutVcf"];
	morphMessage( str );
}

void SetGlobals( Options* ptrMainOptions )
{
	if ( ptrMainOptions->OptMap["noClean"] )
		CLEAN = 0;
	else
		CLEAN = 1;
	WIN = 600;
	MIN_CLIP = 20;
	NMEI = 3;
	if (ptrMainOptions->ArgMap["nNegRef"].compare(".") != 0)
		N_NEG_REF = stoi(ptrMainOptions->ArgMap["nNegRef"]);
	else
		N_NEG_REF = 6000;
}


// set LOG in globals
// put default log as outvcf.log
void initializeLogFile( Options* ptrMainOptions )
{
	string log_name;
	if (ptrMainOptions->ArgMap["Log"].compare(".")==0)
		log_name = ptrMainOptions->ArgMap["OutVcf"] + ".log";
	else
		log_name = ptrMainOptions->ArgMap["Log"];
	LOG.open(log_name.c_str());
	if (!LOG.is_open())
		morphErrorFile(log_name);
}


void LoadSampleList(string & sample_list_name, vector<SingleSampleInfo> & msinfo )
{
	string line;
	ifstream sample_list_file;
	sample_list_file.open( sample_list_name.c_str() );
	if ( !sample_list_file.is_open() ) {
		string str = "ERROR: Unable to open sample list: " + sample_list_name;
		morphError(str, 51);
	}
	while( getline( sample_list_file, line ) ) {
		std::stringstream ss;
		ss << line;
		vector<string> current_info;
		string field;
		while( getline( ss, field, '\t') )
			current_info.push_back(field);
		if ( current_info.size() != 6 ) {
			string str = "[LoadSampleList()] current line:\n" + line + "\n";
			str += "sample list should have 6 columns: sample-name bam depth read-length avr-ins-size var-avr-ins-size.";
			morphError(str, 1);
		}
		SingleSampleInfo new_info;
		new_info.sample_name = current_info[0];
		new_info.bam_name = current_info[1];
		new_info.depth = stof(current_info[2]);
		new_info.avr_read_length = stoi(current_info[3]);
		new_info.avr_ins_size = stoi(current_info[4]);
		new_info.var_avr_ins_size = stoi(current_info[5]);
		msinfo.push_back(new_info);
	}
	sample_list_file.close();
}

void setMainOptionStr( string & arg_str, string & dummy_str )
{
	arg_str = "-Win=600;-Chr=-1;-Range=.;";
	arg_str += "-GenomeFasta=/net/wonderland/home/saichen/reference/archive/hs37d5.fa;";
	arg_str += "-SampleList= ;-OutVcf= ;-Log=.;";
	arg_str += "-SiteVcf=.;-rSingle=.;-rStat=.;-nNegRef=.";
	dummy_str = "--noClean;--siteOnly;--statOnly;--noPower;--noRefExclusion";
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
		string str = "[setChrList] Unable to open fasta index: " + fai_name;
		morphError(str, 10);
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
			string str = "[MultiSampleCalling->setChrList] cannot get correct chr name from fasta index: " + fai_name;
			morphError(str, 12);
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
			string str = "Cannot match chromosome " + ref_name + " in -Chr option to fasta index";
			morphError(str, 11);
		}
	}
	fai.close();
}

void setMeiFileNames( vector<string> & ms_name )
{
	ms_name.resize(NMEI);
	ms_name[0] = PATH + "refs/Alu.fa";
	ms_name[1] = PATH + "refs/L1.fa";
	ms_name[2] = PATH + "refs/SVA.fa";
}

string generateSiteVcf( Options * ptrMainOptions, vector<SingleSampleInfo> & msinfo, 
	vector< vector<string> > & MEnames, vector< vector<string> > & MEseqs )
{
	// generate site list
	string ref_path = PATH + "refs/";
	vector<string> mei_list_name;
	mei_list_name.resize(3);
	mei_list_name[0] = ref_path + "Alu.bed";
	mei_list_name[1] = ref_path + "L1.bed";
	mei_list_name[2] = ref_path + "SVA.bed";
	map<int, map<string, map<int, int> > > candidate_sites; // mtype-> chr -> start : end
	vector<string> chr_list;
	setChrList( ptrMainOptions, chr_list );
	string single_chr = "";
	Sites preliminary_sites( mei_list_name, chr_list );
	preliminary_sites.MakePreliminarySiteList( msinfo, ptrMainOptions->ArgMap["Range"], candidate_sites );

/*
for( int i=0; i<3; i++ ) {
	std::cout << "at " << i << ", size=" << candidate_sites[i].size() << std::endl;
	for( map<int, int>::iterator t=candidate_sites[i]["17"].begin(); t !=candidate_sites[i]["17"].end(); t++ )
		std::cout << "  " << t->first << "->" << t->second << std::endl;
}
*/


// assembly
	AsbSites assembled_sites( MEnames, MEseqs );
	assembled_sites.Assembly( candidate_sites, msinfo );
	string hard_filter_vcf_name = ptrMainOptions->ArgMap["OutVcf"] + ".hardfilter.site.vcf";
	assembled_sites.PrintToVcf( hard_filter_vcf_name );
	// then generate pass site vcf
	string pass_site_vcf_name = ptrMainOptions->ArgMap["OutVcf"] + ".pass.site.vcf";
	generatePassVcf( hard_filter_vcf_name, pass_site_vcf_name );

	// return for futher step
	return pass_site_vcf_name;
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
    if ( path.length() < 3 )
    	morphError("Unable to get Morphling bin path");
    if ( path.length() == 3 )
    	path = "";
    else
    	path = path.substr( 0, path.length() - 3 );
    return path;
}

bool LoadSitesFromVcf( map<string, vector<GtRec> > & can_sites, string & vcf_name, vector< vector<string> > & MEnames )
{
	if (IsEmptyFile(vcf_name)) {
		string str = "Empty file " + vcf_name;
		morphWarning(str);
		return 0;
	}

	// initialize
	map<string, int> counts;
	std::ifstream vcf;
	vcf.open(vcf_name.c_str());
	string line;
	bool header_checked = 0;
	string chr;
	while( getline( vcf, line ) ) {
		if (!header_checked) {
			if ( IsHeaderLine(line) )
				continue;
			header_checked = 1;
		}
		VcfRecord vr( line, 0 ); // no gt
		if ( chr.compare( vr.GetChr() ) != 0) {
			chr = vr.GetChr();
			counts[chr] = 1;
			continue;
		}
		counts[chr]++;
	}
	vcf.close();
	if ( counts.empty() ) {
		string str = "  " + vcf_name + " has no valid record";
		morphWarning(str);
		return 0;
	}
	for( map<string, int>::iterator t = counts.begin(); t != counts.end(); t++ )
		can_sites[t->first].resize( t->second );

	// load	
	vcf.open(vcf_name.c_str());
	header_checked = 0;
	int idx = 0;
	chr = "";
	while( getline( vcf, line ) ) {
		if (!header_checked) {
			if ( IsHeaderLine(line) )
				continue;
			header_checked = 1;
		}
		VcfRecord vr( line, 0 ); // no gt
		if ( chr.compare( vr.GetChr() ) != 0 ) {
			idx = 0;
			chr = vr.GetChr();
		}
		can_sites[chr][idx].position = vr.GetPosition();
		can_sites[chr][idx].mtype = vr.GetMeiType();
		can_sites[chr][idx].subtype = vr.GetMeiSubtypeIndex( MEnames );
		idx++;
	}
	vcf.close();
	return 1;
}

void loadRefFileMap( string & filename, map<string, string> & fv )
{
	ifstream fl;
	fl.open(filename.c_str());
	if (!fl.is_open())
		morphErrorFile(filename);
	string line;
	while(std::getline(fl, line)) {
		std::stringstream ss;
		ss << line;
		string name;
		getline( ss, name, '\t');
		string rf;
		getline( ss, rf, '\t');
		if (name.empty() || rf.empty()) {
			string str = "[loadRefFileMap] Error at line:\n  " + line;
			morphError(str, 20);
		}
		fv[name] = rf;
	}
	fl.close();
}


void printGtVcf( vector<SingleSampleInfo> & msinfo, string & hard_filter_vcf_name,
	map<string, vector<GtRec> > & can_sites, string & out_name )
{
	std::ifstream invcf;
	invcf.open( hard_filter_vcf_name.c_str() );
	if (!invcf.is_open())
		morphErrorFile( hard_filter_vcf_name );
	std::ofstream outvcf;
	outvcf.open( out_name.c_str() );
	if (!outvcf.is_open())
		morphErrorFile( out_name );
	// read & print
	string line;
	string chr;
	vector<GtRec>::iterator psite;
	bool pass_header = 0;
	while( getline(invcf, line) ) {
		if (!pass_header) {
			if (!IsHeaderLine(line))
				morphError("No #CHROM line in vcf");
			if ( line[1] == 'C' ) {
				outvcf << line;
				for(int i=0; i<msinfo.size(); i++)
					outvcf << "\t" << msinfo[i].sample_name;
				outvcf << endl;
				pass_header = 1;
			}
			else
				outvcf << line << endl;;
			continue;
		}
		VcfRecord vr(line, 0);
		if ( chr.compare( vr.GetChr() ) != 0 ) {
			chr = vr.GetChr();
			psite = can_sites[chr].begin();
		}
		if ( vr.GetPosition() != psite->position ) {
			string str = "At " + chr + ":" + std::to_string(vr.GetPosition())  + " can_sites is " + std::to_string(psite->position);
			morphError( str, 41 );
		}
		vector<float> prior;
		bool af_set = setAFfromGL( prior, psite->info ); // ref allele frequency
		if (af_set) { // only print the one with gl info
			vector<int> gt;
			gt.resize( psite->info.size() );
			vector<float> log10_prior;
			log10_prior.resize(3);
			for(int i=0; i<3; i++) // convert to prior
				log10_prior[i] = log10(prior[i]);
			for( int i=0; i<gt.size(); i++ )
				gt[i] = getGtFromGl( log10_prior, psite->info[i].gl );
			// set af & ac
			vr.SetAFinfo( prior,gt );

			// print
			vr.PrintNoEnding( outvcf );
			outvcf << "\tGT:PL:AD:OD";
			for( int i=0; i<psite->info.size(); i++ ) {
				outvcf << "\t";
				PrintLhRecToVcf( outvcf, gt[i], psite->info[i] );
			}
			outvcf << std::endl;
		}
		psite++;
	}
	if (!pass_header)
		morphError("No #CHROM line in vcf");
	invcf.close();
	outvcf.close();
}


void estimateHomPower()
{

}




















