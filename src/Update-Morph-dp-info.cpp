#include <sys/types.h> // getpid
#include <unistd.h> // getpid
#include <sstream>
#include <string>
#include <iostream>
#include "Utilities.h"
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;

// add QUAL to QC.info

void ResetQC( string & dname );
int SetQualThreshold( string & rec, string & qc_dir, string & bed_name, string & ex_name );
int getQualThredFromFile( string & qlog_name );
string GetMeiNameFromIndex( int mt );

string MPATH = "/net/wonderland/home/saichen/Morphling/";

int main(int argc, char * argv[])
{
	ifstream qlist;
	qlist.open( argv[1] );
	CheckInputFileStatus(qlist, argv[0]);
	string line;
	while( getline( qlist, line ) ) {
		ResetQC( line );
	}
	qlist.close();
	return 0;
}

void ResetQC( string & dname )
{
	if ( dname[dname.size()-1] != '/' )
		dname += '/';
	string qc_dir = dname + "QC/";
	string ex_name = MPATH + "refs/ref-MEI.bed";
	
	int min_quality = -1;
	for( int mei_type = 0; mei_type<=2; mei_type++ ) {
		string outRecord = qc_dir + "refLH." + std::to_string(mei_type) + ".report";
		string mei_name = GetMeiNameFromIndex( mei_type );
		string log_name = qc_dir + "QC.log";
		string bed_name = MPATH + "refs/ref-MEI." + std::to_string( mei_type );
		if ( min_quality < 0 )
			min_quality = SetQualThreshold( outRecord, qc_dir, bed_name, ex_name );
		string power_cmd = MPATH + "bin/Evaluate-ctrl-performance.pl -i " + outRecord + " -m " + mei_name + " -o " + log_name + " -r " + bed_name + " -e " + ex_name + " -t " + qc_dir + " -q " + std::to_string(min_quality);
		if ( mei_type > 0 )
			power_cmd += " --append";
		ExecuteCmd( power_cmd );
	}
}


int SetQualThreshold( string & rec, string & qc_dir, string & bed_name, string & ex_name )
{
	int minq = -1;
	pid_t pid = getpid();
	int np = (int)pid;
	string tmp_log = qc_dir + "QC." + std::to_string(np) + ".tmp";
	for( int q=10; q<=30; q+=5 ) {
		string power_cmd = MPATH + "bin/Evaluate-ctrl-performance.pl -i " + rec + " -m " + std::to_string(q) + " -o " + tmp_log + " -r " + bed_name + " -e " + ex_name + " -t " + qc_dir + " -q " + std::to_string(q);
		if ( q != 10 )
			power_cmd += " --append";
		ExecuteCmd( power_cmd );
	}
	minq = getQualThredFromFile( tmp_log );
	if ( minq > 0 ) { // refine between (minq-5, minq)
		for( int q=minq-4; q<minq; q++ ) {
			string power_cmd = MPATH + "bin/Evaluate-ctrl-performance.pl -i " + rec + " -m " + std::to_string(q) + " -o " + tmp_log + " -r " + bed_name + " -e " + ex_name + " -t " + qc_dir + " -q " + std::to_string(q);
			if ( q != minq-4 )
				power_cmd += " --append";
			ExecuteCmd( power_cmd );
		}
		int min2 = getQualThredFromFile( tmp_log );
		if ( min2 > 0 ) // if return -1 (all>150), then keep using minq in round 1
			minq = min2;
	}
	else {
		cerr << "Warning: All Qs have N>150, use Q=30!" << endl;
		minq = 30;
	}
// set qc.info
	string cmd = string("echo ") + "\"Qual" + "\t" + std::to_string(minq) + "\" >> " + qc_dir + "QC.info";
	ExecuteCmd( cmd );
// clear
	string rm_cmd = "rm -f " + tmp_log;
	ExecuteCmd( rm_cmd );
	return minq;
}

// return first q encounter when n < 150
// if all > 150, minq = -1
int getQualThredFromFile( string & qlog_name )
{
	int minq = -1;
	std::ifstream qlog;
	qlog.open( qlog_name.c_str() );
	CheckInputFileStatus( qlog, qlog_name.c_str() );
	bool header = 0;
	string line;
	while( getline( qlog, line ) ) {
		if ( !header ) {
			header = 1;
			continue;
		}
		std::stringstream ss;
		ss << line;
		string field;
		string qstr;
		getline(ss, qstr, '\t');
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		if ( !std::all_of( field.begin(), field.end(), isdigit ) ) {
			cerr << "ERROR: 4th field in " << qlog_name << " is not numeric!" << endl;
			exit(1);
		}
		int ncount = stoi( field );
		if ( ncount < 150 ) {
			minq = stoi( qstr );
			break;
		}
	}
	qlog.close();
	return minq;
}

string GetMeiNameFromIndex( int mt )
{
	string str;
	if ( mt == 0 )
		str = string("ALU");
	else if ( mt == 1)
		str = string("L1");
	else if ( mt == 2)
		str = string("SVA");
	else {
		cerr << "ERROR: [ComputeLHMEI: GetMeiNameFromIndex] mt = " << mt << ", out of range!" << endl;
		exit(1);
	}
	return str;
}
