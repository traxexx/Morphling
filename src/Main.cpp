#include <iostream>
#include <string>
#include "ComputeLHMEI.h"
#include "MultiSampleCalling.h"
#include "Wrappers.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

void RunDiscover( int argc, char * argv[] );
void RunGenotype( int argc, char * argv[] );

int main(int argc, char * argv[])
{
// if need to print help info only
	if (argc <= 1) {
		DisplayUsageInfo();
		return 0;
	}
	
// if only one argument, display corresponding usage info
	string FirstArg = string(argv[1]);
	argc--;
	argv++;
	if ( FirstArg.compare("Discover") == 0 )
		RunDiscover( argc, argv );
	else if ( FirstArg.compare("Genotype") == 0)
		RunGenotype( argc, argv );
	else {
		cerr << "ERROR: Please select mode: Discover or Genotype!" << endl;
		exit(1);
	}
	return 0; // successful execution
}


void RunDiscover( int argc, char * argv[] )
{
// help info
	if ( argc <= 1 ) {
		DisplayDiscoveryUsageInfo();
		return;
	}
	string FirstArg = string(argv[1]);
	 if (FirstArg.compare("-hello") == 0) { // display debug usage info
		DisplayDiscoveryDetailedUsageInfo();
		DisplayDiscoveryDebugUsageInfo();
		return;
	}
	if (FirstArg.compare("-h") == 0) { // display detailed help info
		cout << endl;
		DisplayDiscoveryDetailedUsageInfo();
		return;
	}
	
// get path first
	string Path = GetExePath(); // secured last is '/'
	if ( Path.length() <= 4 ) {
		std::cerr << "ERROR: LHMEI-Discovery is not in $ProgramDir/bin/" << std::endl;
		exit(1);
	}
	Path = Path.substr(0, Path.size() - 4); // remove bin/
	string RefPath = Path + "refs/";
	
// arguments
	string ArgString = "-Win=600;-Step=100;-CtrlChr=20;-Simplify=1;-NonOffset=1;";
	ArgString += "-GenomeFasta=/net/wonderland/home/saichen/reference/archive/hs37d5.fa;";
	ArgString += "-MElist=" + RefPath + "MobileElement.list;-MEcoord=" + RefPath + "MobileElement.coord;-HetIndex=" + RefPath + "hs37d5-chr20-MEI-slice.het-index;";
	ArgString += "-SliceFA=" + RefPath + "slice-chr20-hs37d5.fa;";
	ArgString += "-Mapper=/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-mem;";
	ArgString += "-refPrefix=refStats;-ReadLen=-1;-InsSize=-1;-MeiType=-1;-Depth=-1;";
	string Dummies = string("--verbose;--passOnly;--debug;--keepIntermediates;--includeSingleAnchor;--pseudoChr;--printNonVariant;--printRefStats;--noCtrlVcf;--disableDPfilter;--noRefAllele;--noBreakPoint");
	FirstArg = string(argv[1]);
	if ( FirstArg.compare("Test") == 0 ) { // test mode
		cout << endl;
		cout << "  Running discovery test mode..." << endl;
		argc--;
		argv++;
		ArgString += "-Sample=TestSample;-Bam=" + Path + "usage_test/1612_test.bam;-WorkDir=" + Path + "usage_test/output;-Chr=20;";	
	}
	else {
		ArgString += string("-Sample=.;-Bam= ;-WorkDir= ;-Chr=-1;");
	}

// run discovery
	Options MainOptions( argc, argv, ArgString, Dummies );
	Options * ptrMainOptions = &MainOptions;
	ComputeLHMEI(ptrMainOptions);	
}


void RunGenotype( int argc, char * argv[] )
{
// help info
	if ( argc <= 1 ) {
		DisplayGenotypeUsageInfo();
		return;
	}
	string FirstArg = string(argv[1]);
	if (FirstArg.compare("-h") == 0) { // display detailed help info
		cout << endl;
		DisplayGenotypeDetailedUsageInfo();
		return;
	}
	
// get path first
	string Path = GetExePath(); // secured last is '/'
	if ( Path.length() <= 4 ) {
		std::cerr << "ERROR: LHMEI-Discovery is not in $ProgramDir/bin/" << std::endl;
		exit(1);
	}
	Path = Path.substr(0, Path.size() - 4); // remove bin/
	string RefPath = Path + "refs/";	

// arguments
	string ArgString = string("-Win=600;-Step=100;-CtrlChr=20;-HetIndex=") + RefPath + "hs37d5-chr20-MEI-slice.het-index;-MElist=" + RefPath + "MobileElement.list;";
	string Dummies = string("--Parallel;--verbose;--passOnly;--debug;--keepIntermediates;--includeSingleAnchor;--pseudoChr;--printNonVariant;--printRefStats;--noCtrlVcf;--disableDPfilter;--noRefAllele;--noBreakPoint");

	FirstArg = string( argv[1] );
	if ( FirstArg.compare("Test") == 0 ) { // test mode
		cout << endl;
		cout << "  Running discovery test mode..." << endl;
		argc--;
		argv++;
		ArgString += "-SampleList=" + Path + "usage_test/gt.list;-WorkDir=" + Path + "usage_test/gt_out;-Chr=20;-MeiType=0;"; 
	}
	else {
		ArgString += "-SampleList= ;-WorkDir= ;-Chr=-1;MeiType=-1;";
	}
	Options MainOptions( argc, argv, ArgString, Dummies );
	Options * ptrMainOptions = &MainOptions;
	MultiSampleCalling( ptrMainOptions );	
}

