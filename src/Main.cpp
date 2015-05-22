#include <iostream>
#include <string>
#include "ComputeLHMEI.h"
#include "MultiSampleCalling.h"
#include "Wrappers.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char * argv[])
{
// if need to print help info only
	if (argc <= 1) {
		DisplayUsageInfo();
		DisplayDiscoveryUsageInfo();
		return 0;
	}
	string FirstArg = string(argv[1]);
	 if (FirstArg.compare("-hello") == 0) { // display debug usage info
		DisplayDiscoveryDetailedUsageInfo();
		DisplayDiscoveryDebugUsageInfo();
		return 0;
	}
	if (FirstArg.compare("-h") == 0) { // display detailed help info
		std::cout << std::endl;
		DisplayDiscoveryDetailedUsageInfo();
		return 0;
	}

// decide mode
	int current_mode = 0; // 0 discovery 1 genotype 2 re-genotype
	if ( FirstArg.compare( "Genotype" ) == 0 ) {
		current_mode = 1;
		if ( argc <= 2 ) {
			DisplayGenotypeUsageInfo();
			return 0;
		}
	}
	else if ( FirstArg.compare("ReGenotype") == 0 )
		current_mode = 2;
	if ( current_mode > 0 ) { // get rid of first arg
		argc--;
		argv++;
	}


// get path first
	string Path = GetExePath(); // secured last is '/'
	if ( Path.length() <= 4 ) {
		std::cerr << "ERROR: LHMEI-Discovery is not in $ProgramDir/bin/" << std::endl;
		exit(1);
	}
	Path = Path.substr(0, Path.size() - 4); // remove bin/
	string RefPath = Path + "refs/";
	
// do by mode
	if ( current_mode == 0 ) { // discovery
		string ArgString = "-Win=600;-Step=100;-CtrlChr=20;-Simplify=1;-NonOffset=1;-Chr=-1;";
		ArgString += "-GenomeFasta=/net/wonderland/home/saichen/reference/archive/hs37d5.fa;";
		ArgString += "-MElist=" + RefPath + "MobileElement.list;-MEcoord=" + RefPath + "MobileElement.coord;-HetIndex=" + RefPath + "hs37d5-chr20-MEI-slice.het-index;";
		ArgString += "-SliceFA=" + RefPath + "slice-chr20-hs37d5.fa;";
		ArgString += "-Mapper=/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-mem;";
		ArgString += "-refPrefix=refStats;-ReadLen=-1;-InsSize=-1;-MeiType=-1;-Depth=-1;";
		string Dummies = string("--verbose;--debug;--keepIntermediates;--includeSingleAnchor;--pseudoChr;--printNonVariant;--printRefStats;--noCtrlVcf;--disableDPfilter;--noRefAllele;--noBreakPoint");
		if (FirstArg.compare("Test") == 0) {  // test mode (intermediate file is default kept)
			cout << endl;
			cout << "  Running discovery test mode..." << endl;
			cout << endl;
			argc--;
			argv++;
			ArgString += "-Sample=TestSample;-Bam=" + Path + "usage_test/1612_test.bam;-WorkDir=" + Path + "usage_test/output";
		}
		else {
			ArgString += string("-Sample=.;-Bam= ;-WorkDir= ;");
		}
		Options MainOptions( argc, argv, ArgString, Dummies );
		Options * ptrMainOptions = &MainOptions;
		ComputeLHMEI(ptrMainOptions);
	}
	else if ( current_mode == 1 ) { // genotype
		string ArgString = string("-SampleList= ;-WIN=600;-WorkDir= ;-HetIndex=") + RefPath + "hs37d5-chr20-MEI-slice.het-index;";
		ArgString += "MeiType=-1;-Chr=-1;-MElist=" + RefPath + "MobileElement.list;";
		string Dummies = string("--Parallel;--verbose");
		Options MainOptions( argc, argv, ArgString, Dummies );
		Options * ptrMainOptions = &MainOptions;
		MultiSampleCalling( ptrMainOptions );
	}
/*	else if ( current_mode == 2 ) { // re-genotype single sample
		string ArgString = string("-Sample=.;-Bam= ;-Vcf= ;-CtrlDir= ;-SiteList= ;-WIN=600;");
		string Dummies = string("");
		ReGenotypeSingleSample( ptrMainOptions );
	}		
*/	
	return 0;	
}


