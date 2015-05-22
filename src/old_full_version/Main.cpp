#include <iostream>
#include <string>
#include "ComputeLHMEI.h"

void LHMEI_Help();

int main(int argc, char * argv[])
{
	if (argc <= 1) {
		LHMEI_Help(); return 0;
	}
	
	std::string ArgString;
	std::string Dummies;
	std::string PATH = std::string("/net/wonderland/home/saichen/LHMEI_v4");
	
	ArgString += "-Win=600;-Step=100;-MElist=" + PATH + "/refs/MobileElement.list;-MEcoord=";
	ArgString += PATH + "/refs/MobileElement.coord;-Mapper=bwa;-HetIndex=";
	ArgString += PATH + "/refs/hs37d5-chr20-MEI-slice.het-index;";
	
	Dummies = std::string("--verbose");
	
	std::string FirstArg = std::string(argv[1]);
	if (FirstArg.compare("Test") != 0) {
		ArgString += std::string("-Sample= ;-Bam= ;-WorkDir= ;-GenomeFasta= ;-BinDir= ;");
		Dummies += ";--keepIntermediates";
	}
	else { // test mode (intermediate file is default kept)
		argc--;
		argv++;
		std::cout << "  Running test mode..." << std::endl;
		ArgString += std::string("-Sample=TestSample;-Bam=") + PATH + "/usage_test/1612_test.bam;";
		ArgString += "-GenomeFasta=/net/wonderland/home/saichen/reference/archive/hs37d5.fa;-WorkDir=";
		ArgString += PATH + "/usage_test/output;-BinDir=" + PATH;
	}
	Options MainOptions( argc, argv, ArgString, Dummies );
	
	Options * ptrMainOptions = &MainOptions;
	ComputeLHMEI(ptrMainOptions);

	return 0;	
}


using std::cout;
using std::endl;

void LHMEI_Help ()
{
	cout << endl;
	cout << "LHMEI Options: [] is default" << endl;
	cout << "    Test:  run test commands. Override all options below except for --verbose." << endl;
	cout << endl;
	cout << "  Required:" << endl;
	cout << "    -Sample:  Sample name." << endl;
	cout << "    -Bam:  Raw bam files for LHMEI discovery. MUST be sorted by cooridinate." << endl;
	cout << "    -WorkDir:  Working directory (if not exist, LHMEI will create one)." << endl;
	cout << "    -BinDir: Directory of LHMEI." << endl;
	cout << "    -MElist:  MEI Consensus sequence list. [refs/MobileElement.list]" << endl;
	cout << "    -GenomeFasta:  Reference genome sequence. [hs37d5.fa]" << endl;
	cout << "  With defaults:" << endl;
	cout << "    -Mapper:  Mapper fore genearating -Bam. [bwa]" << endl;
	cout << "    -HetIndex:  Contrl stats of sliced MEI. [/LHMEI-dir/refs/hs37d5-chr20-MEI-slice.het-index]" << endl;
	cout << endl;
	cout << "  Optional:" << endl;
	cout << "    -MEcoord:  Genomic regions of MEs. To speed up. Recommended. [ /refs/MobileElement.coord ]" << endl;
	cout << endl;
}
