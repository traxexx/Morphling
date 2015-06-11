#include "Assembly.h"
#include "Sites.h"
#include "Globals.h"
#include <iostream>
#include <utility>
#include <fstream>
#include <algorithm>    // std::count

using std::cout;
using std::endl;
using std::cerr;

/*
Within MEI inner cluster:
	cluster reads on each side. Use anchor to determine strand.
	For each subtype:
		map disc cluster reads back:
			For long refs (e.g., 6k L1) use 1st disc to decide approximate mapping region.
				repeat doing this for all reads until finding max sw score
			decide inner anchor region.
		map clip reads back around inner anchor
			Map >=40 bp clip along with disc
		especially, if polyA/T, it is +/-. No mapping required.
		use max sw score as target.
	report:
		sub-type
		MEI length
		%missing middle part
		MEI average coverage
		Alt sequence ( use mapping position to locate each base )
*/

void Assembly( Options* ptrMainOptions )
{
// set globals
	SetAsbGlobals( ptrMainOptions );

// load sites & bam list
	Sites allSites( ptrMainOptions->ArgMap["Vcf"], ptrMainOptions->ArgMap["SampleList"], ptrMainOptions->ArgMap["Out"], ptrMainOptions->ArgMap["MElist"] );

// cluster reads in each sites & print
	allSites.AssemblySubtypes();

// report finish
	cout << "Morphling Assembly finished with no error reported. Check final vcf at: " << ptrMainOptions->ArgMap["Out"] << endl;
}


// global sub function
void SetAsbGlobals( Options* ptrMainOptions )
{
	WIN = stoi( ptrMainOptions->ArgMap["Win"] );
	if ( WIN <= 0 ) {
		cerr << "ERROR: win size = " << WIN << ", please specify a valida win size use -Win option." << endl;
		exit(1);
	}
	std::ifstream in_list( ptrMainOptions->ArgMap["SampleList"].c_str() );
	NSAMPLE = std::count(std::istreambuf_iterator<char>(in_list), std::istreambuf_iterator<char>(), '\n');
	in_list.close();
	if ( NSAMPLE <= 0 ) {
		cerr << "ERROR: empty sample list: " << ptrMainOptions->ArgMap["SampleList"] << endl;
		exit(1);
	}
}

























