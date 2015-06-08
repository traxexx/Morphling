#include "Assembly.h"
#include "Sites.h"
#include "iostream"

using std::cout;
using std::endl;

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
// load sites & bam list
	Sites allSites( ptrMainOptions->ArgMap["Vcf"], ptrMainOptions->ArgMap["SampleList"] );

// cluster reads in each sites & print
	allSites.AssemblySubtype( ptrMainOptions );

// close all bam
	allSites.Close();

// report finish
	cout << "Morphling Assembly finished with no error reported. Check final vcf at: " << ptrMainOptions->ArgMap["Out"] << endl;
}



























