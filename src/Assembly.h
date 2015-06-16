#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "SamRecord.h"
#include "Options.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

// for a single sample, generate pre assembly info ( info + seq )
void PreAssemble( Options* ptrMainOptions );

// do assembly based on merged vcf
void Assembly( Options* ptrMainOptions );


// used in assembly
void LoadAssemblySampleList( string & vcf_name, string & sample_list_name, vector<string> & pres );

// set globals used in Assembly (also set pre-assembly globals)
void SetAsbGlobals( Options* ptrMainOptions );

// set pre-assembly globals
void SetPreAsbGlobals( Options* ptrMainOptions );

// check if all pre exist
// if not, exit with error
void CheckPreAssembles( vector<string> & pres);

// check read length & insert size
bool AssemblyDiscPass( SamRecord & sam_rec );

#endif


/* old scheme (too time consuming not used)
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
