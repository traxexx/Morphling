#ifndef DISCOVERY_H
#defing DISCOVERY_H

#include <string>
#include <vector>
#include "OriginalStats.h"

using std::string;
using std::vector;

typedef struct {
	string genotype;
	vector<int> pl;
} variant; 


// main function for doing discovery
void DiscoverMeiHits( const char* proper_prefix, const char* disc_prefix, const char* focus_chr,
	const char* ctrl_proper_prefix, const char* ctrl_disc_prefix, const char* out_prefix );


// print set GL to Vcf	
void PrintMeiHitsAsVcf( vector< MergeCell > & MergeData, const char* vcf_name )

/********** Utility functions *****/


/**** used in print vcf ******/

int GetVecSum( vector<int> & counts ); // get sum of int vec

int GetMeiCountSum( vector<int> & counts ); // get mei read sum of int vec

void SetGTandPL( Varint * variant, vector<int> & GL ); // set varint info from GL


#endif