#ifndef MULTISAMPLECALLING_H
#define MULTISAMPLECALLING_H

#include "Options.h"
#include "FastaUtil.h"
#include "SamFile.h"
#include "SamFileHeader.h"
#include <string>
#include <vector>
#include <map>
#include <utility>
#include "RefStats.h"

using std::string;
using std::pair;
using std::vector;
using std::map;

void MultiSampleCalling( Options * ptrMainOptions );

void ReGenotypeSingleVcf( vector<int> & siteVec, vector<RefSeq*> & REF_SEQ, RefStats & rstats, SamFile & samIn, SamFileHeader & samHeader, string & in_vcf_name, string & out_vcf_name );

// read sample list & parse
void LoadSampleList( string & sample_list_name, vector<vector<string> > & SampleList );


// some global options
void SetGenotypeGlobalOptions( Options * ptrMainOptions );

void SetGenotypeGlobalParameters( Options * ptrMainOptions );

void SetGenotypeReadMapGlobals( string & qinfo_name );

#endif