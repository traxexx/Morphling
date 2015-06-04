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

void ReGenotype( Options * ptrMainOptions );

// read sample list & parse
void LoadSampleList( string & sample_list_name, vector<vector<string> > & SampleList );

// some global options
void SetGenotypeGlobalOptions( Options * ptrMainOptions );

void SetGenotypeGlobalParameters( Options * ptrMainOptions );

void SetGenotypeReadMapGlobals( string & qinfo_name );

// print command to file
void PrintParallelCommand( ofstream & paraFile, vector<string> & subsample, string & site_list_name, Options* ptrMainOptions, string & mt, string & current_chr, string & rgdir );

// re-genotype single sample
void ReGenotypeSingleVcf( vector<RefSeq*> REF_SEQ, vector<int> & siteVec, vector<string> subinfo, string & rg_dir, string & mt, string & current_chr );

// used in ReGenotypeSingleVcf
void implementSingleVcf( vector<int> & siteVec, vector<RefSeq*> & REF_SEQ, RefStats & rstats, SamFile & samIn, SamFileHeader & samHeader, string & in_vcf_name, string & out_vcf_name );

#endif