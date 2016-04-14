#ifndef REFSTATS_H
#define REFSTATS_H

#include "singleSample.h"
#include "lhData.h"
#include "SamFile.h"
#include <map>
#include <vector>

using std::map;
using std::string;
using std::vector;

void GenerateRefStats(
	SingleSampleInfo & msinfo, vector< vector<string> > & MEseq, vector< vector<string> > & MEname, 
	string & vcf_name, string & out_prefix, bool ref_exclusion );

// print het only
void PrintRefStats( vector< vector<float> > & het_stat, string & out_name );

void loadUnliftBed( std::map<int, int> & ref_coord, string & ref_bed_name);

// split neg ctrl base bed by window size
void splitBaseBed( string & base_ctrl_bed, string & split_out, int nex );

void makeRegionBam( map<int, int> & ref_coord, SingleSampleInfo & si, string & remap_bam_prefix);

void remapAndClean( string & remap_bam_prefix );

void ExecuteCmd( string & cmd );

void generatePassVcf( string & raw_vcf_name, string & pass_vcf_name );

#endif