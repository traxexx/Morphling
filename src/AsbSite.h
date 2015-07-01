#ifndef ASBSITE_H
#define ASBSITE_H

#include <string>
#include <vector>
#include <map>
#include "SamFile.h"
#include "SamFileHeader.h"
#include "SamRecord.h"
#include "MeiSeqs.h"

using std::vector;
using std::map;
using std::string;

class AsbSite{
  public:
  	AsbSite(vector< vector<string> > & Seqs, subCluster & sc1, subCluster & sc2, bool strand, string & ref_seq);

// get methods
	bool IsAssembled();
	int GetSVlength();
	float GetSVdepth();
	int GetMissingBaseCount();
	int GetCorrection();
	int GetConsecutiveMiss();
	int GetLeftLength();
	int GetRightLength();
	int GetLeftMost();
	int GetRightMost();
	int GetSampleCount();
	int GetValidSampleCount();

  	
  private:
  	void findMostLikelyCenter( subCluster & sc, vector< vector<string> > & Seqs );
  	void setAssemblyInfo( subCluster & cluster1, subCluster & cluster2 );
  	void setMiddleBase( vector<int> & basecov );
  	int calculateMissingBase( vector<int> & basecov );
	int calculateBaseCount( vector<int> & basecov );
  	void setSampleCount( subCluster & cluster1, subCluster & cluster2 );
  
  	string rSeq;
  	bool assembled;
 	int sample_count;
 	int valid_sample_count;
 	int left_most;
 	int right_most;
	int missing_base;
	int basecount;
	bool corrected;
	int consecutive_miss;
	int left_length;
	int right_length;
};

// set args preAssemble -Win 600 -Bam /net/wonderland/home/saichen/Morphling/usage_test/1612_test.bam -Vcf usage_test/gt_out/final/test.vcf -Out usage_test/assembly/pres/1612.pre --verbose
// set args Assembly -SampleList usage_test/assembly.list -Vcf usage_test/gt_out/final/test.vcf -Out usage_test/assembly/final.test.vcf

#endif

