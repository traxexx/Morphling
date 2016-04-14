#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "SamFile.h"
#include "lhData.h"
#include "singleSample.h"
#include "Globals.h"
#include <string>
#include <vector>

using std::vector;
using std::string;

// class to store ref stats & calculate LH
class Likelihood
{
  public:
  	Likelihood(  vector< vector<string> > & MEseq, SingleSampleInfo & si, string & ref_file_prefix );
  	void SetLhRec( lhRec & lr, string chr, int pos, int mtype, int subtype );

  private:
  	float getLog10Likelihood( int m, int g, vector<int> & wc );
    void setFromRefFile( int m, int g, string & ref_file_name );
    bool setWindowStat( vector<int> & stat, string & chr, int pos, string & mseq );
    int getFlankingCount( string & chr, int pos );

    // outer vector is mtype
  	vector< vector<RefStat> > refStat;
  	vector< vector<int> > nref; // #uniq refs
  	vector< vector<float> > log10_total; // log10(ncount)

    // for re-mapping use only
    SamFile bam;
    SamFileHeader bam_header;
    SamFile alt_bam;
    SamFileHeader alt_bam_header;
    vector< vector<string> > * pMEseq;
    int avr_read_length;
    int max_ins_size;
};

// sum up 2 log 10
float SumTwoLog10s( float l1, float l2 );

// print gt info
void PrintLhRecToVcf( std::ofstream & vcf, int gt, lhRec & lr );

#endif

