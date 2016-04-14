#ifndef GENREF_H
#define GENREF_H

#include "lhData.h"
#include "SamFile.h"
#include "singleSample.h"
#include <map>
#include <vector>
#include <string>

using std::string;
using std::vector;

class genRef
{
  public:
  	genRef( string & ref_bed_name, SingleSampleInfo & si,
      vector<string> & subMEseq, vector<string> & subMEname, bool neg_ctrl );
  	void Print( string & out_name );
  	void Export( vector< vector<float> > & het_stat ); // export to het

  private:
  	void transferRawStatToRefStat();
  	void transferSingleRawStatToRefStat(int index, vector<float> & raw);
    static bool sortRawStat( const vector<float> & x, const vector<float> & y );

    string bam_name;
    int avr_read_length;
    float max_ins_size;
    float min_frac; // minimum fraction of reads in window when count is zero
    vector<string> * psMEseq;
  	vector< vector<float> > rawStat; // fraction. in order.
  	RefStat refStat; // log10 scaled fraction. deduped
};

void loadRefBed( std::map<int, string> & ref_coord, string & ref_bed_name, bool neg_ctrl);

bool equalFloatVector( vector<float> & f1, vector<float> & f2 );

bool SetWindowStatFromBam( vector<int> & stat, int avr_read_length, float max_ins_size,
  string & mei_seq, string & chr, int st, int ed, SamFile & bam, SamFileHeader & bam_header, 
  SamFile & alt_sam, SamFileHeader & alt_sam_header, bool both_strand );

int getSingleTypeIndex( SamRecord & rec, SamRecord & mate_rec, string & mei_seq );

#endif