#ifndef REFSTATS_H
#define REFSTATS_H

#include <vector>
#include <iterator>
#include <map>

#include "OriginalStats.h"

using std::vector;
using std::string;
using std::map;

typedef struct {
	int dups;
	vector<float> log_frac;
} StatCell;

typedef vector<StatCell>::iterator StatCellPtr;

/*** het index related structs ***/

typedef struct {
	int location;
	int lift_over;
} Loc;

typedef struct {
	Loc hom;
	vector< Loc > hets;
} HetRecord;

typedef vector< vector<HetRecord> > AllHetIndex; // read het index only once

/*** end *********/

/** For remove self loci: AdjustUsedLoci ****/

typedef struct {
	unsigned int hom_index; // index in hom stats
	vector< unsigned int > hets; // index in refStats[1];
} SelfHomElement;

typedef struct {
	unsigned int lift;
	unsigned int stat_index; // index in neg stats: 3.18 added: no het will exist in neg
} SelfNegElement;

typedef map< int, SelfNegElement > SelfSingleMap;

/***********/

class RefStats
{
  public:
  	RefStats( string & ctrl_proper_prefix, string & disc_prefix, int mei_type, AllHetIndex & allHetIndex ); // constructor from ctrl stat file
  	~RefStats();
  	
  	void SetRecordGL( MergeCellPtr & merge ); // per record
  // ref GL
	void SetCtrlGLs();
  	void AdjustUsedLoci( OriginalStats* dataOsPtr ); // remove self-contribution from LH
  	void MarkRefLHasDone(); // mark ref as done
  	void ReAdjustSelfsWithLiftOver();
  	void PrintCtrlGLasRecord( string & outRecord, string & ctrl_bam, string & ctrl_fasta ); // generate ctrl lh, report %power, #novel
  	
	void PrintRefStats( string & out_prefix ); // debug function: print to out_prefix.0 neg .1 het, .2 hom

  private:
// inner data
  	OriginalStats * osPtr; // ref counts
  	vector< vector< StatCell> > refStats; // 0 to 2: neg het hom (according to gt)
  
  	void setSingleRefGL( vector< MergeCell >::iterator & merge, int & s_ref );
// build ref stats
	void setStatRef( AllHetIndex & allHetIndex ); // main set function
	void generateIndexBasedGenomeLink();
	void setHomAndHetRef( AllHetIndex & allHetIndex );
	void setNegRef( AllHetIndex & allHetIndex );
  	
  // re-forge genome link
  	void _setOtherIndex( vector<int> & OtherIndex );
  	float _getAvrLogFrac( float hom_log, float neg_log );
	void _destroyGenomeLink( StatCellPtr & ptr );
	vector< StatCellPtr > indexBasedGenomeLink; // for adding
	vector< StatCell > copiedStats; // copied from merge data
	
  // for adjust
	void adjustSingleGL( vector<int> & counts, StatCellPtr & stat_ptr, MergeCellPtr & special_ptr, int s_ref );
	void setSingleRefGLWithExclusion(  vector<int> & counts, StatCellPtr & stat_ptr, MergeCellPtr & special_ptr, int s_ref );
	void adjustMultipleGL( vector<int> & counts, vector< unsigned int> & hets, MergeCellPtr & special_ptr );
	void updateSelfSingleMapKeys( SelfSingleMap & selfMap );
  	map< int, SelfHomElement > selfHom; // 1st is genome index
  	SelfSingleMap selfHet; // remove the het contribution
  	SelfSingleMap selfNeg; // remove the neg contribution
	
  // indicator
	bool ref_lh_set; // is ref lh set?
	
  // parameters
  	const unsigned int CountSize;
	const int min_log; // = log( 10e-6 ), minimal frac if counts equal zero
	const int mei_index;
	const float CountOffSet; // offset for 0 in log_frac
	string REF_CHR;
	StatCellPtr NullStatCellPtr; // indicate null stat cell ptr
	vector<StatCell> Dummy;
	vector< int > RefCounts;
};

void SetAllHetIndex( string & het_index_name, AllHetIndex & allHetIndex );

// utility function
int GetHomCoordIndex( int location );

int GetHetCoordIndex( int location );



#endif