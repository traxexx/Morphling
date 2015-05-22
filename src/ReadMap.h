#include <string>
#include <vector>
#include <map>
#include <utility> // pair
#include "SamFile.h"
#include "FastaUtil.h"

#include "QC.h"
#include "ProperDeck.h"
#include "DiscPair.h"

using std::vector;
using std::map;
using std::string;

/*
1. Read Proper. For each chr, read & print
  print format ( as delimiter):
    MEI? stats
  file: proper + short-clip
		alu
		line
		sva
2. Same. Do disc.
3. For each chr, Load Proper & Read Disc.
	for eacho MEI:
		if (disc-MEI or proper-MEI) print stats
4. read stats & sort & do likelihood (do special for REF_CHR)
*/


typedef vector< vector<int> > ProperMapCell; // proper (2), A (8), L (8), S (8): direction * is_begin * remap

typedef struct {
	bool valid; // see if this region should be excluded
	vector< vector<int> >  stats;
} DiscMapCell;

class ReadMap
{
  public:
  
  // constructor & destructor
    ReadMap( string & mei_list);
    ~ReadMap();
    
	void SetMapFromBam( string & bam, string & disc_bam, string & output_prefix, string & focus_chr, string & level_file_prefix);
	void SetMapFromCtrlBam( string & ctrl_bam, string & output_prefix, string & mei_coord_list );
	
  private:
  // use in set map from bam
  	void setCandidateRegion(  vector< std::pair<int, int> > & candidate_region, int chr_length, string & level_name);
  
// process proper bam
  	void processProperReadsBySection( SamFile & samIn, SamFileHeader & samHeader, string & chr_name, int st, int ed, QC * ProperQC);	
  	void addRetrievedIndexToProperMap( RetrievedIndex & rIndex );
  	void printProperMapBySection(string & output_prefix);
  // related utilities
  	void _setProperMapAddIndex (vector<int> & add_index, RetrievedIndex & rIndex);
  
// process disc bam	
	void processDiscReads( string & discBam );
	void initializeDiscMap( SamFileHeader & samHeader );
	void printDiscMap( string & output_prefix );
  // related util
	bool isSingleReadInCandidate( string & chr_name, const int loc);
	void addLociToDiscMap( Loci & loci );
	int getIndexInDiscMapCellFromLoci( Loci & loci , int which_mei);


/*** data structors ***/
  // proper map
  	vector<ProperMapCell> properMap;
  // disc map
  	map<string, vector< std::pair<int, int> > > discCandidates; // clear after setting discMap;
  	map<string, vector<DiscMapCell> > discMap;
  // ref mei sequence
	vector<RefSeq*> REF_SEQ;
};

void InitializeMeiSeqRef( vector<RefSeq*> & REF_SEQ, string & mei_list );



