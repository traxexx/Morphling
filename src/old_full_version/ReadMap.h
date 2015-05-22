#include <string>
#include <vector>
#include <deque>
#include <iterator>
#include <utility>
#include "SamFile.h"
#include "FastaUtil.h"

#include "MapData.h"
#include "ProcessReadPairs.h"
#include "PL.h"
#include "ControlStats.h"

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


/* deque for adding proper from coord sort bam */
typedef struct {
	int MatePosition;
	std::string ReadID;
	char ClipType; // 0x1 clip exist; 0x2 Sufficient length; 0x3 L(1) R(0); 0x4; Map to Alu; 0x5; LINE1; 0x6 SVA;
	int ClipLength;
} ProperDeckCell;


typedef std::map<std::string, std::vector< GenomeLocationCell > > genomeLocationLink;
typedef std::map<std::string, std::vector< GenomeLocationCell > >::iterator genomeLocationLinkIterator;


class ReadMap
{
  public:
  
  // constructor & destructor
    ReadMap(std::string & sample_name, int & winlen, int & steplen, std::string & ref_chr, std::string & ref_fasta, std::string & mei_list);
    ~ReadMap();
    
	void SetMapFromBam( std::string & bam, std::string & ref_fasta, std::string & work_dir,
		std::string & mei_list, std::string & mei_coord_list );
	void SetControls( std::string & ctrl_bam, std::string & work_dir );
	void SetDataPL ( int & mei_type, std::string & work_dir, std::string & ctrl_dir, std::string & het_index );
	void PrintToVcf( std::string & outVcfName, int & mei_type );
	
	// Read type count from raw bam
	std::map<std::string, int> ChrLenTable; // adjusted by win_len & step_len
	// disc map
	std::map<std::string, std::vector<DiscMapCell> > DiscRawMap;
	std::map<std::string, std::vector< std::vector<Coord> > > MeiCoordRef; // coord
	std::vector<RefSeq*> REF_SEQ;
	
// SHOULD I HAVE THIS ONE BUILD-IN (we have ALU & L1 & SVA)?	 final merge map
	std::map<std::string, std::vector<ReadMapCell> > RawReadMap; // only for read in SetDataPL
	std::vector<MergeCell> MergeMap;
	genomeLocationLink GenomeLocationLink; // only contain wds with mei-reads
	// something extra for windows on ctrl chr (with some excluded ctrl-windows)
	std::vector< SpecialPLcell > SpecialPLs; // for some windows in REF_CHR
	std::vector< MergeCell > specialMergeMap; // for wd with exclusion
	
	// control part
	std::vector<StatCell> MergeCtrlMap;
	std::vector< std::vector<StatCell>::iterator > CtrlLocationLink; // this one contains all windows in chr20
	std::vector< std::vector<StatCell>::iterator > HomCtrl;
	std::vector< std::vector<StatCell>::iterator > NegCtrl;
	std::vector<HetCtrlCell> HetCtrl;
	std::vector<ReadMapCell> RawCtrlMap; // destroyed after merge is set
	
  private:
 	void initializeMeiSeqRef(std::string & mei_list);
 	void initializeChrLenTable(std::string & ref_fasta);
 	void initializeMeiCoordRef( std::string & mei_coord_list);
 	 
  
  // raw-bam process
  	void processProperReads( std::string & bam, std::string & work_dir, std::string & disc_name, std::string & focus_chr );
	void processDiscReads( std::string & disc_name, std::string & work_dir, std::string & focus_chr );
  	
  	void addToDeck (SamRecord & sam_rec);
  	void addPairToProperRawMap (std::vector<ProperMapCell> & properMap, SamRecord & sam_rec);
  	void printProperRawMap(std::vector<ProperMapCell> & properMap, std::string & work_dir, std::string & current_chr);

	void addDiscPairToRawMap (DiscReadPair * rp, SamRecord & sam_rec);
	void addDiscSingleToRawMap(DiscReadPair * rp, SamRecord & sam_rec, bool mate_mapped, bool first_anchor);
	void printDiscRawMap(std::string & work_dir, std::string current_chr);

  // ctrl related
  	void setGenomicIndexToControl();
  	void setCtrlLiftOver(std::string & het_index_name);
  	void setCtrlStat( int & mei_type, std::string & het_index );
  	
  // data related
	void setDataFromRawCounts();
  	void setPerMergeCellPL (std::vector<MergeCell>::iterator & current_merge_it);
  	
  	std::deque < std::vector<ProperDeckCell> > properDeck;
  	int deck_start;
  	std::string LastChr;
  	int ReadMapCellSize_;
  	
  	int win_len;
  	int step_len;
  	std::string REF_CHR;
  	
// for QC
  	int Supplementary_Info;
	int PCR_Duplicates;
	int QC_Fail;
	int Secondary_Alignment; 
	
// constants
	int MinQuality;
	int MinClip;
	int ShortClip;
	std::string SampleName;
};