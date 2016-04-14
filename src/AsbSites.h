#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "SamFile.h"
#include "Aligner.h"
#include "singleSample.h"

/* For each site, read & store mapped read info */

using std::vector;
using std::string;
using std::map;

/* mapping info to each subtype.
   left & right cluster are separated. */
struct SubtypeInfo {
	vector<float> coverage; // score / MATCH / mapping-length
	int evidence; // #supporting reads. NOT adjusted with 2 anchor refinement!
};

/* contain subtype info to one site */
struct AsbSiteListElement {
	int nclip;
	int ndisc;
	int raw_breakp; // estimated breakpoint -> refined breakpoint
	int true_breakp;
	int disc_breakp;
	int extend; // extend of cluster from breakp
//	int n_spanning; // #spanning reads across the breakpoint
//	int depth_sum; // window flanking base count
//	bool doing_rescue;
	int polyA_left;
	int polyA_right;
	int polyT_left;
	int polyT_right;
	vector< vector<SubtypeInfo> > p_cluster; // 4 cluster: plus-left-right | minus ...
} ;

/* all sites */
typedef map<string, map<int, AsbSiteListElement> > AsbSiteList; // chr -> raw_breakp -> each subtype info

/* site-level adjusted assembled info */
struct CompleteSiteInfo {
	bool assembled; // indicate if this site is assembled
	int breakpoint;
	bool imprecise;
	int mei_type;
	float score_sum;
	bool strand;
	int subtype;
	int left_most;
	int right_most;
	int left_anchor_length;
	int right_anchor_length;
	int n_left_poly; // #polyA/T reads in left cluster
	int n_right_poly; // #polyA/T reads in right cluster
	int n_left_anchor;
	int n_right_anchor;
	int svlen;
	int n_missing;
	float avr_mei_depth;
	int n_spanning;
	float avr_flank_depth;
	bool stringent_filter; // FALSE means pass
	int filter;
	// indicate if other faimilies are merged. Used in rescue stage in Assembly()
	bool merged[3];
//	bool merge0; // if ALU is merged?
//	bool merge1; // if L1 is merged?
//	bool merge2; // if SVA is merged?
};

/* all adjusted assembled site info */
typedef map<string, map<int, CompleteSiteInfo> > CompleteSiteList; // chr -> round(position/WIN) -> csi info

class AsbSites
{
public:
	AsbSites( vector< vector<string> > & mei_names, vector< vector<string> > & mei_seqs );
	void Assembly( map<int, map<string, map<int, int> > > & candidate_sites, vector<SingleSampleInfo> & msinfo );
	void PrintToVcf(string & vcf_name);

private:
	AsbSiteList * pCandidateSites; // data structure to record info from samples
	CompleteSiteList compSites; // adjusted from pCandidateSites

	void LoadCandidateSitesFromSites( map<int, map<string, map<int, int> > > & candidate_sites );
	void AddEvidenceFromSingleBam( SingleSampleInfo & si );
	void ImplementBreakpFromSingleBam( SingleSampleInfo & si );

	void assemblySites(); // do assembly and remove overlapped sites
	void rescueCsi( vector<SingleSampleInfo> & msinfo ); // secondary assembly
	void removeUnAssembledSites(); // do this after rescue
	void removeNearbySites(); // final step in assembly

	void constructRescueAsbList();

	void addSectionEvidence( int m, map<int, AsbSiteListElement>::iterator & it, SamFile & bam, SamFileHeader & bam_header, SamFile & alt_bam, SamFileHeader & alt_bam_header);
	void addProperRead( int m, map<int, AsbSiteListElement>::iterator & it, SamRecord & sam_rec );
	void addDiscRead( int m, map<int, AsbSiteListElement>::iterator & it, SamRecord & sam_rec, SamFile & alt_bam, SamFileHeader & alt_bam_header );
	bool addToSubtypeInfo( int m, string & seq, map<int, AsbSiteListElement>::iterator & it, bool is_left_cluster );
	void refineAnchorCluster( SubtypeInfo & sbi );
	bool _realignAndAddSingleSide( int m, vector<SubtypeInfo> & p, map<int, AsbSiteListElement>::iterator & it, string & seq );
	void _updateBaseCoverageWithRead( vector<float> & coverage, Aligner & al );
	void updateTrueBreakp( AsbSiteListElement & ase, int breakp );
	void updateDiscBreakp( AsbSiteListElement & ase, SamRecord & rec );
	void setCsiBreakPoint( CompleteSiteInfo & csi, AsbSiteListElement & ase );

	bool setCompleteSiteInfoFromAsbSiteListElement( int m, CompleteSiteInfo & csi, AsbSiteListElement & ase );
	bool isNewCsiBetter( CompleteSiteInfo & new_csi, map<int, CompleteSiteInfo>::iterator & t );

	void implementSection( SamFile & bam, SamFileHeader & bam_header, CompleteSiteInfo & csi );
	void addToSpanning( SamRecord & rec, CompleteSiteInfo & csi );
	void addToDepth( SamRecord & rec, CompleteSiteInfo & csi );
//	bool isGoodSectionDepth( const char* chr, int st, int ed, SamFile & bam, SamFileHeader & bam_header ); // move to discover phase

	void printHardFilterSiteVcfHeader( std::ofstream & vcf );
	void printSingleCsiToVcf( string & chr, std::ofstream & vcf, CompleteSiteInfo & csi );
	void _setRecordStatics();
	bool _convertMostLikelySubtypeToCsiInfo( CompleteSiteInfo & csi, int mei_type, bool is_plus_strand, int optimized_subtype, AsbSiteListElement & ase );
	int _getLeftMost( vector<float> & coverage );
	int _getRightMost( vector<float> & coverage );
	int _getAnchorLength( vector<float> & coverage );
	int _getMissingBaseCount( vector<float> & coverage );
	int _getSumDepth( vector<float> & coverage );
	void _setCsiFilter( CompleteSiteInfo & csi );
	void _setCsiStringentFilter( CompleteSiteInfo & csi );

	float _sumFloatVector( vector<float> & v ); // because std::accumulate dosen't work!
	int _getCompSiteCount();

	// MEseqs
	vector< vector<string> > * pMEseqs;
	vector< vector<string> > * pMEnames;

	// multi-sample stats
	bool first_add;
	string stat_chr; // chr to get ins size and read length. 1st of smallest m sites
	int sample_count;
	vector<float> sample_depth;
	float total_depth;
	int avr_read_length_across_all;
	int avr_ins_size_across_all;

	// vcf print related
	int MIN_SPAN;
	int MIN_QUALITY; // min anchor quality
	int avr_ins_size; // dynamic
	int avr_read_length;
	int min_read_length;
	string current_bam_name;
	map<int, string> filter_name;
	map<bool, char> strand_name;
};


