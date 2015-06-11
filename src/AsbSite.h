#ifndef ASBSITE_H
#define ASBSITE_H

#include <string>
#include <vector>
#include <map>
#include "SamFile.h"
#include "SamFileHeader.h"
#include "SamRecord.h"

using std::vector;
using std::map;
using std::string;

// for storing mapping info of each read on each mei subtype
struct EviInfo
{
	int Boundary; // if minus, do not consider any mapping results on left; otherwise right.
	int LAlign;
	int RAlign;
	int Score; // SW score
	string Seq; // sequence. if map to '-' strand, revert it
};

// temporarily store discordant 2nd read info, inferred from 1st read
struct DiscInfo
{
	string Chr;
	int Position;
	int MatePosition;
};


// single site assembly
class AsbSite {
 public:
 	AsbSite( string & chr, int position, bool* gtList, vector<SamFile> & BamFiles, vector<SamFileHeader> & BamFileHeaders, vector<string> & MEseqs );

// get methods
	bool IsAssembled();
	int GetPosition();
	int GetSubtype();
	int GetSVlength();
	float GetSVdepth();
	int GetMissingBaseCount();
	int GetLeftMost();
	int GetRightMost();
	bool GetStrand(); // +=TRUE, -=FALSE
	int GetSampleCount();
 
	void Assembly();
 	
 private:
 	void setClusterInfoFromBams(); // store read info
 	void findMostLikelyCluster(); // set cluster & pars
 
 // used in set cluster	
 	void add2ClusterFromSingleBam( int sp, int cluster_start, int cluster_end );
 	void addProper2Cluster( SamRecord & sam_rec );
 	void addDisc2Cluster( SamRecord & sam_rec );
 	void setEviInfoByRemap( string & seq, vector< map<int, vector<EviInfo> > > & evec, bool boundary, bool lbound );
 
 // used in finalize the cluster
	void setAssemblyInfo( map<int, vector<EviInfo> > & cluster1, map<int, vector<EviInfo> > & cluster2 ); 
 
 	string Chr;
 	int Position;
 	bool* GtList;
 	vector<SamFile> * pBamFiles; // point to opened bam files
 	vector<SamFileHeader> * pBamFileHeaders;
 	vector<string> * pMEseqs; // decide by mei type
 	
 	bool assembled;
 	vector< vector< map<int, vector<EviInfo> > > > Clusters; // left(+/-) right (+/-). strand based on MEI sequence (NOT genomic flanking)
 	int sample_count;
 	int subtype;
 	bool plusstr;
 	int left_most;
 	int right_most;
	int missing_base;
	int basecount;
};

string ReverseCompSeq( string & seq );
char GetCompNt( char nt);

#endif

