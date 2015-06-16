#ifndef CLUSTER_H
#define CLUSTER_H

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "SamRecord.h"
#include "MeiSeqs.h"

using std::string;
using std::map;
using std::vector;
using std::ofstream;

const int MinClip = 20;

// temporarily store discordant 2nd read info, inferred from 1st read
struct DiscInfo
{
	string Chr;
	int Position;
	int MatePosition;
};

class Cluster {
  public:
  	Cluster( int mei_index, vector<string> & MEseqs );
  	void AddProper( SamRecord & sam_rec );
  	void AddDisc( SamRecord & sam_rec );
  	void Print( ofstream & out_info, ofstream & seq_info );
  	
  	int CountPolyA( string & seq );
  	int CountPolyT( string & seq );
  	int GetClipLengthFromCigar( string & cigar );
  	
  private:
  	void setEviInfoByRemap( string & seq, vector< subCluster > & evec, bool boundary, bool lbound );
  
  	int mei_index;
  	int npolyA;
  	int npolyT;
  	vector<string> * pMEseqs; // decide by mei type
  	vector< vector< subCluster > > rClusters; // 0~4 -> subtype -> read map
  	vector<string> Seqs; // site -> sample -> find seq by SeqKey
};


#endif
