#ifndef PROPERDECK_H
#define PROPERDECK_H

#include <string>
#include <deque>
#include <vector>
#include "SamFile.h"
#include "FastaUtil.h"

using std::string;
using std::deque;
using std::vector;


/* deque for adding proper from coord sort bam */
typedef struct {
	bool valid; // if valid == 0, do not use
	int MatePosition;
	string ReadID;
	char ClipType; // 0x1 clip exist; 0x2 Sufficient length; 0x3 L(1) R(0); 0x4; Map to Alu; 0x5; LINE1; 0x6 SVA;
	int ClipLen;
	int FirstAlignEnd;
} ProperDeckCell;


/* structure for returning info from proper deck and adding to raw proper map */
typedef struct {
	bool valid; // if not valid, not retrieved
	int min_position;
	int max_position;
	char type; // 1-6: exist enough-length at-begin A L S
	bool mei_on_1st; // if mei clip is on 1st read
} RetrievedIndex;

// handling proper reads
class ProperDeck
{
  public:
  	ProperDeck( vector<RefSeq*> & Ref_Seq );
  	~ProperDeck();
  	void Add( SamRecord & sam_rec );
  	RetrievedIndex RetrieveIndex( SamRecord & sam_rec );
  	
  private:
  	char getSingleClipType ( SamRecord & sam_rec);
  
  	deque< vector<ProperDeckCell> > properDeck;
  	int deck_start;
  	vector<RefSeq*> rSeq;
};

#endif

