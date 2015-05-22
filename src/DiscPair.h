#ifndef DISCPAIR_H
#define DISCPAIR_H

#include "SamRecord.h"
#include "FastaUtil.h"
#include <string>

using std::string;

typedef struct {
	string chr;
	int st;
	int ed;
	bool mateMap; // mate is unmapped(0) or disc (1)
	bool sense_strand;
	int type;
} Loci;

class DiscPair
{
 public:
	DiscPair( bool first_in_candidate, bool second_in_candidate, SamRecord & sam_rec, vector<RefSeq*> & refSeq );
	~DiscPair();

	bool IsSamePair( SamRecord & sam_rec );
	void AddSecondToPair( SamRecord & sam_rec );
	
	bool FirstValid();
	bool SecondValid();
	
	int GetFirstAlignPosition();
	int GetSecondAlignPosition();

	Loci GetFirstLoci();
	Loci GetSecondLoci();
	
	
 private:
 
	bool valid;
	bool in_candidate1;
	bool in_candidate2;
	string chr_name1;
	string chr_name2;
	bool strand1; //+ 1 - 0
	bool strand2;
	int loc1;
	int loc2; // the above 4 for checking if the same read
	string seq1;
	string seq2;
	bool high_qual1;
	bool high_qual2;
	bool is_map1;
	bool is_map2;
	int align1_end;
	int align2_end;
	int em1;
	int em2;
	
	vector<RefSeq*> rSeq;
};

#endif




