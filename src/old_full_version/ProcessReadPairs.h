#ifndef PROCESSREADPAIRS_H
#define PROCESSREADPAIRS_H

#include "SamFile.h"
#include <Cigar.h>

#include "Params.h"
#include "FastaUtil.h"

#include <vector>
#include <string>
#include<iterator>

typedef struct {
	int start;
	int end;
} Coord;

typedef std::vector<Coord>::iterator CordPtr;

typedef struct {
	bool is_empty;
	std::string id;
	std::string seq1;
	bool first_high_qual;
	std::string seq2;
	bool second_high_qual;
	bool first_map;
	bool second_map;
	int first_alignment_end;
	int mc1; // mc tag
	int mc2;
} DiscReadPair;

// remove QO tags for disc sam
void simplifySamRec( SamRecord & sam_rec );

// set EM tag for disc sam. If unmap, EM=0
void setEMtag( SamRecord & sam_rec, std::vector<CordPtr> & CordPtrVec, std::vector<CordPtr> & EndPtrVec );



int getMaxClipLength ( SamRecord & sam_rec);

char getClipType ( SamRecord & sam_rec, std::vector<RefSeq*> & REF_SEQ);

void clearDiscReadPair(DiscReadPair * rp);

void setDiscReadPair (DiscReadPair * rp, SamRecord & sam_rec);

// check: 2nd read length, both qual, both chr
bool sanityCheckDiscReadPair(DiscReadPair * rp, SamRecord & sam_rec);

void rescueSingleMCinReadPair( int & mc, std::string & seq, std::vector<RefSeq*> & REF_SEQ );

#endif
