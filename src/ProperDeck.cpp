#include "ProperDeck.h"
#include "ProperDeckUtility.h"
#include <iterator>
#include <iostream>
#include <utility> // exit, pair
#include "Globals.h"
#include "ReadMapGlobals.h"

using std::cout;
using std::cerr;
using std::endl;

// constructor
ProperDeck::ProperDeck( vector<RefSeq*> & Ref_Seq ):
	deck_start(0),
	rSeq(Ref_Seq)
{
	properDeck.resize( WIN );
}

// destructor
ProperDeck::~ProperDeck() {}

// principle
/*
	<--front		back-->
	0  1  2  ...	3  4 5
	deck length = WIN
*/ 


// add 1st to deck
void ProperDeck::Add( SamRecord & sam_rec )
{
	int dist = sam_rec.get1BasedPosition() - deck_start;
// locate first
	deque< vector<ProperDeckCell> >::iterator it_deck;
	if (dist >= WIN * 2 - 1) { // start new deck
		deck_start = sam_rec.get1BasedPosition();
		properDeck.clear();
		properDeck.resize(WIN);
		it_deck = properDeck.begin();
	}
	else if ( dist >= WIN ) { // pop & push
		int pop_len = dist - WIN + 1;
		deck_start += pop_len;
		for(int i=0; i < pop_len; i++)
			properDeck.pop_front();
		properDeck.resize( WIN );
		it_deck = properDeck.end();
		it_deck--;
	}
	else { // just add
		it_deck = properDeck.begin();
		it_deck += dist;
	}
// add to deck using it_deck
	it_deck->resize( it_deck->size() + 1 );
	vector<ProperDeckCell>::iterator it_vec = it_deck->end();
	it_vec--;
	it_vec->MatePosition = sam_rec.get1BasedMatePosition();
	it_vec->ReadID = sam_rec.getReadName();
	if ( sam_rec.getMapQuality() < MinQuality || sam_rec.getReadLength() < MinReadLen )
		it_vec->valid = 0;
	else
		it_vec->valid = 1;
	it_vec->ClipType = getSingleClipType(sam_rec);
	it_vec->ClipLen = it_vec->ClipType & 1 ? abs(getMaxClipLen(sam_rec)) : 0; // length include short clip
	it_vec->FirstAlignEnd = sam_rec.get1BasedAlignmentEnd();
}

// get mate info from deck
RetrievedIndex ProperDeck::RetrieveIndex( SamRecord & sam_rec )
{
	RetrievedIndex rIndex;
// find mate
	deque< vector<ProperDeckCell> >::iterator it_deck = properDeck.begin();
	int first_loc = sam_rec.get1BasedMatePosition();
	int deck_dist = first_loc - deck_start;
// some additional sanity check
	if (deck_dist < 0 || deck_dist >= WIN) {
//		if (DEBUG_MODE)
//			std::cerr << "Mate not found in deck (possibly low-qual): " << sam_rec.getReadName() << std::endl;
		rIndex.valid = 0;
		return rIndex;
	}
// then search for mate
	it_deck += deck_dist;
  // if nothing in here
	if ( it_deck->size() == 0 ) {
//		if (DEBUG_MODE)
//			std::cerr << "Nothing in designated position: " << sam_rec.getReadName() << ", possibly boundary of sections?" << std::endl;
		rIndex.valid = 0;
		return rIndex;
	
	}
	
	vector<ProperDeckCell>::iterator it_mate = it_deck->end(); // indicator: if == end, not found
  	vector<int> same_1st_vec;
  // 1st loc first
	for( vector<ProperDeckCell>::iterator pit = it_deck->begin(); pit != it_deck->end(); pit++ ) {
		if ( pit->MatePosition == sam_rec.get1BasedPosition() )
			same_1st_vec.push_back( pit - it_deck->begin() );
	}
	if ( same_1st_vec.size() == 1 ) { // that's the only right one
		it_mate = it_deck->begin();
		it_mate += same_1st_vec[0];
	}
	else if (same_1st_vec.size() > 1) { // compare id
		for( vector<int>::iterator nit = same_1st_vec.begin(); nit != same_1st_vec.end(); nit++ ) {
			vector<ProperDeckCell>::iterator it = it_deck->begin();
			it += *nit;
			if ( it->ReadID.compare( sam_rec.getReadName() ) == 0 )
				it_mate = it;
		}
	}
// adjust rIndex
	if (it_mate == it_deck->end()) {
		rIndex.valid = 0;
		return rIndex;
	}
	if (!it_mate->valid) {
		rIndex.valid = 0;
		return rIndex;
	}
	rIndex.valid = 1;
	
// adjust type & mei_on_1st: it_mate is used here
	bool use2 = 0;
	if (it_mate->ClipType & 1) {
		int clip2_len = getMaxClipLen(sam_rec);
		if (abs(clip2_len) > it_mate->ClipLen) { // use 2
			use2 = 1;
			rIndex.type = getSingleClipType( sam_rec );
			rIndex.mei_on_1st = 0;
		}
		else { // use 1: if it's equal this actually favors clip 1 but who cares...
			rIndex.type = it_mate->ClipType;
			rIndex.mei_on_1st = 1;
		}
	}
	else { // no clip in 1, use 2 (yes, 2 may not have clip either)
		use2 = 1;
		rIndex.type = getSingleClipType( sam_rec );
		rIndex.mei_on_1st = 0;
	}
	
// if proper: min = 1st-start, max = 2nd-end
// if clip( no matter short or sufficient ):
//			  min = clip-read-start  max = clip-read-end (all mean mapped)
	if ( rIndex.type == 0 ) { // proper
		rIndex.max_position = sam_rec.get1BasedAlignmentEnd();
		rIndex.min_position	= first_loc;
	}
	else { // clip
		if (!use2) { // use 1
			rIndex.max_position = it_mate->FirstAlignEnd;
			rIndex.min_position = first_loc;
		}
		else { // use 2
			rIndex.max_position = sam_rec.get1BasedAlignmentEnd();
			rIndex.min_position = sam_rec.get1BasedPosition();
		}
	}
	
	return rIndex;
}


/*** inner funcitons **********/
char ProperDeck::getSingleClipType( SamRecord & sam_rec )
{
	int clip_len = getMaxClipLen(sam_rec);
	if (clip_len == 0)
		return 0; // no clip

// change clip len to abs	
	bool isBeginClip;
	if (clip_len < 0) {
		clip_len *= -1;
		isBeginClip = 0;
	}
	else
		isBeginClip = 1;
		
	if (clip_len < ShortClip)
		return 0; // no clip, just return;
	if (clip_len < MinClip)
		return 1; // short clip
		
// adjust clip_str & isBegin
	char clipType = 3;
	string seq = sam_rec.getSequence();
	string clip_str;
	if (isBeginClip) {
		clipType |= 4;
		clip_str = seq.substr(0, clip_len);
	}
	else
		clip_str = seq.substr(seq.length() - clip_len);

// remap
	int init_i = 4; // 8(2^3) A 16 L 32 S
	for(int i=0; i<3; i++) {
		init_i *= 2;
		if ( rSeq[i]->MeiMap(clip_str) )
			clipType |= init_i;
	}
	return clipType;
}


























