#include "DiscPair.h"
#include "ReadMapGlobals.h"

// constructor: add 1st in pair
// before doing this already checked if unmapped
DiscPair::DiscPair( bool first_in_candidate, bool second_in_candidate, SamRecord & sam_rec, vector<RefSeq*> & refSeq )
{
	valid = (first_in_candidate || second_in_candidate) ? 1 : 0;
	if (!valid)
		return;
	in_candidate1 = first_in_candidate;
	in_candidate2 = second_in_candidate;
// add others
	chr_name1 = sam_rec.getReferenceName();
	chr_name2 = sam_rec.getMateReferenceName();
	if (in_candidate2)
		seq1 = sam_rec.getSequence();
	loc1 = sam_rec.get1BasedPosition();
	loc2 = sam_rec.get1BasedMatePosition();
	is_map1 = (sam_rec.getFlag() & 0x4) ? 0 : 1;
	is_map2 = (sam_rec.getFlag() & 0x8) ? 0 : 1;
	if (is_map1) {
		high_qual1 = sam_rec.getMapQuality() >= MinQuality ? 1 : 0;
		int * ptr_em = sam_rec.getIntegerTag("EM");
		if (ptr_em)
			em1 = *ptr_em; // em does not care about qual
		else // null pointer: no em tag. This happens in unmap
			em1 = 0;
		if (high_qual1) {
			strand1 = (sam_rec.getFlag() & 0x10) ? 0 : 1;
			align1_end = sam_rec.get1BasedAlignmentEnd();
		}
		else
			in_candidate1 = 0;
	}
	rSeq = refSeq;
}


// destructor
DiscPair::~DiscPair() {}

// check if 2nd the same pair as this DiscPair (do not check read name)
bool DiscPair::IsSamePair( SamRecord & sam_rec )
{
	if ( chr_name1.compare( sam_rec.getMateReferenceName() ) != 0 )
		return 0;
	if ( chr_name2.compare( sam_rec.getReferenceName() ) != 0 )
		return 0;
	if ( loc1 != sam_rec.get1BasedMatePosition() )
		return 0;
	if ( loc2 != sam_rec.get1BasedPosition() )
		return 0;
// pass all checkes
	return 1;
}

// add second to pair if pass IsSamePair
void DiscPair::AddSecondToPair( SamRecord & sam_rec )
{
	if (!valid)
		return;
// set others
	if (is_map2) {
		int * em_ptr = sam_rec.getIntegerTag("EM");
		if ( !em_ptr )
			em2 = 0;
		else
			em2 = *em_ptr;
		high_qual2 = sam_rec.getMapQuality() >= MinQuality ? 1 : 0;
		if (high_qual2) {
			strand2 = (sam_rec.getFlag() & 0x10) ? 0 : 1;
			align2_end = sam_rec.get1BasedAlignmentEnd();
		}
		else
			in_candidate2 = 0;
	}
	if (in_candidate1)
		seq2 = sam_rec.getSequence();
// re-adjust	
	if ( !in_candidate1 && !in_candidate2) {
		valid = 0;
		return;
	}
  // re-map (rescue): 1 2 4 A L S
	if (in_candidate2) {
		int init_i = 1;
		for(int i=0; i<3; i++) {
			if (!(em1 & init_i)) {
				if ( rSeq[i]->MeiMap(seq1) )
					em1 |= init_i;
			}
			init_i *= 2;
		}
	}
	if (in_candidate1) {
		int init_i = 1;
		for(int i=0; i<3; i++) {
			if (!(em2 & init_i)) {
				if ( rSeq[i]->MeiMap(seq2) )
					em2 |= init_i;
			}
			init_i *= 2;
		}
	}
}


// see if 1st can be used as anchor
bool DiscPair::FirstValid()
{
	if (in_candidate1)
		return 1;
	else
		return 0;
}


// see if 2nd can be used as anchor
bool DiscPair::SecondValid()
{
	if (in_candidate2)
		return 1;
	else
		return 0;
}


int DiscPair::GetFirstAlignPosition()
{
	return loc1;
}

int DiscPair::GetSecondAlignPosition()
{
	return loc2;
}


// get 1st loci for adding to discmap
// 1st anchor
// start = alignment end
// end = alignment_end + 1/2 * read length
// 2nd anchor
// start = alignment start - 1/2 * read length
// end = alignment start
Loci DiscPair::GetFirstLoci()
{
	Loci loci;
	loci.chr = chr_name1;
	loci.sense_strand = strand1;
	if ( strand1 ) {
		loci.st = align1_end;
		loci.ed = align1_end + seq1.length() * 0.5;
	}
	else {
		loci.st = loc1 - seq1.length() * 0.5;
		loci.ed = loc1;
	}
	if ( loci.st < 0 )
		loci.st = 0;
	loci.type = em2;
	loci.mateMap = is_map2;
	return loci;
}


// 2nd loci
Loci DiscPair::GetSecondLoci()
{
	Loci loci;
	loci.chr = chr_name2;
	loci.sense_strand = strand2;
	if (strand2) {
		loci.st = align2_end;
		loci.ed = align2_end + seq2.length() * 0.5;
	}
	else {
		loci.st = loc2 - seq2.length() * 0.5;
		loci.ed = loc2;
	}
	if ( loci.st < 0 )
		loci.st = 0;
	loci.type = em1;
	loci.mateMap = is_map1;
	return loci;
}










