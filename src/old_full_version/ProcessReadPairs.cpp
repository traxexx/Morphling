#include "ProcessReadPairs.h"
#include "Params.h"

// OQ
void simplifySamRec( SamRecord & sam_rec )
{
	bool status = sam_rec.rmTag("OQ", 'Z');
	if (!status)
		std::cerr << "Warning: OQ tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("XA", 'Z');
	if (!status)
		std::cerr << "Warning: XA tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("MD", 'Z');
	if (!status)
		std::cerr << "Warning: MD tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("RG", 'Z');
	if (!status)
		std::cerr << "Warning: RG tag not removed in " << sam_rec.getReadName() << std::endl;			
}

// 1 2 3 A L S
void setEMtag( SamRecord & sam_rec, std::vector<CordPtr> & CordPtrVec, std::vector<CordPtr> & EndPtrVec )
{
	if (sam_rec.getFlag() & 0x4 || sam_rec.getFlag() & 0x8) {
		sam_rec.addIntTag("EM", 0);
		return; // nothing else to do
	}
// check if inside MeiCoord
	int tag = 0;
	int align_tag = 1;
	for( int i=0; i<3; i++ ) {
		if (i != 0)
			align_tag *= 2;
		CordPtr cp = CordPtrVec[i];
		if (cp == EndPtrVec[i])  continue;
		while (sam_rec.get1BasedUnclippedEnd() <= cp->start) {// should be unclipped here for further reads
			cp++;
			if (cp == EndPtrVec[i])  continue;
		}
		CordPtrVec[i] = cp;
		if (sam_rec.get1BasedPosition() <= cp->end && sam_rec.get1BasedAlignmentEnd() >= cp->start)
			tag |= align_tag;
	}
	sam_rec.addIntTag("EM", tag);
}


int getMaxClipLength( SamRecord & sam_rec )
{
	Cigar * myCigar = sam_rec.getCigarInfo();
	int beginClip = myCigar->getNumBeginClips();
	int endClip = myCigar->getNumEndClips();
	if (beginClip >= endClip) {
		return beginClip;
	}
	else {
		return -endClip;
	}	
}


char getClipType (SamRecord & sam_rec, std::vector<RefSeq*> & REF_SEQ)
{
	char clipType = 0;
	bool isBeginClip = 0;
	int maxClipLen = getMaxClipLength(sam_rec);
	if (maxClipLen < 0) {
		maxClipLen *= (-1);
	}
	else
		isBeginClip = 1;
	
	
	if (maxClipLen >= ShortClip) {
		clipType |= 1;
		if (maxClipLen >= MinClip) {
			clipType |= 2;
			std::string clipString;
			std::string seq = sam_rec.getSequence();
			if (isBeginClip) {
				clipType |= 4;
				clipString = seq.substr(0,maxClipLen);
			}
			else {
				clipString = seq.substr(seq.length() - maxClipLen);
			}
			if (REF_SEQ[0]->MeiMap(clipString))
				clipType |= 8;
			if (REF_SEQ[1]->MeiMap(clipString))
				clipType |= 16;
			if (REF_SEQ[2]->MeiMap(clipString))
				clipType |= 32;
		}			
	}
	return clipType;
}

void clearDiscReadPair(DiscReadPair * rp)
{
	rp->is_empty = 1;
	rp->id.clear();
	rp->seq1.clear();
	rp->seq2.clear();
}

void setDiscReadPair (DiscReadPair * rp, SamRecord & sam_rec)
{
	rp->seq1 = sam_rec.getSequence();
	if (rp->seq1.length() >= MinReadLength) {
		rp->is_empty = 0;
		rp->id = sam_rec.getReadName();
		rp->first_high_qual = sam_rec.getMapQuality() >= MinQuality ? 1 : 0;
		if (rp->first_high_qual)
			rp->first_alignment_end = sam_rec.get1BasedAlignmentEnd();
		rp->mc1 = *sam_rec.getIntegerTag("EM");
	}
	else
		rp->seq1.clear();
}


// 2nd read length, unmap, both qual
bool sanityCheckDiscReadPair(DiscReadPair * rp, SamRecord & sam_rec)
{
	rp->seq2 = sam_rec.getSequence();
// 2nd length
	if (rp->seq2.length() < MinReadLength)
		return 0;
// unmap + quality check
	rp->second_high_qual = sam_rec.getMapQuality() >= MinQuality ? 1 : 0;
// check if unmap anchor low qual
	int flag = sam_rec.getFlag();
	rp->first_map = (flag & 0x8) ? 0 : 1;
	rp->second_map = (flag & 0x4) ? 0 : 1;
	
	if ( (!rp->first_map) & (!rp->second_high_qual))
		return 0;
	if ( (!rp->second_map) & (!rp->first_high_qual))
		return 0;
// check both quality	
	if (!(rp->first_high_qual | rp->second_high_qual)) // both low qual
		return 0;
	
	rp->mc2 = *sam_rec.getIntegerTag("EM");
// all pass
	return 1;	
}


void rescueSingleMCinReadPair( int & mc, std::string & seq, std::vector<RefSeq*> & REF_SEQ ) {
	int tag_val = 1;
	for(int i=0; i<3; i++) {
		if (i != 0)
			tag_val *= 2;
		if (mc & tag_val)  continue;
		if (REF_SEQ[i]->MeiMap(seq))
			mc |= tag_val;
	}
}

















