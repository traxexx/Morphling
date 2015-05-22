#include "ProperDeckUtility.h"
#include "Cigar.h"

int getMaxClipLen( SamRecord & sam_rec )
{
	Cigar * myCigar = sam_rec.getCigarInfo();
	int begin_clip = myCigar->getNumBeginClips();
	int end_clip = myCigar->getNumEndClips();
	if (begin_clip >= end_clip)
		return begin_clip;
	else
		return -end_clip;
}

