#include "AnchorRank.h"
#include "GLs.h"
#include "Globals.h"
#include <math.h>


// Here all GL or count size should be CHECKED before doing the function!

/* in the functions below, posterior is represented by:
	10*log(x)
	%fraction is:
	100 * (%frac)
	it's a sign difference with PL...
	*/
	
	
// similar to GetVariantQuality but skipped the sanity check, and assume GL[1] or GL[2] > GL[0]	
int GetVariantPosterior( vector<float> & GL)
{
	float lVariant = SumGL( GL[1], GL[2] );
	float lAll = SumGL( GL[0], lVariant );
	float pVariant = GetProbFromGLs( lVariant, lAll );
	int post = round(pVariant * 10);
	return post;
}

float GetSupportReadFraction( vector<int> & counts, int depth )
{
	int supports = getSumSupportClips( counts ) + getSumSupportDiscs( counts ) + getSumSupportUnmaps( counts );
	float frac = float(supports) / depth;
	return frac;
}


int GetModifiedProperReadFraction( vector<int> & counts, int depth )
{
	int frac = counts[0] * 100 / depth;
	return frac;
}


// get supports
int getSumSupportClips( vector<int> & counts )
{
	int sum = counts[4] + counts[5] + counts[8] + counts[9];
	return sum;
}

int getSumSupportDiscs( vector<int> & counts)
{
	int sum = counts[11] + counts[13];
	return sum;
}

int getSumSupportUnmaps( vector<int> & counts)
{
	int sum = counts[15] + counts[17];
	return sum;
}





