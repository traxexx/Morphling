#include "glCalc.h"
#include "morphError.h"
#include <math.h>
#include <algorithm>

// note gl is log10(gl)
int getGtFromGl( vector<float> & log10_prior, vector<float> & gl )
{
	if (gl.empty()) // missing genotype
		return -1;
	if (gl.size() != 3)
		morphError("[getGtFromGl] gl size not 3");
	vector<float> lpost;
	lpost.resize(3);
	for( int i=0; i<3; i++ )
		lpost[i] = log10_prior[i] + gl[i];
	vector<float>::iterator pmax = std::max_element(lpost.begin(), lpost.end());
	int gt = pmax - lpost.begin();
	if (gt < 0 || gt > 2)
		morphError("[getGtFromGl] abnormal gt");
	return gt;
}

// EM based on gl
// return ref allele frequency
bool setAFfromGL( vector<float> & gf, vector<lhRec> & vl )
{
	gf.resize(3);
	int iter = 0;
	gf[0] = 1.0 / 3;
	gf[1] = gf[0];
	gf[2] = gf[1];
	vector<float> log10_gf;
	log10_gf.resize(3);
	for(int i=0; i<3; i++)
		log10_gf[i] = log10(1.0/3);
	float mse = 1;
	float diff = 0;
	int n = 0;
	for(int i=0; i<vl.size(); i++) {
		if (!vl[i].gl.empty())
			n++;
	}
	if (n==0)
		return 0; // unable to set AF because no gl
//		morphError( "[setLog10AF] no info at site" );
	
	while( mse > 1e-3 && iter < 50 ) {
		float mle_gf[3] = {0, 0, 0};
		vector<float> lpost (3,0);
		mse = 0;
	// for each sample
		for( int i = 0; i < vl.size(); ++i) {
			if ( vl[i].gl.empty()) // skip missing
				continue;
			for( int j = 0; j < 3; j++ )
				lpost[j] = log10_gf[j] + vl[i].gl[j];
			// sum & normalize
			convertToProb( lpost );
			for( int j = 0; j < 3; j++ )
				mle_gf[j] += lpost[j];
		}	
		for( int j = 0; j < 3; j++ )
			mle_gf[j] = mle_gf[j] / n;
		for( int j = 0; j < 3; j++ ) {
			diff = gf[j] - mle_gf[j];
			mse += ( diff * diff );
		}
		for( int j = 0; j < 3; j++ ) {
			gf[j] = mle_gf[j];
			log10_gf[j] = log10(gf[j]);
		}
		iter++;
	}
// final check
	if ( iter >= 50 ) {
		string str = "[ConsensusVcfRecord::EstimateAFviaEM()] EM did not converge. mse = " + std::to_string(mse);
		morphWarning(str);
	}
	return 1;
}


// log10 to prob, normalized
void convertToProb( vector<float>  & lp )
{
	float sum = 0;
	vector<float>::iterator pm = std::max_element( lp.begin(), lp.end() ); // because it's negative
	int max_val = *pm;
	for(int i=0; i<3; i++) {
		lp[i] -= (max_val); // scale
		lp[i] = pow(10, lp[i]); // back to prob
		sum += lp[i]; // get sum
	}
	for(int i=0; i<3; i++)
		lp[i] /= sum; // normalize
}
