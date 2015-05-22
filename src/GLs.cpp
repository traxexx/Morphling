#include "GLs.h"

#include <iostream>
#include <iterator>
#include <math.h> // exp, log, round

using std::cerr;
using std::endl;

// summing all together
float GetGLfromCounts( vector<int> & counts, vector<float> & ref )
{
	float gl = 0;
	vector<float>::iterator ref_it = ref.begin();
	for( vector<int>::iterator data = counts.begin(); data!= counts.end(); data++, ref_it++ ) {
		gl += (*data) * (*ref_it);
	}
	return gl;
}


// sum 2 log-p
float SumGL( float original, float single )
{
	if ( original - single >= 14 )
		return original;
	if ( original - single <= -14 )
		return single;
// do sum
	float mid = ( original + single ) / 2;
	float conjugate = mid + log( exp(original - mid) + exp(single-mid) );
	return conjugate;
}


// minus a log-p from conjugate. If compensate > original, return 1;
float MinusGL( float original, float compensate)
{	
	if ( original - compensate < 0 ) { // need to re-calculate outside
		cerr << "ERROR: compensate - original > 5. Cannot calculate minusGL!" << endl;
		exit(1);
	}
	if ( original - compensate >= 7) { // no compensate
		return original;
	}
	float mid = ( original + compensate ) / 2;
	float conjugate = mid + log( exp(original - mid) - exp(compensate-mid) );
	return conjugate;
}

// sum of vector
int GetVecSum( vector<int> & counts )
{
	int sum = 0;
	for( vector<int>::iterator it = counts.begin(); it != counts.end(); it++ )
		sum += (*it);
	return sum;
}

// # zeros in vec
int GetNumberOfZerosInVec( vector<int> & counts )
{
	int sum = 0;
	for( vector<int>::iterator it = counts.begin(); it != counts.end(); it++ ) {
		if ( *it == 0 )
			sum++;
	}
	return sum;
}

/**** PL & vcf related ***/

// Prob is exp(BaseGL) / exp(AddGL);
float GetProbFromGLs( float BaseGL, float AddGL )
{
	float prob;
	if ( BaseGL - AddGL >= 14 )
		prob = 0.999999;
	else if ( AddGL - BaseGL >= 14 )
		prob = 0.000001;
	else {
		float mid = ( BaseGL + AddGL ) / 2;
		prob = exp( BaseGL - mid ) / exp( AddGL - mid );
	}
	return prob;
}

// pl of p(no variant) / p(all)
int GetVariantQualityFromGL( vector<float> & GL )
{
	int qual;
	if (GL.size() != 3) {
		cerr << "ERROR: GL size != 3 at GetGenotype() " << endl;
		exit(1);
	}
	
// start
	float lVariant = SumGL( GL[1], GL[2] );
	float lAll = SumGL( GL[0], lVariant );
	float pVariant = GetProbFromGLs( lVariant, lAll );
	if ( GL[0] >= GL[1] && GL[0] >= GL[2] ) { // no variant, return -10 * log10 ( variant )
		qual = round( -10 * log10( pVariant ) );
	}
	else { // exist variant, return p( no variant )
		if ( pVariant >= 0.999999 )
			qual = 60; // 10E(-6) as minimum
		else
			qual = round( -10 * log10( 1 - pVariant ) );
	}
	return qual;
}


int GetAlleleDosage( vector<float> & GL )
{
	int dosage;
	if ( GL[0] >= GL[1] ) {
		if ( GL[0] >= GL[2] )
			dosage = 0;
		else // 0 < 2
			dosage = 2;
	}
	else { // 1 > 0
		if ( GL[1] >= GL[2] )
			dosage = 1;
		else // 2 > 1
			dosage = 2;
	}
	return dosage;
}


string GetGenotype( vector<float> & GL )
{
	if (GL.size() != 3) {
		cerr << "ERROR: GL size != 3 at GetGenotype() " << endl;
		exit(1);
	}
	string gt;
	if ( GL[0] >= GL[1] ) {
		if ( GL[0] >= GL[2] )
			gt = string( "0/0" );
		else // 0 < 2
			gt = string( "1/1" );
	}
	else { // 1 > 0
		if ( GL[1] >= GL[2] )
			gt = string( "0/1" );
		else // 2 > 1
			gt = string( "1/1" );
	}
	return gt;
}

void SetPLsFromGL( vector<int> & PL, vector<float> & GL )
{
	float base10 = log(10);
	if ( !PL.empty() ) {
		cerr << "ERROR: PL is required to be empty in SetPLsFromGL!" << endl;
		exit(1);
	}
	PL.resize(3);
	if (GL.size() != 3) {
		cerr << "ERROR: GL size != 3 at SetPLsFromGL! " << endl;
		exit(1);
	}
	for( int i = 0; i <= 2; i++ ) // maybe PL is not that precise due to GL is stored as int
		PL[i] = round( -GL[i] * 10 / base10 );
}


// gt quality = wrong varint gt / all variant gt
int GetGenotypeQuality( vector<float> & GL )
{
	int gq;
	float sum =  SumGL( GL[1], GL[2] );
	if ( sum == GL[1] || sum == GL[2] ) // variant GL is dominating
		gq = 60;
	else { // not dominating
		float p1 = GetProbFromGLs( GL[1], sum );
		if ( GL[1] >= GL[2] ) {// 1 is right gt
			if ( p1 >= 0.999999 )
				gq = 60;
			else
				gq = round( (-10) * log10( 1 - p1 ) );
		}
		else // 2 is right gt
			gq = round( (-10) * log10( p1) );
	}
	return gq;
}

int GetDosageFromGenotype( string & gt )
{
	int ds = 0;
	if ( gt[0] == '.' ) // NA
		ds = -1;
	else if ( gt[0] == '1' )
		ds = 1;
		
	if ( ds < 0 ) // if one allele is NA, no matter what the other is, return NA
		return ds;
		
	if ( gt[2] == '.' )
		ds = -1;
	else if ( gt[2] == '1' )
		ds++;
		
	return ds;	
}

























