#include "ConsensusVcfRecord.h"
#include "Utilities.h"
#include <iostream>
#include <iterator>
#include <algorithm> // all_of
#include <sstream>
#include <math.h>
#include <map>
#include "GLs.h" // GetDosageFromGenotype

using std::to_string;
using std::cout;
using std::cerr;
using std::endl;
using std::map;

ConsensusVcfRecord::ConsensusVcfRecord( int position )
{
	this->position = position;
}

ConsensusVcfRecord::~ConsensusVcfRecord() {}

bool ConsensusVcfRecord::CheckPosition( int & location )
{
	if ( location == position )
		return 1;
	else
		return 0;
}

int ConsensusVcfRecord::GetPosition()
{
	return this->position;
}

int ConsensusVcfRecord::GetRawVariantQuality()
{
	if ( variant_quality <= 0 )
		estimateRawVariantQuality();
	return variant_quality;
}


int ConsensusVcfRecord::GetSumDosage()
{
	int sum = 0;
	for( vector<int>::iterator ptr = this->Dosages.begin(); ptr != this->Dosages.end(); ptr++ ) {
		if ( *ptr >= 0 )
			sum += *ptr;
	}
	return sum;
}



void ConsensusVcfRecord::SetSampleSize( int ns )
{
	PLs.resize( ns );
	DPs.resize( ns, 0 );
	GQs.resize( ns, 0 );
	Dosages.resize( ns, 0 );
}

void ConsensusVcfRecord::SetAltAllele( string & alt )
{
	this->alt_allele = alt;
}

void ConsensusVcfRecord::SetInfoFields( int sp_index, string & info_str )
{
	std::stringstream ss;
	ss << info_str;
	map< string, string > info_fields;
	string field;
	while( getline( ss, field, ';' ) ) {
		std::stringstream flss;
		flss << field;
		string sub1;
		string sub2;
		getline( flss, sub1, '=' );
		getline( flss, sub2, '=' );
		info_fields[ sub1 ] = sub2;
	}
	
	map< string, string >::iterator ptr;
// average variant end
	ptr = info_fields.find( "END" );
	if ( ptr != info_fields.end() ) {
		int current_end = stoi(ptr->second);
		updateVariantEnd( current_end );
		updateCIEND( current_end );
	}
// see if this break is ci_low or ci_low
	ptr = info_fields.find( "BREAKP" );
	if ( ptr != info_fields.end() ) {
		updateCIPOS( stoi(ptr->second) );
	}
}

// field must be GT:DP:GQ:GL
void ConsensusVcfRecord::SetGLFields( int sp_index, string & gl_str )
{
	std::stringstream ss;
	ss << gl_str;
	vector< string > items;
	string field;
	while( getline(ss, field, ':') )
		items.push_back( field );
	if ( items.size() != 4 )
		cerr << "Warning: abnormal GL field: " << gl_str << " at position: " << position << endl;
		
	// parse field
	Dosages[ sp_index ] = GetDosageFromGenotype( items[0] );
	DPs[ sp_index ] = stoi( items[1] );
	GQs[ sp_index ] = stoi( items[2] );
	vector< string > gl_items;
	std::stringstream glss;
	glss << items[3];
	while( getline( glss, field, ',' ) )
		gl_items.push_back( field );
	if ( gl_items.size() != 3 )
		cerr << "Warning: abnormal GL sub field: " << items[3] << " at position: " << position << endl;
	for( int i = 0; i <= 2; i++ )
		PLs[ sp_index ][i] = stoi( gl_items[i] );
// normalize
	int norm = PLs[sp_index][ Dosages[ sp_index ] ];
	if ( norm != 0 ) {
		for( int i = 0; i <= 2; i++ )
			PLs[ sp_index ][i] -= norm;
	}
}

void ConsensusVcfRecord::Print( string & chr, ofstream & out_vcf )
{
// update all empty fields to dot
	out_vcf << chr << "\t" << position << "\t.\t" << ref_allele << "\t" << alt_allele;
	string info_field = generateInfoStr();
	int prvq = getVarintQualityWithAF();
	out_vcf << "\t" << prvq << "\t" << filter << "\t" << info_field << "\t";
  // these fields shouldn't be empty
	out_vcf << "GT:DP:GQ:PL";
	for( int sp = 0; sp < (int)DPs.size(); sp++ ) {
		string format = generateGLStr( sp );
		out_vcf << "\t" << format;
	}
	out_vcf << endl;
}


// convert PL to probability. Let sum(probs) = 1
void ConsensusVcfRecord::setProbViaPLs( vector< vector<float> > prob_gls )
{
	int n = (int)PLs.size();
	prob_gls.resize( n );
	for( int i = 0; i < n; i++ ) {
		prob_gls[i].resize(3);
		float sum = 0;
		for( int j = 0; j < 3; j++ ) {
			prob_gls[i][j] = pow( 10, -(float)PLs[i][j] / 10 );
			sum += prob_gls[i][j];
		}
		for( int j = 0; j < 3; j++ )
			prob_gls[i][j] /= sum;
	}
}

/*** set methods used to polish record ***/
void ConsensusVcfRecord::SetRefAllele( GenomeSequence * gs, string & chr_name )
{
	ref_allele = gs->getBase( chr_name.c_str(), this->position );
}

void ConsensusVcfRecord::EstimateFIC()
{
	float num = 0;
	float denum = 0;
	vector< vector<float> > prob_gls;
	setProbViaPLs( prob_gls );
	float gf[3];
	gf[0] = gt_freq[0];
	gf[1] = gt_freq[1];
	gf[2] = 1 - gf[0] - gf[1];
	float hwe[3];
	hwe[0] = ref_af * ref_af;
	hwe[1] = ( 1 - ref_af ) * ref_af;
	hwe[2] = ( 1 - ref_af ) * ( 1 - ref_af );
	
	for( int sp = 0; sp < (int)Dosages.size(); sp++ ) {
		if ( DPs[sp] <= 0 ) // skip missing
			continue;
		float o_het_sum = gf[1] * prob_gls[sp][1];
		float o_sum = gf[0] * prob_gls[sp][0];
		o_sum += o_het_sum;
		o_sum += gf[2] * prob_gls[sp][2];
		
		float e_het_sum = hwe[1] * prob_gls[sp][1];
		float e_sum = hwe[0] * prob_gls[sp][0];
		e_sum += e_het_sum;
		e_sum += hwe[2] * prob_gls[sp][2];
		
		num += o_het_sum / o_sum;
		denum += e_het_sum / e_sum;
	}
	if ( denum == 0 ) {
		cerr << "Warning: unable to estimate FIC with missing PL info at position: " << position << endl;
		fic = 2;
	}
	else
		fic = 1 - num / denum;
}

void ConsensusVcfRecord::RefineInfoFields()
{

}

void ConsensusVcfRecord::SetFilter()
{
	filter = "PASS";
}


void ConsensusVcfRecord::EstimateAFviaEM()
{
	int iter = 0;
	float gf[3];
	gf[0] = 1.0 / 3;
	gf[1] = gf[0];
	gf[2] = gf[1];
	float mse = 0;
	float diff = 0;
	int n = 0;
	for( vector< int >::iterator it = DPs.begin(); it != DPs.end(); it++ ) {
		if( *it > 0 )
			n++;
	}
	if ( n == 0 ) {
		cerr << "Warning: no valid info at position: " << position << endl;
		ref_af = -1;
		return;
	}
	vector< vector<float> > prob_gls;
	setProbViaPLs( prob_gls );
	
	while( mse > 1e-3 && iter < 50 ) {
		float mle_gf[3] = {0, 0, 0};
		float gf_indi[3];
	// for each sample
		for( int i = 0; i < n; ++i) {
			if ( DPs[i] <= 0 ) // skip missing
				continue;
			for( int j = 0; j < 3; j++ )
				gf_indi[j] = gf[0] * prob_gls[i][j];
			float prob_sum = gf_indi[0] + gf_indi[1] + gf_indi[2];
			for( int j = 0; j < 3; j++ )
				gf_indi[j] /= prob_sum;
			for( int j = 0; j < 3; j++ )
				mle_gf[j] += gf_indi[j];
		}	
		for( int j = 0; j < 3; j++ )
			mle_gf[j] = mle_gf[j] / n;
		for( int j = 0; j < 3; j++ ) {
			diff = gf[j] - mle_gf[j];
			mse += ( diff * diff );
		}
		for( int j = 0; j < 3; j++ )
			gf[j] = mle_gf[j];
		iter++;
	}
// final check
	if ( iter >= 50 )
		cerr << "Warning: EM not converge..." << endl;
	ref_af = gf[0] + gf[1] / 2;
}


string ConsensusVcfRecord::generateInfoStr()
{
	string info;
	int nv = countSampleWithVariant();
	int variant_end;
	if ( nv == 0 ) {
		cerr << "Warning: no variant at position: " << this->position << endl;
		variant_end = 0;
	}
	else
		variant_end = round( (float)variant_end_sum / nv );
	int ci_low = left_most - position;
	int ci_high = right_most - position;
	int ci_end_low = left_most_end - variant_end;
	int ci_end_high = right_most_end - variant_end;
	info = string("SVTYPE=INS;END=") + to_string(variant_end) + ";CIPOS=" + to_string(ci_low) + ", " + to_string(ci_high);
	info += ";CIEND=" + to_string(ci_end_low) + "," + to_string(ci_end_high) + ";SVLEN=" + to_string(sv_length);
	float raf = round( ref_af * 10000 ) * 0.0001;
	float rfic = round( fic * 10000 ) * 0.0001;
	info += ";AF=" + to_string(raf) + ";FIC=" + to_string(rfic); 
	return info;
}


string ConsensusVcfRecord::generateGLStr( int sp )
{
	string format = getGenotypeFromDosage( Dosages[sp] );
	format += ":" + to_string( DPs[sp] ) + ":" + to_string( GQs[sp] ) + ":";
	format += to_string(PLs[sp][0]) + "," + to_string(PLs[sp][1]) + "," + to_string(PLs[sp][2]);
	return format;
}

string ConsensusVcfRecord::getGenotypeFromDosage( int dosage )
{
	string gt;
	if ( dosage == 0 )
		gt = "0/0";
	else if ( dosage == 1 )
		gt = "0/1";
	else if ( dosage == 2 )
		gt = "1/1";
	else
		gt = "./.";
	return gt;
}


void ConsensusVcfRecord::estimateRawVariantQuality()
{
	int sum = 0;
	for( vector<vector<int>>::iterator ptr = this->PLs.begin(); ptr != this->PLs.end(); ptr++ ) {
		if ( (*ptr)[0] >= 0 )
			sum += (*ptr)[0];
	}
	this->variant_quality = sum;
}


int ConsensusVcfRecord::getVarintQualityWithAF()
{
	if ( variant_quality <= 0 )
		estimateRawVariantQuality();
	int ns = (int)DPs.size();
	int prvq = variant_quality - 10*ns*log10( gt_freq[0] );
	return prvq;
}

void ConsensusVcfRecord::updateCIPOS( int pos )
{
	if ( pos < this->left_most )
		this->left_most = pos;
	if ( pos > this->right_most )
		this->right_most = pos;
}

void ConsensusVcfRecord::updateVariantEnd( int current_end )
{
	variant_end_sum += current_end;
}


void ConsensusVcfRecord::updateCIEND( int current_end )
{
	if ( current_end < this->left_most_end )
		this->left_most_end = current_end;
	if ( current_end > this->right_most_end )
		this->right_most_end = current_end;
}

int ConsensusVcfRecord::countSampleWithVariant()
{
	int nv = 0;
	for( vector<int>::iterator ptr = Dosages.begin(); ptr != Dosages.end(); ptr++ ) {
		if ( *ptr > 0 )
			nv++;
	}
	return nv;
}



