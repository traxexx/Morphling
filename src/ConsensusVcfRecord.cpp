#include "ConsensusVcfRecord.h"
#include "Utilities.h"
#include <iostream>
#include <iterator>
#include <algorithm> // all_of
#include <sstream>
#include <math.h>
#include <map>
#include <iomanip> // setprecision

using std::to_string;
using std::cout;
using std::cerr;
using std::endl;
using std::map;

ConsensusVcfRecord::ConsensusVcfRecord( int position )
{
	this->position = position;
	this->ref_allele = '.';
	this->variant_quality = 0;
	this->filter.clear();
	this->gt_freq[0] = 0;
	this->gt_freq[1] = 0;
	this->Dosages.clear();
	this->DPs.clear();
	this->GQs.clear();
	this->GLs.clear();
	this->lanchor = 0;
	this->ranchor = 0;
	this->variant_end_baseline = position;
	this->variant_end_diff_sum = 0;
	this->left_most = 0;
	this->right_most = 0;
	this->left_most_end = 0;
	this->right_most_end = 0;
	this->evidence = 0;
	this->wcount = 0;
	this->ref_af = -1;
	this->hwe = -1;
	this->fic = 2;
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


int ConsensusVcfRecord::GetEvidence()
{
	return evidence;
}

int ConsensusVcfRecord::GetWinCount()
{
	return wcount;
}

void ConsensusVcfRecord::SetSampleSize( int ns )
{
	GLs.resize( ns );
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

// average variant end
	map< string, string >::iterator end_ptr = info_fields.find( "END" );
	if ( end_ptr != info_fields.end() ) { // it is a variant
	// update end
		int current_end = stoi(end_ptr->second);
		updateVariantEnd( current_end );
		updateCIEND( current_end );
	// update position
		map< string, string >::iterator ptr = info_fields.find( "Breakp" );
		if ( ptr != info_fields.end() )
			updateCIPOS( stoi(ptr->second) );
		else
			updateCIPOS( this->position );
	// update evidence
		if ( info_fields.find("CLIP") == info_fields.end() || info_fields.find("DISC") == info_fields.end() || info_fields.find("UNMAP") == info_fields.end()) {
			cerr << "ERROR: [ConsensusVcfRecord::SetInfoFields] can't find CLIP/DISC/UNMAP in INFO field. Are you using vcf from Morphling Discover output?" << endl;
			exit(1);
		}
		evidence = stoi(info_fields["CLIP"]) + stoi(info_fields["DISC"]) + stoi(info_fields["UNMAP"]);
	// update wcount
		if ( info_fields.find("WCOUNT") == info_fields.end() ) {
			cerr << "ERROR: [ConsensusVcfRecord::SetInfoFields] can't find WCOUNT in INFO field. Are you using vcf from Morphling Discover output?" << endl;
			exit(1);
		}
		wcount = stoi( info_fields["WCOUNT"] );
	// update left & right anchor
		if ( info_fields.find("ANC") == info_fields.end() ) {
			cerr << "ERROR: [ConsensusVcfRecord::SetInfoFields] can't find ANC in INFO field. Are you using vcf from Morphling Discover output?" << endl;
			exit(1);
		}
		std::stringstream cbss;
		cbss << info_fields["ANC"];
		getline( cbss, field, ',' );
		this->lanchor += stoi(field);
		getline( cbss, field, ',' );
		this->ranchor += stoi(field);
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
	if ( sp_index < 0 || sp_index >= (int)Dosages.size() ) {
		cerr << "ERROR: [ConsensusVcfRecord::SetGLFields] sp_index = " << sp_index << ", but Dosages size = " << (int)Dosages.size() << endl;
		exit(1);
	}
	Dosages[ sp_index ] = getDosageFromGenotype( items[0] );
	DPs[ sp_index ] = stoi( items[1] );
	int gq = stoi( items[2] );
	GQs[ sp_index ] = gq;
	vector< string > gl_items;
	std::stringstream glss;
	glss << items[3];
	while( getline( glss, field, ',' ) )
		gl_items.push_back( field );
	if ( gl_items.size() != 3 )
		cerr << "Warning: abnormal GL sub field: " << items[3] << " at position: " << position << endl;
	if ( GLs[sp_index].empty() )
		GLs[sp_index].resize(3);
	for( int i = 0; i <= 2; i++ )
		GLs[ sp_index ][i] = stof( gl_items[i] );
// normalize
	float norm = GLs[sp_index][ Dosages[ sp_index ] ];
	if ( norm != 0 ) {
		for( int i = 0; i <= 2; i++ )
			GLs[ sp_index ][i] -= norm;
	}
}

void ConsensusVcfRecord::Print( string & chr, ofstream & out_vcf )
{
// update all empty fields to dot
	out_vcf << chr << "\t" << position << "\t.\t" << ref_allele << "\t" << alt_allele <<  "\t" << getVarintQualityWithAF() << "\t" << filter << "\t";
	printInfoStr( out_vcf );
  // these fields shouldn't be empty
	out_vcf << "\tGT:DP:GQ:PL";
	for( int sp = 0; sp < (int)DPs.size(); sp++ ) {
		out_vcf << "\t";
		printGLStr( sp, out_vcf );
	}
	out_vcf << endl;
}


// convert PL to probability. Let sum(probs) = 1
void ConsensusVcfRecord::setProbViaGLs( vector< vector<float> > & prob_gls )
{
	int n = (int)GLs.size();
	prob_gls.resize( n );
	for( int i = 0; i < n; i++ ) {
		prob_gls[i].resize(3);
		float sum = 0;
		for( int j = 0; j < 3; j++ ) {
			prob_gls[i][j] = exp( GLs[i][j] );
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

// SR>10 || <0.1 --> SR
// DP vote --> DP
void ConsensusVcfRecord::SetFilter( vector<float> & dps, vector<int> & ref_dp_cuts, vector<int> & alt_dp_cuts )
{
// SR
	filter.clear();
	if ( lanchor == 0 || ranchor == 0 )
		filter = "SR";
	else {
		float sratio = (float)lanchor / ranchor;
		if ( sratio>10 || sratio <0.1 )
			filter = "SR";
	}

// DP
	int dpvote = 0;
	for( int i=0; i<(int)Dosages.size(); i++ ) {
		if ( Dosages[i] == 0 ) { // ref
			if ( DPs[i] < ref_dp_cuts[i]/4 || DPs[i] > ref_dp_cuts[i]*3 )
				dpvote -= dps[i];
			else
				dpvote += dps[i];
		}
		else {
			if ( DPs[i] < alt_dp_cuts[i]/4 || DPs[i] > alt_dp_cuts[i]*3 ) { // see if het fail
				if ( DPs[i] < (alt_dp_cuts[i]+ref_dp_cuts[i])/8 || DPs[i] > (alt_dp_cuts[i]+ref_dp_cuts[i])*3/2 )
					dpvote -= dps[i];
				else
					dpvote += dps[i];
			}
			else
				dpvote += dps[i];
		}
	}
	if (dpvote<0) {
		if (!filter.empty())
			filter += '+';
		filter += "DEPTH";
	}

// final
	if ( filter.empty() )	
		filter = "PASS";
}

// also estimate HWE test chi-square & FIC
void ConsensusVcfRecord::EstimateAFviaEM()
{
	int iter = 0;
	float gf[3];
	gf[0] = 1.0 / 3;
	gf[1] = gf[0];
	gf[2] = gf[1];
	float mse = 1;
	float diff = 0;
	int n = 0;
	for( vector< int >::iterator it = DPs.begin(); it != DPs.end(); it++ ) {
		if( *it > 0 )
			n++;
	}
	if ( n == 0 ) {
		cerr << "Warning: [ConsensusVcfRecord::EstimateAFviaEM()] No info at position: " << position << endl;
		ref_af = -1;
		return;
	}
	vector< vector<float> > prob_gls;
	setProbViaGLs( prob_gls );
	
	while( mse > 1e-3 && iter < 50 ) {
		float mle_gf[3] = {0, 0, 0};
		float gf_indi[3];
		mse = 0;
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
		cerr << "Warning: [ConsensusVcfRecord::EstimateAFviaEM()] EM did not converge. mse = " << mse << endl;
	for( int j=0; j<2; j++ ) {
		if( gf[j] == 0 )
			gf[j] = (float)1/(int)Dosages.size()/10 > 0.000001 ? (float)1/(int)Dosages.size()/10 : 0.000001;
	}
	this->ref_af = gf[0] + gf[1] / 2;
	this->gt_freq[0] = gf[0];
	this->gt_freq[1] = gf[1];
	
// then calculate hwe chi-square
	float hwe_gf[3];
	hwe_gf[0] = ref_af * ref_af;
	hwe_gf[1] = ref_af * ( 1 - ref_af );
	hwe_gf[2] = ( 1 - ref_af ) * ( 1 - ref_af );
	float num = 0;
	float denum = 0;
	for( int i = 0; i < n; ++i) {
		if ( DPs[i] <= 0 ) // skip missing
				continue;
		for( int j=0; j<3; j++ ) {
			num += prob_gls[i][j] * hwe_gf[j];
			denum += prob_gls[i][j] * gf[j];
		}
	}
	this->hwe = -2 * log(num / denum);
	
// then FIC
	num = 0;
	denum = 0;
	for( int sp = 0; sp < (int)Dosages.size(); sp++ ) {
		if ( DPs[sp] <= 0 ) // skip missing
			continue;
		float o_het_sum = gf[1] * prob_gls[sp][1];
		float o_sum = gf[0] * prob_gls[sp][0];
		o_sum += o_het_sum;
		o_sum += gf[2] * prob_gls[sp][2];
		
		float e_het_sum = hwe_gf[1] * prob_gls[sp][1];
		float e_sum = hwe_gf[0] * prob_gls[sp][0];
		e_sum += e_het_sum;
		e_sum += hwe_gf[2] * prob_gls[sp][2];
		
		num += o_het_sum / o_sum;
		denum += e_het_sum / e_sum;
	}
	if ( denum == 0 ) {
		cerr << "Warning: [ConsensusVcfRecord::EstimateFIC()] Unable to estimate FIC with missing PL info at position: " << position << endl;
		this->fic = 2;
	}
	else
		this->fic = 1 - num / denum;	
}


void  ConsensusVcfRecord::printInfoStr( ofstream & ovcf )
{
	int nv = countSampleWithVariant();
	if ( nv == 0 ) {
		cerr << "Warning: [ConsensusVcfRecord::generateInfoStr()] No variant at position: " << this->position << endl;
		ovcf << "END=NA;CIPOS=0,0;CIEND=0,0;AF=NA;FIC=NA";
	}
	else {
		int variant_end = round( (float)variant_end_diff_sum / nv ) + variant_end_baseline;
		int ci_low = left_most - position;
		int ci_high = right_most - position;
		int ci_end_low = left_most_end - variant_end;
		int ci_end_high = right_most_end - variant_end;
		ovcf << "END=" << variant_end << ";CIPOS=" << ci_low << "," << ci_high;
		ovcf << ";CIEND=" << ci_end_low <<  "," << ci_end_high;
		if ( lanchor == 0 )
			ovcf << ";SR=0";
		else if ( ranchor == 0 )
			ovcf << ";SR=INF";
		else {
			float sratio = (float)lanchor / ranchor;
			ovcf << ";SR=" << std::setprecision(2) << sratio;
		}
		ovcf << ";AF=" << std::setprecision(4) << 1-ref_af << ";HWE=" << hwe << ";FIC=";
	// check fic
		if ( fic >= -1 && fic <= 1 ) {
			ovcf << fic;
		}
		else
			ovcf << "NA";
	}
}


void ConsensusVcfRecord::printGLStr( int sp, ofstream & ovcf )
{
	ovcf << getGenotypeFromGLandAF( sp ) << ":" << DPs[sp] << ":" << GQs[sp] << ":";
	int PLs[3];
	setPLsFromGL( sp, &PLs[0] );
	ovcf << PLs[0] << "," << PLs[1] << "," << PLs[2];
}

string ConsensusVcfRecord::getGenotypeFromGLandAF( int sp )
{
	string gt;
	float sum = 0;
	for( int i=0; i<3; i++ )
		sum += GLs[sp][i];
	if ( sum == 0 ) // missing
		gt = "./.";
	else { // calculate posterior
		float ngl[3];
		sum = 0;
		for( int i=0; i<3; i++ )
			ngl[i] = exp(GLs[sp][i]);
		ngl[0] *= gt_freq[0];
		ngl[1] *= gt_freq[1];
		ngl[2] *= (1 - gt_freq[0] - gt_freq[1]);
		if ( ngl[2] > ngl[1] ) {
			if ( ngl[2] > ngl[0] )
				gt = "1/1";
			else // 0 >= 2
				gt = "0/0";
		}
		else { // 1 >= 2
			if ( ngl[1] > ngl[0] )
				gt = "0/1";
			else // 0 >= 1
				gt = "0/0";		
		}
	}
	return gt;
}

int ConsensusVcfRecord::getDosageFromGenotype( string & gt )
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


void ConsensusVcfRecord::estimateRawVariantQuality()
{
	int sum = 0;
	for( vector< int >::iterator ptr = GQs.begin(); ptr != GQs.end(); ptr++) {
		if ( (*ptr) >= 0 )
			sum += (*ptr);
	}
	this->variant_quality = sum;
}


int ConsensusVcfRecord::getVarintQualityWithAF()
{
	if ( variant_quality <= 0 )
		estimateRawVariantQuality();
	int ns = (int)DPs.size();
	int prvq;
	if ( gt_freq[0] == 0 ) // use 1/3n as minimum gt freq
		prvq = variant_quality - 10*ns*log10( 0.33 / ns );
	else
		prvq = variant_quality - 10*ns*log10( gt_freq[0] );
	return prvq;
}

void ConsensusVcfRecord::setPLsFromGL( int sp, int* vp)
{
// check if no info
	if ( GLs[sp][0] == 0 && GLs[sp][1] == 0 && GLs[sp][2] == 0) {
		vp[0] = 0;
		vp[1] = 0;
		vp[2] = 0;
		return;
	}
	
// convert ln to log
// GL is already normalized. no need to normalize again
	vp[0] = round( -GLs[sp][0] / 2.302585 * 10);
	vp[1] = round( -GLs[sp][1] / 2.302585 * 10);
	vp[2] = round( -GLs[sp][2] / 2.302585 * 10);	
}

void ConsensusVcfRecord::updateCIPOS( int pos )
{
// left
	if ( left_most == 0 )
		left_most = pos;
	else {
		if ( pos < left_most )
			left_most = pos;
	}
// right
	if ( right_most == 0 )
		right_most = pos;
	else {
		if ( pos > right_most )
			right_most = pos;
	}
}

void ConsensusVcfRecord::updateVariantEnd( int current_end )
{
	variant_end_diff_sum += (current_end - this->variant_end_baseline);
}


void ConsensusVcfRecord::updateCIEND( int current_end )
{
// left
	if ( left_most_end == 0 )
		left_most_end = current_end;
	else {
		if ( current_end < this->left_most_end )
			this->left_most_end = current_end;
	}
// right
	if ( right_most_end == 0 )
		right_most_end = current_end;
	else {
		if ( current_end > this->right_most_end )
			this->right_most_end = current_end;
	}
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



