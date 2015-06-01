#include <sstream>
#include <algorithm> // std::all_of
#include "VcfRecord.h"
#include "Globals.h"
#include "SamFile.h"
#include "GLs.h"
#include "bamUtility.h"
#include "QC.h"
#include "Utilities.h"
#include "ProperDeckUtility.h"
#include "ReadMapGlobals.h"

using std::stringstream;
using std::getline;
using std::cerr;
using std::endl;
using std::to_string;

// constructor
VcfRecord::VcfRecord():
	chr( string("") ),
	position( -1 ),
	ref_allele( string("") ),
	alt_allele( string("") ),
	variant_quality( -1 ),
	filter( string("") ),
	info_field( string("") ),
	format_field( string("") ),
	gl_field( string("") ),
	info_parsed( 0 ),
	gl_parsed( 0 ),
	depth( -1 ),
	dosage( -1 ),
	breakp_refined( 0 ),
	both_end( 0 ),
	variant_end( -1 ),
	ci_low(0),
	ci_high(0),
	clip_count(-1),
	disc_count(-1),
	unmap_count(-1),
	evidence(-1),
	win_count( -1 )
{
	read_counts.clear();
	GL.clear();
}

// destructor
VcfRecord::~VcfRecord() {}

/*****set methods******/

// from file
void VcfRecord::SetFromLine( string & line )
{
	if ( line.size() < 1 ) {// create an empty class
		cerr << "Warning: use an empty line to set vcf record!" << endl;
		return;
	}

	stringstream ss;
	ss << line;
	getline(ss, chr, '\t');
	string position_str;
	getline(ss, position_str, '\t');
	position = stoi( position_str );
	getline(ss, id, '\t'); // ref
	getline(ss, ref_allele, '\t'); // ref
	getline(ss, alt_allele, '\t'); // alt
	string qual_str;
	getline(ss, qual_str, '\t'); // qual
	variant_quality = stoi( qual_str );
	getline(ss, filter, '\t'); // filter
	getline(ss, info_field, '\t'); // info
	getline(ss, format_field, '\t'); // format
	getline(ss, gl_field, '\t'); // gl field
}

/* from input parameters
void VcfRecord::SetFromInputParameters(string & chr_name, int breakp, string & ref_al, string & alt_al, int quality, string & filter_str, string & info_str, string & format_str, string & gl_str)
{
	chr = chr_name;
	position = breakp;
	ref_allele = ref_al;
	alt_allele = alt_al;
	variant_quality = quality;
	filter = filter_str;
	info_field = info_str;
	format_field = format_str;
	gl_field = gl_str;
}
*/

// from merge cell ptr
void VcfRecord::SetFromGenomeLocationCell( GenomeLocationCell & glc )
{
	position = glc.win_index * STEP + WIN / 2;
	win_count = 1;
	UpdateFromMergeCellPtr( glc.ptr );
}

/*** set methods members ***/
void VcfRecord::SetChrName( string & chr_name) {
	if ( !chr.empty() )
		cerr << "Warning: replace existing chr: " << chr << " to " << chr_name << "?" << endl;
	chr = chr_name;
}

void VcfRecord::SetPosition( int center )
{
	position = center;
}

void VcfRecord::SetAltAllele( string & alt_str )
{
	alt_allele = alt_str;
}

// add: info_str=info_val;
void VcfRecord::AddIntegerInfoField( const char * info_str, int info_val )
{
// update info first
	if ( info_field.empty() )
		updateInfoField();
// then add
	info_field += ( ";" + string(info_str) + "=" + std::to_string(info_val) );
}

/** get methods for members **/
string VcfRecord::GetChromosome()
{
	if ( chr.empty() ) {
		cerr << "ERROR: chr name not defined at: " << position << endl;
		exit(1);
	}
	return chr;
}	
	
int VcfRecord::GetPosition()
{
	if ( position < 0 ) {
		cerr << "ERROR: position not defined at: " << position << endl;
		exit(1);
	}
	return position;
}

int VcfRecord::GetVariantQuality()
{
	return variant_quality;
}

bool VcfRecord::GetFilterPassOrNot()
{
	if ( filter.compare("PASS") == 0 )
		return 1;
	else
		return 0;
}

int VcfRecord::GetVariantEnd()
{	
	if ( variant_end < 0 ) {
		cerr << "ERROR: variant end is not set at: " << position << endl;
		exit(1);
	}
	return variant_end;
}

int VcfRecord::GetVariantDepth()
{
	if ( !info_parsed ) {
		parseGLfield();
	}
	if ( depth < 0 ) {
		cerr << "ERROR: dp not defined at: " << chr << ":" << position << endl;
		exit(1);
	}
	return depth;
}

int VcfRecord::GetDosage()
{
	if ( !gl_parsed ) {
		parseGLfield();
	}
	if ( dosage < 0 ) {
		cerr << "ERROR: gt not defined at: " << chr << ": " << position << endl;
		exit(1);
	}
	return dosage;
}


int VcfRecord::GetEvidenceDepth()
{
	if ( evidence < 0 )
		parseInfoField();
	return evidence;
}


/** update methods from merge cell ptr ***/
void VcfRecord::UpdateFromMergeCellPtr( MergeCellPtr & ptr )
{
	GL = ptr->GL;
	variant_quality = GetVariantQualityFromGL( GL );
	variant_end = position + WIN / 2;
	depth = GetVecSum( ptr->counts );
	read_counts = ptr->counts;
// take as parsed
	info_parsed = 1;
	gl_parsed = 1;
// set necessary fields
	if ( dosage < 0 )
		dosage = GetAlleleDosage( ptr->GL );
//	updateGLfield();
	if ( dosage > 0 ) { // MEI hit
		updateEvidenceFromReadCount();
		updateBothEnd();
		if ( win_count < 1 )
				win_count = 1;
		updateFilter();
//		updateInfoField();
	}
	else { // reference
		info_field = string(".");
		filter = string(".");
	}
}

/*** update methods for anchor window in discovery phase ***/
bool VcfRecord::UpdateByRankAnchor( GlcPtr & anchor, GlcPtr & new_anchor )
{
// rule: dosage --> posterior-variant --> %(disc + clip + unmap) --> less %proper -->depth --> use anchor
//bool CompareDosage( vector<int> & anchorGL, vector<int> & newGL );

// update necessary fields first
	win_count++;

	MergeCellPtr anchor_ptr = anchor->ptr;
	MergeCellPtr new_ptr = new_anchor->ptr;
	int use_new = -1; // will change to 0 or 1 during comparison. if equal, then keep -1
// update gq first
	int anchor_gq_sig_higher = 0; // for comparing posterior next
	int new_gq = GetVariantQualityFromGL( new_ptr->GL );
	if ( new_gq > variant_quality ) {
		if ( new_gq - variant_quality >= 30 )
			anchor_gq_sig_higher = -1;
	}
	else if ( new_gq < variant_quality ) {
		if ( variant_quality - new_gq >= 30 )
			anchor_gq_sig_higher = 1;
	}

// compare posterior variant
	if ( anchor_gq_sig_higher > 0 )
		use_new = 0;
	else if ( anchor_gq_sig_higher < 0 )
		use_new = 1;
	else {
		// calculate depth here but do not use as comparison for now
		int new_dp = GetVecSum( new_ptr->counts );
		
	// %support
		float anchor_support_frac = float(evidence) / depth;
		float new_support_frac = GetSupportReadFraction( new_ptr->counts, new_dp );
		if ( anchor_support_frac > new_support_frac )
			use_new = 0;
		else if ( anchor_support_frac < new_support_frac )
			use_new = 1;
		else {
		// %proper
			int anchor_proper = GetModifiedProperReadFraction( anchor_ptr->counts, depth );
			int new_proper = GetModifiedProperReadFraction( new_ptr->counts, new_dp );
			if ( anchor_proper < new_proper )
				use_new = 0;
			else if ( anchor_proper > new_proper )
				use_new = 1;
			else {	
			// depth
				if ( depth > new_dp )
					use_new = 0;
				else if ( new_dp > depth )
					use_new = 1;
			}
		}
	}
	
// set return value & update
	bool bool_new;
	if ( use_new < 0 ) {
		updateRecWithEqualAnchorIndex( new_anchor->win_index );
		bool_new = 0;
	}
	else if ( use_new > 0 ) {
		updateRecWithNewAnchorIndex( new_anchor->win_index );
		UpdateFromMergeCellPtr( new_ptr );
		bool_new = 1;
	}
	else // use old anchor
		bool_new = 0;
			
	return bool_new;
}

void VcfRecord::updateRecWithNewAnchorIndex( int win_index )
{
	position = win_index * STEP + WIN / 2;
}

void VcfRecord::updateRecWithEqualAnchorIndex( int win_index )
{
	int new_position = win_index * STEP + WIN / 2;
	int dist = (new_position - position) / 2;
	ci_low -=  dist;
	ci_high += dist;
	position = ( position + new_position ) / 2;
	variant_end = new_position + WIN / 2;
	updateBothEnd();
}

void VcfRecord::updateBothEnd()
{
	if ( read_counts.size() == 0 ) {
		cerr << "ERROR: cannot update both end info since no read_count vector exists in VcfRecord: " << chr << ": " << position <<"!" << endl;
		exit(1);
	}
// if already both end, return
	if ( both_end )
		return;
	bool left_present = read_counts[4] + read_counts[13] + read_counts[17];
	bool right_present = read_counts[9] + read_counts[11] + read_counts[15];
	both_end = (left_present & right_present);
}

void VcfRecord::updateFilter()
{
	filter.clear();
	if ( !both_end )
		filter = "BOTH_END";
	bool depth_pass = depthQC();
	if( !depth_pass ) {
		if ( !filter.empty() )
			filter += '+';
		filter += "DEPTH";
	}
	if ( evidence < LEVEL ) {
		if ( !filter.empty() )
			filter += '+';
		filter += "+SUP_READS";
	}
	if ( filter.empty() )
		filter = "PASS";
}

void VcfRecord::updateInfoField()
{
	info_field = string( "SVTYPE=INS;END=" ) + to_string(variant_end) + ";CIPOS=" + to_string(ci_low) + "," + to_string(ci_high);
	info_field += ";CLIP=" + to_string(clip_count) + ";DISC=" + to_string(disc_count) + ";UNMAP=" + to_string(unmap_count) + ";WCOUNT=" + to_string(win_count);
}

void VcfRecord::updateEvidenceFromReadCount()
{
	clip_count = getSumSupportClips( read_counts );
	disc_count = getSumSupportDiscs( read_counts );
	unmap_count = getSumSupportUnmaps( read_counts );
	updateEvidence();
}

void VcfRecord::updateEvidence()
{
	if ( clip_count < 0 || unmap_count < 0 || disc_count < 0 ) {
		cerr << "ERROR: unable to get evidence count at info field: " << info_field << endl;
		exit(1);
	}
	evidence = clip_count + unmap_count + disc_count;
}

void VcfRecord::updateGLfield()
{
	format_field = string("GT:DP:GQ:PL");
	vector<int> PL;
	PL.clear();
	SetPLsFromGL( PL, GL );
	int gt_quality = dosage > 0 ? GetGenotypeQuality( GL ) : 0;
	string genotype = GetGenotype( GL );
	gl_field = genotype + ":" + to_string(depth) + ":" + to_string(gt_quality) + ":";
	gl_field += to_string(PL[0]) + "," + to_string(PL[1]) + "," + to_string(PL[2]);
	
}

/*** filters ***/
bool VcfRecord::depthQC()
{
	if ( depth < MIN_READ_IN_WIN || depth > MAX_READ_IN_WIN )
		return 0;
	else
		return 1;
}

/*** inner parsers ***/
void VcfRecord::parseGLfield()
{
	if ( format_field.empty() || gl_field.empty() ) {
		cerr << "ERROR: no format or gl field recorded: " << chr << ": " << position << endl;
		exit(1);
	}
// parse format field
	stringstream ss;
	ss << format_field;
	int gt_index = -1;
	int dp_index = -1;
	int fc = 0;
	string field;
	while( getline(ss, field, ':') ) {
		if ( field.compare("GT") == 0 )
			gt_index = fc;
		else if ( field.compare("DP") == 0 )
			dp_index = fc;
		fc++;
	}
	if ( fc == 0 ) {
		cerr << "ERROR: can't parse gl field: " << gl_field << ", maybe not separated by ':' ?" << endl;
		exit(1);
	}
/* parse gl field
	vector<string> gl_items;
	gl_items.resize( fc );
	vector<string>::iterator it = gl_items.begin();
*/
	stringstream glss;
	glss << gl_field;
	fc = 0;
	while( getline( glss, field, ':' ) ) {
		if ( fc == gt_index ) {
			dosage = GetDosageFromGenotype( field );
		}
		else if ( fc == dp_index ) {
			depth = stoi( field );
		}
		fc++;
	}
// clear
	gl_parsed = 1;
}


void VcfRecord::parseInfoField()
{
	if ( info_field.empty() ) {
		cerr << "ERROR: no info field to parse at: " << chr << ": " << position << endl;
		exit(1);
	}
// parse
	stringstream ss;
	ss << info_field;
	string field;
	while( getline(ss, field, ';') ) {
		stringstream sub_ss;
		sub_ss << field;
		string sub_item;
		getline( sub_ss, sub_item, '=' );
		if ( sub_item.compare("END") == 0 ) {
			string ve_str;
			getline( sub_ss, ve_str, '=' );
			variant_end = stoi( ve_str );
		}
		else if ( sub_item.compare("CIPOS") == 0 ) {
			string int_str;
			getline( sub_ss, int_str, '=' );
			clip_count = stoi( int_str );
		}
		else if ( sub_item.compare("CLIP") == 0 ) {
			string int_str;
			getline( sub_ss, int_str, '=' );
			clip_count = stoi( int_str );		
		}
		else if ( sub_item.compare("DISC") == 0 ) {
			string int_str;
			getline( sub_ss, int_str, '=' );
			disc_count = stoi( int_str );		
		}
		else if ( sub_item.compare("UNMAP") == 0 ) {
			string int_str;
			getline( sub_ss, int_str, '=' );
			unmap_count = stoi( int_str );		
		}
		else if ( sub_item.compare("WCOUNT") == 0 ) {
			string int_str;
			getline( sub_ss, int_str, '=' );
			win_count = stoi( int_str );		
		}
	}
	updateEvidence();
// clear
	info_parsed = 1;
}


/*** refine break point ***/
void VcfRecord::SetBreakPointAndCIFromBam( SamFile & currentSam, SamFileHeader & currentSamHeader )
{
	if ( breakp_refined )
		return;

	if ( depth <= 0 ) {
		if (DEBUG_MODE)
			cerr << "Warning: No reads in section: " << chr << ": " << position << endl;
		return;
	}

	bool section_status = currentSam.SetReadSection( chr.c_str(), position - WIN / 2 - 200, position + WIN / 2 + 200 );
	if ( !section_status ) {
		cerr << "ERROR: Unable to set bam section: " << chr << ": " << position - WIN / 2 - 200 << " - " << position + WIN / 2 + 200 << endl;
		exit(1);
	}
	int start = position - WIN / 2;
	
// build vector of spanning counts
	vector<int> locs;
	vector<int> clips;
	clips.resize( WIN, 0 );
	locs.resize( WIN, 0 );	
	SamRecord sam_rec;
	while( currentSam.ReadRecord( currentSamHeader, sam_rec ) ) {
		bool qc = PassQC( sam_rec );
		if ( !qc )
			continue;
		int st = sam_rec.get1BasedPosition() - start;
		int ed = sam_rec.get1BasedAlignmentEnd() - start;
		ed -= 10;
		if ( ed < 0 )  // in locs?
			continue;
		st += 10;
		if ( st >= WIN )  // in locs?
			continue;
		if ( ed - st < 0 ) // enough length?
			continue;
		int true_st = st >= 0 ? st : 0;
		int true_ed = ed < WIN ? ed : WIN - 1;
		for( int i = true_st; i <= true_ed; i++ ) {
			locs[i]++;
		}
// add clip info
		int cliplen = getMaxClipLen( sam_rec );
		if ( cliplen == 0 )
			continue;
		int index;
		if ( abs(cliplen) < 10 ) //skip short clip
			continue;
		if ( cliplen > 0 ) { // begin clip
			index = st - 10;
		}
		else { // end clip
			index = ed + 10;
		}
		if ( index < 0 || index >= WIN ) // skip those out of WIN region
			continue;
	// add to map
		clips[ index ]++; 
	}

// get breakp
// utilize clip first
	bool no_clip = std::all_of( clips.begin(), clips.end(), [](int i) { return i == 0;} );
	if ( !no_clip ) { // exist clip
		vector<int> cluster = clips;
		for( int i = 24; i < WIN - 25; i++ ) {
			for( int j = 1; j <= 24; j++ )
				cluster[i] += clips[ i - j ];
			for( int j = 1; j <= 24; j++ )
				cluster[i] += clips[ i + j ];
		}
		int max_cluster = *std::max_element( cluster.begin(), cluster.end() );
		int avr_index = GetAvrLocationOfCertainValue( cluster, max_cluster, ci_low, ci_high );
		position = avr_index + start;
		ci_low = avr_index - ci_low > 25 ? ci_low - avr_index : -25;
		ci_high = ci_high - avr_index > 25 ? ci_high - avr_index : 25;
	// set event end at the next point where #span increase
		variant_end = start + WIN;
		for( int inc = avr_index + 1; inc < WIN; inc++ ) {
			if( locs[ inc ] > locs[inc - 1] ) {
				variant_end = inc + start;
				break;
			}
		}
	}
	else {
// if no clip, use span info: average each 3 nearby locs --> get minimum --> if equal values exist, take average
		vector<int> avrs;
		avrs.resize( WIN - 2 );
		for( int i = 1; i <= WIN - 2; i++ ) {
			avrs[i-1] = round( float(locs[i-1] + locs[i] + locs[i+1]) / 3 );
		}
		int min_span = *std::min_element( avrs.begin(), avrs.end() );
		int avr_index = GetAvrLocationOfCertainValue( avrs, min_span,  ci_low, ci_high );
		variant_end = start + WIN;
		ci_low = avr_index - ci_low > 25 ? ci_low - avr_index : -25;
		ci_high = ci_high - avr_index > 25 ? ci_high - avr_index : 25;
 		position = avr_index + start + 1;
	}	
	breakp_refined = 1;
}

void VcfRecord::PrintRecord( ofstream & out_vcf )
{
// set necessary fields
	if ( filter.empty() )
		updateFilter();
	if ( info_field.empty() )
		updateInfoField();
	if ( gl_field.empty() )
		updateGLfield();

// update all empty fields to dot	
	out_vcf << chr << "\t" << position << "\t.\t";
  // ref allele
	if ( ref_allele.empty() )
		out_vcf << ".\t";
	else
		out_vcf << ref_allele << "\t";
  // alt allele
	if ( alt_allele.empty() )
		out_vcf << ".\t";
	else
		out_vcf << alt_allele << "\t";
  // int
	out_vcf << variant_quality << "\t";
  // filter
	if ( filter.empty() )
		out_vcf << ".\t";
	else
		out_vcf << filter << "\t";
  // info
	if ( info_field.empty() )
		out_vcf << ".\t";
	else
		out_vcf << info_field << "\t";
  // these fields shouldn't be empty
	out_vcf << format_field << "\t" << gl_field << "\t" << endl;
}

// if start with '#', then is header line
bool IsInfoLine( string & line )
{
	if ( line[0] == '#' )
		return 1;
	else
		return 0;
}
