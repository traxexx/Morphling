#include "OriginalStats.h"
#include "Utilities.h"
#include "GLs.h"
#include "Globals.h"
#include "bamUtility.h"
#include "VcfRecord.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;
using std::to_string;
using std::stoi;
using std::ifstream;
using std::stringstream;
using std::getline;

// constructor: only set mei index
OriginalStats::OriginalStats( int mei_type, string & sample_name ):
	ClipStart (2),
	DiscStart (10),
	RawCellSize (18),
	mei_index( mei_type ),
	SampleName( sample_name ),
	current_add_start( 0 )
{
	rawStats.clear();
	chrNameHash.clear();
	chrIndexVec.clear();
	if ( mei_index == 0 )
		mei_name = string("Alu");
	else if ( mei_index == 1 )
		mei_name = string("L1");
	else if ( mei_index == 2 )
		mei_name = string("SVA");
	else {
		cerr << "ERROR: mei_index = " << mei_index << ", should be >=0 & <=2 !" << endl;
		exit(1);
	}
}

// destructor
OriginalStats::~OriginalStats() {}


// add
bool OriginalStats::Add( string current_chr, string & proper_name, string & disc_name )
{
	ifstream properFile;
	string proper_zero = string(proper_name) + ".0";
	properFile.open( proper_zero.c_str() );
	if ( !properFile.is_open() ) {
		cerr << "Warning: Can't open " << proper_zero << " when adding to OriginalStats!" << endl;
		return 0; // add fail
	}
// add .0 first
	int current_chr_index = chrNameHash.size();
	if ( chrNameHash.find( current_chr ) != chrNameHash.end() ) {
		cerr << "ERROR: chr" << current_chr << " already exists!" << endl;
		exit(1);
	}
	chrNameHash[ current_chr ] = current_chr_index;
  // new cell for add
	RawCell * add_raw = new RawCell;
	add_raw->chr_index = current_chr_index;
	add_raw->win_index = 0;
	add_raw->counts.resize( RawCellSize, 0);
// for adding the next 2
	current_add_start = rawStats.size();
	string line;
	while( getline( properFile,line ) ) {
		stringstream ss;
		ss << line;
		string str_cord_index;
		getline( ss, str_cord_index, '\t');
		string str_proper;
		getline( ss, str_proper, '\t');
		string str_short;
		getline( ss, str_short, '\t');
		add_raw->win_index = stoi( str_cord_index );
		add_raw->counts[0] = stoi( str_proper);
		add_raw->counts[1] = stoi( str_short);
		rawStats.push_back( *add_raw );
	}
	delete add_raw;
	properFile.close();

// the proper-chr.mei_index
	string clip_name = string(proper_name) + "." + to_string(mei_index + 1);
	appendRawStats( clip_name, ClipStart );
// then disc
	string new_disc_name = disc_name + string(".") + to_string(mei_index);
	appendRawStats( new_disc_name, DiscStart );
	return 1;
}	

// append, used in add clip & disc
void OriginalStats::appendRawStats( string & rec_name, int base )
{
	ifstream recFile;
	recFile.open( rec_name.c_str() );
	CheckInputFileStatus( recFile, rec_name.c_str() );
	vector< RawCell >::iterator rs = rawStats.begin();
	if ( rawStats.size() == 0 ) {
		cerr << "ERROR: no 0x2 read pairs added to rawStats. Does the original bam consist of single-end reads?" << endl;
		exit(1);
	}
	rs += current_add_start;
	bool reach_end = 0;
	string line;
	while( getline( recFile, line ) ) {
		stringstream ss;
		ss << line;
		string str_cord_index;
		getline( ss, str_cord_index, '\t');
		int cord_index = stoi( str_cord_index ); // index of current window to add
		
		if ( cord_index > rs->win_index ) { // find next
			while( cord_index > rs->win_index ) {
				rs++;
				if (rs == rawStats.end()) {
					reach_end = 1;
					break;
				}
			}
			if (reach_end) // no more to add (no available proper any more )
				break;
		}
		// now cord_index <= rs->win_index
		
		if ( rs->win_index == cord_index) { // add
			for( int i = 0; i < 8; i++ ) {
				string field;
				getline( ss, field, '\t' );
				rs->counts[ i + base ] = stoi( field );
			}
			rs++;
			if ( rs == rawStats.end() )
				break; // all added
		}
		else { // cord < ref, no clip in current cord
			continue;
		}
	}
	recFile.close();
}


// construct genome link
void OriginalStats::ReOrganize()
{
// make chrIndexHash
	buildChrIndexVec();
	
// sort rawStats
	std::sort( rawStats.begin(), rawStats.end(), sortRawStats );

	vector<int> dupVec; // record dup
	setDupVecForRawStats( dupVec );
	int empty = 0; // first several empty cells
	vector<RawCell>::iterator rs = rawStats.begin();
	int sum = 0;
	for( vector<int>::iterator it = rs->counts.begin(); it != rs->counts.end(); it++ )
		sum += *it;
	if (sum == 0)
		empty += dupVec[0];
// pre-allocate merge data
	if (empty == 0)
		MergeData.resize( dupVec.size() );
	else
		MergeData.resize( dupVec.size() - 1 );
		
// raw->merge & set genome location map
	vector<RawCell>::iterator raw_ptr = rawStats.begin();
	vector<int>::iterator dup_ptr = dupVec.begin();
	MergeCellPtr merge_ptr = MergeData.begin();
	if (empty != 0) {
		raw_ptr += dupVec[0];
		dup_ptr++;
	}
  // start: clear rawStats when finish copy (to save memory)
  	GenomeLocationCell new_cell;
  	for( ; dup_ptr != dupVec.end(); dup_ptr++, merge_ptr++ ) {
  		merge_ptr->counts = std::move( raw_ptr->counts );
  		merge_ptr->dups = (*dup_ptr);
  		for( int ct = 0; ct < (*dup_ptr); ct++ ) {
  			new_cell.win_index = raw_ptr->win_index;
  			new_cell.ptr = merge_ptr;
  			GenomeLocationMap[ convertChrIndexToName( raw_ptr->chr_index ) ].push_back( new_cell );
  			raw_ptr++;
  		}
  	}

// sort genome location map per chr
	for( map< string, vector< GenomeLocationCell > >::iterator map_it = GenomeLocationMap.begin(); map_it != GenomeLocationMap.end(); map_it++ ) {
		std::sort( map_it->second.begin(), map_it->second.end(), sortGenomeLocations );
	}
	
// clear
	rawStats.clear();
}

// clear the cells with level < LEVEL
void OriginalStats::ClearUnderLevelMergeCells()
{
	int erase_cell = 0;
	for( MergeCellPtr mptr = MergeData.begin(); mptr != MergeData.end(); mptr++ ) {
		int n_support = getSumSupportDiscs( mptr->counts ) + getSumSupportUnmaps(mptr->counts) + getSumSupportClips(mptr->counts);
		if ( n_support < LEVEL ) {
			mptr->counts.clear();
			erase_cell++;
		}
	}
	if (DEBUG_MODE)
		cout << "Removed under level cell = " << erase_cell << endl;
}

// print out: single chr
void OriginalStats::PrintGLasVcf( string & vcf_name, string & bam_name, string & ref_fasta, string & focus_chr )
{
// open bam & load fasta
	if ( !currentSam.IsOpen() )
		OpenBamAndBai( currentSam, currentSamHeader, bam_name );
	string alt_str = string("<INS:ME:") + mei_name + ">";
// do print. Also merge consecutive windows. Print out the one with highest GL
	ofstream outVcf;
	outVcf.open( vcf_name.c_str() );
	CheckOutFileStatus( outVcf, vcf_name.c_str() );
	string dummy = string("");
	if ( GenomeLocationMap.find( focus_chr ) == GenomeLocationMap.end() ) {
		cerr << "Warning: no data for chr: " << focus_chr << ", output empty split vcf..." << endl;
		return;
	}
	GlcPtr anchor; // mark the anchor
	VcfRecord * current_vcf_rec = NULL;
	for( GlcPtr cell_it = GenomeLocationMap[ focus_chr ].begin(); cell_it != GenomeLocationMap[ focus_chr ].end(); cell_it++ ) {
		if (PRINT_NON_VARIANT) { // for debug purpose only
			outVcf << focus_chr << "\t" << cell_it->win_index << "\t.\t.\t<INS:ME:" << mei_name << ">\t";
			MergeCellPtr Merge = cell_it->ptr;
			int clip_count = getSumSupportClips( Merge->counts );
			int disc_count = getSumSupportDiscs( Merge->counts);
			int unmap_count = getSumSupportUnmaps( Merge->counts );
			bool left_present = Merge->counts[4] + Merge->counts[13] + Merge->counts[17];
			bool right_present = Merge->counts[9] + Merge->counts[11] + Merge->counts[15];
			bool both_end = left_present & right_present;
			outVcf << "-1\tPASS\tBOTH_END=" << both_end << "\t" << ";CLIP=" << clip_count << ";DISC=" << disc_count << ";UNMAP=" << unmap_count;
			outVcf << ";";
			for( vector<int>::iterator p = Merge->counts.begin(); p != Merge->counts.end(); p++ ) {
				outVcf << *p << ",";
			}
			int dp = GetVecSum( Merge->counts );
			outVcf << "\tGT:DP:GQ:PL\t";
			if( cell_it->ptr->GL.size() == 3 ) {
				string genotype = GetGenotype( Merge->GL );
				outVcf << genotype << ":" << dp << ":-1:";
				outVcf << Merge->GL[0] << "," << Merge->GL[1] << "," << Merge->GL[2] << endl;
			}
			else {
				outVcf << "0/0:" << dp << ":-1:NA,NA,NA" << endl;
			}			
		}
		if (cell_it->ptr->GL.size() != 3)
			continue;
		if ( cell_it->ptr->GL[0] >= cell_it->ptr->GL[1] && cell_it->ptr->GL[0] >= cell_it->ptr->GL[2] ) { // skip the window with no variant
			continue;
		}
		if ( current_vcf_rec == NULL ) { // new anchor
			anchor = cell_it;
			current_vcf_rec = new VcfRecord;
			current_vcf_rec->SetFromGenomeLocationCell( *anchor );
		}
		else { // compare with previous anchor
			int new_start = cell_it->win_index * STEP;
			if ( new_start - current_vcf_rec->GetVariantEnd() > 0 ) { // not consecutive. start new anchor. print old anchor
				current_vcf_rec->SetChrName( focus_chr );
				current_vcf_rec->SetAltAllele( alt_str );
				if ( REFINE_BREAK_POINT )
					current_vcf_rec->SetBreakPointAndCIFromBam( currentSam, currentSamHeader );
				current_vcf_rec->PrintRecord( outVcf);
				delete current_vcf_rec;
				anchor = cell_it;
				current_vcf_rec = new VcfRecord;
				current_vcf_rec->SetFromGenomeLocationCell( *anchor );
			}
			else { // consecutive. compare with anchor
				bool use_new = current_vcf_rec->UpdateByRankAnchor( anchor, cell_it );
				if ( use_new )
					anchor = cell_it;
			}
		}
	}
// print out last rec
	if ( current_vcf_rec != NULL ) {
		current_vcf_rec->SetChrName( focus_chr );
		current_vcf_rec->SetAltAllele( alt_str );
		if ( REFINE_BREAK_POINT )
			current_vcf_rec->SetBreakPointAndCIFromBam( currentSam, currentSamHeader );
		current_vcf_rec->PrintRecord( outVcf );
	}
	outVcf.close();	
}


/**** inner functions ****/

// build chrIndexVec from chrNameHash for later use
void OriginalStats::buildChrIndexVec()
{
	int vec_size = chrNameHash.size();
	chrIndexVec.resize( vec_size, string("") );
	for( map< string, int >::iterator map_it = chrNameHash.begin(); map_it != chrNameHash.end(); map_it++ ) {
		if ( map_it->second >= vec_size ) {
			cerr << "ERROR: chr index " << map_it->second << " larger than chrNameHash size!" << endl;
			exit(1);
		}
		chrIndexVec[ map_it->second ] = map_it->first;
	}
// sanity check: if any chrIndexVec element is uninitialized
	for( vector< string >::iterator it = chrIndexVec.begin(); it != chrIndexVec.end(); it++ ) {
		if ( it->size() == 0 ) {
			cerr << "ERROR: index " << ( it - chrIndexVec.begin() ) << " in chrIndexVec is not initialized!" << endl;
			exit(1);
		}
	}
}

// set dup vec by counting dups in sorted rawStats (inner)
void OriginalStats::setDupVecForRawStats( vector<int> & dupVec )
{
	if ( rawStats.empty() ) {
		cerr << "ERROR: no available read info loaded in rawStats. Check if proper* and disc* in ctrl_tmps/ are all empty!" << endl;
		exit(1);
	}
	dupVec.clear();
	vector<RawCell>::iterator raw = rawStats.begin();
	vector<RawCell>::iterator prev = raw;
	raw++;
	int dup = 1;
	for( ; raw != rawStats.end(); raw++, prev++ ) {
		if ( prev->counts == raw->counts ) // it is ok since int == method is defined
			dup++;
		else {
			dupVec.push_back( dup );
			dup = 1;
		}
	}
	dupVec.push_back( dup );
}



// sort rawStats: counts > chr > win
bool OriginalStats::sortRawStats( RawCell x, RawCell y )
{
	if ( x.counts.size() != y.counts.size() ) {
		cerr << "ERROR: sortRawStats x y do not have same counts!" << endl;
		exit(1);
	}
	vector<int>::iterator x_it = x.counts.begin();
	vector<int>::iterator y_it = y.counts.begin();
	for( ; x_it != x.counts.end(); x_it++, y_it++ ) {
		if ( *x_it < *y_it )
			return 1;
		else if ( *x_it > *y_it )
			return 0;
	}
// all the same then coord
	if ( x.chr_index < y.chr_index )
		return 1;
	else if ( x.chr_index > y.chr_index )
		return 0;
	else { // compare win index
		if ( x.win_index < y.win_index )
			return 1;
		else if ( x.win_index > y.win_index )
			return 0;
		else {
			cerr << "ERROR: sortRawStats: Same location record appear twice at: chr" << x.chr_index << ": " << x.win_index << endl;
			exit(1);
		}
	}
}


// sort GenomeLocation in each chr
bool OriginalStats::sortGenomeLocations( GenomeLocationCell x, GenomeLocationCell y )
{
	if ( x.win_index < y.win_index )
		return 1;
	else if ( x.win_index > y.win_index )
		return 0;
	else {
		cerr << "ERROR: sortGenomeLocation: Same location record appear twice at: " << x.win_index << endl;
		exit(1);
	}
}


/*** other utility functions ****/

// convert chr name to index
int OriginalStats::convertChrNameToIndex( string chr_name )
{
	if ( chrNameHash.find( chr_name ) == chrNameHash.end() ) {
		cerr << "ERROR: Can't find " << chr_name << " in chrNameHash!" << endl;
		exit(1);
	}
	int chr_index = chrNameHash[ chr_name ];
	return chr_index;
}

// convert chr index back to name
string OriginalStats::convertChrIndexToName( int chr_index )
{
	if ( chr_index >= int( chrIndexVec.size() ) ) {
		cerr << "ERROR: Can't find index " << chr_index << " in chrIndexVec!" << endl;
		exit(1);
	}
	string chr_name = chrIndexVec[ chr_index ];
	return chr_name;
}















