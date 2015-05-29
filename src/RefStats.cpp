#include "RefStats.h"
#include "GLs.h"
#include "Globals.h"

#include <iostream>
#include <math.h> // log round
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "Utilities.h"

using std::cout;
using std::endl;
using std::cerr;

/* A big constructor:
	load stats from files
	sort statst + merge dup
	set genome link
	mei_index 0~2 A L S
*/
RefStats::RefStats( string & proper_prefix, string & disc_prefix, int mei_type, AllHetIndex & allHetIndex ):
	ref_lh_set( 0 ),
	CountSize(18),
	min_log( log(0.000001) ),
	mei_index ( mei_type ),
	CountOffSet(0.3),
	REF_CHR( string("20") )
{
// reserve refstats
	refStats.resize(3);

// set null stat cell ptr
	Dummy.clear();
	NullStatCellPtr = Dummy.end();
// start doing real things	
	string test_name = string( "CtrlRemap" );
	osPtr = new OriginalStats( mei_type, test_name );
	osPtr->Add( REF_CHR, proper_prefix, disc_prefix);
	osPtr->ReOrganize();
	setStatRef( allHetIndex );
}

RefStats::~RefStats() {}

// for each MergeCell, based on counts, set GL
void RefStats::SetRecordGL( MergeCellPtr & merge )
{
// skip the cleared cells
	if ( merge->counts.size() == 0 )
		return;
	if ( merge->counts.size() != CountSize ) {
		cerr << "ERROR: wrong count size at SetPerRecordGL!" << endl;
		exit(1);
	}
	if (merge->GL.size() != 3)
		merge->GL.resize(3);
	for(int s_ref = 0; s_ref <= 2; s_ref++) {
		setSingleRefGL( merge, s_ref );
	}
}

// print ctrl power & novel info: detailed is an old parameter, no use now
void RefStats::PrintCtrlGLasRecord( string & outRecord, string & ctrl_bam, string & ctrl_fasta )
{	
	osPtr->PrintGLasVcf( outRecord, ctrl_bam, ctrl_fasta, REF_CHR );
// clear ref data
	delete osPtr;
}

/*** set ctrl GL ****/
void RefStats::SetCtrlGLs()
{
// drop the stats with #disc < LEVEL
	osPtr->ClearUnderLevelMergeCells();

// set GLs
	for( MergeCellPtr merge_it = osPtr->MergeData.begin(); merge_it != osPtr->MergeData.end(); merge_it++ ) {
		SetRecordGL( merge_it );
	}
	
	AdjustUsedLoci( osPtr );
}


/**** remove self-contribution *****/

// main function: set use_lift = 1 on real data
// when exclude_nearby, exclude window +- 1 around hit
// always initialize special merge when there is a hit (no matter hom/neg)
void RefStats::AdjustUsedLoci(  OriginalStats * dataOsPtr )
{
// initialize special merge data
	int destroy_range = WIN / STEP * 2;
	dataOsPtr->SpecialMergeData.resize( selfHom.size() + (selfHet.size() + selfNeg.size()) * (destroy_range * 2 + 1) );
	MergeCellPtr special_ptr = dataOsPtr->SpecialMergeData.begin();
// add to special
	for( GlcPtr genome = dataOsPtr->GenomeLocationMap[ REF_CHR ].begin(); genome != dataOsPtr->GenomeLocationMap[ REF_CHR ].end(); genome++ ) {
  // skip cleared cell first
  		if ( genome->ptr->counts.size() == 0 )
  			continue;
	
		map< int, SelfHomElement >::iterator hom_map = selfHom.find( genome->win_index );
		MergeCellPtr data_ptr = genome->ptr;
		if ( hom_map != selfHom.end() ) { // in hom
//cout << genome->win_index << "\tHOM" << endl;
			special_ptr->GL = data_ptr->GL;
			special_ptr->dups = 1;
			special_ptr->counts = data_ptr->counts;
		// adjust hom
			StatCellPtr hom_ptr = refStats[2].begin();
			if ( hom_map->second.hom_index >= refStats[2].size() ) {
				cerr << "ERROR: selfHom index: " << hom_map->second.hom_index << " > refStats[2] size: " << refStats[2].size() << endl;
				exit(1);
			}
			hom_ptr += hom_map->second.hom_index;
			adjustSingleGL( data_ptr->counts, hom_ptr, special_ptr, 2);
		// adjust het
			adjustMultipleGL( data_ptr->counts, hom_map->second.hets, special_ptr);
			genome->ptr = special_ptr;
		// move special ptr
			special_ptr++;
		}
		else { // maybe it's neg or het?
			vector<int> search_range;
			for(int i = -destroy_range; i<= destroy_range; i++)
				search_range.push_back(i);
			
			for( vector<int>::iterator range_ptr = search_range.begin(); range_ptr != search_range.end(); range_ptr++ ) {
				SelfSingleMap::iterator neg_map = selfNeg.find( genome->win_index + (*range_ptr) );
				SelfSingleMap::iterator het_map = selfHet.find( genome->win_index + (*range_ptr) );
				if ( neg_map != selfNeg.end() ) { //  in neg
//cout << genome->win_index << "\tNEG" << endl;
					special_ptr->GL = data_ptr->GL;
					special_ptr->dups = 1;
					special_ptr->counts = data_ptr->counts;
					StatCellPtr neg_ptr = refStats[0].begin();
					if ( neg_map->second.stat_index >= refStats[0].size() ) {
						cerr << "ERROR: selfNeg index: " << neg_map->second.stat_index << " > refStats[0] size: " << refStats[0].size() << endl;
						exit(1);
					}
					neg_ptr += neg_map->second.stat_index;
					adjustSingleGL( data_ptr->counts, neg_ptr, special_ptr, 0 );
					genome->ptr = special_ptr;
					special_ptr++;
				}
				else if ( het_map != selfHet.end() ) { // in het
//cout << genome->win_index << "\tHET" << endl;
					special_ptr->GL = data_ptr->GL;
					special_ptr->dups = 1;
					special_ptr->counts = data_ptr->counts;
					StatCellPtr het_ptr = refStats[1].begin();
					if ( het_map->second.stat_index >= refStats[1].size() ) {
						cerr << "ERROR: selfHet index: " << het_map->second.stat_index << " > refStats[1] size: " << refStats[1].size() << endl;
						exit(1);
					}
					het_ptr += het_map->second.stat_index;
					adjustSingleGL( data_ptr->counts, het_ptr, special_ptr, 1 );
					genome->ptr = special_ptr;
					special_ptr++;
				}
			}
		}
	}
// resize special
	if ( special_ptr == dataOsPtr->SpecialMergeData.begin() ) {
		dataOsPtr->SpecialMergeData.clear();
	}
	else {
		special_ptr--;
		int actual_size = special_ptr - dataOsPtr->SpecialMergeData.begin() + 1;
		if ( actual_size < 0 ) {
			cerr << "ERROR: special_ptr out of boundary!" << endl;
			exit(1);
		}
		if ( actual_size > int( dataOsPtr->SpecialMergeData.size() ) ) {
			cerr << "ERROR: Actual special size > reserved size!" << endl;
			exit(1);
		}
		dataOsPtr->SpecialMergeData.resize( actual_size );
	}
	if (DEBUG_MODE)
		cout << "Special Merge Data size = " << dataOsPtr->SpecialMergeData.size() << endl;
}

// mark ref_lh as done
void RefStats::MarkRefLHasDone()
{
	ref_lh_set = 1;
}

// remove slices from selfHom & selfNeg: MUST DO THIS AFTER REF LH IS SET!
void RefStats::ReAdjustSelfsWithLiftOver()
{
	if ( !ref_lh_set ) {
		cerr << "Warning: ref lh is not set. Are you sure you want to skip ref lh?" << endl;
	}
// simply clear selfHom
	selfHom.clear();
// update neg map keys by lift over
	updateSelfSingleMapKeys( selfNeg );
	updateSelfSingleMapKeys( selfHet );
}


// update keys by lift over
// way to solve collision: if find same key, add the new to extra; then append extra to selfNeg
  // why? since no 2 new key should be the same
void RefStats::updateSelfSingleMapKeys( SelfSingleMap & selfMap )
{
	SelfSingleMap extra_map;
	vector<int> old_keys;
	for ( SelfSingleMap::iterator map_it = selfMap.begin(); map_it != selfMap.end(); map_it++ )
		old_keys.push_back( map_it->first );
  // add new keys
	for( vector<int>::iterator val = old_keys.begin(); val != old_keys.end(); val++ ) {
		SelfSingleMap::iterator old_it = selfMap.find( *val );
		int new_key = *val + old_it->second.lift;
		if ( selfMap.find(new_key) != selfMap.end() ) { // collision
			extra_map[ new_key ] = old_it->second;
		}
		else { // no, add directly
			selfMap[ new_key ] = old_it->second;
		}
	}
  // delete old keys
  	for( vector<int>::iterator val = old_keys.begin(); val != old_keys.end(); val++ ) {
  		selfMap.erase( *val );
  	}
	
	if ( DEBUG_MODE ) {
		cout << "ReAdjustSelfsWithLiftOver: Collision extra map size = " << extra_map.size() << endl;
	}
	if ( extra_map.size() == 0 ) // no collision to process
		return;
  // append extra to selfNeg	
	for( SelfSingleMap::iterator extra = extra_map.begin(); extra != extra_map.end(); extra++ ) {
		selfMap[ extra->first ] = extra->second;
	}
}


// if return true, already redirect to special merge data
void RefStats::adjustSingleGL( vector<int> & counts, StatCellPtr & stat_ptr, MergeCellPtr & special_ptr, int s_ref )
{
	float compensate = GetGLfromCounts( counts, stat_ptr->log_frac );
	float original =  special_ptr->GL[ s_ref ] + log( RefCounts[ s_ref ] );
	if ( original - compensate <= 3 ) { // set new GL
		setSingleRefGLWithExclusion( counts, stat_ptr, special_ptr, s_ref);
	}
	else { // minus GL
		float new_gl = MinusGL( original, compensate );
		new_gl -= log(RefCounts[ s_ref ] - 1);
		special_ptr->GL[ s_ref ] = new_gl;	
	}
}


// used in adjustSinlgeGL when need to recalculate
void RefStats::setSingleRefGLWithExclusion( vector<int> & counts, StatCellPtr & stat_ptr, MergeCellPtr & special_ptr, int s_ref )
{
	float conjugate_gl = 0;
	bool skipped = 0;
	bool is_first = 1;
	for( vector<StatCell>::iterator upper_it = refStats[ s_ref ].begin(); upper_it != refStats[s_ref].end(); upper_it++) {
	  // see if skipped
		if (!skipped) {
			if ( upper_it == stat_ptr ) { // skip exclusion
				skipped = 1;
				continue;
			}
		}
		float single_gl = GetGLfromCounts( counts, upper_it->log_frac );
		if ( s_ref == 0 ) // for neg, also multiply #dups as weight
			single_gl += log( upper_it->dups );
// add 2 together
		if ( is_first ) {
			conjugate_gl = single_gl;
			is_first = 0;
		}
		else
			conjugate_gl = SumGL( conjugate_gl, single_gl );
	}
	conjugate_gl -= log( RefCounts[ s_ref ] - 1 ); // minus the exclude
	special_ptr->GL[ s_ref ] = conjugate_gl;	
}

// adjust hom_rec with multiple het
void RefStats::adjustMultipleGL( vector<int> & counts, vector<unsigned int> & hets, MergeCellPtr & special_ptr )
{
	vector<unsigned int>::iterator it_het = hets.begin();
// sum all compensates
	float compensate;
	bool is_first = 1;
	for( ; it_het != hets.end(); it_het++ ) {
		float single_gl = GetGLfromCounts( counts, refStats[1][*it_het].log_frac );
		if ( is_first ) {
			compensate = single_gl;
			is_first = 0;	
		}
		else
			compensate = SumGL( compensate, single_gl );
	}
// see if need to re-calculate
	float original =  special_ptr->GL[1] + log( RefCounts[1] );
	if ( original - compensate > 3 ) { // only do minus
		float new_gl = MinusGL( original, compensate );	
		new_gl -= log( RefCounts[1] - 10 );
		special_ptr->GL[1] = new_gl;
		return;
	}

// add now
	unsigned int dist = 0;
	float conjugate_gl;
	is_first = 1;
	for( vector<StatCell>::iterator upper_it = refStats[1].begin(); upper_it != refStats[1].end(); dist++, upper_it++) {
	// check if exclude
		if (it_het != hets.end()) {
			if ( dist == *it_het ) {
				it_het++;
				continue;
			}
		}
		float single_gl = GetGLfromCounts( counts, upper_it->log_frac );
    // add 2 together
		if ( is_first ) {
			conjugate_gl = single_gl;
			is_first = 0;
		}
		else
			conjugate_gl = SumGL( conjugate_gl, single_gl );		
	}
	conjugate_gl -= log(RefCounts[1] - 10); // counts
	special_ptr->GL[1] = conjugate_gl;
}

/************** end *****************/


/******* other inner functions ******/

// only for a ref in 0~2
void RefStats::setSingleRefGL( MergeCellPtr & merge, int & s_ref )
{
	float conjugate_gl;
	bool is_first = 1;
	for( vector<StatCell>::iterator upper_it = refStats[ s_ref ].begin(); upper_it != refStats[s_ref].end(); upper_it++) {
		float single_gl = GetGLfromCounts( merge->counts, upper_it->log_frac );
		if ( s_ref == 0 ) // for neg, also multiply #dups as weight
			single_gl += log( upper_it->dups );
// add 2 together
		if ( is_first ) {
			conjugate_gl = single_gl;
			is_first = 0;
		}
		else
			conjugate_gl = SumGL( conjugate_gl, single_gl );
	}
	conjugate_gl -= log(RefCounts[ s_ref ]);
	merge->GL[ s_ref ] = conjugate_gl;
}

/**** functions used in build ref stats *****/

// set hom & ctrl by label and re-direct
// method: // re-forge genome link to a index-based version
//			generate hom ref by copying
//			generate het ref by averaging
//			eliminate nearby windows
// 			re-direct neg
void RefStats::setStatRef( AllHetIndex & allHetIndex )
{
	generateIndexBasedGenomeLink();
	setHomAndHetRef( allHetIndex );
	setNegRef( allHetIndex );
// set refCounts
	RefCounts.resize(3);
	RefCounts[2] = refStats[2].size();
	RefCounts[1] = refStats[1].size();
	int sum = 0;
	for( vector< StatCell>::iterator it = refStats[0].begin(); it != refStats[0].end(); it++ ) {
		sum += it->dups;
	}
	RefCounts[0] = sum;
	cout << "Final 0/0 ref size = " << refStats[0].size() << ", counts = " << RefCounts[0] << endl;
}

// copy genome link to here, generate index based
// also convert count to log
void RefStats::generateIndexBasedGenomeLink()
{
// copy merge data	
	copiedStats.resize( osPtr->MergeData.size() );
	MergeCellPtr md_it = osPtr->MergeData.begin();
	for( StatCellPtr stat = copiedStats.begin(); stat != copiedStats.end(); stat++, md_it++ ) {
		stat->dups = md_it->dups;
		int sum = GetVecSum( md_it->counts );
		if (sum == 0) {
			cerr << "ERROR: empty cell labeled as exist!" << endl;
			exit(1);
		}
		stat->log_frac.resize( md_it->counts.size() );
	  // get log_frac
	  	int zeros = GetNumberOfZerosInVec( md_it->counts );
	  	float sum_float = float(sum + CountOffSet * zeros);
		for( int i = 0; i < int(md_it->counts.size()); i++ ) {
			stat->log_frac[i] = md_it->counts[i] > 0 ? log( md_it->counts[i] / sum_float ) : log( CountOffSet / sum_float );
		}
	}
	
// reset genome link
	if ( osPtr->GenomeLocationMap.size() != 1 ) {
		cerr << "ERROR: proper reads in other chromosomes are mixed in ctrl_bam!" << endl;
		exit(1);
	}
	map< string, vector< GenomeLocationCell > >::iterator map_it = osPtr->GenomeLocationMap.begin();
	if ( map_it->first.compare(REF_CHR) != 0 ) {
		cerr << "ERROR: ref-chr " << REF_CHR << " is not in genome-location-map!" << endl;
		exit(1);
	}
	if ( map_it->second.size() == 0 ) {
		cerr << "ERROR: empty ref-chr in osPtr!" << endl;
		exit(1);
	}
	GlcPtr vit = map_it->second.end();
	vit--;
	int ref_chr_length = vit->win_index + 1; // this is not the actual chr lenght, just end of data
	indexBasedGenomeLink.resize( ref_chr_length, NullStatCellPtr );
	MergeCellPtr merge_data_begin = osPtr->MergeData.begin();

	for( GlcPtr os_genome = map_it->second.begin(); os_genome != map_it->second.end(); os_genome++ ) {
		int merge_dist = os_genome->ptr - merge_data_begin;
		if ( merge_dist < 0 || merge_dist >= int( copiedStats.size() ) ) {
			cerr << "ERROR: out of copiedStats boundary when linking!" << endl;
			exit(1);
		}
		if ( os_genome->win_index >= int( indexBasedGenomeLink.size() ) ) {
			cerr << "ERROR: os_genome->win_index ( " << os_genome->win_index << ") >= indexBasedGenomeLink size (" << indexBasedGenomeLink.size() << ")!" << endl;
			exit(1);
		}
		indexBasedGenomeLink[ os_genome->win_index ] = ( copiedStats.begin() + merge_dist );
	}
}

// set hom ref by copying
// also destroy link in hom nearby region
// since hom & het in selfHom / selfHet need to add together, otherwise need to search for hom pair when add selfHet
void RefStats::setHomAndHetRef( AllHetIndex & allHetIndex )
{
// REMOVE OTHER MEI SLICE WINDOWS!
	vector<int> OtherIndex;
	_setOtherIndex( OtherIndex );
// destroy link
	int destroy_range = WIN / STEP * 2; // on one side of hom #windows should be destroyed from neg
	for( int i = 0; i <= 1; i++ ) {
		int other_index = OtherIndex[i];
		for( vector<HetRecord>::iterator het_rec = allHetIndex[ other_index ].begin(); het_rec != allHetIndex[ other_index ].end(); het_rec++ ) {
			vector< StatCellPtr >::iterator ptr = indexBasedGenomeLink.begin();
			if ( het_rec->hom.location >= int( indexBasedGenomeLink.size() ) )	{
				cerr << "Warning: hom location: " << het_rec->hom.location << " > genome link size: " << indexBasedGenomeLink.size()<< ", skipped!" << endl;
				continue;
			}
			ptr += het_rec->hom.location;
		// check begin boundary
			int upper;
			if ( het_rec->hom.location < destroy_range ) {
				ptr = indexBasedGenomeLink.begin();
				upper = het_rec->hom.location;	
			}
			else {
				ptr -= destroy_range;
				upper = destroy_range;
			}
		// check end boundary
			int end_dist = (indexBasedGenomeLink.end() - 1 - ptr);
			if ( end_dist <= 2*destroy_range )
				upper += end_dist;
			else
				upper += destroy_range;
			for( int i = 0; i < upper; i++ )
				_destroyGenomeLink( *ptr );
		}	
	}
	
// then add hom
	refStats[2].resize( allHetIndex[mei_index].size() );
	refStats[1].resize( refStats[2].size() * 10 ); // het max is depend on available hom
	StatCellPtr het_add = refStats[1].begin();
	StatCellPtr hom_add = refStats[2].begin();
// selfHom/Het is a map so no need to have pointer or do resize
	int rec_skipped = 0;
	int het_rec_skipped = 0;
	for( vector<HetRecord>::iterator all_rec = allHetIndex[mei_index].begin(); all_rec != allHetIndex[mei_index].end(); all_rec++ ) {
		int hom_index = all_rec->hom.location;
	  // skip hom at no record
		if ( indexBasedGenomeLink[ hom_index ] == NullStatCellPtr ) {
			rec_skipped++;
//			if (DEBUG_MODE)
//				cout << "Skipped hom index at: " << hom_index << endl;
			continue;
		}
		StatCellPtr hom_genome_ptr = indexBasedGenomeLink[ hom_index ];
	  // add to refStats[0]
		hom_add->dups = 1;
		hom_add->log_frac = hom_genome_ptr->log_frac;
	  // add to selfHom
		selfHom[ hom_index ].hom_index = hom_add - refStats[2].begin();
		selfHom[ hom_index ].hets.clear();
	// add het according to hets Loc vec
		for( vector< Loc >::iterator het_cell = all_rec->hets.begin(); het_cell != all_rec->hets.end(); het_cell++ ) {
			int neg_index = het_cell->location;
			StatCellPtr neg_genome_ptr = indexBasedGenomeLink[ neg_index ];
			if ( neg_genome_ptr == NullStatCellPtr ) {
				het_rec_skipped++;
//				if (DEBUG_MODE)
//					cout << "Skipped neg index at het: " << neg_index << endl;
				continue;
			}
		  // add to refStats[1]
			het_add->dups = 1;
			het_add->log_frac.resize( CountSize );
			for( unsigned int i = 0; i < neg_genome_ptr->log_frac.size(); i++ ) {
				het_add->log_frac[i] = _getAvrLogFrac( hom_genome_ptr->log_frac[i], neg_genome_ptr->log_frac[i] );
			}
			int het_ref_index = het_add - refStats[1].begin();
		  // add to selfHom
		  	selfHom[ hom_index ].hets.push_back( het_ref_index );
		  
		  // add to selfHet
			selfHet[ neg_index ].lift = het_cell->lift_over;
			selfHet[ neg_index ].stat_index = het_ref_index;
			
		  // destroy het link
			_destroyGenomeLink( indexBasedGenomeLink[ neg_index ] );
		  // move het ptr
			het_add++;
		}
	// destroy hom & nearby link
		int destroy_start = hom_index >= 5 ? (hom_index - 5) : 0;
		int destroy_end = hom_index + 5 < int( indexBasedGenomeLink.size() ) ? hom_index + 5 : int( indexBasedGenomeLink.size() ) - 1;
		for( int destroy_index = destroy_start; destroy_index < destroy_end; destroy_index++ )
			_destroyGenomeLink( indexBasedGenomeLink[ destroy_index ] );
	// move hom ptr
		hom_add++;
	}

// sanity check	
	if ( hom_add == refStats[0].begin() ) {
		cerr << "Error: empty 1/1 ref. Stop!" << endl;
		exit(2);	
	}
	if ( het_add == refStats[1].begin() ) {
		cerr << "Error: empty 0/1 ref. Stop!" << endl;
		exit(2);	
	}

	het_add--;
	hom_add--;	
	refStats[2].resize( hom_add - refStats[2].begin() + 1 );
	refStats[1].resize( het_add - refStats[1].begin() + 1);

	cout << "Final 1/1 ref size = " << refStats[2].size() << ", skipped " << rec_skipped << " refs." << endl;
	cout << "Final 0/1 ref size = " << refStats[1].size() << ", skipped " << het_rec_skipped << " refs." << endl;
}


// generate neg by going over and de-dup overlapping windows
void RefStats::setNegRef( AllHetIndex & allHetIndex )
{
// destroy and re-adjust first
	int dist = 0; // dist to genome location of last neg record
	int Times = WIN / STEP;

	for( vector< StatCellPtr >::iterator genome = indexBasedGenomeLink.begin(); genome != indexBasedGenomeLink.end(); genome++, dist++ ) {
		if ( *genome == NullStatCellPtr ) {
			continue;
		}
		StatCellPtr ptr = *genome;
		if ( dist < Times ) // destroy
			_destroyGenomeLink( ptr );
		else // can use this
			dist = 0;
	}
	
// add to selfNeg by looping through genome index
// also add copied-stats index to a temporary map
	int genome_index = 0;
	bool reach_end = 0;
	map<int, int> cpConversion;
	vector<HetRecord>::iterator rec = allHetIndex[mei_index].begin();
	for( vector< StatCellPtr >::iterator genome = indexBasedGenomeLink.begin(); genome != indexBasedGenomeLink.end(); genome++, genome_index++ ) {
		if ( *genome == NullStatCellPtr )
			continue;		
	// to selfNeg
		if ( (*genome)->dups == 1 ) {
			if ( !reach_end ) {
				while( genome_index > rec->hom.location ) {
					rec++;
					if ( rec == allHetIndex[mei_index].end() ) {
						reach_end = 1;
						break;
					}
				}
				rec--;
			}
			selfNeg[genome_index].lift = rec->hom.lift_over;
			int copied_index = *genome - copiedStats.begin();
			if ( copied_index < 0 || copied_index >= int(copiedStats.size()) ) {
				cerr << "ERROR: copied_index out of boundary!" << endl;
				exit(1);
			}
			selfNeg[genome_index].stat_index = copied_index;
			cpConversion[ copied_index ] = 0;
		}
	}
	if ( DEBUG_MODE )
		cout << "cpConversion size = " << cpConversion.size() << endl;
// clear genome index
	indexBasedGenomeLink.clear();


// add to refStats[0] by copying copiedStats
// update cpConversion
	refStats[0].resize( copiedStats.size() );
	vector< StatCell>::iterator ref_it = refStats[0].begin();
	int cpIndex = 0;
	for( StatCellPtr cstat = copiedStats.begin(); cstat != copiedStats.end(); cstat++, cpIndex++ ) {
		if (cstat->log_frac.empty())
			continue;
		ref_it->dups = cstat->dups;
		ref_it->log_frac = std::move( cstat->log_frac );
		if ( cpConversion.find( cpIndex ) != cpConversion.end() )
			cpConversion[ cpIndex ] = ref_it - refStats[0].begin();
		ref_it++;
	}
// sanity check
	if ( ref_it == refStats[0].begin() ) {
		cerr << "ERROR: Empty refStats[0]!" << endl;
		exit(2);
	}
	
	ref_it--;
	refStats[0].resize( ref_it - refStats[0].begin() + 1 );

// remove copied stats
	copiedStats.clear();
	
// update selfNeg
	for( SelfSingleMap::iterator map_it = selfNeg.begin(); map_it != selfNeg.end(); map_it++ ) {
		int cp_index = map_it->second.stat_index;
		if ( cpConversion.find( cp_index ) == cpConversion.end() ) {
			cerr << "ERROR: converting index not found in cpConversion!" << endl;
			exit(1);
		}
		map_it->second.stat_index = cpConversion[ cp_index ];
	}
}


// set other index
void RefStats::_setOtherIndex( vector<int> & OtherIndex )
{
	if (mei_index == 0) {
		OtherIndex = {1,2};
	}
	else if ( mei_index == 1 ) {
		OtherIndex = {0,2};
	}
	else { //2
		OtherIndex = {0,1};
	}
}

// get avr log of log_frac
float RefStats::_getAvrLogFrac( float hom_log, float neg_log )
{
	float avrLog;
	avrLog = log( exp(hom_log) + exp(neg_log) ) - log(2);
/*
	if ( neg_log == min_log )
		avrLog = round(hom_log - log(2));
	else if ( hom_log == min_log )
		avrLog = round(neg_log - log(2));
	else {
		avrLog = round(log( exp(hom_log) + exp(neg_log) ) - log(2)); // here since it's 0-1, no need to care about log = Inf;
	}
	if (avrLog < min_log) // make sure min is min_log
		return min_log;
	else
		return avrLog;
*/
	return avrLog;
}

// use to destroy link
void RefStats::_destroyGenomeLink( StatCellPtr & ptr )
{
// see if it's empty
	if ( ptr == NullStatCellPtr )
		return;

	if (ptr->dups == 0) {
		cerr << "ERROR: dups = 0, something goes wrong in indexBasedGenomeLink";
		if ( ptr->log_frac.empty() )
			cerr << ". Log frac is cleared!" << endl;
		else
			cerr << ". Log frac is not cleared!" << endl;
		exit(1);
	}
	ptr->dups--;
	if ( ptr->dups == 0 ) {
		ptr->log_frac.clear(); // not deleted completely but copied stats will be cleared in the end so doesn't matter
	}
	ptr = NullStatCellPtr;
}


/*************** end *******************/


/*** read het index only once ****/

// store them in AllHetIndex
// when create a new ref stats, create refs from AllHetIndex
// free AllHetIndex when discover phase is done
void SetAllHetIndex( string & het_index_name, AllHetIndex & allHetIndex )
{
	using std::getline;
	allHetIndex.resize(3);
	
// 1st char identifier
	map< char, int > identifier_map;
	identifier_map['A'] = 0;
	identifier_map['L'] = 1;
	identifier_map['S'] = 2;

// read file
	std::ifstream het_index_file;
	het_index_file.open( het_index_name.c_str() );
	CheckInputFileStatus( het_index_file, het_index_name.c_str() );
	string line;
	vector< Loc >::iterator ptr_add_het;
	int dist = 0; // to make sure ptr_add_het do not go out of .end() boundary
	while( getline( het_index_file, line ) ) {
		std::stringstream ss;
		ss << line;
		string ref_name;
		getline( ss, ref_name, '\t');
		string loc_str;
		getline( ss, loc_str, '\t' );
		string ref_type_str;
		getline( ss, ref_type_str, '\t' );
		string lift_over_str;
		getline( ss, lift_over_str, '\t' );
		if ( std::stoi(ref_type_str) == 0 ) { // start a new record
			int mei_type = identifier_map[ ref_name[0] ];
			HetRecord het_record;
			het_record.hom.location = GetHomCoordIndex( stoi( loc_str ) );
			het_record.hom.lift_over = stoi( lift_over_str );
			allHetIndex[mei_type].push_back( het_record );
			vector< HetRecord >::iterator vit = allHetIndex[mei_type].end();
			vit--;
			vit->hets.resize(10);
			ptr_add_het = vit->hets.begin();
			dist = 0;
		}
		else { // het, add according to the last one
			if ( dist >= 10 ) {
				cerr << "ERROR: #het record > designated number!" << endl;
				exit(1);
			}
			ptr_add_het->location = GetHetCoordIndex( stoi( loc_str ) );
			ptr_add_het->lift_over = stoi( lift_over_str );
			ptr_add_het++;
			dist++;
		}
	}
	het_index_file.close();
}


// useful functions

// location is center
int GetHomCoordIndex( int location )
{
	int index;
	if ( location < WIN / 2 )
		index = 0;
	else {
		index = round( float(location - WIN / 2) / STEP );
	}
	return index;
}

// location is left boundary
int GetHetCoordIndex( int location )
{
	return (location / STEP );
}


/**************************  debug functions ************************/
// print out ref stats for debug purpose
void RefStats::PrintRefStats( string & out_prefix )
{
	for(int s_ref = 0; s_ref <= 2; s_ref++) {
		string s_ref_name = out_prefix + "." + std::to_string(s_ref);
		std::ofstream s_ref_file;
		s_ref_file.open( s_ref_name.c_str() );
		CheckOutFileStatus( s_ref_file, s_ref_name.c_str() );
		for( vector< StatCell >::iterator sp_it = refStats[s_ref].begin(); sp_it != refStats[s_ref].end(); sp_it++ ) {
			s_ref_file << sp_it->dups;
			for( vector<float>::iterator it = sp_it->log_frac.begin(); it != sp_it->log_frac.end(); it++ )
				s_ref_file << "\t" << (*it);
			s_ref_file << endl;
		}
		s_ref_file.close();	
	}
}
