#include "PL.h"
#include <algorithm> // min_element
#include <iostream>
#include <map>
#include <iterator>

using std::vector;
using std::map;

int getPerRefPerMergeCellPL( 
	vector<MergeCell>::iterator & current_merge_it, vector< vector<StatCell>::iterator > & ctrl,
	map< vector<StatCell>::iterator, int > & exclude_stats, bool take_count)
{
	int lh;
	int skipped_count = 0;
	map< vector<StatCell>::iterator, int > Exclude_stats = exclude_stats;
	for( vector< vector<StatCell>::iterator >::iterator it = ctrl.begin(); it != ctrl.end(); it++ ) {
		map< vector<StatCell>::iterator, int >::iterator find_it = Exclude_stats.find(*it);
		int count_offset = 0;
		if ( find_it != Exclude_stats.end() ) {
			if ( find_it->first->count == 1 ) {
				skipped_count++; Exclude_stats.erase(find_it);
				continue;
			}
			else {
				if (take_count) { // minus count, do lh still
					count_offset = find_it->second;
					skipped_count += find_it->second;
				}
				else { // skip this record (hom)
					skipped_count++; continue;
				}
			}
		}
		int current_lh = calculatePLfromSingleIntVec( current_merge_it->count_table, (*it)->log_table);
		if (take_count)	lh -= log((*it)->count - count_offset);
		if (it == ctrl.begin()) {
			lh = current_lh; continue;
		}
		lh = calculateSumOfPLs( lh, current_lh );
	}
// decide final size	
	int final_size;
	if (take_count) {
		final_size = 0;
		for( vector< vector<StatCell>::iterator >::iterator it = ctrl.begin(); it != ctrl.end(); it++ )
			final_size += (*it)->count;
		final_size -= skipped_count;
	}
	else
		final_size = ctrl.size() - skipped_count;
	
	if (final_size == 0) {
		std::cerr << "ERROR: no ctrl available in setPerRefPerMergeCellPL()... Something must be wrong!" << std::endl; exit(2);
	}
		
	lh += log( final_size );
	return lh;
}

int getPerRefPerMergeCellHetPL( vector<MergeCell>::iterator & current_merge_it, vector<HetCtrlCell> & HetCtrl, map< vector<StatCell>::iterator, int > & exclude_stats)
{
	int lh;
	int skipped = 0;
	for( vector<HetCtrlCell>::iterator it = HetCtrl.begin(); it != HetCtrl.end(); it++ ) {
	
		map< vector<StatCell>::iterator, int >::iterator find_it = exclude_stats.find( it->hom_trace );
		if (find_it != exclude_stats.end()) {
			skipped++; continue;
		}
		find_it = exclude_stats.find( it->neg_trace );
		if (find_it != exclude_stats.end()) {
			skipped++; continue;
		}
				
		int current_lh = calculatePLfromSingleIntVec( current_merge_it->count_table, it->log_table );
		if (it == HetCtrl.begin()) {
			lh = current_lh; continue;
		}
		lh = calculateSumOfPLs( lh, current_lh );
	}
	
	int final_size = HetCtrl.size() - skipped;
	if (final_size == 0) {
		std::cerr << "ERROR: no available HetCtrl in getPerRefPerMergeCellHetPL()... Something must be wrong!" << std::endl; exit(2);
	}
	lh += ( final_size );
	return lh;
}


int calculatePLfromSingleIntVec( vector<int> & data, vector<int> & ref )
{
	int single_lh = 0;
	vector<int>::iterator ref_it = ref.begin();
	for( vector<int>::iterator data_it = data.begin(); data_it != data.end(); data_it++ ) {
		single_lh += (*data_it) * (*ref_it);
		ref_it++;
	}
	return single_lh;
}

int calculateSumOfPLs(int & val1, int & val2)
{
	if (val1 - val2 >= 5)
		return val2;
	else if (val1 - val2 >= -5) {
		int mid = (val1 + val2) / 2;
		int new_val = mid - log( exp(-val1 + mid) + exp(-val2 + mid) );
		return new_val;
	}
	else return val1;
}

// for print gt in vcf
string getGenotypeFromPL(vector<int> & pl)
{
	std::string gt;
	vector<int>::iterator it = std::min_element(pl.begin(), pl.end()); // remember pl is actually negative
	int dist = it - pl.begin();
	switch(dist) {
		case(0):
			gt = std::string("1/1"); break;
		case(1):
			gt = std::string("0/1"); break;
		case(2):
			gt = std::string("0/0"); break;
		default:
			std::cerr << "ERROR: getGenotypeFromPL() irregular max_element!" << std::endl; exit(2);
	}
	return gt;
}

// for print pl in vcf
void setNormalizedPLfromPL( vector<int> & pl, vector<int> & ref_pl)
{
	if (pl.size() != 3)
		pl.resize(3);
	int min_pl = *std::min_element(ref_pl.begin(), ref_pl.end());
	vector<int>::iterator rit = ref_pl.begin();
	for( vector<int>::iterator pit = pl.begin(); pit != pl.end(); pit++ ) {
		*pit = *rit - min_pl;
		rit++;
	}
}






































