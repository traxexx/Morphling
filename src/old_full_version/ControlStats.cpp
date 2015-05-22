#include "ControlStats.h"
#include <iostream>
#include <math.h> // log
#include <algorithm> // find

bool sortLiftOver( std::pair<int, int> x, std::pair<int, int> y )
{
	if (x.first < y.first)
		return 1;
	else if (x.first > y.first)
		return 0;
	else {
		std::cerr << "ERROR: Lift Over index is equal: " << x.first << "." << std::endl; exit(1);
	}
}

// use to convert Hom & Neg
void convertStatCellVecFromCountToLog( std::vector< std::vector<StatCell>::iterator > & sv )
{
	for( std::vector< std::vector<StatCell>::iterator >::iterator sv_it = sv.begin(); sv_it != sv.end(); sv_it++ ) {
		convertIntVecFromCountToLog( (*sv_it)->log_table );
	}
}

// use in *Vector & building het ctrl
void convertIntVecFromCountToLog( std::vector<int> & vec )
{
	float sum = float(vec.size()) / 2;
	for( std::vector<int>::iterator it = vec.begin(); it != vec.end(); it++ )
		sum += *it;
	for( std::vector<int>::iterator it = vec.begin(); it != vec.end(); it++ ) {
		*it = -round( log((*it +0.5) / sum) );
	}
}

// overlap

bool ExistOverlap( std::vector< std::vector<StatCell>::iterator > & ctrl, std::map< std::vector<StatCell>::iterator, int > & exclude_stats)
{
	for( std::vector< std::vector<StatCell>::iterator >::iterator cit = ctrl.begin(); cit != ctrl.end(); cit++) {
		if ( exclude_stats.find(*cit) != exclude_stats.end() )
			return 1;
	}
	return 0;
}


bool ExistHetOverlap( std::vector<HetCtrlCell> & ctrl, std::map< std::vector<StatCell>::iterator, int > & exclude_stats)
{
	for( std::vector<HetCtrlCell>::iterator het_it = ctrl.begin(); het_it != ctrl.end(); het_it++ ) {
		if ( exclude_stats.find( het_it->hom_trace ) != exclude_stats.end() )
			return 1;
		if ( exclude_stats.find( het_it->neg_trace ) != exclude_stats.end() )
			return 1;			
	}
	return 0;
}



