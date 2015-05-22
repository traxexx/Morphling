#ifndef PL_H
#define PL_H

#include <vector>
#include <iterator>
#include <string>
#include <map>

#include "ControlStats.h"
#include "MapData.h"

using std::string;
using std::vector;
using std::map;

typedef struct {
	int index; // which window
	std::vector<int> PL;
} SpecialPLcell;


int getPerRefPerMergeCellPL( vector<MergeCell>::iterator & current_merge_it, vector< vector<StatCell>::iterator > & ctrl,
	map< vector<StatCell>::iterator, int > & exclude_stats, bool take_count);

int getPerRefPerMergeCellHetPL( vector<MergeCell>::iterator & current_merge_it, vector<HetCtrlCell> & HetCtrl, map< vector<StatCell>::iterator, int > & exclude_stats);

int calculatePLfromSingleIntVec( std::vector<int> & data, std::vector<int> & ref );

int calculateSumOfPLs(int & val1, int & val2);

string getGenotypeFromPL(vector<int> & pl);

void setNormalizedPLfromPL( vector<int> & pl, vector<int> & ref_pl) ;

#endif


			