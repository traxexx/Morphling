#ifndef CONTROLSTATS_H
#define CONTROLSTATS_H


#include <utility>
#include <vector>
#include <iterator>
#include <map>

typedef struct {
	int count;
	std::vector<int> log_table;
} StatCell;

typedef struct {
	std::vector<int> log_table;
	std::vector<StatCell>::iterator hom_trace;
	std::vector<StatCell>::iterator neg_trace;
} HetCtrlCell;

bool sortLiftOver(std::pair<int, int> x, std::pair<int, int> y);

void convertStatCellVecFromCountToLog( std::vector< std::vector<StatCell>::iterator > & sv );
void convertIntVecFromCountToLog( std::vector<int> & vec );

// check if ctrl & exclude stats have overlap
bool ExistOverlap( std::vector< std::vector<StatCell>::iterator > & ctrl, std::map< std::vector<StatCell>::iterator, int > & exclude_stats);

bool ExistHetOverlap( std::vector<HetCtrlCell> & ctrl, std::map< std::vector<StatCell>::iterator, int > & exclude_stats);

#endif

