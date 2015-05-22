#ifndef MAPDATA_H
#define MAPDATA_H

#include <vector>
#include <string>

using std::vector;

/* Raw data structure for generating read type count */
typedef struct {
	vector<int> regular;
	vector<int> Alu;
	vector<int> L1;
	vector<int> SVA;
} ProperMapCell; // win_index, count_table

typedef struct {
	vector<int> Alu;
	vector<int> L1;
	vector<int> SVA;
} DiscMapCell;


typedef struct {
	std::string chr;
	int win_index;
	vector<int> stats;
	bool isMEI;
} ReadMapCell;

/* Final data structure, after sort */
typedef struct {
	vector<int> count_table;
	vector<int> PL; // LH result
} MergeCell;


typedef struct {
	int win_index;
	std::vector<MergeCell>::iterator data;
} GenomeLocationCell;

// read data from output stats
bool setPerChrRawStats( std::vector<ReadMapCell> & raw_data, std::string & regular_name, std::string & clip_name, std::string & disc_name, bool skip_nonMEI);

// set duplicate vec for ReadMapCell vector
void setDupVecForReadMapCellVec( std::vector<int> & dup, std::vector<ReadMapCell> & raw_data );

// sort link (for output purpose)
bool sortGenomeLocationLink( GenomeLocationCell x, GenomeLocationCell y );

// check if disc map unit is empty
bool isEmptyDiscRawMapUnit( std::vector<DiscMapCell> & dv );

#endif