#include "MapData.h"
#include <fstream>
#include <iterator>
#include "Params.h"
#include "Utilities.h"
#include <sstream>
#include <iostream>

bool setPerChrRawStats(std::vector<ReadMapCell> & raw_data, std::string & regular_name, std::string & clip_name, std::string & disc_name, bool skip_nonMEI)
{
// first check if we need to skip this chr
	std::ifstream f1(regular_name.c_str());
	std::ifstream f2(clip_name.c_str());
	if (!f1.good() && !f2.good())
		return 0;
	else {
		f1.close(); f2.close(); }

// if skip_nonMEI is false, then doing construction of raw control map. Do not set chr as well (parent function).
// disc first
	std::ifstream disc_file; disc_file.open(disc_name);
	CheckFileStatus(disc_file);
	std::vector<ReadMapCell>::iterator it_raw = raw_data.begin();
	std::string line, field;
	while(std::getline(disc_file, line)) {
		std::stringstream ss; ss << line;
		std::getline(ss, field, '\t');
		it_raw->isMEI = (field.compare("1") == 0) ? 1 : 0;
		it_raw->stats.resize(FIELD_COUNT,0);
		std::vector<int>::iterator nt = it_raw->stats.begin(); nt += 10;
		for(int i=0; i<8; i++, nt++) {
			std::getline(ss, field, '\t');
			*nt = stoi(field);
		}
		it_raw++;
	}
	disc_file.close();

	std::ifstream clip_file; clip_file.open(clip_name);
	CheckFileStatus(clip_file);
	it_raw = raw_data.begin();
	while(std::getline(clip_file, line)) {
		std::stringstream ss; ss << line;
		std::getline(ss, field, '\t');
		if (skip_nonMEI) {
			if (field.compare("1") != 0) {
				if (it_raw->stats.size() > 0 && (!it_raw->isMEI))
					it_raw->stats.clear();
				it_raw++; continue;
			}
		}
		if (it_raw->stats.size() == 0)
			it_raw->stats.resize(FIELD_COUNT);
		std::vector<int>::iterator nt = it_raw->stats.begin(); nt += 2;
		for(int i=0; i<8; i++, nt++) {
			std::getline(ss, field, '\t');
			*nt = stoi(field);			
		}
		it_raw++;
	}
	clip_file.close();
	
	std::ifstream regular_file; regular_file.open(regular_name);
	CheckFileStatus(regular_file);
	it_raw = raw_data.begin();
	while(std::getline(clip_file, line)) {
		std::stringstream ss; ss << line;
		if (skip_nonMEI) {
			if (it_raw->stats.size() == 0) {
				it_raw++; continue;
			}
		}
		std::getline(ss, field, '\t');
		it_raw->stats[0] = stoi(field);
		std::getline(ss, field, '\t');
		it_raw->stats[1] = stoi(field);
		it_raw++;
	}
	regular_file.close();
	return 1;
}


// dup[0] means how many empty
void setDupVecForReadMapCellVec( std::vector<int> & dup, std::vector<ReadMapCell> & raw_data )
{
	std::vector<ReadMapCell>::iterator raw_it = raw_data.begin();
  	int empty_rmc = 0;
  	for(; raw_it != raw_data.end(); raw_it++) {
  		if (raw_it->stats.size() > 0)
  			break;
  		empty_rmc++;
  	}
  	if (raw_it == raw_data.end()) {
  		std::cerr << "ERROR: Empty Raw Map. Anything goes wrong in Re-map? " << std::endl; exit(1);
  	}
  	dup.resize(1, empty_rmc);
  	
  // set duplicates vector
  	int ct = 1;
  	for(std::vector<ReadMapCell>::iterator prev = raw_it; raw_it != raw_data.end(); raw_it++, prev++) {
  		if (prev->stats == raw_it->stats)
  			ct++;
  		else {
  			dup.push_back( ct );
  			ct = 1;
  		}
  	}
  	dup.push_back(ct);
}


bool sortGenomeLocationLink( GenomeLocationCell x, GenomeLocationCell y )
{
	if (x.win_index < y.win_index)
		return 1;
	else if (x.win_index > y.win_index)
		return 0;
	else {
		std::cerr << "ERROR: sortGenomeLocationLink x == y! win_index = " << x.win_index << std::endl; exit(1);
	}
}


bool isEmptyDiscRawMapUnit( std::vector<DiscMapCell> & dv )
{
	for( std::vector<DiscMapCell>::iterator it = dv.begin(); it != dv.end(); it++ )
		if (!(it->Alu.size() < 1 && it->L1.size() < 1 && it->SVA.size() < 1))
			return 0;
	return 1;
}





