#include <iostream>
#include <stdlib.h>
#include "Utilities.h"
#include <iterator>

using std::vector;

void CheckFileStatus(std::ifstream & infile)
{
	if (!infile.is_open()) {
		std::cerr << "ERROR: Can't open input file " << infile << std::endl;
		exit(1);
	}
}

void CheckOutFileStatus(std::ofstream & outfile)
{
	if (!outfile.is_open()) {
		std::cerr << "ERROR: Can't open output file " << outfile << std::endl;
		exit(1);
	}
}


std::string getMeiTypeString(int & mei_type)
{
	std::string MeiType;
	switch( mei_type ) {
		case 0:
			MeiType = std::string("Alu"); break;
		case 1:
			MeiType = std::string("L1"); break;
		case 2:
			MeiType = std::string("SVA"); break;
		default:
			std::cerr << "ERROR: illegal mei_type in ReadMap::PrintToVcf!" << std::endl; exit(1);
	}
	return MeiType;
}

// use as DP in vcf
int getSumOfVector( vector<int> & count_table )
{
	int sum = 0;
	for( vector<int>::iterator it = count_table.begin(); it != count_table.end(); it++ )
		sum += *it;
	return sum;
}

// use as GQ in vcf
/*
0 proper 1 short [4 5] [8 9] (clip) [12 13] [16 17] (disc, unmap)
*/
int getSumOfMeiReads( vector<int> & count_table )
{
	vector<int> asMei = {4,5,8,9,12,13,16,17};
	int sum = 0;
	for( vector<int>::iterator it = asMei.begin(); it != asMei.end(); it++ )
		sum += count_table[*it];
	return sum;
}


void CheckCmdStatus(std::string & cmd, int & cmd_status)
{
	if (cmd_status != 0) {
		std::cerr << "ERROR: Fail to run: " << cmd << std::endl;
		std::cerr << "    Exit " << cmd_status << std::endl;
		exit(1);
	}
}


void GenerateDoneFile( std::string & work_dir, const char * prefix )
{
	std::string touch_cmd = std::string("touch ") + work_dir + "/" + std::string(prefix) + ".Done";
	int touch_cmd_status = system(touch_cmd.c_str());
	CheckCmdStatus(touch_cmd, touch_cmd_status);
}

bool ExistDoneFile( std::string & work_dir, const char * prefix )
{
	std::string fname = work_dir + "/" + std::string(prefix) + ".Done";
	std::ifstream dfile(fname.c_str());
	if ( dfile.good() ) {
		dfile.close();
		std::cout << "Warning: exist " << fname << ", skip related steps!" << std::endl;
		return 1;
	}
	return 0;
}


























