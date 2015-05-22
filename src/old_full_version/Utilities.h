#include <fstream>
#include <string>
#include <vector>

using std::vector;

void CheckFileStatus(std::ifstream & infile);

void CheckOutFileStatus(std::ofstream & outfile);

std::string getMeiTypeString(int & mei_type);

int getSumOfVector( vector<int> & count_table );

int getSumOfMeiReads( vector<int> & count_table );

void CheckCmdStatus(std::string & cmd, int & cmd_status);

void GenerateDoneFile( std::string & work_dir, const char * prefix );

bool ExistDoneFile( std::string & work_dir, const char * prefix );


