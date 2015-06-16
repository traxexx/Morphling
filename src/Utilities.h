#ifndef UTILITIES_H
#define UTILITIES_H

#include <fstream>
#include <string>
#include <vector>

using std::vector;
using std::string;

void CheckInputFileStatus(std::ifstream & infile, const char* name);

void CheckOutFileStatus(std::ofstream & outfile, const char* name);

std::string getMeiTypeString(int & mei_type);

int getSumOfVector( vector<int> & count_table );

int GetAvrLocationOfCertainValue( vector<int> & vec, int & val, int & min_index, int & max_index );

int getSumOfMeiReads( vector<int> & count_table );

void ExecuteCmd( std::string & cmd );

void GenerateDoneFile( std::string & work_dir, const char * prefix );

bool ExistDoneFile( std::string & work_dir, const char * prefix );

string GetRemapCmd( string & full_mapper, std::string fastq_prefix, std::string ref_fasta, std::string remapSam );

string GetFileBaseName( string & full_name );

string GetFileDirName( string & full_name );

void MakeDirectory( string & dir_name );

string GetExePath();

int GetTabLocation( int search_start, int noccur, string & line );

int CountContinuousChar( string & seq, char c, int n);

#endif