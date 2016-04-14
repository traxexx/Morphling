#include "qcCheck.h"
#include <fstream>
#include <algorithm>
#include "morphError.h"

// int is also counted as valid float number
bool str_is_float( string & str )
{
	if ( str[str.size()-1] == '.' )
		return 0;
	bool success = 1;
	bool dot = 0;
	for( int i=0; i<str.size(); i++ ) {
		if ( str[i] == '.' ) {
			if ( !dot ) {
				dot = 1;
				continue;
			}
			else {
				success = 0;
				break;
			}
		}
		if ( !isdigit( str[i] ) ) {
			success = 0;
			break;
		}
	}
	return success;
}

int GetFileLineCount(string & ref_file_name)
{
	std::ifstream file;
	file.open(ref_file_name.c_str());
	if (!file.is_open())
		morphErrorFile( ref_file_name );
	int lc = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	file.close();
	return lc;
}

bool IsEmptyFile( string & name )
{
	bool status;
	std::ifstream file;
	file.open(name.c_str());
	if (!file.is_open())
		morphErrorFile( name );	
	if (file.peek() == std::ifstream::traits_type::eof())
		status = 1;
	else status = 0;
	file.close();
	return status;
}
