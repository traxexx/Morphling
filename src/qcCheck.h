#ifndef QCCHECK_H
#define QCCHECK_H

#include <string>

using std::string;

// return true if string is a float number
bool str_is_float( string & str );

int GetFileLineCount(string & ref_file_name);

bool IsEmptyFile( string & name );

#endif
