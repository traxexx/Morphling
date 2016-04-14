#ifndef MEISEQ_H
#define MEISEQ_H

#include <string>
#include <vector>

using std::string;
using std::vector;

void SetMeiSeqAndNames( string & filename, vector<string> & me_names, vector<string> & me_seqs );

string GetMeiNameFromIndex( int m );

int GetMeiIndexFromName( string & name );

#endif