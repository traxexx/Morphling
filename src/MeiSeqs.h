#ifndef MEISEQS_H
#define MEISEQS_H

#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::map;

// for storing mapping info of each read on each mei subtype
struct EviInfo
{
	int Boundary; // if ==2, do not consider any mapping results on left; ==1 no right. ==0 no boundary.
	int LAlign;
	int RAlign;
	int Score; // SW score
	int SeqKey;
	int SampleKey; // note sample index. Only used in Sites & AsbSites
};

typedef map<int, vector<EviInfo> > subCluster;


void loadMEsequence( string & me_list_name, vector< vector<string> > & MEseqs, vector< vector<string> > & MEnames );

string ReverseCompSeq( string & seq );

char GetCompNt( char nt);

int GetMEtypeFromAlt( string & field);

#endif
