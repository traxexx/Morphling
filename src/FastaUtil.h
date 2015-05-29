#ifndef FASTAUTIL_H
#define FASTAUTIL_H

#include <iostream>
#include <algorithm>
#include <string>
#include <map>

using std::string;
using std::map;

typedef std::map<string, string> SeqHash;
typedef SeqHash::iterator seq_it;

class RefSeq
{
public:
	RefSeq();
	~RefSeq();
	void ReadSeq( string & meiName );
	void addPolyAtail();
	void printAll();
	bool MeiMap( string & seq );
	int GetSeqHashSize(); // debug function
	
private:
	SeqHash SeqH;
	bool singlePartMap( string & seq );
	void RevComp( string & seq, string & rev );
	char CompNt( char sx );	
};

#endif

