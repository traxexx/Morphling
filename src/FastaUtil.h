#ifndef FASTAUTIL_H
#define FASTAUTIL_H


#include <iostream>
#include <algorithm>
#include <string>
#include <map>

typedef std::map<std::string, std::string> SeqHash;
typedef SeqHash::iterator seq_it;

class RefSeq
{
public:
	void ReadSeq(std::string & meiName);
	void addPolyAtail();
	void printAll();
	bool MeiMap(std::string & seq);
	
private:
	SeqHash SeqH;
	bool singlePartMap( std::string & seq );
	void RevComp(std::string & seq, std::string & rev);
	char CompNt(char & sx);	
};

#endif

