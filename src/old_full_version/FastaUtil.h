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
	void ReadSeq(const char * meiName, float map_ratio, float mismatch_ratio);
	void addPolyAtail();
	void printAll();
	bool MeiMap(std::string & seq);
	
private:
	SeqHash SeqH;
	void RevComp(std::string & seq, std::string & rev);
	char CompNt(char & sx);
	
	float MAP_RATIO_;
	float MISMATCH_RATIO_;
};

#endif

