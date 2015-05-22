#include "FastaUtil.h"
#include "ssw_cpp.h"
#include <fstream>

using namespace std;

int match = 2;

void RefSeq::ReadSeq(const char * meiName, float map_ratio, float mismatch_ratio)
{
	MAP_RATIO_ = map_ratio;
	MISMATCH_RATIO_ = mismatch_ratio;
	string chrName="", line, meiRef="";
	int ix=1;
	bool firstLine = 1;
	ifstream MeiFA;
	MeiFA.open(meiName);
	if(MeiFA.is_open())
	{
		while(MeiFA.good())
		{
			getline(MeiFA, line);
			if (line.size() == 0)
				continue;
			if (line.at(0) == '>')
			{
				if (!firstLine)
				{
					SeqH[chrName] = meiRef;
					ix = 1; chrName = ""; meiRef = "";
				}
				else
					firstLine = 0;
				while((line.at(ix) != ' ') && (line.at(ix) != '\t') && (line.at(ix) != '\n'))
				{
					chrName += line.at(ix); ix++;
				}
			}
			else
			{
				meiRef += line.substr(0, line.size());
			}
		}
		SeqH[chrName] = meiRef;
		MeiFA.close();
	}
	else std::cerr << "Unable to open " << meiName << endl;
}


void RefSeq::addPolyAtail()
{
	std::string polyA = std::string(30,'a');
	std::string tail;
	// since all converted to lower, this transform is unecessary
	for(seq_it iter = SeqH.begin(); iter != SeqH.end(); iter++)
	{
		tail = iter->second.substr(iter->second.length()-21,20);
		transform(tail.begin(), tail.end(), tail.begin(), ::tolower);
		if (tail.compare(0,20,polyA,0,20) != 0)
			iter->second.append(polyA);
	}
}

void RefSeq::printAll()
{
	for(seq_it iter = SeqH.begin(); iter != SeqH.end(); ++iter)
	{
		cout << iter->first << "  " << iter->second << endl;
	}
}

bool RefSeq::MeiMap(string & seq)
{
	int readLen = seq.length();
	int SR = readLen * match * MAP_RATIO_;
	int MM = readLen * MISMATCH_RATIO_;
	
	int SR23 = SR * 0.67;

	for(seq_it iter = SeqH.begin(); iter != SeqH.end(); iter++) {
		int refLen = iter->second.length();
		StripedSmithWaterman::Aligner aligner;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		aligner.Align(seq.c_str(), iter->second.c_str(), refLen, filter, & alignment);
		 if(alignment.sw_score >= SR && alignment.mismatches <= MM)
		 	return 1;
		 else if (alignment.sw_score <= SR23) // no further alignment for low scores
		 	break;
	}
	
	string revSeq(seq.length(),'\0');
	RevComp(seq, revSeq);
	for(seq_it iter = SeqH.begin(); iter != SeqH.end(); iter++) {
		int refLen = iter->second.length();
		StripedSmithWaterman::Aligner aligner;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		aligner.Align(revSeq.c_str(), iter->second.c_str(), refLen, filter, & alignment);
		 if(alignment.sw_score >= SR && alignment.mismatches <= MM)
		 	return 1;
		 else if (alignment.sw_score <= SR23) // no further alignment for low scores
		 	break;
	}
	
	return 0;
}

void RefSeq::RevComp(string & seq, string & rev)
{
        int len = seq.length();
        if (len < 2)
	{
               rev = string("N"); return;
	}
        else
        {
            for (int ix=0; ix<len; ix++)
            	rev[ix] = CompNt(seq[len-ix-1]);
        }
}

char RefSeq::CompNt(char & sx)
{
        char c;
        switch(sx)
        {
                        case ('A'):
                                c = 'T'; break;
                        case ('T'):
                                c = 'A'; break;
                        case ('C'):
                                c = 'G'; break;
                        case ('G'):
                                c = 'C'; break;
                        case ('a'):
                                c = 'T'; break;
                        case ('t'):
                                c = 'A'; break;
                        case ('c'):
                                c = 'G'; break;
                        case ('g'):
                                c = 'C'; break;
                        default:
                                c = 'N'; break; 

        }
        return c;
}

































