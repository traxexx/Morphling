#include "FastaUtil.h"
#include "MapParameters.h"
#include "ssw_cpp.h"
#include <fstream>
#include <vector>

using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

RefSeq::RefSeq()
{
	SeqH.clear();
}

RefSeq::~RefSeq(){}

void RefSeq::ReadSeq( string & meiName )
{
	string chrName;
	string line;
	string meiRef;
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
					meiRef = string("");
				}
				else
					firstLine = 0;
				chrName = line.substr(0, line.size() - 2);
/*				while((line.at(ix) != ' ') && (line.at(ix) != '\t') && (line.at(ix) != '\n'))
				{
					chrName += line.at(ix); ix++;
				}
*/
			}
			else
			{
				meiRef += line.substr(0, line.size());
			}
		}
		SeqH[chrName] = meiRef;
		MeiFA.close();
	}
	else cerr << "Unable to open " << meiName << endl;
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

// function for doing mapping
bool RefSeq::singlePartMap( string & seq )
{
	int readLen = seq.length();
	int SR = readLen * match * MAP_RATIO;
	int MM = readLen * MISMATCH_RATIO;
	
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
		 else if (alignment.sw_score < SR23) // no further alignment for low scores
		 	break;
	}
	
	return 0;
}

// turn on part mapping
bool RefSeq::MeiMap(string & seq)
{
	if ( SeqH.empty() ) // if dummy class, nothing will be mapped
		return 0;

	int read_length = seq.length();
// short read map them all
	if ( read_length <= 30 ) {
		bool map_state = singlePartMap( seq );
		return map_state;
	}
	else { // long read split map. If one part map, put into use
		vector< string > vec_str;
		int part_count = read_length / 30 + 1;
		vec_str.resize( part_count );
		for( int i = 0; i < part_count - 1; i++ ) {
			vec_str[i] = seq.substr( i*30, 30 );
		}
		vec_str[ part_count - 1 ] = seq.substr( read_length - 30 );
// map these parts
		int tail_length = part_count * 30 - read_length;
		int map_count;
		if ( tail_length < 15 && part_count > 2 )
			map_count = part_count - 1;
		else
			map_count = part_count;
		for( int i = 0; i < map_count; i++ ) {
			bool map_state = singlePartMap( vec_str[i] );
			if ( map_state )
				return 1;
		}
		return 0;
	}
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

char RefSeq::CompNt(char sx)
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

/*** debug functions ***/
int RefSeq::GetSeqHashSize()
{
	int hsize = SeqH.size();
	return hsize;
}































