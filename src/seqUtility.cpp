#include "seqUtility.h"
#include "Cigar.h"
#include <ctype.h>
#include <algorithm>

using std::cerr;
using std::endl;


// check scheme:
//  if exists 15 consecutive A
//  if first or last 20 nts in read contain >=15 A
bool ContainRepetative( string & seq, char nt )
{
    string pnt;
    pnt.resize(15, nt);
    
    // check polyN first
    if (seq.find(pnt) != std::string::npos)
        return 1;
    
    // if not, then try to find 75% polyNT
    if (seq.length() < 20)
        return 0;

    // if not A, =1; is A, =0
    int head [3];
    int tail [2];
    for(int i=0;i<3; i++) {
    	if (seq[i] != nt)
    		head[i] = 1;
    	else
    		head[i] = 0;
    }
    if(seq[18]!=nt)
    	tail[0] = 1;
    else
    	tail[1] = 0;

    if(seq[19]!=nt)
    	tail[0] = 1;
    else
    	tail[1] = 0;
    // then slide
    int i1 = 2; // end position of head
    int i2 = 19; // end position of tail
    bool found = 0;
    while( i2 < seq.length() ) {
    	int na = head[0] + head[1] + head[2] + tail[0] + tail[1];
    	i1++;
    	bool notN;
    	if (seq[i1] != nt) {
    		na++;
    		notN = 1;
    	}
    	else
    		notN = 0;
    	if (na <= 5) {
 	   		for(int i=i1+1; i<i1+15; i++) {
 //std::cout << "i=" << i << ", seq=" << seq.length() << std::endl;
    			if (seq[i] != nt)
    				na++;
    			if (na>5)
    				break;
    		}
    	}
    	if (na>5) { // fail, move to next
    		i2++;
    		if (i2 >= seq.length())
    			break;
    		head[0] = head[1];
    		head[1] = head[2];
    		if (notN)
    			head[2] = 1;
    		else
    			head[2] = 0;
    		tail[0] = tail[1];
    		if (seq[i2] != nt)
    			tail[1] = 1;
    		else
    			tail[1] = 0;

    	}
    	else {
    		found = 1;
    		break;
    	}
    }
    if (found)
    	return 1;
    else
    	return 0;

/*    pnt.resize(10, nt);
    if (seq.find(pnt) != std::string::npos) {

    }

    // then check if contain more than 75%
    // check head
    int sum = 0;
    for(int i=0; i<20; i++) {
        if (seq[i]==nt)
            sum++;
    }
    if (sum>=15)
        return 1;
    // check tail
    sum = 0;
    for(int i=0; i<20; i++) {
        if (seq[seq.length()-1-i]==nt)
            sum++;
    }
    if (sum>=15)
        return 1;
*/
    // not found
    return 0;
}

int GetTotalClipLength( string & cigar )
{
	int c1 = -1;
	int c2 = -1;
	int clip_len = 0;

	// begin clip
	for(int i=0; i<cigar.length(); i++) {
		if (cigar[i] == 'S') {
			c1 = i;
			break;
		}
		else if (cigar[i] == 'M')
			break;
	}
	if (c1 != -1) {
		string numb = cigar.substr( 0, c1 );
		if ( !std::all_of( numb.begin(), numb.end(), ::isdigit ) ) {
			cerr << "ERROR: [GetTotalClipLength] unidentified clip string " << cigar << endl;
			exit(10);
		}
		clip_len = stoi(numb);
	}

	// end clip
	if (cigar[cigar.length()-1] == 'S') {
		for(int i=cigar.length()-2; i>=0; i--) {
			if (!isdigit(cigar[i])) {
				c2 = i;
				break;
			}
		}
		int lc2 = cigar.length() - c2 - 2;
		if (lc2 <=0) {
			cerr << "ERROR: [GetTotalClipLength] unidentified clip string length " << cigar << endl;
			exit(10);
		}
		string numb = cigar.substr( c2+1, lc2 );
		if ( !std::all_of( numb.begin(), numb.end(), ::isdigit ) ) {
			cerr << "ERROR: [GetTotalClipLength] unidentified clip string " << cigar << endl;
			exit(10);
		}
		clip_len += stoi(numb);
	}

	return clip_len;
}

// end clip is negative
int GetMaxClipLen( SamRecord & rec )
{
	Cigar * myCigar = rec.getCigarInfo();
	int begin_clip = myCigar->getNumBeginClips();
	int end_clip = myCigar->getNumEndClips();
	if (begin_clip >= end_clip)
		return begin_clip;
	else
		return -end_clip;
}

string GetMaxClipSeq( SamRecord & rec )
{
	int clip_len = GetMaxClipLen( rec );
	if (clip_len==0) {
		cerr << "ERROR: [GetMaxClipSeq] no clip in sequence. Something is wrong!\n" << endl;
		exit(20);
	}

	string seq = rec.getSequence();
	string clip;
	if ( clip_len > 0 )
		clip = seq.substr(0, clip_len);
	else
		clip = seq.substr( seq.length() + clip_len );

	return clip;
}

string RevCompSeq( string & seq ) {
	string rv;
	int n = (int)seq.size();
	rv.resize( n );
	for( int i=n-1; i>=0; i-- )
		rv[n-i-1] = getCompNt( seq[i] );
	return rv;
}

char getCompNt( char nt)
{
	char c;
	switch( nt ) {
		case('A'):
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

