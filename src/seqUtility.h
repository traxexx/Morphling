#ifndef SEQUTILITY_H
#define SEQUTILITY_H

#include "SamRecord.h"
#include <string>

using std::string;

// check if the read contains polyA/T
bool ContainRepetative( string & seq, char nt );

int GetTotalClipLength( string & cigar );

int GetMaxClipLen( SamRecord & rec );

string GetMaxClipSeq( SamRecord & rec );

string RevCompSeq( string & seq );

char getCompNt( char nt);

// one end of the pair is unmap
bool OneEndUnmap( int flag );

#endif