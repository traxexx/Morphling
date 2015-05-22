#include "SetReadMapGlobals.h"
#include "ReadMapGlobals.h"
#include "Globals.h"
#include <math.h>

bool IsRMGset()
{
	if ( AvrReadLen < 0 )
		return 0;
	if ( MinReadLen < 0 )
		return 0;
	if ( AvrInsSize < 0 )
		return 0;
	if ( MinDiscIns < 0)
		return 0;
	
// all filter pass
	return 1;
}

void SetRMGReadLength( int avr_read_len )
{
	AvrReadLen = avr_read_len;
	MinReadLen = avr_read_len / 2;
}

void SetRMGinsSize( int avr_ins_size )
{
	AvrInsSize = avr_ins_size;
	MinDiscIns = avr_ins_size * 3;
}

void SetRMGdepth( float depth )
{
	DEPTH = depth;
	float rc = depth * WIN / AvrReadLen / 2;
	MIN_READ_IN_WIN = round( rc / 4 );
	MAX_READ_IN_WIN = round( rc * 7 / 4 );
}