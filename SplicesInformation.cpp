#include "SplicesInformation.h"

struct _spliceInformation splicesInformation[ HASH_MAX ] ;

void SplicesInformation_Init() 
{
	memset( splicesInformation, -1, sizeof( splicesInformation ) ) ;
}

void SplicesInformation_Reset() 
{
	memset( splicesInformation, -1, sizeof( splicesInformation ) ) ;	
}

// val: non-negative integers. 0-deleted. others: the point merged to.
void SplicesInformation_Add( int pos, int strand, int val ) 
{
	int key = pos % HASH_MAX ;
	//printf( "Add: %d %d: %d\n", pos, strand, val ) ;
	while ( splicesInformation[ key ].pos != -1 && 
		splicesInformation[key].pos != pos && 
		splicesInformation[key].strand != strand )
	{
		++key ;
		if ( key >= HASH_MAX )
			key = 0 ;
	}
	//printf( "Add 2: %d %d %d\n", splicesInformation[key].pos, splicesInformation[key].strand, val ) ;
	splicesInformation[ key ].pos = pos ;
	splicesInformation[ key ].val = val ;
	splicesInformation[ key ].strand = strand ;
	//static int cnt = 0 ;
	//++cnt ;
	//printf( "###%d\n", cnt ) ;
}

// -1: Non-exist. 0-deleted. other postive integers: the point merged to.
int SplicesInformation_Get( int pos, int strand ) 
{
	int key = pos % HASH_MAX ;
	while ( splicesInformation[key].pos != -1  )
	{
		if ( splicesInformation[key].pos == pos && splicesInformation[key].strand == strand )
			return splicesInformation[key].val ;

		++key ;
		if ( key >= HASH_MAX )
			key = 0 ;
	}
	return -1 ;
}

// Almost the same as SpliceInformation_Get, except that it returns the position if that splice sites is not in the table.
int SplicesInformation_GetPos( int pos, int strand )
{
	int ret = SplicesInformation_Get( pos, strand ) ;
	//if ( ret != -1 )
	//	printf( "### %d %d\n", pos, ret ) ;
	if ( ret == -1 )
		return pos ;
	else
		return ret ;
}
