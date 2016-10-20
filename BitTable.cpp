#include "BitTable.h"

BitTable::BitTable()
{
	//tab = new INT64[1] ;
	//size = UNIT_SIZE ;
	tab = NULL ;
	size = 0 ;
}

BitTable::BitTable( int s )
{
	if ( s == 0 )
		s = UNIT_SIZE ;
		
	if ( s & UNIT_MASK )
	{
		asize = s / UNIT_SIZE + 1 ;
		tab = new INT64[ asize ]() ;
	}
	else
	{
		asize = s / UNIT_SIZE ;
		tab = new INT64[ asize ]() ;
	}

	size = s ;
	Reset() ;
}

BitTable::~BitTable()
{
	//printf( "hi %d\n", size ) ;
	if ( tab != NULL )
		delete [] tab ;
	//size = -1 ;
}

void BitTable::Init( int s )
{
	if ( tab != NULL )
		delete [] tab ;

	if ( s == 0 )
		s = UNIT_SIZE ;
		
	if ( s & UNIT_MASK )
	{
		asize = s / UNIT_SIZE + 1 ;
		tab = new INT64[ asize ]() ;
	}
	else
	{
		asize = s / UNIT_SIZE ;
		tab = new INT64[ asize ]() ;
	}

	size = s ;
	Reset() ;
}

void BitTable::Set( int i )
{
	int ind, offset ;
	ind = i / UNIT_SIZE ;
	offset = i & UNIT_MASK ;
//printf( "%d,%d %d,%d: %llx", 311 / UNIT_SIZE, 311 & UNIT_MASK, 279 / UNIT_SIZE, 279 & UNIT_MASK, tab[ind] ) ; 	
	tab[ind] |= ( (INT64)1 << offset ) ;
//printf( " %llx %d\n", tab[ind], UNIT_SIZE ) ;
}

void BitTable::Unset( int i )
{
	int ind, offset ;
	ind = i / UNIT_SIZE ;
	offset = i & UNIT_MASK ;
	
	tab[ind] &= ( (INT64)-1 ^ ( (INT64)1 << offset ) ) ;
}

void BitTable::Flip( int i )
{
	int ind, offset ;
	ind = i / UNIT_SIZE ;
	offset = i & UNIT_MASK ; 	
	
	tab[ind] ^= ( (INT64)1 << offset ) ;
}

bool BitTable::Test( int i ) 
{
	int ind, offset ;
	ind = i / UNIT_SIZE ;
	offset = i & UNIT_MASK ; 
	
	if ( tab[ind] & ( (INT64)1 << offset ) )
		return true ;
	else
		return false ;
}

void BitTable::Not()
{
	int i ;
	for ( i = 0 ; i < asize ; ++i )
	{
		tab[i] = ~tab[i] ;
	}
}

void BitTable::And( const BitTable &in )
{
	if ( asize != in.asize )
		return ;
	int i ;
	for ( i = 0 ; i < asize ; ++i )
		tab[i] &= in.tab[i] ; 
}

void BitTable::Or( const BitTable &in )
{
	if ( asize != in.asize )
		return ;
	int i ;
	for ( i = 0 ; i < asize ; ++i )
		tab[i] |= in.tab[i] ; 
}

int BitTable::Count()
{
	if ( size <= 0 )
		return 0 ;
	INT64 k ;
	int i, j, ret = 0 ;
	for ( i = 0 ; i < asize - 1 ; ++i )
	{
		k = tab[i] ;
		while ( k )
		//for ( j = 0 ; j < UNIT_SIZE ; ++j ) 	
		{
			if ( k & 1 )
				++ret ;
			//printf( "### %d %d %d %d\n", ret, asize, size, k ) ;
			k /= 2 ;
		}
	}
	//printf( "(%d) ", ret ) ;	
	k = tab[ asize - 1 ] ;
	for ( i = 0 ; i < size & UNIT_MASK ; ++i )
	{
		if ( k & 1 )
			++ret ;
		k /= 2 ;
	}
	//for ( i = 0 ; i < asize ; ++i ) 
	//	printf( "%llu ", tab[i] ) ;
	//printf( "(%d)\n", ret ) ;
	return ret ;
}

void BitTable::Reset()
{
	for ( int i = 0 ; i < asize ; ++i )
		tab[i] = 0 ;
}

bool BitTable::IsEqual( const BitTable &in )
{
	int i, k ;
	if ( in.size != size )
		return false ;
	//printf( "== %d %d\n", in.size, size ) ;
	//if ( size & UNIT_MASK )
	//	k = size / UNIT_SIZE + 1 ;
	//else
	//	k = size / UNIT_SIZE ;	
		
	for ( i = 0 ; i < asize ; ++i )	
		if ( tab[i] != in.tab[i] )
			return false ;
	return true ;
}
