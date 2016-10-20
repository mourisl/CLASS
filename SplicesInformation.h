/**
   Use a hash table to hold special splice junctions, like those removed and merged.

   Sept 8, 2012 - Li Song
*/

#ifndef _SPLICES_INFORMATION_HEADER
#define _SPLICES_INFORMATION_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HASH_MAX 400027

struct _spliceInformation
{
	int pos ;
	int strand ;
	int val ;
} ;

void SplicesInformation_Init() ;
void SplicesInformation_Reset() ;
// val: non-negative integers. 0-deleted. others: the point merged to.
void SplicesInformation_Add( int pos, int strand, int val ) ;
// -1: Non-exist. 0-deleted. other postive integers: the point merged to.
int SplicesInformation_Get( int pos, int strand ) ;
// Almost the same as SplicesInformation_Get, except that it returns the position if that splice sites is not in the table.
int SplicesInformation_GetPos( int pos, int strand ) ;
#endif
