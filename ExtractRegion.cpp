/**
Extract the interested regions from the header file based on the information provided by *.as file.
Build a new header file containing those regions.
It stores the new right combinations in temporary file "interest_right_comb.out".


Usage: ./a.out event organism *.as header right_comb

Li Song
July 27, 2012
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX 6000000
#define MAX_POINT 1000
struct _event
{
	int start, end ;
	char chrom[6] ;
} ;
struct _event events[MAX] ;

struct _point
{
	int pos, strand ;
	int support ;
} ;

int eventCnt = 0 ;

char line[1024], line2[1024] ;
struct _point start[MAX_POINT], end[MAX_POINT] ;

struct _chrom
{
	int startInd ;
	char name[6] ;
} ;
struct _chrom chroms[100] ;
int chromCnt = 0 ;

struct _spliceIndex
{
	int start, end ;
	int strand ;
} ;
struct _spliceIndex spliceIndices[MAX_POINT] ;

int Comp( const void *p1, const void *p2 )
{
	struct _event a, b ;
	a = *(struct _event *)p1 ;
	b = *(struct _event *)p2 ;
	
	if ( a.start != b.start )
		return a.start - b.start ;
	else
		return a.end - b.end ;
}

int GetChromStartInd( char *s )
{
	int i ;
	for ( i = 0 ; i < chromCnt ; ++i )
		if ( !strcmp( chroms[i].name, s ) )
			return chroms[i].startInd ;
	return -1 ;
}

int main( int argc, char *argv[] )
{
	if ( argc != 6 )
	{
		printf( "Usage: ./a.out event organism *.as header right_comb\n" ) ;
		return 0 ;
	}

	int i, j, k ;   
	FILE *fpHeader, *fpAs, *fpRightComb, *fpNewRightComb ;
	bool pair ; // Read two lines?
	char id[32], chrom[10], type[50], fid[50] ;
	char prevChrom[10] = "-1" ;
	int tag = -1 ;
	int offset, startCnt, endCnt ;
	int rstart, rend ;
	int siCnt ;

	sprintf( line, "%s_right_comb.out", argv[1] ) ;

	fpAs = fopen( argv[3], "r" ) ;
	fpHeader = fopen( argv[4], "r" ) ;
	fpRightComb = fopen( argv[5], "r" ) ;
	fpNewRightComb = fopen( line, "w" ) ;

	pair = false ;
	if ( strstr( argv[1], "IR" ) || strstr( argv[1], "AE" ) )
		pair = true ;

	// Read in the events    
	while ( fgets( line, sizeof( line ), fpAs ) )
	{
		if ( strstr( line, "_ON" ) || strstr( line, "AE" ) )
		{
			fgets( line2, sizeof( line2 ), fpAs ) ;
			if ( !strstr( line2, argv[2] ) )
				continue ;
		}

		if ( !strstr( line, argv[1] ) || !strstr( line, argv[2] ) )
			continue ;

		sscanf( line, "%s %s %s %s %d %d", id, type, fid, events[ eventCnt ].chrom, &events[ eventCnt ].start, &events[ eventCnt ].end ) ;
		
		if ( strcmp( prevChrom, events[ eventCnt ].chrom ) )
		{
			chroms[ chromCnt ].startInd = eventCnt ;
			strcpy( chroms[ chromCnt ].name, events[ eventCnt ].chrom ) ;
			++chromCnt ;
		}

		if ( tag == -1 )
		{
			strcpy( prevChrom, events[ eventCnt ].chrom ) ;
			tag = 0 ;
		}
		else if ( strcmp( prevChrom, events[ eventCnt ].chrom ) )
		{
			qsort( &events[tag], eventCnt - tag, sizeof( events[0] ), Comp ) ;

			strcpy( prevChrom, events[ eventCnt ].chrom ) ;
			tag = eventCnt ;
		}

		++eventCnt ;
	}
       //for ( i = 0 ; i < eventCnt ; ++i )
	//	printf( "%s %d %d\n", events[i].chrom, events[i].start, events[i].end ) ; 
	// Go through all the regions and output those containing the alternative splicing events.
	k = -1 ; // The tag for the events array
	while ( fscanf( fpHeader, "%s %d %d %d", chrom, &offset, &startCnt, &endCnt ) != EOF )
	{

		//sscanf( line, "%s %d %d %d", chrom, &offset, &startCnt, &endCnt ) ;
		rstart = offset ;
		rend = -1 ;
		for ( i = 0 ; i < startCnt ; ++i )
		{
			fscanf( fpHeader, "%d %d %d", &start[i].pos, &start[i].strand, &start[i].support ) ;
		}

		for ( i = 0 ; i < endCnt ; ++i )
		{
			fscanf( fpHeader, "%d %d %d", &end[i].pos, &end[i].strand, &end[i].support ) ;
			if ( end[i].pos > rend )
				rend = end[i].pos ;
		}	
		
		fscanf( fpHeader, "%d", &siCnt ) ;
		for ( i = 0 ; i < siCnt ; ++i )
			fscanf( fpHeader, "%d %d %d", &spliceIndices[i].start, &spliceIndices[i].end, &spliceIndices[i].strand ) ;
		
		if ( k == -1 || strcmp( events[k].chrom, chrom ) )
			k = GetChromStartInd( chrom ) ;
	
		while ( k < eventCnt && events[k].end < rstart && !strcmp( events[k].chrom, chrom ) )
			++k ;
		
		if ( k < eventCnt && !strcmp( events[k].chrom, chrom ) && events[k].start <= rend )
		{
			printf( "%s %d %d %d\n", chrom, offset, startCnt, endCnt  ) ;
			for ( i = 0 ; i < startCnt ; ++i )
			{
				printf( "%d %d %d\n", start[i].pos, start[i].strand, start[i].support ) ;
			}
			for ( i = 0 ; i < endCnt ; ++i )
			{
				printf( "%d %d %d\n", end[i].pos, end[i].strand, end[i].support ) ;
			}
			printf( "%d\n", siCnt ) ;	
			for ( i = 0 ; i < siCnt ; ++i )
				printf( "%d %d %d\n", spliceIndices[i].start, spliceIndices[i].end, spliceIndices[i].strand ) ;
			
			fgets( line2, sizeof( line2 ), fpRightComb ) ;
			j = atoi( line2 ) ;
			fprintf( fpNewRightComb, "%s", line2 ) ;	
			for ( i = 0 ; i < j ; ++i )
			{
				fgets( line2, sizeof( line2 ), fpRightComb ) ;
				fprintf( fpNewRightComb, "%s", line2 ) ;
			}
		}
		else if ( strcmp( events[k].chrom, chrom ) )
			--k ;
	}	

	fclose( fpHeader ) ;
	fclose( fpRightComb ) ;
	fclose( fpNewRightComb ) ;
	return 0 ;
}
