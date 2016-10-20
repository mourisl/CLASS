#ifndef _LSONG_READS_HEADER
#define _LSONG_READS_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "sam.h"

#define MAX_READ 1000000 // The # of reads for a region
#define MAX_READ_SPLICE 4
//#define MAX_READ_GENE 2000000 // The # of reads for genes

#define ID_LENGTH 100 

#define GENE_READ_MULT 10
#define GENE_READ_HASH_MULT 1 
//#define READS_LENGTH 75

struct _readTree
{
	int index ;
	int timeStamp ;
	char id[ID_LENGTH] ; 
	struct _readTree *left, *right ;
} ;

struct _readFile
{
	FILE *fp ;
	samfile_t *fpsam ;
	bool sam ; // Use sam file API or just FILE pointer.
} ;

struct _read
{
	int start, end ;
	short strand ;
	int splices[MAX_READ_SPLICE][2] ; // The coordinate of the splice junction
	short scnt ; // splice cnt 
	int mateInd ; // -1 if its mate is out of the region	
	bool unique ;
} ;


struct _geneRead
{
	struct _read *reads[GENE_READ_MULT] ;
	struct _readTree *hash[GENE_READ_MULT] ;
} ;

struct _readFile OpenReadFile( char *prefix ) ;

//int GetReadsLength( FILE *fp ) ;
void GetReadsInfo( struct _readFile file, int &readsLen, int &fragLen, int &fragStd, long long &totalReadCnt ) ;

/**
  Extract reads covered the region from the SAM file.
  extent stores the start of the left-most reads and the end of the right-most reads
  @return: The number of reads. And it puts the reads in the array passed from the main program
*/
//int ExtractReads( FILE *fp, char *rchrom, int rstart, int rend, struct _read reads[], int extent[2] ) ; 
int ExtractReads( struct _readFile file, char *rchrom, int rstart, int rend, struct _read reads[], int extent[2] ) ; 

void InitGeneReads( struct _geneRead *geneReads ) ;

//int ExtractGeneReads( FILE *fp, char *rchrom, int rstart, int rend, struct _geneRead geneReads, int extent[2] ) ;
int ExtractGeneReads( struct _readFile file, char *rchrom, int rstart, int rend, struct _geneRead geneReads, int extent[2] ) ;

void ClearReadHash( struct _geneRead &geneReads ) ;
#endif
