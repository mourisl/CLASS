#ifndef _LSONG_SOLVE_REGION_HEADER
#define _LSONG_SOLVE_REGION_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
//#include "sam.h"
#include "lp_lib.h"
#include "Reads.h"
#include "SplicesInformation.h"
#include "Evidence.h"
#include "TTS.h"

#define MAX_LENGTH 2000000
#define MAX_POINT 1000 
#define ALPHA 0.01
#define SPLICE_ALPHA 2 
#define MERGE_DISTANCE 6 // When two points is less than distance, merge one of them.
#define SLACK_PENALTY 20 
#define MAX_NEIGHBOR 20  // The maximal # of neighbor exons.
#define MAX_EXON 2000000
#define LOG_BASE 1.25 
#define LONG_EXON 356 

#define INF 2000000000

struct _point
{ 
	int pos ;
	int strand ; // 0-negative, 1-positive, -1-not sure.
	int support ; // The # of reads that supports this point.
	int merge ; // Which point does this point merged to for each strand
	int strong ; // -1, not a splice point. 0-soft splice point, 1-strong splice point
} ;
 

struct _splice
{
	int pos ;
	int type ; // 0-start(end of an exon), 1-end(start of an exon)
	int otherInd ; // The index of the other end. -1, if the other end is out of the region.  
	int otherPos ; // The coordinates of the other end
	int support ;
	int strand ;
	int uniqSupport, secSupport ;
	int uniqEditDistance, secEditDistance ;
	//int soft ; // The soft boundary position along with the splice site. Like [...)....(...] the two splice sites have the soft boundary. -1: If N/A.
} ;

struct _spliceIndex
{
	int start, end ; // A splice junction begining with "end of an exon" and terminating with "start of an exon". -1 means out of the region. 
	int strand ;
	int startPos, endPos ;
	int support ; // The # of reads support this splice.
} ;

struct _exon
{
	int start, end ;
	int strand ; // -1: unknown. 0:-, 1:+.

	bool hasSoft ; // Whether it has a soft end	
	int *prev, *next ; // The coordinates of the previous exon's end and the next exon's start.
	int pcnt, ncnt ; // How many before and next used. If 0, it must be a soft boundary.		
	
	int psize, nsize ;	
} ;

// Just pair.
struct _pair 
{
	int a, b ;
} ;

/*int SolveRegion( char *chrom,  
		struct _point inputStart[], struct _point inputEnd[], int inputStartCnt, int inputEndCnt, 
		struct _spliceIndex spliceIndices[], int siCnt, 
		int *depth, 
		struct _read reads[], int readCnt,
	        struct _exon allExons[], int &exonCnt,
		bool force ) ; */

// Convert the splices from FindRegions to the arguments for SolveRegion
int SolveRegion_Wrapper( char *chrom, struct _splice splices[], int scnt, int softStart, int softEnd,
		int *depth, struct _read reads[], int readCnt, struct _evidence *evidences, int eviCnt, struct _exon allExons[], int &exonCnt,
		int minSpliceSupportPlus, int minSpliceSupportMinus, bool force ) ;
#endif
