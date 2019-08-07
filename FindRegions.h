/**
  Scan the depth(coverage) file and the splice file to determine an region containing exons.

  Li Song
  Aug 16, 2012
*/
#ifndef _LSONG_FIND_REGIONS_HEADER
#define _LSONG_FIND_REGIONS_HEADER

#include <string.h>
#include "SolveRegion.h"

//#define CUT_THRESHOLD 1
//#define CUT_MERGE 15

//int ExtractExons( FILE *fp ) ;

void FindRegions_Init( char *prefix, char *otherJunctionFile ) ;
void FindRegions_Destroy() ;

int FindRegions_Next( char *chrom, struct _splice splices[], int &scnt, int &start, int &end, struct _evidence *evidences, int &eviCnt ) ; 
#endif
