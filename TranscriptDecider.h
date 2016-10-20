#ifndef _LSONG_TRANSCRIPT_DECIDER_HEADER
#define _LSONG_TRANSCRIPT_DECIDER_HEADER

#include <stdio.h>
#include "SolveRegion.h"
#include "BitTable.h"
#include "Evidence.h"
#include <pthread.h>
#include <unistd.h>

void TranscriptDecider_Init( char *prefix ) ;
void TranscriptDecider_Go( char *chrom, struct _exon exons[], int exonCnt, struct _evidence *inEvidences, int inEviCnt ) ;

#endif
