/**
Read in the transcripts or information from evidence gtf file
*/
#ifndef _LSONG_EVIDENCE_HEADER
#define _LSONG_EVIDENCE_HEADER

#include "SolveRegion.h"
#include <stdlib.h>

#define EVIDENCE_MARGIN 10

struct _evidence
{
	struct _exon *exons ;
	int ecnt ;
} ;

struct _evidence_hash
{
	int p[4] ;
	int val ;
} ;

void Evidence_Init( char *evidenceFile, struct _evidence **evidence ) ;

// Extract evidence from start to end.
int Evidence_Extract( char *chrom, struct _evidence *e ) ;

// The new evidences from the transcripts
void Evidence_Add( int intron1[2], int intron2[2], int val ) ;
// return 0, if does not exist.
int Evidence_Val( int intron1[2], int intron2[2] ) ;

bool WithinEvidenceMargin( int a, int b ) ;

// Reset the hash table
void Evidence_ClearHash() ;

// Processing the exons of evidences
int ExtractEvidenceExons( struct _exon **evidenceExons, struct _evidence *evidences, int eviCnt ) ;
int IsPosInEvidenceExon( int pos, struct _exon *evidenceExons, int eviExonCnt ) ;
#endif
