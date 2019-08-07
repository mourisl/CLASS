#include "TranscriptDecider.h"

#define MAX_TRANSCRIPT 1000003
#define BEST_COVER 1
#define ALL_TRANSCRIPT 0 
#define EVIDENCE_MARGIN 10
#define CUTOFF_RATIO 0 
#define USE_F3 150 
#define USE_DP 200000
#define DP_EXON_ID_STRIDE 6
#define DP_STRIDE_CNT 5
#define MAX_TRANSCRIPT_CONSTRAINT (MAX_READ+1000)

#define SIZEOF_STRIDE ( ( DP_EXON_ID_STRIDE + 1 )*( DP_EXON_ID_STRIDE + 1 )*( DP_EXON_ID_STRIDE + 1 )*( DP_EXON_ID_STRIDE + 1 )*( DP_EXON_ID_STRIDE + 1 ) )


extern int FRAG_LENGTH, FRAG_STD ;

extern int READS_LENGTH ;

extern long long TOTAL_READ_COUNT ;

extern bool VERBOSE ;

extern int NUM_OF_THREADS ;
extern pthread_attr_t pthreadAttr ;

bool USE_SET_COVER ;
double FPKM_FRACTION = 0.05 ;

extern char *gtfIdLabel ;

struct _transcript
{
	int *eid ; // Exon id
	int ecnt ; // exon cnt
	int strand ;
	double abundance ; // The heuristic number of fragments belonging to the transcript
	int geneId, transcriptId ;
	double FPKM ;
} ;

struct _transcriptConstraint
{
	int type ; // 0-consecutive exons, 1-splice junction, 2-evidence
	int info[2] ; // type 0: the exon indices. info[1], the end position.
	
	int leftInfo[5][2], rightInfo[5][2] ; 	
	
	int infoCnt, lcnt, rcnt ;
	int strand ;
	int support ; // The number of reads supporting it.
	int uniqSupport ; // The number of unique reads supporting it
	double abundance ;

	double normAbund ; // normalized abundance 
	//int effectiveLength ; // normAbund * effectiveLength = support
	double effectiveCount ; // How many count when updating the count for a transcript.
} ;

// The exon node in the splice graph
struct _exonNode
{
	int id ; // The index in the exon array.
	int *next, *prev ;
	int ncnt, pcnt ;
	int nsize, psize ; // The # of memory allocated for next and prev
	int farthest ; // The farthest end position of substrancripts starting from this one
	int before ; // The exons before
	bool *used ;
	int subcnt ;
	int weight ; // some dummy weight
	int supportTranscript ; // How many transcripts contains this exon
	int geneId ;
} ;

// The argument passed to the pthread calling
struct _pthreadArgPickTranscripts
{
	char *chrom ;
	struct _exon *exons ;
	int exonCnt ;
	struct _exonNode *nodes ;
	struct _transcriptConstraint *tc ;
	int tcCnt ;
	int from, to ;

	int evidenceIndFrom, evidenceIndTo ;
} ;

int currPickTxptThreadsCnt ; 
pthread_mutex_t pickTxptMutex ;
pthread_mutex_t outputTxptMutex ;
pthread_cond_t idleCond ;
pthread_cond_t clearCond ;

FILE *td_fpDepth ;
struct _readFile td_fpReads ;
int transcriptStart = -1, transcriptEnd = -1, transcriptLastExonStart = -1 ;

//struct _transcript transcripts[MAX_TRANSCRIPT], alltranscripts[MAX_TRANSCRIPT] ;
struct _transcript outputTranscripts[ MAX_TRANSCRIPT ] ;
int outputTxptCnt ;

//BitTable btable[MAX_TRANSCRIPT] ;
double maxAbundance[MAX_EXON] ;

//int tcnt, atcnt ;

//struct _read transcriptReads[MAX_READ_GENE] ;
struct _geneRead transcriptReads ;

int trCnt ;

// abundanceTranscriptConstraints holds the constraints from single-end constraints.
// transcriptConstraints combines the abundanceTranscriptConstraints if the data is paired-end
struct _transcriptConstraint transcriptConstraints[MAX_TRANSCRIPT_CONSTRAINT], 
		alltc[MAX_TRANSCRIPT_CONSTRAINT], abundanceConstraints[MAX_TRANSCRIPT_CONSTRAINT] ;

//int pairExonTc[MAX_TRANSCRIPT_CONSTRAINT] ;
//int tcUsedCnt[MAX_TRANSCRIPT_CONSTRAINT] ; // How many times a transcript constraint used.

//int eidUsed[MAX_TRANSCRIPT] ;
int euCnt = 0 ;
//bool eidUsedFlag[MAX_EXON] ;

int geneId[MAX_EXON] ; // This is actually the set.
int geneTranscriptId[MAX_EXON] ; // The transcript id used for each gene
int geneIdCnt ;
double geneAbundance[MAX_EXON] ;
//int bufferEid[MAX_EXON], toEndDist[MAX_EXON] ;

struct _evidence *evidences ;
int eviTag,	// The beginning of the evidences  
	eviCnt ; // The count of the evidences

struct _dp
{
	int *eid ;
	int ecnt ;
	double cover ;

	int timeStamp ;
	double minAbundance ;
} ;

struct _dpStride
{
   	// To memorize the consecutive exons up to 6 (5 strides) by looking at the strides
   	struct _dp p[SIZEOF_STRIDE] ;
} ;

struct _dpHash
{
	int *eid ;
	int ecnt ;
	double cover ;
	int cnt ; // The length of the header for eid.
	double minAbundance ;
	int timeStamp ;
} ;

// A list of transcriptConstraints Id
struct _tcList 
{
	int *id ;
	int cnt ;
} ;

struct _dpAttribute
{
	struct _dp *f1, **f2, ***f3 ; 
	struct _dpStride *stride ;
	struct _dpHash *hash ;
	struct _dpHash *longHash ; // The hash for long subtranscript 
	int offset ;
	int size ;
	double minAbundance ;
	bool forAbundance ; // Is the DP for calculating most expressed transcript?
	
	bool solveRemaining ; // The main part using evidences finished. 

	struct _tcList *pairExonTc ; // Store the tc starting with a exon and compatible with it and its next. 
	int *pairExonTcIds ; 
	int *bufferEid ;
	//struct _tcList *partPairExonTc ; // The 
	int timeStamp ;
	int *solveCaller ; // count how many times we called SolveSub without using memoir.
} ;

struct _dp SolveSubTranscript( int visit[], int vcnt, const struct _dpAttribute &attr, struct _exonNode nodes[], struct _exon exons[], 
	struct _transcriptConstraint *tc, int tcCnt ) ;

// For set cover
struct _setcover
{
	int choose[1024] ;
	int cnt ;
	double weight ;
} ;

void TranscriptDecider_Init( char *prefix )
{
	char buffer[1024] ;
	
	sprintf( buffer, "%s.sam", prefix ) ;
	td_fpReads = OpenReadFile( prefix ) ; 
	
	//printf( "file pointer: %lld\n", td_fpReads.fp ) ;
	sprintf( buffer, "%s.depth", prefix ) ;
	td_fpDepth = fopen( buffer, "r" ) ;
	//atcnt = 0 ;
	transcriptEnd = -1 ;
	transcriptLastExonStart = -1 ;
	InitGeneReads( &transcriptReads ) ;
	//for ( int i = 0 ; i < MAX_TRANSCRIPT ; ++i )
	//	alltranscripts[i].eid = NULL ;	

	/*if ( evidence != NULL )
		Evidence_Init( evidence, &evidences ) ;
	else
		evidences = NULL ;*/
	memset( geneTranscriptId, 0, sizeof( geneTranscriptId ) ) ;
}

void TranscriptDecider_Destroy()
{
	CloseReadFile( td_fpReads ) ;
	ReleaseGeneReads( &transcriptReads ) ;
	fclose( td_fpDepth ) ;
}

// @return: FPKM of the transcripts
double GetTranscriptAbundance( struct _exon exons[], struct _transcript *transcript )
{
	int i, len = 0 ;
	double factor = ( FRAG_LENGTH == READS_LENGTH ) ? 1000000 : 2000000 ;
	double ret ; 
	for ( i = 0 ; i < transcript->ecnt ; ++i )
		len += exons[ transcript->eid[i] ].end - exons[ transcript->eid[i] ].start + 1 ;
	//printf( "%lf (%d)\n", transcript->abundance, len ) ;
	ret = transcript->abundance / ( len / (double)1000.0 ) / ( TOTAL_READ_COUNT / factor ) ;
	return ret ;
}

void OutputTranscript( char *chrom, struct _exon exons[], struct _transcript *transcript )
{
	//return ;
		//printf( "hi2 %d\n", transcript->eid ) ;
	int i ;
	int a = transcript->eid[0] ;
	int b = transcript->eid[ transcript->ecnt - 1 ] ;
	char strand[2] = "+" ;
	if ( exons[a].strand == 0 )
		strand[0] = '-' ;
	else if ( exons[a].strand == -1 )
		strand[0] = '.' ;
	
	//transcript->geneId = 0 ;
	//transcript->transcriptId = 0 ;
	//transcript->FPKM = 0 ;

	char prefix[512] ;
	prefix[0] = '\0' ;
	if ( gtfIdLabel != NULL )
	{
		sprintf( prefix, "%s_", gtfIdLabel ) ;
	}
	printf( "%s\tCLASS\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; Abundance \"%.6lf\";\n", 
			chrom, exons[a].start, exons[b].end, strand, 
			prefix, chrom, transcript->geneId, 
			prefix, chrom, transcript->geneId, transcript->transcriptId, transcript->FPKM ) ;
	for ( i = 0 ; i < transcript->ecnt ; ++i )
	{
		printf( "%s\tCLASS\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; "
				"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; Abundance \"%.6lf\"\n",
				chrom, exons[ transcript->eid[i] ].start, exons[ transcript->eid[i] ].end, strand,
				prefix, chrom, transcript->geneId, 
				prefix, chrom, transcript->geneId, transcript->transcriptId, 
				i + 1, transcript->FPKM ) ;
	}
	//++transcriptId ;
}

int CompInt( const void *p1, const void *p2 )
{
	return ( *( int * )p1 ) - ( *( int * )p2 ) ;
}

int CompDouble( const void *p1, const void *p2 )
{
	double a = *(double *)p1 ;
	double b = *(double *)p2 ;
	if ( a < b )
		return -1 ;
	else if ( a > b )
		return 1 ;
	else
		return 0 ;
}

// Returnt the nearest exons ending before dest
int EvidenceBinarySearch( int dest, struct _exon *e, int ecnt )
{
	int l, r, m ;
	l = 0 ;
	r = ecnt -1 ;
	while ( l <= r )
	{
		m = ( l + r ) / 2 ;
		if ( e[m].end == dest )
			return m ;
		else if ( e[m].end < dest )
			l = m + 1 ;
		else
			r = m - 1 ;
	}
	return l - 1 ;
}

bool IsInEvidence( int tpair[][2], int strand )
{
	int j, l ;
	int epair[2][2] ;
	j = Evidence_Val( tpair[0], tpair[1] ) ;

	if ( j == 1 )
		return true ;
	else if ( j == 2 )
		return false ;

	for ( j = eviTag ; j < eviCnt && evidences[j].exons[0].start <= tpair[0][0] + EVIDENCE_MARGIN ; ++j )
	{
		if ( evidences[j].exons[0].strand != strand )
			continue ;
		l = EvidenceBinarySearch( tpair[0][0] - EVIDENCE_MARGIN, evidences[j].exons, evidences[j].ecnt ) ;
		if ( l < 0 )
			l = 0 ;
		for ( ; l < evidences[j].ecnt - 2 && 
				evidences[j].exons[l].end <= tpair[0][0] + EVIDENCE_MARGIN ; ++l )
		{
			epair[0][0] = evidences[j].exons[l].end ;
			epair[0][1] = evidences[j].exons[l + 1].start ;
			epair[1][0] = evidences[j].exons[l + 1].end ;
			epair[1][1] = evidences[j].exons[l + 2].start ;

			if ( WithinEvidenceMargin( tpair[0][0], epair[0][0] ) && WithinEvidenceMargin( tpair[0][1], epair[0][1] ) 
					&& WithinEvidenceMargin( tpair[1][0], epair[1][0] ) && WithinEvidenceMargin( tpair[1][1], epair[1][1] ) )
			{
				//printf( "- (%d, %d, %d) %d %d\n", i, atcnt, k, tpair[0][0], tmpCnt ) ;
				Evidence_Add( tpair[0], tpair[1], 1 ) ;
				return true ;
			}
		}
	}
	Evidence_Add( tpair[0], tpair[1], 2 ) ;
	return false ;
}

// Search the first constraints starts after dest.
int BinarySearchTranscriptConstraints( int dest, struct _transcriptConstraint *tc, int tcCnt )
{
	//return 0 ;
	int l, r, m ;
	l = 0 ;
	r = tcCnt - 1 ;
	while ( l <= r )
	{
		m = ( l + r ) / 2 ;
		if ( tc[m].leftInfo[0][0] >= dest )
			r = m - 1 ;
		else
			l = m + 1 ;
	}
	/*if ( l > 0 )
		return l - 1 ;
	else
		return 0 ;*/
	return l ;
}

// Search the last exon index starting before dest
int BinarySearchExons( int dest, int eid[], int ecnt, struct _exon exons[] ) 
{
	//return 0 ;
	int l, r, m ;
	l = 0 ;
	r = ecnt - 1 ;
	while ( l <= r )
	{
		m = ( l + r ) / 2 ;
		
		if ( exons[ eid[m] ].start > dest )
			r = m - 1 ;
		else
			l = m + 1 ;
	}

	if ( r >= 0 ) 
		return r ;
	else
		return 0 ;
}

int TranscriptConstraintDifference( struct _transcriptConstraint *a, struct _transcriptConstraint *b )
{
	/*if ( ( a->type == 0 || b->type == 0 ) && a->type != b->type )
		return a->type - b->type ;
	if ( a->leftInfo[0][0] != b->leftInfo[0][0] && a->type != 0 && b->type != 0 )
		return a->leftInfo[0][0] - b->leftInfo[0][0] ;*/
	if ( a->leftInfo[0][0] != b->leftInfo[0][0] )
		return a->leftInfo[0][0] - b->leftInfo[0][0] ;

	if ( a->type != b->type )
		return a->type - b->type ;
	else if ( a->strand != b->strand )
		return a->strand - b->strand ;
	else
	{
		if ( a->type == 0 )
			return a->info[0] - b->info[0] ;
		else if ( a->type == 1 )
		{
			int i ;
			if ( a->lcnt != b->lcnt )
				return a->lcnt - b->lcnt ;
			if ( a->rcnt != b->rcnt )
				return a->rcnt - b->rcnt ;
			for ( i = 0 ; i < a->lcnt ; ++i )
			{
				if ( a->leftInfo[i][0] != b->leftInfo[i][0] )
					return a->leftInfo[i][0] - b->leftInfo[i][0] ;
				if ( a->leftInfo[i][1] != b->leftInfo[i][1] )
					return a->leftInfo[i][1] - b->leftInfo[i][1] ;
			}
			for ( i = 0 ; i < a->rcnt ; ++i )
			{
				if ( a->rightInfo[i][0] != b->rightInfo[i][0] )
					return a->rightInfo[i][0] - b->rightInfo[i][0] ;
				if ( a->rightInfo[i][1] != b->rightInfo[i][1] )
					return a->rightInfo[i][1] - b->rightInfo[i][1] ;
			}
		}
		else if ( a->type == 2 )
		{
			int i ;
			for ( i = 0 ; i < 2 ; ++i )
			{
				if ( a->leftInfo[i][0] != b->leftInfo[i][0] )
					return a->leftInfo[i][0] - b->leftInfo[i][0] ;
				if ( a->leftInfo[i][1] != b->leftInfo[i][1] )
					return a->leftInfo[i][1] - b->leftInfo[i][1] ;
			}
		}
	}
	if ( a->effectiveCount > b->effectiveCount )
		return 1 ;
	else if ( a->effectiveCount < b->effectiveCount )
		return -1 ;

	return 0 ;

}

int CompTranscriptConstraint( const void *p1, const void *p2 )
{
	struct _transcriptConstraint *a, *b ;
	a = ( ( struct _transcriptConstraint * )p1 ) ;
	b = ( ( struct _transcriptConstraint * )p2 ) ;
	return TranscriptConstraintDifference( a, b ) ;
}

int GetFather( int tag, int *father )
{
	if ( father[tag] != tag )
		return father[tag] = GetFather( father[tag], father ) ; 
	return tag ;
}

/**
  Use the LP to quantificate transcripts.
*/
void Quantification( char *chrom, int *points, int pcnt, struct _exon exons[],
	struct _transcript *transcripts, int tcnt )
{
	return ; // not used anymore
	int i, j, k ;
	int *depth = (int *)malloc( ( pcnt + 10 ) * sizeof( *depth ) ) ;
	int tmp, psum = 0, sum = 0 ;
	char buffer[100] ;
	int position, d, tag = -1 ;
	bool started = false ;

	// Read in the depth file and build the depth for each region.
	while ( 1 )
	{
		fscanf( td_fpDepth, "%s %d %d", buffer, &position, &d ) ;
		tmp = strcmp( buffer, chrom ) ;
		if ( ( !tmp && position > points[pcnt - 1] ) || ( tmp && started ) )
			break ;
		if ( ( tmp && !started ) || ( !tmp && position < points[0] ) )
			continue ;
		started = true ;
		
		psum += d ;
		sum += d ;
		if ( position == points[tag + 1] ) 
		{
			if ( tag == -1 )
			{
				++tag ;
				continue ;
			}
			depth[tag] = psum ;
			psum = 0 ;
			++tag ;
		}
	}
	//printf( "### %d %d\n", tag, pcnt ) ;
	// Build the lp system.
	int intervalCnt = pcnt - 1 ;
	int Ncol, *colno ;
	double *row ;
	int p ;
	lprec *lp ;

	Ncol = tcnt + intervalCnt ;
	colno = (int *)malloc( Ncol * sizeof( *colno ) ) ;
	row = (double *)malloc( Ncol * sizeof( *row ) ) ;
	lp = make_lp( 0, Ncol ) ;
	set_add_rowmode( lp, TRUE ) ;

	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		int k = 0 ;
		for ( j = 0 ; j < tcnt ; ++j )
		{
			for ( p = 0 ; p < transcripts[j].ecnt ; ++p )
			{
				if ( points[i] >= exons[ transcripts[j].eid[p] ].start &&
					points[i + 1] <= exons[transcripts[j].eid[p] ].end )
				{
					break ;
				}
				//printf( "%d: %d %d %d %d\n", j, points[i], points[i + 1], exons[transcripts[j].eid[p] ].start, exons[ transcripts[j].eid[p] ].end ) ;
			}

			if ( p >= transcripts[j].ecnt )
				continue ;
			colno[k] = j + 1 ;
			row[k] = 1 ;
			++k ;
		}
		
		colno[k] = tcnt + i + 1 ;
		row[k] = -1 ;
		add_constraintex(lp, k + 1, row, colno, LE, (double)depth[i] / ( points[i + 1] - points[i] + 1 ) ) ;

		row[k] = 1 ;
		add_constraintex(lp, k + 1, row, colno, GE, (double)depth[i] / ( points[i + 1] - points[i] + 1 ) ) ;
	}

	/*for ( i = 0 ; i < tcnt ; ++i )
	{
		colno[0] = i + 1 ;
		row[i] = 1 ;
		add_constraintex( lp, 1, row, colno, GE, 1 ) ;
	}*/

	for ( i = 0 ; i < tcnt ; ++i )
	{
		tmp = 0 ;
		for ( p = 0 ; p < transcripts[i].ecnt ; ++p )
			tmp += ( exons[ transcripts[i].eid[p] ].end - exons[transcripts[i].eid[p] ].start + 1 ) ;
		colno[i] = i + 1 ;
		row[i] = tmp ;
	}
	add_constraintex( lp, tcnt, row, colno, LE, sum + tcnt ) ;
	add_constraintex( lp, tcnt, row, colno, GE, sum - tcnt ) ;

	set_add_rowmode( lp, FALSE ) ;
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		colno[i] = tcnt + i + 1 ;
		row[i] = 1 ;
	}
	set_obj_fnex( lp, intervalCnt, row, colno ) ;
	set_minim( lp ) ;
	set_verbose( lp, IMPORTANT ) ;
	
	solve( lp ) ;

	get_variables( lp, row ) ;
	//write_LP( lp, stdout ) ;
	for ( i = 0 ; i < tcnt ; ++i )
		transcripts[i].abundance = row[i] ;
	free( colno ) ;
	free( row ) ;
	delete_lp( lp ) ;
	free( depth ) ;
}

int enumTid[MAX_TRANSCRIPT] ;

void EnumerateSetCover( int tag, int tcCnt, int depth, BitTable &destCover, struct _setcover &best,
	BitTable *btable, struct _transcript *transcripts, int tcnt )
{
	int i, j, k ;
	BitTable cover( tcCnt ) ;
	cover.Reset() ;
	for ( i = 0 ; i < depth ; ++i )
		cover.Or( btable[ enumTid[i] ] ) ;
	//printf( "# %d\n", cover.Count() ) ;
	if ( cover.IsEqual( destCover ) )
	{
		double weight = 0 ;
		for ( i = 0 ; i < depth ; ++i )
			weight += transcripts[ enumTid[i] ].abundance ;
		if ( weight < best.weight )
		{
			best.weight = weight ;
			best.cnt = depth ;
			for ( i = 0 ; i < depth ; ++i )
				best.choose[i] = enumTid[i] ;
		}
		//printf( "%d %d %d\n", tcCnt, cover.Count(), destCover.Count() ) ;
		return ;
	}
	//printf( "hi %d %d %d\n", tag, depth, destCover.Count() ) ;
	for ( i = tag ; i < tcnt ; ++i ) //- ( best.cnt - depth ) + 1 ; ++i )
	{
		//printf( "hi %d %d %d\n", i, depth, destCover.Count() ) ;
		enumTid[depth] = i ;
		double weight ;
		cover.Reset() ;
		
		for ( j = 0 ; j <= depth ; ++j )
		{
			cover.Or( btable[ enumTid[j] ] ) ;
			weight += transcripts[ enumTid[j] ].abundance ;
		}
		for ( j = tag + 1 ; j < tcnt ; ++j )
			cover.Or( btable[j] ) ;

		if ( !cover.IsEqual( destCover ) || weight > best.weight )
			continue ;
		EnumerateSetCover( i + 1, tcCnt, depth + 1, destCover, best, btable, transcripts, tcnt ) ;
	}
}

// Test whether a constraints is compatible with the transcript.
// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
int IsConstraintInTranscript( int *eid, int ecnt, const struct _transcriptConstraint &c, struct _exon exons[], int *toEndDist = NULL )
{
	int i, j, k, tag ;
	if ( c.type == 0 )
	{
		for ( i = BinarySearchExons( c.leftInfo[0][0], eid, ecnt, exons ) ; i < ecnt ; ++i )
		{
			if ( eid[i] < c.info[0] )
				break ;
			if ( eid[i] == c.info[0] )
				return 1 ;
		}
		return 0 ;
	}
	else if ( c.type == 1 )
	{
		int tag = 0 ;
		bool first = true, flag ;
		int leftInd, rightInd ;
		if ( c.strand != -1 && exons[ eid[0] ].strand != -1 && c.strand != exons[ eid[0] ].strand )
			return 0 ;
		if ( c.lcnt )
		{
			k = 0 ;
			for ( tag = BinarySearchExons( c.leftInfo[0][0], eid, ecnt, exons ) ; tag < ecnt ; ++tag )
			{
				if ( c.leftInfo[0][0] < exons[ eid[tag] ].start )
					return 0 ;

				if ( c.leftInfo[k][0] >= exons[ eid[tag] ].start
						&& c.leftInfo[k][1] <= exons[ eid[tag] ].end 
						&& ( k == c.lcnt - 1 || c.leftInfo[k][1] == exons[ eid[tag]].end ) 
						&& ( k == 0 || c.leftInfo[k][0] == exons[ eid[tag] ].start ) )
					break ;
			}
			if ( tag >= ecnt )
				return 0 ;
			++tag ;
			for ( k = 1 ; k < c.lcnt && tag < ecnt ; ++k, ++tag )
			{
				if ( c.leftInfo[k][0] >= exons[ eid[tag] ].start
						&& c.leftInfo[k][1] <= exons[ eid[tag] ].end 
						&& ( k == c.lcnt - 1 || c.leftInfo[k][1] == exons[ eid[tag]].end ) 
						&& ( k == 0 || c.leftInfo[k][0] == exons[ eid[tag] ].start ) )
					continue ;
				break ;
			}

			if ( k < c.lcnt && tag >= ecnt )
			{
				return 2 ;
			}
			if ( k < c.lcnt )
				return 0 ;
			leftInd = tag - 1 ;
		}
		//printf( "### %d %d %d %d %d\n", i, j, transcripts[i].ecnt, tc[j].lcnt, tc[j].rcnt ) ;		

		if ( c.rcnt )
		{
			int len = 0 ;
			if ( c.lcnt && c.rightInfo[0][0] > exons[ eid[ ecnt - 1] ].end )
			{
				// Test the insert size.
				len = 0 ;
				for ( i = leftInd + 1 ; i < ecnt ; ++i ) 
					len += ( exons[ eid[i] ].end - exons[ eid[i] ].start + 1 ) ;
				len += ( exons[ eid[ leftInd ] ].end - c.leftInfo[ c.lcnt - 1 ][1] + 1 ) ;
				if ( len >= FRAG_LENGTH + 2 * FRAG_STD - 2 * READS_LENGTH )
					return 0 ;
				return 2 ;
			}
			k = 0 ;
			for ( tag = BinarySearchExons( c.rightInfo[0][0], eid, ecnt, exons ) ; tag < ecnt ; ++tag )
			{
				if ( c.rightInfo[0][0] < exons[ eid[ tag ] ].start )
					return 0 ;

				if ( c.rightInfo[k][0] >= exons[ eid[tag] ].start
						&& c.rightInfo[k][1] <= exons[ eid[tag] ].end 
						&& ( k == c.rcnt - 1 || c.rightInfo[k][1] == exons[ eid[tag]].end ) 
						&& ( k == 0 || c.rightInfo[k][0] == exons[ eid[tag] ].start ) )
					break ;
			}

			if ( tag >= ecnt )
				return 0 ;
			rightInd = tag ;
			++tag ;
			// Test the insert size.
			len = 0 ;
			for ( k = leftInd + 1 ; k <= rightInd - 1 ; ++k )
			{
				len += ( exons[ eid[k] ].end - exons[ eid[k] ].start + 1 ) ;
			}

			if ( leftInd < rightInd )
			{
				len += ( exons[ eid[ leftInd ] ].end - c.leftInfo[ c.lcnt - 1 ][1] + 1 ) ;
				len += ( c.rightInfo[0][0] - exons[ eid[ rightInd ] ].start + 1 ) ;
			}

			if ( len >= FRAG_LENGTH + 2 * FRAG_STD - 2 * READS_LENGTH )
				return 0 ;

			for ( k = 1 ; k < c.rcnt && tag < ecnt ; ++k, ++tag )
			{
				if ( c.rightInfo[k][0] >= exons[ eid[tag] ].start
						&& c.rightInfo[k][1] <= exons[ eid[tag] ].end 
						&& ( k == c.rcnt - 1 || c.rightInfo[k][1] == exons[ eid[tag]].end ) 
						&& ( k == 0 || c.rightInfo[k][0] == exons[ eid[tag] ].start ) )
					continue ;
				break ;
			}

			if ( k < c.rcnt && tag >= ecnt )
				return 2 ;
			if ( k < c.rcnt )
				return 0 ;
			//printf( "success\n" ) ;
		}
		//if ( k >= tc[j].rcnt && tc[j].rcnt )
		//	continue ;
		/*if ( tc[j].leftInfo[0][1] == 2742878 && tc[j].leftInfo[1][0] == 2757641 )
		  {
		// printf( "%d %d %d %d\n", i, j, k, tc[j].lcnt ) ;
		//exit( 1 ) ;
		}*/
		return 1 ;
	}
	else if ( c.type == 2 )
	{
		if ( c.strand != exons[ eid[0] ].strand )
			return 0 ;
		int pairs[2][2] ;
		for ( i = BinarySearchExons( c.leftInfo[0][0], eid, ecnt, exons ) ; i < ecnt - 2 ; ++i )
		{
			if ( c.leftInfo[0][0] < exons[ eid[i] ].end )
				break ;

			pairs[0][0] = exons[ eid[i] ].end ;
			pairs[0][1] = exons[ eid[i + 1] ].start ;
			pairs[1][0] = exons[ eid[i + 1] ].end ;
			pairs[1][1] = exons[ eid[i + 2] ].start ;

			for ( j = 0 ; j < 2 ; ++j )
			{
				if ( pairs[j][0] != c.leftInfo[j][0] || pairs[j][1] != c.leftInfo[j][1] )
					break ;
			}

			if ( j >= 2 )
				return 1 ;
		}

		if ( exons[ eid[ecnt -1 ] ].end == c.leftInfo[0][0] )
			return 2 ;
		if ( ecnt >= 2 && exons[ eid[ecnt - 2] ].end == c.leftInfo[0][0] &&
			exons[ eid[ecnt - 1] ].start == c.leftInfo[0][1] && 
			exons[ eid[ecnt - 1] ].end == c.leftInfo[1][0] )
			return 2 ;
	}
	return 0 ;
}

/**
   Compare two transcripts.
*/
int TranscriptsDifference( const struct _transcript &t1, const struct _transcript &t2 )
{
	int i ;
	int cnt = t1.ecnt < t2.ecnt ? t1.ecnt : t2.ecnt ;	
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( t1.eid[i] != t2.eid[i] )
			return t1.eid[i] - t2.eid[i] ;
	}

	return t1.ecnt - t2.ecnt ;
}

// For qsort
int CompTranscripts( const void *p1, const void *p2 )
{
	return TranscriptsDifference( *( struct _transcript * )p1, *( struct _transcript * )p2 ) ;
}

int CompTranscriptsAbundance( const void *p1, const void *p2 )
{
	if ( ( ( struct _transcript *)p1 )->abundance != ( ( struct _transcript *)p2 )->abundance )
	{
		if ( ( ( struct _transcript *) p1 )->abundance < ( ( struct _transcript *)p2 )->abundance )
			return 1 ;
		else
			return -1 ;
	}
	else
		return 0 ;
}

int GetTranscriptExonWeight( struct _exonNode nodes[], const struct _transcript &t )
{
	int ret = 0 ;
	int i ;
	for ( i = 0 ; i < t.ecnt ; ++i )
		ret += nodes[ t.eid[i] ].weight ;
	return ret ;
}

inline double ComputeScore( double cnt, double a, double A )
{
	return cnt * ( 1 + pow( a / A, 0.25 ) ) ;
}

/**
	Pick the minimum # of transcripts in region [tstart, tend] that satisfy 
	all the constraints.
*/
void PickTranscripts( char *chrom, struct _exonNode nodes[], struct _exon exons[],
	struct _transcript *alltranscripts, const int &atcnt, struct _transcript *transcripts, int &tcnt, 
	struct _transcriptConstraint *tc, const int &tcCnt )
{
	int i, j, k ;
	int tag ;
	int mate ;
	//struct _transcriptConstraint *tc = transcriptConstraints ;
	BitTable *btable = new BitTable[ atcnt ] ; 
	//BitTable *btable = (BitTable *)malloc( sizeof( *btable ) * atcnt ) ;
	//delete [] btable ;
	//int *tcUsedCnt ;
	//tcUsedCnt = ( int * )malloc( sizeof( int ) * tcCnt ) ;
	//memset( tcUsedCnt, 0, sizeof( tcUsedCnt ) ) ;
	for ( i = 0 ; i < atcnt ; ++i )
		btable[i].Init( tcCnt ) ;
	tcnt = 0 ;

#if ALL_TRANSCRIPT	
	tcCnt = 0 ;
	for ( i = 0 ; i < atcnt ; ++i )
	{
		int transcriptId = 1 ;
		alltranscripts[i].transcriptId = transcriptId ;
		OutputTranscript( chrom, exons, &alltranscripts[i] ) ;
	}
	delete [] btable ;
	return ;
#endif
	
	for ( i = 0 ; i < atcnt ; ++i )
	{
		for ( j = 0 ; j < tcCnt ; ++j )
		{	
			if ( IsConstraintInTranscript( alltranscripts[i].eid, alltranscripts[i].ecnt, tc[j], exons) == 1 )
				btable[i].Set( j ) ;
		}
	}
	
	BitTable bAll( tcCnt ) ;	
	for ( i = 0 ; i < atcnt ; ++i )
		bAll.Or( btable[i] ) ;
	//printf( "# %d\n", tcCnt ) ;
	/*for ( i = 0 ; i < tcCnt ; ++i )
	{
		if ( !bAll.Test(i) )
		{
			//tcUsedCnt[i] = MAX_READ ;
			tc[i].support = tc[i].abundance = 0 ;
		}
		printf( "## %d %lf:", i, tc[i].normAbund ) ;
		printf( "%d %d: ", tc[i].lcnt, tc[i].rcnt ) ;
		for ( j = 0 ; j < tc[i].lcnt ; ++j )
			printf( "(%d %d) ", tc[i].leftInfo[j][0], tc[i].leftInfo[j][1] ) ;
		for ( j = 0 ; j < tc[i].rcnt ; ++j )
			printf( "[%d %d] ", tc[i].rightInfo[j][0], tc[i].rightInfo[j][1] ) ;
		printf( "\n" ) ;

	}*/

	/*for ( i = 0 ; i < tcCnt ; ++i )
	{
		if ( tc[i].support == 1 )
		{
			printf( "%d %d: ", tc[i].lcnt, tc[i].rcnt ) ;
			for ( j = 0 ; j < tc[i].lcnt ; ++j )
				printf( "(%d %d) ", tc[i].leftInfo[j][0], tc[i].leftInfo[j][1] ) ;
			for ( j = 0 ; j < tc[i].rcnt ; ++j )
				printf( "(%d %d) ", tc[i].rightInfo[j][0], tc[i].rightInfo[j][1] ) ;
			printf( "\n" ) ;
		}
	}*/
	/*for ( i = 0 ; i < atcnt ; ++i )
	{
		double *abundances = (double *)malloc( sizeof( double ) * tcCnt ) ;
		int aused = 0 ;
		int cut ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( tc[j].support == MAX_READ || tc[j].support == 0 || !btable[i].Test(j))
				continue ;
			abundances[aused] = tc[j].normAbund ; //tc[j].support ;
			++aused ;	
		}

		qsort( abundances, aused, sizeof( abundances[0] ), CompDouble ) ;
		cut = abundances[ int( aused * 0.5 ) ] ;
		int abundance = MAX_READ ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			//if ( btable[i].Test(j) )
			//	printf( "%d %d %d\n", tc[j].abundance, tc[j].support, tc[j].type ) ;
			if ( tc[j].abundance >= cut && tc[j].abundance < abundance && btable[i].Test(j) )
				abundance = tc[j].abundance ;
		}
		if ( abundance == MAX_READ )
			abundance = 1 ;	
		printf( "%d %d : ", i, btable[i].Count() ) ;
		for ( j = 0 ; j < aused ; ++j )
			printf( "%lf ", (double)abundances[j] ) ;
		printf( "\n" ) ;
		alltranscripts[i].transcriptId = i ;
		//++i ;
		OutputTranscript( chrom, exons, &alltranscripts[i] ) ;
		free( abundances ) ;
	}
	return ;*/


	for ( i = 0 ; i < 0 ; ++i )
	{
		OutputTranscript( chrom, exons, &alltranscripts[i] ) ;
		
		printf( "constraints: " ) ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( tc[j].type != 1 )
				continue ;
			printf( "%d", btable[i].Test( j ), tc[j].lcnt, tc[j].rcnt ) ;
		}
		printf( ": %d", btable[i].Count() ) ;
		printf( "\n" ) ;
		/*printf( "evidence: " ) ;
		for ( j = 0 ; j < alltranscripts[i].ecnt - 2 ; ++j )
		{
			int tpair[2][2] ;
			tpair[0][0] = exons[ alltranscripts[i].eid[j] ].end ;
			tpair[0][1] = exons[ alltranscripts[i].eid[j + 1] ].start ;
			tpair[1][0] = exons[ alltranscripts[i].eid[j + 1] ].end  ;
			tpair[1][1] = exons[ alltranscripts[i].eid[j + 2] ].start ;
			printf( "%d", 2 - Evidence_Val( tpair[0], tpair[1] ) ) ;
		}
		printf( ": %d\n", transcriptWeight[i] ) ;*/
	}

	/*if ( 0 )//atcnt <= 20 )
	{
		BitTable destCover( tcCnt ) ;
		for ( i = 0 ; i < atcnt ; ++i )
			destCover.Or( btable[i] ) ;
		for ( i = 1 ; i <= atcnt ; ++i )
		{
			//printf( "hi %d %d %d\n", i, atcnt, destCover.Count() ) ;
			bestWeight = -1 ;
			EnumerateSetCover( 0, tcCnt, 0, i, destCover ) ;
			if ( bestWeight != -1 )
			{
				for ( j = 0 ; j < i ; ++j )
					transcripts[j] = alltranscripts[ bestEnumTid[j] ] ;
				tcnt = i ;
				break ;
			}
		}
	}
	else*/ if ( 1 )
	{
		double max[BEST_COVER + 1] ;
		int  maxtag[BEST_COVER + 1] ; // TODO: Make it parameter
		double maxcnt[BEST_COVER+1] ;
		int maxExonWeight[BEST_COVER + 1] ;
		int maxAllCnt[BEST_COVER + 1] ;
		double cover ;
		double bestCnt ;
		double maxAbundance = -1 ; // The abundance of the most-abundant transcript.
		bool solveRemaining = false ;
		int remainingLen = 0 ; // The transcript remaining must be shorter than this.

		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = MAX_READ ;
			int tag = -1 ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				if ( btable[i].Test(j) && tc[j].type == 1 && tc[j].abundance > 0 )
				{
					if ( tc[j].abundance < value )
					{
						value = tc[j].abundance ;
						tag = j ;
					}
				}
			}
			if ( value == MAX_READ )
				value = 1e-6 ;
			/*printf( "%lf:", value ) ;
			for ( j = 0 ; j < tc[tag].lcnt ; ++j )
				printf( "(%d %d) ", tc[tag].leftInfo[j][0], tc[tag].leftInfo[j][1] ) ;
			for ( j = 0 ; j < tc[tag].rcnt ; ++j )
				printf( "[%d %d] ", tc[tag].rightInfo[j][0], tc[tag].rightInfo[j][1] ) ;
			printf( "\n" ) ;
			for ( j = 0 ; j < alltranscripts[i].ecnt ; ++j )
				printf( "<%d %d> ", exons[ alltranscripts[i].eid[j] ].start, exons[ alltranscripts[i].eid[j] ].end ) ;
			printf( "\n" ) ;*/
			if ( value > maxAbundance )
				maxAbundance = value ;
		}
		//printf( "maxAbundance = %lf %d\n", maxAbundance, atcnt ) ;
		//printf( "Choose: " ) ;
		while ( 1 ) 
		{

			//printf( "## %lf\n", maxAbundance ) ;
			/*printf( "=============\n" ) ;
			  for ( i = 0 ; i < atcnt ; ++i )
			  {
			  printf( "%d %d %d %d %d\n", alltranscripts[i].transcriptId, btable[i].Count(), transcriptWeight[ alltranscripts[i].transcriptId ], geneWeight[ alltranscripts[i].geneId ], alltranscripts[i].ecnt ) ;
			  }
			  printf( "=============\n" ) ;*/

			for ( i = 0 ; i < BEST_COVER ; ++i )
			{
				max[i] = -1 ;
				maxtag[i] = -1 ;
				maxcnt[i] = -1 ;
				maxExonWeight[i] = -1 ;
				maxAllCnt[i] = -1 ;
			}

			for ( i = 0 ; i < atcnt ; ++i )
			{
				double penalty ;
				double value = MAX_READ ;
				double cnt ;
				int allCnt ;
				double evidence = 0 ;

				//penalty = 1 ;

				cnt = 0 ;
				allCnt = 0 ;
				
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					/*if ( tc[j].abundance > 0 && btable[i].Test(j) )
					  {
					  if ( tc[j].abundance < value )
					  value = tc[j].abundance ;
					  ++cnt ;
					  }
					  else if ( tc[j].abundance == 0 && tc[j].support != MAX_READ && btable[i].Test(j) )
					  {
					  value = 1 ;
					//++cnt ;
					}*/
					
					if ( btable[i].Test(j) )
					{
						/*if ( tc[j].abundance > 0 )
						{
							if ( tc[j].abundance < value )
								value = tc[j].abundance ;
						}
						else if ( tc[j].abundance == 0 && tc[j].support != MAX_READ ) 
						{
							value = 1 ;
						}*/
						//if ( tc[j].support < value )
						//	value = tc[j].support ;
						if ( tc[j].normAbund < value && tc[j].type == 1 )
						{
							value = tc[j].normAbund ;
						}
						if ( tc[j].abundance > 0 && tc[j].type == 1 )
						{
							//cnt += 1 ; /// (double)( tcUsedCnt[j] + 1 )  ;
							//if ( tc[j].lcnt != 0 && tc[j].rcnt != 0 )
							//	cnt += 1 ;

							cnt += tc[j].effectiveCount ;
						}
						else if ( tc[j].type == 2 )
						{
							//TODO: partial evidence
							/*for ( k = 0 ; k < tcCnt ; ++k )
							{
								if ( tc[k].type != 1 || tc[k].abundance <= 0 )
									continue ;
								if ( ( tc[k].leftInfo[0][0] > tc[j].leftInfo[1][1] || 
											tc[k].leftInfo[ tc[k].lcnt - 1 ][1] < tc[j].leftInfo[0][0] ) &&
										( tc[k].rightInfo[0][0] > tc[j].leftInfo[1][1] || 
										  tc[k].rightInfo[ tc[k].rcnt - 1 ][1] < tc[j].leftInfo[0][0] )
								   )
								   continue ;
								if ( btable[i].Test(k) )
									break ;
							}
							if ( k < tcCnt )*/
							evidence += tc[j].effectiveCount ;
						}
						//else if ( tc[j].abundance > 0 && tc[j].support == MAX_READ )
						//	++cnt ;
						/*else if ( tc[j].abundance == 0 && tc[j].support != MAX_READ )
						{
							cnt = 0 ;
							break ;	
						}*/

						if ( tc[j].type == 1 || tc[j].type == 2 )
							allCnt += tc[j].effectiveCount ;
					}
				}
				if ( value == MAX_READ )
					value = 1e-6 ;
				//int tmp = btable[i].Count() ;/// penalty ;
				penalty = 1 ;
				/*if ( value > maxAbundance )
				{
					printf( "ERROR: %lf %lf\n", value, maxAbundance ) ;
					cnt = 0 ;
					evidence = 0 ;
					value = MAX_READ ;
					for ( j = 0 ; j < tcCnt ; ++j )
					{
						if ( btable[i].Test(j) )
						{
							printf( "hi %d %d %d %lf\n", i, j, tc[j].type, tc[j].normAbund ) ;
							if ( tc[j].normAbund < value )
								value = tc[j].normAbund ;
							if ( tc[j].abundance > 0 && tc[j].type != 2 )
							{
								cnt += 1 ; /// (double)( tcUsedCnt[j] + 1 )  ;
							}
							else if ( tc[j].type == 2 )
								evidence += 1 ;	
						}
					}
					printf( "%lf %lf %lf\n", cnt, evidence, value ) ;
					exit( 1 ) ;
				}*/

				double tmp ;
				int exonWeight = GetTranscriptExonWeight( nodes, alltranscripts[i] ) ;
				if ( !USE_SET_COVER )
				{
					/*if ( !solveRemaining )
						tmp = ( cnt + evidence ) / ( ( 2 - value / maxAbundance ) ) ; /// penalty ;
					else
						tmp = ( cnt + evidence / (double)tcCnt ) / ( ( 2 - value / maxAbundance ) ) ; /// penalty ;*/
					if ( !solveRemaining )
						tmp = ComputeScore( cnt + evidence, value, maxAbundance ) ; 
					else
						tmp = ComputeScore( cnt + evidence / (double)tcCnt, value, maxAbundance ) ;
					/*if ( !solveRemaining )
						tmp = ( cnt + evidence ) / ( ( 2 - sqrt( value / maxAbundance ) ) ) ; /// penalty ;
					else
						tmp = ( cnt + evidence / (double)tcCnt ) / ( ( 2 - sqrt( value / maxAbundance ) ) ) ; /// penalty ;*/
				}
				else
				{
					if ( alltranscripts[i].ecnt > 2 )
						penalty = (double)( alltranscripts[i].ecnt - 2 - evidence ) 
							/ ( alltranscripts[i].ecnt - 2 ) + 0.1 ;
					else
						penalty = 1.1 ;	 
					tmp = cnt / penalty ; 	
				}
				//double tmp = cnt ;  
				//int tmp = cnt ;
				for ( k = BEST_COVER - 1 ; k >= 0 ; --k )
				{
					if ( tmp > max[k] )
					{
						max[k + 1] = max[k] ;
						maxtag[k + 1] = maxtag[k] ;
						maxcnt[k + 1] = maxcnt[k] ;
						maxExonWeight[k + 1] = maxExonWeight[k] ;
						maxAllCnt[k + 1] = maxAllCnt[k] ;
					}
					else if ( tmp == max[k] && ( cnt + evidence ) > maxcnt[k] )
					{	
						max[k + 1] = max[k] ;
						maxtag[k + 1] = maxtag[k] ;
						maxcnt[k + 1] = maxcnt[k] ;
						maxExonWeight[k + 1] = maxExonWeight[k] ;
						maxAllCnt[k + 1] = maxAllCnt[k] ;
					}
					else if ( tmp == max[k] && ( cnt + evidence ) == maxcnt[k] && allCnt > maxAllCnt[k] )
					{
						max[k + 1] = max[k] ;
						maxtag[k + 1] = maxtag[k] ;
						maxcnt[k + 1] = maxcnt[k] ;
						maxExonWeight[k + 1] = maxExonWeight[k] ;
						maxAllCnt[k + 1] = maxAllCnt[k] ;
					}
					else
						break ;
				}
				max[k + 1] = tmp ;
				maxtag[k + 1] = i ;
				maxcnt[k + 1] = cnt + evidence ;
				maxExonWeight[k + 1] = exonWeight ;
				maxAllCnt[k + 1] = allCnt ;
				if ( k + 1 == 0 )
					bestCnt = cnt ;
			} // end for atcnt
			
			if ( bestCnt <= 1e-6 && ( solveRemaining || evidences == NULL ) )
				break ;
			else if ( bestCnt <= 1e-6 )
			{
				//break ;
				solveRemaining = true ;
				remainingLen = MAX_EXON ;
				for ( i = 0 ; i < tcnt ; ++i )
				{
					if ( transcripts[i].ecnt == 1 )
						continue ;
					if ( transcripts[i].ecnt < remainingLen )
						remainingLen = transcripts[i].ecnt ;
				}
				continue ;
			}
			//OutputTranscript( chrom, exons, &alltranscripts[maxtag] ) ;
			//BitTable btag( tcCnt ) ;
			//BitTable btag ;
			//btag.Init( tcCnt ) ;
			for ( i = 0 ; i < BEST_COVER ; ++i )
			{
				if ( max[i] <= 1e-6 )
					break ;
				//if ( !solveRemaining || alltranscripts[maxtag[i]].ecnt < remainingLen )
				//{
				transcripts[ tcnt ] = alltranscripts[ maxtag[i] ] ;
					//transcripts[ tcnt ].transcriptId = transcriptId ;
					//++transcriptId ;
				transcripts[tcnt].abundance = 0 ;
				for ( j = 0 ; j < transcripts[tcnt].ecnt ; ++j )
					++nodes[ transcripts[tcnt].eid[j] ].weight ;	
				//}
				//printf( "%d\n", tcnt ) ;
				// Update the abundance
				double value = MAX_READ ;

				// NOTES: Actually, the transcript's abundance is 1 if it containts some depleted constraints.
				// However, if we only substract abundance of each constraints by 1, in the next iteration, we still
				// will choose the same transcript. So, it is safe to subtract the minimal non-zero abundance from each
				// compatible constraints.
				k = -1 ;
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( btable[ maxtag[i] ].Test(j) && tc[j].abundance > 0 )
					{
						if ( tc[j].abundance < value )
						{
							value = tc[j].abundance ;
							k = j ;
						}
					}
					else if ( btable[ maxtag[i] ].Test(j) && tc[j].abundance == 0 )
						value = 0 ;

					/*if ( tc[j].abundance > 0 && btable[ maxtag[i] ].Test(j) )
					{
						if ( tc[j].abundance < value )
							value = tc[j].abundance ;
					}
					else if ( tc[j].abundance == 0 && tc[j].support != MAX_READ && btable[maxtag[i]].Test(j) )
					{
						value = 0 ;
					}*/
				}
				if ( value < 1e-6 )
				{
					//--tcnt ;
					value = MAX_READ ;
					for ( j = 0 ; j < tcCnt ; ++j )
						if ( btable[ maxtag[i] ].Test(j) && tc[j].abundance > 0 )
						{
							if ( tc[j].abundance < value )
							{
								value = tc[j].abundance ;
								k = j ;
							}
						}
					//if ( value == MAX_READ )
					//	value = 1 ;
					//value = 1 ; 
				}

				//if ( value <= tc[k].normAbund * 0.5 )
				//	--tcnt ;
				/*double cnt = 0 ;
				int a = MAX_READ ;
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( btable[ maxtag[i] ].Test(j) && tc[j].abundance > 0 )
					{
						cnt += 1 ;/// (double)( tcUsedCnt[j] + 1 )  ;
						//printf( "%d\n", tcUsedCnt[j] ) ;
					}
					if ( btable[ maxtag[i] ].Test(j) && tc[j].support < a )
						a = tc[j].support ; 
				}
				printf( "%lf %lf\n" , cnt, max[i] ) ;
				printf( "\t" ) ;
				for ( j = 0 ; j < transcripts[tcnt - 1].ecnt ; ++j )
					printf( "%d ", transcripts[tcnt - 1].eid[j] ) ;
				printf( "\n" ) ;*/
				
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( btable[ maxtag[i] ].Test(j) && tc[j].abundance > 0 )
					{
						if ( tc[j].type == 0 ) // type-0 constraints
							tc[j].abundance = 0 ;
						else if ( tc[j].type == 1 )
						{
							tc[j].abundance -= 1 * value ;
							double factor = 1 ;
							if ( tc[j].lcnt > 0 && tc[j].rcnt > 0 )
								factor = 2 ;
							transcripts[ tcnt ].abundance += ( tc[j].support * value / tc[j].normAbund * factor ) ; 
							
							if ( solveRemaining || USE_SET_COVER )
								tc[j].abundance = 0 ;
						}
						if ( tc[j].abundance < 0 )
							tc[j].abundance = 0 ;
					}

					/*if ( btable[ maxtag[i] ].Test(j) )
					{
						++tcUsedCnt[j] ;
					}*/
				}
				//if ( value == MAX_READ )
				//	value = 1 ;

				//btag.Or( btable[ maxtag[i] ] ) ;
				//btable[ maxtag[i] ].Reset() ;
				
				++tcnt ;
				
				if ( tcnt >= MAX_TRANSCRIPT )
				{
					qsort( transcripts, tcnt, sizeof( struct _transcript ), CompTranscripts ) ;
					int cnt ;
					cnt = 1 ;
					for ( i = 1 ; i < tcnt ; ++i )
					{
						if ( TranscriptsDifference( transcripts[i], transcripts[i - 1] ) )
						{
							transcripts[cnt] = transcripts[i] ;
							++cnt ;
						}
						else
							transcripts[cnt - 1].abundance += transcripts[i].abundance ;
						//else
						//	free( transcripts[i].eid ) ;
					}					
					tcnt = cnt ;
				}
				if ( tcnt >= MAX_TRANSCRIPT )
				{
					delete [] btable ;
					return ;
				}
			} // end for BEST_COVER

			//for ( i = 0 ; i < BEST_COVER ; ++i )
			//	printf( "%d ", alltranscripts[ maxtag[i] ].transcriptId ) ;

			//btag.Not() ;
			//for ( i = 0 ; i < atcnt ; ++i )
			//	btable[i].And( btag ) ;
			/*for ( i = 0 ; i < tcCnt ; ++i )
				if ( tcUsedCnt[i] < 2 && tc[i].abundance > 0 )
					break ;
			if ( i >= tcCnt )
				break ;*/
		} // end for while(1)
	} // end for if(1)
	
	//Quantification( chrom, points, pcnt, exons ) ;
	//printf( "\n==============================\n\n" ) ;
	delete [] btable ;
	// Remove the repeated transcripts
	qsort( transcripts, tcnt, sizeof( struct _transcript ), CompTranscripts ) ;
	if ( tcnt == 0 )
		return ;
	int cnt ;
	cnt = 1 ;
	for ( i = 1 ; i < tcnt ; ++i )
	{
		if ( TranscriptsDifference( transcripts[i], transcripts[i - 1] ) )
		{
			transcripts[cnt] = transcripts[i] ;
			++cnt ;
		}
		else
			transcripts[cnt - 1].abundance += transcripts[i].abundance ;
	}
	tcnt = cnt ;
}


// Build the transcript constraints based on the exons[from..to] spanning tstart to tend.
// forceSingle means to force to separate the paired reads and use the single reads. 
bool IsTwoExonsInConstraint( int *eid, struct _transcriptConstraint &c, struct _exon exons[] )
{
	int i ;
	if ( c.strand != -1 && exons[ eid[0] ].strand != -1 && c.strand != exons[ eid[0] ].strand )
		return false ;
	for ( i = 0 ; i < c.lcnt - 1 ; ++i )
	{
		if ( c.leftInfo[i][0] >= exons[eid[0]].start && c.leftInfo[i][1] == exons[ eid[0] ].end )
			break ;
	}
	if ( i < c.lcnt - 1 )
	{
		if ( c.leftInfo[i + 1][0] == exons[eid[1]].start && c.leftInfo[i + 1][1] <= exons[ eid[1] ].end )
			return true ;
	}
	
	for ( i = 0 ; i < c.rcnt - 1 ; ++i )
	{
		if ( c.rightInfo[i][0] >= exons[eid[0]].start && c.rightInfo[i][1] == exons[ eid[0] ].end )
			break ;
	}
	if ( i < c.rcnt - 1 )
	{
		if ( c.rightInfo[i + 1][0] == exons[eid[1]].start && c.rightInfo[i + 1][1] <= exons[ eid[1] ].end )
			return true ;
	}
	return false ;
}

int BuildTranscriptConstraints( char *chrom, struct _exonNode nodes[], struct _exon exons[], 
			int tstart, int tend, int from, int to, bool forceSingle  )
{
	int i, j, k ;
	int tag ;
	int tcCnt = 0 ;
	int mate ;
	struct _transcriptConstraint *tc ;
	if ( forceSingle )
		tc = abundanceConstraints ;
	else
		tc = transcriptConstraints ;
	int extent[2] ;
	//printf( "file pointer: %lld\n", td_fpReads.fp ) ;
	if ( tend == -1 )
		return -1 ;

	if ( forceSingle ) // Later, the non-abundance can directly use the result.
	{ 
		trCnt = ExtractGeneReads( td_fpReads, chrom, tstart, tend, transcriptReads, extent ) ;
	}

	int *tmpPoints, *points ;	
	int pcnt = 0 ;
	tmpPoints = ( int * )malloc( 3 * ( to - from + 1 ) * sizeof( int ) ) ;
	points = ( int * )malloc( 3 * ( to - from + 1 ) * sizeof( int ) ) ;
	for ( i = from ; i <= to ; ++i )
	{
		tmpPoints[2 * ( i - from ) ] = exons[ i ].start ;
		tmpPoints[2 * ( i - from ) + 1] = exons[ i ].end ;
	}
	k = 2 * ( to - from + 1 ) ;
	qsort( tmpPoints, k, sizeof( int ), CompInt ) ;

	if ( to - from + 1 > 0 )
	{	
		points[0] = tmpPoints[0] ;
		//printf( "%d\n", points[0] ) ;
		pcnt = 1 ;
		for ( i = 1 ; i < k ; ++i )
		{
			if ( tmpPoints[i] == tmpPoints[i - 1] )
				continue ;
			points[pcnt] = tmpPoints[i] ;
			++pcnt ;
		}
	}
	else
		pcnt = 0 ;
	//for ( i = 0 ; i < pcnt ; ++i )
	//	printf( "%d\n", points[i] ) ;
	
	// All the exons must show up.
	if ( !forceSingle )
	{	
		for ( i = from ; i <= to ; ++i )
		{
			//for ( k = tcCnt - 1 ; k >= 0 && tc[k].info[0] > eidUsed[j] ; --k )
			//	;
			//if ( k >= 0 && tc[k].info[0] == transcripts[i].eid[j] )
			//	continue ;
			//for ( k = tcCnt ; k >= 0 && tc[k].info[0] > transcripts[i].eid[j] ; --k )
			//	tc[k + 1] = tc[k] ;

			alltc[tcCnt].type = 0 ;
			alltc[tcCnt].info[0] = i ;
			alltc[tcCnt].infoCnt = 1 ;
			alltc[tcCnt].strand = exons[ i ].strand ;
			alltc[tcCnt].support = MAX_READ ;

			alltc[tcCnt].leftInfo[0][0] = exons[i].start ;
			alltc[tcCnt].leftInfo[0][1] = exons[i].end ;
			++tcCnt ;
		}
	}

	// For the reads.
	/*for ( i = 0 ; i < trCnt ; ++i )
	{
		if ( reads[i].scnt && reads[i].splices[0][0] == 2742878 )
			printf( "%d %d\n", reads[i].splices[0][1], reads[i].mateInd ) ;
	}*/
	for ( i = 0 ; i < trCnt ; ++i )
	{
		int a, b, s, e ;
		struct _read &read = transcriptReads.reads[i / MAX_READ][i % MAX_READ] ;
		
		// TODO: Remove this line later.
		//if ( forceSingle )
		//	read.mateInd = -1 ;
		if ( !forceSingle && read.mateInd < i && read.mateInd != -1  )
			continue ;
		alltc[ tcCnt ].type = 1 ;
		alltc[ tcCnt ].lcnt = alltc[ tcCnt ].rcnt = 0 ;
		alltc[ tcCnt ].strand = read.strand ; //-1 ;
		k = 0 ;
	//printf( "begin %d %d %d %d\n", tcCnt, i, trCnt, j ) ;
		for ( j = 0 ; j <= read.scnt ; ++j )
		{
			if ( j == 0 )
				a = read.start ;
			else
				a = read.splices[j - 1][1] ;

			if ( j == read.scnt )
				b = read.end ;	
			else
				b = read.splices[j][0] ;
			
			if ( a == b )
			{
				if ( j == 0 )
					--a ;
				else if ( j == read.scnt )
					++b ;
			}
			
			if ( a < points[0] || b > points[ pcnt - 1 ] )
			{
				k = 0 ;
				break ;
			}
			for ( s = 0 ; s < pcnt && a >= points[s] ; ++s )
				;
			--s ;

			if ( ( j > 0 && points[s] != a ) )
			{
				k = 0 ;
				break ;
			}	
			//TODO: Look into the effect of allowing e=s. 
			for ( e = s ; e < pcnt && points[e] < b ; ++e ) // the first points[e] greater or equal b. 
				;
			if ( j < read.scnt && points[e] != b )
			{
				k = 0 ;
				break ;
			}
			
			//exit( 1 ) ;
			alltc[ tcCnt ].leftInfo[k][0] = points[s] ;
			alltc[ tcCnt ].leftInfo[k][1] = points[e] ;
			++k ;
		}
		alltc[ tcCnt ].lcnt = k ;	
		/*if ( k > 1 && !forceSingle )
		{
			printf( "%d %d\n", read.start, read.end ) ;
			exit( 1 ) ;
		}*/
	//printf( "end0 %d\n", tcCnt ) ;
		k = 0 ;

		if ( !forceSingle )
		{
			int matetag = read.mateInd / MAX_READ ;
			int mateoffset = read.mateInd - matetag * MAX_READ ;

			for ( j = 0 ; read.mateInd != -1 && j <= transcriptReads.reads[matetag][mateoffset].scnt  ; ++j )
			{
				struct _read &read = transcriptReads.reads[matetag][mateoffset] ;
				if ( j == 0 )
					a = read.start ;
				else
					a = read.splices[j - 1][1] ;

				if ( j == read.scnt )
					b = read.end ;	
				else
					b = read.splices[j][0] ;

				if ( a == b )
				{
					if ( j == 0 )
						--a ;
					else if ( j == read.scnt )
						++b ;
				}

				if ( a < points[0] || b > points[ pcnt - 1 ] )
				{
					k = 0 ;
					break ;
				}
				
				for ( s = 0 ; s < pcnt && a >= points[s] ; ++s )
					;
				--s ;
				if ( j > 0 && points[s] != a )
				{
					k = 0 ;
					break ;
				}

				for ( e = s + 1 ; e < pcnt && points[e] < b ; ++e )
					;
				if ( j < read.scnt && points[e] != b )
				{
					k = 0 ;
					break ;
				}
				alltc[ tcCnt ].rightInfo[k][0] = points[s] ;
				alltc[ tcCnt ].rightInfo[k][1] = points[e] ;
				++k ;
			}
		}
		alltc[ tcCnt ].rcnt = k ;
		alltc[ tcCnt ].support = 1 ;
		alltc[ tcCnt ].uniqSupport = (read.unique ? 1 : 0) ;
		
		if ( alltc[ tcCnt ].lcnt == 0 && alltc[ tcCnt ].rcnt != 0 )
		{
			// This happens when we cut the left soft boundary, the left read might also be cut and 
			// that part will be throw away.
			alltc[ tcCnt ].lcnt = alltc[tcCnt].rcnt ;
			memcpy( alltc[tcCnt].leftInfo, alltc[tcCnt].rightInfo, sizeof( alltc[tcCnt].rightInfo ) ) ;
			alltc[ tcCnt ].rcnt = 0 ;
		}

		if ( alltc[ tcCnt ].lcnt > 0 && alltc[ tcCnt ].rcnt > 0 )
		{
			if ( alltc[tcCnt].leftInfo[0][0] > alltc[tcCnt].rightInfo[0][0] )
			{
				// Should not happen
				int tmp[100][2] ;
				int tmp2 ;
				for ( j = 0 ; j < alltc[tcCnt].lcnt ; ++j )
				{
					tmp[j][0] = alltc[tcCnt].leftInfo[j][0] ;
					tmp[j][1] = alltc[tcCnt].leftInfo[j][1] ;
				}
				
				
				for ( j = 0 ; j < alltc[tcCnt].rcnt ; ++j )
				{
					alltc[tcCnt].leftInfo[j][0] = alltc[tcCnt].rightInfo[j][0] ;
					alltc[tcCnt].leftInfo[j][1] = alltc[tcCnt].rightInfo[j][1] ;
				}


				for ( j = 0 ; j < alltc[tcCnt].lcnt ; ++j )
				{
					alltc[tcCnt].rightInfo[j][0] = tmp[j][0] ;
					alltc[tcCnt].rightInfo[j][1] = tmp[j][1] ;
				}

				tmp2 = alltc[tcCnt].lcnt ;
				alltc[tcCnt].lcnt = alltc[tcCnt].rcnt ;
				alltc[tcCnt].rcnt = tmp2 ;
			}
		}
		if ( alltc[ tcCnt ].lcnt || alltc[ tcCnt ].rcnt )
		{
			++tcCnt ;
		}

		if ( tcCnt >= MAX_TRANSCRIPT_CONSTRAINT )
		{
			// If there are too many constraints, merge them now.
			qsort( alltc, tcCnt, sizeof( *alltc ), CompTranscriptConstraint ) ;
			k = 1 ;
			for ( j = 1 ; j < tcCnt ; ++j )
			{
				if ( TranscriptConstraintDifference( &alltc[j], &alltc[j - 1] ) )
				{
					alltc[k] = alltc[j] ;
					++k ;
				}
				else
				{
					alltc[k - 1].support += alltc[j].support ;
					alltc[k - 1].uniqSupport += alltc[j].uniqSupport ;
				}
			}
			//printf( "Clean: %d->%d\n", tcCnt, k ) ;
			tcCnt = k ;
		}
		if ( tcCnt >= MAX_TRANSCRIPT_CONSTRAINT )
			break ;
	}
	//printf( "end1\n" ) ;
	//free( tmpPoints ) ;
	//free( points ) ;
	/*for ( i = 0 ; i < tcCnt ; ++i )
	{
		printf( "%d %d %d ", alltc[i].type, alltc[i].lcnt, alltc[i].rcnt ) ;
		for ( j = 0 ; j < alltc[i].lcnt ; ++j )
			printf( "(%d %d)", alltc[i].leftInfo[j][0], alltc[i].leftInfo[j][1] ) ;
		for ( j = 0 ; j < alltc[i].rcnt ; ++j )
			printf( "(%d %d)", alltc[i].rightInfo[j][0], alltc[i].rightInfo[j][1] ) ;
		printf( "\n" ) ;
	}*/
	// Build the evidence constraints
	if ( !forceSingle && evidences != NULL )
	{
		int a, b, c ;
		int pairs[2][2] ;
		for ( ; eviTag < eviCnt && evidences[ eviTag ].exons[ evidences[eviTag].ecnt - 1 ].end < tstart ; ++eviTag )
			;
		for ( a = from ; a <= to ; ++a )	
		{
			for ( i = 0 ; i < nodes[a].ncnt ; ++i )
			{
				b = nodes[a].next[i] ;
				for ( j = 0 ; j < nodes[b].ncnt ; ++j )
				{
					c = nodes[b].next[j] ;

					pairs[0][0] = exons[a].end ;
					pairs[0][1] = exons[b].start ;
					pairs[1][0] = exons[b].end ;
					pairs[1][1] = exons[c].start ;

					if ( IsInEvidence( pairs, exons[a].strand ) )
					{
						alltc[ tcCnt ].type = 2 ;
						alltc[ tcCnt ].leftInfo[0][0] = pairs[0][0] ;
						alltc[ tcCnt ].leftInfo[0][1] = pairs[0][1] ;
						alltc[ tcCnt ].leftInfo[1][0] = pairs[1][0] ;
						alltc[ tcCnt ].leftInfo[1][1] = pairs[1][1] ;
						alltc[ tcCnt ].support = 1 ;
						alltc[ tcCnt ].strand = exons[a].strand ;
						++tcCnt ;
						continue ;
					}
				}
			}
		}
	}


	if ( !tcCnt )
	{
		free( tmpPoints ) ;
		free( points ) ;
		return 0 ;
	}
	qsort( alltc, tcCnt, sizeof( *alltc ), CompTranscriptConstraint ) ;

	k = 1 ;
	tc[0] = alltc[0] ;
	for ( i = 1 ; i < tcCnt ; ++i )
	{
		if ( TranscriptConstraintDifference( &alltc[i], &alltc[i - 1] ) )
		{
		//	printf( "	support = %d\n", tc[k - 1].support ) ;
			tc[k] = alltc[i] ;
			++k ;
		}
		else
		{
			tc[k - 1].support += alltc[i].support ;
			tc[k - 1].uniqSupport += alltc[i].uniqSupport ;
		}
	}
	tcCnt = k ;

	/*for ( i = 0 ; i < tcCnt ; ++i )
	{
		printf( "%d %d %d ", tc[i].type, tc[i].lcnt, tc[i].rcnt ) ;
		for ( j = 0 ; j < tc[i].lcnt ; ++j )
			printf( "(%d %d)", tc[i].leftInfo[j][0], tc[i].leftInfo[j][1] ) ;
		for ( j = 0 ; j < tc[i].rcnt ; ++j )
			printf( "(%d %d)", tc[i].rightInfo[j][0], tc[i].rightInfo[j][1] ) ;
		if ( i > 0 )
			printf( "%d", (bool)TranscriptConstraintDifference( &tc[i], &tc[i - 1] ) ) ;
		printf( "\n" ) ;
	}
	exit( 1 ) ;*/
		//printf( "### %d\n", tcCnt ) ;
	// Calculate the end of type1 constraints 
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		if ( tc[j].type == 0 )
		{
			tc[j].info[1] = exons[ tc[j].info[0] ].end ;
			//if ( !forceSingle )
			//	printf( "hi %d %d\n", j, tcCnt ) ;
		}
		else if ( tc[j].type == 1 )
		{	
			int a = -1, b = -1 ;
			if ( tc[j].lcnt )
				a = tc[j].leftInfo[ tc[j].lcnt - 1][1] ;			
			if ( tc[j].rcnt )
				b = tc[j].rightInfo[ tc[j].rcnt - 1 ][1] ;
			tc[j].info[1] = a > b ? a : b ;
		}
		else if ( tc[j].type == 2 )
		{
			tc[j].info[1] = tc[j].leftInfo[1][1] ;
		}
	}

	for ( j = 0 ; j < tcCnt ; ++j )
	{
		tc[j].effectiveCount = 1 ;
	}
	/*if ( forceSingle )
	{
		printf( "# %d\n", tcCnt ) ;
		for ( i = 0 ; i < tcCnt ; ++i )
		{
			printf( "%d: ", tc[i].support ) ;
			printf( "%d %d: ", tc[i].lcnt, tc[i].rcnt ) ;
			for ( j = 0 ; j < tc[i].lcnt ; ++j )
				printf( "(%d %d) ", tc[i].leftInfo[j][0], tc[i].leftInfo[j][1] ) ;
			for ( j = 0 ; j < tc[i].rcnt ; ++j )
				printf( "(%d %d) ", tc[i].rightInfo[j][0], tc[i].rightInfo[j][1] ) ;
			printf( "\n" ) ;
		}
		printf( "\n" ) ;
	}*/

	//printf( "%d %d %d\n", tcnt, tcCnt, k ) ;	
	if ( forceSingle )
	{
		// Caculate the normalized abundance
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int a = -1, b = -1 ;
			for ( i = 0 ; i < pcnt ; ++i )
			{
				if ( points[i] == tc[j].leftInfo[0][0] )
					break ;
			}

			a = i ;
			for ( i = a + 1 ; i < pcnt ; ++i )
			{
				if ( points[i] == tc[j].leftInfo[ tc[j].lcnt - 1 ][1] )
					break ;
			}
			b = i ;
			
			if ( a >= pcnt || b >= pcnt )
			{
				/*printf( "ERROR: %d %d %d (%d, %d)\n", a, b,pcnt, j, tcCnt ) ;
				for ( i = 0 ; i < pcnt ; ++i )
					printf( "%d ", points[i] ) ;
				printf( "\n%d %d\n", tc[j].leftInfo[0][0], tc[j].leftInfo[ tc[j].lcnt - 1 ][1] ) ;
				exit( 1 ) ;*/
				tc[j].support = -1 ;
				continue ;
			}

			if ( b == a + 1 )
			{
				if ( points[b] - READS_LENGTH + 2 - points[a] > 0 )
				{
					tc[j].normAbund = tc[j].support / (double)( points[b] - READS_LENGTH + 2 - points[a] ) ;
				}
				else
				{
					tc[j].normAbund = tc[j].support / (double)( READS_LENGTH - ( points[b] - points[a] ) );
				}
				/*for ( k = 0 ; k < tc[j].lcnt ; ++k )
				  printf( "(%d %d) ", tc[j].leftInfo[k][0], tc[j].leftInfo[k][1] ) ;
				  printf( ": %lf: %d\n", tc[j].normAbund, tc[j].strand ) ;*/
				continue ;
			}

			int midlength = 0 ; // The concatenated length from point[a+1] to point[b-1]+1
			int length ;
			if ( tc[j].lcnt == 1 )
			{
				for ( i = a + 2 ; i < b ; ++i )
					midlength += ( points[i] - points[i - 1] ) ;
				//++midlength ;
				midlength += 2 ;
			}
			else
			{
				for ( i = 1 ; i < tc[j].lcnt - 1 ; ++i )
					midlength += ( tc[j].leftInfo[i][1] - tc[j].leftInfo[i][0] + 1 ) ;
				for ( i = a + 2 ; i < pcnt && points[i] <= tc[j].leftInfo[0][1] ; ++i )
				{
					midlength += ( points[i] - points[i - 1] ) ;
				}

				for ( i = b - 1 ; points[i - 1] >= tc[j].leftInfo[ tc[j].lcnt - 1 ][0] ; --i )
					midlength += ( points[i] - points[i - 1] ) ;
					
				midlength += 3 ;
			}
			
			int left = points[a + 1] - ( READS_LENGTH - midlength ) ;
			int overb1 = points[b-1] + 1 ;
			if ( left < points[a] )
			{
				overb1 += points[a] - left ;
				left = points[a] ;
			}
			int right = points[a + 1] - 1 ;
			
			
			/*if ( points[b - 1] + 1 + ( right - left ) > points[b] )
			{
				right = left + ( points[b] - points[b - 1] - 1 ) ;
			}*/
			if ( right - left > points[b] - overb1 )
			{
				right = points[b] - overb1 + left ;
			}
			
			if ( right - left + 1 <= 0 )
			{
				// Caused by the actual read length for the constraints is longer than the normal read length.
				// Or the anchor size of a splice junction is only 1.
				right = left ;
			}
			tc[j].normAbund = tc[j].support / (double)( right - left + 1 ) ;
			
			
			/*for ( k = 0 ; k < tc[j].lcnt ; ++k )
				printf( "(%d %d) ", tc[j].leftInfo[k][0], tc[j].leftInfo[k][1] ) ;
			printf( ": %lf: %d\n", tc[j].normAbund, tc[j].strand ) ;*/
			//for ( k = 0 ; k < tc[j].rcnt ; ++k )
			//	printf( "R(%d %d)\n", tc[j].rightInfo[k][0], tc[j].rightInfo[k][1] ) ;
			/*if ( right < left )
			{
				for ( i = a ; i <= b ; ++i )
					printf( "%d ", points[i] ) ;
				printf( "\n" ) ;
				for ( i = 0 ; i < tc[j].lcnt ; ++i )
					printf( "(%d %d) ", tc[j].leftInfo[i][0], tc[j].leftInfo[i][1] ) ; 
				printf( "\n" ) ;
				printf( "err %d: %d %d %d: %d\n", j, tc[j].support, left, right, tc[j].lcnt ) ;
				exit( 1 ) ;
			}*/
		}
	}

	if ( 0 ) //!forceSingle )
	{
		// Build the penalty constraints.
		// Find candidate part in the splice graph.  
		// Candidate part: a non-branch path, and the exons outside the path sinks into the path
		for ( i = from + 1 ; i <= to - 1 ; ++i )
		{
			// Find the start nodes for the candidate part
			if ( nodes[i].pcnt < 2 ) 
				continue ;
			for ( j = 0 ; j < nodes[i].pcnt ; ++j )				
			{
				//This is a wrong test, the next might be another node in the previous node of the start node
				if ( nodes[ nodes[i].prev[j] ].ncnt > 1 )
					break ;
			}
			if ( j < nodes[i].pcnt )
				continue ;

			int start = i ;
			//printf( "hi1 %d %d\n", exons[start].start, exons[start].pcnt ) ;
			for ( k = i ; k < to ; ++k )
			{
				// It seems we can allow edges to come into the part
				/*if ( nodes[k].pcnt > 1 )
				  break ;*/
				if ( nodes[k].ncnt > 1 )
					break ;
			}
			if ( k >= to )
				continue ;
			for ( j = 0 ; j < nodes[k].ncnt ; ++j )
			{
				if ( nodes[ nodes[k].next[j] ].pcnt > 1 )
					break ;
			}
			if ( j < nodes[k].ncnt )
				continue ;
			int end = k ;
			//printf( "hi2 %d %d\n", exons[start].start, exons[k].end ) ;

			// The candidate part is from 'start' to 'end'.
			// Now we need to check whether they are well supported
			// by other constraints.
			// Goal: a constraint can cover the entering and leaving edge at the same time for the candidate part,
			//       there union should cover all the edges entering and leaving the part.
			int *inExons = ( int * )malloc( sizeof( int ) * nodes[start].pcnt ) ;
			int *outExons = ( int * )malloc( sizeof( int ) * nodes[end].ncnt ) ;
			int inCnt = nodes[start].pcnt ;
			int outCnt = nodes[end].ncnt ;
			bool *covered = ( bool * )malloc( sizeof( bool ) * inCnt * outCnt ) ;
			
			memset( covered, false, sizeof( bool ) * inCnt * outCnt ) ;
			for ( int in = 0 ; in < inCnt ; ++in )
			{
				for ( int out = 0 ; out < outCnt ; ++out )
				{
					for ( k = 0 ; k < tcCnt ; ++k ) // TODO: not start from 0.
					{
						if ( tc[k].type != 1 )
							continue ;
						if ( tc[k].lcnt + tc[k].rcnt < 3 )
							continue ;
						int eid[2] ;
						eid[0] = nodes[start].prev[in] ;
						eid[1] = start ;
						if ( !IsTwoExonsInConstraint( eid, tc[k], exons ) )
							continue ;
						eid[0] = end ;
						eid[1] = nodes[ end ].next[out] ;
						if ( IsTwoExonsInConstraint( eid, tc[k], exons ) )
							break ;
					}

					if ( k < tcCnt )
					{
						//printf( "hi2 %d %d %d: %d %d(%d %d) %d\n", exons[start].start, in, out, k, tc[k].lcnt, tc[k].leftInfo[0][0], tc[k].leftInfo[0][1], tc[k].rcnt ) ;
						covered[ in + out * inCnt ] = true ;
					}
				}
			}
			
			bool isCandidate = true ;
			for ( int in = 0 ; in < inCnt ; ++in )
			{
				int out ;
				for ( out = 0 ; out < outCnt ; ++out )
					if ( covered[ in + out * inCnt ] )
						break ;
				if ( out >= outCnt )
				{
					isCandidate = false ; 
					break ;
				}
			}
			for ( int out = 0 ; out < outCnt ; ++out )
			{
				int in ;
				for ( in = 0 ; in < inCnt ; ++in )
					if ( covered[ in + out * inCnt ] )
						break ;
				if ( in >= inCnt )
				{
					isCandidate = false ;
					break ;
				}
			}

			if ( isCandidate )
			{
				for ( int in = 0 ; in < inCnt ; ++in )
				{
					for ( int out = 0 ; out < outCnt ; ++out )
					{
						if ( !covered[ in + out * inCnt ] )
						{
							tc[ tcCnt ].type = 3 ;
							tc[ tcCnt ].lcnt = 2 ;
							tc[ tcCnt ].leftInfo[0][0] = exons[ nodes[start].prev[ in ] ].start ;
							tc[ tcCnt ].leftInfo[0][1] = exons[ nodes[start].prev[ in ] ].end ;
							tc[ tcCnt ].leftInfo[1][0] = exons[ start ].start ;
							tc[ tcCnt ].leftInfo[1][1] = exons[ start ].end ;

							tc[ tcCnt ].rcnt = 2 ;
							tc[ tcCnt ].rightInfo[0][0] = exons[end].start ;
							tc[ tcCnt ].rightInfo[0][1] = exons[end].end ;
							tc[ tcCnt ].rightInfo[1][0] = exons[ nodes[end].next[out] ].start ;
							tc[ tcCnt ].rightInfo[1][1] = exons[ nodes[end].next[out] ].end ;
						
							tc[tcCnt].strand = exons[start].strand ;
							tc[tcCnt].effectiveCount = 0 ;
							//printf( "Add: %d %d\n", in, out ) ;
							++tcCnt ;
						}
					}
				}
			}

			free( inExons ) ;
			free( outExons ) ;
			free( covered ) ;
		}
	}

	free( tmpPoints ) ;
	free( points ) ;
	return tcCnt ;
}

void DetermineTranscriptConstraintsAbundance( int abdCnt, int tcCnt )
{
	int i, j, k ;
	//int lsupport, rsupport ;
	double labund, rabund ;
	struct _transcriptConstraint *tc = transcriptConstraints ;
	struct _transcriptConstraint *ac = abundanceConstraints ;
	
	//for ( i = 0 ; i < tcCnt ; ++i )
	//{
	//	printf( "	%d %d\n", tc[i].lcnt, tc[i].rcnt ) ;
	//}

	for ( i = 0 ; i < tcCnt ; ++i )	
	{
		if ( tc[i].type == 0 || tc[i].type == 2 || tc[i].effectiveCount < 1 )
		{
			tc[i].support = MAX_READ ;
			tc[i].normAbund = MAX_READ ;
			//printf( "Abundance %d: %d\n", i, tc[i].support) ;
			if ( tc[i].type == 3 )
			{
				tc[i].normAbund = 1e-6 ;
				tc[i].type = 1 ;
			}
			continue ;
		}
		/*if ( !tc[i].lcnt || !tc[i].rcnt )
		{
			//printf( "Abundance %d: %d\n", i, tc[i].support) ;
			continue ;
		}*/
		//lsupport = rsupport = MAX_READ ;
		labund = rabund = 0 ; //MAX_READ ;
		for ( j = 0 ; j < abdCnt ; ++j )
		{
			if ( tc[i].type != ac[j].type )
				continue ;
			if ( tc[i].strand != -1 && ac[j].strand != -1 && tc[i].strand != ac[j].strand )
				continue ;
			if ( tc[i].lcnt != ac[j].lcnt )
				continue ;
			for ( k = 0 ; k < tc[i].lcnt ; ++k )
				if ( tc[i].leftInfo[k][0] != ac[j].leftInfo[k][0]
					|| tc[i].leftInfo[k][1] != ac[j].leftInfo[k][1] )
					break ;
			if ( k >= tc[i].lcnt )
			{
				//lsupport = ac[j].support ;
				labund += ac[j].normAbund ;
				//break ;
			}
		}

		for ( j = 0 ; j < abdCnt ; ++j )
		{
			if ( tc[i].rcnt != ac[j].lcnt )
				continue ;
			if ( tc[i].strand != -1 && ac[j].strand != -1 && tc[i].strand != ac[j].strand )
				continue ;
			for ( k = 0 ; k < tc[i].rcnt ; ++k )
				if ( tc[i].rightInfo[k][0] != ac[j].leftInfo[k][0]
					|| tc[i].rightInfo[k][1] != ac[j].leftInfo[k][1] )
					break ;
			if ( k >= tc[i].rcnt )
			{
				//rsupport = ac[j].support ;
				rabund += ac[j].normAbund ;
				//break ;
			}
		}
		//printf( "Abundance %d: %d <= %d %d\n", i, tc[i].support, lsupport, rsupport ) ;
		//tc[i].support = lsupport < rsupport ? lsupport : rsupport ;
		if ( labund == 0 )
			labund = MAX_READ ;
		if ( rabund == 0 )
			rabund = MAX_READ ;
		tc[i].normAbund = labund < rabund ? labund : rabund ;
		//printf( "DetermineAbundance: %d %lf\n", i, tc[i].normAbund ) ;
	}

	for ( i = 0 ; i < tcCnt ; ++i )
	{
		//tc[i].abundance = tc[i].support ;
		if ( tc[i].type != 1 )
			tc[i].abundance = tc[i].support ;
		else
			tc[i].abundance = tc[i].normAbund ;
	}

}

void EnumerateTranscript( int tag, int visit[], int vcnt, struct _exonNode nodes[], struct _exon exons[], int tcCnt,
	struct _transcript *alltranscripts, int &atcnt ) 
{
	visit[ vcnt ] = tag ;
	++vcnt ;
	int i, j, k ;
	//printf( "## %d %d\n", vcnt, tag ) ;
	/*if ( tag == 351 )
	{
		printf( "###### hi %d %d %d\n", exons[tag].start, exons[tag].end, exons[tag].ncnt ) ;
		exit( 1 ) ;
	}*/	
	/*printf( "###### %d (%d %d):\n", tag, exons[tag].start, exons[tag].end ) ;
	for ( i = 0 ; i < exons[tag].ncnt ; ++i )
		printf( "%d ", exons[tag].next[i] ) ;
	printf( "\n" ) ;*/
	if ( !nodes[tag].ncnt )//&& vcnt > 1 )
	{
		struct _transcript transcript ;
		transcript.ecnt = vcnt ;
		transcript.strand = exons[tag].strand ;
		
		transcript.eid = ( int * )malloc( vcnt * sizeof( int ) ) ;
		transcript.geneId = geneId[tag] ;

		for ( i = 0 ; i < vcnt ; ++i )
		{
			transcript.eid[i] = visit[i] ;
		}
		//transcript.transcriptId = transcriptId ;
		//++transcriptId ;
		
		/*if ( atcnt < MAX_TRANSCRIPT )
		{
			alltranscripts[ atcnt ] = transcript ;
			++atcnt ;
		}*/
		/*if ( exons[tag].start == 1451284 )
		{
			printf( "#### %d %d %d\n", nodes[tag].ncnt, exons[tag].ncnt, vcnt ) ;
			OutputTranscript( "test", exons, &transcript ) ;
			printf( "finish output\n" ) ;
		}*/

		//InsertTranscript( exons, tcCnt, &transcript ) ;
		alltranscripts[ atcnt ] = transcript ;
		++atcnt ;
		
		/*if ( exons[tag].start == 1451284 )
		{
			printf( "##2## %d %d\n", nodes[tag].ncnt, exons[tag].ncnt ) ;
			OutputTranscript( "test", exons, &transcript ) ;
		}*/
	}
	//else if ( !nodes[tag].ncnt )
	//{
	//	return ;
	//}
	else
	{
		// It has next exon, enumerate them.
		for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
			EnumerateTranscript( nodes[tag].next[i], visit, vcnt, 
				nodes, exons, tcCnt, alltranscripts, atcnt ) ;
	}
}


// Search the splice graph and mark all the exons in the gene.
void SearchGeneExons( struct _exonNode nodes[], struct _exon exons[], int tag, int id  ) 
{
	if ( geneId[tag] != -1 )
		return ;
	//printf( "%d %d: %d %d\n", exons[tag].start, exons[tag].end, tag, id ) ;
	geneId[tag] = id ;
	int i ;
	// Find a nice method here.
	if ( exons[tag].end > transcriptEnd )
		transcriptEnd = exons[tag].end ;
	if ( exons[tag].start > transcriptLastExonStart )
		transcriptLastExonStart = exons[tag].start ;
	for ( i = 0 ; i < nodes[tag].pcnt ; ++i )
		SearchGeneExons( nodes, exons, nodes[tag].prev[i], id ) ;
	for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
		SearchGeneExons( nodes, exons, nodes[tag].next[i], id ) ;
		
}

int CoverTheEdges( int tag, struct _exonNode nodes[], int f[], int next[] )
{
	if ( f[tag] != -1 )
		return f[tag] ;	
	int i, tmp ;
	int max = -1, maxtag ;

	if ( nodes[tag].ncnt == 0 )
	{
		next[tag] = -1 ;
		f[tag] = 0 ;
		return 0 ;
	}

	for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
	{
		tmp = CoverTheEdges( nodes[tag].next[i], nodes, f, next ) ;
		if ( tmp + ( nodes[tag].used[i] ? 0 : 1 ) > max )
		{
			max = tmp + ( nodes[tag].used[i] ? 0 : 1 ) ;
			maxtag = i ;
		}
	}
	f[tag] = max ;
	next[tag] = maxtag ;
	return max ;
}

void SearchSubTranscript( int tag, int parents[], int pcnt, struct _dp &pdp, int visit[], int vcnt, 
	int extends[], int extendCnt, const struct _dpAttribute &attr, struct _exonNode nodes[], struct _exon exons[], 
	struct _transcriptConstraint *tc, int tcCnt, int tcStartInd ) 
{
	int i, j, jj, k, cnt, allCnt ;
	double cover = 0 ;
	struct _dp visitdp ;	
	int ecnt ;
	int *eid = attr.bufferEid ;
	int *pairExonTc = attr.pairExonTcIds ;

	visit[ vcnt ] = tag ;
	++vcnt ;
	visitdp.cover = MAX_READ ;
	visitdp.ecnt = 0 ;
	if ( nodes[tag].ncnt && extends != NULL && exons[tag].end >= tc[ extends[ extendCnt - 1 ] ].info[1] )
	{
		// Solve the subtranscript beginning with visit.
		// Now we got the optimal transcript for visit. 
		visitdp = SolveSubTranscript( visit, vcnt, attr, nodes, exons, tc, tcCnt ) ;

		ecnt = pcnt + visitdp.ecnt ;
		//eid = ( int * )malloc( sizeof( int ) * ecnt ) ;
		for ( i = 0 ; i < pcnt ; ++i )
			eid[i] = parents[i] ;
		for ( j = 0 ; j < visitdp.ecnt ; ++j, ++i )
			eid[i] = visitdp.eid[j] ;
		if ( visitdp.cover != -1 )
			free( visitdp.eid ) ;
	}
	else if ( nodes[tag].ncnt && extends != NULL )//&& exons[tag].end < transcriptConstraints[ extends[ extendCnt - 1 ] ].info[1] )
	{
		int tmp ;
		ecnt = pcnt + vcnt ;
		eid = ( int * )malloc( sizeof( int ) * ( ecnt + 1 ) ) ; 
		for ( i = 0 ; i < pcnt ; ++i )
			eid[i] = parents[i] ;
		for ( j = 0 ; j < vcnt ; ++j, ++i )
			eid[i] = visit[j] ;
		for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
		{
			eid[ecnt] = nodes[tag].next[i] ;
			for ( j = extendCnt - 1 ; 
				( tmp = IsConstraintInTranscript( eid, ecnt + 1, tc[ extends[j]], exons ) ) == 0 ||
				 tc[ extends[j] ].info[1] > nodes[ nodes[tag].next[i] ].farthest ; --j )
			{
				;
			}
			/*if ( tmp == 1 && transcriptConstraints[ extends[j] ].type == 1 && transcriptConstraints[ extends[j] ].abundance <= attr.minAbundance )
			{
				if ( pdp.cover == -1 )
				{	
					pdp.cover = 0 ;
					pdp.ecnt = 0 ;
				}
				//printf( "stop early" ) ;
				free( eid ) ;
				return ;
			}*/
			SearchSubTranscript( nodes[tag].next[i], parents, pcnt, pdp, visit, vcnt, extends, j + 1, attr,
				nodes, exons, tc, tcCnt, tcStartInd ) ;
		}
		free( eid ) ;
		return ;
	}
	else //if ( !nodes[tag].ncnt )
	{
		ecnt = pcnt + vcnt ;
		//eid = ( int * )malloc( sizeof( int ) * ecnt ) ;
		for ( i = 0 ; i < pcnt ; ++i )
			eid[i] = parents[i] ;
		for ( j = 0 ; j < vcnt ; ++j, ++i )
			eid[i] = visit[j] ;
		//ecnt = pcnt + vcnt ;
	}

	// Get the constraints between parent and visit.
	if ( !attr.forAbundance )
	{
		double cnt = 0 ;
		int exonTcCnt = 0 ;
		allCnt = 0 ;

		if ( nodes[tag].ncnt )
			k = pcnt - 1 ;
		else
			k = pcnt + vcnt - 1 ;
		// Extract the transcript constraints
		if ( nodes[tag].ncnt )
		{
			for ( i = 0 ; i <= k ; ++i )
			{
				int tag = ( eid[i] - attr.offset ) * attr.size + ( eid[i + 1] - attr.offset ) ;
				for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
				{
					pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
					++exonTcCnt ;
				}
			}
		}
		else
		{
			int tag ;
			for ( i = 0 ; i < k ; ++i )
			{
				tag = ( eid[i] - attr.offset ) * attr.size + ( eid[i + 1] - attr.offset ) ;
				for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
				{
					pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
					++exonTcCnt ;
				}
			}
			tag = ( eid[i] - attr.offset ) * attr.size + eid[i] - attr.offset ;
			for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
			{
				pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
				++exonTcCnt ;
			}
		}


		//for ( j = BinarySearchTranscriptConstraints( exons[ eid[0] ].start, transcriptConstraints, tcCnt ) ; j < visitdp.tcTag ; ++j )
		//for ( j = 0 ; j < visitdp.tcTag ; ++j )
		//for ( j = 0 ; j < tcCnt ; ++j )
		//for ( j = tcStartInd ; j < tcCnt ; ++j )
		for ( jj = 0 ; jj < exonTcCnt ; ++jj )
		{
			j = pairExonTc[jj] ;
			//if ( transcriptConstraints[j].abundance == 0 )
			//	continue ;
			if ( tc[j].leftInfo[0][0] > exons[ eid[k] ].end )
				break ;

			if ( IsConstraintInTranscript( eid, ecnt, tc[j], exons ) == 1 )
			{
				if ( tc[j].abundance > 0 )
				{
					if ( tc[j].type == 1 )
						cnt += tc[j].effectiveCount ;
					else if ( tc[j].type == 2 ) 
					{
						//printf( "%d %d hi\n", eid[0], j ) ;
						if ( !attr.solveRemaining )
							cnt += tc[j].effectiveCount ;
						else
							cnt += ( tc[j].effectiveCount ) / (double)tcCnt ;
					}
				}
				if ( tc[j].normAbund <= attr.minAbundance )
				{
					cnt = -1 ;
					break ;
				}

				if ( tc[j].type == 1 || tc[j].type == 2 )
					allCnt += tc[j].effectiveCount ;
			}
		}

		// The reason why I shift is to differentiate cover no constraints and invalid.
		if ( nodes[tag].ncnt == 0 )
		{
			++cnt ; // shift it
		}
		else
		{
			if ( cnt == -1 || visitdp.cover == 0 )
				cnt = 0 ;
			else
				cnt = cnt + visitdp.cover ;
		}

		cover = cnt ;
		/*if ( transcriptId == 125 && attr.minAbundance == 9 )
		  {
		  printf( "process: %d (%d %d) %d %d\n\t", cnt, j, tcCnt, transcriptConstraints[j].support, attr.minAbundance ) ;
		  for ( i = 0 ; i < ecnt ; ++i )
		  printf( "%d ", eid[i] ) ;
		  printf( "\n" ) ;
		  }*/
	}
	else
	{
		int exonTcCnt = 0 ;
		cover = visitdp.cover ;
		//cover = MAX_READ ;
		if ( nodes[tag].ncnt )
			k = pcnt - 1 ;
		else
			k = pcnt + vcnt - 1 ;
		
		if ( nodes[tag].ncnt )
		{
			for ( i = 0 ; i <= k ; ++i )
			{
				int tag = ( eid[i] - attr.offset ) * attr.size + ( eid[i + 1] - attr.offset ) ;
				for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
				{
					pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
					++exonTcCnt ;
				}
			}
		}
		else
		{
			int tag ;
			for ( i = 0 ; i < k ; ++i )
			{
				tag = ( eid[i] - attr.offset ) * attr.size + ( eid[i + 1] - attr.offset ) ;
				for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
				{
					pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
					++exonTcCnt ;
				}
			}
			tag = ( eid[i] - attr.offset ) * attr.size + eid[i] - attr.offset ;
			for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
			{
				pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
				++exonTcCnt ;
			}
		}

		//for ( j = BinarySearchTranscriptConstraints( exons[ eid[0] ].start, transcriptConstraints, tcCnt ) ; j < visitdp.tcTag ; ++j )
		//for ( j = 0 ; j < visitdp.tcTag ; ++j )
		//for ( j = 0 ; j < tcCnt ; ++j )
		//for ( j = tcStartInd ; j < tcCnt ; ++j )
		for ( jj = 0 ; jj < exonTcCnt ; ++jj )
		{
			j = pairExonTc[jj] ;
			//if ( transcriptConstraints[j].abundance == 0 )
			//	continue ;
			if ( tc[j].leftInfo[0][0] > exons[ eid[k] ].end )
				break ;
			/*if ( transcriptConstraints[j].abundance == 375 )
			  {
			  printf( "%d %d: %d\n", 375, transcriptConstraints[j].info[1], tag ) ;
			  }*/
			if ( tc[j].type == 1 && IsConstraintInTranscript( eid, ecnt, tc[j], exons ) == 1 &&
					tc[j].normAbund < cover  )
			{
				cover = tc[j].normAbund ;
			}
		}
		if ( cover == MAX_READ )
			cover = 1e-6 ;
		//cover = cnt ;
	}
	
	//printf( "cnt = %d\n", cnt ) ;
	// Compare with the optimal transcript of parent.
	if ( cover > pdp.cover ) 
	{
		//printf( "%d\n", exons[tag].start ) ;
		if ( pdp.cover != -1 )
			free( pdp.eid ) ;
		//pdp.eid = eid ;
		pdp.eid = (int *)malloc( sizeof( int ) * ecnt ) ;
		memcpy( pdp.eid, eid, sizeof( int ) * ecnt ) ;
		pdp.ecnt = ecnt ;
		//pdp.cover = cnt ;
		pdp.cover = cover ;
	}
	else if ( cover == pdp.cover )
	{
		int tmp = 0 ;
		if ( pdp.cover == -1 )
			pdp.ecnt = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
			if ( ( tc[j].type == 1 || tc[j].type == 2 ) &&
				IsConstraintInTranscript( eid, ecnt, tc[j], exons ) == 1 )
			{
				++tmp ;
			}

		if ( allCnt > tmp )
		{
			if ( pdp.cover != -1 )
				free( pdp.eid ) ;
			//pdp.eid = eid ;
			pdp.eid = (int *)malloc( sizeof( int ) * ecnt ) ;
			memcpy( pdp.eid, eid, sizeof( int ) * ecnt ) ;
			pdp.ecnt = ecnt ;
			pdp.cover = cover ;
		}
	}
	/*else
	{
		//if ( visitdp.cover != -1 )
		free( eid ) ;
	}*/
}

int CompTcEnd( const void *p1, const void *p2 )
{
	//const struct _transcriptConstraint &tc1 = transcriptConstraints[*( int *)p1] ;
	//const struct _transcriptConstraint &tc2 = transcriptConstraints[*( int *)p2] ;
	//return tc1.info[1] - tc2.info[1] ;
	//return *( int *)p1 - *( int * )p2 ;
	return ( ( struct _pair *)p1 )->a - ( ( struct _pair *)p2 )->a ;
}

struct _dp SolveSubTranscript( int visit[], int vcnt, const struct _dpAttribute &attr, struct _exonNode nodes[], struct _exon exons[], 
	struct _transcriptConstraint *tc, int tcCnt )
{
	int i, j, jj, k ;
	//struct _transcriptConstraint *tc = transcriptConstraints ;
	int cnt = 0, tmp ; 
	//int newend = exons[ visit[ vcnt - 1 ] ].end ;
	int vcover ;
	struct _dp visitDp, ret ;
	int tag ;
	int tcStartInd ;
	bool belowMin ; 
	int *pairExonTc = attr.pairExonTcIds ;
	//if ( attr.forAbundance == true )
	//	printf( "" ) ;

	++(*attr.solveCaller) ;
	visitDp.ecnt = ret.ecnt = 0 ;
	visitDp.eid = ret.eid = NULL ;
	
	if ( vcnt == 1 && !nodes[visit[0]].ncnt && !nodes[visit[0]].pcnt )
	{
		// Single-exon transcript
		visitDp.cover = -1 ;
		visitDp.ecnt = 0 ;
		int nextv[2] ;
		SearchSubTranscript( visit[0], visit, 0, visitDp, nextv, 0, NULL, 0, 
			attr, nodes, exons, tc, tcCnt, 0 ) ;
		return visitDp ;
	}
	// If the visit sequence is stored.
	if ( vcnt == 1 )
	{
		int a = visit[0] - attr.offset ;
		if ( attr.f1[a].cover != -1 && ( attr.f1[a].timeStamp == attr.timeStamp || 
			( attr.f1[a].minAbundance <= attr.minAbundance && attr.f1[a].cover == 0 ) ) )
		{
			//printf( "works1\n" ) ;
			ret = attr.f1[a] ;
			ret.eid = (int *)malloc( sizeof( int ) * attr.f1[a].ecnt  ) ;
			memcpy( ret.eid, attr.f1[a].eid, sizeof( int ) * attr.f1[a].ecnt ) ;
			return ret ;
		}
	}
	else if ( vcnt == 2 )
	{
		int a = visit[0] - attr.offset ;
		int b = visit[1] - attr.offset ;
		if ( attr.f2[a][b].cover != -1 && ( attr.f2[a][b].timeStamp == attr.timeStamp || 
			( attr.f2[a][b].minAbundance <= attr.minAbundance && attr.f2[a][b].cover == 0 ) ) )
		{
			//printf( "works2\n" ) ;
			ret = attr.f2[a][b] ;
			ret.eid = (int *)malloc( sizeof( int ) * attr.f2[a][b].ecnt  ) ;
			memcpy( ret.eid, attr.f2[a][b].eid, sizeof( int ) * attr.f2[a][b].ecnt ) ;
			return ret ;
		}

	}
	else if ( vcnt == 3 && attr.size <= USE_F3 )
	{
		int a = visit[0] - attr.offset ;
		int b = visit[1] - attr.offset ;
		int c = visit[2] - attr.offset ;

		if ( attr.f3[a][b][c].cover != -1 && ( attr.f3[a][b][c].timeStamp == attr.timeStamp || 
			( attr.f3[a][b][c].minAbundance <= attr.minAbundance && attr.f3[a][b][c].cover == 0 ) ) )
		{
			//printf( "works3\n" ) ;
			ret = attr.f3[a][b][c] ;
			ret.eid = (int *)malloc( sizeof( int ) * attr.f3[a][b][c].ecnt  ) ;
			memcpy( ret.eid, attr.f3[a][b][c].eid, sizeof( int ) * attr.f3[a][b][c].ecnt ) ;
			return ret ;
		}
	}
	else /*if ( vcnt <= 8 )*/
	{
		int maxStride = -1 ;
		for ( i = 1 ; i < vcnt ; ++i )
		{
			if ( visit[i] - visit[i - 1] > maxStride )
				maxStride = visit[i] - visit[i - 1] ;
		}

		if ( maxStride <= DP_EXON_ID_STRIDE && vcnt <= DP_STRIDE_CNT + 1 )
		{
			int strides[DP_STRIDE_CNT] ;
			int key = 0 ;
			int a = visit[0] - attr.offset ; 
			for ( i = 0 ; i < vcnt - 1 ; ++i )
				strides[i] = visit[i + 1] - visit[i] ; 
			for ( ; i < DP_STRIDE_CNT ; ++i )
				strides[i] = 0 ;
			for ( i = 0 ; i < DP_STRIDE_CNT ; ++i )
				key = key * ( DP_EXON_ID_STRIDE + 1 ) + strides[i] ;
			
			if ( attr.stride[a].p[key].cover != -1 && ( attr.stride[a].p[key].timeStamp == attr.timeStamp ||
				( attr.stride[a].p[key].minAbundance <= attr.minAbundance && attr.stride[a].p[key].cover == 0 ) ) )
			{
				//printf( "stride works\n" ) ;
				ret = attr.stride[a].p[key] ;
				ret.eid = ( int * )malloc( sizeof( int ) * attr.stride[a].p[key].ecnt ) ;
				memcpy( ret.eid, attr.stride[a].p[key].eid, sizeof( int ) * attr.stride[a].p[key].ecnt ) ;
				return ret ;
			}
		}

		int key = 0 ;
		for ( i = 0 ; i < vcnt ; ++i )
		{
			key = ( key * 17 + ( visit[i] - attr.offset ) ) % MAX_TRANSCRIPT ;
		}

		if ( attr.hash[key].cover != -1 && attr.hash[key].cnt == vcnt && 
			( attr.hash[key].timeStamp == attr.timeStamp || 
				( attr.hash[key].minAbundance <= attr.minAbundance && attr.hash[key].cover == 0 ) ) )
		{
			for ( i = 0 ; i < vcnt ; ++i )
			{
				if ( attr.hash[key].eid[i] != visit[i] )
					break ;
			}

			if ( i >= vcnt )
			{
				//printf( "works hash\n" ) ;
				ret.eid = (int *)malloc( sizeof( int ) * attr.hash[key].ecnt ) ;
				memcpy( ret.eid, attr.hash[key].eid, sizeof( int ) * attr.hash[key].ecnt ) ;
				ret.ecnt = attr.hash[key].ecnt ;
				ret.cover = attr.hash[key].cover ;
				return ret ;
			}
		}
	}
	

	int *nextv = ( int * )malloc( sizeof( int ) * attr.size ) ;
	int *extends = ( int *)malloc( sizeof( int ) * ( tcCnt + 1 ) ) ; // The possible compatible transcript constraints ids sorted by their ends.
	int exonTcCnt = 0 ;
	/*printf( "%d: ", vcnt ) ;
	for ( i = 0 ; i < vcnt ; ++i )
		printf( "%d ", visit[i] ) ;
	printf( "\n" ) ;*/
	//extends[0] = exons[ visit[ vcnt - 1 ] ].end ; 
	k = 0 ;
	//tcStartInd = BinarySearchTranscriptConstraints( exons[ visit[0] ].start, tc, tcCnt ) ;
	belowMin = false ;
	
	for ( i = 0 ; i < vcnt - 1 ; ++i )
	{
		tag = ( visit[i] - attr.offset ) * attr.size + ( visit[i + 1] - attr.offset ) ;
		for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
		{
			pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
			++exonTcCnt ;
		}
	}
	tag = ( visit[i] - attr.offset ) * attr.size + visit[i] - attr.offset ;
	for ( j = 0 ; j < attr.pairExonTc[tag].cnt ; ++j )
	{
		pairExonTc[ exonTcCnt ] = attr.pairExonTc[tag].id[j] ;
		++exonTcCnt ;
	}
	//for ( j = tcStartInd ; j < tcCnt ; ++j )
		//for ( j = 0 ; j < tcCnt ; ++j )
	for ( jj = 0 ; jj < exonTcCnt ; ++jj )
	{
		j = pairExonTc[jj] ;

		if ( tc[j].leftInfo[0][0] > exons[ visit[ vcnt - 1 ] ].end )
			break ;

		if ( tc[j].type == 0 && tc[j].info[0] == visit[ vcnt - 1 ] )
		{
			extends[k] = j ;
			++k ;
			continue ;
		}
		if ( tc[j].normAbund <= attr.minAbundance || tc[j].abundance == 0 || tc[j].type == 0 )//( tcUsed[j] )
			continue ;


		tmp = IsConstraintInTranscript( visit, vcnt, tc[j], exons ) ;

		if ( tmp == 1 && tc[j].normAbund <= attr.minAbundance )
		{
			belowMin = true ;
			break ;
		}


		if ( tmp == 2 )
		{
			/*int a = -1, b = -1, c;
			  if ( tc[j].lcnt )
			  a = tc[j].leftInfo[ tc[j].lcnt - 1][1] ;			
			  if ( tc[j].rcnt )
			  b = tc[j].rightInfo[ tc[j].rcnt - 1 ][1] ;
			  c = a > b ? a : b ;*/
			/*if ( vcnt < 30 )
			  {
			  for ( i = 0 ; i < vcnt ; ++i )
			  printf( "(%d %d) ", exons[ visit[i] ].start, exons[ visit[i] ].end ) ;
			  printf( "| " ) ;
			  for ( i = 0 ; i < tc[j].lcnt ; ++i )
			  printf( "[%d %d] ", tc[j].leftInfo[i][0], tc[j].leftInfo[i][1] ) ;
			  printf( "| " ) ;
			  for ( i = 0 ; i < tc[j].rcnt ; ++i )
			  printf( "[%d %d] ", tc[j].rightInfo[i][0], tc[j].rightInfo[i][1] ) ;
			  printf( "|| %d\n", tc[j].info[1] ) ;
			  }*/
			extends[k] = j ;
			++k ;	
		}
	}
	//printf( "Extends count: %d\n", k ) ;

	if ( k > 0 )
	{
		// Sort the extends based on their ends
		struct _pair *extendsPair = (struct _pair *)malloc( sizeof( *extendsPair ) * k ) ;
		for ( i = 0 ; i < k ; ++i )
		{
			extendsPair[i].a = tc[ extends[i] ].info[1] ;
			extendsPair[i].b = extends[i] ;
		}
		qsort( extendsPair, k, sizeof( struct _pair ), CompTcEnd ) ;
		for ( i = 0 ; i < k ; ++i )
			extends[i] = extendsPair[i].b ;

		free( extendsPair ) ;
	}
	/*printf( "##\n" ) ;
	  for ( i = 0 ; i < k ; ++i )
	  printf( "%d\n", tc[ extends[i] ].type ) ;
	  printf( "==\n" ) ;*/

	//printf( ": %d\n", exons[ visit[ vcnt - 1 ] ].end ) ;
	//printf( "## %d %d\n", exons[ visit[0] ].start, newend ) ;
	visitDp.cover = -1 ;
	visitDp.ecnt = 0 ;

	if ( !belowMin )
	{
		tag = visit[ vcnt - 1 ] ;
		for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
		{
			visit[vcnt] = nodes[tag].next[i] ;
			for ( j = k - 1 ; 
					IsConstraintInTranscript( visit, vcnt + 1, tc[extends[j]], exons ) == 0 ||
					tc[ extends[j] ].info[1] > nodes[ nodes[tag].next[i] ].farthest ; --j )
				;
			//printf( "hi %d %d %d: %d\n", visit[ vcnt - 1 ], tc[ extends[0] ].type, tc[ extends[0] ].info[0], j ) ;
			//printf( "## %d => %d\n", newend, extends[j] ) ;
			//printf( "### %d\n", j ) ;
			//if ( exons[tag].end > 49390571 )
			//	printf( "	%d: %d %d %d\n", vcnt, exons[ visit[0] ].end, exons[ nodes[tag].next[i] ].end, tc[ extends[j] ].info[1] ) ;

			SearchSubTranscript( nodes[tag].next[i], visit, vcnt, visitDp, nextv, 0, extends, j + 1, 
					attr, nodes, exons, tc, tcCnt, tcStartInd ) ;	
		}
	}
	else
	{
		visitDp.cover = 0 ;
		visitDp.eid = (int *)malloc( 0 ) ;
	}
	
	free( nextv ) ;
	free( extends ) ;
	// Memorize the result
	ret = visitDp ;
	if ( vcnt == 1 )
	{	
		int a = visit[0] - attr.offset ;
		if ( attr.f1[a].timeStamp != 0 )
			free( attr.f1[a].eid ) ;
		attr.f1[a].eid = (int *)malloc( sizeof( int ) * ret.ecnt  ) ;
		memcpy( attr.f1[a].eid, ret.eid, sizeof( int ) * ret.ecnt ) ;
		attr.f1[a].ecnt = ret.ecnt ;
		attr.f1[a].cover = ret.cover ;
		attr.f1[a].minAbundance = attr.minAbundance ;
		attr.f1[a].timeStamp = attr.timeStamp ;
	}
	else if ( vcnt == 2 )
	{
		int a = visit[0] - attr.offset ;
		int b = visit[1] - attr.offset ;

		if ( attr.f2[a][b].timeStamp != 0 )
			free( attr.f2[a][b].eid ) ;
		
		attr.f2[a][b].eid = (int *)malloc( sizeof( int ) * ret.ecnt  ) ;
		memcpy( attr.f2[a][b].eid, ret.eid, sizeof( int ) * ret.ecnt ) ;
		attr.f2[a][b].ecnt = ret.ecnt ;
		attr.f2[a][b].cover = ret.cover ;
		attr.f2[a][b].minAbundance = attr.minAbundance ;
		attr.f2[a][b].timeStamp = attr.timeStamp ;
	}
	else if ( vcnt == 3 && attr.size <= USE_F3 )
	{
		int a = visit[0] - attr.offset ;
		int b = visit[1] - attr.offset ;
		int c = visit[2] - attr.offset ;

		if ( attr.f3[a][b][c].timeStamp != 0 )
			free( attr.f3[a][b][c].eid ) ;
		attr.f3[a][b][c].eid = (int *)malloc( sizeof( int ) * ret.ecnt  ) ;
		memcpy( attr.f3[a][b][c].eid, ret.eid, sizeof( int ) * ret.ecnt ) ;
		attr.f3[a][b][c].ecnt = ret.ecnt ;
		attr.f3[a][b][c].cover = ret.cover ;
		attr.f3[a][b][c].minAbundance = attr.minAbundance ;
		attr.f3[a][b][c].timeStamp = attr.timeStamp ;
	}
	else /*if ( vcnt <= 8 )*/
	{
		int maxStride = -1 ;
		for ( i = 1 ; i < vcnt ; ++i )
		{
			if ( visit[i] - visit[i - 1] > maxStride )
				maxStride = visit[i] - visit[i - 1] ;
		}

		if ( maxStride <= DP_EXON_ID_STRIDE && vcnt <= DP_STRIDE_CNT + 1 )
		{
			// use stride
			int strides[DP_STRIDE_CNT] ;
			int key = 0 ;
			int a = visit[0] - attr.offset ; 
			for ( i = 0 ; i < vcnt - 1 ; ++i )
				strides[i] = visit[i + 1] - visit[i] ; 
			for ( ; i < DP_STRIDE_CNT ; ++i )
				strides[i] = 0 ;
			for ( i = 0 ; i < DP_STRIDE_CNT ; ++i )
				key = key * ( DP_EXON_ID_STRIDE + 1 ) + strides[i] ;

			if (  attr.stride[a].p[key].timeStamp != 0 )
				free( attr.stride[a].p[key].eid ) ;
			attr.stride[a].p[key].eid = ( int * )malloc( sizeof( int ) * ret.ecnt ) ;
			memcpy( attr.stride[a].p[key].eid, ret.eid, sizeof( int ) * ret.ecnt ) ;
			attr.stride[a].p[key].ecnt = ret.ecnt ;
			attr.stride[a].p[key].cover = ret.cover ;
			attr.stride[a].p[key].minAbundance = attr.minAbundance ;
			attr.stride[a].p[key].timeStamp = attr.timeStamp ;
		}
		else
		{
			// try hash
			int key = 0 ;
			bool update = false ;
			for ( i = 0 ; i < vcnt ; ++i )
			{
				key = ( key * 17 + ( visit[i] - attr.offset ) ) % MAX_TRANSCRIPT ;
			}
			
			if ( attr.hash[key].timeStamp != attr.timeStamp || attr.hash[key].cnt != vcnt )
				update = true ;
			else
			{
				for ( i = 0 ; i < vcnt ; ++i )
					if ( visit[i] != attr.hash[key].eid[i] )
						break ;
				if ( i < vcnt )
					update= true ;
			}

			if ( update ) //attr.hash[key].cover == -1  )
			{
				if ( attr.hash[key].timeStamp != 0 )
					free( attr.hash[key].eid ) ;

				attr.hash[key].eid = ( int * )malloc( sizeof( int ) * ret.ecnt ) ;
				memcpy( attr.hash[key].eid, ret.eid, sizeof( int ) * ret.ecnt ) ;
				attr.hash[key].cnt = vcnt ;
				attr.hash[key].ecnt = ret.ecnt ;
				attr.hash[key].cover = ret.cover ;
				attr.hash[key].minAbundance = attr.minAbundance ;
				attr.hash[key].timeStamp = attr.timeStamp ;
			}
		}
	}
	return ret ;
}

void CleanDpMemo( const struct _dpAttribute &attr )
{
	int i, j, k ;
	for ( i = 0 ; i < attr.size ; ++i )
	{
		if ( attr.f1[i].cover != -1 )
		{
			free( attr.f1[i].eid ) ;
			attr.f1[i].cover = -1 ;
		}
		for ( j = 0 ; j < attr.size ; ++j )
		{
			if ( attr.f2[i][j].cover != -1 )
			{
				free( attr.f2[i][j].eid ) ;
				attr.f2[i][j].cover = -1 ;
			}
			if ( attr.size <= USE_F3 )
			{
				for ( k = 0 ; k < attr.size ; ++k )
				{
					if ( attr.f3[i][j][k].cover != -1 )
					{
						free( attr.f3[i][j][k].eid ) ;
						attr.f3[i][j][k].cover = -1 ;
					}
				}
			}
		}

		for ( j = 0 ; j < SIZEOF_STRIDE ; ++j )
			if ( attr.stride[i].p[j].cover != -1 )
			{
				free( attr.stride[i].p[j].eid ) ;
				attr.stride[i].p[j].cover = -1 ;
			}
	}
	for ( i = 0 ; i < MAX_TRANSCRIPT ; ++i )
		if ( attr.hash[i].cover != -1 )
		{
			free( attr.hash[i].eid ) ;
			attr.hash[i].cover = -1 ;
		}
}

void PickTranscriptsByDP( char *chrom, int from, int to, struct _exonNode nodes[], struct _exon exons[], 
	struct _transcript *transcripts, int &tcnt,
	struct _transcriptConstraint *tc, const int &tcCnt )
{
	int i, j, k ;
	struct _dp best, tmp, tmp2 ;
	struct _transcript transcript ;
	struct _dpAttribute attr ;
	//struct _transcriptConstraint *tc = transcriptConstraints ;
	int head[5] ;
	double maxAbundance = -1 ;
	tcnt = 0 ;

	//memset( tcUsed, false, sizeof( bool ) * tcCnt ) ;
	attr.size = to - from + 1 ;
	attr.offset = from ;
	attr.f1 = ( struct _dp * )malloc( sizeof( struct _dp ) * ( to - from + 1 ) ) ;

	attr.f2 = ( struct _dp ** )malloc( sizeof( struct _dp * ) * ( to - from + 1 ) ) ;
	for ( i = 0 ; i < to - from + 1 ; ++i )
		attr.f2[i] = ( struct _dp * )malloc( sizeof( struct _dp ) * ( to - from + 1 ) ) ;
	if ( attr.size <= USE_F3 )
	{
		attr.f3 = ( struct _dp *** )malloc( sizeof( struct _dp ** ) * ( to - from + 1 ) ) ;	
		for ( i = 0 ; i < to - from + 1 ; ++i )
			attr.f3[i] = ( struct _dp ** )malloc( sizeof( struct _dp* ) * ( to - from + 1 ) ) ;
		for ( i = 0 ; i < to - from + 1 ; ++i )
		{
			for ( j = 0 ; j < to - from + 1 ; ++j )
				attr.f3[i][j] = ( struct _dp * )malloc( sizeof( struct _dp ) * ( to - from + 1 ) ) ;
		}
	}

	attr.stride = ( struct _dpStride * )malloc( sizeof( struct _dpStride ) * ( to - from + 1 ) ) ;
	attr.hash = ( struct _dpHash * )malloc( sizeof( struct _dpHash ) * MAX_TRANSCRIPT ) ;
	
	for ( i = 0 ; i < to - from + 1 ; ++i )
	{
		attr.f1[i].cover = -1 ;
		attr.f1[i].timeStamp = 0 ;
		for ( j = 0 ; j < to - from + 1 ; ++j )
		{
			attr.f2[i][j].cover = -1 ;
			attr.f2[i][j].timeStamp = 0 ;
			if ( attr.size <= USE_F3 )
				for ( k = 0 ; k < to - from + 1 ; ++k )
				{
					attr.f3[i][j][k].cover = -1 ;
					attr.f3[i][j][k].timeStamp = 0 ;
				}
		}

		for ( j = 0 ; j < SIZEOF_STRIDE ; ++j )
		{
			attr.stride[i].p[j].cover = -1 ;
			attr.stride[i].p[j].timeStamp = 0 ;
		}
	}
	for ( i = 0 ; i < MAX_TRANSCRIPT ; ++i )
	{
		attr.hash[i].cover = -1 ;
		attr.hash[i].timeStamp = 0 ;
	}
	attr.solveRemaining = false ;
	
	// Get the constraints for paired exon
	attr.pairExonTc = ( struct _tcList * )malloc( sizeof( struct _tcList ) * 
			( to - from + 1 ) * ( to - from + 1 ) ) ;
	attr.pairExonTcIds = ( int * )malloc( sizeof( int ) * tcCnt ) ;
	attr.bufferEid = (int * )malloc( sizeof( int ) * ( to - from + 1 ) ) ;

	int *pairExonTc = attr.pairExonTcIds ;
	memset( attr.pairExonTc, 0, sizeof( struct _tcList ) * ( to - from + 1 ) * ( to - from + 1 ) ) ;
	for ( i = from ; i <= to ; ++i )
	{
		int cnt ;
		int startTcInd = BinarySearchTranscriptConstraints( exons[i].start, tc, tcCnt ) ;
		for ( j = 0 ; j < nodes[i].ncnt ; ++j )
		{
			cnt = 0 ;
			for ( k = startTcInd ; k < tcCnt ; ++k )
			{
				int eid[3] ;
				if ( tc[k].leftInfo[0][0] > exons[i].end )
					break ;
				eid[0] = i ;
				eid[1] = nodes[i].next[j] ;
				if ( IsConstraintInTranscript( eid, 2, tc[k], exons ) != 0 )
				{
					pairExonTc[cnt] = k ;
					++cnt ;
				}
			}

			if ( cnt > 0 )
			{
				int tag = ( i - from ) * ( to - from + 1 ) + nodes[i].next[j] - from ;
				attr.pairExonTc[tag].cnt = cnt ;
				attr.pairExonTc[tag].id = ( int * )malloc( sizeof( int ) * cnt ) ;
				memcpy( attr.pairExonTc[tag].id, pairExonTc, sizeof( int ) * cnt ) ;
			}
		}
		// use the diagonal element to store the single-exon compatible constraint
		cnt = 0 ;
		for ( k = startTcInd ; k < tcCnt ; ++k )
		{
			int eid[3] ;
			if ( tc[k].leftInfo[0][0] > exons[i].end )
				break ;
			eid[0] = i ;
			if ( IsConstraintInTranscript( eid, 1, tc[k], exons ) != 0 )
			{
				pairExonTc[cnt] = k ;
				++cnt ;
			}
		}

		if ( cnt > 0 )
		{
			int tag = ( i - from ) * ( to - from + 1 ) + i - from ;
			attr.pairExonTc[tag].cnt = cnt ;
			attr.pairExonTc[tag].id = ( int * )malloc( sizeof( int ) * cnt ) ;
			memcpy( attr.pairExonTc[tag].id, pairExonTc, sizeof( int ) * cnt ) ;
		}
	}

	// Calculate the max abundance.
	attr.forAbundance = true ;
	attr.minAbundance = 0 ;
	attr.timeStamp = 1 ;
	maxAbundance = -1 ;
	int solveCaller = 0 ;
	attr.solveCaller = &solveCaller ;
	//printf( "before find maxAbundance\n" ) ;
	for ( i = from ; i <= to ; ++i )
	{
		if ( nodes[i].pcnt )//|| !nodes[i].ncnt )
			continue ;
		head[0] = i ;
		//printf( "hi1: %d\n", i ) ;
		tmp = SolveSubTranscript( head, 1, attr, nodes, exons, tc, tcCnt ) ;
		//printf( "hi2: %d\n", i ) ;
		if ( tmp.cover > maxAbundance ) // Actually, it is abundance
			maxAbundance = tmp.cover ;
		if ( tmp.cover != -1 )
			free( tmp.eid ) ;
	}
	//printf( "after find maxAbundance: %lf %d\n", maxAbundance, solveCaller ) ;
	attr.forAbundance = false ;
	//CleanDpMemo( attr ) ;
	++attr.timeStamp ;

	//printf( "### %d\n", attr.size) ;
	//printf( "%lf\n", maxAbundance ) ;
	//for ( i = 0 ; i < tcCnt ; ++i )
	//	printf( "%d ", tc[i].support ) ;
	//printf( "\n" ) ;
	//printf( "%lf %d\n", maxAbundance, tcCnt ) ;
	//exit( 1 ) ;		
	while ( 1 )
	{
		attr.minAbundance = 0 ;
		//printf( "## %lf\n", maxAbundance ) ;
		// Calculate the transcript with most abundance.
		best.cover = 0 ;
		double bestScore = 0 ;
		double bestAbundance ;
		double cover ;

		while ( 1 ) 
		{
			tmp2.cover = 0 ;
			for ( i = from ; i <= to ; ++i )
			{
				if ( nodes[i].pcnt )// || !nodes[i].ncnt )
					continue ;
				head[0] = i ;
				//printf( "hi %lf: from exon %d\n", attr.minAbundance, i ) ;
				tmp = SolveSubTranscript( head, 1, attr, nodes, exons, tc, tcCnt ) ;
				--tmp.cover ;	
				//printf( "hi1 %d %lf\n", tmp.ecnt, tmp.cover ) ;
				if ( tmp.cover > tmp2.cover )
				{
					if ( tmp2.cover != 0 )
						free( tmp2.eid ) ;
					tmp2 = tmp ;
				}
				else
				{
					//if ( tmp.ecnt != 0 )
						free( tmp.eid ) ; 
				}
			}

			// Get the non-evidence constraints
			if ( evidences != NULL && tmp2.cover != 0 )
			{
				cover = 0 ;
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( tc[j].type == 2 )
						continue ;
					if ( tc[j].type == 1 && tc[j].abundance > 0 && IsConstraintInTranscript( tmp2.eid, tmp2.ecnt, tc[j], exons ) == 1 )
						++cover ;
				}
			}
			else
				cover = tmp2.cover ;
			if ( cover == 0 || ComputeScore( tmp2.cover, maxAbundance, maxAbundance ) < bestScore ) //tmp2.cover / double(2 - 1.0) < bestScore )
				break ;
			// Get the abundance of tmp2.
			double abundance = MAX_READ ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				//if ( tc[j].abundance == 0 )
				//	continue ;
				if ( tc[j].normAbund < abundance && IsConstraintInTranscript( tmp2.eid, tmp2.ecnt, tc[j], exons ) == 1 )
					abundance = tc[j].normAbund ;
			}
			bool onlyCoverExons = false  ;
			if ( abundance == MAX_READ )
			{
				abundance = 1e-6 ;
				onlyCoverExons = true ;
			}
			//printf( "%d (%lf => %lf)\n", abundance, tmp2.cover / (double)( 2 - abundance / maxAbundance ), bestScore ) ;
			//printf( "	%lf: %d %d %lf\n", tmp2.cover / (double)( 2 - abundance / maxAbundance ), tmp2.cover, abundance, maxAbundance ) ;
			if ( ComputeScore( tmp2.cover, abundance, maxAbundance) > bestScore )
			{
				bestScore = ComputeScore( tmp2.cover, abundance, maxAbundance ) ;
				if ( best.cover != 0 )
					free( best.eid ) ;
				best = tmp2 ;
				bestAbundance = abundance ;
			}
			else
			{
				//if ( tmp2.ecnt != 0 )
					free( tmp2.eid ) ;
			}
			//printf( "	%d %d\n", attr.minAbundance, abundance ) ;	
			if ( onlyCoverExons )
				break ;
			attr.minAbundance = abundance ;
			//CleanDpMemo( attr ) ;		
			++attr.timeStamp ;
			if ( attr.timeStamp == 0 )
				++attr.timeStamp ;
		} //end of while tmp2
		/*double ttt = 0 ;
		if ( evidences != NULL )
		{
			cover = 0 ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				//if ( tc[j].type == 2 )
				//
				//	continue ;
				if ( tc[j].abundance > 0 && IsConstraintInTranscript( best.eid, best.ecnt, tc[j], exons ) == 1 )
				{
					if ( tc[j].type != 2 )
						++cover ;
					else		
						++ttt ;
				}
			}
		}
		else
			cover = best.cover ;*/
		
		if ( best.cover == 0 && ( attr.solveRemaining || evidences == NULL ) )
			break ;
		else if ( best.cover == 0 )
		{
			attr.solveRemaining = true ;
			continue ;
		}
		transcripts[tcnt].ecnt = best.ecnt ;
		transcripts[tcnt].eid = best.eid ;
		//transcripts[tcnt].eid = ( int * )malloc( sizeof( int ) * best.ecnt ) ;
		//memcpy( transcripts[tcnt].eid, best.eid, sizeof( int ) * best.ecnt ) ;
		//transcripts.transcriptId = transcriptId ;
		//++transcriptId ;
		transcripts[tcnt].geneId = geneId[ best.eid[0] ] ;
		transcripts[tcnt].abundance = 0 ;
		for ( i = 0 ; i < transcripts[tcnt].ecnt ; ++i )
			++nodes[ transcripts[tcnt].eid[i] ].weight ;
		//OutputTranscript( chrom, exons, &transcript ) ;
		
		/*printf( "%lf %lf\n",  best.cover, bestScore ) ;
		printf( "\t" ) ;
		for ( i = 0 ; i < best.ecnt ; ++i )
			printf( "%d ", best.eid[i] ) ;
		printf( "\n" ) ;*/

		/*if ( transcriptId == 26 )
		  {
		  int eid[17]={76,78,79, 80, 81, 82, 83, 85, 87, 88, 89, 91, 92, 93, 94, 95, 96}; 
		  double cnt = 0 ;
		  for ( j = 0 ; j < tcCnt; ++j )
		  {
		  if ( tc[j].abundance > 0 && IsConstraintInTranscript( eid, 17, tc[j], exons ) == 1 )
		  cnt += 1 / (double)( tcUsedCnt[j] + 1 ) ;
		  }
		  printf( "26: %lf => ", cnt ) ;
		  }*/
		double abundance = MAX_READ ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( IsConstraintInTranscript( best.eid, best.ecnt, tc[j], exons ) == 1 )
			{
				if ( tc[j].abundance > 0 && tc[j].abundance < abundance )
					abundance = tc[j].abundance ;
				//++tcUsedCnt[j] ;
			}
		}
		if ( abundance == MAX_READ )
			abundance = 1e-6 ;

		// Remove the constraints it covers
		double factor = 1 ;
		if ( tcCnt > 3000 )
			factor *= ( ( tcCnt - 3000.0 ) / 1000.0 + 1 ) ;
		if ( solveCaller > 5000 )
			factor *= ( ( factor - 5000.0 ) / 2500.0 ) ;
		factor = factor * factor ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( IsConstraintInTranscript( best.eid, best.ecnt, tc[j], exons ) == 1 && tc[j].abundance > 0 )
			{
				if ( tc[j].type == 0 )
					tc[j].abundance = 0 ;
				else if ( tc[j].type == 1 )
				{
					double time2 = 1 ;
					if ( tc[j].lcnt > 0 && tc[j].rcnt > 0 )
						time2 = 2 ;
					if ( tc[j].abundance > factor * abundance )
						transcripts[ tcnt ].abundance += ( tc[j].support * ( abundance / tc[j].normAbund ) * factor * time2 ) ;
					else
						transcripts[ tcnt ].abundance += ( tc[j].support * ( tc[j].abundance / tc[j].normAbund ) * time2 ) ;

					tc[j].abundance -= factor * abundance ;	
					
					if ( attr.solveRemaining || USE_SET_COVER )
						tc[j].abundance = 0 ;
				}
				if ( tc[j].abundance < 0 )
					tc[j].abundance = 0 ;
				//tc[j].abundance = 0 ;
			}
			//tcUsed[j] = true ;
		}
		
		/*for ( i = 0 ; i < tcCnt ; ++i )
			if ( tcUsedCnt[i] == 0 && tc[i].abundance > 0 )
				break ;
		if ( i >= tcCnt )
		{
			break ;
		}*/
		
		++tcnt ;
		if ( tcnt >= MAX_TRANSCRIPT )
		{
			qsort( transcripts, tcnt, sizeof( struct _transcript ), CompTranscripts ) ;
		
			int cnt ;
			cnt = 1 ;
			for ( i = 1 ; i < tcnt ; ++i )
			{
				if ( TranscriptsDifference( transcripts[i], transcripts[i - 1] ) )
				{
					transcripts[cnt] = transcripts[i] ;
					++cnt ;
				}
				else
				{
					transcripts[cnt - 1].abundance += transcripts[i].abundance ;
					free( transcripts[i].eid ) ;
				}
			}
			tcnt = cnt ;		
		}
		if ( tcnt >= MAX_TRANSCRIPT )
		{
			//CleanDpMemo( attr ) ;
			break ;
		}		
		
		//CleanDpMemo( attr ) ;	
		//printf( "Find one: %d %d %d\n", best.ecnt, tcnt, tcCnt ) ;	
		//free( transcript.eid ) ;
		++attr.timeStamp ;
		if ( attr.timeStamp == 0 )
			++attr.timeStamp ;
	} // end of while best
	CleanDpMemo( attr ) ;
	free( attr.f1 ) ;
	for ( i = 0 ; i < attr.size ; ++i )
		free( attr.f2[i] ) ;
	free( attr.f2 ) ;
	if ( attr.size <= USE_F3 )
	{
		for ( i = 0 ; i < attr.size ; ++i )
			for ( j = 0 ; j < attr.size ; ++j )
				free( attr.f3[i][j] ) ;
		for ( i = 0 ; i < attr.size ; ++i )
			free( attr.f3[i] ) ;
		free( attr.f3 ) ;
	}
	free( attr.stride ) ;
	free( attr.hash ) ;

	for ( i = 0 ; i < attr.size * attr.size ; ++i )
		if ( attr.pairExonTc[i].cnt > 0 )
			free( attr.pairExonTc[i].id ) ;
	free( attr.pairExonTc ) ;
	free( attr.pairExonTcIds ) ;
	free( attr.bufferEid ) ;

	//for ( i = 0 ; i < tcnt ; ++i )
	//	alltranscripts[i] = transcripts[i] ;
	//atcnt = tcnt ;

	qsort( transcripts, tcnt, sizeof( struct _transcript ), CompTranscripts ) ;
	if ( tcnt == 0 )
		return ;
	int cnt ;
	cnt = 1 ;
	for ( i = 1 ; i < tcnt ; ++i )
	{
		if ( TranscriptsDifference( transcripts[i], transcripts[i - 1] ) )
		{
			transcripts[cnt] = transcripts[i] ;
			++cnt ;
		}
		else
			transcripts[cnt - 1].abundance += transcripts[i].abundance ;
	}
	tcnt = cnt ;
	/*if ( transcriptId > 26 )
		exit( 1 ) ;*/
}

int FindFarthestBefore( int tag, struct _exonNode nodes[], struct _exon exons[] ) 
{
	if ( nodes[tag].before != -1 )
		return nodes[tag].before ;
	if ( !nodes[tag].pcnt )
		return nodes[tag].before = exons[tag].start ;
	int tmp, min = INF ;
	int i ;
	for ( i = 0 ; i < nodes[tag].pcnt ; ++i )
	{
		tmp = FindFarthestBefore( nodes[tag].prev[i], nodes, exons ) ;
		if ( tmp < min )
			min = tmp ;
	}
	return nodes[tag].before = min ;
}

int FindFarthest( int tag, struct _exonNode nodes[], struct _exon exons[] )
{
	if ( nodes[tag].farthest != -1 )
		return nodes[tag].farthest ;
	if ( !nodes[tag].ncnt )
		return nodes[tag].farthest = exons[tag].end ;
	int tmp, max = -1 ;
	int i ;
	for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
	{
		tmp = FindFarthest( nodes[tag].next[i], nodes, exons ) ;
		if ( tmp > max )
			max = tmp ;
	}
	return nodes[tag].farthest = max ;
}

int SubTranscriptCount( int tag, struct _exonNode nodes[] )
{
	if ( nodes[tag].subcnt != -1 )
		return nodes[tag].subcnt ;
	if ( !nodes[tag].ncnt )	
		return nodes[tag].subcnt = 1 ;
	int tmp = 1, i ;
	for ( i = 0 ; i < nodes[tag].ncnt ; ++i )
		tmp += SubTranscriptCount( nodes[tag].next[i], nodes ) ;
	return nodes[tag].subcnt = tmp ;
}

// Test whether two exons are connected, if so return the shortest "exon-bp" distance between them.
// Otherwise, return INF ; 
int AreExonsConnected( int f, int t, int *reach, int offset, int size, struct _exon exons[], struct _exonNode nodes[] ) 
{
	int tag = ( f - offset ) * size + t - offset ;
	int i ;
	int tmp ;
	int len = exons[f].end - exons[f].start + 1 ;

	if ( reach[tag] != -1 )
		return reach[tag] ;
	if ( f == t )
		return reach[tag] = len ;
	if ( f > t )
		return reach[tag] = INF ;
	for ( i = 0 ; i < nodes[f].ncnt ; ++i )
	{
		if ( tmp = AreExonsConnected( nodes[f].next[i], t, reach, offset, size, exons, nodes ) )
		{
			if ( tmp != INF && ( reach[tag] == -1 || tmp + len < reach[tag] ) )
				reach[tag] = tmp + len ;
		}
	}
	if ( reach[tag] == -1 )
		reach[tag] = INF ;
	return reach[tag] ;
}

// Count the number of possible path between two exons where the path length is bounded 
INT64 PathsBetweenExons( int f, int t, INT64 *ways, int *reach, int offset, int size, struct _exon exons[], struct _exonNode nodes[] ) 
{
	int tag = ( f - offset ) * size + t - offset ;
	int i ;
	int lenf = exons[f].end - exons[f].start + 1, lent = exons[t].end - exons[t].start + 1 ;
	INT64 tmp ;

	if ( ways[tag] != -1 )
		return ways[tag] ;
	if ( f == t )
		return ways[tag] = 1 ;
	if ( f > t )
		return ways[tag] = 0 ;
	ways[tag] = 0 ;
	for ( i = 0 ; i < nodes[f].ncnt ; ++i )
	{
		tmp = PathsBetweenExons( nodes[f].next[i], t, ways, reach, offset, size, exons, nodes ) ;
		if ( AreExonsConnected( nodes[f].next[i], t, reach, offset, size, exons, nodes ) - lenf - lent > 
			FRAG_LENGTH + 2 * FRAG_STD - 2 * READS_LENGTH )
			tmp = 0 ;
		ways[tag] += tmp ;
	}
	return ways[tag] ;
}

// Test whether each constraints can be satisfies, if not, set its support to -1.
// TODO: Handle the case that the splice in the reads can not create a splice junction and the its two ends still exist in exons.
int VerifyConstraints( struct _exonNode nodes[], struct _exon exons[], int from, int to, int tcCnt )
{
	int i, j, k, l ;
	int a, b ;
	int *reach ;
	int newTcCnt = tcCnt ;
	struct _transcriptConstraint *tc = transcriptConstraints ;
	reach = ( int * )malloc( sizeof( int ) * ( to - from + 1 ) * ( to - from + 1 ) ) ;
	memset( reach, -1, sizeof( int ) * ( to - from + 1 ) * ( to - from + 1 ) ) ;

	for ( j = 0 ; j < tcCnt ; ++j )
	{
		if ( tc[j].type == 0 )
		{
			//if ( !nodes[ tc[j].info[0] ].pcnt && !nodes[tc[j].info[0] ].ncnt )
			//{
			//	tc[j].support = -1 ;
			//}
			continue ;
		}
		else if ( tc[j].type == 2 )
			continue ;
		// Each sub constraints should compatible to some exons
		k = from ;
		for ( i = 0 ; i < tc[j].lcnt ; ++i )
		{
			for ( ; k <= to ; ++k )
			{
				//if ( !nodes[k].pcnt && !nodes[k].ncnt )
				//	continue ;

				if ( exons[k].start <= tc[j].leftInfo[i][0] && exons[k].end >= tc[j].leftInfo[i][1] && 
					( i == tc[j].lcnt - 1 || exons[k].end == tc[j].leftInfo[i][1] ) &&
					( i == 0 || exons[k].start == tc[j].leftInfo[i][0] ) )
				{
					if ( i < tc[j].lcnt - 1 )
					{
						for ( l = 0 ; l < nodes[k].ncnt ; ++l )
						{
							if ( exons[ nodes[k].next[l] ].start == tc[j].leftInfo[i + 1][0] )
								break ;
						}
						if ( l >= nodes[k].ncnt )
							continue ;
					}
					break ;
				}
			}

			if ( k > to )
				break ;
		}
		if ( i < tc[j].lcnt )
		{
			//tcUsedCnt[j] = 1 ;
			//tc[j].abundance = tc[j].support = 0 ;
			tc[j].lcnt = 0 ;
			//tc[j].support = -1 ;
			//continue ;
			//continue ;
			if ( tc[j].rcnt == 0 )
			{
				tc[j].support = -1 ;
				continue ;
			}
		}

		k = from ;
		for ( i = 0 ; i < tc[j].rcnt ; ++i )
		{
			for ( ; k <= to ; ++k )
			{
				if ( !nodes[k].pcnt && !nodes[k].ncnt )
					continue ;
				if ( exons[k].start <= tc[j].rightInfo[i][0] && exons[k].end >= tc[j].rightInfo[i][1] && 
						( i == tc[j].rcnt - 1 || exons[k].end == tc[j].rightInfo[i][1] ) &&
						( i == 0 || exons[k].start == tc[j].rightInfo[i][0] ) )
				{
					if ( i < tc[j].rcnt - 1 )
					{
						for ( l = 0 ; l < nodes[k].ncnt ; ++l )
						{
							if ( exons[ nodes[k].next[l] ].start == tc[j].rightInfo[i + 1][0] )
								break ;
						}
						if ( l >= nodes[k].ncnt )
							continue ;
					}
					break ;
				}
			}

			if ( k > to )
				break ;
		}
		if ( i < tc[j].rcnt )
		{
			//tcUsedCnt[j] = 1 ;
			//tc[j].abundance = tc[j].support = 0 ;
			//continue ;
			tc[j].rcnt = 0 ;

			//tc[j].support = -1 ;
			//continue ;
			if ( tc[j].lcnt == 0 )
			{
				tc[j].support = -1 ;
				continue ;
			}
		}
		else if ( tc[j].lcnt == 0 )
		{
			// If the left part is empty.
			tc[j].lcnt = tc[j].rcnt ;
			for ( k = 0 ; k < tc[j].rcnt ; ++k )
			{
				tc[j].leftInfo[k][0] = tc[j].rightInfo[k][0] ;
				tc[j].leftInfo[k][1] = tc[j].rightInfo[k][1] ;
			}
			tc[j].rcnt = 0 ;
		}

		// Now test the connectivity
		if ( tc[j].lcnt > 0 && tc[j].rcnt > 0 && tc[j].leftInfo[ tc[j].lcnt - 1 ][1] < tc[j].rightInfo[0][0] )
		{
			int dist ;
			int flank ;
			for ( a = from ; a <= to ; ++a )
			{
				if ( !nodes[a].pcnt && !nodes[a].ncnt )
					continue ;
				if ( exons[a].start <= tc[j].leftInfo[ tc[j].lcnt - 1 ][0] && exons[a].end >= tc[j].leftInfo[ tc[j].lcnt - 1 ][1] &&
						( tc[j].lcnt == 1 || exons[a].start == tc[j].leftInfo[ tc[j].lcnt - 1 ][0] ) )	
				{
					for ( b = from ; b <= to ; ++b )
					{
						if ( !nodes[b].pcnt && !nodes[b].ncnt )
							continue ;
						if ( exons[b].start <= tc[j].rightInfo[0][0] && exons[b].end >= tc[j].rightInfo[0][1] &&
								( tc[j].rcnt == 1 || exons[b].end == tc[j].rightInfo[0][1] ) )
						{
							if ( a == b )
								break ;
							//i = a < b ? a : b ;
							//k = a > b ? a : b ;
							i = a ;
							k = b ;
							flank = 0 ;
							dist = AreExonsConnected( i, k, reach, from, to - from + 1, exons, nodes ) ;
							
							if ( i < k )
							{
								flank = exons[a].end - tc[j].leftInfo[ tc[j].lcnt - 1][1] ;
								flank += tc[j].rightInfo[0][0] - exons[b].start ;
							}

							if ( dist - ( exons[a].end - exons[a].start + 1 ) - ( exons[b].end - exons[b].start + 1 ) + flank 
								< FRAG_LENGTH + 2 * FRAG_STD - 2 * READS_LENGTH  )
							{
								break ;
							}
						}
					}
					if ( b <= to )
						break ;
				}
			}
			if ( a > to )
			{
				//tcUsedCnt[j] = 1 ;
				//tc[j].abundance = tc[j].support = 0 ;
				//continue ;
				tc[j].support = -1 ;
				//continue ;
			}
		} // if for connectivity
	}

	// Remove the unsatisfiable ones
	k = 0 ; 
	for ( i = 0 ; i < newTcCnt ; ++i )
	{
		if ( tc[i].support == -1 )
		{
			continue ;
		}
		tc[k] = tc[i] ;
		++k ;
	}
	newTcCnt = k ;
	if ( k == 0 )
	{
		free( reach ) ;
		return 0 ;
	}
	qsort( tc, newTcCnt, sizeof( *tc ), CompTranscriptConstraint ) ;
	/*printf( "%d\n", points[0] ) ;
	for ( i = 0 ; i < tcCnt ; ++i )
	{
		printf( "%d %d %d ", alltc[i].type, alltc[i].lcnt, alltc[i].rcnt ) ;
		for ( j = 0 ; j < alltc[i].lcnt ; ++j )
			printf( "(%d %d)", alltc[i].leftInfo[j][0], alltc[i].leftInfo[j][1] ) ;
		for ( j = 0 ; j < alltc[i].rcnt ; ++j )
			printf( "(%d %d)", alltc[i].rightInfo[j][0], alltc[i].rightInfo[j][1] ) ;
		printf( "\n" ) ;
	}*/

	k = 1 ;
	for ( i = 1 ; i < newTcCnt ; ++i )
	{
		if ( TranscriptConstraintDifference( &tc[i], &tc[i - 1] ) )
		{
		//	printf( "	support = %d\n", tc[k - 1].support ) ;
			tc[k] = tc[i] ;
			++k ;
		}
	}
	tcCnt = k ;

	// Split a paired end reads, if the combination of the insertion is too much.
	//printf( "Split Reads: %d\n", tcCnt ) ;
	/*newTcCnt = tcCnt ;
	INT64 *ways ;
	ways = (INT64 *)malloc( sizeof( *ways ) * ( to - from + 1 ) * ( to - from + 1 ) ) ;
	memset( ways, -1, sizeof( *ways ) * (to - from + 1 ) * ( to - from + 1 ) ) ;
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		INT64 comb = 0 ;
		if ( tc[j].type != 1 || tc[j].lcnt == 0 || tc[j].rcnt == 0 || tc[j].leftInfo[ tc[j].lcnt - 1 ][1] >= tc[j].rightInfo[0][0] )
			continue ;

		for ( a = from ; a <= to ; ++a )
		{
			if ( !nodes[a].pcnt && !nodes[a].ncnt )
				continue ;
			if ( exons[a].start <= tc[j].leftInfo[ tc[j].lcnt - 1 ][0] && exons[a].end >= tc[j].leftInfo[ tc[j].lcnt - 1 ][1] &&
					( tc[j].lcnt == 1 || exons[a].start == tc[j].leftInfo[ tc[j].lcnt - 1 ][0] ) )	
			{
				for ( b = from ; b <= to ; ++b )
				{
					if ( !nodes[b].pcnt && !nodes[b].ncnt )
						continue ;
					if ( exons[b].start <= tc[j].rightInfo[0][0] && exons[b].end >= tc[j].rightInfo[0][1] &&
							( tc[j].rcnt == 1 || exons[b].end == tc[j].rightInfo[0][1] ) )
					{
						//if ( a == b )
						//	break ;
						//i = a < b ? a : b ;
						//k = a > b ? a : b ;
						i = a ;
						k = b ;
						//printf( "%d %d\n", i, k ) ;
						comb += PathsBetweenExons( i, k, ways, reach, from, to - from + 1, 
							exons, nodes ) ;

					} // if b is good
				} // for b
			} // if exon[a] is good
		} // for a
		printf( "comb: %lld\n", comb ) ;
		//if ( comb == 0 )
		//	exit( 1 ) ;
		if ( comb >= 10 )
		{
			// The insertion is too complicated, and we should split this read
			// NOTE: When this happens, the gene must be complicated, and we definitely will use
			// 	DP. So the enumeration method will not be affected. 
			tc[j].effectiveCount = 0.5 ;
			tc[newTcCnt] = tc[j] ;

			tc[newTcCnt].lcnt = tc[j].rcnt ;
			tc[newTcCnt].rcnt = 0 ;
			for ( i = 0 ; i < tc[j].rcnt ; ++i )
			{
				tc[newTcCnt].leftInfo[i][0] = tc[j].rightInfo[i][0] ;
				tc[newTcCnt].leftInfo[i][1] = tc[j].rightInfo[i][1] ;
			}

			++newTcCnt ;
			tc[j].rcnt = 0 ;
		}
	}
	free( ways ) ;
	
	qsort( tc, newTcCnt, sizeof( *tc ), CompTranscriptConstraint ) ;
	k = 1 ;
	for ( i = 1 ; i < newTcCnt ; ++i )
	{
		if ( TranscriptConstraintDifference( &tc[i], &tc[i - 1] ) )
		{
		//	printf( "	support = %d\n", tc[k - 1].support ) ;
			tc[k] = tc[i] ;
			++k ;
		}
	}
	tcCnt = k ;*/
	//for ( i = 0 ; i < tcCnt ; ++i )
	//	printf( "%d: (%d %d)\n", tc[i].type, tc[i].lcnt, tc[i].rcnt) ;
	//printf( "Finished Split Reads: %d\n", tcCnt ) ;

	// Calculate the end of type1 constraints 
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		if ( tc[j].type == 0 )
		{
			tc[j].info[1] = exons[ tc[j].info[0] ].end ;
		}
		else if ( tc[j].type == 1 )
		{	
			int a = -1, b = -1 ;
			if ( tc[j].lcnt )
				a = tc[j].leftInfo[ tc[j].lcnt - 1][1] ;			
			if ( tc[j].rcnt )
				b = tc[j].rightInfo[ tc[j].rcnt - 1 ][1] ;
			tc[j].info[1] = a > b ? a : b ;
		}
	}

	free( reach ) ;
	return tcCnt ;
}

// When there are too many constraints, we will clean up them.
// We might lose some information.
int CleanTranscriptConstraints( int abdCnt, int tcCnt )
{
	int i, j, k ;
	int ret = 0 ;
	// Keep the constraints not type 1.

	return ret ;	
}

void RefineTranscripts( char *chrom, int from, int to, struct _exonNode nodes[], struct _exon exons[], 
	struct _transcript *transcripts, int tcnt,
	struct _transcriptConstraint *tc, const int &tcCnt, 
	struct _pair *geneMerge, int geneMergeCnt, int evidenceIndFrom, int evidenceIndTo )
{
	int i, j, k ;
	int mainStrand ;
	//struct _transcriptConstraint *tc = transcriptConstraints ;
	/*for ( i = 0 ; i < tcnt ; ++i )
	{
		btable[i].Init( tcCnt ) ;
		//printf( "%d: %d %d %d\n", tcCnt, transcripts[i].eid[0], transcripts[i].eid[1], transcripts[i].eid[2] ) ;
		for ( j = 0 ; j < tcCnt ; ++j )
			if ( IsConstraintInTranscript( transcripts[i].eid, transcripts[i].ecnt, tc[j], exons ) == 1 )
			{
				btable[i].Set( j ) ;
			}
	}*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		for ( j = 0 ; j < transcripts[i].ecnt ; ++j )
			nodes[ transcripts[i].eid[j] ].geneId = transcripts[i].geneId ;
	}

	// Remove shadow transcripts
	int plusCnt = 0, minusCnt = 0 ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].strand == 1 )
			++plusCnt ;
		if ( transcripts[i].strand == 0 )
			++minusCnt ;
	}
	if ( plusCnt > minusCnt )
		mainStrand = 1 ;
	else
		mainStrand = 0 ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].ecnt != 2 || transcripts[i].strand == mainStrand )
			continue ;
		for ( j = 0 ; j < tcnt ; ++j )
		{
			if ( transcripts[j].ecnt <= 2 || transcripts[j].strand == transcripts[i].strand )
				continue ;
			for ( k = 0 ; k < transcripts[j].ecnt - 1 ; ++k )
			{
				if ( exons[ transcripts[j].eid[k] ].end >= exons[ transcripts[i].eid[0] ].end - 20 
					&& exons[ transcripts[j].eid[k] ].end <= exons[ transcripts[i].eid[0] ].end + 20 
					|| ( exons[transcripts[j].eid[k + 1]].start >= exons[ transcripts[i].eid[1] ].start - 20 
						&& exons[transcripts[j].eid[k + 1]].start <= exons[ transcripts[i].eid[1] ].start + 20 ) )
				{
					transcripts[i].geneId = -1 ;
					j = tcnt ;
					break ;
				}
			}
		}
		
		/*if ( transcripts[i].ecnt < 2 || transcripts[i].strand == mainStrand )
			continue ;
		for ( j = 0 ; j < tcnt ; ++j )
		{
			if ( transcripts[j].strand == transcripts[i].strand )
				continue ;
			for ( k = 0 ; k < transcripts[j].ecnt - 1 ; ++k )
			{
				if ( exons[ transcripts[j].eid[k] ].end >= exons[ transcripts[i].eid[0] ].end - 20 
					&& exons[ transcripts[j].eid[k] ].end <= exons[ transcripts[i].eid[0] ].end + 20 
					&& exons[transcripts[j].eid[k + 1]].start >= exons[ transcripts[i].eid[1] ].start - 20 
					&& exons[transcripts[j].eid[k + 1]].start <= exons[ transcripts[i].eid[1] ].start + 20 )
				{
					transcripts[i].geneId = -1 ;
					j = tcnt ;
					break ;
				}
			}
		}*/
	}
		
	if ( 1 ) //tcnt > 1 )
	{
		for ( i = 0 ; i < tcnt ; ++i )
		{
			if ( transcripts[i].ecnt != 2 )
				continue ;
			if ( exons[ transcripts[i].eid[0] ].end - exons[ transcripts[i].eid[0] ].start + 1 < 110
					&& exons[ transcripts[i].eid[1] ].end - exons[ transcripts[i].eid[1] ].start + 1 < 110 )
			{
				int cnt = 0 ;
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( tc[j].type != 1 )
						continue ;
					if ( IsConstraintInTranscript( transcripts[i].eid, transcripts[i].ecnt, tc[j], exons ) )
					{
						cnt += tc[j].support ;
					}
				}

				if ( cnt < 10 )
					transcripts[i].geneId = -1 ;
			}
		}
	}

	/*if ( 1 )
	{
		for ( i = 0 ; i < tcnt ; ++i )
		{
			int len = 0 ;
			for ( j = 0 ; j < transcripts[i].ecnt ; ++j )
				len += exons[ transcripts[i].eid[j] ].end - exons[ transcripts[i].eid[j] ].start + 1 ;
			if ( len < 200 )
				transcripts[i].geneId = -1 ;
		}
	}*/
	// Compute the FPKM 
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 )
			continue ;
		transcripts[i].FPKM = GetTranscriptAbundance( exons, &transcripts[i] ) ;
		geneAbundance[ transcripts[i].geneId ] = -1 ;
	}
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 )
			continue ;
		if ( transcripts[i].FPKM > geneAbundance[ transcripts[i].geneId ] )
			geneAbundance[ transcripts[i].geneId ] = transcripts[i].FPKM ;
	}

	// Possible wrong TSS or TES exon
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 )
			continue ;
		for ( j = 0 ; j < transcripts[i].ecnt ; ++j )
		{
			++nodes[ transcripts[i].eid[j] ].supportTranscript ;
		}
	}

	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 || transcripts[i].FPKM > 0.3 * geneAbundance[ transcripts[i].geneId ] ) 
			continue ;

		if ( nodes[ transcripts[i].eid[0] ].supportTranscript == 1 && 
			exons[ transcripts[i].eid[0] ].end - exons[ transcripts[i].eid[0]].start + 1 < 70 )
		{
			for ( j = 0 ; j < tcnt ; ++j )
			{
				if ( transcripts[j].geneId != transcripts[i].geneId )
					continue ;
				if ( nodes[ transcripts[j].eid[0] ].supportTranscript >= 2 )
				{
					transcripts[i].geneId = -1 ;
					break ;
				}
			}
		}

		k = transcripts[i].ecnt - 1 ;
		if ( nodes[ transcripts[i].eid[k] ].supportTranscript == 1 && 
			exons[ transcripts[i].eid[k] ].end - exons[ transcripts[i].eid[k] ].start + 1 < 70 )
		{
			for ( j = 0 ; j < tcnt ; ++j )
			{
				if ( transcripts[j].geneId != transcripts[i].geneId )
					continue ;
				if ( nodes[ transcripts[j].eid[ transcripts[j].ecnt - 1] ].supportTranscript >= 2 )
				{
					transcripts[i].geneId = -1 ;
					break ;
				}
			}
		}
	}


	// Remove the transcripts containing constraints, like gene merges.
	int *tid = (int *)malloc( sizeof( int ) * tcnt ) ;
	int tidCnt ;
	for ( i = 0 ; i < geneMergeCnt ; ++i )
	{
		tidCnt = 0 ;
		double maxFPKM = -1 ;
		int maxFPKMTag = 0 ;
		double totalAbundance = 0 ;
		int totalExonLen = 0 ;
		//printf( "try merge: %d %d\n", exons[ geneMerge[i].a].end, exons[ geneMerge[i].b ].start ) ;

		for ( j = 0 ; j < tcnt ; ++j )
		{
			if ( transcripts[j].geneId == -1 )
				continue ;
			for ( k = 0 ; k < transcripts[j].ecnt - 1 ; ++k )
			{
				if ( transcripts[j].eid[k] == geneMerge[i].a && transcripts[j].eid[k + 1] == geneMerge[i].b )
				{
					tid[ tidCnt ] = j ;
					++tidCnt ;
					if ( transcripts[j].FPKM > maxFPKM ) 
					{
						maxFPKM = transcripts[j].FPKM ;
						maxFPKMTag = j ;
					}
					break ;
				}
			}
		}

		maxFPKM = -1 ;
		for ( j = 0 ; j < tcnt ; ++j )
		{
			if ( transcripts[j].geneId != transcripts[ maxFPKMTag ].geneId )
				continue ;
			totalAbundance += transcripts[j].abundance ;
		}
		for ( j = from ; j <= to ; ++j )
		{
			if ( nodes[j].geneId != transcripts[ maxFPKMTag ].geneId )
				continue ;
			totalExonLen += exons[j].end - exons[j].start + 1 ;
		}

		if ( totalAbundance * READS_LENGTH / (double)totalExonLen < 20 )
		{
			//printf( "No merge: %lf / %d = %lf\n", totalAbundance, totalExonLen, totalAbundance / totalExonLen ) ;
			continue ;
		}

		//printf( "merge: %d %d\n", exons[ geneMerge[i].a].end, exons[ geneMerge[i].b ].start ) ;
		for ( j = 0 ; j < tidCnt ; ++j )
		{
			if ( tid[j] != maxFPKMTag )
			{
				transcripts[ tid[j] ].geneId = -1 ;
			}
		}
	}
	free( tid ) ;

	for ( i = 0 ; i < tcnt ; ++i )
	{
		/*if ( exons[transcripts[i].eid[0]].start == 108769308 )
		{
			//printf( "hi %d %lf\n", transcripts[i].geneId, geneAbundance[ transcripts[i].geneId ] ) ;
			OutputTranscript( chrom, exons, &transcripts[i] ) ;	
		}*/
		if ( transcripts[i].geneId != -1 && transcripts[i].FPKM < geneAbundance[ transcripts[i].geneId ] * FPKM_FRACTION )
		{
			// Check whether it is in the evidence. If not, then we remove it.
			//if ( exons[ transcripts[i].eid[0] ].end == 45340498 )
			//	printf( "hi %d %d: %d %d\n", evidenceIndFrom, evidenceIndTo, evidences[evidenceIndFrom].exons[0].end, exons[ transcripts[i].eid[0] ].end ) ;
			if ( evidences != NULL )
			{
				for ( j = evidenceIndFrom ; j <= evidenceIndTo ; ++j )
				{
					if ( evidences[j].ecnt != transcripts[i].ecnt || transcripts[i].ecnt == 1 )
						continue ;
					for ( k = 0 ; k < transcripts[i].ecnt - 1 ; ++k )
					{
						if ( evidences[j].exons[k].end != exons[ transcripts[i].eid[ k ] ].end ||
								evidences[j].exons[k + 1].start != exons[ transcripts[i].eid[ k + 1 ] ].start )
						{	
							break ;
						}
					}
					if ( k >= transcripts[i].ecnt - 1 )
						break ;

				}
				if ( j <= evidenceIndTo )
				{
					continue ;
				}
				else
					transcripts[i].geneId = -1 ;
			}
			else
				transcripts[i].geneId = -1 ;
		}
	}
	
	/*==================================================================
	Remove transcripts that seems duplicated
	====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 )
			continue ;
		int cnt = 0, uniqCnt = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( tc[j].type != 1 )
				continue ;
			if ( IsConstraintInTranscript( transcripts[i].eid, transcripts[i].ecnt, tc[j], exons ) )
			{
				++cnt ;
				if ( tc[j].uniqSupport >= 0.05 * tc[j].support )
					++uniqCnt ;
			}
		}
		//printf( "%d %d\n", uniqCnt, cnt ) ;
		if ( uniqCnt <= 2 && uniqCnt <= 0.1 * cnt )
		{
			transcripts[i].geneId = -1 ;
		}
	}

	/*for ( j = 0 ; j < tcCnt ; ++j )
	{
		if ( tc[j].type != 1 )
			continue ;
		for ( i = 0 ; i < tc[j].lcnt ; ++i )
			printf( "(%d %d) ", tc[j].leftInfo[i][0], tc[j].leftInfo[i][1] ) ;
		for ( i = 0 ; i < tc[j].rcnt ; ++i )
			printf( "[%d %d] ", tc[j].rightInfo[i][0], tc[j].rightInfo[i][1] ) ;
		printf( ": %d %d strand: %d\n", tc[j].uniqSupport, tc[j].support, tc[j].strand ) ;
	}*/

	//int cnt = 0 ;

	pthread_mutex_lock( &outputTxptMutex ) ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].geneId == -1 )
			continue ;
		//transcripts[i].transcriptId = transcriptId ;
		//++transcriptId ;
		//++cnt ;
		//OutputTranscript( chrom, exons, &transcripts[i] ) ;	
		outputTranscripts[ outputTxptCnt ] = transcripts[i] ;
		outputTranscripts[ outputTxptCnt ].eid = ( int * )malloc( sizeof( int ) * transcripts[i].ecnt ) ;
		memcpy( outputTranscripts[ outputTxptCnt ].eid, transcripts[i].eid, sizeof( int ) * transcripts[i].ecnt ) ;
		++outputTxptCnt ;
	}
	pthread_mutex_unlock( &outputTxptMutex ) ;

	//if ( cnt > 100 )
	//printf( "RefineTranscripts: %d\n", cnt ) ;
	//fflush( stdout ) ;

	return ;
}

// The wrapper for pthread calling
void *PickTranscripts_Thread( void *arg )
{
	int i, j, k ;
	int *visit ;
	struct _pair *geneMerge ;
	struct _transcript *transcripts ;
	int geneMergeCnt = 0 ;
	bool useDP ;
	int tcnt = 0, atcnt = 0 ;

	char *chrom = ( ( struct _pthreadArgPickTranscripts *)arg )->chrom ;
	struct _exon *exons = ( ( struct _pthreadArgPickTranscripts *)arg )->exons ;
	int exonCnt = ( ( struct _pthreadArgPickTranscripts *)arg )->exonCnt ;
	struct _exonNode *nodes = ( ( struct _pthreadArgPickTranscripts *)arg )->nodes ;
	struct _transcriptConstraint *tc = ( ( struct _pthreadArgPickTranscripts *)arg )->tc ;		
	const int tcCnt = ( ( struct _pthreadArgPickTranscripts *)arg )->tcCnt ;
	int from = ( ( struct _pthreadArgPickTranscripts *)arg )->from ;
	int to = ( ( struct _pthreadArgPickTranscripts *) arg )->to ;

	geneMerge = ( struct _pair * )malloc( sizeof( struct _pair ) * ( to - from + 1 ) * ( to - from + 1 ) ) ;	
	transcripts = ( struct _transcript *)malloc( sizeof( struct _transcript ) * MAX_TRANSCRIPT ) ;
	
	for ( i = from ; i <= to ; ++i )
	{
		//if ( nodes[i].pcnt )
		//	continue ;
		FindFarthestBefore( i, nodes, exons ) ;
		FindFarthest( i, nodes, exons ) ;
	}

	// Find the junction that is a possible gene merge
	for ( i = from ; i <= to ; ++i )
	{
		if ( nodes[i].ncnt >= 2 )
		{
			for ( j = 0 ; j < nodes[i].ncnt ; ++j )
			{
				k = nodes[i].next[j] ;
				// The junction between node i and k.
				if ( nodes[ k ].pcnt >= 2 && nodes[i].ncnt >= 2 )
				{
					bool isGeneMerge = true ;
					int n, p ;
					for ( n = 0 ; n < nodes[i].ncnt ; ++n )
					{
						if ( nodes[i].next[n] == k )
							continue ;
						for ( p = 0 ; p < nodes[k].pcnt ; ++p )
						{
							if ( nodes[k].prev[p] == i )
								continue ;

							if ( nodes[ nodes[i].next[n] ].farthest >= nodes[ nodes[k].prev[p] ].before )
							{
								isGeneMerge = false ;
								break ;
							}
						}
					}

					if ( isGeneMerge == true )
					{
						/*for ( p = 0 ; p < nodes[k].pcnt ; ++p )
							if ( nodes[k].prev[p] == i )
								break ;
						for ( n = j ; n < nodes[i].ncnt ; ++n )
							nodes[i].next[n] = nodes[i].next[n + 1] ;
						--nodes[i].ncnt ;
						
						for ( ; p < nodes[k].pcnt ; ++p )
							nodes[k].prev[p] = nodes[k].prev[p + 1] ;
						--nodes[k].pcnt ;
						--j ;*/
						geneMerge[ geneMergeCnt ].a = i ;
						geneMerge[ geneMergeCnt ].b = k ;
						++geneMergeCnt ;
						//printf( "Gene merge: %d %d\n", exons[i].end, exons[k].start ) ;
					}
				}
			}
		}
	}

	k = 0 ;
	for ( i = from ; i <= to ; ++i )
	{
		if ( nodes[i].pcnt )
			continue ;
		k += SubTranscriptCount( i, nodes ) ;
	}
	atcnt = k ;
	useDP = false ;
	if ( k <= USE_DP )
	{
		for ( i = from ; i <= to ; ++i )
			if ( nodes[i].subcnt > USE_DP )
			{
				useDP = true ;
				break ;
			}
	}
	else
		useDP = true ;

	if ( !useDP )
	{
		// If there are too many constraints
		if ( tcCnt * k > 600000 || tcCnt * k < 0 )
			useDP = true ;
	}

	if ( USE_SET_COVER )
		useDP = false ;
	if ( !useDP )
	{
		// If the possible # of transcript is not much
		visit = ( int * )malloc( sizeof( int ) * ( exonCnt ) ) ;
		struct _transcript *alltranscripts = (struct _transcript *)malloc( sizeof( *alltranscripts ) * (atcnt + 10 ) ) ;
		//printf( "estimated atcnt:%d\n", atcnt ) ;
		atcnt = 0 ;
		for ( i = from ; i <= to ; ++ i )
		{
			if ( nodes[i].pcnt )
				continue ;	
			EnumerateTranscript( i, visit, 0, nodes, exons, tcCnt, alltranscripts, atcnt ) ;
		}
				
		//printf( "actual atcnt: %d\n", atcnt ) ;
	//printf( "before pick %d\n", atcnt ) ;
		PickTranscripts( chrom, nodes, exons, alltranscripts, atcnt, transcripts, tcnt, tc, tcCnt ) ;
		RefineTranscripts( chrom, from, to, nodes, exons, 
			transcripts, tcnt, tc, tcCnt, geneMerge, geneMergeCnt,
			( ( struct _pthreadArgPickTranscripts *)arg )->evidenceIndFrom, 
			( ( struct _pthreadArgPickTranscripts *)arg )->evidenceIndTo ) ;
	//printf( "after pick\n") ;
		free( visit ) ;

		for ( i = 0 ; i < atcnt; ++i )
			free( alltranscripts[i].eid ) ;
		free( alltranscripts ) ;
	}
	else
	{
		PickTranscriptsByDP( chrom, from, to, nodes, exons, transcripts, tcnt, tc, tcCnt ) ;
		RefineTranscripts( chrom, from, to, nodes, exons, 
			transcripts, tcnt, tc, tcCnt, geneMerge, geneMergeCnt,
			( ( struct _pthreadArgPickTranscripts *)arg )->evidenceIndFrom, 
			( ( struct _pthreadArgPickTranscripts *)arg )->evidenceIndTo ) ;

		for ( i = 0 ; i < tcnt ; ++i )
			free( transcripts[i].eid ) ;
	}
	free( geneMerge ) ;	
	free( transcripts ) ;
	free( tc ) ;
	free( arg ) ;
	//pthread_exit
	pthread_mutex_lock( &pickTxptMutex ) ;
	--currPickTxptThreadsCnt ;
	if ( currPickTxptThreadsCnt == NUM_OF_THREADS - 1 )
	{
		pthread_cond_signal( &idleCond ) ;
	}
	if ( currPickTxptThreadsCnt == 0 )
	{
		pthread_cond_signal( &clearCond ) ;
	}
	pthread_mutex_unlock( &pickTxptMutex ) ;

	pthread_exit( NULL ) ;
}

// Solve the transcripts with exons "from"..."to"
void SolveGene( char *chrom, struct _exonNode nodes[], struct _exon exons[], int exonCnt, int from, int to )
{
	int i, j, k ;
	//int *f, *next ;
	int tcCnt, abdCnt  ;
	struct _pthreadArgPickTranscripts *pArg = ( struct _pthreadArgPickTranscripts *)malloc( sizeof( *pArg ) ) ;

	//f = ( int * )malloc( sizeof( int ) * exonCnt ) ;
	//next = ( int * )malloc( sizeof( int ) * exonCnt ) ;
	//char buffer[1024] ;	
	//fgets( buffer, sizeof( buffer ), td_fpReads.fp ) ;
	//printf( "%s", buffer ) ;
	abdCnt = BuildTranscriptConstraints( chrom, nodes, exons, exons[from].start, transcriptEnd, from, to, true  ) ;
	tcCnt = BuildTranscriptConstraints( chrom, nodes, exons, exons[from].start, transcriptEnd, from, to, false  ) ;
	tcCnt = VerifyConstraints( nodes, exons, from, to, tcCnt ) ;
	/*for ( i = 0 ; i < tcCnt ; ++i )
		printf( "%d ", transcriptConstraints[i].type ) ;
	printf( ":%d\n", tcCnt ) ;*/
	DetermineTranscriptConstraintsAbundance( abdCnt, tcCnt ) ;

	int evidenceIndFrom = -1, evidenceIndTo = -1 ;
	if ( evidences != NULL )
	{
		evidenceIndFrom = eviTag ;
		for ( evidenceIndTo = evidenceIndFrom ; evidenceIndTo < eviCnt ; ++evidenceIndTo )
		{
			if ( evidences[ evidenceIndTo ].exons[0].start > transcriptEnd )
			{
				break ;
			}
		}
		//printf( "hi3: %d %d (%d %d %d) %d %d\n", evidences[evidenceIndFrom].exons[0].start, evidences[ evidenceIndTo ].exons[0].start,
		//	eviCnt, evidenceIndFrom, evidenceIndTo, exons[from].start, transcriptEnd ) ;
	}
	//printf( "%d\n",  ) ;
	//memset( tcUsedCnt, 0, sizeof( int ) * tcCnt ) ;

	/*printf( "%d\n", trCnt ) ;	
	for ( i = 0 ; i < tcCnt ; ++i )
		printf( "(%d %d) ", transcriptConstraints[i].type, transcriptConstraints[i].support ) ;
	printf( "\n" ) ;*/
	/*k = 0 ; 
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		if ( transcriptConstraints[j].support <= 1 )
			continue ;
		transcriptConstraints[k] = transcriptConstraints[j] ;
		++k ;
		if ( transcriptConstraints[k].type == 0 )
			tcUsedCnt[k] = MAX_READ ;
	}
	tcCnt = k ;*/
		// Put the single exons
	/*for ( i = from ; i <= to ; ++i )
	{
		if ( !nodes[i].pcnt && !nodes[i].ncnt )
		{
			transcripts[tcnt].eid = (int *)malloc( sizeof( int ) ) ;
			transcripts[tcnt].eid[0] = i ;
			transcripts[tcnt ].ecnt = 1 ;
			transcripts[tcnt].geneId = geneId[i] ;
			transcripts[tcnt].strand = exons[i].strand ;
			++tcnt ;
		}
	}*/

	//printf( "Transcript Counts: %d %d %d %d\n", to - from + 1, k, tcCnt, trCnt ) ;

	if ( VERBOSE )
	{
		printf( "# Assembling from %d to %d, which have %d exons, %d constraints.\n", 
			exons[from].start, transcriptEnd, to - from + 1, tcCnt ) ;
		fflush( stdout ) ;
	}

	pArg->chrom = chrom ;
	pArg->exons = exons ;
	pArg->exonCnt = exonCnt ;
	pArg->nodes = nodes ;
	pArg->tc = ( struct _transcriptConstraint * )malloc( sizeof( struct _transcriptConstraint ) * tcCnt ) ;
	memcpy( pArg->tc, transcriptConstraints, sizeof( struct _transcriptConstraint ) * tcCnt ) ;
	pArg->tcCnt = tcCnt ;
	pArg->from = from ;
	pArg->to = to ;
	pArg->evidenceIndFrom = evidenceIndFrom ;
	pArg->evidenceIndTo = evidenceIndTo ;

	// Test available thread here
	pthread_mutex_lock( &pickTxptMutex ) ;
	if ( currPickTxptThreadsCnt >= NUM_OF_THREADS )
		pthread_cond_wait( &idleCond, &pickTxptMutex ) ;
	++currPickTxptThreadsCnt ;
	pthread_t thread ;
	pthread_create( &thread, &pthreadAttr, PickTranscripts_Thread, (void *)pArg ) ;
	pthread_mutex_unlock( &pickTxptMutex ) ;

	//PickTranscripts_Thread( (void *)pArg ) ;

	//if ( to - from + 1 == 102 )
	//	exit(1) ;

	/*k = 0 ;
	  for ( i = 0 ; i < tcCnt ; ++i )
	  {
	  if ( transcriptConstraints[i].type != 1 )
	  continue ;
	  k += transcriptConstraints[i].support ;
	  }
	  printf( "%d\n", k ) ;*/


	//printf( "SolveGene: %d %d\n", from, to ) ;
	/*while ( 1 )
	  {
	  memset( f, -1, sizeof( int ) * exonCnt ) ;
	  for ( i = from ; i <= to ; ++i )
	  {
	  if ( nodes[i].pcnt )
	  continue ;
	  if ( CoverTheEdges( i, nodes, f, next ) )
	  break ;
	  }
	//printf( "### %d %d\n", i, to ) ;	
	if ( i > to )
	break ;
	int p = i ;
	transcript.ecnt = 0 ;
	transcript.strand = exons[i].strand ;
	transcript.geneId = geneId[i] ;
	transcript.transcriptId = transcriptId ;
	++transcriptId ;
	while ( 1 )
	{
	transcript.eid[ transcript.ecnt ] = p ;
	++transcript.ecnt ;
	if ( next[p] != -1 )
	{
	//printf( "%d %d\n", p, next[p] ) ;
	nodes[p].used[ next[p] ] = true ;
	p = nodes[p].next[ next[p] ] ;
	}
	else
	break ;
	}
	//printf( "=== %d %d\n", i, to ) ;	
	OutputTranscript( chrom, exons, &transcript ) ;
	}*/

	//free( f ) ;
	//free( next ) ;
}

void TranscriptDecider_Go( char *chrom, struct _exon exons[], int exonCnt, struct _evidence *inEvidences, int inEviCnt )
{
	int i, j, k ;
	int pos ;
	int visit[500] ;
	//memset( eidUsedFlag, false, sizeof( eidUsedFlag ) ) ;
	memset( geneId, -1, sizeof( geneId ) ) ;
	memset( geneTranscriptId, 0, sizeof( geneTranscriptId ) ) ;
	
	euCnt = 0 ;
	geneIdCnt = 0 ;
	struct _exonNode *nodes ;
	evidences = inEvidences ;
	eviCnt = inEviCnt ;
	//memset( maxAbundance, -1, sizeof( maxAbundance ) ) ;

	//GetWeight() ;
	//if ( evidences )
	//	eviCnt = Evidence_Extract( chrom, evidences ) ;
	//printf( "Evidence #: %d\n", eviCnt ) ;
	//eviCnt = 0 ;
	eviTag = 0 ;

	// First, use the merged position for  the prev[] and next[] in exons from positions.
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		for ( j = 0 ; j < exons[i].pcnt ; ++j )
			exons[i].prev[j] = SplicesInformation_GetPos( exons[i].prev[j], exons[i].strand ) ;
		// we now may have repeated prev and next.
		int l ;
		for ( l = 0, k = 0 ; l < exons[i].pcnt ; ++l )
		{
			for ( j = 0 ; j < k ; ++j )
				if ( exons[i].prev[j] == exons[i].prev[l] )	
					break ;
			if ( j >= k )
			{
				exons[i].prev[k] = exons[i].prev[l] ;
				++k ;
			}
		}
		exons[i].pcnt = k ;
		for ( j = 0 ; j < exons[i].ncnt ; ++j )
			exons[i].next[j] = SplicesInformation_GetPos( exons[i].next[j], exons[i].strand ) ;
		for ( l = 0, k = 0 ; l < exons[i].ncnt ; ++l )
		{
			for ( j = 0 ; j < k ; ++j )
				if ( exons[i].next[j] == exons[i].next[l] )	
					break ;
			if ( j >= k )
			{
				exons[i].next[k] = exons[i].next[l] ;
				++k ;
			}
		}
		exons[i].ncnt = k ;

	}
	
	/*for ( i = 0 ; i < exonCnt ; ++i )
	{
		printf( "%d %d\n\t", exons[i].start, exons[i].end ) ;
		for ( j = 0 ; j < exons[i].pcnt ; ++j )
			printf( "%d ", exons[i].prev[j] ) ;
		printf( "\n\t" ) ;
		for ( j = 0 ; j < exons[i].ncnt ; ++j )
			printf( "%d ", exons[i].next[j] ) ;
		printf( "\n") ;
	}
	fflush( stdout ) ;*/
	// Convert the exons into the nodes of the splice graph
	nodes = ( struct _exonNode *)malloc( sizeof( struct _exonNode ) * exonCnt ) ;
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		nodes[i].next = ( int * )malloc( sizeof( int ) * ( exons[i].ncnt + 3 )  ) ; // TODO: Make it variable
		nodes[i].used = ( bool * )malloc( sizeof( bool ) * ( exons[i].ncnt + 3 )  ) ;
		nodes[i].prev = ( int * )malloc( sizeof( int ) * ( exons[i].pcnt + 3 ) ) ;
		nodes[i].ncnt = 0 ;
		nodes[i].pcnt = 0 ;
		nodes[i].nsize = exons[i].ncnt + 3 ;
		nodes[i].psize = exons[i].pcnt + 3 ;
		nodes[i].farthest = -1 ;
		nodes[i].before = -1 ;
		nodes[i].subcnt = -1 ;
		nodes[i].weight = 0 ;
		nodes[i].supportTranscript = 0 ;
	}
		
	// Build the adjacent list.
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		for ( k = 0 ; k < exons[i].ncnt ; ++k )
		{
		/*if ( exons[i].end == 3047451 )
		{
			printf( "# hi %d %d: %d %d %d\n", exons[i].end, exons[i].next[k], exons[i + 1].start, exons[i + 2].start, exons[i+3].start  ) ; 
		}*/
			for ( j = i + 1 ; j < exonCnt && exons[j].start <= exons[i].next[k] ; ++j )
			{
		/*if ( exons[i].end == 3047451 )
		{
			printf( "# hi %d %d: %d %d %d\n", exons[i].end, exons[i].next[k], exons[j].start, exons[i].strand, exons[j].strand  ) ; 
		}*/
				if ( exons[j].strand == exons[i].strand && exons[i].next[k] == exons[j].start )
				{
					while ( nodes[i].ncnt >= nodes[i].nsize )
					{
						char *buffer = (char *)malloc( sizeof( int ) * nodes[i].nsize + 4 ) ;
						memcpy( buffer, nodes[i].next, sizeof( int ) * nodes[i].nsize ) ;
						free( nodes[i].next ) ;
						nodes[i].nsize += ( exons[i].ncnt + 2 ) ;
						nodes[i].next = ( int * )malloc( sizeof( int ) * nodes[i].nsize ) ;
						memcpy( nodes[i].next, buffer, sizeof( int ) * nodes[i].ncnt ) ;
						free( buffer ) ;

						free( nodes[i].used ) ;
						nodes[i].used = ( bool * )malloc( sizeof( bool ) * nodes[i].nsize ) ;
					}
					nodes[i].next[ nodes[i].ncnt ] = j ;
					++nodes[i].ncnt ;

					while ( nodes[j].pcnt >= nodes[j].psize )
					{
						char *buffer = (char *)malloc( sizeof( int ) * nodes[j].psize + 4 ) ;
						memcpy( buffer, nodes[j].prev, sizeof( int ) * nodes[j].psize ) ;
						free( nodes[j].prev ) ;
						nodes[j].psize += ( exons[j].pcnt + 2 ) ;
						nodes[j].prev = ( int * )malloc( sizeof( int ) * nodes[j].psize ) ;
						memcpy( nodes[j].prev, buffer, sizeof( int ) * nodes[j].pcnt ) ;
						free( buffer ) ;
					}
					nodes[j].prev[ nodes[j].pcnt ] = i ;
					++nodes[j].pcnt ;
				}
			}
		}
		//if ( exons[i].end == 9572851 )
		//	printf( "### %d %d\n", nodes[i].ncnt, exons[i].next[0] ) ;
		/*if ( nodes[i].ncnt > exons[i].ncnt * 20 )
		{
			printf( "### WARNING! %d %d\n", nodes[i].ncnt, exons[i].ncnt ) ; exit( 1 ) ;
		}*/
		//if ( nodes[i].ncnt )
		//	memset( nodes[i].used, false, sizeof( bool ) * nodes[i].ncnt ) ;
	}


	/*for ( i = 0 ; i < exonCnt ; ++i )
	{
		FindFarthestBefore( i, nodes, exons ) ;
		FindFarthest( i, nodes, exons ) ;
	}
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		if ( nodes[i].ncnt >= 2 )
		{
			for ( j = 0 ; j < nodes[i].ncnt ; ++j )
			{
				k = nodes[i].next[j] ;
				// The junction between node i and k.
				if ( nodes[ k ].pcnt >= 2 && nodes[i].ncnt >= 2 )
				{
					bool geneMerge = true ;
					int n, p ;
					for ( n = 0 ; n < nodes[i].ncnt ; ++n )
					{
						if ( nodes[i].next[n] == k )
							continue ;
						for ( p = 0 ; p < nodes[k].pcnt ; ++p )
						{
							if ( nodes[k].prev[p] == i )
								continue ;

							if ( nodes[ nodes[i].next[n] ].farthest >= nodes[ nodes[k].prev[p] ].before )
							{
								geneMerge = false ;
								break ;
							}
						}
					}

					if ( geneMerge == true )
					{
						for ( p = 0 ; p < nodes[k].pcnt ; ++p )
							if ( nodes[k].prev[p] == i )
								break ;
						for ( n = j ; n < nodes[i].ncnt ; ++n )
							nodes[i].next[n] = nodes[i].next[n + 1] ;
						--nodes[i].ncnt ;
						
						for ( ; p < nodes[k].pcnt ; ++p )
							nodes[k].prev[p] = nodes[k].prev[p + 1] ;
						--nodes[k].pcnt ;
						--j ;
						printf( "Gene merge: %d %d\n", exons[i].end, exons[k].start ) ;
					}
				}
			}
		}
	}*/

	//printf( "Begin\n" ) ;
	//printf( "%d\n", exonCnt ) ;
	transcriptEnd = -1 ;
	transcriptLastExonStart = -1 ;
	outputTxptCnt = 0 ;

	// prepare the mutex for pthread
	pthread_mutex_init( &pickTxptMutex, NULL ) ;
	pthread_mutex_init( &outputTxptMutex, NULL ) ;

	pthread_cond_init( &idleCond, NULL ) ;
	pthread_cond_init( &clearCond, NULL ) ;

	currPickTxptThreadsCnt = 0 ;

	k = -1 ;
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		if ( geneId[i] != -1 )
			continue ;
		if ( exons[i].start > transcriptLastExonStart ) //(transcriptEnd )
		{
			//printf( "%d %d: %d %d\n", transcriptLastExonStart, transcriptEnd, exons[i].start, exons[i].end ) ;
			if ( exons[i].start < transcriptEnd )
			{
				// Adjust the transcriptEnd, if they overlap
				if ( transcriptEnd > exons[i].end ) // The first exon of next gene is totall within the last exon
					transcriptEnd = ( exons[i].start + exons[i].end ) / 2 ;
				else
					transcriptEnd = ( transcriptEnd + exons[i].start ) / 2 ;
			}
			if ( k != -1 && k <= i - 1 )
				SolveGene( chrom, nodes, exons, exonCnt, k, i - 1 ) ;
			k = i ;
		}

		SearchGeneExons( nodes, exons, i, geneIdCnt ) ;
		//printf( "==%d %d\n", i, transcriptEnd  ) ;
		maxAbundance[ geneIdCnt ] = -1 ;
		++geneIdCnt ;
	}
	if ( exonCnt != 0 )
		SolveGene( chrom, nodes, exons, exonCnt, k, i - 1 ) ;

	// Wait for all the threads finishing
	pthread_mutex_lock( &pickTxptMutex ) ;
	if ( currPickTxptThreadsCnt > 0 )
		pthread_cond_wait( &clearCond, &pickTxptMutex ) ;
	pthread_mutex_unlock( &pickTxptMutex ) ;
	/*while ( 1 )
	{
		pthread_mutex_lock( &pickTxptMutex ) ;
		if ( currPickTxptThreadsCnt == 0 )
		{
			pthread_mutex_unlock( &pickTxptMutex ) ;
			break ;
		}
		pthread_mutex_unlock( &pickTxptMutex ) ;
		usleep( 500000 ) ;
	}*/

	// Output the transcripts.
	qsort( outputTranscripts, outputTxptCnt, sizeof( struct _transcript ), CompTranscripts ) ;
	for ( i = 0 ; i < outputTxptCnt ; ++i )
	{
		outputTranscripts[i].transcriptId = geneTranscriptId[ outputTranscripts[i].geneId ] ;
		++geneTranscriptId[ outputTranscripts[i].geneId ] ;
		OutputTranscript( chrom, exons, &outputTranscripts[i] ) ;
		free( outputTranscripts[i].eid ) ;
	}

	/*for ( i = 0 ; i < exonCnt ; ++i )
	{
		//if ( exons[i].pcnt != 0 )
		if ( eidUsedFlag[i] )
			continue ;
			
		if ( exons[i].start > transcriptEnd )
		{
			//printf( "### %d %d %d\n", transcriptStart, transcriptEnd, tcnt ) ;
			PickTranscripts( chrom, exons, transcriptStart, transcriptEnd ) ;
			
			//printf( "=== %d %d %d\n", transcriptStart, transcriptEnd, tcnt ) ;
			euCnt = 0 ;
			transcriptEnd = -1 ;
			atcnt = 0 ;
			transcriptStart = exons[i].start ;
		}
		//printf( "### %d %d %d\n", i, exons[i].start, exons[i].end ) ;
		EnumerateTranscript( chrom, i, visit, 0, exons, exonCnt ) ;
	}
	PickTranscripts( chrom, exons, transcriptStart, transcriptEnd ) ;*/
	transcriptEnd = -1 ;
	transcriptLastExonStart = -1 ;
	pthread_mutex_destroy( &pickTxptMutex ) ;
	pthread_mutex_destroy( &outputTxptMutex ) ;
	pthread_cond_destroy( &idleCond ) ;
	pthread_cond_destroy( &clearCond ) ;

	ClearReadHash( transcriptReads ) ;
	
	for ( i = 0 ; i < exonCnt ; ++i )
	{
		free( nodes[i].next ) ;
		free( nodes[i].prev ) ;
		free( nodes[i].used ) ;
	}
	free( nodes ) ;
}

