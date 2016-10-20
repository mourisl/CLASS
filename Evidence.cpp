#include "Evidence.h"

#define MAX_EVIDENCE 4254587

char evi_buffer[5000] ;
FILE *fpEvidence ;
struct _exon tmpExons[500] ;
int  tmpExonCnt ;
int hashUsed ;

int eviExonTag ;
struct _evidence_hash evidenceHash[MAX_EVIDENCE] ;

void Evidence_Init( char *evidenceFile, struct _evidence **evidence ) 
{
	fpEvidence = NULL ;
	fpEvidence = fopen( evidenceFile, "r" ) ;
	if ( fpEvidence == NULL )
	{
		fprintf( stderr, "Can not find file %s.\n", evidenceFile ) ;
		*evidence = NULL ;
		//return ;
		exit( 1 ) ;
	}
	*evidence = ( struct _evidence *)malloc( sizeof( struct _evidence ) * MAX_EVIDENCE ) ;
	int i ;
	for ( i = 0 ; i < MAX_EVIDENCE ; ++i )
	{
		(*evidence)[i].exons = NULL ;
		(*evidence)[i].ecnt = 0 ;
	}
	evi_buffer[0] = '\0' ;
}


// Test whether |a-b|<=EVIDENCE_MARGIN
bool WithinEvidenceMargin( int a, int b )
{
	return ( a - b <= EVIDENCE_MARGIN && a - b >= - EVIDENCE_MARGIN ) ;
}

int Comp_Evidence( const void *p1, const void *p2 )
{
	struct _evidence *a, *b ;
	a = ( struct _evidence *)p1 ;
	b = ( struct _evidence *)p2 ;
	return a->exons[0].start - b->exons[0].start ;
} 

// Extract evidence of a chromosome
int Evidence_Extract( char *inChrom, struct _evidence *e ) 
{
	char chrom[50], software[50], type[50], score[50], strand[10] ;
	int start, end ;
	char transcript_id[500] = "" ;
	char stmp[500], buffer[500] ;
	bool first = true ;
	char *p ;
	int ret = 0 ;
	int i ;

	tmpExonCnt = 0 ;
	if ( evi_buffer[0] != '\0' )
	{
		sscanf( evi_buffer, "%s %s %s %d %d %s %s", 
			chrom, software, type, &start, &end, score, strand ) ;		
		if ( !strcmp( inChrom, chrom ) && !strcmp( type, "exon" ) )
		{
			first = false ;
			tmpExons[ tmpExonCnt ].start = start ;
			tmpExons[ tmpExonCnt ].end = end ;
			tmpExons[ tmpExonCnt ].strand = strand[0] == '+' ? 1 : 0 ;
			++tmpExonCnt ;

			if ( p = strstr( evi_buffer, "transcript_id" ) )
			{
				p += 15 ;
				for ( i = 0 ; *(p + i) && *(p + i) != '\"' ; ++i )	
					;
				strncpy( stmp, p, i ) ;
				stmp[i] = '\0' ;

				strcpy( transcript_id, stmp ) ;
			}
		}
	}
	while ( fgets( evi_buffer, sizeof( evi_buffer ), fpEvidence ) != NULL )
	{
		sscanf( evi_buffer, "%s %s %s %d %d %s %s", 
			chrom, software, type, &start, &end, score, strand ) ;	
		if ( strcmp( inChrom, chrom ) && first )
			continue ;
		else if ( strcmp( inChrom, chrom ) && !first )
			break ;

		if ( strcmp( type, "exon" ) )
			continue ;
		first = false ;
		bool newTranscript = false ;
		if ( p = strstr( evi_buffer, "transcript_id" ) )
		{
			p += 15 ;
			for ( i = 0 ; *(p + i) && *(p + i) != '\"' ; ++i )	
				;
			strncpy( stmp, p, i ) ;
			stmp[i] = '\0' ;
			
			if ( strcmp( transcript_id, stmp ) )
			{
				// Begin a new transcript
				strcpy( transcript_id, stmp ) ;
				newTranscript = true ;
			}
		}
		
		if ( newTranscript == false )
		{
			if ( tmpExonCnt >= 500 )
				continue ;

			tmpExons[ tmpExonCnt ].start = start ;
			tmpExons[ tmpExonCnt ].end = end ;
			if ( strand[0] == '+' )
				tmpExons[tmpExonCnt].strand = 1 ;
			else if ( strand[0] == '-' )
				tmpExons[tmpExonCnt].strand = 0 ;
			else
				tmpExons[tmpExonCnt].strand = -1 ;
			++tmpExonCnt ;
		}
		else
		{
			if ( tmpExonCnt >= 3 )
			{
				e[ ret ].ecnt = tmpExonCnt ;
				if ( e[ret].exons != NULL )
					free( e[ret].exons ) ;
				e[ret].exons = ( struct _exon *)malloc( sizeof( *( e[ret].exons ) ) * tmpExonCnt ) ;
				for ( i = 0 ; i < tmpExonCnt ; ++i )
					e[ret].exons[i] = tmpExons[i] ;
				/*if ( e[ret].exons[0].start == 67637 && e[ret].exons[ e[ret].ecnt - 1 ].end == 30827626 )
				  {
				  printf( "ERROR %s\n", evi_buffer ) ; exit( 1 ) ;
				  }*/ 	
				++ret ;
			}
			tmpExonCnt = 0 ;

			tmpExons[ tmpExonCnt ].start = start ;
			tmpExons[ tmpExonCnt ].end = end ;
			if ( strand[0] == '+' )
				tmpExons[tmpExonCnt].strand = 1 ;
			else if ( strand[0] == '-' )
				tmpExons[tmpExonCnt].strand = 0 ;
			else
				tmpExons[tmpExonCnt].strand = -1 ;
			++tmpExonCnt ;

		}
	} // end of while

	if ( tmpExonCnt >= 3 )
	{
		e[ ret ].ecnt = tmpExonCnt ;
		if ( e[ret].exons != NULL )
			free( e[ret].exons ) ;
		e[ret].exons = ( struct _exon *)malloc( sizeof( *( e[ret].exons ) ) * tmpExonCnt ) ;
		for ( i = 0 ; i < tmpExonCnt ; ++i )
			e[ret].exons[i] = tmpExons[i] ;
		++ret ;
		tmpExonCnt = 0 ;
	}
	else
		tmpExonCnt = 0 ;
	qsort( e, ret, sizeof( *e ), Comp_Evidence ) ;	
	memset( evidenceHash, -1, sizeof( evidenceHash ) ) ;
	//printf( "Find %d evidences in %s", ret, inChrom ) ;
	//exit( 1 ) ;
	return ret ;
}

void Evidence_ClearHash()
{
	hashUsed = 0 ;
	memset( evidenceHash, -1, sizeof( evidenceHash ) ) ;
}

// Add the evidence into the table. Already make sure it does not exist.
void Evidence_Add( int intron1[2], int intron2[2], int val )
{
	if ( hashUsed > MAX_EVIDENCE / 8 )
		Evidence_ClearHash() ;

	int i, key = 0 ;
	int p[4] = { intron1[0], intron1[1], intron2[0], intron2[1] } ;
	for ( i = 0 ; i < 4 ; ++i )
	{
		key = key * 101 + p[i] ;
		if ( key > MAX_EVIDENCE )
			key %= MAX_EVIDENCE ;
	}
	//if ( key < 0 )
	//	key += MAX_EVIDENCE ;
	while ( evidenceHash[key].p[0] != -1 )
	{
		++key ;
		if ( key > MAX_EVIDENCE )
			key = 0 ;
	}
	++hashUsed ;
	for ( i = 0 ; i < 4 ; ++i )
		evidenceHash[key].p[i] = p[i] ;
	evidenceHash[key].val = val ;
}


int Evidence_Val( int intron1[2], int intron2[2] ) 
{
	int i, key = 0 ;
	int p[4] = { intron1[0], intron1[1], intron2[0], intron2[1] } ;
	for ( i = 0 ; i < 4 ; ++i )
	{
		key = key * 101 + p[i] ;
		if ( key > MAX_EVIDENCE )
			key %= MAX_EVIDENCE ;
	}
	//if ( key < 0 )
	//	key += MAX_EVIDENCE ;

	while ( evidenceHash[key].p[0] != -1 )
	{
		for ( i = 0 ; i < 4 ; ++i )
			if ( evidenceHash[key].p[i] != p[i] )
				break ;
		if ( i >= 4 )
			return evidenceHash[key].val ;
		++key ;
		if ( key > MAX_EVIDENCE )
			key = 0 ;
	}
	return 0 ;
}

int Comp_EvidenceExon( const void *p1, const void *p2 )
{
	struct _exon *a = ( struct _exon * )p1 ;
	struct _exon *b = ( struct _exon * )p2 ;
	if ( a->start != b->start )
		return a->start - b->start ;
	else
		return b->end - a->end ;
}

int ExtractEvidenceExons( struct _exon **evidenceExons, struct _evidence *evidences, int eviCnt )
{
	int maxExonCnt = 0 ;
	eviExonTag = 0 ;
	for ( int i = 0 ; i < eviCnt ; ++i )	
	{
		maxExonCnt += evidences[i].ecnt ;			
	}
	
	if ( maxExonCnt == 0 )
	{
		*evidenceExons = NULL ;
		return 0 ;
	}

	if ( *evidenceExons != NULL )
		free( *evidenceExons ) ;
	*evidenceExons = ( struct _exon * )malloc( maxExonCnt * sizeof( struct _exon ) ) ;
	int cnt = 0 ;
	for ( int i = 0 ; i < eviCnt ; ++i )
	{
		for ( int j = 0 ; j < evidences[i].ecnt ; ++j ) // Only consider the internal exons
		{
			(*evidenceExons)[cnt] = evidences[i].exons[j] ;
			if ( j == 0 )
				(*evidenceExons)[cnt].pcnt = 0 ;
			else
				(*evidenceExons)[cnt].pcnt = 1 ;
			
			if ( j == evidences[i].ecnt - 1 )
				(*evidenceExons)[cnt].ncnt = 0 ;
			else
				(*evidenceExons)[cnt].ncnt = 1 ;
				
			++cnt ;
		}
	}

	if ( cnt == 0 )
	{
		*evidenceExons = NULL ;
		return 0 ;
	}
	
	qsort( *evidenceExons, cnt, sizeof( struct _exon ), Comp_EvidenceExon ) ;		
	int ret = 0 ;
	// Build exon bundles. clutering overlapping exons.
	for ( int i = 1 ; i < cnt ; ++i )
	{
		if ( (*evidenceExons)[i].start > (*evidenceExons)[ret].end )
		{
			++ret ;
			(*evidenceExons)[ret] = (*evidenceExons)[i] ;
			continue ;
		}
		// Below is overlap case
		if ( (*evidenceExons)[i].pcnt == 0 )
			continue ;

		if ( (*evidenceExons)[ret].pcnt == 0 )
		{
			(*evidenceExons)[ret].start = (*evidenceExons)[i].start ;
			(*evidenceExons)[ret].pcnt = (*evidenceExons)[i].pcnt ;
		}

		while ( ret >= 0 && (*evidenceExons)[ret].ncnt == 0 )
		{
			--ret ;
		}
		
		if ( ret < 0 )
		{
			++ret ;
			(*evidenceExons)[ret] = (*evidenceExons)[i] ;
			continue ;
		}

		if ( (*evidenceExons)[i].end > (*evidenceExons)[ret].end )
		{
			(*evidenceExons)[ret].end = (*evidenceExons)[i].end ;
			(*evidenceExons)[ret].ncnt = (*evidenceExons)[i].ncnt ;
		}
	}
	return ret + 1 ;
}


// Already make sure the evidences holds the same position
int IsPosInEvidenceExon( int pos, struct _exon *evidenceExons, int eviExonCnt )
{
	if ( eviExonTag >= eviExonCnt )
		return false ;
	if ( pos >= evidenceExons[eviExonTag].start && pos <= evidenceExons[eviExonTag].end )
		return true ;
	if ( pos < evidenceExons[eviExonTag].start )
		return false ;

	for ( ; eviExonTag < eviExonCnt ; ++eviExonTag )
	{
		if ( evidenceExons[ eviExonTag ].end < pos )
			continue ;
		break ;
	}

	if ( eviExonTag >= eviExonCnt )
		return false ;

	if ( pos >= evidenceExons[eviExonTag].start && pos <= evidenceExons[eviExonTag].end )
		return true ;
	if ( pos < evidenceExons[eviExonTag].start )
		return false ;
	
	// Should never reach here.
	return false ;
}





