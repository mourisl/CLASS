#include "Reads.h"

#define LINE_SIZE 10000
#define HASH_MAX 1000003 

#define OVERLAP_READS_MAX 100000

char line[LINE_SIZE] ;
struct _readTree hash[HASH_MAX] ; // Dynamic allocated

struct _read overlapReads[OVERLAP_READS_MAX] ;
char overlapReadsId[OVERLAP_READS_MAX][ID_LENGTH] ; // Record the reads that overlap more than one region.
int orCnt = 0 ;

extern int READS_LENGTH ;
extern bool VAR_RD_LEN ;

void CigarToString( bam1_core_t *c, uint32_t *in_cigar, char *out_cigar )
{
	int k, op, l ;
	char opcode = '\0' ;

	*out_cigar = '\0' ;
	for ( k = 0 ; k < c->n_cigar ; ++k )
	{
		op = in_cigar[k] & BAM_CIGAR_MASK ;
		l = in_cigar[k] >> BAM_CIGAR_SHIFT ;
		switch (op)
		{
			case BAM_CMATCH: opcode = 'M' ; break ;
			case BAM_CINS: opcode = 'I' ; break ;
			case BAM_CDEL: opcode = 'D' ; break ;
			case BAM_CREF_SKIP: opcode = 'N' ; break ;
			case BAM_CSOFT_CLIP: opcode = 'S' ; break ;
			case BAM_CHARD_CLIP: opcode = 'H' ; break ;
			case BAM_CPAD: opcode = 'P' ; break ;
			default: opcode = 'M' ; break ;
		}
		sprintf( out_cigar + strlen( out_cigar ), "%d%c", l, opcode ) ;
	}
}

struct _readFile OpenReadFile( char *prefix ) 
{
	char buffer[1024] ;
	FILE *fp = NULL ;
	struct _readFile ret ;
	sprintf( buffer, "%s.bam", prefix ) ;
	fp = fopen( buffer, "r" ) ;
	
	if ( fp  )
	{
		ret.fpsam = samopen( buffer, "rb", 0 ) ;
		if ( !ret.fpsam->header ) 
		{
			fprintf( stderr, "Can not open %s.\n", buffer ) ;
			exit( 1 ) ;
		}
		ret.sam = true ;
		fclose( fp ) ;
		return ret ;
	}

	sprintf( buffer, "%s.sam", prefix ) ;
	fp = NULL ;
	fp = fopen( buffer, "r" ) ;
	if ( fp )
	{
		ret.fp = fopen( buffer, "r" ) ;
		ret.sam = false ;
		fclose( fp ) ;
		return ret ;
	}

	fprintf( stderr, "Can not find alignment file %s.bam or %s.sam\n", prefix, prefix ) ;
	exit( 1 ) ;
}

extern int CompInt( const void *p1, const void *p2 ) ;

void GetReadsInfo( struct _readFile file, int &readsLen, int &fragLen, int &fragStd, long long &totalReadCnt ) 
{	
	char readid[ID_LENGTH], chrom[50], mapq[10], cigar[1000], mateChrom[50], matePos[10], seq[1000] ;
	int start, mstart, tlen, tag ;
	int i, j, k, len, cnt ;
	const int LEN_CASE = 1000000, CASE = 1000000  ;
	int lens[LEN_CASE] ;
	int *offset ;
	totalReadCnt = 0 ;
	
	offset = ( int * )malloc( sizeof( int ) * CASE ) ;
	k = 0 ;
	readsLen = -1 ;
	cnt = 0 ;

	if ( file.sam )
	{
		bam1_t *read ;
		read = bam_init1() ;
		for ( i = 0 ; i < CASE && samread( file.fpsam, read ) > 0 ; ++i )
		{
			if ( i < LEN_CASE )
			{
				len = read->core.l_qseq ;
				lens[k] = len ;
				++k ;
				//if ( len > readsLen )
				//	readsLen = len ;
			}

			if ( read->core.tid == read->core.mtid )
				offset[i] = read->core.mpos - read->core.pos ;
			else
				offset[i] = -1 ;

			bam_destroy1( read ) ;
			read = bam_init1() ;
			++totalReadCnt ;
			++cnt ;
		}

		while ( samread( file.fpsam, read ) > 0 )
			++totalReadCnt ;
	}
	else
	{
		for ( i = 0 ; i < CASE && fgets( line, sizeof( line ), file.fp ) ; ++i )
		{
			sscanf( line, "%s %d %s %d %s %s %s %d %d %s", readid, &tag, chrom, &start, mapq, cigar, mateChrom, &mstart, &tlen, seq ) ;
			if ( i < LEN_CASE )
			{
				len = strlen( seq ) ;
				if ( len > readsLen )
					readsLen = len ;
				lens[k] = len ;
				++k ;
			}

			if ( mateChrom[0] == '=' )
			{
				offset[i] = mstart - start ;
			}
			else
				offset[i] = -1 ;
			++totalReadCnt ;
			++cnt ;
		}
		
		while ( fgets( line, sizeof( line ), file.fp ) )
			++totalReadCnt ;
		rewind( file.fp ) ;
	}
	/*int max = 0 ;
	for ( i = 0 ; i < k ; ++i )
		if ( lens[i][0] >= max )
		{
			max = lens[i][0] ;
			readsLen = lens[i][0] ;
		}*/
	qsort( lens, k, sizeof( int ), CompInt ) ;
	if ( !VAR_RD_LEN )
		readsLen = lens[ k / 2 ] ;
	else
		readsLen = lens[k - 1] ;

	long long sum = 0 ;
	long long sumsq = 0 ;
	k = 0 ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( offset[i] < 0 )
			continue ;
		if ( offset[i] + readsLen > 10 * readsLen )
			continue ;
		//sum += offset[i] + readsLen ;
		//sumsq += ( offset[i] + readsLen ) * ( offset[i] + readsLen ) ;
		offset[k] = offset[i] ;
		++k ;
	}
	cnt = k ;
	//printf( "%d\n", offset[10] ) ;
	qsort( offset, k, sizeof( int ), CompInt ) ;
	for ( i = 0 ; i < cnt * 0.9 ; ++i )
	{
		sum += offset[i] + readsLen ;
		sumsq += ( offset[i] + readsLen ) * ( offset[i] + readsLen ) ;
	}
	k = cnt * 0.9 ;
	
	if ( k != 0 )
	{
		fragLen = sum / k ;
		fragStd = (int)sqrt( sumsq / k - fragLen * fragLen ) ;
	}
	else
	{
		fragLen = readsLen ;
		fragStd = 0 ;
	}

	for ( i = 0 ; i < HASH_MAX ; ++i )
	{
		hash[i].timeStamp = 0 ;
		hash[i].left = hash[i].right = NULL ; 
	}
	free( offset ) ;
	//printf( "%d %d %d\n", readsLen, fragLen, fragStd ) ; exit( 1 ) ;
}

/*bool IntervalOverlap( int a0, int a1, int b0, int b1 )
{
	if ( a1 <= b0 || a0 >= b1 )
		return false ;
	return true ;
}*/

// Test whether mate pair a, b compatible with each other
bool AreMatesCompatible( const struct _read &a, const struct _read &b )
{
	int i, j ;
	if ( b.end < a.end )
		return false ;
	if ( a.end <= b.start )
		return true ;
	int sa[20][2], sb[20][2] ; // Converts the splice junctions to sub exons.
	sa[0][0] = a.start ;
	for ( i = 0 ; i < a.scnt ; ++i )
	{
		sa[i][1] = a.splices[i][0] ;
		sa[i + 1][0] = a.splices[i][1] ;
	}
	sa[i][1] = a.end ;
	
	sb[0][0] = b.start ;
	for ( i = 0 ; i < b.scnt ; ++i )
	{
		sb[i][1] = b.splices[i][0] ;
		sb[i + 1][0] = b.splices[i][1] ;
	}
	sb[i][1] = a.end ;

	// Now the mates must overlap.
	// Each intron should not overlap with any exon from its mate.
	for ( i = 0 ; i < a.scnt ; ++i )
	{
		for ( j = 0 ; j <= b.scnt ; ++j )
		{
			if ( !( sa[i + 1][0] <= sb[j][0] || sa[i][1] >= sb[j][1] ) )
				return false ;
		}
	}

	for ( j = 0 ; j < b.scnt ; ++j )
	{
		for ( i = 0 ; i <= a.scnt ; ++i )
		{
			if ( !( sb[j + 1][0] <= sa[i][0] || sb[j][1] >= sa[i][1] ) )
				return false ;
		}
	}
	
	return true ;
}

int SearchReadTree( struct _readTree *p, char *id, int timeStamp )
{
	//printf( "# %s %s %s\n", hash[384].id[0] != 'H' ? "NULL" : hash[384].id , 
	//		hash[384].left != NULL ? "NULL" : hash[384].left->id ,
	//		hash[384].right != NULL ? "NULL" : hash[384].right->id ) ; 
			
	if ( p == NULL || p->timeStamp != timeStamp )
	{
		return -1 ;
	}
	int cmp = strcmp( id, p->id ) ;
	if ( cmp == 0 )
		return p->index ;
	else if ( cmp < 0 )
		return SearchReadTree( p->left, id, timeStamp ) ;
	else
		return SearchReadTree( p->right, id, timeStamp ) ;
}

void InsertReadTree( struct _readTree *p, char *id, int index, int timeStamp )
{
	if ( p->timeStamp != timeStamp )
	{
		strcpy( p->id, id ) ;
		p->index = index ;
		p->timeStamp = timeStamp ;
		return ;
	}

	int cmp = strcmp( id, p->id ) ;
	if ( cmp < 0 )
	{
		if ( p->left == NULL )
		{
			p->left = ( struct _readTree *)malloc( sizeof( struct _readTree ) ) ;
			strcpy( p->left->id, id ) ;
			p->left->index = index ;
			p->left->timeStamp = timeStamp ;
			p->left->left = p->left->right = NULL ;
		}
		else
			InsertReadTree( p->left, id, index, timeStamp ) ;
	}
	else
	{
		if ( p->right == NULL )
		{
			p->right = ( struct _readTree *)malloc( sizeof( struct _readTree ) ) ;
			strcpy( p->right->id, id ) ;
			p->right->index = index ;
			p->right->timeStamp = timeStamp ;
			p->right->left = p->right->right = NULL ;
		}
		else
			InsertReadTree( p->right, id, index, timeStamp ) ;
	}
} 

void ClearReadTree( struct _readTree *p )
{
	if ( p == NULL )
		return ;
	ClearReadTree( p->left ) ;
	ClearReadTree( p->right ) ;
	free( p ) ;	 	
}

// Search and insert
int SearchHash( char *id, int index, int timeStamp )
{
	int i, key = 0 ;
	int ret ;
	for ( i = 0 ; id[i] ; ++i )
		key = ( key * 31 + ( id[i] - '0' ) ) % HASH_MAX ; 
	
	if ( key < 0 )
		key += HASH_MAX ;
	/*
	   while ( 1 )
	   {
	   if ( hash[key] == -1 || !strcmp( hashId[key], id ) )
	   break ;
	   ++key ;
	   if ( key >= HASH_MAX )
	   key = 0 ;
	   }
	   return key ;
	 */
	/*printf( "%s %s %s\n", hash[384].id[0] != 'H' ? "NULL" : hash[384].id , 
			hash[384].left == NULL ? "NULL" : hash[384].left->id ,
			hash[384].right == NULL ? "NULL" : hash[384].right->id ) ; */
	ret = SearchReadTree( &hash[key], id, timeStamp ) ;
	//printf( "%d\n", ret ) ; fflush( stdout ) ;
	if ( ret == -1 )
	{
		//printf( "repeat %d %s %s %d %d\n", key, hash[key].id, id, hash[key].pos, mpos ) ;
		InsertReadTree( &hash[key], id, index, timeStamp ) ;
	}
	return ret ;
}	

int SearchGeneHash( struct _geneRead geneReads, char *id, int index, int timeStamp )
{
	int i, key = 0 ;
	int tag, offset ;
	int ret ;
	for ( i = 0 ; id[i] ; ++i )
		key = ( key * 31 + ( id[i] - '0' ) ) % ( GENE_READ_HASH_MULT * HASH_MAX ) ; 
	if ( key < 0 )
		key += GENE_READ_HASH_MULT * HASH_MAX ;

	tag = key / HASH_MAX ; 
	offset = key - tag * HASH_MAX ;
	
	/*while ( 1 )
	{
		if ( geneReads.hash[tag][offset] == -1 || !strcmp( &geneReads.hashId[tag][offset * ID_LENGTH], id ) )
			break ;
		++offset ;
		if ( offset >= HASH_MAX )
		{
			++tag ;
			offset = 0 ;	
			if ( tag >= GENE_READ_MULT )
				tag = 0 ;
		}
	}*/
	//if ( timeStamp == 119 )
	// printf( "%d %d: %d\n", tag, offset, key ) ;
	ret = SearchReadTree( &geneReads.hash[tag][offset], id, timeStamp ) ;
	if ( ret == -1 )
		InsertReadTree( &geneReads.hash[tag][offset], id, index, timeStamp ) ;
	return ret ;
}

int ExtractReads( struct _readFile file, char *rchrom, int rstart, int rend, struct _read reads[], int extent[2] ) 
{
	char readid[ID_LENGTH], chrom[50], mapq[10], cigar[1000], mateChrom[50] ;
	int start, mstart, flag ; // read start and mate read start
	bool first = true ;
	int num, len ;
	int i, j, k ; 
	int rcnt = 0 ;
	int tmp ;
	bool secondary = false ;
	bam1_t *b = NULL ;
	if ( rend < rstart )
		return 0 ;
	reads[0].scnt = 0 ;
	extent[1] = -1 ;
	static int timeStamp = 0 ;
	++timeStamp ;
	
	//orCnt = 0 ;
	//memset( hash, -1, HASH_MAX * sizeof( hash[0] ) ) ;
	/*for ( i = 0 ; i < HASH_MAX ; ++i )
	{
		if ( hash[i].index != -1 )
		{
			ClearReadTree( hash[i].left ) ;
			ClearReadTree( hash[i].right ) ;
			hash[i].index = -1 ;
		}
	}*/
	//printf( "hi2\n" ) ;

	// Firstly, check the overlaping reads
	for ( i = 0 ; i < orCnt ; ++i )
	{
		if ( overlapReads[i].end < rstart || overlapReads[i].start > rend )
		{
			overlapReads[i].start = -1 ;
			continue ;
		}
		//k = SearchHash( ) ;
		reads[rcnt] = overlapReads[i] ;
		
		if ( overlapReads[i].mateInd == -1 || overlapReads[i].mateInd > rend )
			reads[rcnt].mateInd = -1 ;
		else
		{
			k = SearchHash( overlapReadsId[i], rcnt, timeStamp ) ; // reuse mateInd as mate position
			if ( k == -1 )
			{
				reads[ rcnt ].mateInd = -1 ;
			}			
			else if ( reads[k].start == overlapReads[i].mateInd && AreMatesCompatible( reads[k], reads[rcnt] ) )
			{
				reads[k].mateInd = rcnt ;
				reads[ rcnt ].mateInd = k ; 

				if ( reads[rcnt].strand == -1 )
					reads[ rcnt ].strand = reads[k].strand ;
				if ( reads[ k ].strand == -1 )
					reads[k].strand = reads[rcnt].strand ;
			}
			else
				reads[rcnt].mateInd = -1 ;
		}
		++rcnt ;
		
		if ( rcnt >= MAX_READ )
			return rcnt ;
	}
	// Clean the overlapping reads
	k = 0 ;
	for ( i = 0 ; i < orCnt ; ++i )
	{
		if ( overlapReads[i].start == -1 )
			continue ;
		overlapReads[k] = overlapReads[i] ;
		strcpy( overlapReadsId[k], overlapReadsId[i] ) ;
		++k ;
	}
	orCnt = k ;
	
	while ( 1 )
	{
		//if ( rstart == 9393714 ) //12618480 )
		//	printf( "start: %d\n", rcnt ) ;
		if ( file.sam )
		{
			if ( b )
				bam_destroy1( b ) ;
			b = bam_init1() ;
			if ( samread( file.fpsam, b ) <= 0 )
				break ;
			if ( strlen( bam1_qname(b) ) > ID_LENGTH )
				strcpy( readid, "-1" ) ;
			else
				strcpy( readid, bam1_qname( b ) ) ;
				
			if ( b->core.tid != -1 )
				strcpy( chrom, file.fpsam->header->target_name[ b->core.tid ] ) ;
			else
				strcpy( chrom, "-1" ) ;

			CigarToString( &(b->core), bam1_cigar(b), cigar ) ;
			start = b->core.pos + 1 ;	
			mstart = b->core.mpos + 1 ;
			flag = b->core.flag  ;
			
			if ( b->core.mtid == b->core.tid )
				mateChrom[0] = '=' ;	
			else
				mateChrom[0] = '*' ;	
		
			if ( !VAR_RD_LEN && b->core.l_qseq != READS_LENGTH )
				continue ;
		}
		else
		{
			if ( !fgets( line, sizeof( line ), file.fp ) ) 
				break ;
			int tlen ;
			char seq[1000] ;
			sscanf( line, "%s %d %s %d %s %s %s %d %d %s", readid, &flag, chrom, &start, mapq, cigar, mateChrom, &mstart, &tlen, seq ) ;
			if ( !VAR_RD_LEN && strlen( seq ) != READS_LENGTH )
				continue ;
		}
		//printf( "= %d %s: (%d %d) %d\n", rcnt, readid, rstart, rend, start ) ;	
		tmp = strcmp( chrom, rchrom ) ;		
		if ( first && tmp )
			continue ;
		else if ( !first && tmp )
			break ;
		else if ( !first && start > rend )
			break ;
		
		num = 0 ;
		len = 0 ;
		reads[ rcnt ].scnt = 0 ;
		for ( i = 0 ; cigar[i] ; ++i )
		{
			if ( cigar[i] >= '0' && cigar[i] <= '9' )
				num = num * 10 + cigar[i] - '0' ;
			else if ( cigar[i] == 'I' || cigar[i] == 'S' || cigar[i] == 'H'
				|| cigar[i] == 'P' )
			{
				num = 0 ;
			}
			else if ( cigar[i] == 'N' )
			{
				if ( reads[ rcnt ].scnt >= MAX_READ_SPLICE ) 
					break ;
					
				reads[ rcnt ].splices[ reads[ rcnt ].scnt ][0] = start + len - 1 ;
				reads[ rcnt ].splices[ reads[ rcnt ].scnt ][1] = start + len + num ;
				++reads[ rcnt ].scnt ;
				
				len += num ;
				num = 0 ;
			}
			else
			{
				len += num ;
				num = 0 ;
			}
		}

		//printf( "# %d %s: (%d %d) %d %d\n", rcnt, readid, rstart, rend, start, flag & 4 ) ;		
		if ( start + len - 1 < rstart )
			continue ;
		
		//if (rstart == 9908278 )
		//	printf( "%s %d: %d\n", readid, flag, rcnt ) ;	
		first = false ;
	
		reads[ rcnt ].start = start ;
		reads[ rcnt ].end = start + len - 1 ;

		if ( reads[rcnt].end > extent[1] )
			extent[1] = reads[rcnt].end ;
		
		if ( file.sam )
		{
			if ( bam_aux_get( b, "XS" ) )
			{
				if ( bam_aux2A( bam_aux_get( b, "XS" ) ) == '-' )
					reads[rcnt].strand = 0 ;
				else
					reads[rcnt].strand = 1 ;
			}
			else
			{
				if ( reads[rcnt].scnt > 0 )
					reads[rcnt].strand = 1 ;
				else
					reads[rcnt].strand = -1 ;
			}

			secondary = false ;
			if ( bam_aux_get( b, "NH" ) )
			{
				if ( bam_aux2i( bam_aux_get( b, "NH" ) ) >= 10 && ( flag & 0x100 ) )
					secondary = true ;
			}
		}
		else
		{
			if ( strstr( line, "XS:A:-" ) )
				reads[ rcnt ].strand = 0 ;
			else if ( reads[ rcnt ].scnt > 0 )
				reads[ rcnt ].strand = 1 ;
			else
				reads[ rcnt ].strand = -1 ;
		}
		
		if ( ( flag & 4 ) || secondary )
			--rcnt ;	
		else if ( mateChrom[0] != '=' || mstart > rend )
			reads[ rcnt ].mateInd = -1 ;
		else
		{
			k = SearchHash( readid, rcnt, timeStamp ) ;

			if ( k == -1 )
			{
				reads[ rcnt ].mateInd = -1 ;
			}			
			else if ( reads[k].mateInd == -1 && reads[k].start == mstart && AreMatesCompatible( reads[k], reads[rcnt] ) )
			{
				reads[k].mateInd = rcnt ;
				reads[ rcnt ].mateInd = k ; 

				if ( reads[rcnt].strand == -1 )
					reads[ rcnt ].strand = reads[k].strand ;
				if ( reads[ k ].strand == -1 )
					reads[k].strand = reads[rcnt].strand ;
			}
			else
				reads[k].mateInd = -1 ;
		}
		if ( (flag & 4) == 0 && !secondary && reads[rcnt].end > rend + 6 && orCnt < OVERLAP_READS_MAX )
		{
			overlapReads[orCnt] = reads[rcnt] ;
			if ( mateChrom[0] != '=' || mstart > rend )
				overlapReads[orCnt].mateInd = -1 ;
			else
				overlapReads[orCnt].mateInd = mstart ;
			strcpy( overlapReadsId[orCnt], readid ) ;
			++orCnt ;
		}

		++rcnt ;	
		
		if ( rcnt >= MAX_READ )
			break ;			
	}
	extent[0] = reads[0].start ;
	return rcnt ;	
}

void InitGeneReads( struct _geneRead *geneReads )
{
	int i, j ;
	for ( i = 0 ; i < GENE_READ_MULT ; ++i )
	{
		geneReads->reads[i] = ( struct _read *)malloc( sizeof( struct _read ) * MAX_READ ) ;
	}

	for ( i = 0 ; i < GENE_READ_HASH_MULT ; ++i )
	{
		geneReads->hash[i] = ( struct _readTree * )malloc( sizeof( struct _readTree ) * HASH_MAX ) ;
		for ( j = 0 ; j < HASH_MAX ; ++j )
		{
			geneReads->hash[i][j].timeStamp = 0 ;
			geneReads->hash[i][j].left = geneReads->hash[i][j].right = NULL ;
		}
	}	
}


int ExtractGeneReads( struct _readFile file, char *rchrom, int rstart, int rend, struct _geneRead geneReads, int extent[2] )
{
	char readid[ID_LENGTH], samtag[1000], chrom[50], mapq[10], cigar[1000], mateChrom[10] ;
	int flag, start, mstart ; // read start and mate read start
	bool first = true ;
	int num, len ;
	int tag, offset ;
	int i, j, k ; 
	int rcnt = 0 ;
	int tmp ;
	bool secondary = false ;

	bam1_t *b = NULL ;
	if ( rend < rstart )
		return 0 ;
	static int timeStamp = 0 ;
	++timeStamp ;
	//reads[0].scnt = 0 ;
	extent[1] = -1 ;
	/*for ( i = 0 ; i < GENE_READ_MULT ; ++i )
	{
		memset( geneReads.hash[i], -1, sizeof( int ) * HASH_MAX ) ;
	}*/
	//printf( "%d %d\n", rstart, timeStamp ) ;

	while ( 1 )
	{
		//if ( rstart == 9393714 ) //12618480 )
		//	printf( "start: %d\n", rcnt ) ;
		if ( file.sam )
		{
			if ( b )
				bam_destroy1( b ) ;
			b = bam_init1() ;
			if ( samread( file.fpsam, b ) <= 0 )
				break ;
			strcpy( readid, bam1_qname( b ) ) ;
			if ( b->core.tid >= 0 )
				strcpy( chrom, file.fpsam->header->target_name[ b->core.tid ] ) ;
			else
				strcpy( chrom, "-1" ) ;
			CigarToString( &(b->core), bam1_cigar(b), cigar ) ;
			start = b->core.pos + 1 ;	
			mstart = b->core.mpos + 1 ;
			flag = b->core.flag  ;
			if ( b->core.mtid == b->core.tid )
				mateChrom[0] = '=' ;		
			
			if ( !VAR_RD_LEN && b->core.l_qseq != READS_LENGTH )
				continue ;
		}
		else
		{
			if ( !fgets( line, sizeof( line ), file.fp ) ) 
				break ;
			int tlen ;
			char seq[1000] ;
			sscanf( line, "%s %d %s %d %s %s %s %d %d %s", readid, &flag, chrom, &start, mapq, cigar, mateChrom, &mstart, &tlen, seq ) ;
			if ( !VAR_RD_LEN && strlen( seq ) != READS_LENGTH )
				continue ;
		}

		tmp = strcmp( chrom, rchrom ) ;		
		if ( first && tmp )
			continue ;
		else if ( !first && tmp )
			break ;
		else if ( !first && start > rend )
			break ;
		if ( rcnt >= MAX_READ * GENE_READ_MULT )
			continue ;	
		num = 0 ;
		len = 0 ;
		tag = rcnt / MAX_READ ;
		offset = rcnt - tag * MAX_READ ;

		geneReads.reads[tag][offset].scnt = 0 ;
		for ( i = 0 ; cigar[i] ; ++i )
		{
			if ( cigar[i] >= '0' && cigar[i] <= '9' )
				num = num * 10 + cigar[i] - '0' ;
			else if ( cigar[i] == 'I' || cigar[i] == 'S' || cigar[i] == 'H'
				|| cigar[i] == 'P' )
			{
				num = 0 ;
			}
			else if ( cigar[i] == 'N' )
			{
				if ( geneReads.reads[tag][offset].scnt >= MAX_READ_SPLICE ) 
					break ;
				geneReads.reads[tag][offset].splices[ geneReads.reads[tag][offset].scnt ][0] = start + len - 1 ;
				geneReads.reads[tag][offset].splices[ geneReads.reads[tag][offset].scnt ][1] = start + len + num ;
				++geneReads.reads[tag][offset].scnt ;
				
				len += num ;
				num = 0 ;
			}
			else
			{
				len += num ;
				num = 0 ;
			}
		}

		if ( start + len - 1 < rstart )
			continue ;
		first = false ;
		geneReads.reads[tag][offset].start = start ;
		geneReads.reads[tag][offset].end = start + len - 1 ;
		

		if ( geneReads.reads[tag][offset].end > extent[1] )
			extent[1] = geneReads.reads[tag][offset].end ;
		// TODO: Introduce the undecided strand.
		if ( file.sam )
		{
			if ( bam_aux_get( b, "XS" ) )
			{
				if ( bam_aux2A( bam_aux_get( b, "XS" ) ) == '-' )
					geneReads.reads[tag][offset].strand = 0 ;
				else
					geneReads.reads[tag][offset].strand = 1 ;
			}
			else
			{
				if ( geneReads.reads[tag][offset].scnt > 0 )
					geneReads.reads[tag][offset].strand = 1 ;
				else
					geneReads.reads[tag][offset].strand = -1 ;
			}
			secondary = false ;
			geneReads.reads[tag][offset].unique = true ;
			if ( bam_aux_get( b, "NH" ) )
			{
				if ( bam_aux2i( bam_aux_get( b, "NH" ) ) >= 10 && ( flag & 0x100 ) )
				{
					secondary = true ;
				}

				if ( bam_aux2i( bam_aux_get( b, "NH" ) ) > 1 )
					geneReads.reads[tag][offset].unique = false ;
			}
		}
		else
		{
			if ( strstr( line, "XS:A:-" ) )
				geneReads.reads[tag][offset].strand = 0 ;
			else if ( geneReads.reads[tag][offset].scnt > 0 )
				geneReads.reads[tag][offset].strand = 1 ;
			else
				geneReads.reads[tag][offset].strand = -1 ;
		}	
		//if ( len >= 50000 )
		//{
		//	printf( "%s\n", readid ) ;
		//}	
		if ( ( flag & 4 ) || secondary )
			--rcnt ;	
		else if ( mateChrom[0] != '=' || mstart > rend )
			geneReads.reads[tag][offset].mateInd = -1 ;
		else
		{
			//if ( rstart == 106068489 )
			//	printf( "\thi1\n" ) ;
			k = SearchGeneHash( geneReads, readid, rcnt, timeStamp ) ;
			//if ( rstart == 9393714) //12618480 )
			//	printf( "%d\n", rcnt ) ;
			//k = -1 ;
			if ( k == -1 )
			{
				geneReads.reads[tag][offset].mateInd = -1 ;
			}			
			else//if ( geneReads.reads[htag][hoffset].start == mstart )
			{
				int tmptag, tmpoffset ;
				tmptag = k / MAX_READ ;
				tmpoffset = k - MAX_READ * tmptag ;
				if ( geneReads.reads[tmptag][tmpoffset].mateInd == -1 && geneReads.reads[tmptag][tmpoffset].start == mstart &&
					AreMatesCompatible( geneReads.reads[tmptag][tmpoffset], geneReads.reads[tag][offset] ) )
				{
					//if ( start - mstart > 50000 )
					//{
					//	printf( "%d\n", flag ) ;
					//	printf( "%s\n", readid ) ;
					//}
					geneReads.reads[tmptag][tmpoffset].mateInd = rcnt ;
					geneReads.reads[tag][offset].mateInd = k ; 

					if ( geneReads.reads[tag][offset].strand == -1 )
						geneReads.reads[tag][offset].strand = geneReads.reads[tmptag][tmpoffset].strand ;
					if ( geneReads.reads[tmptag][tmpoffset].strand == -1 )
						geneReads.reads[tmptag][tmpoffset].strand = geneReads.reads[tag][offset].strand ;
				}
				else
					geneReads.reads[tag][offset].mateInd = -1 ;
			}
		}
		++rcnt ;		
	}
	extent[0] = geneReads.reads[0][0].start ;
	return rcnt ;	
}

void ClearReadHash( struct _geneRead &geneReads )
{
	int i, j ;
	for ( i = 0 ; i < GENE_READ_HASH_MULT ; ++i )
	{
		for ( j = 0 ; j < HASH_MAX ; ++j )
		{
			ClearReadTree( geneReads.hash[i][j].left ) ;
			ClearReadTree( geneReads.hash[i][j].right ) ;
			geneReads.hash[i][j].left = geneReads.hash[i][j].right = NULL ;
		}
	}

	for ( j = 0 ; j < HASH_MAX ; ++j )
	{
		ClearReadTree( hash[j].left ) ;	
		ClearReadTree( hash[j].right ) ;
		hash[j].left = hash[j].right = NULL ;
	}
	orCnt = 0 ;
}


