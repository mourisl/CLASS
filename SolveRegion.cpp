/**
	Given a region and its coverage, splice site information, find the best possible regionExons that can explain the information by Linear Programming method.
*/

#include "SolveRegion.h" 


extern int READS_LENGTH ;
extern int FRAG_LENGTH, FRAG_STD ;

double psum[MAX_LENGTH] ;
char buffer[MAX_LENGTH], buffer2[MAX_LENGTH] ;
int OFFSET ;

struct _point inputStart[MAX_POINT], inputEnd[MAX_POINT] ;
struct _spliceIndex spliceIndices[MAX_POINT * MAX_POINT] ;

//struct _point start[MAX_POINT], end[MAX_POINT] ; // The points after merging

struct _point *allStart, *allEnd ;
int allStartCnt, allEndCnt ;

//int point[MAX_POINT] ;
//double interval[MAX_POINT] ;
// The enumerated regionExons.
struct _enumExon
{
	int start, end ;
	int startInd, endInd ;
	int strand ;
	
	bool inEvidence ;
} ;
struct _enumExon regionExons[MAX_POINT * MAX_POINT] ;
bool exonHasSoft[MAX_POINT * MAX_POINT][2] ; // Does a exon have soft boundary. 0-left, 1-right?

//int ecnt = 0 ;
//bool eused[MAX_POINT * MAX_POINT] ;
//int startCover[MAX_POINT][3], endCover[MAX_POINT][3] ;
//int intervalCnt, startCnt, endCnt ;

struct _exonfrag
{
	int eid ; // exon id
	int length ; 
	int previd ; // exon frag index of the previous fragment of the same exon.

	int type ; // 0-normal, 1-at the end, which can have higher difference between its neighbor.
} ;

//struct _exonfrag exonfrag[MAX_POINT * MAX_POINT] ;
//int efCnt ;

// Data structure used to find the constraint between splice junctions
// TODO: Can be more efficient
//bool graph[MAX_POINT][MAX_POINT] ; 
//bool visited[MAX_POINT][2] ;
//bool beyond[MAX_POINT][2][2] ; // Whether one side of the exon is beyond the region. 0-left end, 1-right end. The last "2" is the strand.
//int leftEf[MAX_POINT], rightEf[MAX_POINT] ;
//int leftEfCnt, rightEfCnt ;
//int efid[MAX_POINT][2] ; // The first and last ef index of an exon. 0-left end, 1-right end

// Data structure for traversing "exon path" when testing contraint points.
// Means: is there a path between those two exon?
//int exonPath[MAX_POINT][MAX_POINT] ;

//int exonDist[2][MAX_POINT][MAX_POINT] ; // The minimal and maximal distance between two regionExons. 0-min, 1-max

// Constraints from the mate pair reads.
struct _constraintPoint
{
	int startPos, endPos ;
	int strand ;
	int type ; // 0-must belong to the same exon, 1-they are connected through several regionExons. 2-the Pos is the index of an interval	
		   // 3-the chunk [startPos, endPos] must be from the same exon. 4-reachability between two mate pairs.	
	int support ; // # of supporting reads.
	int otherPos[2] ;
	int precise[2] ;
	bool valid ;
} ;

struct _constraintPoint constraintPoints[MAX_POINT * MAX_POINT] ;
//int cpCnt = 0 ;

// For returning
//bool bestEused[MAX_POINT * MAX_POINT] ;
//double bestScore  ;
//int bestEcnt ;

//int mainStrand ; // What is the main strand in this region

extern FILE *fpPolyA ;

// The attributes for best combination.
struct _bestCombAttr
{
	int bestEcnt ;
	int bestEvidenceExonCnt ;
	int bestLongExonCnt ;
	int bestBalancedExonCnt ; 
} ;

struct _solveRegionData 
{
	bool *graph ;
	bool *visited ;
	bool *beyond ; 
	struct _exonfrag *exonfrag ;
	int *efid ;
	int *leftEf, *rightEf ;
	struct _enumExon *regionExons ;
	bool *eused ;
	bool *exonHasSoft ;
	int *exonPath ;
	int *exonDist ;
	struct _constraintPoint *constraintPoints ;
	struct _point *start, *end ;
	int *point ;
	double *interval ;
	int *startCover, *endCover ;
	struct _spliceIndex *spliceIndices ;

	struct _bestCombAttr bestCombAttr ;
} ;

struct _pthreadArgSolveRegion
{
	struct _enumExon *regionExons ;
	int ecnt ;

	bool *exonHasSoft ;
	int OFFSET ;
	int mainStrand ;

	struct _point *start, *end ;
	int *startCover, *endCover ;
	int *point ;
	double *interval ;
	int startCnt, endCnt, intervalCnt ;

	struct _spliceIndex *spliceIndices ;
	int siCnt ;

	struct _constraintPoint *constraintPoints ;
	int cpCnt ;

	struct _exon *allExons ;
	int *exonCnt ;

	int readCnt ;
} ;

extern pthread_attr_t pthreadAttr ;
extern int NUM_OF_THREADS ;

int currSolveRegionThreadsCnt ;
pthread_mutex_t solveRegionMutex ;
pthread_mutex_t allExonsMutex ;
pthread_cond_t idleSolveRegionCond ;
pthread_cond_t clearSolveRegionCond ;

/**
  Test whether s1 and s2 can belong to the same strand
*/
bool IsSameStrand( int s1, int s2 )
{
	if ( s1 >= 0 && s2 >= 0 )
		return s1 == s2 ;
	if ( s1 == -1 || s2 == -1 )
		return true ;
	if ( ( s1 == -2 && s2 == 1 ) 
		|| ( s1 == 1 && s2 == -2 ) )
		return true ;
	if ( ( s1 == -3 && s2 == 0 ) 
		|| ( s1 == 0 && s2 == -3 ) )
		return true ;
	return false ;
}

/**
Suppose a, b are the reads supports the junctions.
Test whether a significantly great b.
*/
bool SignificantGreater( int a, int b )
{
	/*if ( a <= b )
		return false ;
	if ( b <= 20 && a > b + 100 )
		return true ;
	if ( b <= 50 && a > b + 200 )
		return true ;
	if ( b <= 100 && a > b + 400 )
		return true ;
	if ( a > 5 * b )
		return true ;*/
	if ( a > 100 * b )
		return true ;
	if ( b <= 100 && a - 6 * sqrt( a ) > b + 6 * sqrt( b ) )
		return true ;
	return false ;
}


/**
  Find the representative of the set that the point is merged to.
*/
int GetMergeFather( struct _point points[], int tag )
{
	if ( points[tag].merge != tag )
		return points[tag].merge = GetMergeFather( points, points[tag].merge ) ;
	return tag ;
}

/**
	Recursively determine the exon fragments on the left side and right side of a "complete" splice.
	tag: the current exon.
	direction: 0-from left to right(right end of an exon), 1-from right to left(left end of an exon)
	@return : false if an exon's boundary has its pair beyond the region.
*/
bool FindSpliceExon( int tag, int direction, int efCnt, int ecnt, 
	int &leftEfCnt, int &rightEfCnt, struct _solveRegionData &data )
{
	bool *eused = data.eused ;
	bool *graph = data.graph ;
	bool *visited = data.visited ;
	bool *beyond = data.beyond ;
	struct _enumExon *regionExons = data.regionExons ;
	int *efid = data.efid ;
	int *leftEf = data.leftEf ;
	int *rightEf = data.rightEf ;

	if ( visited[tag * 2 + direction] || !eused[tag] )
		return true ;
	visited[tag * 2 + direction] = true ;	
	int i ;
	if ( direction == 0 )
	{
		if ( regionExons[tag].strand >= 0 && beyond[ tag * 4 + 2 + regionExons[tag].strand ] )
			return false ;

		leftEf[ leftEfCnt ] = efid[tag * 2 + 1] ;
		++leftEfCnt ;
		
		for ( i = tag + 1 ; i < ecnt ; ++i )
		{
			if ( graph[tag * ecnt + i] )
				if ( !FindSpliceExon( i, 1, efCnt, ecnt, leftEfCnt, rightEfCnt, data ) )
					return false ;
		}
	}
	else
	{
		if ( regionExons[tag].strand >= 0 && beyond[ tag * 4 + 2 * 0 + regionExons[tag].strand ] )
			return false ;

		rightEf[ rightEfCnt ] = efid[tag * 2 + 0] ;
		++rightEfCnt ;
		for ( i = 0 ; i < tag ; ++i )
			if ( graph[i * ecnt + tag] )
				if ( !FindSpliceExon( i, 1, efCnt, ecnt, leftEfCnt, rightEfCnt, data ) )
					return false ;
	}
	return true ;
}

double LP( int ecnt, int intervalCnt, int startCnt, int endCnt, struct _solveRegionData &data, bool useExonFrag = true )
{
	int i, j, k ;
	int p ;
	int a, b ;
	int miCnt ;
	int Ncol ;
	int *colno ;
	double *row ;
	lprec *lp ;	
	bool right = false ; // right combination?	
	int spliceExonfrag[2][MAX_POINT] ; // The indices of the exonfrags that spanning the splice.
	int spliceSupport[2], scCnt[2] ;
	//long long intervalCover[MAX_POINT] ; // Use bit to represent which regionExons covered this interval.
	int *cid ; // The current id(exonfrag id) for corresponded exon.
	int varCnt ; // The number of variables in the LP that will be used.
	double ret ;	
	int efCnt ;
	double sum = 0 ;
	struct _enumExon *regionExons = data.regionExons ;
	bool *exonHasSoft = data.exonHasSoft ;
	bool *eused = data.eused ;
	struct _exonfrag *exonfrag = data.exonfrag ;
	int *efid = data.efid ;
	int *point = data.point ;
	double *interval = data.interval ;
	bool *visited = data.visited ;
	int *leftEf = data.leftEf ;
	int *rightEf = data.rightEf ;
	struct _point *start = data.start ;
	struct _point *end = data.end ;

	cid = (int *)malloc( sizeof( int ) * ecnt ) ;
	efCnt = 0 ;
	memset( cid, -1, sizeof( int ) * ecnt ) ;

	/*for ( i = 0 ; i <= intervalCnt ; ++i )
	{
		printf( "%d ", point[i] ) ;
	}
	printf( "\n" ) ;*/
	/*for ( i = 0 ; i < ecnt ; ++i )
		printf( "%d", eused[i] ) ;
	printf( "\n" ) ;*/
	for ( i = 0 ; i < intervalCnt ; ++i )
		sum += interval[i] ;
	
	if ( useExonFrag )
	{	
		for ( i = 0 ; i < intervalCnt ; ++i )
		{
			//k = 0 ;
			//intervalCover[i] = 0 ;
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( eused[j] && regionExons[j].start <= point[i + 1] && regionExons[j].end > point[i]  )
				{
					//colno[k] = j + 1 ;
					//row[k] = 1 ;
					/*a = regionExons[j].start > point[i] ? regionExons[j].start : point[i] ;
					  b = regionExons[j].end < point[i + 1] ? regionExons[j].end : point[i + 1] ;
					  row[k] = ( b - a + 1 ) / (double)( regionExons[j].end - regionExons[j].start + 1 ) ;  */
					//printf( "### %d %d %lf\n", a, b, row[k] ) ;				
					//++k ;

					exonfrag[ efCnt ].eid = j ;
					exonfrag[ efCnt ].length = point[i + 1] - point[i] ;
					exonfrag[ efCnt ].previd = cid[j] ;
					
					if ( ( regionExons[j].start == point[i] + 1 && exonHasSoft[j * 2 + 0] ) || 
						( regionExons[j].end == point[i + 1] && exonHasSoft[j * 2 + 1] ) )
					{
						exonfrag[ efCnt ].type = 1 ;
					}
					else
						exonfrag[ efCnt ].type = 0 ;

					if ( cid[j] == -1 )
						efid[j * 2 + 0] = efCnt ;
					cid[j] = efCnt ;
					efid[j * 2 + 1] = efCnt ;
					++efCnt ;
				} 
			}
		}
	}
	else
	{
		for ( i = 0 ; i < ecnt ; ++i )
		{
			exonfrag[ efCnt ].eid = i ;
			exonfrag[ efCnt ].length = regionExons[i].end - regionExons[i].start + 1 ;
			exonfrag[ efCnt ].previd = -1 ;
			efid[i * 2 + 0] = efid[i * 2 + 1] = efCnt ; //i ;
			cid[i] = efCnt ;
			++efCnt ;
		}
	}

	Ncol = efCnt + intervalCnt + efCnt + intervalCnt ; // The last two is for the slack variable.
	colno = (int *)malloc( Ncol * sizeof( *colno ) ) ;	
	row = (double *)malloc( Ncol * sizeof( *row ) ) ;	
	lp = make_lp( 0, Ncol ) ;
	set_add_rowmode(lp, TRUE) ;
	varCnt = efCnt + intervalCnt ;

	//printf( "#### %d\n", efCnt ) ;
	p = 0 ;
	for ( i = 0 ; i < intervalCnt ; ++i )
	{	
		k = 0 ;
		for ( j = 0 ; j < ecnt ; ++j )
		{
			if ( eused[j] && regionExons[j].start <= point[i + 1] && regionExons[j].end > point[i]  )
			{
				if ( useExonFrag )
					colno[k] = p + 1 ;
				else
					colno[k] = j + 1 ;
				row[k] = 1 ;
				++k ;
				++p ;
			}
		}

		colno[k] = efCnt + i + 1 ;
		row[k] = -1 ; 
		//printf( "###%d %d %d\n", i, mergeInterval[i].read, mergeInterval[i].length) ;
		
		add_constraintex(lp, k + 1, row, colno, LE, (double)interval[i] / ( point[i + 1] - point[i] ) ) ;		

		row[k] = 1 ;
		add_constraintex(lp, k + 1, row, colno, GE, (double)interval[i] / ( point[i + 1] - point[i] ) ) ;
	}
	
	/**
	 Constraint from the splice junctions.
	*/
	// For the start points
	p = 0 ;
	for ( i = 0 ; i < 0/*intervalCnt*/ ; ++i )
	{
		scCnt[0] = scCnt[1] = 0 ;
		for ( j = 0 ; j < ecnt ; ++j )
		{
			if ( eused[j] && regionExons[j].start <= point[i + 1] && regionExons[j].end > point[i]  )
			{
				if ( regionExons[j].start == point[i] + 1 )
				{
					a = regionExons[j].strand ;
					spliceExonfrag[a][ scCnt[a] ] = p ;
					++scCnt[a] ;
					spliceSupport[a] = start[ regionExons[j].startInd ].support ;
				}
				++p ;
			}
		}

		for ( k = 0 ; k <= 1 ; ++k )
		{
			if ( !spliceSupport[k] || !scCnt[k] )
				continue ;
			for ( j = 0 ; j < scCnt[k] ; ++j )
			{
				colno[j] = spliceExonfrag[k][j] + 1 ;
				row[j] = 1 ;	
			} 
			//add_constraintex( lp, scCnt[k], row, colno, GE, sqrt( spliceSupport[k] * 0.1 )) ; //- SPLICE_ALPHA * ( point[i + 1] - point[i] ) / 2 ) ;
			add_constraintex( lp, scCnt[k], row, colno, LE, sqrt( spliceSupport[k] * 10 ) ) ; //+ SPLICE_ALPHA * ( point[i + 1] - point[i] ) / 2 ) ;
		}
	}	
	
	// For the end points
	p = 0 ;
	for ( i = 0 ; i < 0/*intervalCnt*/ ; ++i )
	{
		scCnt[0] = scCnt[1] = 0 ;
		for ( j = 0 ; j < ecnt ; ++j )
		{
			if ( eused[j] && regionExons[j].start <= point[i + 1] && regionExons[j].end > point[i]  )
			{
				if ( regionExons[j].end == point[i + 1] )
				{
					a = regionExons[j].strand  ;
					spliceExonfrag[a][ scCnt[a] ] = p ;
					++scCnt[a] ;
					spliceSupport[a] = end[ regionExons[j].endInd ].support ;
				}
				++p ;
			}
		}

		for ( k = 0 ; k <= 1 ; ++k )
		{
			if ( !spliceSupport[k] || !scCnt[k] )
				continue ;
			for ( j = 0 ; j < scCnt[k] ; ++j )
			{
				colno[j] = spliceExonfrag[k][j] + 1 ;
				row[j] = 1 ;	
			} 
			//add_constraintex( lp, scCnt[k], row, colno, GE, spliceSupport[k] * 0.1 ) ; //- SPLICE_ALPHA * ( point[i + 1] - point[i] ) / 2 ) ;
			//add_constraintex( lp, scCnt[k], row, colno, LE, spliceSupport[k] * 10 ) ; //+ SPLICE_ALPHA * ( point[i + 1] - point[i] ) / 2 ) ;
		}
	}		
	
	/**
	 Every exonfrag should be greater than 1. ( TODO: ? )
	*/
	for ( j = 0 ; j < efCnt ; ++j )
	{
		colno[0] = j + 1 ;
		row[0] = 1 ;
		add_constraintex( lp, 1, row, colno, GE, 1 ) ; //( regionExons[j].end - regionExons[j].start + 1 ) ) ;
	}	
	
	/**
	 Continuity: The relationship between adjacent exonfrag.
	*/
	for ( j = 0 ; j < ecnt ; ++j )
	{
		if ( !eused[j] )
			continue ;
		p = cid[j] ;
		while ( 1 )
		{
			k = exonfrag[p].previd ; 
			if ( k == -1 )
				break ;
			colno[0] = k + 1 ;
			colno[1] = p + 1 ;

			row[0] = 1 ;
			row[1] = -1 ;
		
			/*if ( regionExons[j].strand == 0 )
				row[0] = pow( 1.002, ( exonfrag[p].length + exonfrag[k].length ) / 2 ) ;
			else
				row[1] = -pow( 1.002, ( exonfrag[p].length + exonfrag[k].length ) / 2 ) ;*/

			/*colno[2] = varCnt + 1 ;
			row[2] = 1 ;	
			add_constraintex( lp, 3, row, colno, GE, -( exonfrag[p].length + exonfrag[k].length ) * (double)ALPHA / 2 ) ;

			//row[0] = ALPHA ;
			//row[1] = -1 ;
			row[2] = -1 ;	
			add_constraintex( lp, 3, row, colno, LE, ( exonfrag[p].length + exonfrag[k].length ) * (double)ALPHA / 2 ) ;
			++varCnt ;*/
			//row[0] = 1 ;
			//row[1] = -1 ;
			
			//add_constraintex( lp, 2, row, colno, EQ, 0 ) ;
		
			
			colno[2] = varCnt + 1 ;
			row[2] = ( exonfrag[p].length + exonfrag[k].length ) / 2 ;
			//add_constraintex( lp, 3, row, colno, GE, 0 ) ;
			//if ( exonfrag[k].type == 0 && exonfrag[p].type == 0 )
			//{
			double tmp = ( exonfrag[p].length + exonfrag[k].length ) * (double)ALPHA / 2 ;
			if ( tmp > 10 )
				tmp = 10 ;
			add_constraintex( lp, 3, row, colno, GE, -tmp ) ;
			//}
			//else
			//{
				
			//	add_constraintex( lp, 3, row, colno, GE, -( exonfrag[p].length + exonfrag[k].length ) * (double)1 / 2 ) ;
			//}

			row[2] = -row[2] ;
			//add_constraintex( lp, 3, row, colno, LE, 0 ) ;
			//if ( exonfrag[k].type == 0 && exonfrag[p].type == 0 )
			//{
				add_constraintex( lp, 3, row, colno, LE, tmp ) ;
			//}
			//else
			//{
			//	add_constraintex( lp, 3, row, colno, LE, ( exonfrag[p].length + exonfrag[k].length ) * (double)1 / 2 ) ;
			//}
			
			++varCnt ;
			p = k ;
		}
	}
	
	/**
		The relationship from splice junctions that the regionExons on the left side should euqal to right side
	*/
	memset( visited, false, sizeof( visited[0] ) * ecnt * 2 ) ;
	for ( i = 0 ; i < 0/*ecnt*/ ; ++i )
	{
		if ( !eused[i] )
			continue ;
		int leftEfCnt = 0, rightEfCnt = 0 ;
		//printf( "## %d %d %d\n", i, visited[3][0], graph[3][5] ) ;
		if ( !FindSpliceExon( i, 0, efCnt, ecnt, leftEfCnt, rightEfCnt, data ) )
			continue ;
		if ( rightEfCnt )
		{
		//	printf( "## %d %d\n", leftEf[0], right) ;
			for ( i = 0 ; i < leftEfCnt ; ++i )
			{
				colno[i] = leftEf[i] + 1 ;
				row[i] = 1 ;
			}	
			
			for ( ; i < leftEfCnt + rightEfCnt ; ++i )
			{
				colno[i] = rightEf[i - leftEfCnt] + 1 ;
				row[i] = -1 ;
			}
			
			/*colno[i] = varCnt + 1 ;
			row[i] = 1 ;	
			add_constraintex( lp, leftEfCnt + rightEfCnt + 1, row, colno, GE, 
						-( exonfrag[ leftEf[0] ].length + exonfrag[ rightEf[0] ].length ) * (double)ALPHA / 2 ) ;
			row[i] = -1 ;
			add_constraintex( lp, leftEfCnt + rightEfCnt + 1, row, colno, LE, 
		    			( exonfrag[ leftEf[0] ].length + exonfrag[ rightEf[0] ].length ) * (double)ALPHA / 2 ) ;						
			++varCnt ;*/

			colno[i] = varCnt + 1 ;
			row[i] = ( exonfrag[ leftEf[0] ].length + exonfrag[ rightEf[0] ].length ) / 2 ;
			add_constraintex( lp, leftEfCnt + rightEfCnt + 1, row, colno, GE, 
						-( exonfrag[ leftEf[0] ].length + exonfrag[ rightEf[0] ].length ) * (double)ALPHA / 2 ) ;
			row[i] = -row[i] ;
			//add_constraintex( lp, leftEfCnt + rightEfCnt + 1, row, colno, LE, 0 ) ;
			add_constraintex( lp, leftEfCnt + rightEfCnt + 1, row, colno, LE, 
		    			( exonfrag[ leftEf[0] ].length + exonfrag[ rightEf[0] ].length ) * (double)ALPHA / 2 ) ;						
			++varCnt ;
		}
	}
	//if ( OFFSET == 156914811 )
	//	write_LP( lp, stdout ) ;
	/**
	 The sum of them should be equal to the total number of reads.
	*/
	k = 0 ;
	for ( j = 0 ; j < efCnt ; ++j )
	{
		colno[j] = j + 1 ;
		row[j] = exonfrag[j].length ;
	}
	add_constraintex( lp, efCnt, row, colno, GE, sum - efCnt ) ;
	add_constraintex( lp, efCnt, row, colno, LE, sum + efCnt ) ;

	/*for ( i = 0 ; i < ecnt ; ++i )
		set_int( lp, i + 1, TRUE ) ;*/

	set_add_rowmode(lp, FALSE ) ;

	for ( k = 0 ; k < intervalCnt ; ++k )
	{
		colno[k] = efCnt + k + 1 ;
		row[k] = 1 ;	
	}
	for ( ; k < varCnt - efCnt ; ++k )
	{
		colno[k] = k + efCnt + 1 ;
		row[k] = SLACK_PENALTY ; 
	}
	set_obj_fnex( lp, varCnt - efCnt, row, colno ) ;

	set_minim( lp ) ;
	set_verbose(lp, CRITICAL);	
	if ( solve( lp ) != 0 )
	{
		free( colno ) ;
		free( row ) ;
		free( cid ) ;
		delete_lp( lp ) ;
		return -1 ;
	}
	//printf( "LP_sum=%lf\n", sum ) ;
	//write_LP(lp, stdout);
	/*printf( "|" ) ;
	for ( i = 0 ; i < ecnt ; ++i )
	{
		if ( eused[i] )
			printf( "%d ", i ) ;
	}
	printf( "\n" ) ;*/
	
	k = 0 ;

	ret = get_objective( lp ) ;
	/*if ( right )
		printf("#Objective value: %lf\n", ret ) ;
	else	
		printf("Objective value: %lf\n", ret ) ;
	if ( right )
		printf("#Objective value: " ) ;
	else	
		printf("Objective value: " ) ;*/

	/*get_variables(lp, row) ;
	  for ( i = 0 ; i < efCnt + intervalCnt ; ++i )
	  printf( "%lf ", row[i] ) ;
	printf( "\n" ) ;*/
	free( cid ) ;
	free( colno ) ;
	free( row ) ;
	delete_lp( lp ) ;
	//printf( "\n\n" ) ;
	return ret ;
}

/**
  Test whether there is a path from exon s to exon t.
*/
int TestExonPath( int s, int t, int ecnt, struct _solveRegionData &data )
{
	int *exonPath = data.exonPath ;
	bool *eused = data.eused ;
	bool *graph = data.graph ;
	int index = s * ecnt + t ;

	if ( exonPath[index] != -1 )
		return exonPath[index] ;
	if ( s == t )
		return exonPath[index] = 1 ;
	int i ;
	for ( i = s + 1 ; i <= t ; ++i )
	{
		if ( eused[i] && graph[s * ecnt + i] && TestExonPath( i, t, ecnt, data ) )
			return exonPath[index] = 1 ;	
	}
	return exonPath[index] = 0 ;
}

/**
Verify whether an enumeration satisfy the constraint from the mate-pair reads.
@parm: initialTest, remove the constraints that contradict with all the exons.
*/
bool VerifyConstraintPoints( int ecnt, int cpCnt, int mainStrand, struct _solveRegionData &data, bool initialTest = false  )
{
	int i, j, k ;
	struct _enumExon *regionExons = data.regionExons ;
	bool *eused = data.eused ;
	int *exonPath = data.exonPath ;
	struct _constraintPoint *constraintPoints = data.constraintPoints ;

	memset( exonPath, -1, sizeof( int ) * ecnt * ecnt ) ;

	if ( initialTest )
	{
		for ( i = 0 ; i < cpCnt ; ++i )
			constraintPoints[i].valid = true ;
		for ( i = 0 ; i < ecnt ; ++i )
			eused[i] = true ;

		for ( i = 0 ; i < cpCnt ; ++i )
			if ( constraintPoints[i].type == 0 ||
				constraintPoints[i].type == 3 ||
				constraintPoints[i].type == 4 )
			{
				if ( constraintPoints[i].support < 10 )
					constraintPoints[i].valid = false ;
			}
	}

	for ( i = 0 ; i < cpCnt ; ++i )
	{
		if ( !constraintPoints[i].valid )
			continue ;

		if ( constraintPoints[i].type == 0 )
		{
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( !eused[j] )
					continue ;
				
				if ( regionExons[j].start == constraintPoints[i].startPos && 
					regionExons[j].end == constraintPoints[i].endPos &&
					regionExons[j].strand == constraintPoints[i].strand )
					break ;
			}
			
			if ( j >= ecnt )
			{
				if ( !initialTest )
					return false ;
				else
					constraintPoints[i].valid = false ;	
			}
		}
		/*else if ( constraintPoints[i].type == 1 )
		{
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( !eused[j] || regionExons[j].strand != constraintPoints[i].strand || 
					regionExons[j].start != constraintPoints[i].startPos )
					continue ;
				for ( k = j ; k < ecnt ; ++k )
				{
					if ( !eused[k] || regionExons[k].strand != constraintPoints[i].strand || 
							regionExons[k].end != constraintPoints[i].endPos )
						continue ;
					if ( TestExonPath( j, k ) )
					{
						break ;	
					}
				}

				if ( k < ecnt )
					break ;
			}

			if ( j >= ecnt )
			{
				if ( !initialTest )
					return false ;
				else
					constraintPoints[i].valid = false ;	
			}
		}
		else if ( constraintPoints[i].type == 2 )
		{
			if ( constraintPoints[i].strand != -1 ) // TODO: because of noise
				continue ;
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( !eused[j] || regionExons[j].start > point[ constraintPoints[i].startPos ] + 1 || 
					regionExons[j].end < point[ constraintPoints[i].startPos + 1] )
					continue ;
				for ( k = j ; k < ecnt ; ++k )
				{
					if ( !eused[k] || regionExons[k].start > point[ constraintPoints[i].endPos ] + 1 ||
						regionExons[k].end < point[ constraintPoints[i].endPos + 1 ] )
						continue ;
					if ( regionExons[j].strand != regionExons[k].strand || 
					     ( constraintPoints[i].strand != -1 && regionExons[j].strand != constraintPoints[i].strand ) )
					     continue ;
					if ( TestExonPath( j, k ) )
						break ;
				}
				if ( k < ecnt )
					break ;
			}
			if ( j >= ecnt )
			{
				printf( "=.= %d %d %d\n" , constraintPoints[i].startPos, constraintPoints[i].endPos, constraintPoints[i].strand ) ;
				return false ;
			}
		}*/
		else if ( constraintPoints[i].type == 3 )
		{
			if ( constraintPoints[i].strand != -1 && constraintPoints[i].strand != mainStrand ) // TODO: because of noise
				continue ;
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( eused[j] && regionExons[j].start <= constraintPoints[i].startPos 
					&& regionExons[j].end >= constraintPoints[i].endPos 
					&& ( constraintPoints[i].strand == -1 || constraintPoints[i].strand == regionExons[j].strand ) )
					break ;
			}

			if ( j >= ecnt )
			{
				if ( !initialTest )
					return false ;
				else
					constraintPoints[i].valid = false ;	
			}
			/*{
				printf( "+.+ %d %d %d %d %d %d\n" , constraintPoints[i].startPos + OFFSET, constraintPoints[i].endPos + OFFSET, constraintPoints[i].strand, regionExons[1].start + OFFSET, regionExons[1].end + OFFSET,
						ecnt ) ;
				return false ;
			}*/
		}
		else if ( constraintPoints[i].type == 4 )
		{
			if ( constraintPoints[i].strand != -1 && constraintPoints[i].strand != mainStrand ) // TODO: because of noise
				continue ;
			for ( j = 0 ; j < ecnt ; ++j )
			{
				if ( !eused[j] )
					continue ;
				if ( regionExons[j].start > constraintPoints[i].startPos ||
					regionExons[j].end < constraintPoints[i].endPos )
					continue ;
				if ( constraintPoints[i].precise[0] && regionExons[j].start != constraintPoints[i].startPos )
					continue ;

				for ( k = j ; k < ecnt ; ++k )
				{
					if ( !eused[k] )
						continue ;
					if ( regionExons[k].start > constraintPoints[i].otherPos[0] ||
						regionExons[k].end < constraintPoints[i].otherPos[1] )
						continue ;
					if ( constraintPoints[i].precise[1] && regionExons[k].end != constraintPoints[i].otherPos[1] )
						continue ;
					//printf( "### %d %d %d %d: %d %d\n", i, j, k, TestExonPath(j , k), constraintPoints[i].startPos, constraintPoints[i].endPos ) ;
					if ( TestExonPath( j, k, ecnt, data ) )
						break ;
				}

				if ( k < ecnt )
					break ;
			}
			if ( j >= ecnt )
			{
				if ( !initialTest )
					return false ;
				else
					constraintPoints[i].valid = false ;	
			}
			/*{
				printf( "=.= %d %d %d\n" , constraintPoints[i].startPos, constraintPoints[i].endPos, constraintPoints[i].strand ) ;
				return false ;
			}*/
		}
	}

	return true ;
}

/**
  The (inner boundaries) distance between two regionExons/
*/
void FindExonDist( int s, int t, int ecnt, struct _solveRegionData &data )
{
	struct _enumExon *regionExons = data.regionExons ;
	int *exonDist = data.exonDist ;
	bool *graph = data.graph ;
	bool *eused = data.eused ;
	int index = s * ecnt + t ;

	if ( exonDist[ index * 2 + 0] != -1 )
		return ;
	if ( s == t || graph[s * ecnt + t] )
	{
		exonDist[index * 2 + 0] = 0 ;
		exonDist[index * 2 + 1] = 0 ;
		return ;
	}
	int i, k ;
	for ( i = s + 1 ; i < t ; ++i )
	{
		int index2 = i * ecnt + t ;
		if ( !eused[i] || !graph[s * ecnt + i] || !TestExonPath( i, t, ecnt, data ) )
			continue ;
		FindExonDist( i, t, ecnt, data ) ;
		k = exonDist[index2 * 2 + 0] ;
		if ( exonDist[ index * 2 + 0] == -1 || 
			exonDist[ index2 * 2 + 0] + regionExons[i].end - regionExons[i].start + 1 < exonDist[index * 2 + 0] )
			exonDist[ index * 2 + 0] = exonDist[index2 * 2 + 0] + regionExons[i].end - regionExons[i].start + 1 ;

		k = exonDist[index2 * 2 + 1] ;
		if ( exonDist[index * 2 + 1] == -1 || 
			exonDist[ index2 * 2 + 1] + regionExons[i].end - regionExons[i].start + 1 > exonDist[index * 2 + 1] )
			exonDist[ index * 2 + 1] = exonDist[ index2 * 2 + 1] + regionExons[i].end - regionExons[i].start + 1 ;
	}
}

/**
  Find the distance between two intervals from the type-2 constraint points.
  return: The weight of pair intervals whose distance does not include the normal range. 
*/
double IntervalDistance( int readCnt, int ecnt, int cpCnt, struct _solveRegionData &data )
{
	int i, j, k ;
	double ret = 0 ;
	int min, max ;
	int rmin, rmax ;
	int scale = 0 ;

	struct _enumExon *regionExons = data.regionExons ;
	struct _constraintPoint *constraintPoints = data.constraintPoints ;
	int *exonDist = data.exonDist ;
	int *point = data.point ;
	bool *eused = data.eused ;

	memset( exonDist, -1, sizeof( exonDist[0] ) * ecnt * ecnt * 2  ) ;
	for ( i = 0 ; i < cpCnt ; ++i )
	{
		if ( constraintPoints[i].type != 2 )
			continue ;
		++scale ;
		rmin = MAX_LENGTH ;
		rmax = -1 ;
		for ( j = 0 ; j < ecnt ; ++j )
		{
			if ( !eused[j] || regionExons[j].start > point[ constraintPoints[i].startPos ] + 1 || 
					regionExons[j].end < point[ constraintPoints[i].startPos + 1] )
				continue ;
			for ( k = j ; k < ecnt ; ++k )
			{
				if ( !eused[k] || regionExons[k].start > point[ constraintPoints[i].endPos ] + 1 ||
						regionExons[k].end < point[ constraintPoints[i].endPos + 1 ] )
					continue ;
				if ( regionExons[j].strand != regionExons[k].strand || 
						( constraintPoints[i].strand != -1 && regionExons[j].strand != constraintPoints[i].strand ) )
					continue ;
				if ( k != j )
				{
					FindExonDist( j, k, ecnt, data ) ;
					int index = j * ecnt + k ;
					min = exonDist[2 * index + 0] + regionExons[j].end - point[ constraintPoints[i].startPos + 1 ] + 1 + 
						point[ constraintPoints[i].endPos ] - regionExons[k].start + 2 ;
					max = exonDist[2 * index + 1] + regionExons[j].end - point[ constraintPoints[i].startPos ] + 
						point[ constraintPoints[i].endPos + 1 ] - regionExons[k].start + 1 ;
				}
				else
				{
					min = point[ constraintPoints[i].endPos ] - point[ constraintPoints[i].startPos + 1 ] + 2 ;
					max = point[ constraintPoints[i].endPos + 1 ] - point[ constraintPoints[i].startPos ] ; 
				}

				if ( min < rmin )
					rmin = min ;
				if ( max > rmax )
					rmax = max ;
			}
		}
		if ( rmax < FRAG_LENGTH - 2 * FRAG_STD || rmin > FRAG_LENGTH + 2 * FRAG_STD ) // TODO: add parameter here.
		{
			//printf( "=.=: %d %d %d %d\n", constraintPoints[i].startPos, constraintPoints[i].endPos, rmin, rmax ) ;
			ret += constraintPoints[i].support / (double)readCnt ;	
		}		
	} // end for i
	return ret * scale ;
}

void EnumerateExon( int depth, int readCnt, int ecnt, int intervalCnt, int startCnt, int endCnt, int cpCnt, int mainStrand, 
	double &bestScore, bool *bestEused,
	struct _solveRegionData &data, bool *fixedEused = NULL )
{
	double lpResult, score ;
	struct _enumExon *regionExons = data.regionExons ;
	bool *eused = data.eused ;
	int *startCover = data.startCover ;
	int *endCover = data.endCover ;
	struct _point *start = data.start, *end = data.end ;

	if ( depth == ecnt )
	{
		int i, j, k ;
	
		// Validate this combination
		for ( i = 0 ; i < startCnt ; ++i )
			if ( !startCover[i * 3 + 0] || !startCover[i * 3 + 1] || !startCover[i * 3 + 2] )
				return ;
				
		for ( i = 0 ; i < endCnt ; ++i )
			if ( !endCover[i * 3 + 0] || !endCover[i * 3 + 1] || !endCover[i * 3 + 2] )
				return ;
		// Every interval should be covered by some regionExons.
		/*for ( i = 0 ; i < intervalCnt ; ++i )
		{
			for ( j = 0 ; j < ecnt ; ++j )
				if ( eused[j] && regionExons[j].start <= point[i] + 1 && regionExons[j].end >= point[i + 1] )
					break ;
			if ( j >= ecnt )
				return ;
		}*/
		/*if ( OFFSET == 156906608 )
		{
			for ( i = 0 ; i < ecnt ; ++i )
				if ( eused[i] )
					printf( "%d ", i ) ;
			printf( "\n" ) ;
			//exit( 1 ) ;
		}*/

		if ( !VerifyConstraintPoints( ecnt, cpCnt, mainStrand, data ) )
			return ;

		lpResult = LP( ecnt, intervalCnt, startCnt, endCnt, data ) ;
		if ( lpResult == -1 )
			return ;

		//exit( 0 ) ;
		int eusedCnt = 0 ;
		int &bestEcnt = data.bestCombAttr.bestEcnt ;
		int &bestEvidenceExonCnt = data.bestCombAttr.bestEvidenceExonCnt ;
		int &bestLongExonCnt = data.bestCombAttr.bestLongExonCnt ;
		int &bestBalancedExonCnt = data.bestCombAttr.bestBalancedExonCnt ; 
		
		int evidenceExonCnt = 0 ;
		int longExonCnt = 0 ;
		int balancedExonCnt = 0 ;
		for ( i = 0 ; i < ecnt ; ++i )
		{
			if ( eused[i] )
			{
				int l ;
				++eusedCnt ;
				//if ( regionExons[i].end - regionExons[i].start >= LONG_EXON )
				//	++longExonCnt ;

				if ( regionExons[i].inEvidence )
					++evidenceExonCnt ;

				for ( j = 0 ; j < startCnt ; ++j )
					if ( start[j].pos == regionExons[i].start && IsSameStrand( start[j].strand, regionExons[i].strand ) )
						break ;
				for ( k = 0 ; k < endCnt ; ++k )
					if ( end[k].pos == regionExons[i].end && IsSameStrand( end[k].strand, regionExons[k].strand ) )
						break ;
				
				if ( start[j].support <= end[k].support * 10 && end[k].support <= start[j].support * 10 )
					++balancedExonCnt ;

				for ( l = j + 1 ; l < startCnt ; ++l )
				{
					if ( start[l].pos < regionExons[i].end && IsSameStrand( start[l].strand, regionExons[i].strand ) )
						break ;
				}

				if ( l < startCnt )
				{
					++longExonCnt ;
					continue ;
				}

				for ( l = k - 1 ; l >= 0 ; --l )
				{
					if ( end[l].pos > regionExons[i].start && IsSameStrand( end[l].strand, regionExons[i].strand ) )
						break ;
				}
				if ( l >= 0 )
				{
					++longExonCnt ;
					continue ;
				}
			}
		}
		score = lpResult + IntervalDistance( readCnt, ecnt, cpCnt, data ) ; 

		/*for ( i = 0 ; i < ecnt ; ++i )
			printf( "%d ", eused[i] ) ;
		printf( ":" ) ;
		printf( "%lf %lf\n", lpResult, score ) ;*/		
		/*if ( OFFSET == 156914811 )
		{
			for ( i = 0 ; i < ecnt ; ++i )
				if ( eused[i] )
					printf( "%d ", i ) ;
			printf( ": %lf %d %d %d %d\n", lpResult + IntervalDistance( readCnt ), eusedCnt, evidenceExonCnt, balancedExonCnt, longExonCnt ) ;
			//exit( 1 ) ;
		}*/
		//printf( "score=%lf(%lf) bestScore=%lf\n", score, lpResult, bestScore) ;
		if ( ( bestScore > 0.01 * intervalCnt && score < bestScore ) || 
			( score <= 0.01 * intervalCnt && 
				( eusedCnt < bestEcnt ||
				  ( eusedCnt == bestEcnt && evidenceExonCnt < bestEvidenceExonCnt ) ||
				  ( eusedCnt == bestEcnt && evidenceExonCnt == bestEvidenceExonCnt && balancedExonCnt > bestBalancedExonCnt ) ||
				  ( eusedCnt == bestEcnt && evidenceExonCnt == bestEvidenceExonCnt && 
				  	balancedExonCnt == bestBalancedExonCnt && longExonCnt < bestLongExonCnt ) ) ) )
		{
			k = 0 ;
			bestScore = score ;
			for ( i = 0 ; i < ecnt ; ++i )
				bestEused[i] = eused[i] ;
			bestEcnt = eusedCnt ;
			bestEvidenceExonCnt = evidenceExonCnt ;
			bestBalancedExonCnt = balancedExonCnt ;
			bestLongExonCnt = longExonCnt ;
		}
		else if ( score <= bestScore + 1e-6 )
		{
			// both of them are smaller than 0.01 * intervalCnt
			if ( eusedCnt < bestEcnt ||
				( eusedCnt == bestEcnt && evidenceExonCnt < bestEvidenceExonCnt ) ||
				  ( eusedCnt == bestEcnt && evidenceExonCnt == bestEvidenceExonCnt && balancedExonCnt > bestBalancedExonCnt ) ||
				  ( eusedCnt == bestEcnt && evidenceExonCnt == bestEvidenceExonCnt && 
				  	balancedExonCnt == bestBalancedExonCnt && longExonCnt < bestLongExonCnt ) ) 
			{
				bestScore = score ;
				for ( i = 0 ; i < ecnt ; ++i )
					bestEused[i] = eused[i] ;
				bestEcnt = eusedCnt ;
				bestEvidenceExonCnt = evidenceExonCnt ;
				bestBalancedExonCnt = balancedExonCnt ;
				bestLongExonCnt = longExonCnt ;
			}
		}

		return ;
	}

	if ( fixedEused != NULL && fixedEused[depth] )
	{
		EnumerateExon( depth + 1, readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, mainStrand, 
			bestScore, bestEused, data, fixedEused ) ;
		return ;
	}

	eused[depth] = true ;
	int startStrand, endStrand ;
	startStrand = start[ regionExons[depth].startInd ].strand ;
	endStrand = end[ regionExons[depth].endInd ].strand ;
	if ( startStrand >= -1 ) 
		++startCover[ regionExons[depth].startInd * 3 + startStrand + 1] ;
	if ( endStrand >= -1 )
		++endCover[ regionExons[depth].endInd * 3 + endStrand + 1] ;
	EnumerateExon( depth + 1, readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, mainStrand, 
		bestScore, bestEused, data, fixedEused ) ;
	
	eused[depth] = false ;
	if ( startStrand >= -1 ) 
		--startCover[ regionExons[depth].startInd * 3 + startStrand + 1] ;
	if ( endStrand >= -1 )
		--endCover[ regionExons[depth].endInd * 3 + endStrand + 1] ;
	EnumerateExon( depth + 1, readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, mainStrand, 
		bestScore, bestEused, data, fixedEused ) ;
}

int Comp( const void *p1, const void *p2 )
{
	return *(int *)p1 - *(int *)p2 ;
}

int CompPoint( const void *p1, const void *p2 )
{
	return ( (struct _point *)p1 )->pos - ( (struct _point *)p2 )->pos ;
}
// Collect some statists from the interval
void IntervalStatistic( int intervalCnt, double *interval, int *point, double *ret, int *depth, int len )
{
	int i, j, k ;
	int sum ;
	double avg[MAX_POINT], dtmp ;
	ret[0] = ret[1] = -1 ;
	// Here ret[0] is the highest average coverage for an interval, ret[1] is the highest variance.
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		avg[i] = interval[i] / (double)( point[i + 1] - point[i] ) ;
		if ( avg[i] > ret[0] )
			ret[0] = avg[i] ;
	}

	k = 0 ;
	dtmp = 0 ;
	for ( i = 0 ; i < len ; ++i )
	{
		dtmp += ( depth[i] - avg[k] ) * ( depth[i] - avg[k] ) ; 

		if ( i == point[k + 1] )
		{
			dtmp /= (double)( point[k + 1] - point[k] ) ;
			if ( dtmp > ret[1] )
				ret[1] = dtmp ;
			dtmp = 0 ;
			++k ;
		}
	}
	
	// ret[2] is the minimal length of the intervals
	ret[2] = MAX_LENGTH ;
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		if ( point[i + 1] - point[i] < ret[2] )
			ret[2] = point[i + 1] - point[i] ;
	}
}

void BuildExonGraph( int ecnt, int siCnt, struct _solveRegionData &data )
{
	int i, j, k ;
	struct _enumExon *regionExons = data.regionExons ;
	struct _spliceIndex *spliceIndices = data.spliceIndices ;
	bool *graph = data.graph ;
	bool *beyond = data.beyond ;
	
	memset( graph, false, sizeof( graph[0] ) * ecnt * ecnt ) ;
	//memset( beyond, false, sizeof( beyond[0] ) * ecnt * ecnt * 4 ) ;
	for ( i = 0 ; i < ecnt - 1 ; ++i )
	{
		for ( j = i + 1 ; j < ecnt ; ++j )
		{
			if ( regionExons[i].strand != regionExons[j].strand )
				continue ;
			for ( k = 0 ; k < siCnt ; ++k )
			{
				if ( regionExons[i].strand == spliceIndices[k].strand && regionExons[i].endInd == spliceIndices[k].end 
					&& regionExons[j].startInd == spliceIndices[k].start )
				{
					graph[i * ecnt + j] = graph[j * ecnt + i] = true ;
				}
					
			}
		}
	}
	/*for ( i = 0 ; i < ecnt ; ++i )
	{
		for ( k = 0 ; k < siCnt ; ++k )
		{
			if ( regionExons[i].strand == spliceIndices[k].strand && regionExons[i].startInd == spliceIndices[k].start && spliceIndices[k].end == -1 )
				beyond[i * 4 + 2 * 0 + regionExons[i].strand ] = true ;

			if ( regionExons[i].strand == spliceIndices[k].strand && regionExons[i].endInd == spliceIndices[k].end && spliceIndices[k].start == -1 )
				beyond[i * 4 + 2 * 1 + regionExons[i].strand ] = true ;
		}
	}*/
}

/**
  Test whether the exons are in the evidences.
*/
void AreExonsInEvidences( int ecnt, struct _evidence *evidences, int eviCnt )
{
	if ( evidences == NULL )
		return ;

	static int eviTag = 0 ;
	static int prevEviCnt = 0 ;
	int i, j, k ;

	if ( prevEviCnt != eviCnt ) // New chromosome. TODO: Make this more robust
		eviTag = 0 ;
	prevEviCnt = eviCnt ;
	
	while ( eviTag < eviCnt && ( evidences[eviTag].exons[0].start > regionExons[ ecnt - 1 ].start + OFFSET + EVIDENCE_MARGIN 
		|| evidences[eviTag].exons[ evidences[eviTag].ecnt - 1 ].start < regionExons[0].start + OFFSET - EVIDENCE_MARGIN ) )
		++eviTag ;
	
	for ( i = 0 ; i < ecnt ; ++i )
	{
		for ( j = eviTag ; j < eviCnt ; ++j )
		{
			if ( regionExons[i].inEvidence == true || evidences[j].exons[0].start > regionExons[i].start + OFFSET + EVIDENCE_MARGIN )
				break ;

			if ( evidences[j].exons[0].strand != regionExons[i].strand )
				continue ;

			for ( k = 0 ; k < evidences[j].ecnt ; ++k )
				if ( WithinEvidenceMargin( evidences[j].exons[k].start, regionExons[i].start + OFFSET ) 
					&& WithinEvidenceMargin( evidences[j].exons[k].end, regionExons[i].end + OFFSET ) )
				{
					regionExons[i].inEvidence = true ;		
					break ;
				}
				else if ( regionExons[i].start + OFFSET < evidences[j].exons[k].start - EVIDENCE_MARGIN )
					break ;
		}
	}
}

/**
  Merge the points based on strand(-:0, +:1)
  cnt: |points|
*/
void MergePoints( struct _point points[], int cnt, int strand, int factor )
{
	int i, j ;
	for ( j = 0 ; j < cnt && ( points[j].strand != strand || points[j].merge != j ) ; ++j )
		;

	for ( i = j + 1 ; i < cnt ; ++i )
	{
		if ( points[i].strand != strand || points[i].merge != i )//&& points[i].strand != -1 )
			continue ;

		if ( points[i].pos - points[j].pos > factor * MERGE_DISTANCE ||
			( points[i].support > 50 && points[j].support > 50 ) ) 
		{
			j = i ;
			continue ;
		}
		
		if ( points[i].support >= points[j].support )
		{
			// merge j to i.
			points[j].merge = i ;
			
			j = i ;	
			// TODO: Do we need to put j's support to i?
			points[i].support += points[j].support ;
		}
		else
		{
			// merge i to j.
			points[i].merge = j ;
			
			points[j].support += points[i].support ;
		}
	}
}

/**
  Update the splice indices according to the merge information.
*/
void UpdateSpliceIndices( struct _spliceIndex spliceIndices[], int siCnt, struct _point start[], struct _point end[], int startCnt, int endCnt )
{
	int newStartInd[MAX_POINT], newEndInd[MAX_POINT] ;
	int i, k ;	
	
	// Get the indices after merging. These are the indices used in the finalized "start", "end" array.
	k = 0 ;
	for ( i = 0 ; i < startCnt ; ++i )
	{
		if ( start[i].merge != i )
			continue ;
		newStartInd[i] = k ;
		++k ;
	}
	
	k = 0 ;
	for ( i = 0 ; i < endCnt ; ++i )
	{
		if ( end[i].merge != i )
			continue ;
		newEndInd[i] = k ;
		++k ;
	}

	// Update the splice indices
	for ( i = 0 ; i < siCnt ; ++i )
	{
		if ( spliceIndices[i].start != -1 )
		{
			spliceIndices[i].start = newStartInd[ GetMergeFather( start, spliceIndices[i].start ) ] ;
			//spliceIndices[i].startPos = start[ spliceIndices[i].start ].pos ;
		}
		if ( spliceIndices[i].end != -1 )
		{
			spliceIndices[i].end = newEndInd[ GetMergeFather( end, spliceIndices[i].end ) ] ;
			//spliceIndices[i].endPos = end[ spliceIndices[i].start ].pos ;
		}
	}
}

/**
  Collect the constraints from the reads
*/
void GetConstraintPoints( struct _point start[], int startCnt, struct _point end[], int endCnt,
	int point[], int intervalCnt,
	struct _read reads[], int readCnt, 
	struct _constraintPoint constraintPoints[], int &cpCnt )
{
	int i, j, k, mate, l, tag ;
	bool flag ;
	cpCnt = 0 ;
	//printf( "## %d\n", readCnt ) ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		if ( !reads[i].scnt )
			continue ;
		for ( j = 0 ; j + 1 < reads[i].scnt ; ++j )
		{
			// make sure the splice site is in the points after merging.
			for ( l = 0 ; l < startCnt ; ++l ) // TODO: Use binary search
				if ( start[l].strand == reads[i].strand && 
					start[l].pos == reads[i].splices[j][1] )
					break ;
			if ( l >= startCnt )
				continue ;

			for ( l = 0 ; l < endCnt ; ++l ) // TODO: Use binary search
				if ( end[l].strand == reads[i].strand &&
					end[l].pos == reads[i].splices[j + 1][0] )
					break ;
			if ( l >= endCnt )
				continue ;

			//fprintf( stderr, "## hi %d %d %d\n", start[0].pos + OFFSET, end[ endCnt - 1 ].pos + OFFSET, reads[i].start  ) ;
			
			// Removing repeated constraint points.
			flag = false ;
			for ( l = cpCnt - 1 ; l >= 0 && constraintPoints[l].startPos >= reads[i].start - OFFSET ; --l )
				if ( constraintPoints[l].startPos == reads[i].splices[j][1] - OFFSET &&
					constraintPoints[l].endPos == reads[i].splices[j + 1][0] - OFFSET &&
					constraintPoints[l].strand == reads[i].strand && constraintPoints[l].type == 0 )
				{
					flag = true ;
					constraintPoints[l].type = 0 ;
					++constraintPoints[l].support ;
					break ;
				}

			if ( flag )
				continue ;

			constraintPoints[ cpCnt ].startPos = reads[i].splices[j][1] - OFFSET ;
			constraintPoints[ cpCnt ].endPos = reads[i].splices[j + 1][0] - OFFSET;
			constraintPoints[ cpCnt ].strand = reads[i].strand ;			
			constraintPoints[ cpCnt ].type = 0 ;
			constraintPoints[ cpCnt ].support = 1 ;
			++cpCnt ;
		}

		mate = reads[i].mateInd ;
		if ( mate <= i || !reads[mate].scnt || reads[i].splices[ reads[i].scnt - 1 ][1] >= reads[mate].splices[0][0])
			continue ;
			
		//fprintf( stdout, "## hi %d %d %d %d %d\n", start[0].pos + OFFSET, end[ endCnt - 1 ].pos + OFFSET, reads[i].splices[0][0], reads[i].splices[0][1], reads[i].scnt  ) ;
		for ( l = 0 ; l < startCnt ; ++l ) // TODO: Use binary search
			if ( start[l].strand == reads[i].strand && 
					start[l].pos == reads[i].splices[ reads[i].scnt - 1 ][1] - OFFSET )
				break ;
		if ( l >= startCnt )
			continue ;

		for ( l = 0 ; l < endCnt ; ++l ) // TODO: Use binary search
			if ( end[l].strand == reads[i].strand && 
				end[l].pos == reads[mate].splices[0][0] - OFFSET )
				break ;
		if ( l >= endCnt )
			continue ;

		// We only need to add the last splice point of the first read and the first splice point of the second read.
		// Other constraint can be deduced from the type-0 constraint points.
		flag = false ;
		for ( l = cpCnt - 1 ; l >= 0 && constraintPoints[l].startPos >= reads[i].start - OFFSET ; --l )
			if ( constraintPoints[l].startPos == reads[i].splices[ reads[i].scnt - 1 ][1] - OFFSET &&
					constraintPoints[l].endPos == reads[ mate ].splices[0][0] - OFFSET &&
					constraintPoints[l].strand == reads[i].strand )
			{
				flag = true ;
				//constraintPoints[l].type = 0 ;
				break ;
			}
		if ( flag )
			continue ;
		
		constraintPoints[ cpCnt ].startPos = reads[i].splices[j][1] - OFFSET ;
		constraintPoints[ cpCnt ].endPos = reads[mate].splices[0][0] - OFFSET;
		constraintPoints[ cpCnt ].strand = reads[i].strand ;			
		constraintPoints[ cpCnt ].type = 1 ;
		++cpCnt ;
	}	

	// Type-2 constraint. Used to measure the distance between two intervals if there are mate pair on it.
	tag = cpCnt ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		mate = reads[i].mateInd ;	
		if ( mate < i ) //|| reads[i].end >= reads[mate].start ) //|| reads[i].strand != -1 ) //|| reads[i].start - OFFSET < 0 || reads[i].end > point[ intervalCnt + 1 ] )
			continue ;
		int a, b ;
		a = reads[i].start - OFFSET ;
		b = reads[mate].end - OFFSET ;

		for ( j = 0 ; j < intervalCnt ; ++j )
		{
			if ( point[j] < a && a <= point[j + 1] )
				break ;
		}
		
		for ( l = j ; l < intervalCnt ; ++l )
		{
			if ( point[l] < b && b <= point[l + 1] )
				break ;
		}

		if ( j >= intervalCnt || l >= intervalCnt )
			continue ;

		for ( k = tag ; k < cpCnt ; ++k )
			if ( constraintPoints[k].strand == reads[i].strand && constraintPoints[k].startPos == j && constraintPoints[k].endPos == l )
			{
				++constraintPoints[k].support ;
				break ;
			}
		if ( k >= cpCnt )
		{
			constraintPoints[ cpCnt ].startPos = j ;
			constraintPoints[ cpCnt ].endPos = l ;
			constraintPoints[ cpCnt ].strand = reads[i].strand ;
			constraintPoints[ cpCnt ].type = 2 ;
			constraintPoints[ cpCnt ].support = 1 ;	
			++cpCnt ;
		}
	}
	tag = cpCnt ;
	// type-3
	for ( i = 0 ; i < readCnt ; ++i )
	{
		int rstart, rend, s ;
		int l ;
		for ( j = 0 ; j < reads[i].scnt + 1 ; ++j )
		{
			if ( j == 0 )
				rstart = reads[i].start ;
			else
				rstart = reads[i].splices[j - 1][1] ;

			if ( j == reads[i].scnt )
				rend = reads[i].end ;
			else
				rend = reads[i].splices[j][0] ;

			s = -1 ;
			// Ignore the part before and after the region.
			for ( k = 1 ; k < intervalCnt ; ++k )
			{
				if ( rstart - OFFSET <= point[k] && rend - OFFSET > point[k] )
				{
					if ( s == -1 )
						s = k ;
				}
				else if ( rstart - OFFSET < point[k] )
					break ;
			}

			if ( s != -1 )
			{
				l = k - 1 ;
				for ( k = tag ; k < cpCnt ; ++k  )
				{
					if ( constraintPoints[k].startPos <= point[s - 1] + 1 &&
						constraintPoints[k].endPos >= point[l + 1] )
					{
						++constraintPoints[k].support ;
						break ;
					}
					else if ( constraintPoints[k].startPos > point[s - 1] + 1 &&
						constraintPoints[k].endPos < point[l + 1] )
					{
						constraintPoints[k].startPos = point[s - 1] + 1 ;
						constraintPoints[k].endPos = point[l + 1] ;
						++constraintPoints[k].support ;
						break ;
					}

					/*if ( constraintPoints[k].startPos == point[s - 1] + 1 &&
						constraintPoints[k].endPos == point[l + 1] )
					{
						++constraintPoints[k].support ;
						break ;
					}*/
				}

				if ( k >= cpCnt )
				{
					constraintPoints[ cpCnt ].type = 3 ;
					constraintPoints[ cpCnt ].strand = reads[i].strand ;
					constraintPoints[ cpCnt ].startPos = point[s - 1] + 1 ;
					constraintPoints[ cpCnt ].endPos = point[l + 1] ;
					constraintPoints[ cpCnt ].support = 1 ;
					++cpCnt ;
				}
			}
		}
	}
	// Type4: The reachability between two mate pairs.
	tag = cpCnt ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		mate = reads[i].mateInd ;	
		if ( mate < i ) //|| reads[i].end >= reads[mate].start ) //|| reads[i].strand != -1 ) //|| reads[i].start - OFFSET < 0 || reads[i].end > point[ intervalCnt + 1 ] )
			continue ;
		int a, b, c, d ;
		int tmp ;
		int precise[2] = { 0, 0 };
		int m, n ;
		if ( reads[i].scnt > 0 )
		{
			a = reads[i].splices[ reads[i].scnt - 1 ][1] ;
			if ( SplicesInformation_Get( a, reads[i].strand ) != -1 )
				continue ;
			precise[0] = 1 ;
			a -= OFFSET ;
		}
		else
			a = reads[i].start - OFFSET ;
		b = reads[i].end - OFFSET ;
		c = reads[mate].start - OFFSET ;
		if ( reads[mate].scnt > 0 )
		{
			d = reads[mate].splices[0][0] ;
			if ( SplicesInformation_Get( d, reads[i].strand ) != -1 )
				continue ;
			precise[1] = 1 ;
			d -= OFFSET ;
		}
		else
			d = reads[mate].end - OFFSET ;

		if ( b < a )
		{
			a = tmp ; a = b ; b = tmp ;
		}

		for ( j = 0 ; j < intervalCnt ; ++j )
		{
			if ( point[j] < a && a <= point[j + 1] )
				break ;
		}
		
		for ( l = j ; l < intervalCnt ; ++l )
		{
			if ( point[l] < b && b <= point[l + 1] )
				break ;
		}

		for ( m = l ; m < intervalCnt ; ++m )
		{
			if ( point[m] < c && c <= point[m + 1] )
				break ;
		}

		for ( n = m ; n < intervalCnt ; ++n )
		{
			if ( point[n] < d && d <= point[n + 1] )
				break ;
		}
		
		if ( j >= intervalCnt || l >= intervalCnt || m >= intervalCnt || n >= intervalCnt )
			continue ;

		for ( k = tag ; k < cpCnt ; ++k )
			if ( constraintPoints[k].strand == reads[i].strand && constraintPoints[k].startPos == point[j] + 1 && constraintPoints[k].endPos == point[l + 1] &&
				constraintPoints[k].otherPos[0] == point[m] + 1 && constraintPoints[k].otherPos[1] == point[n + 1] &&
				constraintPoints[k].precise[0] == precise[0] && constraintPoints[k].precise[1] == precise[1] )
			{
				++constraintPoints[k].support ;
				break ;
			}
		if ( k >= cpCnt )
		{
			constraintPoints[ cpCnt ].startPos = point[j] + 1 ;
			constraintPoints[ cpCnt ].endPos = point[l + 1] ;
			constraintPoints[ cpCnt ].otherPos[0] = point[m] + 1 ;
			constraintPoints[ cpCnt ].otherPos[1] = point[n + 1] ;
			constraintPoints[ cpCnt ].precise[0] = precise[0] ;
			constraintPoints[ cpCnt ].precise[1] = precise[1] ;
			constraintPoints[ cpCnt ].strand = reads[i].strand ;
			constraintPoints[ cpCnt ].type = 4 ;
			constraintPoints[ cpCnt ].support = 1 ;
			++cpCnt ;
		}
	}
	//for ( i = 0 ; i < cpCnt ; ++i )
	//	printf( "%d %d %d\n", constraintPoints[i].type, constraintPoints[i].startPos + OFFSET, constraintPoints[i].endPos + OFFSET ) ;
}

/**
	Directly build the exons by pairing the nearest splice sites.
	@return: the Objective value for this combination.
*/
double DirectBuildExons( int readCnt, int ecnt, int intervalCnt, int startCnt, int endCnt, int cpCnt, struct _solveRegionData &data )
{
	int i, j, k ;
	int startStrand, endStrand ;
	struct _enumExon *regionExons = data.regionExons ;
	bool *eused = data.eused ;
	int *startCover = data.startCover ;
	int *endCover = data.endCover ;
	struct _point *start = data.start, *end = data.end ;
	memset( eused, false, sizeof( eused[0] ) * ecnt) ;
	memset( startCover, 0, sizeof( startCover[0] ) * ecnt * 3 ) ;
	memset( endCover, 0, sizeof( endCover[0] ) * ecnt * 3 ) ;
//	printf( "%d\n", OFFSET ) ;

	for ( i = 0 ; i < ecnt ; ++i )
	{
		if ( regionExons[i].inEvidence )
		{
			eused[i] = true ;
			startStrand = start[ regionExons[i].startInd ].strand ;
			endStrand = end[ regionExons[i].endInd ].strand ;

			if ( startStrand >= -1 )
				startCover[ regionExons[i].startInd * 3 + startStrand + 1 ] = 1 ;
			if ( endStrand >= -1 )
				endCover[ regionExons[i].endInd * 3 + endStrand + 1 ] = 1 ;
		}
	}
	for ( i = 0 ; i < startCnt ; ++i )
	{
		if ( start[i].strand < -1 ||  startCover[i * 3 + start[i].strand + 1] )
			continue ;
		// First, try to pair with the splice sites of similar support
		for ( j = 0 ; j < endCnt ; ++j )
			if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) &&
				end[j].pos - start[i].pos + 1 < LONG_EXON && 
				end[j].support <= start[i].support * 10 && start[i].support <= end[j].support * 10 )
			{
				for ( k = 0 ; k < ecnt ; ++k )
				{
					if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
					{
						startCover[i * 3 + start[i].strand + 1 ] = 1 ;
						endCover[j * 3 + end[j].strand + 1] = 1 ;
				//printf( "Short: %d %d %d\n", i, j, k ) ;
						eused[k] = true ;
						break ;
					}
				}
				break ;
			}
		// If did not find a suitable pairing site
		if ( startCover[i * 3 + start[i].strand + 1 ] == 0 )
		{
			for ( j = 0 ; j < endCnt ; ++j )
				if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
				{
					for ( k = 0 ; k < ecnt ; ++k )
					{
						if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
						{
							startCover[i * 3 + start[i].strand + 1 ] = 1 ;
							endCover[j * 3 + end[j].strand + 1] = 1 ;
							//printf( "Short: %d %d %d\n", i, j, k ) ;
							eused[k] = true ;
							break ;
						}
					}
					break ;
				}
		}
	}

	for ( j = 0 ; j < endCnt ; ++j )
	{
		if ( end[j].strand < -1 ||  endCover[j * 3 + end[j].strand + 1] )
			continue ;
		
		// First, try to pair with the splice sites of similar support
		for ( i = startCnt - 1 ; i >= 0 ; --i )
			if ( end[j].pos > start[i].pos && IsSameStrand( start[i].strand, end[j].strand ) &&
				end[j].pos - start[i].pos + 1 < LONG_EXON && 
				end[j].support <= start[i].support * 10 && start[i].support <= end[j].support * 10 )
			{
				for ( k = 0 ; k < ecnt ; ++k )
				{
					if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
					{
						//printf( "Short: %d %d\n", i, j ) ;
						startCover[i * 3 + start[i].strand + 1 ] = 1 ;
						endCover[j * 3 + end[j].strand + 1] = 1 ;
						eused[k] = true ;
						break ;
					}
				}
				break ;
			}

		// If did not find a suitable pairing site
		if ( endCover[j * 3 + end[j].strand + 1] == 0 )
		{
			for ( i = startCnt - 1 ; i >= 0 ; --i )
				if ( end[j].pos > start[i].pos && IsSameStrand( start[i].strand, end[j].strand ) )
				{
					for ( k = 0 ; k < ecnt ; ++k )
					{
						if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
						{
							//printf( "Short: %d %d\n", i, j ) ;
							startCover[i * 3 + start[i].strand + 1 ] = 1 ;
							endCover[j * 3 + end[j].strand + 1] = 1 ;
							eused[k] = true ;
							break ;
						}
					}
					break ;
				}
		}
	}

	// Two long exons to cover the region
	for ( j = endCnt - 1 ; j >= 0 ; --j )
	{
		i = 0 ;
		if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
		{
			for ( k = 0 ; k < ecnt ; ++k )
			{
				if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
				{
					//printf( "Long: %d %d %d\n", i, j, k ) ;
					eused[k] = true ;
					break ;
				}
			}

			break ;
		}
	}

	for ( i = 0 ; i < startCnt ; ++i )
	{
		j = endCnt - 1 ;
		if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
		{
			for ( k = 0 ; k < ecnt ; ++k )
			{
				if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
				{
					//printf( "Long: %d %d %d\n", i, j, k ) ;
					eused[k] = true ;
					break ;
				}
			}

			break ;
		}
	}
	//return 0 ;
	return LP( ecnt, intervalCnt, startCnt, endCnt, data, true )  + IntervalDistance( readCnt, ecnt, cpCnt, data ) ;
}


void AddAllExons( int OFFSET, int ecnt, int siCnt, bool *bestEused, 
	struct _enumExon *regionExons, struct _spliceIndex *spliceIndices,
	struct _exon *allExons, int &exonCnt ) 
{
	int i, j, k ;
	for ( i = 0 ; i < ecnt ; ++i )
	{
		if ( bestEused[i] )
		{
			allExons[ exonCnt ].start = regionExons[i].start + OFFSET ;
			allExons[ exonCnt ].end = regionExons[i].end + OFFSET ;
			allExons[ exonCnt ].strand = regionExons[i].strand ;
			allExons[ exonCnt ].pcnt = allExons[ exonCnt ].ncnt = 0 ;
			/*if ( end[ regionExons[i].endInd ].strand < -1 )//&& end[ regionExons[i].endInd ].strong == 0 ) 
			{
				printf( "##: polya %d %d\n", regionExons[i].end + OFFSET, end[ regionExons[i].endInd ].strand + 3  ) ;	
			}
			if ( start[regionExons[i].startInd].strand < -1 )//&& start[ regionExons[i].startInd ].strong == 0 )
			{
				printf( "##: polya %d %d\n", regionExons[i].start + OFFSET, start[regionExons[i].startInd].strand + 3  ) ;	
			}*/

			/*if ( regionExons[i].start + OFFSET == 4461426 )
			  {
			  for ( j = 0 ; j < siCnt ; ++j )
			  {
			  printf( "# %d %d %d %d\n", spliceIndices[j].endPos, spliceIndices[j].startPos, spliceIndices[j].strand, regionExons[i].strand ) ;
			  }
			  exit( 1 ) ;
			  }*/	

			for ( j = 0 ; j < siCnt ; ++j )
			{
				if ( IsSameStrand( spliceIndices[j].strand, regionExons[i].strand ) &&
						spliceIndices[j].start == regionExons[i].startInd )
				{
					for ( k = 0 ; k < allExons[ exonCnt ].pcnt ; ++k )
						if ( allExons[ exonCnt ].prev[k] == spliceIndices[j].endPos )
							break ;
					if ( k < allExons[ exonCnt ].pcnt )
						continue ;
					if ( k >= allExons[ exonCnt ].psize ) // No need to use while here
					{
						if ( allExons[ exonCnt ].psize == -1 )
						{
							allExons[ exonCnt ].psize = 40 ;
							allExons[ exonCnt ].prev = (int *)malloc( sizeof( int ) * allExons[ exonCnt ].psize ) ;
						}		
						else
						{
							memcpy( buffer2, allExons[ exonCnt ].prev, sizeof( int ) * allExons[ exonCnt ].psize ) ;
							allExons[ exonCnt ].psize *= 2 ;
							free( allExons[ exonCnt ].prev ) ;

							allExons[ exonCnt ].prev = (int *)malloc( sizeof( int ) * allExons[ exonCnt ].psize ) ;
							memcpy( allExons[ exonCnt ].prev, buffer2, sizeof( int ) * allExons[ exonCnt ].psize / 2 ) ;
						}
					}

					allExons[ exonCnt ].prev[k] = spliceIndices[j].endPos ;
					++allExons[exonCnt].pcnt ;
				}	
			}	

			for ( j = 0 ; j < siCnt ; ++j )
			{
				//printf( "== %d %d %d %d %d\n", j, spliceIndices[j].endPos, spliceIndices[j].startPos,
				//		spliceIndices[j].end, regionExons[i].endInd ) ;
				if ( IsSameStrand( spliceIndices[j].strand, regionExons[i].strand ) &&
						spliceIndices[j].end == regionExons[i].endInd )
				{
					for ( k = 0 ; k < allExons[ exonCnt ].ncnt ; ++k )
						if ( allExons[ exonCnt ].next[k] == spliceIndices[j].startPos )
							break ;
					if ( k < allExons[ exonCnt ].ncnt )
						continue ;
					if ( k >= allExons[ exonCnt ].nsize ) // No need to use while here
					{
						if ( allExons[ exonCnt ].nsize == -1 )
						{
							allExons[ exonCnt ].nsize = 40 ;
							allExons[ exonCnt ].next = (int *)malloc( sizeof( int ) * allExons[ exonCnt ].nsize ) ;
						}		
						else
						{
							memcpy( buffer2, allExons[ exonCnt ].next, sizeof( int ) * allExons[ exonCnt ].nsize ) ;
							allExons[ exonCnt ].nsize *= 2 ;
							free( allExons[ exonCnt ].next ) ;

							allExons[ exonCnt ].next = (int *)malloc( sizeof( int ) * allExons[ exonCnt ].nsize ) ;
							memcpy( allExons[ exonCnt ].next, buffer2, sizeof( int ) * allExons[ exonCnt ].nsize / 2 ) ;
							//printf( "############### %d\n", allExons[ exonCnt ].nsize ) ;
						}
					}
					allExons[ exonCnt ].next[k] = spliceIndices[j].startPos ;
					//printf( "===%d %d %d\n", exonCnt, allExons[ exonCnt ].strand, allExons[ exonCnt ].next[k] ) ;
					++allExons[exonCnt].ncnt ;
				}
			}
			//printf( "### %d %d %d\n", regionExons[ i ].start + OFFSET, regionExons[ i ].end + OFFSET, regionExons[i].strand ) ;
			//allExons[ exonCnt ].next[0], allExons[ exonCnt ].next[1] ) ;
			/*if ( allExons[ exonCnt ].pcnt >= MAX_NEIGHBOR || allExons[ exonCnt ].ncnt >= MAX_NEIGHBOR )
			  {
			  printf( "###ERROR### : %d %d (%d %d)\n", allExons[ exonCnt ].pcnt, allExons[ exonCnt ].ncnt,
			  allExons[exonCnt].start, allExons[exonCnt].end ) ;
			  exit( 1 ) ;
			  }*/

			/*printf( "%d: (%d %d):", exonCnt, allExons[ exonCnt ].start, allExons[ exonCnt ].end ) ;
			  for ( k = 0 ; k < allExons[ exonCnt ].ncnt ; ++k )
			  printf( "%d ", allExons[ exonCnt ].next[k] ) ;				
			  printf( "\n" ) ;*/
			//printf( "### %d \n", allExons[ exonCnt ].ncnt, allExons[ exonCnt ].next[0], allExons[ exonCnt ].next[1] ) ;


			++exonCnt ;
		}
	}
}

// The function used for pthread
void *SolveRegion_Thread( void *arg )
{
	int i, j, k ;
	struct _pthreadArgSolveRegion *pArg = ( struct _pthreadArgSolveRegion *)arg ;
	double bestScore ;

	// Read out the numbers from the argument bundle
	struct _enumExon *regionExons = pArg->regionExons ;
	int ecnt = pArg->ecnt ;
	bool *exonHasSoft = pArg->exonHasSoft ;
	struct _point *start = pArg->start ;
	int *startCover = pArg->startCover ;
	int startCnt = pArg->startCnt ;
	struct _point *end = pArg->end ;
	int endCnt = pArg->endCnt ;
	int *endCover = pArg->endCover ;
	struct _spliceIndex *spliceIndices = pArg->spliceIndices ;
	int siCnt = pArg->siCnt ;
	struct _constraintPoint *constraintPoints = pArg->constraintPoints ;
	int cpCnt = pArg->cpCnt ;
	struct _exon *allExons = pArg->allExons ;
	int &exonCnt = *( pArg->exonCnt ) ; 
	int OFFSET = pArg->OFFSET ;
	int mainStrand = pArg->mainStrand ;
	int *point = pArg->point ;
	double *interval = pArg->interval ;
	int intervalCnt = pArg->intervalCnt ;
	int readCnt = pArg->readCnt ;
	
	// Allocate  memories for data structures.
	bool *graph = ( bool *)malloc( sizeof( bool ) * ecnt * ecnt ) ;
	bool *eused = ( bool *)malloc( sizeof( bool ) * ecnt ) ;
	bool *bestEused = ( bool * )malloc( sizeof( bool ) * ecnt ) ;
	int *exonPath = ( int * )malloc( sizeof( int ) * ecnt * ecnt ) ;
	int *exonDist = ( int * )malloc( sizeof( int ) * ecnt * ecnt * 2 ) ;
	bool *visited = ( bool *)malloc( sizeof( bool ) * ecnt * 2 ) ;
	bool *beyond = ( bool *)malloc( sizeof( bool ) * ecnt * 4 ) ;
	struct _exonfrag *exonfrag = ( struct _exonfrag *)malloc( sizeof( struct _exonfrag ) * ecnt * intervalCnt) ;
	int *leftEf = ( int * )malloc( sizeof( int ) * ecnt ) ;
	int *rightEf = ( int * )malloc( sizeof( int ) * ecnt ) ;
	int *efid = ( int *)malloc( sizeof( int ) * ecnt * 2 ) ;	

	// Assigne the pointers in the "data"
	struct _solveRegionData data ;
	data.eused = eused ;
	data.graph = graph ;
	data.visited = visited ;
	data.beyond = beyond ;
	data.exonfrag = exonfrag ;
	data.efid = efid ;
	data.leftEf = leftEf ;
	data.rightEf = rightEf ;
	data.regionExons = regionExons ;
	data.exonHasSoft = exonHasSoft ;
	data.exonPath = exonPath ;
	data.exonDist = exonDist ;
	data.constraintPoints = constraintPoints ;
	data.start = start ;
	data.end = end ;
	data.startCover = startCover ;
	data.endCover = endCover ;
	data.spliceIndices = spliceIndices ;
	data.point = point ;
	data.interval = interval ;

	/*data.bestCombAttr.bestEcnt = MAX_LENGTH ;
	data.bestCombAttr.bestEvidenceExonCnt = -1 ;
	data.bestCombAttr.bestLongExonCnt = -1 ;
	data.bestCombAttr.bestBalancedExonCnt = -1 ;*/

	BuildExonGraph( ecnt, siCnt, data ) ;	
	VerifyConstraintPoints( ecnt, cpCnt, mainStrand, data, true ) ;
	
	memset( eused, false, sizeof( eused[0] ) * ecnt ) ;
	//memset( bestEused, false, sizeof( bestEused )) ;
	for ( i = 0 ; i < ecnt ; ++i )
		bestEused[i] = true ;
		
	//printf( "enter %d\n", ecnt ) ;	
	//printf( "############ %s %d %d %d %lf %lf %d\n", chrom, start[0].pos + OFFSET, 
	//		end[ endCnt - 1 ].pos + OFFSET, sum, stat[0], stat[1], (int)stat[2] ) ;
	//printf( "%d\n", ecnt ) ;
	if ( ecnt <= 2 )
	{
		for ( i = 0 ; i < ecnt ; ++i )
			bestEused[i] = true ;
	}
	else if ( ecnt < 16 )
	{
		bestScore = MAX_LENGTH ;
		EnumerateExon( 0, readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, mainStrand, 
			bestScore, bestEused, data ) ;
		//printf( "%d %lf\n", start[0].pos + OFFSET, bestScore ) ;

		if ( bestScore > 1 * intervalCnt )
		{
			DirectBuildExons( readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, data ) ;
			for ( i = 0 ; i < ecnt ; ++i )
				bestEused[i] = eused[i] ;
		}
	}
	else //if ( ecnt < 32 || force )
	{
		//return 0 ;
		//printf( "### %d %d\n", OFFSET, ecnt ) ;
		int mintag = -1 ;
		double prevObj, minObj, obj ;
		bool useExonFrag = ecnt < 50 ;
		int longExons[2] ;
		// Use greedy method to 
		for ( i = 0 ; i < ecnt ; ++i )
			eused[i] = true ;
		prevObj = LP( ecnt, intervalCnt, startCnt, endCnt, data, useExonFrag ) ;
		if ( prevObj != -1 )
			prevObj += IntervalDistance( readCnt, ecnt, cpCnt, data ) ;
		else
			prevObj = -1 ;
		memset( startCover, 0, sizeof( startCover[0] ) * ecnt * 3 ) ;
		memset( endCover, 0, sizeof( endCover[0] ) * ecnt * 3  ) ;
		
		for ( i = 0 ; i < ecnt; ++i )
		{
			int startStrand = start[ regionExons[i].startInd ].strand ;
			int endStrand = end[ regionExons[i].endInd ].strand ;

			if ( startStrand >= -1 )
				++startCover[ regionExons[i].startInd * 3 + startStrand + 1] ;
			else
				startCover[ regionExons[i].startInd * 3 + 0] = MAX_POINT ;

			if ( endStrand >= -1 )
				++endCover[ regionExons[i].endInd * 3 + endStrand + 1] ;
			else
				endCover[ regionExons[i].endInd * 3 + 0] = MAX_POINT ;
		}

		// compute the two long exons, they can not be removed.
		for ( j = endCnt - 1 ; j >= 0 ; --j )
		{
			i = 0 ;
			if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
			{
				for ( k = 0 ; k < ecnt ; ++k )
				{
					if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
					{
						//printf( "Long: %d %d %d\n", i, j, k ) ;
						longExons[0] = k ;
						break ;
					}
				}

				break ;
			}
		}

		for ( i = 0 ; i < startCnt ; ++i )
		{
			j = endCnt - 1 ;
			if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
			{
				for ( k = 0 ; k < ecnt ; ++k )
				{
					if ( regionExons[k].startInd == i && regionExons[k].endInd == j )
					{
						//printf( "Long: %d %d %d\n", i, j, k ) ;
						longExons[1] = k ;
						break ;
					}
				}

				break ;
			}
		}

		while ( 1 )
		{
			minObj = -1 ;
			int startStrand, endStrand ;
			for ( i = 0 ; i < ecnt ; ++i )
			{
				if ( regionExons[i].inEvidence )//|| i == longExons[0] || i == longExons[1] )
					continue ;
				startStrand = start[ regionExons[i].startInd ].strand ;
				endStrand = end[ regionExons[i].endInd ].strand ;
				if ( startStrand < -1 )
					startStrand = -1 ;
				if ( endStrand < -1 )
					endStrand = -1 ;

				if ( eused[i] && startCover[ regionExons[i].startInd * 3 + startStrand + 1] > 1 
						&& endCover[ regionExons[i].endInd * 3 + endStrand + 1] > 1 )
				{
					eused[i] = false ;
					if ( VerifyConstraintPoints( ecnt, cpCnt, mainStrand, data ) == false )
					{
						eused[i] = true ;
						continue ;
					}
					obj = LP( ecnt, intervalCnt, startCnt, endCnt, data, useExonFrag ) ;
					//printf( "#%lf %lf\n", obj, IntervalDistance( readCnt, ecnt, cpCnt, data ) ) ;
					if ( obj != -1 )
					{
						obj += IntervalDistance( readCnt, ecnt, cpCnt, data ) ;
						if ( minObj == -1 || obj < minObj ||
								( obj == minObj && 
								  regionExons[i].end - regionExons[i].start > regionExons[mintag].end - regionExons[mintag].start ) )
						{
							minObj = obj ;
							mintag = i ;
						}
					}
					eused[i] = true ;	
				}
			}

			//printf( "%.8lf %.8lf %d %d\n", prevObj, minObj, regionExons[mintag].start + OFFSET, regionExons[mintag].end + OFFSET ) ;
			//if ( minObj == -1 || ( prevObj >= 0.01 * intervalCnt && minObj > prevObj + 1 ) || 
			//	( prevObj < 0.01 * intervalCnt && minObj >= 0.01 * intervalCnt) )
			if ( minObj == -1 || minObj > prevObj + 1e-6 )
			{
				minObj = prevObj ;
				break ;
			}

			prevObj = minObj ;
			eused[ mintag ] = false ;
			startStrand = start[ regionExons[ mintag ].startInd ].strand ;
			endStrand = end[ regionExons[ mintag ].endInd ].strand ;
			if ( startStrand < -1 )
				startStrand = -1 ;
			if ( endStrand < -1 )
				endStrand = -1 ;

			--startCover[ regionExons[mintag].startInd * 3 + startStrand + 1] ;
			--endCover[ regionExons[mintag].endInd * 3 + endStrand + 1] ;
		} // end while

		int greedyChoose = 0, directChoose = 0 ;
		for ( i = 0 ; i < ecnt ; ++i )
		{
			bestEused[i] = eused[i] ; 
			//if ( eused[i] )
			//{
			//printf( "%d %d\n", regionExons[i].strand ) ;
			//}

			if ( eused[i] )
				++greedyChoose ;
		}
		//printf( "%d\n", greedyChoose ) ;
		if ( 1 ) //greedyChoose >= 16 )
		{
			//double directBuildScore = -1 ;
			double directBuildScore = DirectBuildExons( readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, data ) ;
			for ( i = 0 ; i < ecnt ; ++i )
				if ( eused[i] )
					++directChoose ;
			if ( minObj == -1 || minObj > 1 * intervalCnt || directBuildScore < minObj + 1e-6 || greedyChoose >= 2 * directChoose )
			{
				for ( i = 0 ; i < ecnt ; ++i )
				{
					bestEused[i] = eused[i] ;
				}
			}
			//printf( "Complicated: %lf %lf: %d\n", minObj, directBuildScore, greedyChoose ) ;
		}
		else
		{
			//if ( OFFSET == 156906608 )
			//	printf( "%d\n", VerifyConstraintPoints() ) ;
			bool  *fixedEused = (bool *)malloc( sizeof( bool ) * ecnt ) ;
			memset( fixedEused, false ,sizeof( bool ) * ecnt ) ;
			for ( i = 0 ; i < ecnt ; ++i ) 
			{
				if ( !eused[i] )
					fixedEused[i] = true ;

			}
			//if ( OFFSET == 156906608 )
			//	printf( "%d\n", VerifyConstraintPoints() ) ;
			memset( startCover, -1, sizeof( int ) * ecnt * 3 ) ;
			memset( endCover, -1, sizeof( int ) * ecnt * 3 ) ;
			memset( eused, false, sizeof( bool ) * ecnt ) ;


			//for ( i = 0 ; i < ecnt ; ++i )
			//	printf( "( %d %d)\n", regionExons[i].start + OFFSET, regionExons[i].end + OFFSET ) ;
			for ( k = 0 ; k < ecnt ; ++k )
			{
				if ( fixedEused[k] )
					continue ;

				i = regionExons[k].startInd ;
				j = regionExons[k].endInd ;

				if ( start[i].strand >= -1 )
					startCover[i * 3 + start[i].strand + 1] = 0 ;
				if ( end[j].strand >= -1 )
					endCover[j * 3 + end[j].strand + 1] = 0 ;
			}

			bestScore = MAX_LENGTH ;
			//EnumerateExon( 0, readCnt, fixedEused ) ;
			EnumerateExon( 0, readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, mainStrand, 
				bestScore, bestEused, data, fixedEused ) ;
			//printf( "%d %lf\n", start[0].pos + OFFSET, bestScore ) ;
			if ( bestScore > 1 * intervalCnt )
			{
				DirectBuildExons( readCnt, ecnt, intervalCnt, startCnt, endCnt, cpCnt, data ) ;
				for ( i = 0 ; i < ecnt ; ++i )
					bestEused[i] = eused[i] ;
			}
			free( fixedEused ) ;
		}
		/*int tmp = 0 ;
		  for ( i = 0 ; i < ecnt ; ++i )
		  {
		  if ( bestEused[i] )
		  ++tmp ;
		  }
		  prntf( "%d\n", tmp ) ; */
		//return 0 ;
	}

	/*for ( i = 0 ; i < ecnt ; ++i )
	  {
	  if ( bestEused[i] )
	  printf( "%d %d\n", regionExons[i].start + OFFSET, regionExons[i].end + OFFSET ) ;
	  }*/
	//eused[0] = eused[2] = eused[3] = eused[5] = true ;
	//LP() ;
	//fclose( fpComb ) ;
	//printf( "%d %d\n", startCnt, endCnt ) ;
	
	//Add lock here
	pthread_mutex_lock( &allExonsMutex ) ;
	//printf( "OFFSET=%d\n", OFFSET ) ;
	AddAllExons( OFFSET, ecnt, siCnt, bestEused, regionExons, spliceIndices, allExons, exonCnt ) ;
	pthread_mutex_unlock( &allExonsMutex ) ;

	// free the memory.
	free( bestEused ) ;

	free( eused ) ;
	free( graph ) ;
	free( visited ) ;
	free( beyond ) ;
	free( exonfrag ) ;
	free( efid ) ;
	free( leftEf ) ;
	free( rightEf ) ;
	free( regionExons ) ;
	free( exonHasSoft ) ;
	free( exonPath ) ;
	free( exonDist );
	free( constraintPoints ) ;
	free( start ) ;
	free( end ) ;
	free( startCover ) ;
	free( endCover ) ;
	free( spliceIndices ) ;
	free( point ) ;
	free( interval ) ;
	
	free( pArg ) ;

	pthread_mutex_lock( &solveRegionMutex ) ;
	--currSolveRegionThreadsCnt ;
	if ( currSolveRegionThreadsCnt == NUM_OF_THREADS - 1 )
	{
		pthread_cond_signal( &idleSolveRegionCond ) ;
	}
	if ( currSolveRegionThreadsCnt == 0 )
	{
		pthread_cond_signal( &clearSolveRegionCond ) ;
	}
	pthread_mutex_unlock( &solveRegionMutex ) ;

	pthread_exit( NULL ) ;
}

int SolveRegion( char *chrom,  
		struct _point inputStart[], struct _point inputEnd[], int inputStartCnt, int inputEndCnt, 
		struct _spliceIndex spliceIndices[], int siCnt, 
		int *depth, 
		struct _read reads[], int readCnt,
		struct _evidence *evidences, int eviCnt,
		struct _exon allExons[], int &exonCnt, int mainStrand, bool force )  
{
	int ecnt ;
	int startCnt, endCnt, intervalCnt ;
	int cpCnt ;
	struct _point start[MAX_POINT], end[MAX_POINT] ; // The points after merging
	int startCover[MAX_POINT][3], endCover[MAX_POINT][3] ;
	double sum ;

	int i, j, k ;
	int a, b ;
	int max = -1 ;
	double stat[10] ; // An array for statistics
	int len ; // The length of the region
	int mergeFactor = 0 ;
	int point[MAX_POINT] ;
	double interval[MAX_POINT] ;
	
	//printf( "### %s %d %d\n", chrom, inputStart[0].pos, inputEnd[ inputEndCnt - 1 ].pos ) ;
	//fpComb = fopen( "one_right_comb.out", "r" ) ;
	allStart = inputStart ;
	allEnd = inputEnd ;
	allStartCnt = inputStartCnt ;
	allEndCnt = inputEndCnt ;
	// TODO: Make this more efficient.
	//qsort( allStart, startCnt, sizeof( allStart[0] ), CompPoint ) ;
	//qsort( allEnd, endCnt, sizeof( allEnd[0] ), CompPoint ) ;

	memcpy( buffer, spliceIndices, sizeof( struct _spliceIndex ) * siCnt ) ;
	while ( 1 )
	{
		
		++mergeFactor ;
		if ( mergeFactor > 3 )
			break ;
		memcpy( spliceIndices, buffer, sizeof( struct _spliceIndex ) * siCnt ) ;
		// Merge the start and end points.
		// Consider the two strands respectively.
		for ( k = 0 ; k <= 1 ; ++k )
		{
			MergePoints( inputStart, inputStartCnt, k, mergeFactor ) ;
			MergePoints( inputEnd, inputEndCnt, k, mergeFactor ) ;		
		}
		//printf( "Merge Factor: %d\n", mergeFactor ) ;
		UpdateSpliceIndices( spliceIndices, siCnt, inputStart, inputEnd, inputStartCnt, inputEndCnt ) ;

		// Extract the merged points
		k = 0 ;
		for ( i = 0 ; i < inputStartCnt ; ++i )
		{
			if ( inputStart[i].merge != i )
			{
				SplicesInformation_Add( inputStart[i].pos, 
						inputStart[i].strand,
						inputStart[ GetMergeFather( inputStart, inputStart[i].merge) ].pos ) ;
				continue ;
			}
			start[k] = inputStart[i] ;
			//start[k].strand = 0 ;
			++k ;
		}
		startCnt = k ;

		k = 0 ;
		for ( i = 0 ; i < inputEndCnt ; ++i )
		{
			if ( inputEnd[i].merge != i )
			{
				SplicesInformation_Add( inputEnd[i].pos, 
						inputEnd[i].strand,
						inputEnd[ GetMergeFather( inputEnd, inputEnd[i].merge ) ].pos ) ;
				continue ;
			}
			end[k] = inputEnd[i] ;
			++k ;
		}
		endCnt = k ;	
		// Finish merging

		// Each splice site should have at least one "pair".
		int softEnd[2] = { -1, -1 }, softStart[2] = {-1, -1} ;
		for ( i = 0 ; i < startCnt ; ++i )
		{
			for ( j = 0 ; j < endCnt ; ++j )
			{
				if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand ) )
					break ;
			}

			if ( j >= endCnt && start[i].strong != 0 )
			{
				for ( k = 0 ; k < startCnt ; ++k )
					if ( start[k].strand != start[i].strand )
						break ;

				if ( i < k && k < startCnt )
					softEnd[ start[i].strand ] = start[i + 1].pos + READS_LENGTH - 1 ;
				else
					softEnd[ start[i].strand ] = start[i].pos + READS_LENGTH ; // TODO: Make it formal.
				//if ( softEnd[ start[i].strand ] > end[ endCnt - 1].pos )
				//	softEnd[]
			}
		}

		for ( j = endCnt - 1 ; j >= 0 ; --j )
		{
			for ( i = 0 ; i < startCnt ; ++i )
			{
				if ( end[j].pos > start[i].pos && IsSameStrand( end[j].strand, start[i].strand) )
					break ;
			}

			if ( i >= startCnt && end[j].strong != 0 )
			{

				for ( k = endCnt - 1 ; k >= 0 ; --k )
					if ( end[k].strand != end[j].strand )
						break ;
				if ( j > k && k >= 0 )
					softStart[ end[j].strand ] = end[j - 1].pos - READS_LENGTH + 1 ;
				else
					softStart[ end[j].strand ] = end[j].pos - READS_LENGTH ;
			}
		}
		
		// Speical case, splice junction in a region looks like: [[ )). Then the soft boundary above may make it look like [[ ] ()), which builds a gap in the middle contradicts the definition of region.
		if ( startCnt && endCnt && start[ startCnt - 1 ].pos < end[0].pos && 
				start[0].strand >= -1 && end[0].strand >= -1 && !IsSameStrand( start[0].strand,  end[0].strand ) )
		{
			for ( i = 1 ; i < startCnt ; ++i )
				if ( start[i].strand != start[i - 1].strand )
					break ;
			for ( j = 1 ; j < endCnt ; ++j )
				if ( end[j].strand != end[j - 1].strand )
					break ;
			if ( i >= startCnt && j >= endCnt )
			{
				softEnd[ start[0].strand ] = end[0].pos ;
				softStart[ end[0].strand ] = start[ startCnt - 1].pos ;
			}
		}

		// Insert the soft boundary back.
		for ( k = 0 ; k <= 1 ; ++k )
		{
			if ( softStart[k] == -1 )
				continue ;
			for ( i = startCnt - 1 ; i >= 0 && start[i].pos > softStart[k] ; --i )
				start[i + 1] = start[i] ;
			start[i + 1].pos = softStart[k] ;
			start[i + 1].strand = k ;
			start[i + 1].support = 0 ;
			start[i + 1].strong = 0 ;

			// Update the spliceIndices.
			for ( j = 0 ; j < siCnt ; ++j )
				if ( spliceIndices[j].start >= i + 1 )
					++spliceIndices[j].start ;

			++startCnt ;
		}

		for ( k = 0 ; k <= 1 ; ++k )
		{
			if ( softEnd[k] == -1 )
				continue ;
			for ( i = endCnt - 1 ; i >= 0 && end[i].pos > softEnd[k] ; --i )
				end[i + 1] = end[i] ;
			end[i + 1].pos = softEnd[k] ;
			end[i + 1].strand = k ;
			end[i + 1].support = 0 ;
			end[i + 1].strong = 0 ;

			for ( j = 0 ; j < siCnt ; ++j )
				if ( spliceIndices[j].end >= i + 1 )
					++spliceIndices[j].end ;
			++endCnt ;
		}

		len = end[endCnt - 1].pos - start[0].pos + 1 ;
		//printf( "## %d\n", len ) ;
		k = 0 ;
		j = 0 ;
		if ( start[0].pos < inputStart[0].pos )
		{
			for ( k = 0 ; k < inputStart[0].pos - start[0].pos ; ++k )
				psum[k] = 0 ;
		}

		j = 0 ;
		if ( start[0].pos >= inputStart[0].pos )
			j = start[0].pos - inputStart[0].pos ;

		if ( depth[j] != 0 )	
			psum[k] = sqrt( depth[j] ) ;//log( (double)depth[j] ) / log(LOG_BASE);
		else
			psum[k] = 0 ;
		for ( j = j + 1 ; j < len ; ++j )
		{
			if ( depth[j] != 0 )
				psum[k + j] = psum[k + j - 1] + sqrt( depth[j] ) ; //log( (double)depth[j] ) / log(LOG_BASE) ;
			else
				psum[k + j] = psum[k + j - 1] ;
		}
		sum = psum[len - 1] ;	
		OFFSET = start[0].pos ;
		for ( i = 0 ; i < startCnt ; ++i )
		{
			start[i].pos -= OFFSET ; 
			point[i] = start[i].pos - 1 ; // Pay attention to the -1 here. Intervals sequence looks like [a,b],[b+1,c]...
		}

		for ( i = 0 ; i < endCnt ; ++i )
		{
			end[i].pos -= OFFSET ;	
			point[i + startCnt] = end[i].pos ;
		}
		qsort( point, startCnt + endCnt, sizeof( point[0] ), Comp ) ;		

		intervalCnt = startCnt + endCnt - 1 ;	
		
		
		interval[0] = psum[ point[1] ] ; 
		for ( i = 2 ; i < intervalCnt + 1 ; ++i )
			interval[i - 1] = ( psum[ point[i] ] - psum[ point[i - 1] ] ) ; /// (double)( point[i] - point[i - 1] + 1 ) ;
		
		memset( startCover, -1, sizeof( startCover ) ) ;
		memset( endCover, -1, sizeof( endCover ) ) ;
		ecnt = 0 ;
		for ( i = 0 ; i < startCnt ; ++i )
			for ( j = 0 ; j < endCnt ; ++j )
			{
				if ( end[j].pos <= start[i].pos || !IsSameStrand( end[j].strand, start[i].strand ) )
					continue ;

				regionExons[ ecnt ].start = start[i].pos ;
				regionExons[ ecnt ].end = end[j].pos ;
				regionExons[ ecnt ].startInd = i ;
				regionExons[ ecnt ].endInd = j ;
				regionExons[ ecnt ].strand = start[i].strand >= end[j].strand ? start[i].strand : end[j].strand ;
				regionExons[ ecnt ].inEvidence = false ;

				if ( !start[i].strong )
					exonHasSoft[ ecnt ][0] = true ;
				else
					exonHasSoft[ ecnt ][0] = false ;

				if ( !end[j].strong )
					exonHasSoft[ ecnt ][1] = true ;
				else
					exonHasSoft[ ecnt ][1] = false ;
				if ( start[i].strand >= -1 )
					startCover[i][ start[i].strand + 1] = 0 ;
				if ( end[j].strand >= -1 )
					endCover[j][ end[j].strand + 1] = 0 ;
				++ecnt ;
			}
		//printf( "In SolveRegion: %d %d %d\n", ecnt, startCnt, endCnt ) ;
		if ( ecnt >= 32 && !force )
			return -1 ;
		if ( force && ecnt >= 32 )
		{
			continue ;
		}
		break ;
	}	
	/*printf( "//// %d %d\n", startCnt, endCnt ) ;
	for ( i = 0 ; i < startCnt ; ++i )
		printf( "%d: %d %d\n", i, start[i].pos + OFFSET, start[i].strand ) ;

	for ( i = 0 ; i < endCnt ; ++i )
		printf( "# %d: %d %d\n", i, end[i].pos + OFFSET, end[i].strand ) ;

	for ( i = 0 ; i < siCnt ; ++i )
		printf( "si %d: %d %d (%d %d)\n", i, spliceIndices[i].start, spliceIndices[i].end, spliceIndices[i].endPos, spliceIndices[i].startPos ) ;*/
	
	// Very simple region
	if ( ecnt <= 2 )
	{
		bool bestEused[2] = {true, true} ;
		pthread_mutex_lock( &allExonsMutex ) ;
		//printf( "OFFSET=%d\n", OFFSET ) ;
		AddAllExons( OFFSET, ecnt, siCnt, bestEused, regionExons, spliceIndices, allExons, exonCnt ) ;
		pthread_mutex_unlock( &allExonsMutex ) ;
		return 0 ;
	}

	AreExonsInEvidences( ecnt, evidences, eviCnt ) ;
	//for ( i = 0 ; i < startCnt ; ++i )
	//	printf( "%d %d %d\n", startCover[i][0], startCover[i][1], startCover[i][2] ) ;

	/*for ( i = 0 ; i < startCnt ; ++i )
		printf( "[ %d %d(%d)\n", start[i].pos + OFFSET, start[i].strand, start[i].strong ) ;
	for ( i = 0 ; i < endCnt ; ++i )
		printf( "] %d %d(%d)\n", end[i].pos + OFFSET, end[i].strand, end[i].strong ) ;
	printf( "%d\n", ecnt ) ;
	  for ( i = 0 ; i < ecnt ; ++i )
	  printf( "possible exons: %d %d %d\n", regionExons[i].start + OFFSET, regionExons[i].end + OFFSET, regionExons[i].strand ) ;*/


	//if ( ecnt >= 30 )
	//	printf( "%s %d %d %d\n", buffer, start[0] + OFFSET, end[ endCnt - 1 ] + OFFSET, ecnt ) ;
	//return 0 ;

	//IntervalStatistic( stat, depth, len ) ;
	//ReadRightCombination( fpComb ) ;
	//BuildExonGraph( spliceIndices, siCnt ) ;	
	GetConstraintPoints( start, startCnt, end, endCnt, point, intervalCnt, reads, readCnt, constraintPoints, cpCnt ) ;
	//VerifyConstraintPoints( true ) ; // The verification will be done by the threads 

	// Use pthread here
	struct _pthreadArgSolveRegion *arg = (struct _pthreadArgSolveRegion *)malloc( sizeof( *arg ) ) ;
	arg->regionExons = ( struct _enumExon *)malloc( sizeof( struct _enumExon ) * ecnt ) ;
	memcpy( arg->regionExons, regionExons, sizeof( regionExons[0] ) * ecnt ) ;
	arg->ecnt = ecnt ;
	arg->exonHasSoft = ( bool * )malloc( sizeof( bool ) * 2 * ecnt ) ;
	memcpy( arg->exonHasSoft, exonHasSoft, sizeof( bool ) * 2 * ecnt ) ;

	arg->point = ( int *)malloc( sizeof( point[0] ) * ( intervalCnt + 2 ) ) ;
	memcpy( arg->point, point, sizeof( point[0] ) * ( intervalCnt + 2 ) ) ;
	arg->interval = (double *)malloc( sizeof( interval[0] ) * intervalCnt ) ;
	memcpy( arg->interval, interval, sizeof( interval[0] ) * intervalCnt ) ;
	arg->intervalCnt = intervalCnt ;
	arg->start = ( struct _point *)malloc( sizeof( start[0] ) * startCnt ) ;
	memcpy( arg->start, start, sizeof( start[0] ) * startCnt ) ;
	arg->startCnt = startCnt ;
	arg->end = ( struct _point *)malloc( sizeof( end[0] ) * endCnt ) ;
	memcpy( arg->end, end, sizeof( end[0] ) * endCnt ) ;
	arg->endCnt = endCnt ;

	arg->startCover = ( int *)malloc( sizeof( int ) * ecnt * 3 ) ;
	for ( i = 0 ; i < startCnt ; ++i )
	{
		for ( j = 0 ; j < 3 ; ++j )
			arg->startCover[i * 3 + j] = startCover[i][j] ; 
	}
	arg->endCover = ( int *)malloc( sizeof( int ) * ecnt * 3 ) ;
	for ( i = 0 ; i < endCnt ; ++i )
	{
		for ( j = 0 ; j < 3 ; ++j )
			arg->endCover[i * 3 + j] = endCover[i][j] ; 
	}

	arg->spliceIndices = ( struct _spliceIndex * )malloc( sizeof( spliceIndices[0] ) * siCnt ) ;
	memcpy( arg->spliceIndices, spliceIndices, sizeof( spliceIndices[0] ) * siCnt ) ;
	arg->siCnt = siCnt ;

	arg->constraintPoints = ( struct _constraintPoint *)malloc( sizeof( constraintPoints[0] ) * cpCnt ) ;
	memcpy( arg->constraintPoints, constraintPoints, sizeof( constraintPoints[0] ) * cpCnt ) ;
	arg->cpCnt = cpCnt ;

	arg->allExons = allExons ;
	arg->exonCnt = &exonCnt ;

	arg->OFFSET = OFFSET ;
	arg->mainStrand = mainStrand ;
	arg->readCnt = readCnt ;
	//return 0 ;

	// wait for idle thread
	pthread_mutex_lock( &solveRegionMutex ) ;
	if ( currSolveRegionThreadsCnt >= NUM_OF_THREADS )
		pthread_cond_wait( &idleSolveRegionCond, &solveRegionMutex ) ;
	++currSolveRegionThreadsCnt ;
	pthread_t thread ;
	pthread_create( &thread, &pthreadAttr, SolveRegion_Thread, (void *)arg ) ;
	
	// SERIAL
	//pthread_cond_wait( &clearSolveRegionCond, &solveRegionMutex ) ;
	
	pthread_mutex_unlock( &solveRegionMutex ) ;
	//SolveRegion_Thread( (void *)arg ) ;
	return 0 ;
}

int SolveRegion_Wrapper( char *chrom, struct _splice splices[], int scnt, int softStart, int softEnd, int *depth, 
		struct _read reads[], int readCnt, struct _evidence *evidences, int eviCnt, struct _exon allExons[], int &exonCnt,
 		int minSpliceSupportPlus, int minSpliceSupportMinus, bool force )
{
	int i, j, k ;
	int startCnt, endCnt ;
	int leftCnt, rightCnt ;
	struct _point *start, *end ;
	int siCnt ;
	int map[2 * MAX_POINT] ;
	int mainStrand ;

	static struct _TTSsite *tts ;
	static int ttsCnt ;
	static int ttsTag = 0 ; // Current TTS site.	
	static char prevTTSChrom[10] = "-1" ;

	startCnt = 0 ;
	endCnt = 0 ;

	start = inputStart ;
	end = inputEnd ;
	
	// Adjust the soft boundary
	// due to the weird coverage on 3',5' UTR
	if ( softStart != -1 )
	{
		int max = -1 ;
		int newStart = -1 ;
		for ( i = softStart ; i < splices[0].pos ; ++i )
			if ( depth[i - softStart ] > max )
				max = depth[ i - softStart ] ;
		
		if ( splices[0].pos - READS_LENGTH > softStart )
			newStart = splices[0].pos - READS_LENGTH ;
		for ( i = softStart ; i < splices[0].pos - READS_LENGTH ; ++i )
			if (  4 * depth[i - softStart] / (double)max >= (double)max / depth[i - softStart] )
			{
				newStart = i ;
				break ;
			}
		//if ( softStart == 198608105 )
		//{
		//	printf( "hi %d %d: %d-%d\n", max, newStart, softStart, splices[0].pos ) ;
		//}
		if ( newStart != -1 )
		{
			softStart = newStart ;
			// Don't forget to shift the depth array.
			depth = depth + newStart - softStart ;
		}
	}
	if ( softEnd != -1 )
	{
		int offset = softStart != -1 ? softStart : splices[0].pos ; 
		int max = -1 ;
		int newEnd = -1 ;
		for ( i = splices[scnt - 1].pos ; i <= softEnd ; ++i )
		{
			if ( depth[i - offset ] > max )
				max = depth[ i - offset ] ;
		}
		if ( softEnd > splices[scnt - 1].pos + READS_LENGTH )
			newEnd = splices[scnt - 1].pos + READS_LENGTH ;
		for ( i = softEnd ; i > splices[scnt - 1].pos + READS_LENGTH ; --i )
			if (  4 * depth[i - offset] / (double)max >= (double)max / depth[i - offset] )
			{
				newEnd = i ;
				break ;
			}
		if ( newEnd != -1 )
		{
			softEnd = newEnd ;
		}
	}

	if ( softStart != -1 )
	{
		start[0].pos = softStart ;
		start[0].strand = -1 ; //splices[0].strand ;
		start[0].support = 0 ;	
		start[0].strong = 0 ;
		++startCnt ;
	}
	
	tts = TTS_GetCurrentTTS( ttsCnt ) ;
	if ( ttsCnt > 0 ) //&& ( softEnd != -1 || softStart != -1 ) )
	{
		// Adjust ttsTAg ;
		if ( strcmp( prevTTSChrom, chrom ) )
		{
			ttsTag = 0 ;
			strcpy( prevTTSChrom, chrom ) ;
		}

		//printf( "hi1\n" ) ;  
		int rstart = -1, rend ;

		if ( softStart != -1 )
			rstart = softStart ;
		else
		{
			for ( i = 0 ; i < scnt ; ++i )
				if ( splices[i].type == 1 && splices[i].pos != -1 && ( rstart == -1 || splices[i].pos < rstart ) )
					rstart = splices[i].pos ;
		}

		if ( softEnd != -1 )
			rend = softEnd ;
		else
		{
			for ( i = 0 ; i < scnt ; ++i )
				if ( splices[i].type == 0 && splices[i].pos != -1 && splices[i].pos > rend )
					rend = splices[i].pos ;
		}

		/*if ( rstart == -1 )
		{
			printf( "%d\n", rstart ) ;
			for ( i = 0 ; i < scnt ; ++i )
				printf( "%d %d\n", splices[i].type, splices[i].pos ) ;
			exit( 1 ) ;
		}*/

		for ( ; ttsTag < ttsCnt ; ++ttsTag )
		{
			if ( tts[ttsTag].pos > rend )
				break ;
			if ( tts[ttsTag].pos < rstart )
				continue ;
			
			if ( tts[ttsTag].sigPos < rstart || tts[ttsTag].sigPos > rend )
				continue ;
			// pos += 30 ;
			//printf( "%s (%d %d): %d\n", chrom, rstart, rend, tts[ ttsTag ].pos ) ;
			//printf( "#polya: %s %d\n", chrom, pos ) ;
			if ( /*softEnd != -1 &&*/ tts[ ttsTag ].strand == 1 )
			{
				end[ endCnt ].pos = tts[ttsTag].pos ; 
				end[ endCnt ].strand = -2 ;//1 ;
				end[ endCnt ].support = 0 ;
				end[ endCnt ].strong = 0 ;
				++endCnt ;
			}
			else if ( /*softStart != -1 &&*/ tts[ttsTag].strand == 0  )
			{
				start[ startCnt ].pos = tts[ttsTag].pos ;
				start[ startCnt ].strand = -3 ; //0 ;
				start[ startCnt ].support = 0 ;
				start[ startCnt ].strong = 0 ;
				++startCnt ;
			}
		}
		//printf( "hi2 %d %d: %d %d, %d %d\n", rstart, rend, startCnt, start[0].pos, endCnt, end[0].pos ) ;
	}
	
	/*for ( i = 0 ; i < scnt ; ++i )
		if ( splices[i].type == 0 )
			break ;
	for ( j = i ; j < scnt ; ++j )
		if ( splices[j].type == 1 )
			break ;*/
	if ( scnt >= 2 ) //&& ( i == scnt - 1 || ( i == 1 && j >= scnt ) ) )
	{	
		// Test the unbalanced exon
		int leftCnt, rightCnt, insertPos ;
		
		/*if ( i == scnt - 1 )
			insertPos = ( splices[ scnt - 2 ].pos + splices[ scnt - 1 ].pos ) / 2 ;
		else if ( i == 1 )
		{
			insertPos = ( splices[0].pos + splices[1].pos ) / 2 ;
		}*/

		leftCnt = 0 ;
		rightCnt = 0 ;
		for ( i = 0 ; i < scnt ; ++i )
			if ( splices[i].type == 0 && splices[i].otherInd >= scnt  )
				rightCnt += splices[i].support ;
			else if ( splices[i].type == 1 && splices[i].otherInd < 0 )
				leftCnt += splices[i].support ;
		//if ( rightCnt > 100 * leftCnt && softStart == -1 )
		if ( SignificantGreater( rightCnt, leftCnt ) && softStart == -1 )
		{
			insertPos = ( splices[0].pos + splices[1].pos ) / 2 ;
			start[ startCnt ].pos = insertPos ;
			start[ startCnt ].strand = -1 ;
			start[ startCnt ].support = 0 ;
			start[ startCnt ].strong = 0 ;
			++startCnt ;
		}

		//if ( leftCnt > 100 * rightCnt && softEnd == -1 )
		if ( SignificantGreater( leftCnt, rightCnt ) && softEnd == -1  )
		{
			insertPos = ( splices[ scnt - 2 ].pos + splices[ scnt - 1 ].pos ) / 2 ;
			end[ endCnt ].pos = insertPos ;
			end[ endCnt ].strand = -1 ;
			end[ endCnt ].support = 0 ;
			end[ endCnt ].strong = 0 ;
			++endCnt ;
		}
	}

	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].pos == -1 )
			continue ;
		if ( splices[i].otherInd >= 0 && splices[i].otherInd < scnt && 
			( ( splices[i].strand == 1 && splices[i].support < minSpliceSupportPlus ) ||
			  ( splices[i].strand == 0 && splices[i].support < minSpliceSupportMinus ) ) )
			continue ;
		
		if ( splices[i].type == 1 )
		{
			for ( j = 0 ; j < startCnt ; ++j )
				if ( start[j].pos == splices[i].pos && IsSameStrand( start[j].strand, splices[i].strand ) )
					break ;
			if ( j < startCnt )
			{
				map[i] = j ;
				if ( start[j].strong == 0 )
				{
					start[ j ].pos = splices[i].pos ;
					start[ j ].strand = splices[i].strand ;
					start[ j ].support = splices[i].support ;
					start[ j ].strong = 1 ;
				}
				else
					start[j].support += splices[i].support ;
				continue ;
			}	

			for ( j = startCnt - 1 ; j >= 0 ; --j )
			{
				if ( start[j].pos < splices[i].pos )
					break ;
				else
					start[j + 1] = start[j] ;
			}
			start[ j + 1 ].pos = splices[i].pos ;
			start[ j + 1 ].strand = splices[i].strand ;
			start[ j + 1 ].support = splices[i].support ;
			start[ j + 1 ].strong = 1 ;
			// NOTE: There is no need to adjust other splices' map, because splices are sorted.
			map[i] = j + 1 ; // startCnt
			++startCnt ;
		}
		else if ( splices[i].type == 0 )
		{
			for ( j = 0 ; j < endCnt ; ++j )
				if ( end[j].pos == splices[i].pos && IsSameStrand( end[j].strand, splices[i].strand ) )
					break ;
			if ( j < endCnt )
			{
				map[i] = j ;
				if ( end[j].strong == 0 )
				{
					end[ j ].pos = splices[i].pos ;
					end[ j ].strand = splices[i].strand ;
					end[ j ].support = splices[i].support ;	
					end[ j ].strong = 1 ;	
				}
				else
					end[j].support += splices[i].support ;
				continue ;
			}

			for ( j = endCnt - 1 ; j >= 0 ; --j )
			{
				if ( end[j].pos < splices[i].pos )
					break ;
				else
					end[j + 1] = end[j] ;
			}
			end[ j + 1 ].pos = splices[i].pos ;
			end[ j + 1 ].strand = splices[i].strand ;
			end[ j + 1 ].support = splices[i].support ;	
			end[ j + 1 ].strong = 1 ;	
			map[i] = j + 1 ;	
			++endCnt ;					
		}
	}



	if ( softEnd != -1 )
	{
		end[ endCnt ].pos = softEnd ;
		end[ endCnt ].strand = -1 ; //splices[ k - 1 ].strand ;
		end[ endCnt ].support = -1 ;
		end[ endCnt ].strong = 0 ;
		++endCnt ;
	}
	/*for ( i = 0 ; i < startCnt ; ++i )
		printf( "%d %d %d\n", inputStart[i].pos, inputStart[i].strand, inputStart[i].strong ) ;
	for ( i = 0 ; i < endCnt ; ++i )
		printf( "# %d %d %d\n", inputEnd[i].pos, inputEnd[i].strand, inputEnd[i].strong ) ;
	printf( "===\n" ) ;*/
	if ( startCnt == 0 && endCnt == 0 )
		return 0 ;
	
	for ( i = 0 ; i < startCnt ; ++i )
		start[i].merge = i ;
	for ( i = 0 ; i < endCnt ; ++i )
		end[i].merge = i ;

	//printf( "### %d %d %d %d\n", start[0].pos, end[ endCnt - 1].pos, start[0].strand, end[ endCnt - 1].strand ) ;	
	// Extract the splice indices.
	siCnt = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].pos == -1 )
			continue ;
		if ( splices[i].otherInd >= 0 && splices[i].otherInd < scnt && 
			( ( splices[i].strand == 1 && splices[i].support < minSpliceSupportPlus ) ||
			  ( splices[i].strand == 0 && splices[i].support < minSpliceSupportMinus ) ) )
			continue ;

		if ( splices[i].type == 0 && splices[i].otherInd < scnt )
		{
			spliceIndices[ siCnt ].start = map[ splices[i].otherInd  ] ;
			spliceIndices[ siCnt ].end = map[i] ;
			spliceIndices[ siCnt ].strand = splices[i].strand ;
			spliceIndices[ siCnt ].startPos = splices[i].otherPos ;
			spliceIndices[ siCnt ].endPos = splices[i].pos ;
			++siCnt ;
		}
		else if ( splices[i].type == 0 && splices[i].otherInd >= scnt )
		{
			spliceIndices[ siCnt ].start = -1 ;			
			spliceIndices[ siCnt ].end = map[i] ;
			spliceIndices[ siCnt ].strand = splices[i].strand ;
			spliceIndices[ siCnt ].startPos = splices[i].otherPos ;
			spliceIndices[ siCnt ].endPos = splices[i].pos ;
			++siCnt ;
		}
		else if ( splices[i].type == 1 && splices[i].otherInd < 0 )
		{
			spliceIndices[ siCnt ].start = map[i] ;
			spliceIndices[ siCnt ].end = -1 ;
			spliceIndices[ siCnt ].strand = splices[i].strand ; 
			spliceIndices[ siCnt ].startPos = splices[i].pos ;
			spliceIndices[ siCnt ].endPos = splices[i].otherPos ;
			++siCnt ;
		} 
	}

	// Find out the main strand.
	mainStrand = -1 ;
	int plusSupport = 0, minusSupport = 0 ; 
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].strand == 1 )
			plusSupport += splices[i].support ;
		if ( splices[i].strand == 0 )
			minusSupport += splices[i].support ;
	}

	if ( 100 * plusSupport > 95 * ( plusSupport + minusSupport ) )
		mainStrand = 1 ;
	if ( 100 * minusSupport > 95 * ( plusSupport + minusSupport ) )
		mainStrand = 0 ;

	//printf( "### %d %d %d: %d\n", start[0].pos, end[ endCnt - 1 ].pos, readCnt, READS_LENGTH ) ;
	//if ( readCnt * READS_LENGTH / ( end - start + 1 ) < 1 )
	//	return 0 ;

	/*for ( i = 0 ; i < siCnt ; ++i )
	{
		printf( "si %d: %d %d (%d %d)\n", i, spliceIndices[i].start, spliceIndices[i].end, spliceIndices[i].endPos, spliceIndices[i].startPos ) ;
	}*/
	return SolveRegion( chrom, start, end, startCnt, endCnt, spliceIndices, 
			siCnt, depth, reads, readCnt, evidences, eviCnt, allExons, exonCnt, mainStrand, force ) ;	
}
