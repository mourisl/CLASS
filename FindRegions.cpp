/**
  Find a region with alternative splice event based on the depth(coverage) file and the junc file generated from Tophat.
*/

#include "FindRegions.h"

#define MAX_CHROM 100
#define ISLAND_COVER 10 

#define ALL_CLEAN 0 

extern int READS_LENGTH ;
extern int FRAG_LENGTH, FRAG_STD ;
extern bool VERBOSE ;

int BIN_SIZE = 50 ;
int BIN_RANGE = 10000 ;
int BIN_USE = 30 ; // The number of bins on one side of the interval that we will use

double IR_ALPHA = 6.0 ;

#define MAX_CURRENT_ISLANDS 10000
#define BIN_ALPHA 6.0 

struct _chrom
{
	char name[10] ;
	//int end ; 
} ;

FILE *fpDepth, *fpSplice, *fpDepthForGroup ;

struct _splice *splices, *tmpSplices ; 
int scnt ;

// A gene with single exon. It is determined in the IsInExon function.
int singleExon[MAX_EXON][2], tmpSingleExon[MAX_EXON][2] ;
int seCnt = 0, tmpSeCnt = 0 ;
double noiseDepth, noiseSqDepth ;
int noiseCnt ; // the single exons regarded as noise. Computed in IsSingleExon

struct _chrom chroms[MAX_CHROM] ;
int chromCnt ;
int curChromInd ;

struct _readFile fpReads, fpReadsForGroup ;
struct _read cutReads[MAX_READ] ;
int cutReadCnt = 0 ;
int prevReadRegion[2] = { -1000, -1000 };

// The average coverage for the islands.
struct _islandCover
{
	double avgCover[MAX_CURRENT_ISLANDS] ;
	int pos[MAX_CURRENT_ISLANDS] ;
	int used ;
} ;
struct _islandCover prevIslandCover,  // The islands before and in the interval.
		futureIslandCover ;	
double islandDepth[MAX_LENGTH] ;

bool IsInExon( char *chrom, int spliceInd, bool isFuture ) ;

// The variables holding the exons of evidences.
struct _exon *evidenceExons ;
int eviExonCnt ;

// Some additional splice information used for finding regions.
struct _spliceInfo
{
	int geneId ;
	int regionId ;

	// The coverage from this splice to the next splice
	double depthSum ;
	double depthSqSum ;
	int depthCnt ;

	double avgSoftDepth[2] ;
	double softDepthSum[2] ;
	int soft[2] ;

	double threshold ; // The threshold decide whether a noise or a signal for the affiliated soft

	int possibleIR ; // if the splice is ] [, is it a possible intron retention? 0-no, 1-yes, 2-might be a forced IR
} ;

struct _spliceInfo *spliceInfo ;

extern int CompInt( const void *p1, const void *p2 ) ;
extern int CompDouble( const void *p1, const void *p2 ) ;
// Get the median of the n numbers in array a.
double GetMedian( double *a, int n )
{
	qsort( a, n, sizeof( a[0] ), CompDouble ) ;
	return a[n / 2] ;	
}

// Test whether there is a paired-end read spanning the cut junction site
bool NoPairReadsSpan( char *chrom, int cutStart, int cutEnd, int start, int end, bool isFuture, double CUT_THRESHOLD )
{
	if ( cutEnd - cutStart > FRAG_LENGTH + 2 * FRAG_STD - 2 * READS_LENGTH ) 
		return true ;
	if ( isFuture )
		return true ;
	int rstart, rend ;
	int i, k, offset[2] = {0,0} ;
	bool flag ;
	/*if ( start == 124839485 + 1 && end == 124839978 - 1 )
	{
		printf( "### Yoooo %d %d\n", cutStart, cutEnd ) ;
		exit( 1 ) ;	
	}*/
	rstart = cutStart - FRAG_LENGTH ;
	if ( rstart < start )
	{
		rstart = start ;
		offset[0] = MERGE_DISTANCE - 1 ;	
	}
	rend = cutEnd + FRAG_LENGTH ; 
	if ( rend > end )
	{
		rend = end ;
		offset[1] = MERGE_DISTANCE - 1 ;
	}
	//printf( "### %d %d %d\n", rstart, rend, prevReadRegion[1] ) ;
	cutReadCnt = ExtractReads( fpReads, chrom, rstart,
		rend, cutReads, prevReadRegion ) ;

	int support = 0 ;
	flag = false ;
	for ( i = 0 ; i < cutReadCnt ; ++i )
	{
		k = cutReads[i].mateInd ;
		if ( flag && cutReads[i].start == cutReads[i - 1].start )
			continue ;
		flag = false ;
		
		if ( k == -1 || k < i )
			continue ;

		if ( cutReads[i].end > rstart + offset[0] && cutReads[k].start + offset[1] < rend &&
		     cutReads[i].end < cutStart && cutReads[k].start > cutEnd &&  
		     cutReads[i].end < rend && cutReads[k].start > rstart &&cutReads[k].end - cutReads[i].start + 1 <= FRAG_LENGTH + 2 * FRAG_STD )
		{
//if ( start == 124839486 && end == 124839977 )
//	printf( "###### hi %d %d %d %d, %d %d,(%d,%d) (%d,%d)\n", cutStart, cutEnd, rstart, rend, i, k, cutReads[i].start, cutReads[i].end, cutReads[k].start, cutReads[k].end ) ;
			++support ;
			flag = true ;
			if ( log( (double)support ) / log( LOG_BASE ) >= CUT_THRESHOLD )
				return false ; 
		}
	}
	return true ;	
}

void GetFutureIslands( char *chrom, int spliceInd )
{
	fpos_t fpDepthPos, fpReadsPos ;
	int soft[2] ;
	int start = splices[ spliceInd ].pos ;
	int tmpPrevReadRegion[2] = { prevReadRegion[0], prevReadRegion[1] };

	fgetpos( fpDepth, &fpDepthPos ) ;
	//fgetpos( fpDepth, &fpReadsPos ) ;
	
	futureIslandCover.used = 0 ;
	while ( spliceInd < scnt - 1 )
	{
		IsInExon( chrom, spliceInd, true ) ;
		if ( futureIslandCover.used >= BIN_USE || splices[ spliceInd ].pos - start >= BIN_RANGE )
			break ;
		++spliceInd ;
	}
	prevReadRegion[0] = tmpPrevReadRegion[0] ;
	prevReadRegion[1] = tmpPrevReadRegion[1] ;

	fsetpos( fpDepth, &fpDepthPos ) ;
	//fsetpos( fpDepth, &fpReadsPos ) ;
}

// According to the bins to decide the bin_alpha
// type - 0: noise covered all the intron. 1: noise creates alternative TSS or TES. 2: Single exon gene
double GetBinAlpha( int l, int r, int type, int pos = -1 )
{
	//return 0 ;
	//printf( "%d %d %d\n", l, r, pos ) ;
	//return 2.0 ;

	if ( type == 0 )
	{
		//if ( l == BIN_USE || r == BIN_USE )
		//	return 3.0 ;
		//return 2.0 ;
		return IR_ALPHA ;
	}
	else if ( type == 1 )
	{
		//if ( l == BIN_USE || r == BIN_USE )
		//	return 4.0 ;
		if ( l >= 0.9 * BIN_USE || r >= 0.9 * BIN_USE )
			return 4.0 ;
		//if ( l >= 0.75 * BIN_USE || r >= 0.75 * BIN_USE )
		//	return 4.0 ;	
		if ( l >= 0.5 * BIN_USE || r >= 0.5 * BIN_USE )
			return 3.0 ;
		//if ( l >= 0.25 * BIN_USE || r >= 0.25 * BIN_USE )
		return 2.0 ;
	}
	else if ( type == 2 )
	{
		if ( l == BIN_USE || r == BIN_USE )
			return 6.0 ;
		if ( l >= 0.9 * BIN_USE || r >= 0.9 * BIN_USE )
			return 5.0 ;
		if ( l >= 0.75 * BIN_USE || r >= 0.75 * BIN_USE )
			return 4.0 ;	
		if ( l >= 0.5 * BIN_USE || r >= 0.5 * BIN_USE )
			return 3.0 ;
		//if ( l >= 0.25 * BIN_USE || r >= 0.25 * BIN_USE )
		return 2.0 ;
		//return 0 ;
	}
	//if ( l >= 3 || r >= 3 )
	//	return 2.0 ;
	return 0 ;
}

// Test whether an island can become a single exon.
bool IsSingleExon( int from, int to, const struct _islandCover &prevIslandCover, const struct _islandCover &currIslandCover, 
		const struct _islandCover &futureIslandCover )
{
	//if ( to - from + 1 < READS_LENGTH )
	//	return false ;
	//return true ;
	int i, j, k ;
	int a, b ;
	for ( i = 0 ; i < currIslandCover.used ; ++i )
	{
		if ( currIslandCover.pos[i] >= from )
			break ;
	}
	a = i - 1 ;
	for ( ; i < currIslandCover.used ; ++i )
	{
		if ( currIslandCover.pos[i] > to )
			break ;
	}
	b = i ;
	int l = 0, r = 0 ;
	double islandCoverMean, islandCoverVar ;
	double islandSum, islandSqSum ;
	double avgCover ;

	//printf( "%lf %d\n", islandSum, islandCnt ) ;
	//printf( "%d\n", start - 1 ) ;
	islandSum = 0 ;
	islandSqSum = 0 ;
	for ( i = a ; i >= 0 ; --i )
	{
		if ( l >= BIN_USE || currIslandCover.pos[i] < from - BIN_RANGE )
			break ;
		++l ;
		islandSum += currIslandCover.avgCover[i] ;
		islandSqSum += currIslandCover.avgCover[i] * currIslandCover.avgCover[i] ;
		//printf( "island: %lf\n", islandCover->avgCover[i] ) ;
	}
	for ( i = prevIslandCover.used - 1 ; i >= 0 ; --i )
	{
		if ( l >= BIN_USE || prevIslandCover.pos[i] < from - BIN_RANGE  )
			break ;
		++l ;
		islandSum += prevIslandCover.avgCover[i] ;
		islandSqSum += prevIslandCover.avgCover[i] * prevIslandCover.avgCover[i];
		//printf( "island: %lf\n", islandCover->avgCover[i] ) ;
	}

	r = 0 ;
	for ( i = b ; i < currIslandCover.used ; ++i )
	{
		if ( r >= BIN_USE || currIslandCover.pos[i] > to + BIN_RANGE )
			break ;
		++r ;
		islandSum += currIslandCover.avgCover[i] ;
		islandSqSum += currIslandCover.avgCover[i] * currIslandCover.avgCover[i] ;
		//printf( "island: %lf\n", islandCover->avgCover[i] ) ;
	}
	for ( i = 0 ; i < futureIslandCover.used ; ++i )
	{
		if ( r >= BIN_USE || futureIslandCover.pos[i] > to + BIN_RANGE )
			break ;
		++r ;
		islandSum += futureIslandCover.avgCover[i] ;
		islandSqSum += futureIslandCover.avgCover[i] * futureIslandCover.avgCover[i] ;
		//printf( "island: %lf\n", islandCover->avgCover[i] ) ;
	}
	if ( l + r > 4 )
	{ 
		islandCoverMean = islandSum / (double)( l + r ) ;
		islandCoverVar = sqrt( islandSqSum / (double)(l + r) - islandCoverMean * islandCoverMean )  ;
	}
	else
	{
		islandCoverMean = 10 ;
		islandCoverVar = 0 ;
	}
	
	avgCover = 0 ;
	for ( i = a + 1 ; i < b ; ++i )	
	{
		avgCover += currIslandCover.avgCover[i] ;
	}
	if ( b == a )
	{
		return false ;
	}
	avgCover /= ( b - a ) ; // Not accurate because of the last bin.
	//printf( "(%d %d): %lf %lf %lf\n", from, to, avgCover, islandCoverMean, islandCoverVar * GetBinAlpha( l, r, 2 ) ) ;
	if ( avgCover > islandCoverMean + islandCoverVar * GetBinAlpha( l, r, 2 ) )
	{
		//printf( "hi %d %d\n", from, to  ) ;
		return true ;
	}
	else
	{
		for ( i = a + 1 ; i < b ; ++i )
		{
			noiseDepth += currIslandCover.avgCover[i] ;
			noiseSqDepth += currIslandCover.avgCover[i] * currIslandCover.avgCover[i] ;
			++noiseCnt ;
		}
		return false ;
	}
}

/**
  Test wether a interval [start, end] is covered by some exons. Note: This interval is inside of a pair of splice sites. 
  If it can be cut, then it is possible to look like [ )...( ] the positions of )...( is returned in soft[2].
  spliceInd is the left splice indices for that interval.
*/
bool IsInExon( char *chrom, int spliceInd, bool isFuture )
{
	int i, j, k, tmp ;
	char inChrom[50] ;
	int pos = -1, idepth ;
	int soft[2] = {-1, -1} ;

	double depthSum = 0 ;
	const int nextSpliceInd = spliceInd + 1 ;
	int start, end ;
	
	if ( spliceInd < 0 )
		start = 0 ;
	else
		start = splices[ spliceInd ].pos + 1 ;
	if ( nextSpliceInd >= scnt )
		end = INF ;
	else
		end = splices[ nextSpliceInd ].pos - 1 ;
	
	if ( end < start && !isFuture )	
		return true ;

	static char newChrom[50] = "" ;
	static int newPos, newDepth ;
	double depth ;
	int prevPos = start - 1 ;
	int cutStart = -1, cutEnd ;
	int binFill = 0, binPos ;
	double binCover = 0 ;

	int seStart, seEnd ; // Note: this start and end means the start and end of an exon.
	bool first = true ;
	bool covered ;
	bool cut = false ;
	
	int fret = 1 ;

	struct _islandCover currIslandCover ;
	currIslandCover.used = 0 ;

	// Record the information of the islands in this interval.

	double currCover = 0 ;
	double leftIntervalCover = 0, rightIntervalCover = 0 ;
	double leftIntervalAdjust = 0, rightIntervalAdjust = 0 ; // The total coverage around the splice junction, which is likely to be the mismatch.

	bool differentStrands ;

	tmpSeCnt = 0 ;

	// Some constant
	// When to ignore a gap
	int CUT_MERGE ; 
	if ( spliceInd >= 0 && nextSpliceInd < scnt )
	 	CUT_MERGE = ( splices[spliceInd].strand == splices[nextSpliceInd].strand && 
				splices[spliceInd].type == 1 && splices[nextSpliceInd].type == 0 && 
				splices[nextSpliceInd].pos - splices[spliceInd].pos < LONG_EXON ) ? 50 : 15 ;
	else
		CUT_MERGE = 15 ;

	// When to regard this depth is actually empty.
	double CUT_THRESHOLD ;
	if ( spliceInd >= 0 && nextSpliceInd < scnt && 
		splices[spliceInd].type == 0 && splices[nextSpliceInd].type == 1 )
	{
		CUT_THRESHOLD = CUT_THRESHOLD = log( 5.0 ) / log( LOG_BASE ) ;
	}		
	else
		CUT_THRESHOLD = 1 ;

	//if ( !isFuture )
	//	printf( "IsInExon: %s %d %d\n", chrom, start, end ) ;
	/*if ( end < start + CUT_MERGE )
	{
		if ( spliceInd >= 0 )
		{
			spliceInfo[ spliceInd ].depthSum = 0 ;
			spliceInfo[ spliceInd ].possibleIR = false ;
		}
		return true ;
	}*/

	if ( !isFuture && nextSpliceInd < scnt )
	{
		GetFutureIslands( chrom, nextSpliceInd ) ;
	}

	if ( spliceInd >= 0 )
		differentStrands = ( splices[ spliceInd ].strand != splices[ nextSpliceInd  ].strand && 
				splices[ spliceInd ].type != splices[ nextSpliceInd ].type ) ; 
	else
		differentStrands = false ;

	soft[0] = soft[1] = seStart = seEnd = -1 ;	
	//seCnt = 0 ;

	if ( isFuture )
		first = false ;

	while ( 1 )
	{
		if ( fret == EOF )
			break ;

		if ( newChrom[0] && !strcmp( newChrom, chrom ) && nextSpliceInd == 0 && !isFuture)
		{
			strcpy( inChrom, newChrom ) ;
			pos = newPos ;
			idepth = newDepth ;
			newChrom[0] = '\0' ;
		}
		else
		{
			if ( newChrom[0] && strcmp( newChrom, chrom ) && !isFuture )
			{
				// This case happens when the remaining splice sites of the chromosome is filter in depth
				// or a chromosome with no splice sites.
				int chromId = GetChromIdFromName( chrom ) ;
				int newChromId = GetChromIdFromName( newChrom ) ;
				if (  chromId < newChromId )
				{
					// the remaining splice sites in this chromsome is filtered in depth file.
					strcpy( inChrom, chrom ) ;
					idepth = 100000 ;
					pos = end + 1 ;
				}
				else //if (chromId > newChromId)
				{
					// there is a chromosome with no splice sites
					char prevChrom[50] ;
					strcpy( prevChrom, newChrom ) ;
					newChrom[0] = '\0' ;
					
					while ( 1 )
					{
						fret = fscanf( fpDepth, "%s %d %d", inChrom, &pos, &idepth ) ;
						if ( fret == EOF )
							break ;
						if ( !strcmp( prevChrom, inChrom ) )
							continue ;
						else
						{
							int inChromId = GetChromIdFromName( inChrom ) ;
							if ( chromId < inChromId )
							{
								strcpy( inChrom, chrom ) ;
								idepth = 100000 ;
								pos = end + 1 ;
								
								break ;
							}
							else if ( chromId == inChromId )
								break ;
							else // chromId > inChromId
							{
								strcpy( prevChrom, inChrom ) ;
								continue ;	
							}
						}
					}
				}
			}
			else
			{
				fret = fscanf( fpDepth, "%s %d %d", inChrom, &pos, &idepth ) ;
			}
		}

		if ( fret == EOF )
		{
			pos = end + 1 ;
			idepth = 100000 ;
		} 
		/*if ( fret != EOF && first && strcmp( inChrom, chrom ) )
			continue ;
		if ( fret != EOF && ( !first || isFuture ) && strcmp( inChrom, chrom ) )
		{
			if ( !isFuture )
			{
				strcpy( newChrom, inChrom ) ;
				newPos = pos ;
				newDepth = idepth ;
			}

			idepth = 100000 ;
			pos = end + 1 ;
			//break ;
		}*/

		if ( fret != EOF && strcmp( inChrom, chrom ) )
		{
			if ( !isFuture )
			{
				// Let the logic at the beginning in this while(1) to decide whether to stop or keep reading.
				// It might happen that the first depth position we get is from next chromosome, 
				// e.g. the first chromosome is filtered in depth file.
				strcpy( newChrom, inChrom ) ;
				newPos = pos ;
				newDepth = idepth ;
				continue ; 
			}
			else
			{
				idepth = 100000 ;
				pos = end + 1 ;
			}
		}

		//if ( !isFuture )
		//	printf( "FindRegions: %s %d %d\n", inChrom, pos, idepth ) ;

		if ( idepth )
			depth = log( (double)idepth ) / log( LOG_BASE ) ; 
		else
			depth = 0 ;

		if ( pos < start )
			continue ;
		else if ( pos > end && ( nextSpliceInd < scnt || depth >= CUT_THRESHOLD ) ) // Now pos should be as far as the splice site.
		{
			if ( prevPos < pos - 1 && cutStart == -1 )
				cutStart = prevPos + 1 ;

			if ( cutStart != -1 )
			{
				// There may be a gap before the end.
				cutEnd = pos - 1 ;
				if ( cutEnd - cutStart + 1 > CUT_MERGE && ( differentStrands || NoPairReadsSpan( chrom, cutStart, cutEnd, start, end, isFuture, CUT_THRESHOLD ) ) )
				{
					if ( soft[0] == -1 )
					{
						soft[0] = cutStart - 1 ;
						if ( cutStart != start )
						{
							leftIntervalCover = ( currCover - leftIntervalAdjust ) / 
								(double)( cutStart - start - 3 * MERGE_DISTANCE + 1 ) ;	
							if ( spliceInd >= 0 )
								spliceInfo[ spliceInd ].softDepthSum[1] = currCover - leftIntervalAdjust ;
						}
					}

					
					soft[1] = cutEnd + 1 ;	
					if ( seStart != -1 )
					{
						if ( !isFuture )
						{
							tmpSingleExon[ tmpSeCnt ][0] = seStart ;
							tmpSingleExon[ tmpSeCnt ][1] = cutStart - 1 ;
							++tmpSeCnt ;
						}

						for ( i = seStart ; i + BIN_SIZE < cutStart - 1 && 
								currIslandCover.used < MAX_CURRENT_ISLANDS ; i += BIN_SIZE )
						{
							binCover = 0 ;
							for ( j = i ; j < i + BIN_SIZE ; ++j )
								binCover += islandDepth[j - seStart] ;				
							currIslandCover.avgCover[ currIslandCover.used ] = binCover / (double)BIN_SIZE ; 
							//currIslandCover.avgCover[ currIslandCover.used ] = GetMedian( &islandDepth[i - seStart], BIN_SIZE ) ;
							currIslandCover.pos[ currIslandCover.used ] = i ;
							++currIslandCover.used ;
						}

						if ( cutStart - 1 - i >= BIN_SIZE / 2 && currIslandCover.used < MAX_CURRENT_ISLANDS )
						{
							binCover = 0 ;
							for ( j = i ; j < cutStart - 1 ; ++j )
								binCover += islandDepth[j - seStart] ;

							currIslandCover.avgCover[ currIslandCover.used ] = binCover / (double)BIN_SIZE ; 
							//currIslandCover.avgCover[ currIslandCover.used ] = GetMedian( &islandDepth[i - seStart], cutStart - 1 - i ) ;
							currIslandCover.pos[ currIslandCover.used ] = i ;
							++currIslandCover.used ;
						}
					}
					cut = true ;
					currCover = 0 ;
					binCover = 0 ;
				}
			}
			rightIntervalCover= ( currCover - rightIntervalAdjust )/ (double)( end - soft[1] + 1 - 3 * MERGE_DISTANCE + 1 ) ;
			//if ( soft[1] == -1 && rightIntervalCover > 0 && nextSpliceInd == 0 )
			//	soft[1] = start ;	
			
			if ( nextSpliceInd < scnt )
				spliceInfo[ nextSpliceInd ].softDepthSum[0] = currCover - rightIntervalAdjust ;
			break ;
		}

		if ( pos > prevPos + 1 )
		{	
			//return false ;
			if ( cutStart == -1 )
				cutStart = prevPos + 1 ;
			cutEnd = pos - 1 ;
			//cut = true ;
		}	
		if ( pos < start - 1 + 3 * MERGE_DISTANCE )
			leftIntervalAdjust += depth ;
		if ( pos > end + 1 - 3 * MERGE_DISTANCE )
			rightIntervalAdjust += depth ;

		depthSum += depth ;
		if ( spliceInd >= 0 && nextSpliceInd < scnt && 
			splices[ spliceInd ].type == 0 && splices[ nextSpliceInd ].type == 1 && 
			( pos < start - 1 + 3 * MERGE_DISTANCE || pos > end + 1 - 3 * MERGE_DISTANCE ) )
			depthSum -= depth ;
		
		if ( depth < CUT_THRESHOLD && ( evidenceExons == NULL  //|| ( spliceInd >= 0 && nextSpliceInd < scnt && splices[spliceInd].type == 0 && splices[nextSpliceInd].type == 1 )
						|| idepth == 0 || !IsPosInEvidenceExon( pos, evidenceExons, eviExonCnt ) ) )
		{
			if ( cutStart == -1 )
			{
				//if ( !isFuture ) printf( "cutstart %d %d\n", pos, idepth ) ;
				cutStart = pos ;
			}
			cutEnd = pos ;
			//cut = true ;
		}
		else
		{
			if ( first && nextSpliceInd == 0 )
			{
				seStart = pos ;
				soft[1] = pos ;
				//printf( "%d\n", cutStart ) ;
			}
			else if ( cutStart != -1 && cutEnd - cutStart + 1 > CUT_MERGE && ( differentStrands || NoPairReadsSpan( chrom, cutStart, cutEnd, start, end, isFuture, CUT_THRESHOLD ) ) ) 
			{
				//printf( "### %d %d %d\n", cutStart, cutEnd, differentStrands ) ;
				// The two soft boundary should be in the same exon.
				if ( soft[0] == -1 )
				{
					soft[0] = cutStart - 1 ;
					if ( cutStart != start )
					{
						leftIntervalCover = ( currCover - leftIntervalAdjust ) / (double)( cutStart - start ) ;	
						if ( spliceInd >= 0 )
							spliceInfo[ spliceInd ].softDepthSum[1] = currCover - leftIntervalAdjust ;
					}
				}
				soft[1] = cutEnd + 1 ;

				if ( seStart != -1 )
				{
					if ( !isFuture )
					{
						tmpSingleExon[ tmpSeCnt ][0] = seStart ;
						tmpSingleExon[ tmpSeCnt ][1] = cutStart - 1 ;
						++tmpSeCnt ;
					}				
					for ( i = seStart ; i + BIN_SIZE < cutStart - 1 && currIslandCover.used < MAX_CURRENT_ISLANDS ; i += BIN_SIZE )
					{
						binCover = 0 ;
						for ( j = i ; j < i + BIN_SIZE ; ++j )
							binCover += islandDepth[j - seStart] ;				
						currIslandCover.avgCover[ currIslandCover.used ] = binCover / (double)BIN_SIZE ; 
						//currIslandCover.avgCover[ currIslandCover.used ] = GetMedian( &islandDepth[i - seStart], BIN_SIZE ) ;
						currIslandCover.pos[ currIslandCover.used ] = i ;
						++currIslandCover.used ;
					}

					if ( cutStart - 1 - i >= BIN_SIZE / 2 && currIslandCover.used < MAX_CURRENT_ISLANDS )
					{
						binCover = 0 ;
						for ( j = i ; j < cutStart - 1 ; ++j )
							binCover += islandDepth[j - seStart] ;
						
						currIslandCover.avgCover[ currIslandCover.used ] = binCover / (double)BIN_SIZE ; 
						//currIslandCover.avgCover[ currIslandCover.used ] = GetMedian( &islandDepth[i - seStart], cutStart - 1 - i ) ;
						currIslandCover.pos[ currIslandCover.used ] = i ;
						++currIslandCover.used ;
					}

//printf( "island %d %d %d %d\n", start, currIslandCover.used, pos - binFill + 1, binFill ) ;
					binFill = 0 ;
				}
				seStart = cutEnd + 1 ;
				cut = true ;
				currCover = 0 ;
			}
			cutStart = -1 ;
			if ( depth >= CUT_THRESHOLD )
				currCover += depth ;

		}
		//if ( seStart != -1 && pos - seStart < MAX_LENGTH )
		//{
		
		tmp = start ;
		if ( seStart != -1 )
			tmp = seStart ;	
		if ( ( seStart != -1 && pos - tmp < MAX_LENGTH ) || 
			( !cut && spliceInd >= 0 && splices[spliceInd].type == 1 && splices[nextSpliceInd ].type == 0 && 
			  pos - tmp < MAX_LENGTH ) ) // Either it is in a island or a intron-retention portion.
		{
			for ( i = ( prevPos + 1 - tmp > 0 ? prevPos + 1 : tmp ) ; i < pos ; ++i )
				islandDepth[i - tmp] = 0 ;
			//printf( "%d %llf\n", pos - tmp, depth ) ;
			islandDepth[pos - tmp] = depth ;
		}
		//}

		first = false ;
		prevPos = pos ;
	}

	//soft[1] = cutEnd + 1 ;
	//if ( soft[1] == -1 && nextSpliceInd == 0 && )

	if ( soft[0] < start - 1 + 3 * MERGE_DISTANCE )
		soft[0] = -1 ;
	if ( soft[1] > end + 1 - 3 * MERGE_DISTANCE  )
		soft[1] = -1 ;
	
	if ( spliceInd >= 0 && !isFuture )
	{
		spliceInfo[ spliceInd ].soft[1] = soft[0] ;
		spliceInfo[ spliceInd ].avgSoftDepth[1] = leftIntervalCover ;
	}
	if ( nextSpliceInd < scnt && !isFuture )
	{
		spliceInfo[ nextSpliceInd ].soft[0] = soft[1] ;
		spliceInfo[ nextSpliceInd ].avgSoftDepth[0] = rightIntervalCover ;
	}	
	//if ( soft[0] != -1 )
	//	printf( "left: %d %lf\n", start, leftIntervalCover ) ;
	//if ( soft[1] != -1 )
	//	printf( "right: %d %lf\n", end, rightIntervalCover ) ;	
	/*if ( cut == true )
	{
	printf( "### %d %d %d %d\n", start, end, soft[0], soft[1] ) ;
	}*/

	if ( spliceInd >= 0 && !isFuture )
	{
		//if ( !isFuture && start == 155166029 )
		//	printf( "hi %d %d\n", start, end ) ;
		spliceInfo[ spliceInd ].depthSum = depthSum ;
		spliceInfo[ spliceInd ].possibleIR = 0 ;
		if ( !cut && spliceInd >= 0 && nextSpliceInd < scnt && 
				splices[ spliceInd ].type == 0 && splices[nextSpliceInd].type == 1 && end - start + 1 > 6 * MERGE_DISTANCE - 2 )
		{
			spliceInfo[ spliceInd ].possibleIR = 1 ;
			cut = true ;
		}
		else if ( !cut && spliceInd >= 0 && nextSpliceInd < scnt && 
				splices[ spliceInd ].type == 0 && splices[nextSpliceInd].type == 1 && end - start + 1 <= 6 * MERGE_DISTANCE - 2 )
		{
			spliceInfo[ spliceInd ].depthSum = currCover ; //( leftIntervalAdjust + rightIntervalAdjust ) / 2 ;
			spliceInfo[ spliceInd ].possibleIR = 1 ;
			cut = true ;
		}

		if ( spliceInfo[ spliceInd ].possibleIR == 1 )
		{
			// Test whether the coverage of the intron is much higher than the splice junction
			if ( 1 ) //splices[ spliceInd ].otherInd == nextSpliceInd )
			{
				int len ;
				double avg ;
				if ( splices[ nextSpliceInd ].pos - splices[ spliceInd ].pos - 1 > 6 * MERGE_DISTANCE - 2 ) 
					len = splices[ nextSpliceInd ].pos - splices[ spliceInd ].pos - 1 - ( 6 * MERGE_DISTANCE - 2 ) ;
				else
					len = splices[ nextSpliceInd ].pos - splices[ spliceInd ].pos - 1 ;
				avg = spliceInfo[ spliceInd ].depthSum / (double)len ;
				if ( pow( LOG_BASE, avg ) > 4 * ( splices[ spliceInd ].support + splices[ nextSpliceInd ].support ) )
					//&& splices[ spliceInd ].otherInd == nextSpliceInd )
				{
					//printf( "%d: %d %d: %d %d: %d\n", end - start - 1, splices[ spliceInd ].support, splices[ nextSpliceInd ].support,
					//	splices[ spliceInd ].pos, splices[ nextSpliceInd ].pos, spliceInd ) ;
					//cut = false ;
					//spliceInfo[ spliceInd ].possibleIR = 0 ;
					spliceInfo[ spliceInd ].possibleIR = 2 ;
				}
			}

			// Test whether this is actaully from a 3', 5' UTR.
			if ( 0 ) //end - start + 1 >=  3 * BIN_SIZE && end - start + 1 < MAX_LENGTH )
			{
				bool flag = true ;
				// Test increasing from left to right _-^
				double sum, firstSum, prevSum ;
				sum = 0 ;
				for ( i = end ; i > end - BIN_SIZE ; --i )
				{
					sum += islandDepth[i - start] ;
				}
				firstSum = sum ;
				prevSum = sum ;
				for ( ; i >= start + BIN_SIZE - 1 ; i -= BIN_SIZE )
				{
					sum = 0 ;
					for ( j = i ; j > i - BIN_SIZE ; --j )
						sum += islandDepth[j - start] ;
					if ( sum >= prevSum )
					{
						flag = false ;
						break ;
					}
					prevSum = sum ;
				}

				if ( flag && sum / BIN_SIZE <= firstSum / BIN_SIZE - 8 )
				{
					cut = true ;
					spliceInfo[ spliceInd ].possibleIR = 0 ;
					sum = 0 ;
					for ( i = ( start + end ) / 2 ; i < end - 3 * MERGE_DISTANCE ; ++i )
					{
						sum += islandDepth[i - start] ;
					}
					//printf( "Success: %d\n", splices[ nextSpliceInd ].pos ) ;		
					spliceInfo[ nextSpliceInd ].soft[0] = ( start + end ) / 2 ;
					spliceInfo[ nextSpliceInd ].avgSoftDepth[0] = sum / ( i - (start + end) / 2) ;
				}

				// Test decreasing from right to left ^-_
				flag = true ;
				sum = 0 ;
				for ( i = start ; i < start + BIN_SIZE ; ++i )
				{
					sum += islandDepth[i - start] ;
				}
				firstSum = sum ;
				prevSum = sum ;
				for ( ; i <= end - BIN_SIZE + 1 ; i += BIN_SIZE )
				{
					sum = 0 ;
					for ( j = i ; j < i + BIN_SIZE ; ++j )
						sum += islandDepth[j - start] ;
					if ( sum >= prevSum )
					{
						flag = false ;
						break ;
					}
					prevSum = sum ;
				}

				if ( flag && sum / BIN_SIZE <= firstSum / BIN_SIZE - 8 )
				{
					cut = true ;
					spliceInfo[ spliceInd ].possibleIR = 0 ;
					sum = 0 ;
					for ( i = start + 3 * MERGE_DISTANCE ; i <= ( start + end ) / 2 ; ++i )
					{
						sum += islandDepth[i - start] ;
					}
					
					spliceInfo[spliceInd].soft[1] = ( start + end ) / 2 ;
					spliceInfo[spliceInd].avgSoftDepth[1] = sum / ( i - start - 3 * MERGE_DISTANCE) ;
				}

			}
		}
	}

	if ( !isFuture )
	{
		// Get the single exons
		for ( i = 0 ; i < tmpSeCnt && !isFuture ; ++i )
		{
			if ( IsSingleExon( tmpSingleExon[i][0], tmpSingleExon[i][1], prevIslandCover, currIslandCover, futureIslandCover ) )
			{
				//printf( "## %d %d %d %d\n", start, end, tmpSingleExon[i][0], tmpSingleExon[i][1] ) ;
				singleExon[ seCnt ][0] = tmpSingleExon[i][0] ;
				singleExon[ seCnt ][1] = tmpSingleExon[i][1] ;
				++seCnt ;
			}
		}


		// Update the prevIsland by currIsland
		/*for ( i = 0 ; i < prevIslandCover.used ; ++i )
			printf( "#%d ", prevIslandCover.pos[i] ) ;
		printf( "\n" ) ;
		for ( i = 0 ; i < currIslandCover.used ; ++i )
			printf( "%d ", currIslandCover.pos[i] ) ;
		printf( "\n" ) ;*/
		if ( currIslandCover.used >= BIN_USE )
		{
			for ( k = 0, i = currIslandCover.used - BIN_USE ; i < currIslandCover.used ; ++i, ++k )
			{
				prevIslandCover.avgCover[k] = currIslandCover.avgCover[i] ;
				prevIslandCover.pos[k] = currIslandCover.pos[i] ;
			}
			prevIslandCover.used = k ;
		}
		else
		{
			int tmp = BIN_USE - currIslandCover.used ; 
			for ( k = 0, i = ( prevIslandCover.used - tmp < 0 ? 0 : prevIslandCover.used - tmp ) ; i < prevIslandCover.used ; ++i, ++k )
			{
				prevIslandCover.avgCover[k] = prevIslandCover.avgCover[i] ;
				prevIslandCover.pos[k] = prevIslandCover.pos[i] ;
			}
			for ( i = 0 ; i < currIslandCover.used ; ++i, ++k )
			{
				prevIslandCover.avgCover[k] = currIslandCover.avgCover[i] ;
				prevIslandCover.pos[k] = currIslandCover.pos[i] ;
			}
			prevIslandCover.used = k ;
		}
	}
	else
	{
		// Update the futureIsland by currIsland
		for ( i = 0, k = futureIslandCover.used ; i < currIslandCover.used ; ++i, ++k )
		{
			if ( k >= BIN_USE )
				break ;
			futureIslandCover.avgCover[k] = currIslandCover.avgCover[i] ;
			futureIslandCover.pos[k] = currIslandCover.pos[i] ;
		}
		futureIslandCover.used = k ;
	}
	
	/*if ( start == 31246180 )
	{
		printf( "%d: %d %d\n", cut, start, end ) ;
		exit( 1 ) ;
	}*/

#if ALL_CLEAN

	if ( spliceInd >= 0 && splices[spliceInd].type == 0 && splices[nextSpliceInd].type == 1 )
	{
		soft[0] = soft[1] = -1 ;
		cut = true ;	
	}
#endif
	return !cut ;
}

// Decide the regions from splice startInd to splice endInd
double DetermineRegions( char *chrom, int startInd, int endInd ) 
{
	int i, j, k ;
	int idepth ;
	double depth ;
	double avgDepth, stdevDepth ;
	char inChrom[50] ;
	int pos ;
	double depthSum = 0 ;
	int intronLen = 0 ;
	const double alpha = IR_ALPHA ;
	double noiseThreshold = 0 ;
	double threshold ; // TODO: different threshold for intron retention and 3'5' start site
	int exonLen = 0 ;
	double exonDepthSum = 0, exonAvgDepth ;
	const double exonAlpha = 1.5 ;

	int start, end ;
	start = splices[ startInd ].pos ;
	end = splices[ endInd ].pos ;

	// Determine the possibleIR=2 case
	for ( i = startInd ; i < endInd ; ++i )
	{
		if ( spliceInfo[i].possibleIR != 2 )
			continue ;
		int betterBound = 0 ;
		int len ;
		if ( splices[i + 1].pos - splices[i].pos - 1 > 6 * MERGE_DISTANCE - 2 ) 
			len = splices[i + 1].pos - splices[i].pos - 1 - ( 6 * MERGE_DISTANCE - 2 ) ;
		else
			len = splices[i + 1].pos - splices[i].pos - 1 ;
		for ( j = i - 1 ; j >= startInd ; --j )
		{
			if ( spliceInfo[j].regionId != spliceInfo[i].regionId )
				break ;
			if ( splices[j].strand == splices[i].strand && splices[j].type == splices[i].type && 
				splices[j].support > pow( LOG_BASE, spliceInfo[i].depthSum / len )  )
				//splices[j].support > 4 * ( splices[i].support + splices[i + 1].support ) )
			{
				++betterBound ;
				break ;
			}
		}
		for ( j = i + 2 ; j <= endInd ; ++j )
		{
			if ( spliceInfo[j].regionId != spliceInfo[i + 1].regionId )
				break ;
			if ( splices[j].strand == splices[i + 1].strand && splices[j].type == splices[i + 1].type && 
				splices[j].support > pow( LOG_BASE, spliceInfo[i].depthSum / len )  )
			{
				++betterBound ;
				break ;
			}
		}
		if ( betterBound > 0 )
		{
			// We can split the two regions
			//printf( "possibleIR2->1 %d %d\n", splices[i].pos, splices[i + 1].pos ) ;
			spliceInfo[i].possibleIR = 1 ;
			continue ;
		}
		// merge the two regions
		//printf( "possibleIR2->0 %d %d\n", splices[i].pos, splices[i + 1].pos ) ;
		spliceInfo[i].possibleIR = 0 ;
		k = spliceInfo[i + 1].regionId ;
		for ( j = i + 1 ; j <= endInd ; ++j )
		{
			if ( spliceInfo[j].regionId != k )
				break ;
			spliceInfo[j].regionId = spliceInfo[i].regionId ;
		}
	}

	// Get the total coverage of the possible introns.
	for ( i = startInd ; i < endInd ; ++i )
	{
		//if ( spliceInfo[i].soft[1] != -1 || spliceInfo[i + 1].soft[0] != -1 
		//	|| ( splices[i].type == 0 && splices[i + 1].type == 1 ) )
		if ( spliceInfo[i].regionId != spliceInfo[i + 1].regionId )
		{
			depthSum += spliceInfo[i].depthSum ;
			if ( splices[i + 1].pos != splices[i].pos )
			{
				intronLen += ( splices[i + 1].pos - splices[i].pos - 1 ) ;
			}

			if ( i != startInd && splices[i].type == 0 && spliceInfo[i].soft[0] != -1 )
			{
				depthSum -= spliceInfo[i].softDepthSum[0] ;
				intronLen -= ( splices[i].pos - spliceInfo[i].soft[0] ) ;

				//exonDepthSum += spliceInfo[i].softDepthSum[0] ;
				//exonLen += ( splices[i].pos - spliceInfo[i].soft[0] ) ;
			}

			if ( i != endInd && splices[i].type == 1 && spliceInfo[i].soft[1] != -1 )
			{
				depthSum -= spliceInfo[i].softDepthSum[1] ;
				intronLen -= ( spliceInfo[i].soft[1] - splices[i].pos ) ;
				
				//exonDepthSum += spliceInfo[i].softDepthSum[1] ;
				//exonLen += ( splices[i].pos - spliceInfo[i].soft[1] ) ;
			}
		}
		else
		{
			exonDepthSum += spliceInfo[i].depthSum ;
			exonLen += ( splices[i + 1].pos - splices[i].pos ) ;
			//printf( "exon: %d %lf %lf\n", splices[i].pos, spliceInfo[i].depthSum / ( splices[i + 1].pos - splices[i].pos ), spliceInfo[i].depthSum ) ;
		}
	}
	avgDepth = depthSum / intronLen ;
	stdevDepth = sqrt( avgDepth ) ; // Suppose a Poisson model

	if ( exonLen != 0 )
		exonAvgDepth = exonDepthSum / exonLen ;
	else
		exonAvgDepth = 10 * avgDepth ; 

	threshold = avgDepth + alpha * stdevDepth ;
	if ( noiseCnt >= scnt )
	{
		double noiseAvg = noiseDepth / noiseCnt ;
		double noiseStdev = noiseSqDepth / noiseCnt - noiseAvg * noiseAvg ;
		noiseStdev = sqrt( noiseStdev ) ;
		
		noiseThreshold = noiseAvg + 2 * noiseStdev ;
	}
	else //if ( noiseCnt >= scnt / 2 )
	{
		double noiseAvg = noiseDepth / noiseCnt ;
		double noiseStdev = noiseSqDepth / noiseCnt - noiseAvg * noiseAvg ;
		noiseStdev = sqrt( noiseStdev ) ;

		noiseThreshold = noiseAvg + 2 * noiseStdev ;
		noiseThreshold = sqrt( noiseThreshold ) ;
		if ( noiseCnt == 0 )
			noiseThreshold = 0 ;
	}

	if ( noiseThreshold > threshold )
		threshold = noiseThreshold ;
	if ( startInd == endInd )
	{
		if ( splices[i].type == 0 )
			exonAvgDepth = spliceInfo[i].avgSoftDepth[0] ;
		else
			exonAvgDepth = spliceInfo[i].avgSoftDepth[1] ;
		threshold = exonAvgDepth ;
	}
	
	// Find the regions requiring forced soft boundary
	double partialExonDepthSum = 0 ;
	int partialExonLen = 0 ;
	double prevExonDepthSum = -1 ;
	int prevExonLen = -1 ;
	for ( i = startInd ; i <= endInd ; )
	{
		for ( j = i + 1 ; j <= endInd ; ++j )
			if ( spliceInfo[i].regionId != spliceInfo[j].regionId )
				break ;
		int startCnt = 0, endCnt = 0 ; // start-start of an exon
		for ( k = i ; k < j ; ++k )
			if ( splices[k].type == 0 && spliceInfo[k].regionId != spliceInfo[ splices[k].otherInd ].regionId )
				endCnt += splices[k].support ;
			else if ( splices[k].type == 1 && spliceInfo[k].regionId != spliceInfo[ splices[k].otherInd ].regionId )
				startCnt += splices[k].support ;
		if ( startCnt > 100 * endCnt )
		{
			if ( spliceInfo[j - 1].soft[1] != -1 )
				spliceInfo[j - 1].avgSoftDepth[1] = MAX_LENGTH ;
		}
		if ( endCnt > 100 * startCnt )
		{
			if ( spliceInfo[i].soft[0] != -1 )
				spliceInfo[i].avgSoftDepth[0] = MAX_LENGTH ;
		}

		if ( startCnt > 0 && 8 * startCnt < pow( LOG_BASE, exonAvgDepth ) )
		{
			if ( spliceInfo[i].soft[0] != -1 )
				spliceInfo[i].avgSoftDepth[0] = MAX_LENGTH ;
		}

		if ( endCnt > 0 && 8 * endCnt < pow( LOG_BASE, exonAvgDepth ) )
		{
			if ( spliceInfo[j - 1].soft[1] != -1 )
				spliceInfo[j - 1].avgSoftDepth[1] = MAX_LENGTH ;
		}


		// Check whether the left portion of the gene and the right portion of the gene has 
		// siginifcan different coverage.
		/*int currentExonDepthSum = 0 ;
		int currentExonLen = 0 ;
		for ( k = i ; k < j - 1 ; ++k )
		{
			currentExonDepthSum += spliceInfo[k].depthSum ;
			currentExonLen += ( splices[k + 1].pos - splices[k].pos ) ;
		}

		if ( exonLen != 0 && partialExonLen != 0 &&
			i > startInd && j - 1 < endInd && prevExonLen > 0 )
		{
			if ( 4 * partialExonDepthSum / partialExonLen < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) && 
					4 * prevExonDepthSum / prevExonLen < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) && 
					4 * spliceInfo[ startInd ].avgSoftDepth[0] < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) ) 
			{
				spliceInfo[i].avgSoftDepth[0] = MAX_LENGTH ;
				if ( spliceInfo[i].soft[0] == -1 )
					spliceInfo[i].soft[0] = splices[i].pos - 1 ;
			}
		}

		partialExonDepthSum += currentExonDepthSum ;
		partialExonLen += currentExonLen ;
		prevExonDepthSum = currentExonDepthSum ;
		prevExonLen = currentExonLen ;*/

		i = j ;
	}

	/*partialExonDepthSum = 0 ;
	partialExonLen = 0 ;
	prevExonDepthSum = -1 ;
	prevExonLen = -1 ;
	for ( i = endInd ; i >= startInd ; )
	{
		for ( j = i - 1 ; j >= startInd ; --j )
			if ( spliceInfo[i].regionId != spliceInfo[j].regionId )
				break ;
		int currentExonDepthSum = 0 ;
		int currentExonLen = 0 ;
		for ( k = j + 1 ; k < i ; ++k )
		{
			currentExonDepthSum += spliceInfo[k].depthSum ;
			currentExonLen += ( splices[k + 1].pos - splices[k].pos ) ;
		}

		if ( exonLen != 0 && partialExonLen != 0 &&
			i < endInd && j + 1 > startInd && prevExonLen > 0 )
		{
			if ( 4 * partialExonDepthSum / partialExonLen < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) && 
					4 * prevExonDepthSum / prevExonLen < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) &&
					4 * spliceInfo[ endInd ].avgSoftDepth[1] < ( exonDepthSum - partialExonDepthSum ) / ( exonLen - partialExonLen ) ) 
			{
				spliceInfo[i].avgSoftDepth[1] = MAX_LENGTH ;
				if ( spliceInfo[i].soft[1] == -1 )
					spliceInfo[i].soft[1] = splices[i].pos + 1 ;
			}
		}

		partialExonDepthSum += currentExonDepthSum ;
		partialExonLen += currentExonLen ;
		prevExonDepthSum = currentExonDepthSum ;
		prevExonLen = currentExonLen ;

		i = j ;
	}*/
	/*for ( i = startInd ; i < endInd ; ++i )
	{
		if ( spliceInfo[i].regionId == spliceInfo[i + 1].regionId )
			continue ;
		int cnt = splices[i].support + splices[i + 1].support ;
		if ( cnt > 10 )
			continue ;
		if ( 4 * cnt < pow( LOG_BASE, exonAvgDepth ) )
		{
			if ( spliceInfo[i].soft[1] != -1 )
				spliceInfo[i].avgSoftDepth[1] = MAX_LENGTH ;
			if ( spliceInfo[i + 1].soft[0] != -1 )
				spliceInfo[i + 1].avgSoftDepth[0] = MAX_LENGTH ;
		}
	}*/


	// Test the soft ends and IR.
	for ( i = startInd ; i <= endInd ; ++i )
	{
		if ( splices[i].type == 0 && spliceInfo[i].soft[1] != -1 )
		{
			//printf( "%d: %lf | %lf %lf\n", splices[i].pos, spliceInfo[i].avgSoftDepth[1], threshold, exonAvgDepth  ) ;
			if ( spliceInfo[i].avgSoftDepth[1] <= threshold && 
				( spliceInfo[i].avgSoftDepth[1] < exonAlpha * exonAvgDepth || spliceInfo[i].avgSoftDepth[1] < noiseThreshold ) )
				spliceInfo[i].soft[1] = -1 ;
			else
			{
				depthSum -= spliceInfo[i].softDepthSum[1] ;
				intronLen -= ( spliceInfo[i].soft[1] - splices[i].pos ) ;
			}
		}

		if ( splices[i].type == 1 && spliceInfo[i].soft[0] != -1 )
		{
			//printf( "%d: %lf | %lf %lf\n", splices[i].pos, spliceInfo[i].avgSoftDepth[0], threshold, exonAvgDepth  ) ;
			if ( spliceInfo[i].avgSoftDepth[0] <= threshold && 
				( spliceInfo[i].avgSoftDepth[0] < exonAlpha * exonAvgDepth || spliceInfo[i].avgSoftDepth[0] < noiseThreshold ) )
			{
				spliceInfo[i].soft[0] = -1 ;
			}
			else
			{
				depthSum -= spliceInfo[i].softDepthSum[0] ;
				intronLen -= ( splices[i].pos - spliceInfo[i].soft[0] ) ;
			}
		}
	}
	
	// TODO: revise the threshold
	avgDepth = depthSum / intronLen ;
	stdevDepth = sqrt( avgDepth ) ; // Suppose a Poisson model

	threshold = avgDepth + alpha * stdevDepth ;
	if ( noiseThreshold > threshold )
		threshold = noiseThreshold ;

	for ( i = startInd ; i <= endInd ; ++i )
	{
		spliceInfo[i].threshold = threshold ;
		if ( threshold > exonAlpha * exonAvgDepth && exonAvgDepth > noiseThreshold )
		{
			threshold = exonAlpha * exonAvgDepth ;
		}
	}
	//printf( "==%d %d==\n", splices[startInd].pos, splices[endInd].pos ) ;
	//printf( "%lf %d\n", exonAvgDepth, exonLen ) ;
	for ( i = startInd ; i <= endInd ; ++i )
	{
		if ( spliceInfo[i].possibleIR != 0 )
		{
			
			bool forceIR = false ;
			int len ;
			if ( splices[i + 1].pos - splices[i].pos - 1 > 6 * MERGE_DISTANCE - 2 ) 
				len = splices[i + 1].pos - splices[i].pos - 1 - ( 6 * MERGE_DISTANCE - 2 ) ;
			else
				len = splices[i + 1].pos - splices[i].pos - 1 ;
			
			if ( len == 0 )
				forceIR = true ;

			/*printf( "possible IR: %d - %d: intronAvg: %lf threshold: %lf intronLen:%d exonAvg: %lf noise: %lf %d(%d)\n", splices[i].pos, splices[i + 1].pos, 
				spliceInfo[i].depthSum / len, threshold, intronLen, exonAvgDepth,
				noiseThreshold, noiseCnt, scnt ) ;*/

			double otherDepth ;
			

			if ( i - 2 >= startInd && i + 1 <= endInd && spliceInfo[i].regionId == spliceInfo[i - 1].regionId && 
					spliceInfo[i].regionId != spliceInfo[i + 1].regionId &&
					spliceInfo[i - 1].regionId != spliceInfo[i - 2].regionId &&
					splices[i - 1].support > splices[i].support * 100 &&
					splices[i].type == 0 && splices[i - 1].type == 1 )
			{
				forceIR = true ;
			}
			
			if ( i >= startInd && i + 3 <= endInd && spliceInfo[i + 1].regionId == spliceInfo[i + 2].regionId && 
					spliceInfo[i + 1].regionId != spliceInfo[i].regionId &&
					spliceInfo[i + 2].regionId != spliceInfo[i + 3].regionId &&
					splices[i + 2].support > splices[i + 1].support * 100 &&
					splices[i + 1].type == 1 && splices[i + 2].type == 0 )
			{
				forceIR = true ;
			}
		
			/*if ( ( splices[i].support + splices[i + 1].support ) * 4 < spliceInfo[i].depthSum / len && 
				spliceInfo[i].depthSum / len > noiseThreshold )
			{
				forceIR = true ;
			}*/

			if ( forceIR || spliceInfo[i].depthSum / len > threshold || 
				( spliceInfo[i].depthSum / len >= exonAlpha * exonAvgDepth && spliceInfo[i].depthSum / len > noiseThreshold ) )
			{
				//printf( "IR: %d - %d: %lf %lf\n", splices[i].pos, splices[i + 1].pos, 
				//	spliceInfo[i].depthSum / len, threshold ) ;
				// Merge two groups
				k = spliceInfo[i + 1].regionId ;
				for ( j = i + 1 ; j <= endInd ; ++j )
				{
					if ( spliceInfo[j].regionId != k )
						break ;
					spliceInfo[j].regionId = spliceInfo[i].regionId ;
				}
			}
		}
	}
}


void FindRegions_Init( char *prefix, char *otherJunctionFile )
{
	char buffer[1024] ;
	char name[10] ;
	int size ; 
	int i ;
	//FILE *fpSizes = fopen( sizes, "r" ) ;
	sprintf( buffer, "%s.depth", prefix) ;
	fpDepth = fopen( buffer, "r" ) ;
	fpDepthForGroup = fopen( buffer, "r" ) ;	

	if ( otherJunctionFile == NULL )
		sprintf( buffer, "%s.splice", prefix) ;
	else
		strcpy( buffer, otherJunctionFile ) ;

	fpSplice = NULL ;
	fpSplice = fopen( buffer, "r" ) ;
	if ( fpSplice == NULL )
	{
		printf( "Could not find splice junction file %s.\n", buffer ) ;
		exit( 1 ) ;
	}	

	sprintf( buffer, "%s.sam", prefix) ;
	fpReads = OpenReadFile( prefix ) ;
	fpReadsForGroup = OpenReadFile( prefix ) ;
	
	splices = ( struct _splice *)malloc( sizeof( struct _splice ) * 1000003 ) ; 
	tmpSplices = ( struct _splice *)malloc( sizeof( struct _splice ) * 1000003 ) ; 
	scnt = 0 ;
	seCnt = 0 ;
	
	// Get the information about each chromosome from the depth file.
	//printf( "### %s %s\n", depth, tophat ) ;
	chromCnt = INF ;
	curChromInd = -1 ;
	/*chromCnt = 0 ;
	while ( fscanf( fpSizes, "%s %d", name, &size ) != EOF )
	{
		strcpy( chroms[ chromCnt ].name, name ) ;
		chroms[ chromCnt ].end = size ;
		++chromCnt ;
	}
	fclose( fpSizes ) ;*/
	
	prevIslandCover.used = 0 ;
	spliceInfo = NULL ;

	evidenceExons = NULL ;
	eviExonCnt = 0 ;
}

// Groups are just preliminary regions
int GetSpliceGroups( char *chrom, struct _pair groups[] )
{
	if ( scnt == 0 )
		return 0 ;
	int i, j, k ;
	FILE *fpDepthBackup = fpDepth ;
	struct _readFile fpReadsBackup = fpReads ;
	int gid = 0 ;
	fpDepth = fpDepthForGroup ;
	fpReads = fpReadsForGroup ;

	// Reuse splices.soft as group id.
	spliceInfo[0].regionId = gid ;
	groups[0].a = 0 ;

	prevIslandCover.used = 0 ;
	futureIslandCover.used = 0 ;

	IsInExon( chrom, -1, false ) ;	
	for ( i = 1 ; i < scnt ; ++i )
	{
		if ( IsInExon( chrom, i - 1, false ) )
		{
			spliceInfo[i].regionId = gid ;
		}
		else
		{
			++gid ;
			groups[gid - 1].b = i - 1 ;
			groups[gid].a = i ;
			spliceInfo[i].regionId = gid ;
		}
	}
	IsInExon( chrom, scnt - 1, false ) ;
	groups[gid].b = i - 1 ; 


	fpDepth = fpDepthBackup ;
	fpReads = fpReadsBackup ;
	return gid + 1 ;
}

int RegroupSplices( struct _pair groups[] )
{
	if ( scnt == 0 )
		return 0 ;
	int i ;
	int gid = spliceInfo[0].regionId ;
	groups[ spliceInfo[0].regionId ].a = 0 ;

	for ( i = 1 ; i < scnt ; ++i )
	{
		if ( spliceInfo[i].regionId != spliceInfo[i - 1].regionId )
		{
			groups[ spliceInfo[i - 1].regionId ].b = i - 1 ;
			groups[ spliceInfo[i].regionId ].a = i ;
		}
		if ( spliceInfo[i].regionId > gid )
			gid = spliceInfo[i].regionId ;
	}
	groups[ spliceInfo[i - 1].regionId ].b = i - 1 ;
	return gid + 1 ;
}

// Give gene id to the groups.
// diffPM: support for plus splice sites - support minus splice sites. Short for difference between plus and minus.
void SearchGene( int tag, int gene[], int &geneCnt, struct _pair groups[], bool gvisited[], int diffPM )
{
	if ( gvisited[tag] )
		return ;
	gvisited[tag] = true ;

	int i, j, k ;
	int minusCnt = 0, plusCnt = 0 ;
	/*for ( i = groups[tag].a ; i <= groups[tag].b ; ++i )
	{
		if ( splices[i].strand )
			++plusCnt ;
		else
			++minusCnt ;
	}
	if ( plusCnt && minusCnt )
	{
		// A mixture region.
		if ( geneCnt == 0 )
		{
			gene[ geneCnt ] = tag ;
			++geneCnt ;
		}
		return ;
	}*/
	
	for ( i = groups[tag].a ; i <= groups[tag].b ; ++i )
	{
		if ( splices[i].strand == 1 )
			diffPM += splices[i].support ;
		else if ( splices[i].strand == 0 )
			diffPM -= splices[i].support ;
	}

	gene[ geneCnt ] = tag ;
	++geneCnt ;
	for ( i = groups[tag].a ; i <= groups[tag].b ; ++i )
	{
		// noisy splice sites.
		if ( spliceInfo[i].possibleIR != 2 && 
			( ( diffPM < 0 && splices[i].strand == 1 ) ||
			( diffPM > 0 && splices[i].strand == 0 ) ) )
			continue ;
		SearchGene( spliceInfo[ splices[i].otherInd ].regionId, gene, geneCnt, groups, gvisited, diffPM ) ;
	}

	// Force go to next group if there is a possible IR.
	/*if ( spliceInfo[ groups[tag].b ].possibleIR )
	{
		SearchGene( spliceInfo[ groups[tag].b + 1 ].regionId, gene, geneCnt, groups, gvisited, diffPM ) ;
	}*/
}

void TestCleanDataSet( char *chrom, struct _pair *groups, int gcnt, int *groupGeneId )
{
	fpos_t fpos ;
	int i, j, k, tmp ;
	int noiseIntron = 0 ;
	static bool tested = false ;
	static int defaultBinUse = BIN_USE, defaultBinRange = BIN_RANGE ;	
	if ( tested )
		return ;
	char buffer[10] ;
	int pos, cover, prevpos ;
	//int *depth = (int *)malloc( sizeof( int ) * 5000000 ) ;
	tested = true ;
	fgetpos( fpDepth, &fpos ) ;
	rewind( fpDepth ) ;
	j = 0 ;
	for ( i = 0 ; i < gcnt - 1 ; ++i )
	{
		//if ( groupGeneId[i] != groupGeneId[i + 1] )
		//	continue ;
		k = splices[ groups[i + 1].a ].pos - splices[groups[i].b].pos - 1 ; 
		j = 0 ;
		tmp = 0 ;
		prevpos = splices[groups[i].b].pos ;
		bool started = false ;
		//memset( depth, 0, sizeof( int ) * k ) ;
		while ( 1 )
		{
			if ( fscanf( fpDepth, "%s %d %d", buffer, &pos, &cover ) == EOF )
				break ;

			if ( pos <= splices[ groups[i].b ].pos )
				continue ;	
			if ( pos >= splices[ groups[i + 1].a ].pos )
				break ;

			if ( pos > prevpos + 1 || cover == 0 )
			{
				if ( !started )
					started = true ;
				else
				{
					j += tmp ;
					tmp = 0 ;
				}
			}


			if ( started && cover >= 2 )
				++tmp ;
			prevpos = pos ;	
			//depth[ pos - splices[groups[i + 1].a].pos -1 ] = cover ;
		}
		

		if ( (double)j / (double)k >= 0.1 )
			++noiseIntron ;
		//printf( "%d %d\n", j, k ) ;
		//for ( j = 0 ; j < k ; ++k )
	}		

	//printf( "%d %d\n", noiseIntron, gcnt ) ;	
	//exit( 1 ) ;
	if ( noiseIntron <= 0.02 * gcnt )
	{
		BIN_USE = 0 ;
		BIN_RANGE = 0 ;
	}
	else
	{
		BIN_USE = defaultBinUse ;
		BIN_RANGE = defaultBinRange ;
	}
	fsetpos( fpDepth, &fpos ) ;
}


int CompSpliceSupport( const void *p1, const void *p2 )
{
	int a = *(int *)p1 ;
	int b = *(int *)p2 ;
	int oa = splices[a].otherInd ;
	int ob = splices[b].otherInd ;
	if ( splices[a].support != splices[b].support )
		return splices[a].support - splices[b].support ;
	/*else if ( splices[a].uniqSupport != splices[b].uniqSupport )
		return splices[a].uniqSupport - splices[b].uniqSupport ;
	else if ( splices[a].secSupport != splices[b].secSupport )
		return splices[a].secSupport - splices[b].secSupport ;*/
	else
	{
		int lena = splices[a].pos < splices[oa].pos ? ( splices[oa].pos - splices[a].pos ) : 
			( splices[a].pos - splices[oa].pos ) ;
		int lenb = splices[b].pos < splices[ob].pos ? ( splices[ob].pos - splices[b].pos ) :
			( splices[b].pos - splices[ob].pos ) ;
		if ( lenb - lena >= -25 && lenb - lena <= 25 )
		{
			if ( splices[a].uniqSupport != splices[b].uniqSupport )
				return splices[a].uniqSupport - splices[b].uniqSupport ;
			else if ( splices[a].secSupport != splices[b].secSupport )
				return splices[a].secSupport - splices[b].secSupport ;
		}
		else
			return lenb - lena ;
	}
}

void PassSoftBoundary( int tag )
{
	int i, k ;		
	for ( k = tag + 1 ; k < scnt && spliceInfo[k].regionId == spliceInfo[tag].regionId ; ++k )
	{
		if ( splices[k].support == -1 ) //|| splices[k].pos == splices[tag].pos )
			continue ;
		spliceInfo[k].soft[0] = spliceInfo[tag].soft[0] ;
		if ( splices[k].type == 0 && spliceInfo[k].soft[0] == -1 )
		{
			spliceInfo[k].soft[0] = ( splices[k].pos - READS_LENGTH < splices[tag].pos ) ? 
				( splices[k].pos - READS_LENGTH ) : splices[tag].pos ;
			if ( spliceInfo[k].soft[0] <= 0 )
				spliceInfo[k].soft[0] = 1 ;
		}
		else if ( splices[k].type == 1 && spliceInfo[k].soft[0] == -1 && splices[k].pos != splices[tag].pos )
		{
			// Test whether the portion between tag and k belongs to exons.
			int len = 0 ;
			double sum = 0 ;
			for ( i = tag ; i < k ; ++i )
			{
				if ( splices[i].pos == splices[i + 1].pos )
					continue ;

				if ( splices[i + 1].pos - splices[i].pos - 1 > 6 * MERGE_DISTANCE - 2 ) 
					len += splices[i + 1].pos - splices[i].pos - 1 - ( 6 * MERGE_DISTANCE - 2 ) ;
				else
					len += splices[i + 1].pos - splices[i].pos - 1 ;
				sum += spliceInfo[i].depthSum ; 
			}

			if ( sum / len > spliceInfo[tag].threshold )
			{
				//printf( "%d: %lf %lf\n", splices[k].pos, sum / len, spliceInfo[tag].threshold ) ;
				spliceInfo[k].soft[0] = ( splices[k].pos - READS_LENGTH < splices[tag].pos ) ? 
					( splices[k].pos - READS_LENGTH ) : splices[tag].pos ;
				if ( spliceInfo[k].soft[0] <= 0 )
					spliceInfo[k].soft[0] = 1 ;
			}
		}
		break ;
	}

	for ( k = tag - 1 ; k >= 0 && spliceInfo[k].regionId == spliceInfo[tag].regionId ; --k )
	{
		if ( splices[k].support == -1 ) //|| splices[k].pos == splices[tag].pos )
			continue ;
		spliceInfo[k].soft[1] = spliceInfo[tag].soft[1] ;
		if ( splices[k].type == 1 && spliceInfo[k].soft[1] == -1 )
			spliceInfo[k].soft[1] = ( splices[k].pos + READS_LENGTH > splices[tag].pos ) ?
				( splices[k].pos + READS_LENGTH ) : splices[tag].pos ;
		else if ( splices[k].type == 0 && spliceInfo[k].soft[1] == -1 && splices[k].pos != splices[tag].pos )
		{
			int len = 0 ;
			double sum = 0 ;
			for ( i = k ; i < tag ; ++i )
			{
				if ( splices[i].pos == splices[i + 1].pos )
					continue ;

				if ( splices[i + 1].pos - splices[i].pos - 1 > 6 * MERGE_DISTANCE - 2 ) 
					len += splices[i + 1].pos - splices[i].pos - 1 - ( 6 * MERGE_DISTANCE - 2 ) ;
				else
					len += splices[i + 1].pos - splices[i].pos - 1 ;
				sum += spliceInfo[i].depthSum ; 
			}

			if ( sum / len > spliceInfo[tag].threshold )
			{
				spliceInfo[k].soft[1] = ( splices[k].pos + READS_LENGTH > splices[tag].pos ) ? 
					( splices[k].pos + READS_LENGTH ) : splices[tag].pos ;
			}
		}

		break ;
	}
}


void PreprocessSplices( char *chrom, struct _evidence *evidences, int eviCnt )
{
	int i, j, k, l ;
	int gcnt ;
	int geneCnt ;
	bool *gvisited ; // Is a group visited
	struct _pair *groups = ( struct _pair * )malloc( sizeof( struct _pair ) * scnt ) ;
	int *gene ;
	int *groupGeneId ; // gene id for each group
	int *spliceIds ; // the indices of splices
	int geneId = 0 ;
	int eviTag = 0 ;
	int solvedSpliceId ;

	// Firstly, remove the negative splice junctions, or splices spanning too much or the edit distance is too much.
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].otherInd < i || 
			( splices[i].support > 0 && ( splices[ splices[i].otherInd ].pos - splices[i].pos <= 100000 || 
			 splices[i].otherInd - i <= 400 )  
			 && splices[i].support * 2 > splices[i].uniqEditDistance + splices[i].secEditDistance ) ) 
			continue ;
		if ( evidences != NULL )
		{	
			while ( eviTag < eviCnt && evidences[eviTag].exons[ evidences[eviTag].ecnt - 1].end < splices[i].pos )
				++eviTag ;

			for ( k = eviTag ; k < eviCnt && evidences[k].exons[0].start < splices[i].pos ; ++k )
			{
				if ( evidences[k].exons[0].strand != splices[i].strand )
					continue ;
				int e ;
				for ( e = 0 ; e < evidences[k].ecnt - 1 ; ++e )
					if ( WithinEvidenceMargin( splices[i].pos, evidences[k].exons[e].end ) &&
							WithinEvidenceMargin( splices[ splices[i].otherInd ].pos, evidences[k].exons[e + 1].start ) )
						break ;
				if ( e < evidences[k].ecnt - 1 )
				{
					/*if ( splices[i].pos == 14631400 )
					{
printf( "%d %d %d %d\n", splices[i].strand, evidences[k].exons[0].strand, evidences[k].exons[0].start, evidences[k].exons[0].end ) ;
exit( 1 ) ;
					}*/	
					k = -1 ;
					break ;
				}
			}	

			if ( k == -1 && splices[i].support < 0 && 
				splices[i].uniqSupport >= 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) ) 
			{
				splices[i].support = 1 ;
				splices[ splices[i].otherInd ].support = 1 ;
			}
		}
		else if ( splices[i].support > 0 )
		{
			//printf( "revmoed %d %d: %d\n", splices[i].pos, splices[ splices[i].otherInd ].pos, splices[i].support ) ;
			splices[i].support = -1 ;
			splices[ splices[i].otherInd ].support = -1 ;
		}
	}
	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support <= 0 )
		{
			continue ;
		}
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;
	eviTag = 0 ;
	
	// If a site is associate with too many introns, we need to filter some of those, since no matter what
	// the algorithm will have difficulty of picking the right path at this region.
	for ( i = 0 ; i < scnt ;  )
	{
		// sites in the range of [i,j) have the same coordinate.
		for ( j = i + 1 ; j < scnt ; ++j )	
			if ( splices[j].pos != splices[i].pos )
				break ;
		if ( j - i < 10 )
		{
			i = j ;
			continue ;	
		}

		int maxUniqSupport = -1 ;
		for ( k = i ; k < j ; ++k )		
			if ( splices[k].uniqSupport > maxUniqSupport )
				maxUniqSupport = splices[k].uniqSupport ;

		for ( k = i ; k < j ; ++k )
		{
			if ( splices[k].support <= 0.01 * maxUniqSupport )
			{
				splices[k].support = -1 ;
				splices[ splices[k].otherInd ].support = -1 ;
			}
			int span = splices[ splices[k].otherInd].pos - splices[k].pos ;
			if ( span < 0 )
				span = -span ;
			if ( span >= 100000 )
			{
				splices[k].support = -1 ;
				splices[ splices[k].otherInd ].support = -1 ;
			}
		}

		i = j ;
	}
	
	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support <= 0 )
		{
			continue ;
		}
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;


	// Then we remove a splice that are in a bunch of opposite strands splice sites, and support by no more than 2 reads.
	for ( i = 0 ; i < scnt ; ++i )
	{
		int l, r, lcnt, rcnt ;
		bool remove = false ;
		if ( splices[i].otherInd < i || splices[i].support > 2 || splices[i].otherInd - i == 1 )
			continue ;

		lcnt = 0 ;
		rcnt = 0 ;
		for ( l = i - 1 ; l >= 0 && lcnt < 3 ; --l )
		{
			if ( splices[l].support <= 2 )
				continue ;
			if ( splices[l].strand == splices[i].strand )
				break ;
			//if ( splices[l].pos != splices[l + 1].pos )
			++lcnt ;
		}
		for ( r = i + 1 ; r < scnt && rcnt < 3 ; ++r )
		{
			if ( splices[r].support <= 2 )
				continue ;
			if ( splices[r].strand == splices[i].strand )
				break ;
			//if ( splices[r].pos != splices[r- 1].pos )
			++rcnt ;
		}
		
		if ( lcnt >= 3 && rcnt >= 3 )
			remove = true ;

		lcnt = 0 ;
		rcnt = 0 ;
		
		for ( l = splices[i].otherInd - 1 ; l >= 0 && lcnt < 3 ; --l )
		{
			if ( splices[l].support <= 2 )
				continue ;
			if ( splices[l].strand == splices[i].strand )
				break ;
			//if ( splices[l].pos != splices[l + 1].pos )
			++lcnt ;
		}
		for ( r = splices[i].otherInd + 1 ; r < scnt && rcnt < 3 ; ++r )
		{
			if ( splices[r].support <= 2 )
				continue ;
			if ( splices[r].strand == splices[i].strand )
				break ;
			//if ( splices[r].pos != splices[r- 1].pos )
			++rcnt ;
		}
		if ( lcnt >= 3 && rcnt >= 3 )
			remove = true ;

		if ( evidences != NULL && remove )
		{	
			while ( eviTag < eviCnt && evidences[eviTag].exons[ evidences[eviTag].ecnt - 1].end < splices[i].pos )
				++eviTag ;

			for ( k = eviTag ; k < eviCnt && evidences[k].exons[0].start < splices[i].pos ; ++k )
			{
				if ( evidences[k].exons[0].strand != splices[i].strand )
					continue ;
				int e ;
				for ( e = 0 ; e < evidences[k].ecnt - 1 ; ++e )
					if ( WithinEvidenceMargin( splices[i].pos, evidences[k].exons[e].end ) &&
							WithinEvidenceMargin( splices[ splices[i].otherInd ].pos, evidences[k].exons[e + 1].start ) )
						break ;
				if ( e < evidences[k].ecnt - 1 )
				{
					/*if ( splices[i].pos == 14631400 )
					{
printf( "%d %d %d %d\n", splices[i].strand, evidences[k].exons[0].strand, evidences[k].exons[0].start, evidences[k].exons[0].end ) ;
exit( 1 ) ;
					}*/	
					k = -1 ;
					break ;
				}
			}	

			if ( k == -1  )
			{
				continue ;
			}
		}

		if ( remove )
		{
			//printf( "Removed: %d %d\n", splices[i].pos, splices[ splices[i].otherInd ].pos ) ;
			splices[i].support = -1 ;
			splices[ splices[i].otherInd ].support = -1 ;
		}
	}

	// Remove the splice sites with no strand information if the support is weak
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support <= 0 || splices[i].strand >= 0 )
			continue ;

		if ( splices[i].uniqSupport < 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) )
		{
			splices[i].support = -1 ;
			splices[ splices[i].otherInd ].support = -1 ;
		}
	}

	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support <= 0 )
			continue ;
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;
	eviTag = 0 ;

	// We get the rough regions and the single exon transcripts
	gcnt = GetSpliceGroups( chrom, groups ) ;

	/*for ( i = 0 ; i < scnt ; ++i )
	{
		printf( "regionId %d: %d %d possible IR: %d soft: %d %d support: %d other: %d\n", 
			i, splices[i].pos, spliceInfo[i].regionId, spliceInfo[i].possibleIR, 
			spliceInfo[i].soft[0], spliceInfo[i].soft[1],
			splices[i].support,
			splices[splices[i].otherInd].pos ) ;
	}
	fflush( stdout ) ;*/
	//printf( "noiseCnt: %d spliceCnt: %d\n", noiseCnt, scnt ) ;

	gvisited = ( bool * )malloc( sizeof( bool ) * gcnt ) ;
	gene = ( int * )malloc( sizeof( int ) * gcnt ) ;
	groupGeneId = ( int * )malloc( sizeof( int ) * gcnt ) ;
	spliceIds = ( int *)malloc( sizeof( int ) * scnt ) ;
	memset( gvisited, false, sizeof( bool ) * gcnt ) ;
	memset( groupGeneId, -1, sizeof( int ) * gcnt ) ;
	solvedSpliceId = -1 ;
	/*for ( i = 0 ; i < scnt ; ++i )
	{
		printf( "0 soft[1] %d %d %d\n", i, splices[i].pos, spliceInfo[i].soft[1] ) ;
	}*/

	for ( i = 0 ; i < gcnt ; ++i )
	{
		int mins, maxs ;
		if ( gvisited[i] || groups[i].a < solvedSpliceId )
			continue ;
		geneCnt = 0 ;
		SearchGene( i, gene, geneCnt, groups, gvisited, 0 ) ;

		// Fill in the strand information for non-canonical splice sites.
		int prevMainStrand = 1 ;
		for ( j = 0 ; j < geneCnt ; ++j )
		{
			int strandSupport = 0 ;
			for ( k = groups[ gene[j] ].a ; k <= groups[ gene[j] ].b ; ++k )
			{ 
				if ( splices[k].strand >= 0 )
					strandSupport += ( 2 * splices[k].strand - 1 ) * splices[k].support ;
			}
			groupGeneId[ gene[j] ] = strandSupport ; // reuse this array to hold temporary strand information
		}

		for ( j = 0 ; j < geneCnt ; ++j )
		{
			int groupStrand ;
			if ( groupGeneId[ gene[j] ] != 0 )
				groupStrand = groupGeneId[ gene[j] ] > 0 ? 1 : 0 ;
			else
			{
				for ( k = 1 ; k < geneCnt ; ++k )
				{
					if ( j - k < 0 && j + k >= geneCnt )
						break ;
					if ( j + k < geneCnt )
					{
						if ( groupGeneId[ gene[j + k ] ] != 0 )
							groupStrand = groupGeneId[ gene[j + k] ] > 0 ? 1 : 0 ;
					}
					if ( j - k >= 0 )
					{
						if ( groupGeneId[ gene[j - k ] ] != 0 )
							groupStrand = groupGeneId[ gene[j - k] ] > 0 ? 1 : 0 ;
					}
				}
			}
			// If all the introns in the gene are non-canonical, we just assign them to plus strand.
			if ( groupStrand == -1 )
				groupStrand = 1 ;
			for ( k = groups[ gene[j] ].a ; k <= groups[ gene[j] ].b ; ++k )
			{
				if ( splices[k].strand == -1 )
					splices[k].strand = groupStrand ;	
			}
		}

		mins = scnt + 10 ;
		maxs = -1 ;
		for ( j = 0 ; j < geneCnt ; ++j )
		{
			if ( groups[ gene[j] ].a < mins )
				mins = groups[ gene[j] ].a ;
			if ( groups[ gene[j] ].b > maxs )
				maxs = groups[ gene[j] ].b ;
		}
		//printf( "%d-%d: 15076 %d\n", mins, maxs, spliceInfo[15076].possibleIR ) ;
		DetermineRegions( chrom, mins, maxs ) ;
		solvedSpliceId = maxs ;
	}
	memset( groupGeneId, -1, sizeof( int ) * gcnt ) ;

	/*for ( i = 0 ; i < scnt ; ++i )
	{
		printf( "- regionId: %d %d\n", splices[i].pos, spliceInfo[i].regionId ) ;
	}*/

	gcnt = RegroupSplices( groups ) ;
	memset( gvisited, false, sizeof( bool ) * gcnt ) ;

	/**======================================================================
	Remove the splice junctions possibly from gene merge.
	=========================================================================*/
	geneId = 0 ;
	// Mark the splice junctions looks like duplicated
	int *tmpSpliceSupport ;
	tmpSpliceSupport = (int *)malloc( sizeof( int ) * scnt ) ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		//if ( splices[i].otherInd < i )
		//	continue ;
		tmpSpliceSupport[i] = splices[i].support ;
		if ( splices[i].uniqSupport < 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) ) 
		{
			//splices[i].support = 0 ; //splices[i].uniqSupport / 5 ;
			splices[i].support = splices[i].uniqSupport - splices[i].secSupport * 0.01 ;
			if ( splices[i].support < 0 )
				splices[i].support = 0 ;
		}
	}

	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( gvisited[ spliceInfo[i].regionId ] )
			continue ; 
		int cnt = 0 ;
		geneCnt = 0 ;
		SearchGene( spliceInfo[i].regionId, gene, geneCnt, groups, gvisited, 0 ) ;
		
		// We found a gene.
		int maxSupport = -1 ;
		for ( j = 0 ; j < geneCnt ; ++j )
		{
			for ( k = groups[ gene[j] ].a ; k <= groups[ gene[j] ].b ; ++k )
			{
				if ( splices[k].support > maxSupport )
					maxSupport = splices[k].support ;
				spliceIds[ cnt ] = k ;
				spliceInfo[k].geneId = geneId ;
				//printf( "%d %d: %d %d\n", k, splices[k].pos, splices[ splices[k].otherInd ].pos, splices[k].support ) ;
				++cnt ;
			}
			groupGeneId[ gene[j] ] = geneId ;
		}
		++geneId ;
		//printf( "maxSupport %d: %d\n", splices[i].pos, maxSupport ) ;
	
		if ( evidences != NULL )	
			while ( eviTag < eviCnt && evidences[eviTag].exons[ evidences[eviTag].ecnt - 1].end < splices[spliceIds[0]].pos )
				++eviTag ;
		int plusCnt, minusCnt ; // The support for the splices on different strands.
		int mainStrand ;

		qsort( spliceIds, cnt, sizeof( int ), CompSpliceSupport ) ;

		int maxAvgSupport = 0 ;
		if ( cnt > 100 )
		{
			maxAvgSupport = ( splices[ spliceIds[ cnt - 1] ].support + 
				splices[ spliceIds[ cnt - 3 ] ].support + splices[ spliceIds[ cnt - 5 ] ].support ) / 3 ;
		}
		int threshold ;
		
		if ( maxAvgSupport < 100000 )
		{
			threshold = maxSupport * 0.01 + 1 ;
			if ( threshold > 10 )
				threshold = 10 ;
			if ( threshold == 0 )
				threshold = 1 ;
		}
		/*else if ( maxAvgSupport < 200000 )
		{
			threshold = 0.001 * maxSupport ;
		}*/
		else
		{
			threshold = 0.01 * maxSupport ;
		}
		//printf( "%d %d %d\n", maxAvgSupport, maxSupport, threshold ) ;
		plusCnt = 0 ;
		minusCnt = 0 ;
		for ( l = 0 ; l < cnt ; ++l )
		{
			if ( splices[ spliceIds[l] ].strand == 0 )
				minusCnt += splices[ spliceIds[l] ].support ;
			else
				plusCnt += splices[ spliceIds[l] ].support ;
		}
		if ( minusCnt > plusCnt )
			mainStrand = 0 ;
		else
			mainStrand = 1 ;

		//printf( "searching\n" ) ;	
		for ( l = 0 ; l < cnt ; ++l )
		{
			int groupIdJ, groupIdO ;
			bool leftOK = true, rightOK = true;
			j = spliceIds[l] ;
			if ( splices[j].support >= threshold )
				break ;
			int other = splices[j].otherInd ;
			
			groupIdJ = spliceInfo[j].regionId ;
			groupIdO = spliceInfo[ other ].regionId ;

			if ( splices[j].support == -1 || groupGeneId[ groupIdJ ] != groupGeneId[ groupIdO ] )
				continue ;
			
			for ( k = groups[ groupIdJ ].a ; k <= groups[ groupIdJ ].b ; ++k )
			{
				if ( k != j && splices[k].support != -1 && splices[k].type == splices[j].type &&
					splices[k].strand == splices[j].strand && 
					spliceInfo[k].regionId != spliceInfo[ splices[k].otherInd ].regionId )
					break ;
			}
			if ( k > groups[ groupIdJ ].b )
				leftOK = false ;

			for ( k = groups[ groupIdO ].a ; k <= groups[ groupIdO ].b ; ++k )
			{
				if ( k != other && splices[k].support != -1 && splices[k].type == splices[other].type &&
					splices[k].strand == splices[other].strand &&
					spliceInfo[k].regionId != spliceInfo[ splices[k].otherInd ].regionId )
					break ;
			}
			if ( k > groups[ groupIdO ].b )
				rightOK = false ;

			if ( splices[j].strand != mainStrand )
			{
				leftOK = true ;
				rightOK = true ;
			}

			if ( !leftOK && !rightOK )
				continue ;

			if ( evidences != NULL )
			{		
				for ( k = eviTag ; k < eviCnt && evidences[k].exons[0].start < splices[j].pos ; ++k )
				{
					if ( evidences[k].exons[0].strand != splices[j].strand )
						continue ;
					int e ;
					for ( e = 0 ; e < evidences[k].ecnt - 1 ; ++e )
						if ( WithinEvidenceMargin( splices[j].pos, evidences[k].exons[e].end ) &&
							WithinEvidenceMargin( splices[other].pos, evidences[k].exons[e + 1].start ) )
							break ;
					if ( e < evidences[k].ecnt - 1 &&
							splices[i].uniqSupport >= 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) ) 
					{
						k = -1 ;
						break ;
					}
				}	
				if ( k == -1 )
					continue ;
			}
			//printf( "Cleaned Splice: %d-%d %d %d\n", splices[j].pos, splices[other].pos, splices[j].support,threshold ) ;

			// Pass the soft boundaries to other existing splices
			// The logic is: soft[0] is attached to the first splice of a region, 
			// and soft[1] attached to the last splice of a region.
			PassSoftBoundary( j ) ;
			splices[j].support = -1 ;

			PassSoftBoundary( other ) ;
			splices[other].support = -1 ;
		}
	}
	/*for ( i = 0 ; i < scnt ; ++i )
		printf( "%d: %d region: %d soft: %d %d other:(%d %d) support: %d\n", i, splices[i].pos, spliceInfo[i].regionId, spliceInfo[i].soft[0], spliceInfo[i].soft[1], splices[i].otherInd, splices[i].otherPos, splices[i].support ) ;
	fflush( stdout ) ;*/
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support >= 0 )
			splices[i].support = tmpSpliceSupport[i] ;
	}
	free( tmpSpliceSupport ) ;
	//TestCleanDataSet( chrom, groups, gcnt, groupGeneId ) ;
	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support == -1 )
		{
			//printf( "- %d\n", splices[i].pos ) ;
			continue ;
		}
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;
	// Remove the noise splice in a region
	/*for ( i = 0 ; i < scnt ;  )
	  {
	  for ( j = i + 1 ; j < scnt ; ++j )
	  if ( spliceInfo[j].regionId != spliceInfo[i].regionId )
	  break ;


	  i = j ;
	}*/
	/**=================================================================
	Remove splice junction in a false region
	================================================================== */
	for ( i = 0 ; i < scnt ; ++i )
	{
		//if ( splices[i].support != 1 )
		//	continue ;
			
		if ( ( ( i == 0 || spliceInfo[i].regionId != spliceInfo[i - 1].regionId ) && splices[i].type == 0 && spliceInfo[i].soft[0] == -1 ) ||
			( ( i == scnt -1 || spliceInfo[i].regionId != spliceInfo[i + 1].regionId ) && splices[i].type == 1 && spliceInfo[i].soft[1] == -1 ) )
		/*if ( ( i == 0 || spliceInfo[i].regionId != spliceInfo[i - 1].regionId ) && ( i == scnt - 1 || spliceInfo[i].regionId != spliceInfo[i + 1].regionId ) ) 
		{
			if ( ( splices[i].type == 0 && spliceInfo[i].soft[0] == -1 ) ||
				( splices[i].type == 1 && spliceInfo[i].soft[1] == -1 ) )
			{
				int other = splices[i].otherInd ;
					
				for ( k = other + 1 ; k < scnt && spliceInfo[k].regionId == spliceInfo[other].regionId ; ++k )
				{
					if ( splices[k].support == -1 )
						continue ;
					spliceInfo[k].soft[0] = spliceInfo[other].soft[0] ;
					if ( splices[k].type == 0 && spliceInfo[k].soft[0] == -1 )
						spliceInfo[k].soft[0] = ( splices[k].pos - READS_LENGTH < splices[other].pos ) ? 
							( splices[k].pos - READS_LENGTH ) : splices[other].pos ;
					break ;
				}

				for ( k = other - 1 ; k >= 0 && spliceInfo[k].regionId == spliceInfo[other].regionId ; --k )
				{
					if ( splices[k].support == -1 )
						continue ;
					spliceInfo[k].soft[1] = spliceInfo[other].soft[1] ;
					if ( splices[k].type == 1 && spliceInfo[k].soft[1] == -1 )
						spliceInfo[k].soft[1] = ( splices[k].pos + READS_LENGTH > splices[other].pos ) ?
							( splices[k].pos + READS_LENGTH ) : splices[other].pos ;
					break ;
				}
				splices[i].support = -1 ;
				splices[other].support = -1 ;
			}

		}*/
		{
			//printf( "PreprocessSplices: %d-%d\n", splices[i].pos, splices[ splices[i].otherInd ].pos ) ;
			PassSoftBoundary( i ) ;
			splices[i].support = -1 ;

			int other = splices[i].otherInd ;
			if ( splices[other].support != -1 )
			{
				PassSoftBoundary( other ) ;	
				splices[other].support = -1 ;
			}
		} // end of if this splce junction should be removed.
	}
	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support == -1 )
		{
			//printf( "- %d\n", splices[i].pos ) ;
			continue ;
		}
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;

	/*=============================================================================
	Remove paralogous( or duplicated ) gene

	If a gene contains too many splices that have no unique support reads, then 
	the gene is suppose to be duplicated. And this junction should not be saved by
	evidences, since the evidences can have those duplicate genes.
	=============================================================================*/
	// Recalculate the gene id for each splice sites
	/*gcnt = RegroupSplices( groups ) ;
	memset( gvisited, false, sizeof( bool ) * gcnt ) ;

	geneId = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( gvisited[ spliceInfo[i].regionId ] )
			continue ;
		geneCnt = 0 ;
		SearchGene( spliceInfo[i].regionId, gene, geneCnt, groups, gvisited, 0 ) ;
		for ( j = 0 ; j < geneCnt ; ++j )
		{
			for ( k = groups[ gene[j] ].a ; k <= groups[ gene[j] ].b ; ++k )
			{	
				spliceInfo[k].geneId = geneId ;
			}
		}
		++geneId ;
	}

	int *junctionsInGene = (int *)malloc( sizeof( int ) * ( geneId + 10 ) );
	int *uniqJunctionsInGene = ( int * )malloc( sizeof( int ) * ( geneId + 10 ) ) ;

	memset( junctionsInGene, 0, sizeof( int ) * ( geneId + 10 ) ) ;
	memset( uniqJunctionsInGene, 0, sizeof( int ) * ( geneId + 10 ) ) ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].otherInd < i )
			continue ;
		++junctionsInGene[ spliceInfo[i].geneId ] ;
		if ( 0 ) //splices[i].uniqSupport >= 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) )
		{
			//printf( "unique junction: %d-%d\n", splices[i].pos, splices[ splices[i].otherInd ].pos ) ;
			++uniqJunctionsInGene[ spliceInfo[i].geneId ] ;
			//splices[i].support = -1 ;
			//splices[ splices[i].otherInd ].support = -1 ;
		}
		if ( splices[i].uniqSupport < 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) )
		{
			PassSoftBoundary( i ) ;
			splices[i].support = -1 ;
			PassSoftBoundary( splices[i].otherInd ) ;
			splices[ splices[i].otherInd ].support = -1 ;
		}
	}

	for ( i = 0 ; i < scnt ; ++i )
	{
		//printf( "geneId:%d pos: %d : %d <-%d\n", 
		//	spliceInfo[i].geneId, splices[i].pos, splices[i].support, splices[ splices[i].otherInd ].pos  ) ;
		if ( uniqJunctionsInGene[ spliceInfo[i].geneId ] < 0.1 * junctionsInGene[ spliceInfo[i].geneId ] && 
			uniqJunctionsInGene[ spliceInfo[i].geneId ] <= 2 )
			splices[i].support = -1 ;
	}*/

	/*for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].otherInd < i || splices[ splices[i].otherInd ].pos - splices[i].pos <= 200000 )
			continue ;
		if ( splices[i].uniqSupport < 0.05 * ( splices[i].uniqSupport + splices[i].secSupport ) )
		{
			PassSoftBoundary( i ) ;
			splices[i].support = -1 ;
			PassSoftBoundary( splices[i].otherInd ) ;
			splices[ splices[i].otherInd ].support = -1 ;
		}
	}

	k = 0 ;
	for ( i = 0 ; i < scnt ; ++i )
	{
		if ( splices[i].support == -1 )
		{
			continue ;
		}
		splices[k] = splices[i] ;
		spliceInfo[k] = spliceInfo[i] ;
		splices[ splices[k].otherInd ].otherInd = k ;
		++k ;
	}
	scnt = k ;*/
	/*free( junctionsInGene ) ;
	free( uniqJunctionsInGene ) ;*/

	free( gvisited ) ;	
	free( groups ) ;
	free( gene ) ;
	free( groupGeneId ) ;
	free( spliceIds ) ;
	/*printf( "scnt=%d\n", scnt ) ;
	for ( i = 0 ; i < scnt ; ++i )
		printf( "%d: %d region: %d soft: %d %d other:(%d %d)\n", i, splices[i].pos, spliceInfo[i].regionId, spliceInfo[i].soft[0], spliceInfo[i].soft[1], splices[i].otherInd, splices[i].otherPos ) ;
	fflush( stdout ) ;*/
}


// Return 0: if hits the EOF
int GetSplices( char *chrom, struct _evidence *evidences, int &eviCnt )
{
	static char inChrom[50] = "-1", inStrand[10] ;
	static int inStart, inEnd, inSupport, inUniqSupport, inSecondarySupport,
		inUniqEditDistance, inSecondaryEditDistance ;
	char curChrom[50] = "-1" ;
	int i, j, k ;
	int a, b ;

	scnt = 0 ;
	++curChromInd ;
	if ( curChromInd >= chromCnt )
		return 0 ;
	if ( curChromInd != 0 )
	{
		strcpy( curChrom, inChrom ) ;

		splices[0].otherInd = 1 ;
		splices[1].otherInd = 0 ;

		splices[0].pos = inStart ;
		splices[0].type = 0 ;
		//splices[0].strand = inStrand[0] == '-' ? 0 : 1 ;
		if ( inStrand[0] == '-' )
			splices[0].strand = 0 ;
		else if ( inStrand[0] == '+' )
			splices[0].strand = 1 ;
		else
			splices[0].strand = -1 ;

		splices[0].support = inSupport ;
		splices[0].otherPos = inEnd ;
		splices[0].uniqSupport = inUniqSupport ;
		splices[0].secSupport = inSecondarySupport ;
		splices[0].uniqEditDistance = inUniqEditDistance ;
		splices[0].secEditDistance = inSecondaryEditDistance ;

		splices[1].pos = inEnd ;
		splices[1].type = 1 ;
		splices[1].strand = splices[0].strand ;
		splices[1].support = inSupport ;
		splices[1].otherPos = inStart ;
		splices[1].uniqSupport = inUniqSupport ;
		splices[1].secSupport = inSecondarySupport ;
		splices[1].uniqEditDistance = inUniqEditDistance ;
		splices[1].secEditDistance = inSecondaryEditDistance ;
		
		scnt = 2 ;
	}
	
	while ( 1 ) 
	{
		if ( fscanf( fpSplice, "%s %d %d %d %s %d %d %d %d", inChrom, &inStart, &inEnd, 
			&inSupport, inStrand, &inUniqSupport, &inSecondarySupport,
			&inUniqEditDistance, &inSecondaryEditDistance ) == EOF )
		//if ( fscanf( fpSplice, "%s %d %d %d %s", inChrom, &inStart, &inEnd, 
		//	&inSupport, inStrand ) == EOF )
		{
			chromCnt = -1 ;
			break ;
		}

		/*if ( !strcmp( inChrom, "chrM" ) || !strcmp( inChrom, "chrm" ) )
		{
			continue ;
		}*/	
		//printf( "%s %s %d\n", curChrom, inChrom, scnt ) ;
		if ( strcmp( inChrom, curChrom ) )
		{ 
			if ( !strcmp( curChrom, "-1" ) )
			{
				strcpy( curChrom, inChrom ) ;
			}
			else
				break ;
		}

		// Determine whether this is should be a real splice junction.
		if ( ( inUniqSupport == 0 && inSecondarySupport < 5 ) ||
			inUniqSupport + inSecondarySupport <= 0 )
		{
			inSupport = -( inUniqSupport + inSecondarySupport ) ;
		}
		
		if ( inEnd - inStart - 1 >= 200000 && inSupport > 0 )
		{
			if ( inUniqSupport < 0.05 * ( inUniqSupport + inSecondarySupport ) )
				inSupport = -( inUniqSupport + inSecondarySupport ) ;
			else if ( inSecondarySupport > 1 )
				inSupport = inUniqSupport + 1 ;
		}

		for ( i = scnt - 1 ; splices[i].pos > inStart && i >= 0 ; --i )
		{
			if ( splices[i].otherInd != -1 )
				++splices[ splices[i].otherInd ].otherInd ; 
			splices[i + 1] = splices[i] ;
		}
		splices[i + 1].pos = inStart ;
		splices[i + 1].type = 0 ;
		//splices[i + 1].strand = inStrand[0] == '-' ? 0 : 1 ;
		if ( inStrand[0] == '-' )
			splices[i + 1].strand = 0 ;
		else if ( inStrand[0] == '+' )
			splices[i + 1].strand = 1 ;
		else
			splices[i + 1].strand = -1 ;
		splices[i + 1].support = inSupport ;
		splices[i + 1].uniqSupport = inUniqSupport ;
		splices[i + 1].secSupport = inSecondarySupport ;
		splices[i + 1].uniqEditDistance = inUniqEditDistance ;
		splices[i + 1].secEditDistance = inSecondaryEditDistance ;
		a = i + 1 ;
		++scnt ;

		for ( i = scnt - 1 ; splices[i].pos > inEnd && i >= 0 ; --i )
		{
			if ( splices[i].otherInd != -1 )
				++splices[ splices[i].otherInd ].otherInd ;
			splices[i + 1] = splices[i] ;
		}
		splices[i + 1].pos = inEnd ;
		splices[i + 1].type = 1 ;
		splices[i + 1].strand = splices[a].strand ; //inStrand[0] == '-' ? 0 : 1 ;
		splices[i + 1].support = inSupport ;
		splices[i + 1].uniqSupport = inUniqSupport ;
		splices[i + 1].secSupport = inSecondarySupport ;
		splices[i + 1].uniqEditDistance = inUniqEditDistance ;
		splices[i + 1].secEditDistance = inSecondaryEditDistance ;
		b = i + 1 ;
		++scnt ;

		splices[a].otherInd = b ;
		splices[b].otherInd = a ;
		splices[a].otherPos = splices[b].pos ;
		splices[b].otherPos = splices[a].pos ;
	}
	if ( spliceInfo != NULL )
		free( spliceInfo ) ;
	spliceInfo = ( struct _spliceInfo *)malloc( sizeof( struct _spliceInfo ) * ( scnt + 10 ) ) ;
	memset( spliceInfo, -1, sizeof( struct _spliceInfo ) * ( scnt + 10 ) ) ;


	strcpy( chrom, curChrom ) ;
	if ( evidences )
	{	
		eviCnt = Evidence_Extract( chrom, evidences ) ;
		/*for ( i = 0 ; i < eviCnt ; ++i )
		{
			char strand[2] = "+" ;
			if ( evidences[i].exons[0].strand == 0 )
				strand[0] = '-' ;
			else if ( evidences[i].exons[0].strand == -1 )
				strand[0] = '.' ;

			printf( "%s\tCLASS\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"adrenal.%d\"; transcript_id \"adrenal.%d\";\n", 
					chrom, evidences[i].exons[0].start, evidences[i].exons[ evidences[i].ecnt - 1].end,
					strand, i, i ) ;

			for ( j = 0 ; j < evidences[i].ecnt ; ++j )
			{
				printf( "%s\tCLASS\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"adrenal.%d\"; transcript_id \"adrenal.%d\"; exon_number \"%d\";\n",
						chrom, evidences[i].exons[j].start, evidences[i].exons[j].end, strand,
						i, i, j ) ;
			}
		}
		exit( 1 ) ;*/
		eviExonCnt = ExtractEvidenceExons( &evidenceExons, evidences, eviCnt ) ;
	}
	else
		eviCnt = -1 ;
	TTS_Read( chrom ) ;
	if ( VERBOSE )
	{
		printf( "# Preprocessing %s\n", chrom ) ;
		fflush( stdout ) ;
	}
	
	seCnt = 0 ;
	noiseCnt = 0 ;
	noiseDepth = 0 ;
	noiseSqDepth = 0 ;
	PreprocessSplices( chrom, evidences, eviCnt ) ;
	//exit( 1 ) ;
	return 1 ;
}

/**
  Find a possible next region.

  return 0: if EOF.
*/
int FindPossibleNext( char *chrom, struct _splice ret_splices[], int &spliceCnt, int &start, int &end, struct _evidence *evidences, int &eviCnt ) 
{
	int i, j, k ;
	static int usedScnt = 0 ;
	static int usedSingleExon = 0 ;
	static char curChrom[50] = "" ;
	if ( usedScnt >= scnt && usedSingleExon >= seCnt && strcmp( curChrom, "-2" ) )
	{
		strcpy( curChrom, "-2" ) ;
		start = 0 ;
		end = 0 ;
		spliceCnt = 0 ;
		strcpy( chrom, curChrom ) ;
		return 1 ;
	}
	else if ( usedScnt >= scnt && usedSingleExon >= seCnt )
	{
		scnt = 0 ;
		while ( scnt == 0 )
		{
			if ( !GetSplices( curChrom, evidences, eviCnt ) )
				return 0 ;
		}
		usedScnt = 0 ;
		usedSingleExon = 0 ;
	}

	strcpy( chrom, curChrom ) ;
	if ( usedSingleExon < seCnt && ( usedScnt >= scnt || singleExon[ usedSingleExon ][0] < splices[ usedScnt ].pos ) )
	{
		start = singleExon[ usedSingleExon ][0] ;
		end = singleExon[ usedSingleExon ][1] ;
		spliceCnt = 0 ;
		++usedSingleExon ;
		return 1 ;
	}

	for ( k = usedScnt ; k < scnt ; ++k )
	{
		if ( spliceInfo[k].regionId != spliceInfo[ usedScnt ].regionId )
			break ;
	}
	
	start = spliceInfo[ usedScnt ].soft[0] ;
	//if ( start == -1 )
	//	start = splices[ usedScnt ].pos ;
	end = spliceInfo[ k - 1 ].soft[1] ;
	//if ( end == -1 )
	//	end = splices[ k - 1 ].pos ;
	//printf( "possible region: %d %d %d\n", splices[ usedScnt ].pos, splices[k - 1].pos, end ) ;
	
	spliceCnt = k - usedScnt ;
	memcpy( ret_splices, &splices[ usedScnt ], sizeof( splices[0] ) * spliceCnt ) ;

	for ( i = 0 ; i < spliceCnt ; ++i )
		ret_splices[i].otherInd -= usedScnt ;
	usedScnt = k ;
	
	return 1 ;
}

int FindRegions_Next( char *chrom, struct _splice splices[], int &scnt, int &start, int &end, struct _evidence *evidences, int &eviCnt ) 
{
	int i ;
	while ( 1 )
	{
		if ( !FindPossibleNext( chrom, tmpSplices, scnt, start, end, evidences, eviCnt ) )
		{
			if ( spliceInfo != NULL )
				free( spliceInfo ) ;
			return 0 ;
		}
		/*printf( "### %s %d %d %d\n", chrom, startCnt, endCnt, siCnt ) ;	
		  for ( i = 0 ; i < startCnt ; ++i )
		  printf( "%d\n", start[i].pos ) ;
		  for ( i = 0 ; i < endCnt ; ++i )
		  printf( "%d\n", end[i].pos ) ;*/
		//if ( scnt == 0 && start == -1 && end == -1 )
		//	continue ;
		/*for ( i = 0 ; i < endCnt ; ++i )
		  if ( end[i].pos > start[0].pos )
		  break ;
		  if ( i > 0 )
		  continue ;

		  for ( i = startCnt - 1 ; i >= 0 ; --i )
		  if ( start[i].pos < end[ endCnt - 1 ].pos )
		  break ;
		  if ( i < startCnt - 1 )
		  continue ;*/
		/*for ( i = 0 ; i < scnt ; ++i )
			printf( "%d %d\n", splices[i].pos, splices[i].type ) ;
		printf( "%d (%d %d)\n", scnt, start, end ) ;*/
		
		if ( scnt < MAX_POINT )
		{
			for ( i = 0 ; i < scnt ; ++i )
				splices[i] = tmpSplices[i] ;
		}
		else
		{
			// We only keep top-MAX_POINT splices.
			int *cnt = ( int * )malloc( sizeof( int ) * scnt ) ;
			int k, cut ;
			for ( i = 0 ; i < scnt ; ++i )
				cnt[i] = tmpSplices[i].support ;
			qsort( cnt, scnt, sizeof( int ), CompInt ) ;
			cut = cnt[scnt - MAX_POINT] ;
			for ( i = 0, k = 0 ; i < scnt ; ++i )
				if ( tmpSplices[i].support > cut )
				{
					splices[k] = tmpSplices[i] ;
					tmpSplices[i].pos = k ; // Reuse pos for ind
					++k ;
				}
				
			for ( i = 0, k = 0 ; i < scnt ; ++i )
				if ( tmpSplices[i].support > cut )
					splices[k].otherInd = tmpSplices[ tmpSplices[i].otherInd ].pos ;
			scnt = k ;
			free( cnt ) ;
		}
		int s0=-1, s1 = -1, e0=-1, e1 = -1 ;	
		int startCnt = 0, endCnt = 0 ;
		for ( i = 0 ; i < scnt ; ++i )
		{
			if ( splices[i].pos == -1 )
				continue ;
			if ( splices[i].type == 1 )
			{
				if ( s0 == -1 )
					s0 = splices[i].pos ;
				s1 = splices[i].pos ;
				++startCnt ;
			}
			else
			{
				if ( e0 == -1 )
					e0 = splices[i].pos ;
				e1 = splices[i].pos ;
				++endCnt ;
			}
		}

		if ( start != -1 )
		{
			s0 = start ;
			if ( s1 == -1 )
				s1 = start ;
			++startCnt ;
		}
		if ( end != -1 )
		{
			if ( e0 == -1 )
				e0 = end ;
			e1 = end ;
			++endCnt ;
		}
		//printf( "%d %d (%d %d) (%d %d)\n", startCnt, endCnt, s0, s1, e0, e1 ) ;	
		if ( startCnt == 0 || endCnt == 0 )
			continue ;

		//if ( startCnt == 1 && endCnt == 1 && start != -1 && end != -1 )
		//	continue ;

		if ( e0 < s0 || e1 < s1 )
			continue ;
		//if ( scnt == 0 )//&& start == -1 && end == -1 )
		//	continue ;
		//if ( scnt == 1 && start == -1 && end == -1 )
		//	continue ;
		//if ( startCnt ==1 && endCnt == 1 && start[0].strong == 0 && end[0].strong == 0  )
		//	continue ;
		//if ( end[0].pos <= start[0].pos || start[ startCnt - 1 ].pos >= end[ endCnt - 1 ].pos )
		//	continue ;

		//printf( "=== %d %d %d\n", scnt, start, end ) ;	
		
		return 1 ;
	}
}
