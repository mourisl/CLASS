/**
	Driver for the test program
	It go through all the regions, generate the input for the LPAssembly and run it.
	Suppose the order of the chromosomes of the headers file and the depth file are the same.
	Usage: ./a.out headers reads_file depth_file right_comb_file > output
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "SolveRegion.h"
#include "Reads.h"
#include "FindRegions.h"
#include "TranscriptDecider.h"

#define TRANSCRIPT_GO 1

struct _read reads[ MAX_READ ] ;
int depth[MAX_LENGTH] ;
struct _exon exons[MAX_EXON] ;
int exonCnt = 0 ;
struct _splice region_splices[MAX_POINT] ;
int READS_LENGTH ;
int FRAG_LENGTH, FRAG_STD ;
long long TOTAL_READ_COUNT ;
bool VAR_RD_LEN ;

extern int BIN_RANGE, BIN_USE ;
extern bool USE_SET_COVER ;
extern double IR_ALPHA ;
extern double FPKM_FRACTION ;

bool VERBOSE ;
int NUM_OF_THREADS = 1 ;
pthread_attr_t pthreadAttr ;
char *gtfIdLabel = NULL ;
//FILE *fpPolyA = NULL ;

extern int currSolveRegionThreadsCnt ;
extern pthread_mutex_t solveRegionMutex ;
extern pthread_mutex_t allExonsMutex ;
extern pthread_cond_t idleSolveRegionCond ;
extern pthread_cond_t clearSolveRegionCond ;

double timeval_subtract (struct timeval *result,
		struct timeval *x,
		struct timeval *y)
// Subtract the `struct timeval' values X and Y,
// storing the result in RESULT.
// Return 1 if the difference is negative, otherwise 0.
{
	// Perform the carry for the later subtraction by updating y.
	if (x->tv_usec < y->tv_usec)
	{
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}

	if (x->tv_usec - y->tv_usec > 1000000)
	{
		int nsec = (y->tv_usec - x->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	// Compute the time remaining to wait.
	// tv_usec is certainly positive.
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_usec = x->tv_usec - y->tv_usec;

	// Return 1 if result is negative.
	return result->tv_sec + (double)result->tv_usec / 1000000.0 ;
}


int Comp_Int_Desc( const void *p1, const void *p2 )
{
	// Into descending order.
	return ( *( int * )p2 ) - ( *( int * )p1 ) ;   
}

int Comp_Exons( const void *p1, const void *p2 )
{
	struct _exon *a, *b ;
	a = ( struct _exon *)p1 ;
	b = ( struct _exon *)p2 ;

	if ( a->start != b->start )
		return a->start - b->start ;
	if ( a->end != b->end )
		return a->end - b->end ;
	return a->strand - b->strand ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int offset, startCnt, endCnt ;
	int position, d ;
	char chrom[50], preChrom[50] = "-1" ;
	char buffer[2048] ;
	FILE *fpHeader, *fpDepth, *fpOut, *fpComb, *fpOneComb ;
	bool started ;
	int tmp ;
	int siCnt, readCnt ;
	int useEvidence, usePolyA, useOtherJunctionFile ;
	int softStart, softEnd, scnt ;
	struct _readFile fpReads ;
	struct _evidence *evidences ;
	int eviCnt = 0 ;

	struct timeval beginTime, endTime, elapsedTime ;
	gettimeofday(&beginTime,NULL) ;

	if ( argc == 1 )
	{
		//printf( "Usage: ./a.out exon_file reads_file depth_file evidence_file splice_file chroms.size > output\n" ) ;
		printf( "Usage: ./a.out prefix_of_the_sam_depth_splice [options]> output\n"
				"Options:\n"
				"\t-p number_of_threads: specify the number of worker threads (default:1)\n"
				"\t-j junction: the path to the splice junction file\n"
				"\t-e evidence: the path to the evidence files\n"
				//"\t-anno annotation: the path to the annoation file. CLASS will only consider the transcript in specified in the set"
				"\t--var-rd-len: extensive variable read lengths, i.e. reads after trimming (default: auto-detect)\n"
				"\t--set-cover: use set cover to build transcripts from splicing graph\n"
				"\t-F f: do not report the transcripts whose abundance level is lower than f*|most expressed transcript| in a gene\n"
				"\t-l label: add a prefix and a \"_\" to the ids in the GTF file (default: not used)\n"
				"\t--verbose: also output the procedure of CLASS\n" 
				) ;
		return 1 ; 
	}
	//fpHeader = fopen( argv[1], "r" ) ;
	// Find the reads length

	useEvidence = -1 ;
	usePolyA = -1 ;
	useOtherJunctionFile = -1 ;
	VAR_RD_LEN = false ;
	USE_SET_COVER = false ;
	VERBOSE = false ;
	pthread_attr_init( &pthreadAttr ) ;
	pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_DETACHED ) ;

	currSolveRegionThreadsCnt = 0 ;
	pthread_mutex_init( &solveRegionMutex, NULL ) ;
	pthread_mutex_init( &allExonsMutex, NULL ) ;
	pthread_cond_init( &idleSolveRegionCond, NULL ) ;
	pthread_cond_init( &clearSolveRegionCond, NULL ) ;

	// Output the command lines
	printf( "#CLASS_v2.1.7\n#" ) ;
	for ( i = 0 ; i < argc - 1 ; ++i )
	{
		printf( "%s ", argv[i] ) ;
	}
	printf( "%s\n", argv[i] ) ;
	
	//pthread_attr_setstacksize(...) ;
	int extent[2] ;
	//readCnt = ExtractReads( fpReads, "chrM", 1, 16571, reads, extent ) ;
	for ( i = 2 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-e" ) )
		{
			useEvidence = i + 1 ;
			/*FILE *fp = NULL ;
			  fp = fopen( argv[i + 1], "r" ) ;
			  if ( fp == NULL )
			  {
			  printf( "Could not find evidence file %s.\n", argv[i + 1] ) ;
			  return 0 ;
			  }*/
			i += 1 ;
		}
		else if ( !strcmp( argv[i], "--bin" ) )
		{
			BIN_RANGE = atoi( argv[i + 1 ] ) ;
			BIN_USE = atoi( argv[i + 2] ) ;
			//printf( "### %d %d\n", BIN_RANGE, BIN_USE ) ;
			i += 2 ;
		}
		else if ( !strcmp( argv[i], "-j" ) )
		{
			useOtherJunctionFile = i + 1 ;
			i += 1 ;		
		}
		else if ( !strcmp( argv[i], "--set-cover" ) )
		{
			USE_SET_COVER = true ;
		}
		/*else if ( !strcmp( argv[i], "-polya" ) )
		  {
		  usePolyA = i + 1 ; 
		  i += 1 ;
		  }*/
		else if ( !strcmp( argv[i], "--iralpha" ) )
		{
			IR_ALPHA = atof( argv[i + 1] ) ;
			i += 1 ;
		}
		else if ( !strcmp( argv[i], "--var-rd-len" ) )
		{
			VAR_RD_LEN = true ;
		}
		else if ( !strcmp( argv[i], "--verbose" ) )
		{
			VERBOSE = true ;
		}
		else if ( !strcmp( argv[i], "-F" ) )
		{
			FPKM_FRACTION = atof( argv[i + 1] ) ;
			i += 1 ;
		}
		else if ( !strcmp( argv[i], "-p" ) )
		{
			NUM_OF_THREADS = atoi( argv[i + 1] ) ;
			i += 1 ;
		}
		else if ( !strcmp( argv[i], "-l" ) )
		{
			gtfIdLabel = argv[i + 1] ;
			i += 1 ;
		}
		else
		{
			fprintf( stderr, "Unknown parameter %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}


	fpReads = OpenReadFile( argv[1] ) ;
	//rewind( fpReads ) ;
	sprintf( buffer, "%s.depth", argv[1] ) ;
	fpDepth = NULL ;
	fpDepth = fopen( buffer, "r" ) ;
	if ( fpDepth == NULL )
	{
		fprintf( stderr, "Could not find depth file %s\n", buffer ) ;
		return 0 ;
	}

	if ( useEvidence != -1 )
		Evidence_Init( argv[useEvidence], &evidences ) ;
	else
		evidences = NULL ;	

	if ( usePolyA != -1 )
	{
		TTS_Init( argv[ usePolyA ] ) ;
	}

	if ( useOtherJunctionFile == -1 )
		FindRegions_Init( argv[1], NULL ) ;
	else
		FindRegions_Init( argv[1], argv[ useOtherJunctionFile ] ) ;

	SplicesInformation_Init() ;

	if ( VERBOSE )
		printf( "# Getting alignments information.\n" ) ;
	GetReadsInfo( fpReads, READS_LENGTH, FRAG_LENGTH, FRAG_STD, TOTAL_READ_COUNT ) ;
	if ( VERBOSE )
		printf( "# Alignments information of %d reads: Read_Length=%d Fragment_Length=%d Fragment_Stddev=%d\n",
				TOTAL_READ_COUNT, READS_LENGTH, FRAG_LENGTH, FRAG_STD ) ;
	CloseReadFile( fpReads ) ;

	fpReads = OpenReadFile( argv[1] ) ;
	//memset( exons, -1, sizeof( exons ) ) ;	
	for ( i = 0 ; i < MAX_EXON ; ++i )
	{
		exons[i].nsize = 1 ;
		exons[i].psize = 1 ;
		exons[i].next = ( int * )malloc( sizeof( int ) * 1 ) ;
		exons[i].prev = ( int * )malloc( sizeof( int ) * 1 ) ;
	}
	TranscriptDecider_Init( argv[1] ) ;
	//TODO: get the evidences for the next chrom.
	//printf( "hi\n" ) ;
	position = -1 ; // Indicator whether we have read the depth file or not.
	while ( FindRegions_Next( chrom, region_splices, scnt, softStart, softEnd, evidences, eviCnt )  )//&& !TRANSCRIPT_GO )
	{
		/*strcpy( chrom, "chrM ") ;
		  scnt = 0 ;
		  softEnd = 16571 ;
		  softStart = 1 ;*/
		//printf( "Found a new region:\n") ;
		int start = softStart, end = softEnd ;
		if ( start == -1 )
		{
			for ( i = 0 ; i < scnt ; ++i )
			{
				if ( region_splices[i].pos != -1 )
				{
					start = region_splices[i].pos ;
					break ;
				}
			}
		}		
		if ( end == -1 )
		{
			for ( i = scnt - 1 ; i >= 0 ; --i )
			{
				if ( region_splices[i].pos != -1 )
				{
					end = region_splices[i].pos ;
					break ;
				}
			}
		}
		//printf( "== %d %d %d\n", scnt, start, end ) ;
		if ( strcmp( chrom, preChrom ) )
		{
			ClearOverlapReadsBuffer() ;	
			if ( strcmp( preChrom, "-1" ) && strcmp( preChrom, "-2" ) )
			{
				// TODO: Wait for all the solveRegionThread finish
				pthread_mutex_lock( &solveRegionMutex ) ;
				if ( currSolveRegionThreadsCnt > 0 )
					pthread_cond_wait( &clearSolveRegionCond, &solveRegionMutex ) ;
				pthread_mutex_unlock( &solveRegionMutex ) ;
				//printf( "%d\n", exonCnt ) ;
				// Sort the exons
				qsort( exons, exonCnt, sizeof( exons[0] ), Comp_Exons ) ;	
				
				//	printf( "### %s %d\n", preChrom, exonCnt ) ;
				if ( VERBOSE )
				{
					printf( "# Found %d exons in %s.\n", exonCnt, preChrom ) ;
					for ( i = 0 ; i < exonCnt ; ++i )
					{
						printf( "# Exon %d: %s %d %d\n", i, preChrom, exons[i].start, exons[i].end ) ;
						/*printf( "\t# of splice sites before %d: ", exons[i].pcnt ) ;
						for ( j = 0 ; j < exons[i].pcnt ; ++j )
							printf( "%d ", exons[i].prev[j] ) ;
						printf( "\n" ) ;
						printf( "\t# of splice sites after %d: ", exons[i].ncnt ) ;
						for ( j = 0 ; j < exons[i].ncnt ; ++j )
							printf( "%d ", exons[i].next[j] ) ;
						printf( "\n") ;*/
					}

				}
				TranscriptDecider_Go( preChrom, exons, exonCnt, evidences, eviCnt ) ;
				for ( i = 0 ; i < exonCnt ; ++i )
				{
					if ( exons[i].nsize > 1 )
					{
						free( exons[i].next ) ;
						exons[i].next = (int *)malloc( sizeof( int ) * 1 ) ;
						exons[i].nsize = 1 ;
					}

					if ( exons[i].psize > 1 )
					{
						free( exons[i].prev ) ;
						exons[i].prev = (int *)malloc( sizeof( int ) * 1 ) ;
						exons[i].psize = 1 ;
					}
				}

				currSolveRegionThreadsCnt = 0 ;
			}
			strcpy( preChrom, chrom ) ;
			SplicesInformation_Reset() ;
			exonCnt = 0 ;
		}
		if ( !strcmp( chrom, "-2" ) )
			continue ;

		if ( VERBOSE )
		{
			printf( "# Solving region: %s %d %d (%d exons so far)\n", chrom, start, end, exonCnt ) ; 	
			fflush( stdout ) ;
		}
		//fpOneComb = fopen( "one_right_comb.out", "w" ) ;	
		//for ( i = 0 ; i < startCnt ; ++i )
		//	inputStart[i].merge = i ;
		//for ( i = 0 ; i < endCnt ; ++i )
		//	inputEnd[i].merge = i ;

		/*printf( "%s %d %d\n", chrom, startCnt, endCnt ) ;
		  for ( i = 0 ; i < startCnt ; ++i )
		  printf( "%d\n", inputStart[i].pos ) ;
		  for ( i = 0 ; i < endCnt ; ++i )
		  printf( "%d\n", inputEnd[i].pos ) ;
		  for ( i = 0 ; i < siCnt ; ++i )
		  printf( "si %d: %d %d\n", i, spliceIndices[i].startPos, spliceIndices[i].endPos ) ;*/
		started = false ;
		i = 0 ;
		if ( end - start + 10 < MAX_LENGTH )
			memset( depth, 0, sizeof( depth[0] ) * ( end - start + 10 ) ) ;	
		else
			memset( depth, 0, sizeof( depth ) ) ;
		while ( 1 )
		{
			if ( position == -1 )
				if ( fscanf( fpDepth, "%s %d %d", buffer, &position, &d ) == EOF )
					break ;
			
			//printf( "main: %s %d %d\n", buffer, position, d ) ;
			tmp = strcmp( buffer, chrom ) ;
			if ( ( !tmp && position > end ) ||
					( tmp && started ) )
			{
				// End for a region
				break ;
			}
			if ( !tmp && position < start )
			{
				position = -1 ;
				continue ;
			}

			/*if ( ( tmp && !started ) ||
					( !tmp && position < start ) )
			{
				// Before the region
				position = -1 ;
				continue ;
			}*/
			if ( tmp )
			{
				int chromId = GetChromIdFromName( chrom ) ;
				int depthChromId = GetChromIdFromName( buffer ) ;

				if ( depthChromId > chromId )
					break ;
				else // depthChromId < chromId 
				{
					char prevBuffer[2048] ;
					bool skipped = false ;
					strcpy( prevBuffer, buffer ) ;
					while ( 1 )
					{
						if ( fscanf( fpDepth, "%s %d %d", buffer, &position, &d ) == EOF )
						{
							skipped = true ;
							break ;
						}

						if ( !strcmp( prevBuffer, buffer ) )
							continue ;
						strcpy( prevBuffer, buffer ) ;
						depthChromId = GetChromIdFromName( buffer ) ;

						if ( depthChromId < chromId )
							continue ;
						else if ( depthChromId > chromId )
						{
							skipped = true ;
							break ;
						}
						else
							break ;
					}
					if ( skipped )
						break ;
					else 
						continue ; // The position will be kept when jumped back to the start of the while-loop
				}
			}

			started = true ;
			if ( position - start < MAX_LENGTH )
				depth[ position - start ] = d ;
			position = -1 ;
		}
		int extent[2] ;
		readCnt = ExtractReads( fpReads, chrom, start, end, reads, extent ) ;
		if ( VERBOSE )
		{
			printf( "# Found %d reads in this region.\n", readCnt ) ; 	
			fflush( stdout ) ;
		}
		//printf( "readcnt: %d\n", readCnt ) ;
		//if ( start > 9908278 )
		//	return 0 ;
		//printf( "%d %d %d\n", start, end, readCnt ) ;
		//readCnt = 0 ;
		// Output the "right" combination

		/*fgets( buffer, sizeof( buffer ), fpComb ) ;
		  sscanf( buffer, "%d", &j ) ;	
		  fputs( buffer, fpOneComb ) ;
		  for ( i = 0 ; i < j ; ++i )
		  {
		  fgets( buffer, sizeof( buffer ), fpComb ) ;
		  fputs( buffer, fpOneComb ) ;
		  }	
		  fclose( fpOneComb ) ;*/
		int tmp = exonCnt ;
		if ( scnt == 0 && start != -1 && end != -1 )
		{
			// Add lock here
			pthread_mutex_lock( &allExonsMutex ) ;
			exons[ exonCnt ].start = start ;
			exons[ exonCnt ].end = end ;
			exons[ exonCnt ].strand = -1 ;
			exons[ exonCnt ].pcnt = exons[ exonCnt ].ncnt = 0 ;
			//printf( "Single-exon transcript: %d %d %d\n", exonCnt, start, end ) ; //exit( 1 ) ;
			++exonCnt ;
			pthread_mutex_unlock( &allExonsMutex ) ;
		}
		else //if ( 0 )
		{
			int tplus = 1, tminus = 1 ;
			//printf( "Solve region\n") ;
			while ( SolveRegion_Wrapper( chrom, region_splices, scnt, softStart, softEnd, 
						depth, reads, readCnt, evidences, eviCnt, exons, exonCnt, tplus, tminus, false ) == -1 ) 
			{
				int *supports = ( int * )malloc( sizeof( int ) * scnt ) ;

				//printf( "wow\n" ) ; exit( 1 ) ;
				for ( i = 0, k = 0 ; i < scnt ; ++i )
				{
					if ( region_splices[i].strand == 1 && region_splices[i].otherInd >= 0 && region_splices[i].otherInd < scnt 
							&& region_splices[i].pos != -1 && region_splices[i].support >= tplus )
					{
						supports[k] = region_splices[i].support ;
						++k ;
					}
				}
				qsort( supports, k, sizeof( int ), Comp_Int_Desc ) ;

				for ( i = 0 ; i < k - 1 ; ++i )
				{
					if ( supports[i] >= 10 *  supports[i + 1] )
						break ;
				}

				if ( i < k - 1 )
					tplus = supports[i] ;
				else if ( tplus <= tminus && tplus < 10 )
					++tplus ;

				for ( i = 0, k = 0 ; i < scnt ; ++i )
				{
					if ( region_splices[i].strand == 0 && region_splices[i].otherInd >= 0 && region_splices[i].otherInd < scnt 
							&& region_splices[i].pos != -1 && region_splices[i].support >= tminus )
					{
						supports[k] = region_splices[i].support ;
						++k ;
					}
				}
				qsort( supports, k, sizeof( int ), Comp_Int_Desc ) ;

				for ( i = 0 ; i < k - 1 ; ++i )
				{
					//printf( "%d ", supports[i] ) ;
					if ( supports[i] >= 10 *  supports[i + 1] )
						break ;
				}
				if ( i < k - 1 )
				{
					tminus = supports[i] ;
				}
				else if ( tminus <= tplus && tminus < 10 )
					++tminus ; 

				if ( tplus >= 10 && tminus >= 10 )
				{
					// Too much splice junction whose other end is out of the region.
					// and t overflows.
					SolveRegion_Wrapper( chrom, region_splices, scnt, softStart, softEnd, 
							depth, reads, readCnt, evidences, eviCnt, exons, exonCnt, tplus, tminus, true ) ;
					free( supports ) ;
					break ;
				}
				free( supports ) ;
			}
			//printf( "End solving region\n" ) ;
			//for ( i = 0 ; i < exonCnt ; ++i )
			//	printf( "%d %d %d\n", exons[i].start, exons[i].end, exons[i].strand ) ;
		}
		if ( VERBOSE )
		{
			printf( "# Finish solving region: %s %d %d (%d exons so far)\n", chrom, start, end, exonCnt ) ; 	
			fflush( stdout ) ;
		}
		//if ( exons[ exonCnt - 1 ].start >= 10000000 )
		//	break ;
		/*if ( tmp == exonCnt )
		  {
		  printf( "### WRONG\n") ;
		  exit( 1 ) ;
		  }*/	
		//printf( "Finish a region\n" ) ;
	}
	//#if !TRANSCRIPT_GO
	/*FILE *fpExon = fopen( "tmp.exons", "w" ) ;

	  fwrite( &exonCnt, sizeof( exonCnt ), 1, fpExon ) ;
	  for ( i = 0 ; i < exonCnt ; ++i )
	  fwrite( &exons[i], sizeof( exons[0] ), 1, fpExon ) ;
	  fclose( fpExon ) ;
	//printf( "%d\n", exonCnt ) ;*/
	//#endif

	/*FILE *fpExon = fopen( "tmp.exons", "r" ) ;
	  fread( &exonCnt, sizeof( exonCnt ), 1, fpExon ) ;
	  for ( i = 0 ; i < exonCnt ; ++i )
	  fread( &exons[i], sizeof( exons[0] ), 1, fpExon ) ;
	  fclose( fpExon ) ;
	  strcpy( preChrom, "chr12" ) ;
	  printf( "hi %d\n", exonCnt ) ;*/
	//exit( 1 ) ;
	//printf( "hi %d\n", exonCnt ) ; fflush( stdout ) ;
	//			for ( i = 0 ; i < exonCnt ; ++i )
	//				printf( "%d %d %d\n", exons[i].start, exons[i].end, exons[i].strand ) ;

	if ( strcmp( preChrom, "-2" ) )
	{
		if ( VERBOSE )
		{
			printf( "# Found %d exons in %s.\n", exonCnt, preChrom ) ;
			for ( i = 0 ; i < exonCnt ; ++i )
			{
				printf( "# Exon %d: %s %d %d\n", i, preChrom, exons[i].start, exons[i].end ) ;
				/*printf( "\t# of splice sites before %d: ", exons[i].pcnt ) ;
				for ( j = 0 ; j < exons[i].pcnt ; ++j )
					printf( "%d ", exons[i].prev[j] ) ;
				printf( "\n" ) ;
				printf( "\t# of splice sites after %d: ", exons[i].ncnt ) ;
				for ( j = 0 ; j < exons[i].ncnt ; ++j )
					printf( "%d ", exons[i].next[j] ) ;
				printf( "\n") ;*/
			}
		}
		fflush( stdout ) ;

		pthread_mutex_lock( &solveRegionMutex ) ;
		if ( currSolveRegionThreadsCnt > 0 )
			pthread_cond_wait( &clearSolveRegionCond, &solveRegionMutex ) ;
		pthread_mutex_unlock( &solveRegionMutex ) ;
		//sort the exons
		qsort( exons, exonCnt, sizeof( exons[0] ), Comp_Exons ) ;	

		TranscriptDecider_Go( preChrom, exons, exonCnt, evidences, eviCnt ) ;
	}
	
	for ( i = 0 ; i < MAX_EXON ; ++i )
	{
		free( exons[i].next ) ;
		free( exons[i].prev ) ;
	}
	//printf( "### %s %d\n", preChrom, exonCnt ) ;
	//fprintf( stderr, "Finished\n" ) ; fflush( stderr ) ;	
	//fclose( fpHeader ) ;
	//fclose( fpReads ) ;
	CloseReadFile( fpReads ) ;

	//if ( fpPolyA )
	//	fclose( fpPolyA ) ;

	pthread_mutex_destroy( &solveRegionMutex ) ;
	pthread_mutex_destroy( &allExonsMutex ) ;
	pthread_cond_destroy( &idleSolveRegionCond ) ;
	pthread_cond_destroy( &clearSolveRegionCond ) ;
	fclose( fpDepth ) ;

	FindRegions_Destroy() ;
	TranscriptDecider_Destroy() ;

	gettimeofday( &endTime, NULL ) ;	

	if ( VERBOSE )
	{
		printf( "# CLASS finishes successfully with %lf seconds.\n", timeval_subtract( &elapsedTime, &endTime, &beginTime ) ) ;
	}
	return 0 ;

}
