#include <stdio.h>
#include <stdint.h>

#include <string.h>
#include <vector>
#include <map>
#include <fstream>

char usage[] = 
	"./addXS ref.fa [OPTIONS] < in.sam > out.sam\n"
		"From bam to bam: samtools view -h in.bam | ./addXS ref.fa | samtools view -bS - > out.bam\n" ;
char samLine[65537] ;
char seq[65537] ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2,
			-1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, 3,
			-1, -1, -1, -1, -1, -1 } ;
char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

struct _pair
{
	int a, b ;
} ;

class BitSequence 
{
private:
	int len ;
	//std::vector<uint32_t> sequence ;
	uint32_t *sequence ;

public:
	BitSequence() { len = 0 ; sequence = NULL ;} 
	BitSequence( int l )
	{
		len = 0 ;
		sequence = new uint32_t[ l / 16 + 1 ] ;
	}
	
	~BitSequence() 
	{
		//if ( sequence != NULL )
		//	delete[] sequence ;
	}
	
	void SetLength( int l )
	{
		len = 0 ;
		if ( sequence != NULL )
			delete[] sequence ;
		sequence = new uint32_t[ l / 16 + 1 ] ;
	}

	int GetLength()
	{
		return len ;
	}

	void Append( char c )
	{
		if ( ( len & 15 ) == 0 )
		{
			sequence[ len / 16 ] = 0 ;
		}
		++len ;
		//printf( "%d %c\n", len, c ) ;
		Set( c, len - 1 ) ;
	}
	
	// pos is 0-based coordinate
	// notice that the order within one 32 bit butcket is reversed
	void Set( char c, int pos ) 
	{
		if ( pos >= len )		
			return ;

		if ( c >= 'a' && c <= 'z' )
		{
			c = c - 'a' + 'A' ;
		}
		if ( c == 'N' )
			c = 'A' ;

		int ind = pos >> 4 ;
		int offset = pos & 15 ;
		int mask = ( (int)( nucToNum[c - 'A'] & 3 ) ) << ( 2 * offset ) ;
		sequence[ind] = sequence[ind] | mask ;
		//if ( c != 'A' )
		//	printf( "%d: %c %c %d %d : %d\n", pos, c, Get(pos), ind, offset, mask ) ;
		//Print() ;
	}

	char Get( int pos )
	{
		if ( pos >= len )
			return 'N' ;

		int ind = pos >> 4 ;
		int offset = pos & 15 ;
		//printf( "%d: %d\n", pos, sequence[ind] ) ;
		return numToNuc[ ( ( sequence[ind] >> ( 2 * offset ) ) & 3 ) ] ;
	}

	void Print()
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			printf( "%c", Get( i ) ) ;
		printf( "\n" ) ;
	}

	void Print( FILE *fp, int start, int end, bool rc )
	{	
		if ( !rc )
		{
			for ( int i = start ; i <= end ; ++i )
				fprintf( fp, "%c", Get( i ) ) ;
		}
		else
		{
			for ( int i = end ; i >= start ; --i )
			{
				char c = Get( i ) ;
				if ( c == 'A' )
					c = 'T' ;
				else if ( c == 'C' )
					c = 'G' ; 
				else if ( c == 'G' )
					c = 'C' ; 
				else //if ( c == 'T' )
					c = 'A' ; 
				fprintf( fp, "%c", c ) ;
			} 
		}
	}
} ;


int main( int argc, char *argv[] )
{
	int i, j, k ;
	
	std::vector<BitSequence> genome ;
	std::map<std::string, int> chrToChrId ;

	char readid[200], chrom[50], mapq[10], cigar[1000], mateChrom[50] ;
	int start, mstart, flag, tlen ; // read start and mate read start
	int chrCnt = 0 ;
	int chrId ;

	if ( argc < 2 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}
	

	std::ifstream fpRef ;	
	fpRef.open( argv[1] ) ;
	std::string line ;
	
	int motifStrand[1024] ;
	
	/*motifStrand[ 0b10110010 ] = 1 ; // GT/AG
	motifStrand[ 0b01110001 ] = 0 ; // CT/AC
	motifStrand[ 0b10010010 ] = 1 ;// GC/AG
	motifStrand[ 0b01111001 ] = 0 ; // CT/GC
	motifStrand[ 0b00110001 ] = 1 ;	// AT/AC
	motifStrand[ 0b10110011 ] = 0 ; // GT/AT*/
	memset( motifStrand, -1, sizeof( motifStrand )) ;
	motifStrand[ 0xb2 ] = 1 ; // GT/AG
	motifStrand[ 0x71 ] = 0 ; // CT/AC
	motifStrand[ 0x92 ] = 1 ;// GC/AG
	motifStrand[ 0x79 ] = 0 ; // CT/GC
	motifStrand[ 0x31 ] = 1 ;	// AT/AC
	motifStrand[ 0xb3 ] = 0 ; // GT/AT
	 
	k = 0 ;
	while ( getline( fpRef, line ) )
	{
		//printf( "%s\n", line.c_str() ) ;
		int len = line.length() ;
		if ( line[0] == '>' )
		{
			//char *s = strdup( line.c_str() + 1 ) ;
			if ( chrCnt > 0 )
			{
				genome[ chrCnt - 1 ].SetLength( k ) ;
			}
			for ( i = 1 ; i < len ; ++i )
				if ( line[i] == ' ' || line[i] == '\t' )
					break ;
			chrToChrId[ line.substr( 1, i - 1 ) ] = chrCnt ;
			++chrCnt ;
			
			BitSequence bs ;
			genome.push_back( bs ) ;
			
			k = 0 ;
		}	
		else
		{
			k += len ;
		}			
	}
	genome[ chrCnt - 1 ].SetLength( k ) ;
	fpRef.close() ;
	
	fpRef.open( argv[1] ) ;
	while ( getline( fpRef, line ) )
	{
		//printf( "%s\n", line.c_str() ) ;
		int len = line.length() ;
		if ( line[0] == '>' )
		{
			//char *s = strdup( line.c_str() + 1 ) ;
			for ( i = 1 ; i < len ; ++i )
				if ( line[i] == ' ' || line[i] == '\t' )
					break ;
			chrId = chrToChrId[ line.substr( 1, i - 1 ) ] ;
		}	
		else
		{
			BitSequence &bs = genome[chrId] ;
			for ( i = 0 ; i < len ; ++i )
			{
				if ( ( line[i] >= 'A' && line[i] <= 'Z' ) ||
					( line[i] >= 'a' && line[i] <= 'z' ) )			
					bs.Append( line[i] ) ;
			}
		}			
	}	
	fpRef.close() ;
	
	
	
	while ( fgets( samLine, sizeof( samLine ), stdin ) != NULL )
	{
		if ( samLine[0] == '@' )
		{
			printf( "%s", samLine ) ;
			continue ;
		}
		
		sscanf( samLine, "%s %d %s %d %s %s %s %d %d %s", readid, &flag, chrom, &start, mapq, cigar, mateChrom, &mstart, &tlen, seq ) ;
		struct _pair segments[1000] ;
		int segCnt = 0 ;
		int num = 0 ;
		int len = 0 ;
	
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
				segments[segCnt].a = start ;
				segments[segCnt].b = start + len - 1 ;
				++segCnt ;
				
				start = start + len + num ;
				len = 0 ;
				num = 0 ;
			}
			else
			{
				len += num ;
				num = 0 ;
			}
		}
		
		if ( len > 0 )
		{
			segments[segCnt].a = start ;
			segments[segCnt].b = start + len - 1 ;
			++segCnt ;
		}
			
		if ( segCnt <= 1 || strstr( samLine, "XS:A:" ) != NULL )
		{
			printf( "%s", samLine ) ;
			continue ;
		}
			
		std::string schr( chrom ) ;
		int chrId = chrToChrId[schr] ;
		int strand = -1 ;
		
		len = strlen( samLine ) ;
		samLine[len - 1] = '\0' ;
		int motif = 0 ;
		for ( i = 0 ; i < segCnt - 1 ; ++i )
		{
			char m[4] ;
			m[0] = genome[ chrId ].Get( segments[i].b + 1 - 1  ) ;
			m[1] = genome[ chrId ].Get( segments[i].b + 2 - 1 ) ;
			m[2] = genome[ chrId ].Get( segments[i + 1].a - 2 - 1 ) ;
			m[3] = genome[ chrId ].Get( segments[i + 1].a - 1 - 1 ) ;
			motif = 0 ;
			for ( j = 0 ; j < 4 ; ++j )
				motif = ( motif << 2 ) + ( nucToNum[ m[j] - 'A' ] & 3 );
			strand = motifStrand[ motif ] ;
			if ( strand != -1 )
				break ;
		}
		if ( strand == 1 )
		{
			printf( "%s\tXS:A:+\n", samLine ) ;
		}
		else if ( strand == 0 )
		{
			printf( "%s\tXS:A:-\n", samLine ) ;
		}
		else
		{
			printf( "%s\tXS:A:?\n", samLine ) ;
		}
	}
	
	return 0 ;
}
