#include "TTS.h"

#define MAX_TTS 500000

struct _TTSsite *tts = NULL ;
int ttsCnt ;
FILE *fpTTS ;

void TTS_Init( char *file ) 
{
	fpTTS =fopen( file, "r" ) ;
	if ( !fpTTS )
		printf( "Could not find TTS file %s. Program will ignore it.\n", file ) ;
	tts = ( struct _TTSsite * )malloc( sizeof( struct _TTSsite ) * MAX_TTS ) ;
}

void TTS_Read( char *chrom )
{
	static char polyaChrom[10] = "", strand[3] = "" ;
	static int pos = -1, code = -1, sigPos ;
	if ( tts == NULL )
	{
		ttsCnt = -1 ;
		return ;
	}
	ttsCnt = 0 ;
	if ( !strcmp( polyaChrom, chrom ) )
	{
		tts[ ttsCnt ].sigPos = sigPos ;
		tts[ ttsCnt ].pos = pos ; 
		tts[ ttsCnt ].strand = ( strand[0] == '+' ? 1 : 0 ) ;
		tts[ ttsCnt ].code = code ;
		++ttsCnt ;
	}

	while ( fscanf( fpTTS, "%s %d %d %s %d", polyaChrom, &sigPos, &pos, strand, &code  ) != EOF )
	{
		if ( strcmp( polyaChrom, chrom ) )
			break ;

		tts[ ttsCnt ].sigPos = sigPos ;
		tts[ ttsCnt ].pos = pos ; 
		tts[ ttsCnt ].strand = ( strand[0] == '+' ? 1 : 0 ) ;
		tts[ ttsCnt ].code = code ;
		++ttsCnt ;
	}
}

struct _TTSsite *TTS_GetCurrentTTS( int &cnt )
{
	cnt = ttsCnt ;
	return tts ;
}
