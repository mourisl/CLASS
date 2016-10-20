/**
  Store the TTS sites reading from file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _TTSsite
{
	int sigPos ;
	int pos ;
	int strand ;
	int code ; // code for the signal
} ;

void TTS_Init( char *file ) ;
void TTS_Read( char *chrom ) ;
struct _TTSsite *TTS_GetCurrentTTS( int &cnt ) ;
