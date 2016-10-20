/*
  Clean BAM file 
  Usage: ./a.out input [option] >output. 
  	./a.out -h for help.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "sam.h"

#define LINE_SIZE 4097
#define QUEUE_SIZE 1000001
#define HASH_MAX 1000003

#define ABS(x) ((x)<0?-(x):(x))

char line[LINE_SIZE] ;
char strand ; // Extract XS field
bool secondary ;
int mateStart ;

struct _list
{
	bam1_t *b ;
	int start, end ;
	bool keep ;
	struct _list *next ;
} ;

struct _tree
{
	bam1_t *b ;
	int start, end ;
	//bool keep ;
	struct _list *list ;
	struct _tree *left, *right ; 
} ;

int qHead, qTail, qsHead, qsTail ;

samfile_t *fpout ;

void PrintHelp()
{
	printf( 
		"Clean BAM file.\n" 
		"Usage: ./a.out input\n"
	      ) ;
}

void cigar2string( bam1_core_t *c, uint32_t *in_cigar, char *out_cigar )
{
	int k, op, l ;
	char opcode ;

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
		}
		sprintf( out_cigar + strlen( out_cigar ), "%d%c", l, opcode ) ;
	}
}

// Return the pointer to the new node that replacing original t.
struct _tree *RemoveNode( struct _tree *t, struct _tree *parent, bool isLeft )
{
	struct _tree *p, *q, *prev ;
	//if ( t->keep )
	//	samwrite( fpout, t->b ) ;
	//free( t->list ) ;
	//bam_destroy1( t->b ) ;
	p = t->right ;
	if ( p == NULL )
	{
		q = t->left ;
		if ( q == NULL )
		{
			if ( isLeft )
				parent->left = NULL ;
			else
				parent->right = NULL ;
			free( t ) ;
			return NULL ;
		}
		else
		{
			*t = *q ;
			free( q ) ;
			return t ;
		}
	}

	q = p ;
	while ( q->left != NULL )
	{
		prev = q ;
		q = q->left ;
	}
	if ( q == p )
	{
		p->left = t->left ;
		*t = *p ;
		free( p ) ;
		return t ;
	}
	else
	{
		t->b = q->b ;
		//t->keep = q->keep ;
		t->start = q->start ;
		t->end = q->end ;
		t->list = q->list ;
		prev->left = q->right ;	
		free( q ) ;
		return t ;
	}
}

void SearchTree( bam1_t *b, int start, int end, bool keep, struct _tree *t, struct _tree *parent, bool leftChild, struct _list *list )
{
	int cmp ;
	while ( t != NULL && ( t->end < start || t->b->core.tid != b->core.tid ) )
	{
		t = RemoveNode( t, parent, leftChild ) ;
	}

	if ( t != NULL )
	{
		cmp = strcmp( bam1_qname( t->b ) , bam1_qname( b ) ) ;
		if ( cmp < 0 )
			SearchTree( b, start, end, keep, t->left, t, true, list ) ;
		else if ( cmp > 0 )
			SearchTree( b, start, end, keep, t->right, t, false, list ) ;
		else if ( ( t->b->core.flag & ( 0x40 + 0x80 ) ) == ( b->core.flag & ( 0x40 + 0x80 ) ) )
		{
			//printf( "hi\n" ) ;
			SearchTree( b, start, end, keep, t->left, t, true, list ) ;
		}
		else
		{
			if ( start == t->start || end == t->end )
			{
				if ( b->core.mtid != b->core.tid )
					keep = false ;
				else if ( b->core.mtid == b->core.tid && t->b->core.mtid == t->b->core.tid )
				{
					if ( ABS( b->core.pos - b->core.mpos ) < ABS( t->b->core.pos - t->b->core.mpos ) )
						//t->keep = false ;
						t->list->keep = false ;
					else
						keep = false ;
				}
			}
			else
			{
				// just overlap
				keep = false ;
				//t->keep = false ;
				t->list->keep = false ;
			}
			SearchTree( b, start, end, keep, t->left, t, true, list ) ;
		}
	}
	else
	{
		t = ( struct _tree * )malloc( sizeof( struct _tree ) ) ;
		t->b = b ;
		t->start = start ;
		t->end = end ;
		//t->keep = keep ;
		t->left = t->right = NULL ;
		t->list = list ;
		t->list->keep = keep ;

		if ( leftChild )
			parent->left = t ;
		else
			parent->right = t ;
	}
}

void OutputTree( struct _tree *t ) 
{
	if ( t == NULL )
		return ;
	//printf( "Remaining\n" ) ;	
	//if ( t->keep )
		samwrite( fpout, t->b ) ;

	OutputTree( t->left ) ;
	OutputTree( t->right ) ;
}

struct _list * OutputList( struct _list *p, int cut, bam1_t *b )
{
	struct _list *q ;
	while ( p != NULL && ( cut == -1 || p->end < cut || p->b->core.tid != b->core.tid ) )
	{
		if ( p->keep )
			samwrite( fpout, p->b ) ;
		bam_destroy1( p->b ) ;
		q = p ;
		p = p->next ;
		free( q ) ; // b is freed by RemoveNode.
	}
	return p ;
}

int main( int argc, char *argv[] ) 
{
	FILE *fp ;
 	samfile_t *fpsam ;
	bam1_t *b = NULL ;
	struct _tree root ;
	bool useSam = true, hasSplice ;
	char readId[LINE_SIZE] ;
	char cigar[LINE_SIZE] ;
 	int i, len, num ;
 	int start,end ;
	char prevChrome[103] ;
	struct _list *head, *tail, *node ;
	strcpy( prevChrome, "" ) ;

	// processing the argument list
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-h" ) )
		{
			PrintHelp() ;
			return 0 ;
		}
	}
	if ( argc == 1 )
	{
		PrintHelp() ;
		return 0 ;
	}

	len = strlen( argv[1] ) ;
	if ( !( fpsam = samopen( argv[1], "rb", 0 ) ) )
	{
		printf( "Can not open file %s\n", argv[1] ) ;
		return 0 ;
	}

	root.left = NULL ;
	qHead = qTail = qsHead = qsTail = 0 ;
	head = ( struct _list *)malloc( sizeof( *head ) ) ;
	tail = head ;
	head->next = NULL ;

	fpout = samopen( "-", "wb", fpsam->header ) ;
	while ( 1 )
	{
		int flag = 0 ;
		//if ( b )
		//	bam_destroy1( b ) ;
		b = bam_init1() ;
		if ( samread( fpsam, b ) <= 0 )
			break ;
		cigar2string( &(b->core), bam1_cigar( b ), cigar ) ;
		strcpy( readId, bam1_qname( b ) ) ;	
		flag = b->core.flag ;	
		if ( bam_aux_get( b, "NH" ) )
		{	
			if ( bam_aux2i( bam_aux_get(b, "NH" ) ) > 2 )
				secondary = true ;
			else
				secondary = false ;
		}
		else
			secondary = false ;


		mateStart = b->core.mpos + 1 ;

		/*if ( flag & 0x100 )
		  secondary = true ;
		  else 
		  secondary = false ;*/
		if ( bam_aux_get( b, "XS" ) )
			strand = bam_aux2A( bam_aux_get( b, "XS" ) ) ;
		else
			strand = '+' ;
		start = b->core.pos + 1 ;
		
		hasSplice = false ;
		len = num = 0 ;
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
				len += num ;
				num = 0 ;
				hasSplice = true ;
			}
			else
			{
				len += num ;
				num = 0 ;
			}
		}
		end = start + len - 1 ;
		//printf( "# %d %d\n", start, end ) ;
		node = ( struct _list *)malloc( sizeof( *node ) ) ;
		node->b = b ;
		node->start = start ;
		node->end = end ;
		node->next = NULL ;

		tail->next = node ;
		tail = node ;
		//printf( "%d\n", head->next->start ) ;
		head->next = OutputList( head->next, start, b ) ;
		SearchTree( b, start, end, true, root.left, &root, true, node ) ;
	}
	//OutputTree( root.left ) ;
	OutputList( head->next, -1, b ) ;
	samclose( fpsam ) ;
	samclose( fpout ) ;
	//fprintf( stderr, "The number of junctions: %d\n", junctionCnt ) ;
	return 0 ;
}

