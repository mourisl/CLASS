/**
  Use an 64bit integer to represent bit table.

  Li Song
  July 24, 2012
*/


#ifndef _BIT_TABLE_HEADER
#define _BIT_TABLE_HEADER

#include <stdio.h>

typedef unsigned long long int INT64 ;
#define UNIT_SIZE (sizeof( INT64 ) * 8 )
#define UNIT_MASK ((INT64)63) 

class BitTable
{
private:
        INT64 *tab ; // The bits
        int size ; // The size of the table
	int asize ; // The size of the array (size/64).
public:
        BitTable() ; // Intialize a bit table.
        BitTable( int s ) ; // Initalize a bit table with size s.
        ~BitTable() ;
	
	void Init( int s ) ; // Initialize a bit table with size s.
	void Reset() ; // Make every value 0.

        void Set( int i ) ; // Set the ith bit.
        void Unset( int i ) ; // Unset the ith bit
        void Flip( int i ) ; // Flip the ith bit. Same as xor 1.        
        bool Test( int i ) ; // Test the ith bit.

	void Not() ; // Not all the bits.
	void And( const BitTable &in ) ; // Do the "and" on each bits.
	void Or( const BitTable &in ) ; // Do the "or" on each bits.

	int Count() ; // Count how many 1.

        bool IsEqual( const BitTable &in ) ; // Test wether two bit tables equal.   
} ;


#endif
