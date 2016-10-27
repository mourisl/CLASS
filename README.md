CLASS - Constraint-based Local Assembly and Selection of Splice variants
========================================================================

Described in: 
 
Current (dynamic-programming) version:
Song, L., Sabunciyan, S., and Florea, L. (2016). CLASS2: accurate and efficient splice variant annotation from RNA-seq reads. Nucleic Acids Research 2016; doi: 10.1093/nar/gkw158

Previous ('set cover') version:
Song, L. and Florea, L. (2013). CLASS: Constrained Transcript Assembly of RNA-seq Reads.
Third Annual RECOMB Satellite Workshop on Massively Parallel Sequencing - RECOMB-SEQ 2013. BMC Bioinformatics 14(Suppl. 5), S14.

Copyright (C) 2012-2013, and GNU GPL, by Li Song, Liliana Florea

Includes portions copyright from:

lp_solve - Copyright (C) 2005, and GNU LGPL, by Michel Berkelaar, Kjell Eikland, Peter Notebaert
SAMtools - Copyright (C) 2008-2009, Genome Research Ltd, Heng Li

### Content:
I.   What is CLASS?
II.  Install
III.  Usage
IV. Input/Output
V.  Example
VI.   Terms of use
VII.  Support

### I.   What is CLASS?

CLASS is a tool for assembling transcripts from short RNA-seq reads
aligned to a reference genome. It works in three stages. Stage 1 makes
use of linear programming to determine a set of exons.  Stage 2 builds a
splice graph representation of a gene, by connecting the exons (vertices)
via introns (edges) extracted from spliced read alignments. Stage 3
selects a subset of the candidate transcripts encoded in the graph,
according to the constraints derived from mate pairs and spliced alignments
and, optionally, using knowledge about gene structure extracted from known annotation and/or alignments of cDNA sequences.

### II. Install
Run: % sh build.sh

The software is compiled in place

### III.  Usage
	Usage: perl run_class.pl [options]
	Options:
		-a alignment_file (REQUIRED): the path to the alignemtn file(in SAM or BAM format)
		-o output_file: the file storing the output of CLASS (default: ./alignment_file_wo_extension.gtf)
		-p number_of_threads: specify the number of worker threads (default:1)
		-F f: do not report the transcripts whose abundance level is lower than f*|most expressed transcript| in a gene
		-l label: add a prefix and a "_" to the ids in the GTF file (default: not used)
		-j junction: the path to the splice junction file
		-e evidence: the path to the evidence files
		--var_rd_len: extensive variable read lengths, i.e. reads after trimming (default: no)
		--set-cover: use set cover to build transcripts from splicing graph (default: no)
		--verbose: also output the procedure of CLASS (default: no)
		--wd tempoaray_file_directory: the directory storing the temporary files (default: ./class_tmp)
		--clean: whehter to remove the temporary files in -wd (default: no)


### IV. Input/Output

The primary input to CLASS is a set of short read alignments in BAM format
and sorted by chromosome and position, for instance one produced with
the program Tophat2 (http://tophat.cbcb.umd.edu).  

Given an alignment input x.bam, CLASS produces two intermediate data files,
x.depth and x.splice in the temporary working directory.

  * The format of the x.depth file, generate by samtools, is:
    chrom_id position #_of_reads_on_the_position

  * The format of the x.splice file, generate by 'junc', is:
    chrom_id start_intron_position end_intron_position #_of_supporting_reads strand\ 
     #_unique_support_read #_multi-aligned_support_read\ 
     sum_of_unique_support_reads_edit_distance\
     sum_of_multi-aligned_support_reads_edit_distance

    NOTE: When using the '-a' argument, the value #_of_supporting_reads can be
    negative, indicating that this splice junction is invalid.

Lastly, to produce a set of transcripts, the program 'class' takes as
input a BAM/SAM file, the depth file generated by 'samtools' and the splice
junction file generated by 'junc'. The final output, consisting of predicted
transcripts, is in standard GTF format.

NOTE: 'samtools' can only read a BAM file, while 'junc' and 'class' can
read both SAM and BAM formats.

IMPORTANT:
Please make sure to enable the appropriate access to the directory
containing the temporary workind girectory. Also, since the depth file will be created in
that directory and will have one row for each base, please ensure that
sufficient space is available in the work directory.

### V. Example

We will use the file ./Sample/sample.bam as an example.
Running "% perl run_class.pl -a ./Sample/sample.bam" will produce a 
set of transcripts in the file ./sample.gtf . 
 
### VI. Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### VII. Support

Contact us at: lsong10@jhu.edu, florea@jhu.edu

