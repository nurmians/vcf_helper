# VCF tools
Python scripts for DNA sequencing data VCF (Variant Call File) manipulation, sorting and comparison.

# sort_vcf.py


A tool for sorting VCF (Variant Call File) data based on contig and coordinate 
numbering by using the UNIX sort command. Input can be provided from a file or 
STDIN and the sorted output is directed to STDOUT. Contigs outside of chomosomes
1-22, X and Y will be sorted by coordinate and added to the end of the file in the
order they are encountered in the input.

sort_vcf can accept gzipped (and bgzipped) input, but does not output gzipped data,
because VCF files are required to be zipped with bgzip (block gzip), which can only
process complete files (not pipes or streams).

Another useful feature of this script is the ability to remove all special contigs
(-s) and to make contig naming consistent between different files ("chr1" vs. "1").

<pre>
Usage:
  sort_vcf [options] [&lt;file&gt;]

Examples:
  sort_vcf file.vcf
  sort_vcf file.vcf.gz > file_sorted.vcf
  gzip -dc file.vcf.gz | sort_vcf -s -n > file_sorted.vcf

Options:
  -c --chr-prefix       Use chr prefix for standard contigs (chr1-chr22,X,Y) in 
                        the output [default].
  -n --no-prefix        Output standard contigs without the chr prefix.
  -t --tmpdir=DIR       Directory for storing temporary files
  -f --filter=STR       Discard data rows without STR
  -s --standard         Output only standard chomosomes (chr1-chr22,X,Y)
  -a --add-contig=CON   Add contigs to std set of contigs, separate by commas      
  -p --processes=N      Maximum number of worker processes to use [default:1]
  -d --disable-sort     Do not sort, only filter (-f) and remove contigs (-s)
  -v --verbose          More output
  -i --indels-only      Output only indels
  -k --keep-duplicates  Do not remove duplicates (same chrom, coord & alt)
  -b --bed=FILE         Remove calls outside of specified ranges
  -r --remove           Remove all header and comment lines
</pre>

# comp_vcf.py

A tool for comparing VCFs (Variant Call Files) and reporting common and unique
variants called in each file. The VCF files must be sorted to the exact same
chromosomal and coordinate order. If only one input file is specified comp_vcf
will reports statistics on that file. 

The output contains two matrices: 
(i) Comparison matrix shows the number of matching variant calls between every 
pair of files included in the comparison. 
(ii) Count matrix shows for each file, in how many other 
files in the comparison each call is found in. Calls shown on row "IN 1"
are unique to that file and not found in any of the other files in the
comparison. Calls found in one other file are shown on row "IN 2", and so on.

Standart contig names with "chr" prefix can be compared with contigs without
the prefix, e.g. contig "1" matches "chr1".

<h5>Compvcf unpacks multiple ALTs for comparisons (if -c not specified)</h5> 
<pre>"C > A,GT,ATG" => "C > A", "C > GT", "C > ATG"</pre>
<h5>Compvcf unpacks MNVs for comparisons</h5>
<pre>"CT > GA" => "C > G", "T > A"</pre>

<pre>
Usage:
  comp_vcf [options] &lt;file&gt;...

Examples:
  comp_vcf file.vcf  
  comp_vcf --names file1,file2,file3 ~/file1.vcf /data/file2.vcf ~/file3.vcf
  

Options:
  -r --report          Include statistics of individual files in output
  -s --standard        Ignore non-standard chromosomes (chr1-22,X,Y)
  -n --names=NAMES     Column name for each input file (separate by commas)
  -e --baseratios      Calculate and print (with --report) base substitution 
                       ratios.
  -b --bed=FILE        Bed file for specifying genomic intervals. If a 
                       bed-file is provided, everything outside the bed
                       ranges is ignored.
  -g --ignore=FILE     File specifying ranges to ignore in the comparison
                       (Bed or VCF file format).
  -i --ignore-indels   Skip (ignore) indels.
  -I --ignore-snvs     Skip (ignore) SNVs (single nucleotide variants).
  -c --coordinates     Do comparison based on contig and coordinate only
  -d --different       Output rows of the first input file that are not 
                       found in the second file.
  -D --all-different   Output all unmatched rows.
  -m --matching        Output rows of the first input file that are 
                       found in all other input files.
  -M --all-matching    Output all matching rows from all input files.
                       May create duplicate (identical) rows.
  -w --swap            Swap the first and second input files.
  -f --filter=STR      Ignore data rows without STR.
  -h --header          Include header row from first input file to
                       output.
  -a --add-id=STR      Add STR to ID column of each outputted data row.
  -A --add-info        Generate and insert at the beginning of each row
                       a string indicating (Y/N) in which input files 
                       the variant (identical row) is present.
  -p --pretty          Cut outputted lines to 80 chars.
  -P --no-prefix       Output contigs without "chr" prefix. Default is with
                       prefix.
  -o --outform=FORM    Include characters to display output. T=total SNVs,
                       M=comparison matrix, C=common counts, P=Precision 
                       and recall, A=All (default:All but P).
  -t --tabix=CONTIG    Process only contig CONTIG. Requires tabix and .tbi
                       file.
  -q --quiet           Suppress warnings.
  -R --skip-ref        Skip rows where no ALT has been specified.

  </pre>
  
<h3>Comparing two files</h3>

<pre>
comp_vcf 1m.vcf 1mf.vcf

##This list shows the number of total SNV calls in each file (T)
Total number of calls in each file:
1m.vcf: 3
1mf.vcf: 2

##This matrix shows the number of equal calls between two files (M)
Comparison matrix:
        1m.vcf  1mf.vcf
1m.vcf: 3       2
1mf.vcf: 2      2

##This matrix shows in how many files each call made was found (C)
Count matrix:
        1m.vcf  1mf.vcf
IN 1: 1 0
IN 2: 2 2
</pre>

  
Comparing and outputting matching calls in chr17:
<pre>
comp_vcf -t chr17 -m g2.vcf.gz g2f.vcf.gz
</pre>
  
Comparing and outputting calls found only in '1m.vcf':
<pre>
comp_vcf -d -w 1mf.vcf 1m.vcf
</pre>
  
  
  
  
