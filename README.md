# VCF helper
Python scripts for DNA sequencing data VCF (Variant Call File) manipulation, sorting and comparison.


# sort_vcf.py


A tool for sorting VCF (Variant Call File) into numerical order based on the contig
numbering. This tool uses the UNIX sort command. Input can be provided from a file
or STDIN and the sorted output is directed to STDOUT. Contigs otside of chomosomes
1-22, X and Y will be sorted by coordinate and added to the end of the file in the
order they are encountered in the input.

sort_vcf can accept gzipped (and bgzipped) input, but does not output gzipped data,
because VCF files are required to be zipped with bgzip (block gzip), which can only
process complete files (not pipes or streams).

Usage:
  sort_vcf [options] [<file>]

Examples:
  sort_vcf file.vcf
  sort_vcf file.vcf.gz > file_sorted.vcf
  gzip -dc file.vcf.gz | sort_vcf -s -n > file_sorted.vcf

Options:
  -c --chr-prefix      Use chr prefix for standard contigs (chr1-chr22,X,Y) in
                       the output [default].
  -n --no-prefix       Output standard contigs without the chr prefix.
  -t --tmpdir=DIR      Directory for storing temporary files
  -f --filter=STR      Discard data rows without STR
  -s --standard        Output only standard chomosomes (chr1-chr22,X,Y)
  -a --add-contig=CON  Add contigs to std set of contigs, separate by commas
  -p --processes=N     Maximum number of worker processes to use [default:1]
  -d --disable-sort    Do not sort, only filter (-f) and remove contigs (-s)
  -v --verbose         Report progress to stderr
  -i --indels-only     Output only indels


# comp_vcf.py

A tool for comparing VCFs (Variant Call Files) and reporting common and unique
called variants in each file. The VCF files must be sorted to the exact same
chromosomal and coordinate order. If only one input file is specified comp_vcf
will reports statistics on that file.

The output contains two matrices: (i) Comparison matrix shows the number of
matching variant calls (rows) between eavery two files included in the
comparison. (ii) Count matrix shows for each file in how many of the other
files in the comaprison each call is found in. Calls shown on row "IN 1"
are unique to that file and not found in any of the other files in the
comparison. Calls found in one other file are shown on row "IN 2", and so on.

The non-standard contigs in the files in the comparison need to be in the
same order in each file. Using the flag -s ignores all calls in non-standard
contigs to eliminate this problem.

Usage:
  comp_vcf [options] <file>...

Examples:
  comp_vcf file.vcf
  comp_vcf --names file1,file2,file3 ~/file1.vcf /data/file2.vcf ~/file3.vcf


Options:
  -r --report          Include statistics of individual files in output
  -s --standard        Ignore non-standard chromosomes (chr1-22,X,Y)
  -n --names=NAMES     Column name for each input file (separate by commas)
  -b --baseratios      Calculate and print (with --report) base substitution
                       ratios
  -i --ignore_indels   Skip indels
  -c --coordinates     Do comparison basen on contig and coordinate only
  
  
