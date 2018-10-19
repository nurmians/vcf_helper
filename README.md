# VCF helper
Python scripts for DNA sequencing data VCF (Variant Call File) manipulation, sorting and comparison.


# sort_vcf.py

A tool for sorting VCF (Variant Call File) into numerical order based on the contig
numbering. This tool uses the UNIX sort command. Input can be provided from a file
or STDIN and the sorted output is directed to STDOUT. Contigs otside of chomosomes
1-22, X and Y will be sorted by coordinate and added to the end of the file in the
order they are encountered in the input.

Examples:
  sort_vcf.py file.vcf  
  sort_vcf.py file.vcf.gz
  gzip -dc file.vcf.gz | sort_vcf.py -s -n > file_sorted.vcf


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

Examples:
  comp_vcf.py file.vcf  
  comp_vcf.py --names file1,file2,file3 ~/file1.vcf /data/file2.vcf ~/file3.vcf
  
  
