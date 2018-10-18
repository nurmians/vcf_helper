#!/usr/bin/python2.7

"""

A tool for sorting VCF (Variant Call File) into numerical order based on the contig
numbering. This tool uses the UNIX sort command. Input can be provided from a file
or STDIN and the sorted output is directed to STDOUT. Contigs otside of chomosomes
1-22, X and Y will be sorted by coordinate and added to the end of the file in the
order they are encountered in the input.

Usage:
  sort_vcf [options] [<file>]  

Examples:
  sort_vcf file.vcf  
  sort_vcf file.vcf.gz
  gzip -dc file.vcf.gz | sort_vcf -s -n > file_sorted.vcf
  

Options:
  -c --chr-prefix      Use chr prefix for standard contigs (chr1-chr22,X,Y) in 
                       the output [default].
  -n --no-prefix       Output standard contigs without the chr prefix.
  -t --tmpdir=DIR      Directory for storing temporary files
  -f --filter=STR      Discard data rows without STR
  -s --standard        Output only standard chomosomes (chr1-chr22,X,Y)

"""

import docopt
import subprocess, sys, os, datetime, time, re
import string
import random
import gzip

# Order to sort contigs by
STD_CONTIGS = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

def error( aMsg):
  sys.stderr.write( "ERROR: "+aMsg+"\n" )
  sys.exit(-1)

def warning( aMsg):
  sys.stderr.write( "WARNING: "+aMsg+"\n" )

def info( aMsg):
  sys.stderr.write( "INFO: "+aMsg+"\n" )

CONTIG_PATTTERN = re.compile("(##contig=<ID=)(?:chr)*(.*?)([,>].*$)")

def EditContig( aLine, aChrPrefix=True, aOnlyStandard=True):

    if not aLine.startswith( "##contig"): return aLine

    m = CONTIG_PATTTERN.search( aLine)
    if m == None: return aLine
    is_standard = m.group( 2) in STD_CONTIGS
    if aOnlyStandard and not is_standard: return ""

    return m.group( 1) + ("chr" if (aChrPrefix and is_standard) else "") + m.group( 2) + m.group( 3) + "\n"


if __name__ == '__main__':

    args = docopt.docopt(__doc__)

    max_workers = 3

    prefix = "chr"    
    if bool(args['--no-prefix']): prefix = ""
    only_standard = False
    if bool(args['--standard']): only_standard = True
    filename = args['<file>']    
    filter_str = args['--filter']    

    if (not filename or len( filename) == 0):
      if sys.stdin.isatty(): error('No filename or input provided.')
      input_handle = sys.stdin
    else:

      try:
        if filename.endswith(".gz"): input_handle = gzip.open( aFile,'r')
        else: input_handle = open( filename, "r")
      except:
        error("Could not open file '%s'" % aFile)
  

    tmp_dir = args['--tmpdir']
    if (not tmp_dir or len( tmp_dir) == 0): tmp_dir = "."
    if not tmp_dir.endswith("/"): tmp_dir = tmp_dir + "/"


    tmp_files = {}
    non_standard_contigs = []
    # Insert a random str into the temp filenames in case sortvcf is run in parallel
    # for multiple files in the same tmp folder
    while True:
      random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range( 8))
      tmpfilename_template = tmp_dir + "contig_%s_" + random_str + ".tmp"
      if not os.path.isfile( tmpfilename_template % "1"): break #Make sure it does not exist

    header = []

    # Process input
    for line in input_handle:
      cols = line.strip().split("\t")
      contig = ""

      if cols[ 0].startswith("#"):
        # Reformat header contig lines
        line = EditContig( line, len(prefix)>0, only_standard)
        if len( line) > 0: header.append( line)
        continue
      else:
        if cols[ 0].startswith("chr"): 
          contig = cols[ 0][3:]
          line = line[3:] # TODO: better ways to do this
        else: contig = cols[ 0]

        if len( cols) < 2: 
          #warning( "Skipping line (no cols): '%s'" % line)
          continue 

        if len( contig) == 0: 
          warning( "Skipping line: '%s'" % line)
          continue       

        if contig not in STD_CONTIGS and contig not in non_standard_contigs:
          if only_standard: continue 
          non_standard_contigs.append( contig)  


      tmpfilename = tmpfilename_template % contig
  
      # Create tmp file
      if tmpfilename not in tmp_files:
        try: tmp_files[ tmpfilename] = open( tmpfilename, "w")
        except: error("Could not create file '%s'." % tmpfilename)

      #Write to tmp file
      if filter_str and len( filter_str):
        if line.find( filter_str) < 0: continue

      if contig in STD_CONTIGS:
        # Use prefix for standard contigs (if specified)
        tmp_files[ tmpfilename].write( prefix + line)
      else:
        # No prefix
        tmp_files[ tmpfilename].write( line)

    # Input processed
    # Close tmp files
    for fn, fh in tmp_files.items():
      try: fh.close()
      except: warning( "Could not close file: '%s'." % fn)

    # Sort tmp files based on coordinate
    # Use UNIX sort command
    for fn, fh in tmp_files.items():
      #Sort by chromosomal coordinate
      #worker_cmd = "sort -t $'\t' -g -k 2 -o %s %s" % ((fn+".sorted"), fn)
      worker_cmd = "sort -t '\t' -g -k 2 -o %s %s" % ((fn+".sorted"), fn)
      #print "CMD:", worker_cmd
      worker = subprocess.Popen(worker_cmd, shell=True)
      #workers = [subprocess.Popen(worker_cmd) for w in range(max_workers)]
      #for w in workers: w.wait()
      worker.wait()

      try: os.remove( fn)
      except: warning("Could not remove tmp file '%s'" % fn)        

    # Print header
    htf = tmpfilename_template % "header"
    if len( header) == 0: warning( "File has no header.")
    else: sys.stdout.write( "".join( header))

    # Print data rows
    for cont in STD_CONTIGS+non_standard_contigs:
      tf = tmpfilename_template % cont
      tfs = (tmpfilename_template % cont) + ".sorted"
      if os.path.isfile( tfs):
        with open( tfs, "r") as tfsh:
          for line in tfsh:
            sys.stdout.write( line)

        try: os.remove( tfs)
        except: warning("Could not remove tmp file '%s'" % tfs)
      else:
        if not cont in STD_CONTIGS: warning("File '%s' does not exits." % tfs)


    sys.exit( 0)

