#!/usr/bin/python2.7

"""

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


"""

import docopt
import sys, os, datetime, time, re
import string
import random
import gzip
import atexit

import multiprocessing
import subprocess
import shlex
from multiprocessing.pool import ThreadPool

from sys import platform



# Order to sort contigs by
STD_CONTIGS = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
STD_CONTIGS_DICT = { c : True for c in STD_CONTIGS } 

def error( aMsg):
  sys.stderr.write( "ERROR: "+aMsg+"\n" )
  UNEXPECTED_EXIT = False
  sys.exit(-1)

def warning( aMsg):
  sys.stderr.write( "WARNING: "+aMsg+"\n" )

def info( aMsg):
  sys.stderr.write( "INFO: "+aMsg+"\n" )

CONTIG_PATTTERN = re.compile("(##contig=<ID=)(?:chr)*(.*?)([,>].*$)")
UNEXPECTED_EXIT = True
REMOVABLE_TMP_FILES =[]

if platform != "linux" and platform != "linux2":
    error("This program requires the Linux sort algorithm.")
    sys.exit( -1)

def ExitHandler():

    #global UNEXPECTED_EXIT, REMOVABLE_TMP_FILES
    if UNEXPECTED_EXIT:
        try: 
          warning( "Unexpected exit.")
          n_tmpfiles = len( REMOVABLE_TMP_FILES)
          if n_tmpfiles > 0: info( "Removing remaining %i temporary file%s." % ( n_tmpfiles, ("" if n_tmpfiles == 1 else "s")))
        except: pass

        for tmpfile in REMOVABLE_TMP_FILES:
          try: os.remove( tmpfile)
          except: warning("Could not remove temporary file '%s'." % tmpfile)   


def EditContig( aLine, aChrPrefix=True, aOnlyStandard=True):

    if not aLine.startswith( "##contig"): return aLine

    m = CONTIG_PATTTERN.search( aLine)
    if m == None: return aLine
    is_standard = m.group( 2) in STD_CONTIGS
    if aOnlyStandard and not is_standard: return ""

    return m.group( 1) + ("chr" if (aChrPrefix and is_standard) else "") + m.group( 2) + m.group( 3) + "\n"



BED = {}

def ReadBed( aFilename):

  n_contigs = 0
  n_ranges = 0
  try:
    with open( aFilename) as f:
      for line in f:
        if line.startswith("#"): continue #Comments
        if line.startswith("CHROM"): continue #Header
        cols = line.strip().split("\t")
  
        chromosome = cols[ 0]
        start = int(cols[ 1])
        end = int(cols[ 2])
        if chromosome.startswith("chr") and chromosome[3:] in STD_CONTIGS: chromosome = chromosome[3:]
  
        if chromosome not in BED: 
          BED[ chromosome] = []
          n_contigs += 1
        BED[ chromosome].append([start,end])
        n_ranges += 1
  except Exception as ex:
    error( "Failed to read bed file '%s'. (%s)\n" % (aFilename, str(ex)))

  if len(BED) < 1: error( "Bedfile '%s' contained %i contigs and %i ranges." % (aFilename, n_contigs, n_ranges))
  info( "Bedfile '%s' contained %i contigs and %i ranges." % (aFilename, n_contigs, n_ranges)) 

  # Sort ranges for each chromosome
  for k in BED.keys():
    BED[ k].sort(key=lambda x: x[ 0])



def IsBedUsed():
  return len( BED) > 0

BED_WARNED = {}

def IsInBedRange( aChromosome, aCoordinate):

  if len( BED) == 0: return True

  if aChromosome.startswith("chr") and aChromosome[3:] in STD_CONTIGS: aChromosome = aChromosome[3:]
  if aChromosome not in BED:
    if aChromosome not in BED_WARNED:
      warning("Contig '%s' not in bed file." % aChromosome)
      BED_WARNED[ aChromosome] = True
    return False
  for r in BED[aChromosome]:
    if r[ 0] > aCoordinate: break #start
    if aCoordinate <= r[ 1]: return True #end
  return False




def call_proc( cmd):
    """ This runs in a separate thread. """
    return subprocess.call( shlex.split( cmd))  # This will block until cmd finishes
    #p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #out, err = p.communicate()
    #return (out, err)

if __name__ == '__main__':

    args = docopt.docopt(__doc__)
    
    atexit.register( ExitHandler)

    prefix = "chr"    
    if bool(args['--no-prefix']): prefix = ""
    only_standard = False
    if bool(args['--standard']): only_standard = True
    filename = args['<file>']    
    filter_str = args['--filter']    

    disable_sort = False
    if bool(args['--disable-sort']): disable_sort = True

    indels_only = False
    if bool(args['--indels-only']): indels_only = True

    verbose = False
    if bool(args['--verbose']): verbose = True

    keep_duplicates = False
    if bool(args['--keep-duplicates']): keep_duplicates = True

    bedfile = args['--bed']
    if bedfile != None and len( bedfile) > 0: ReadBed( bedfile)

    remove_comments = False
    remove_comments = args['--remove']

    try:
      n_cpus = args['--processes']
      if n_cpus == None: n_cpus = 1
      else: n_cpus = int( n_cpus)
      if n_cpus < 1: n_cpus = 1
    except:
      warning( "Invalid number of processes specified '%s'." % n_cpus)
      n_cpus = 1


    if (not filename or len( filename) == 0):
      if sys.stdin.isatty(): error('No filename or input provided.')
      input_handle = sys.stdin
    else:

      try:
        if filename.endswith(".gz"): input_handle = gzip.open( filename,'r')
        else: input_handle = open( filename, "r")
      except Exception as ex:
        error("Could not open file '%s'. (%s)" % (filename,str(ex)))
  

    tmp_dir = args['--tmpdir']
    if (not tmp_dir or len( tmp_dir) == 0): tmp_dir = "."
    if not tmp_dir.endswith("/"): tmp_dir = tmp_dir + "/"


    extra_contigs = args['--add-contig']
    if extra_contigs and len( extra_contigs) > 0: 
      STD_CONTIGS += map(str.strip, extra_contigs.split( ","))

    tmp_files = {}
    non_standard_contigs = []
    # Insert a random str into the temp filenames in case sortvcf is run in parallel
    # for multiple files in the same tmp folder
    while True:
      random_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range( 8))
      tmpfilename_template = tmp_dir + "contig_%s_" + random_str + ".tmp"
      if not os.path.isfile( tmpfilename_template % "1"): break #Make sure it does not exist

    header = []
    is_dup = False
    n_duplicates = 0
    prev = "?"
    line_num = 0
    in_bed_range = 0    
    in_bed_checked = 0
    header_ended = False
    comments_warned = False

    # Process input
    # Sort lines into separate files based on the contig
    if verbose: info( "Splitting input file by contig...")

    for line in input_handle:
      line_num += 1
      cols = line.strip().split("\t")
      contig = ""

      if cols[ 0].startswith("#"):
        # Reformat header contig lines
        if header_ended and not comments_warned and not remove_comments: 
          warning("Encountered comment line outside of header on line %i." % line_num)
          comments_warned = True
        line = EditContig( line, len(prefix)>0, only_standard)
        if len( line) > 0: header.append( line)
        continue
      else:        
        if cols[ 0].startswith("chr") and cols[ 0][3:] in STD_CONTIGS_DICT: 
          contig = cols[ 0][3:]
          line = line[3:] # TODO: better ways to do this
        else: contig = cols[ 0]

        if len( cols) < 2: 
          #warning( "Skipping line (no cols): '%s'" % line)
          continue 

        if len( contig) == 0: 
          warning( "Skipping line %i: '%s'" % (line_num, line[:50]))
          continue       

        header_ended = True

        # Standard contig check
        if contig not in STD_CONTIGS_DICT and contig not in non_standard_contigs:
          if only_standard: continue 
          non_standard_contigs.append( contig)  

        # Duplicate line removal
        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
        if not keep_duplicates:
          cur = "".join(cols[:5]) # join up to ALT
          if cur == prev:
            n_duplicates += 1
            if verbose or n_duplicates <= 10:
              warning( "Skipping duplicate line %i: '%s'" % (line_num, line[:50]))
              if not verbose and n_duplicates == 10: warning("Only first 10 duplicate lines warned.")
            continue
          prev = cur        

        # BED file filtering
        if IsBedUsed():
          in_bed_checked += 1
          try: int_coord = int( cols[ 1])            
          except: error("Could not convert coordinate '%s' to integer on line %i." % (cols[ 1], line_num))
          if not IsInBedRange( contig, int_coord): continue
          in_bed_range += 1

      tmpfilename = tmpfilename_template % contig
  
      # Create tmp file
      if tmpfilename not in tmp_files:
        try: 
          tmp_files[ tmpfilename] = open( tmpfilename, "w")
          REMOVABLE_TMP_FILES.append( tmpfilename)
        except: 
          error("Could not create file '%s'." % tmpfilename)
      
      # Filter lines
      if filter_str and len( filter_str):
        if line.find( filter_str) < 0: continue

      if indels_only:
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
        line_cols = line.split( "\t", 5)
        if len( line_cols[ 3]) == 1 and len( line_cols[ 4]) == 1: continue        

      #Write to tmp file
      if contig in STD_CONTIGS:
        # Use prefix for standard contigs (if specified)
        tmp_files[ tmpfilename].write( prefix + line)
      else:
        # No prefix
        tmp_files[ tmpfilename].write( line)

    if n_duplicates > 0: info( "Skipped %i duplicate line%s." % (n_duplicates, "" if n_duplicates == 1 else "s"))

    # Input processed
    # Close tmp files
    for fn, fh in tmp_files.items():
      try: fh.close()
      except: warning( "Could not close file: '%s'." % fn)

    if verbose: info( "Input file split.")


    if not disable_sort:

      if verbose: info( "Sorting each contig.")

      # Sort tmp files based on coordinate
      # Use UNIX sort command
      pool = ThreadPool( n_cpus)

      for fn, fh in tmp_files.items():
        #Sort by chromosomal coordinate
        #worker_cmd = "sort -t $'\t' -g -k 2 -o %s %s" % ((fn+".sorted"), fn)
        worker_cmd = "sort -t '\t' -T %s -g -k 2 -o %s %s" % (tmp_dir, (fn+".sorted"), fn)
        pool.apply_async( call_proc, (worker_cmd,))

        REMOVABLE_TMP_FILES.append( fn+".sorted")

        #print "CMD:", worker_cmd
        #worker = subprocess.Popen(worker_cmd, shell=True)
        #workers = [subprocess.Popen(worker_cmd) for w in range(max_workers)]
        #for w in workers: w.wait()
        #worker.wait()

      # Close the pool and wait for each running task to complete
      pool.close()
      pool.join()

      # Remove temporary files
      for fn, fh in tmp_files.items():
        try: 
          os.remove( fn)
          REMOVABLE_TMP_FILES.remove( fn)          
        except: 
          warning("Could not remove tmp file '%s'" % fn)   

    else:
      # Sorting disabled, only rename files as if already sorted
      if verbose: info( "Skipping sort.")

      for fn, fh in tmp_files.items():
        try: 
          os.rename( fn, fn+".sorted")
          REMOVABLE_TMP_FILES.remove( fn)          
          REMOVABLE_TMP_FILES.append( fn+".sorted")
        except: warning("Could not rename tmp file '%s'" % fn)   

    # Print header
    if not remove_comments:
      if verbose: info( "Printing header")
  
      htf = tmpfilename_template % "header"
      if len( header) == 0: warning( "File has no header.")
      else: sys.stdout.write( "".join( header))
  
    # Print data rows
    if verbose: info( "Printing each contig...")

    for cont in STD_CONTIGS+non_standard_contigs:
      
      tf = tmpfilename_template % cont
      tfs = (tmpfilename_template % cont) + ".sorted"
      if os.path.isfile( tfs):

        if verbose: info( "Contig: %s" % cont)

        with open( tfs, "r") as tfsh:
          for line in tfsh:
            sys.stdout.write( line)

        try: 
          os.remove( tfs)
          REMOVABLE_TMP_FILES.remove( tfs)  
        except: warning("Could not remove tmp file '%s'" % tfs)
      else:
        if not cont in STD_CONTIGS: warning("File '%s' does not exits." % tfs)


    if IsBedUsed():
      info( "%.1f%% of sorted calls found in bed range." % (float(in_bed_range)/in_bed_checked*100.0))

    UNEXPECTED_EXIT = False
    if verbose: info( "All Done.")
    sys.exit( 0)

