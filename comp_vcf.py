#!/usr/bin/python2.7

"""


A tool for comparing VCFs (Variant Call Files) and reporting common and unique
called variants in each file. The VCF files must be sorted to the exact same
chromosomal and coordinate order. If only one input file is specified comp_vcf
will reports statistics on that file. 

The output contains two matrices: (i) Comparison matrix shows the number of 
matching variant calls (rows) between every pair of files included in the 
comparison. (ii) Count matrix shows for each file in how many of the other 
files in the comaprison each call is found in. Calls shown on row "IN 1"
are unique to that file and not found in any of the other files in the
comparison. Calls found in one other file are shown on row "IN 2", and so on.

The non-standard contigs in the files in the comparison need to be in the 
same order in each file. Using the flag -s ignores all calls in non-standard 
(special) contigs to eliminate this problem.

Standart contig names with "chr" prefix can be compared with contigs without
the prefix, e.g. contig "1" matches "chr1".

Usage:
  comp_vcf [options] <file>...

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

"""

import docopt
import subprocess, sys, os, datetime, time, re
#from subprocess import Popen, PIPE
import gzip
import itertools
import datetime


# Order to sort contigs by
STD_CONTIGS = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
STD_CONTIGS_DICT = {"1":True,"2":True,"3":True,"4":True,"5":True,"6":True,"7":True,"8":True,"9":True,"10":True,"11":True,"12":True,"13":True,"14":True,"15":True,"16":True,"17":True,"18":True,"19":True,"20":True,"21":True,"22":True,"X":True,"Y":True}

QUIET = False

def error( aMsg):
  sys.stderr.write( "ERROR: "+aMsg+"\n" )
  sys.exit(-1)

def warning( aMsg):
  sys.stderr.write( "WARNING: "+aMsg+"\n" )

def info( aMsg):
  sys.stderr.write( "INFO: "+aMsg+"\n" )

def debug( aMsg):
  sys.stderr.write( "DEBUG: "+aMsg+"\n" )

def progress( aMsg):
  sys.stderr.write( "PROGRESS: "+aMsg+"\n" )

#CONTIG_PATTTERN = re.compile("(##contig=<ID=)(?:chr)*(.*?)([,>].*$)")
STD_BASES = ["A","C","G","T"]
USE_CHR_PREFIX_IN_OUTPUT = True

# Ranges to include
BED = {}
# Ranges to ignore
IGN = {}

#def EditContig( aLine, aChrPrefix=True, aOnlyStandard=True):
#
#    if not aLine.startswith( "##contig"): return aLine
#
#    m = CONTIG_PATTTERN.search( aLine)
#    if m == None: return aLine
#    if aOnlyStandard and m.group( 2) not in STD_CONTIGS: return ""
#
#    return m.group( 1) + ("chr" if aChrPrefix  else "") + m.group( 2) + m.group( 3)
#

def PrintReport( aFile, aName=None, aPrintBaseRatios=False):

  subs = {  "A->C" :0,
            "A->G" :0,
            "A->T" :0,
            "C->A" :0,
            "C->G" :0,
            "C->T" :0,
            "G->A" :0,
            "G->C" :0,
            "G->T" :0,
            "T->A" :0,
            "T->C" :0,
            "T->G" :0 
  }

  stats = { "sub" :0,          
            "indel":0,
            "del"  :0,
            "ins"  :0,
            "in_bed_range":0,
            "not_in_bed_range":0,
            "n_ignored":0,
            "total":0
  }

  contigs = {}
  contig_arr = []

  ##CHROM  POS     ID      REF     ALT     
  try:
    if aFile.endswith(".gz"): file = gzip.open( aFile,'r')
    else: file = open( aFile, "r")

  except:
    error("Could not open file '%s'" % aFile)

  for line in file:
    if line.startswith("#"): continue

    cols = line.split("\t", 5)
    stats["total"] += 1
    contig = cols[ 0]
    if contig not in contigs:
      contigs[ contig] = 1
      contig_arr.append( contig)
    else:
      contigs[ contig] += 1

    if IsIgnUsed() and IsInBedRange( cols[ 0], cols[ 1], IGN): 
      n_ignored += 1
      continue

    if IsInBedRange( cols[ 0], cols[ 1], BED): stats["in_bed_range"] += 1
    else: stats["not_in_bed_range"] += 1
      #sys.stderr.write( "Not in range: %s %s\n" % (cols[ 0], cols[ 1]))

    ref = cols[ 3]
    alt = cols[ 4]
    if len( ref) == 1 and len( alt) == 1:
      stats["sub"] += 1
      if ref in STD_BASES and alt in STD_BASES:
        subs["%s->%s" % (ref, alt)] += 1
    else:
      stats["indel"] += 1
      if len( ref) > len( alt): stats["del"] += 1
      elif len( alt) > len( ref): stats["ins"] += 1

  try: file.close()
  except: pass

  if aName != None and aName != False:    
    print "\nStatistics for '%s':" % aName
  else:
    print "\nStatistics for file '%s':" % aFile

  print "n_Contigs: ", len( contig_arr)
  #print "Contigs: ", contig_arr
  print "Calls per contig: ", ["%s:%i" % (k,contigs[ k]) for k in contig_arr]

  for k in stats.keys():
    if not IsBedUsed() and k == "in_bed_range": continue
    if not IsBedUsed() and k == "not_in_bed_range": continue
    print "%s: %s" % (k, stats[ k])

  if IsBedUsed() and stats["in_bed_range"]+stats["not_in_bed_range"]>0:
    print "bed_in_range_percentage : %.1f%%" %  (float(stats["in_bed_range"])/(stats["in_bed_range"]+stats["not_in_bed_range"])*100.0)


  if aPrintBaseRatios:
    n_subs = int( stats["sub"])
    for k in sorted( subs.keys()):
      print "%s: %3.1f%%" % (k, int(subs[ k]) / float(n_subs)*100.0 )

# Combinations of indels and SNVs can cause coordinates
# to become out of order. Keep a sorted buffer of rows
# before outputting them 
OUTBUF=[]
BUFSIZE=25
BUF_MAX=[]
BUF_INDEX=-1

def buffered_init():
  global BUF_INDEX, OUTBUF, BUF_MAX
  BUF_INDEX += 1
  OUTBUF.append( list())
  BUF_MAX.append( 999999999999)
  return BUF_INDEX

def buffered_yield( bi, aOutCols):
  global OUTBUF, BUFSIZE, BUF_MAX
  lbuf=len(OUTBUF[bi])
  if lbuf > 0 and OUTBUF[bi][ 0][ 0] != aOutCols[ 0]: # New contig    
    buf = buffered_flush( bi)
    while True:
      out = buf.next()
      if out == None: break
      yield out    
    lbuf=0
  OUTBUF[bi].append( aOutCols)
  if aOutCols[ 1] < BUF_MAX[bi] and lbuf>0: OUTBUF[bi].sort(key=lambda x: x[ 1])
  BUF_MAX[bi] = OUTBUF[bi][ lbuf] #un-updated length

  while len( OUTBUF[bi]) >= BUFSIZE: yield OUTBUF[bi].pop( 0) # FIFO
  yield None
  
def buffered_flush( bi ):
  global OUTBUF, BUFSIZE, BUF_MAX
  if len( OUTBUF[bi]) > 1: 
    OUTBUF[bi].sort(key=lambda x: x[ 1])
    #info("buffered_flush: \n%s" % ("  \n".join([str(o[0:5]) for o in OUTBUF[bi]])))
  while len( OUTBUF[bi]) > 0: 
    #info("flush_yield: \n%s" % (str(OUTBUF[ 0][0:5])))
    yield OUTBUF[bi].pop( 0) #FIFO
  BUF_MAX[bi]=999999999999
  yield None
    

#Generator 
def VcfIter( vcf_filename, aOnlyStandardContigs=False, aIgnoreIndels=False, aIgnoreSnvs=False, aFilter=None,
             aOutputHeader=False, aExtraHeaderCol=False, aCoordinatesOnly=False, aSkipRef=False, aTabixContig=None):

    global QUIET

    buf_index = buffered_init()

    try:

      if aTabixContig != None: 
        if not os.path.isfile( vcf_filename + ".tbi"): error("Index file '%s' not found." % (vcf_filename + ".tbi"))
        process = subprocess.Popen(['tabix', '-f', vcf_filename, aTabixContig], stdout=subprocess.PIPE)
        file = process.stdout
      elif vcf_filename.endswith(".gz"): file = gzip.open( vcf_filename,'r')
      else: file = open( vcf_filename, "r")
  
    except:

      if aTabixContig != None:          
        error( "Could not process file '%s' with tabix. (%s)" % (vcf_filename, str( ex)))
      else:
        error("Could not open file '%s'" % aFile)

    linenum = 0
    prev_cols = [""]*5
    dupwarned = 0
    duplines = []
    n_errors = 0

    for line in file:

        linenum += 1
        if not line or len( line) == 0: continue

        if aOutputHeader and line.startswith("#"): 
          if aExtraHeaderCol and line.startswith("#CHROM"): 
            sys.stdout.write( "Present in file\t")          
          sys.stdout.write( line)
          continue

        if line.startswith("#"): continue        
        elif aFilter != None and line.find( aFilter) < 0: continue                
        ##CHROM  POS     ID      REF     ALT  QUAL  PASS   INFO
        try:
          raw_cols = line.split("\t", 5)
          if linenum == 1 and raw_cols[ 1] == "Start": continue # Annovar header             
          raw_cols[ 1] = int( raw_cols[ 1])

          # Remove 'chr 'to be able to compare "chr1" to "1"
          raw_cols[ 0] = FormatContig( raw_cols[ 0])  

          if aCoordinatesOnly: 
            yield raw_cols
            continue

          # Possibly multiple calls in ALT col
          # e.g. CAAA -> C,CA,CAA,CAAAA,CAAAAA
          # Split as if separate rows
          # TODO: read in and sort all lines on the same coordinate
          #       because they can be in different order in different files
          # TODO: ability to compare
          #       chr1    1175123 T       TA
          #       to
          #       chr1    1175123 TA      T,TAA 
          multicalls = sorted( map( str.strip, raw_cols[ 4].split(",")))

          for mcall_index in range( 0, len( multicalls)):

            # Create copy of row with only one ALT
            multi_call = raw_cols[:]
            alt = multicalls[ mcall_index]
            if aSkipRef and (alt == "." or alt == "" or alt == " " or alt == "-"): continue # No ALT
            multi_call[ 4] = alt

            multicols = []        
            mlen = len( multi_call[ 3])

            if mlen > 1 and mlen < 6 and mlen == len( multi_call[ 4]):
              # MNVs  e.g. AC -> GG
              # Split to multiple SNVs (rows)
              # Do not split larger MNVs than 5 bases
              for m in range( 0, mlen):
                split_cols = multi_call[:] #copy
                split_cols[ 1] += m
                split_cols[ 3] = split_cols[ 3][ m]
                split_cols[ 4] = split_cols[ 4][ m]
                multicols.append( split_cols)
            else:
              # SNVs, indels and deletions
              multicols.append( multi_call)
  
            # sort in place according to coordinate
            #multicols.sort(key=lambda x: x[ 1])

            for cols in multicols:

              if cols[ 3] == cols[ 4]: continue # Skip matching REF and ALT
              if aOnlyStandardContigs and cols[ 0] not in STD_CONTIGS: continue
              #two_to_two = (len( cols[ 3]) == 2 and len( cols[ 4]) == 2) #indel GG TT changed into two G T substions later
              if aIgnoreIndels and (len( cols[ 3]) != 1 or len( cols[ 4]) != 1): continue #Skip indels
              if aIgnoreSnvs and len( cols[ 3]) == len( cols[ 4]): continue
              if not IsInBedRange( cols[ 0], cols[1], BED): continue
              if IsIgnUsed() and IsInBedRange( cols[ 0], cols[ 1], IGN): continue
  
              if cols[:5] == prev_cols:
                  # Note: 
                  # Duplicate lines can occur like this when unpacking MNVs
                  # chr1    24032048        .       GC      AT
                  # chr1    24032049        .       C       T
                  #info("dupline: %s" % cols[:5])
                  duplines.append( linenum)
                  dupwarned += 1
                  continue # Skip line     
              prev_cols = cols[:5]

              # Buffered yield
              buf = buffered_yield( buf_index, cols)
              #info( "BUF %s: %s" % (vcf_filename, str(cols[ 0:5])))                  
              while True:
                bufout = buf.next()                
                if bufout == None: break
                #info( "OUT %s: %s" % (vcf_filename, str(bufout[ 0:5])))                  
                yield bufout
              #yield cols


        except Exception as ex:
          if n_errors > 100:
            error("Excessive number of errors (>100), please check the format of file '%s'." % vcf_filename)
          warning( "Bad line format on line %i in file '%s'. (%s)" % (linenum, vcf_filename, str(ex)))
          n_errors += 1
          continue

    try: file.close()
    except: pass

    # Buffered yield
    buf = buffered_flush( buf_index) 
    while True:
      bufout = buf.next()
      if bufout == None: break
      #info( "END %s: %s" % (vcf_filename, str(bufout[ 0:5])))      
      yield bufout

    if dupwarned > 0 and QUIET == False:
      warning("Skipped %i duplicate lines in file '%s'." % (dupwarned, vcf_filename))      
      warning("Duplicated lines in '%s': %s%s" % (os.path.basename(vcf_filename),(",".join( map( str, duplines)[:10])),("" if duplines < 10 else ",...")))
   
    yield None



NON_STD_CONTIGS = []
MAX_INDEX = 999999

def FormatContig( aContig):
  if aContig.startswith("chr") and aContig[3:] in STD_CONTIGS_DICT: aContig = aContig[3:]
  return aContig

def ToggleContigFormat( aContig):
  if aContig.startswith("chr") and aContig[3:] in STD_CONTIGS_DICT: aContig = aContig[3:]
  elif not aContig.startswith("chr"): aContig = ("chr" + aContig.upper())
  return aContig

def GetContigIndices( aContigs):

  n_contigs = len( aContigs)
  indices = [MAX_INDEX]*n_contigs

  for c in range( n_contigs):
    contig = aContigs[ c]
    if contig == None: 
      indices[ c] = MAX_INDEX
      continue
    
    contig = FormatContig( contig)

    try:      
      indices[ c] = STD_CONTIGS.index( contig)
    except ValueError: #non-std contigs
      # NOTE: Indices for non-std contigs are added in the order they are encountered
      # if they are not in the same order in every comapred file, matching lines will
      # be missed.
      if contig not in NON_STD_CONTIGS: NON_STD_CONTIGS.append( contig)
      indices[ c] = NON_STD_CONTIGS.index( contig) + len( STD_CONTIGS)

  return indices


def MinIndex( aDataCols ):

  n_datacols = len( aDataCols)
  chr_order = GetContigIndices( [(None if x == None else x[ 0]) for x in aDataCols])
  min_chr = min( chr_order)
  if min_chr == MAX_INDEX: return None

  #[ expression for item in list if conditional ]
  min_indices = [x for x in range( n_datacols) if chr_order[ x] == min_chr]
  coordinates = [int( aDataCols[ x][ 1]) for x in min_indices]
  min_coord = min( coordinates)
  
  # Check and sort by ALT col value  
  if coordinates.count( min_coord) > 1:
    min_candidate_indices = [x for x in range( n_datacols) if (int( aDataCols[ x][ 1]) == min_coord and chr_order[ x] == min_chr)]
    alts = [aDataCols[ x][ 4] for x in min_candidate_indices]
    first_alt = sorted( alts)[ 0]
    alt_index = alts.index( first_alt)
    dat_col_index = min_candidate_indices[ alt_index]
    return dat_col_index


  #info("coords: %s" % str(map( str, coordinates)))
  #info("min: %s" % str(min_coord))

  for i in min_indices:
    if aDataCols[ i][ 1] == min_coord:
      #print i, "MIN:", aDataCols[ i][0:2], "vs.", aDataCols[ 1 if i == 0 else 0][0:2]
      return i 

  info( "Datacols:\n" + aDataCols)
  error( "Minimum datacol not found.")
  return None


def Equal( aRow1, aRow2, aCoordinatesOnly=False):
  ##CHROM  POS     ID      REF     ALT  
  if aRow1 == None or aRow2 == None: return False
  #return aRow1[ 0] == aRow2[ 0] and aRow1[ 1] == aRow2[ 1] and aRow1[ 4] == aRow2[ 4]
  if aRow1[ 0] == aRow2[ 0] and aRow1[ 1] == aRow2[ 1]:
    if coordinates_only: return True
    return aRow1[ 4] == aRow2[ 4]    
  return False


def OutputRow( aData, aCurbLen=0):

  global USE_CHR_PREFIX_IN_OUTPUT

  chrom = str(aData[ 0])
  if chrom in STD_CONTIGS_DICT and USE_CHR_PREFIX_IN_OUTPUT: aData[ 0] = "chr" + chrom

  out = "\t".join( map( str, aData))
  if aCurbLen != None and aCurbLen > 0 and len( out) > aCurbLen: 
    out = out[:aCurbLen] + "\n"

  #sys.stdout.write( out) 
  try: sys.stdout.write( out) 
  except IOError: sys.exit( 0) 


def CompareFiles( aFiles, aNames=[], aOnlyStandardContigs=False, aFilter=None, aIgnoreIndels=False, aIgnoreSnvs=False,
                  aCoordinatesOnly=False, aPrintMatching=0, aPrintDifferent=0, aIdAddition=None, aCurbLinesTo=None,
                  aOutputHeader=False, aExtraHeaderCol=False, aUsePrefix=True, aOutForm="TMCP", aTabixContig=None,
                  aSkipRef=False):

  global QUIET

  n_files = len( aFiles)

  if len( aNames) == 0: aNames = aFiles[:]
  elif len( aNames) > n_files: aNames = aNames[ 0:n_files]
  elif len( aNames) < n_files: aNames += aFiles[ len( aNames):]  

  # Print header
  #sys.stdout.write( "COMMON\t")
  #for name in aNames:
  #  sys.stdout.write("%s\t" % name)         

  #for i in range( 2, n_files):
  #  sys.stdout.write( "IN %i\t" % i)
  #sys.stdout.write( "\n")

  # NOTE: Do not init matrices this way
  #       rows will be references to a single list
  # comp_matrix = [[0]*n_files]*n_files
  # count_matrix = [[0]*n_files]*n_files
  #
  # Init 2D array (matrix)
  comp_matrix = [[0]*n_files for i in range(n_files)]
  count_matrix = [[0]*n_files for i in range(n_files)]

  file_iters = []
  file_data = []
  finished = [False]*n_files

  if aOutputHeader and aExtraHeaderCol:    
    #sys.stdout.write("# Input files:\t")
    for name in aNames: sys.stdout.write("%s\t" % name)
    #sys.stdout.write("\n")

  # Create line generators
  for filename in aFiles:

    # Data is returned as a list of column values 
    file_iters.append( VcfIter( filename, aOnlyStandardContigs, aIgnoreIndels, aIgnoreSnvs, aFilter, 
                                aOutputHeader=aOutputHeader, aExtraHeaderCol=False, aCoordinatesOnly=aCoordinatesOnly,
                                aSkipRef=aSkipRef, aTabixContig=aTabixContig))    
    data = file_iters[ -1].next()

    if aTabixContig != None:
      if data == None:
        info("Retrying with contig '%s'." % ToggleContigFormat( aTabixContig))
        file_iters[ -1] = TabixIter(filename, ToggleContigFormat( aTabixContig), aIgnoreIndels, aIgnoreSnvs, aFilter, aOutputHeader=aOutputHeader)
        data = file_iters[ -1].next()
      if data == None:
        error( "File '%s' has no data rows for contigs '%s' or '%s'." % (filename, aTabixContig, ToggleContigFormat( aTabixContig))) 

    elif data == None:
      error( "File '%s' has no data rows." % filename)      




    aOutputHeader = False # only first file prints header
    file_data.append( data)

  total = [1]*n_files
  n_outputted_rows = 0
  
  # Manipulate ID column?
  set_addition = (aIdAddition != None and aIdAddition == "SET")  
  id_addition = None
  if aIdAddition != None: id_addition = str( aIdAddition)

  progress_counter = 0
  start_time = datetime.datetime.now()

  while True:

    # Find one of the iterators with the smallest contig & coordinate
    # (multiple iterators can be on the same contig & coordinate)
    # and compare it to all other iterators.
    mi = MinIndex( file_data)
    if mi == None: break # All Done?

    # Compare with self also
    rr = [Equal( file_data[ mi], file_data[ x], aCoordinatesOnly) for x in range( n_files)]
    #print "RR:", rr # DEBUG

    # Save results
    trues = [i for i, x in enumerate(rr) if x == True]
    n_trues = len( trues)    

    if set_addition and any( rr): 
      id_addition = "\t".join( ["TRUE" if x == True else "FALSE" for x in rr])      
      # Find first True index in rr
      first_true = next(i for i,v in enumerate(rr) if v)
      file_data[ first_true].insert( 0, id_addition)

      OutputRow( file_data[ first_true], aCurbLinesTo)
      n_outputted_rows += 1

    else: # set_addition presence indicator should not be used with -m or -d

      if aPrintMatching > 0 and (all( rr) == True):
        if id_addition != None: file_data[ 0][ 2] = id_addition
        OutputRow( file_data[ 0], aCurbLinesTo)
        n_outputted_rows += 1

        if aPrintMatching > 1:
          for i in range( 1, n_files):
            if id_addition != None: file_data[ i][ 2] = id_addition
            OutputRow( file_data[ i], aCurbLinesTo)
          n_outputted_rows += 1
      
      # Only two files in comparison  
      #['--different']): different = 1
      #['--all-different']): different = 2      
      if aPrintDifferent > 0:
        if rr[ 0] == True and rr[ 1] == False:
          if id_addition != None: file_data[ 0][ 2] = id_addition
          if aPrintDifferent == 2: sys.stdout.write("1:")
          OutputRow( file_data[ 0], aCurbLinesTo)
          n_outputted_rows += 1
        if aPrintDifferent == 2 and rr[ 1] == True and rr[ 0] == False: 
          if id_addition != None: file_data[ 1][ 2] = id_addition
          sys.stdout.write("2:")
          OutputRow( file_data[ 1], aCurbLinesTo)
          n_outputted_rows += 1

    # report progress
    progress_counter += 1
    if progress_counter > 100000:
      progress_counter = 0
      cur_time = datetime.datetime.now()
      elapsed_time = cur_time - start_time
      minutes_diff = elapsed_time.total_seconds() / 60.0
      if minutes_diff > 0.5: # Report ~twice a minute
        start_time = cur_time
        if len( file_data) > 0 and len( file_data[ 0]) >= 2:
          contig = file_data[ 0][ 0]
          if len(contig) <= 2: contig = "chr" + contig
          progress( "%s %s" % (contig, file_data[ 0][ 1]))

    # Self comparison (diagonal)
    # set diagonal true
    for t in trues:
    #  print "self:", t
      comp_matrix[ t][ t] += 1

    # +1 to all relevant matrix cells
    combinations = itertools.combinations( trues, 2)
    for row, col in combinations:
      comp_matrix[ row][ col] += 1
      # Mirrored element
      comp_matrix[ col][ row] += 1 
  
    for a in range( n_files):

      if finished[ a] == True: continue

      # Save result
      if rr[ a] == True:

        count_matrix[ a][ n_trues-1] += 1

        # Advance
        file_data[ a] = file_iters[ a].next()
        # If file end reached
        if file_data[ a] == None:  finished[ a] = True        
        else: total[ a] += 1

    if all( finished) == True: break
  
  if n_outputted_rows == 0: pass #info( "Nothing to output.")
  elif (aPrintMatching > 0 or aPrintDifferent > 0) and not QUIET: info( "Outputted %i rows." % n_outputted_rows)

  # With flags --different or --matching
  # Do not print statistics
  if aPrintMatching or aPrintDifferent: return

  if "T" in aOutForm:
    print ""
    if aCurbLinesTo == None:
      print "##This list shows the n of total SNV calls in each file (T)"
    print "Total number of calls in each file:"
    for i in range( n_files):
      print aNames[ i] + ": " + str( total[ i])

  if "M" in aOutForm:
    print ""
    if aCurbLinesTo == None:
      print "##This matrix shows the n of equal calls between two files (M)"
    print "Comparison matrix:"
    print "\t","\t".join( aNames)
    for i in range( n_files):
      print aNames[ i] + ":", "\t".join( map( str, comp_matrix[ i]))

  if "C" in aOutForm:
    print "" 
    if aCurbLinesTo == None:
      print "##This matrix shows in how many files each call made was found (C)"
    print "Count matrix:"
    #print count_matrix 
    print "\t","\t".join( aNames)
    for i in range( n_files):
      print ("IN %i:" % (i+1)), "\t".join( map( str, [count_matrix[ x][ i] for x in range( n_files)]))

  if "P" in aOutForm and n_files == 2:
    print ""
    if aCurbLinesTo == None:
      print "##Precision and Recall are calculated with 2nd file as ground truth (P)"
    FP = count_matrix[ 0][ 0]
    FN = count_matrix[ 1][ 0]
    TP = count_matrix[ 1][ 1]
    if TP == 0: 
      precision = 0.0
      recall = 0.0
    else:
      precision = TP/float( (TP+FP))
      recall = TP/float( (TP+FN))
    f_score = 2*((precision*recall)/(precision+recall))    
    print "[TP: %i, FP: %i, FN: %i]" % (TP,FP,FN)
    print "Precision: %.3f" % precision
    print "Recall: %.3f" % recall
    print "F-score: %.3f" % f_score
  elif "P" in aOutForm and n_files != 2:
    warning("Precision and recall only available when comparing exactly two files.")


  if aCurbLinesTo == None and not QUIET:
    print "All Done."
  else:
    print ""


def ReadBed( aFilename, aDict=BED, aQuiet=False):

  n_contigs = 0
  n_ranges = 0
  warned_end_coordinate = False
  is_vcf = aFilename.lower().endswith(".vcf") or aFilename.lower().endswith(".vcf.gz")
  line_num = 0

  try:

    if aFilename.endswith(".gz"): f = gzip.open( aFilename,'r')
    else: f = open( aFilename, "r")

    for line in f:
      line_num += 1
      if line.startswith("#"): continue #Comments
      if line.startswith("CHROM"): continue #Header
      cols = line.strip().split("\t")
  
      chromosome = cols[ 0]
      start = int(cols[ 1])
      try:
        if is_vcf: end = start
        elif cols[ 2] == ".": end = start
        else: end = int(cols[ 2])
      except:
        end = start
        if not warned_end_coordinate:
          warning("End coordinate column not found in file '%s' (Using start == end)." % aFilename)
          warned_end_coordinate = True
            
      chromosome = FormatContig( chromosome)
  
      if chromosome not in aDict: 
        aDict[ chromosome] = []
        n_contigs += 1
      aDict[ chromosome].append([start,end])
      n_ranges += 1

  except Exception as ex:
    error( "Failed to read bed file '%s'. Line %i: %s\n" % (aFilename, line_num, str(ex)))

  try: f.close()
  except: pass

  if len(aDict) < 1: error( "Bedfile '%s' contained %i contigs and %i ranges." % (aFilename, n_contigs, n_ranges))
  if not aQuiet and not QUIET: info( "Bedfile '%s' contained %i contigs and %i ranges." % (aFilename, n_contigs, n_ranges)) 

  # Sort ranges for each chromosome
  for k in aDict.keys():
    aDict[ k].sort(key=lambda x: x[ 0])





def IsBedUsed():
  return len( BED) > 0

def IsIgnUsed():
  return len( IGN) > 0  

BED_WARNED = {}

def IsInBedRange( aChromosome, aCoordinate, aDict=BED):  

  global QUIET, BED_WARNED

  if len( aDict) == 0: return True

  try:
    aCoordinate = int( aCoordinate)
  except:
    raise Exception("Could not convert '%s' to integer." % aCoordinate)
  
  aChromosome = FormatContig( aChromosome) # e.g. chr19 -> 19

  if aChromosome not in aDict:
    if aChromosome not in BED_WARNED:
      if not QUIET: warning("Contig '%s' not in bed file." % aChromosome)
      BED_WARNED[ aChromosome] = True
    return False

  for r in aDict[ aChromosome]:

    #DEBUG
    #if aChromosome == "19":
    #  between = (r[ 0] < aCoordinate) and (r[ 1] >= aCoordinate)
    #  
    #  if between: 
    #    sys.stderr.write("BED_chr19: ")
    #    sys.stderr.write("%i" % aCoordinate)        
    #    sys.stderr.write(" between: %i-%i\n" % (r[ 0], r[ 1]))
    #  #sys.stderr.write(" start: %i" % (r[ 0] < aCoordinate))
    #  #sys.stderr.write(" end: %i\n" % (r[ 1] >= aCoordinate))

    if r[ 0] > aCoordinate: break #start
    if r[ 1] >= aCoordinate: return True #end



  return False



if __name__ == '__main__':

    args = docopt.docopt(__doc__)

    report = False    
    if bool(args['--report']): report = True

    if bool(args['--no-prefix']): USE_CHR_PREFIX_IN_OUTPUT = False

    only_standard = False
    if bool(args['--standard']): only_standard = True

    coordinates_only = False
    if bool(args['--coordinates']): coordinates_only = True

    ignore_indels = False
    if bool(args['--ignore-indels']): ignore_indels = True

    ignore_snvs = False
    if bool(args['--ignore-snvs']): ignore_snvs = True

    # Sanity check
    if ignore_indels and ignore_snvs:
      error( "Cannot ignore all entries in comparison.")      

    baseratios = False
    if bool(args['--baseratios']): baseratios = True

    different = 0
    if bool(args['--different']): different = 1
    if bool(args['--all-different']): different = 2

    pretty = None
    if bool(args['--pretty']): pretty = 80   

    add_id = args['--add-id']
    extra_header_col = False
    if bool(args['--add-info']): 
      add_id = "SET"
      extra_header_col = True
      if pretty:
        warning("--pretty cannot be used with --add-info.")
        pretty = None

    matching = 0
    if bool(args['--matching']): matching = 1
    if bool(args['--all-matching']): matching = 2

    filenames = args['<file>']    
    n_files = len( filenames)
    
    if bool(args['--swap']): filenames = filenames[::-1]

    if bool(args['--quiet']): QUIET = True

    skip_ref = False
    if bool(args['--skip-ref']): skip_ref = True

    header = False
    if bool(args['--header']): header = True
    

    col_names = args['--names']
    if col_names == None or len( col_names) == 0: col_names = []
    else: 
      col_names = col_names.split(",")
      if n_files != len( col_names): warning( "Number of names should equal number of input files (%ivs.%i)" % (len(col_names),n_files))

    filterstr = args['--filter']

    tabix_str = args['--tabix']

    outform = args['--outform']
    if outform == None or len( outform) == 0: outform = "TMC"
    else: 
      outform = outform.upper()
      if "A" in outform: outform = "TMCP"

    bedfile = args['--bed']
    if bedfile != None and len( bedfile) > 0: ReadBed( bedfile, BED, aQuiet=(pretty!=None))

    ignfile = args['--ignore']
    if ignfile != None and len( ignfile) > 0: ReadBed( ignfile, IGN, aQuiet=(pretty!=None))

    # Check that all files exist
    non_existent_files = []

    if different > 0 and n_files != 2:
      error("Option '--different' requires exactly two input files.")
    if matching > 0 and n_files > 2:      
      info("Outputting SNVs found in all %i input files." % n_files)
      #error("Option '--matching' requires exactly two input files.")

    for f in filenames:      
      if not os.path.isfile( f): 
        non_existent_files.append( f)

    if len( non_existent_files) > 0:
      if len( non_existent_files) > 1: error( "Files '%s' do not exists." % ", ".join( non_existent_files))
      else: error( "File '%s' does not exist." % non_existent_files[ 0])

    # Single input file, report stats only
    if n_files == 1: 
      report = True      
      PrintReport( filenames[ 0], baseratios)
      sys.exit( 0)


    # Multiple input files
    #for f1 in range( n_files):
    #  for f2 in range( f1+1, n_files):
    if report:
      for f in range( n_files):      
        PrintReport( filenames[ f], None if f >= len( col_names) else col_names[ f], baseratios)
    
    try:

      CompareFiles( filenames, col_names, only_standard, filterstr, ignore_indels, ignore_snvs,
                    aCoordinatesOnly=coordinates_only, aPrintMatching=matching, aPrintDifferent=different,
                    aIdAddition=add_id, aCurbLinesTo=pretty, aOutputHeader=header, aExtraHeaderCol=extra_header_col, 
                    aOutForm=outform, aTabixContig=tabix_str, aSkipRef=skip_ref)
    except Exception as ex:
      pass
    
    #print "Processed filenames:", filenames
    sys.exit(0)
   


