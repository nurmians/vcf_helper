#!/usr/bin/python2.7

"""


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
  -i --indels          Compare indels only (not implemented yet)

"""

import docopt
import subprocess, sys, os, datetime, time, re
import gzip
import itertools

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
STD_BASES = ["A","C","G","T"]

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

  if aName != None:    
    print "\nStatistics for '%s':" % aName
  else:
    print "\nStatistics for file '%s':" % aFile

  print "n_Contigs: ", len( contig_arr)
  #print "Contigs: ", contig_arr
  print "Calls per contig: ", ["%s:%i" % (k,contigs[ k]) for k in contig_arr]

  for k in stats.keys():
    print "%s: %s" % (k, stats[ k])

  if aPrintBaseRatios:
    n_subs = int( stats["sub"])
    for k in sorted( subs.keys()):
      print "%s: %3.1f%%" % (k, int(subs[ k]) / float(n_subs)*100.0 )



#Generator
def VcfIter( vcf_filename, aOnlyStandardContigs=False):

    try:
      if vcf_filename.endswith(".gz"): file = gzip.open( vcf_filename,'r')
      else: file = open( vcf_filename, "r")
  
    except:
      error("Could not open file '%s'" % aFile)

    linenum = 0

    for line in file:

        linenum += 1
        if not line or len( line) == 0: continue
        elif line.startswith("#"): continue
        ##CHROM  POS     ID      REF     ALT  
        try:
          cols = line.split("\t", 5)
          cols[ 1] = int( cols[ 1])
          if cols[ 0].startswith("chr"): 
            cols[ 0] = cols[ 0][3:]      

          if aOnlyStandardContigs and cols[ 0] not in STD_CONTIGS: continue

        except:
          warning( "Bad line format on line %i in file '%s'" % (linenum, vcf_filename))
          continue
    
        yield cols


    try: file.close()
    except: pass
    yield None


NON_STD_CONTIGS = []
MAX_INDEX = 999999

def GetContigIndices( aContigs):

  n_contigs = len( aContigs)
  indices = [MAX_INDEX]*n_contigs

  for c in range( n_contigs):
    contig = aContigs[ c]
    if contig == None: 
      indices[ c] = MAX_INDEX
      continue

    if contig.startswith("chr"): contig = contig[3:]
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
  coordinates = [aDataCols[ x][ 1] for x in min_indices]
  min_coord = min( coordinates)

  for i in min_indices:
    if aDataCols[ i][ 1] == min_coord:
      #print i, "MIN:", aDataCols[ i][0:2], "vs.", aDataCols[ 1 if i == 0 else 0][0:2]
      return i 

  info( "Datacols:\n" + aDataCols)
  error( "Minimum datacol not found.")
  return None


def Equal( aRow1, aRow2):
  ##CHROM  POS     ID      REF     ALT  
  if aRow1 == None or aRow2 == None: return False
  return aRow1[ 0] == aRow2[ 0] and aRow1[ 1] == aRow2[ 1] and aRow1[ 4] == aRow2[ 4]


def CompareFiles( aFiles, aNames=[], aOnlyStandardContigs=False):

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
  # rows will be references to a single array
  # comp_matrix = [[0]*n_files]*n_files
  # count_matrix = [[0]*n_files]*n_files
  #
  # Init 2D array (matrix)
  comp_matrix = [[0]*n_files for i in range(n_files)]
  count_matrix = [[0]*n_files for i in range(n_files)]

  file_iters = []
  file_data = []
  finished = [False]*n_files  

  for filename in aFiles:

    file_iters.append( VcfIter( filename, aOnlyStandardContigs))
    data = file_iters[ -1].next()
    if data == None:
      error( "File '%s' has no data rows.")
    file_data.append( data)


  total = [1]*n_files

  while True:

    mi = MinIndex( file_data)
    if mi == None: break # All Done?

    # Compare with self also
    rr = [Equal( file_data[ mi], file_data[ x]) for x in range( n_files)]
    #print "RR:", rr

    #Save results
    trues = [i for i, x in enumerate(rr) if x == True]
    n_trues = len( trues)
    #print "TRues", trues 

    # Self comparison (diagonal)
    for t in trues: 
    #  print "self:", t
      comp_matrix[ t][ t] += 1

    combinations = itertools.combinations( trues, 2)
    for row, col in combinations:
      comp_matrix[ row][ col] += 1
      # Mirrored element
      comp_matrix[ col][ row] += 1 
      pass     

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

  print ""
  print "Total number of calls in each file:"
  for i in range( n_files):
    print aNames[ i] + ": " + str( total[ i])

  print ""
  print "##This matrix shows the number of equal calls between two files"
  print "Comparison matrix:"
  print "\t","\t".join( aNames)
  for i in range( n_files):
    print aNames[ i] + ":", "\t".join( map( str, comp_matrix[ i]))

  print ""
  print "##This matrix shows in how many files each call made was found"
  print "Count matrix:"
  #print count_matrix
  print "\t","\t".join( aNames)
  for i in range( n_files):
    print ("IN %i:" % (i+1)), "\t".join( map( str, [count_matrix[ x][ i] for x in range( n_files)]))

  print "All Done."














if __name__ == '__main__':

    args = docopt.docopt(__doc__)

    report = False    
    if bool(args['--report']): report = True

    only_standard = False
    if bool(args['--standard']): only_standard = True

    #TODO
    only_indels = False
    if bool(args['--indels']): only_indels = True

    baseratios = False
    if bool(args['--baseratios']): baseratios = True

    filenames = args['<file>']    
    n_files = len( filenames)

    col_names = args['--names']
    if col_names == None or len( col_names) == 0: col_names = []
    else: 
      col_names = col_names.split(",")
      if n_files != len( col_names): warning( "Number of names should equal number of input files (%ivs.%i)" % (len(names),n_files))

    # Check that all files exists
    non_existent_files = []

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

    CompareFiles( filenames, col_names, only_standard)







    #print "Filenames:", filenames
    sys.exit(0)
   

