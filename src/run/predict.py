from zbio import gtf, bam, ribo, stat, exp, tools, orf, fa, interval, io
import math, time, itertools
from os.path import isfile

def help():
  return "Predict ORFs from riboseq data"
def set_parser(parser):
  #### basic input options ####
  parser.add_argument("-t", type=strlist, dest="tisbampaths", default=[], help="TIS enriched riboseq bam files, comma seperated")
  parser.add_argument("-b", type=strlist, dest="ribobampaths", default=[], help="Ordinary riboseq bam files, comma seperated")
  parser.add_argument("-g", type=str, dest="genepath", required=True, help='Gene annotation file for TIS background estimation and ORF prediction')
  parser.add_argument("-f", type=str, dest="genomefapath", required=True, help="Genome fasta file")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output result file")
  #### alt input options ####
  parser.add_argument("-i", type=str, dest="input", help="Only test input candidate ORFs, format: transID start stop (0 based, half open)")
  parser.add_argument("--geneformat", type=str, default='auto', help="Gene annotation file format (gtf, bed, gpd, gff, default: auto)")
  parser.add_argument("--tispara", type=strlist, help="Input offset parameter files for -t bam files")
  parser.add_argument("--ribopara", type=strlist, help="Input offset parameter files for -b bam files")
  parser.add_argument("--nparts", type=int, default=10, help="Group transcript according to TIS reads density quantile (default: 10)")
  parser.add_argument("-a", type=str, dest="agenepath", help="Gene file for ORF prediction instead of -g gene file")
  parser.add_argument("-e", type=str, dest="estpath", default='tisBackground.txt', help="Output TIS background estimation result (default: tisBackground.txt)")
  parser.add_argument("-s", type=str, dest="inestpath", help="Input background estimation result file instead of instant estimation")
  parser.add_argument("--alt", action="store_true", help="Use alternative start codons (all codons with 1 base different from ATG)")
  parser.add_argument("--altcodons", type=strlist, help="Use provided alternative start codons, comma seperated, eg. CTG,GTG,ACG")
  #parser.add_argument("--addchr", action="store_true", help="Auto add 'chr' for input gene annotation file")
  parser.add_argument("--tis2ribo", action="store_true", help="Add TIS bam counts to ribo, if specified or -b not provided")
  parser.add_argument("--harr", action="store_true", help="The data is treated with harringtonine (instead of LTM)")
  parser.add_argument("--harrwidth", type=int, default=15, help="Flanking region for harr data, in codons (default: 15)")
  parser.add_argument("--enrichtest", action="store_true", help="Use enrich test instead of frame test")
  parser.add_argument("--nocompatible", action="store_true", help="Do not require reads compatible with transcript splice junctions")
  #parser.add_argument("--maxtisnum", type=int, default=10, help="max number of TIS sites for one transcript")
  parser.add_argument("--minaalen", type=int, default=6, help="Min amino acid length of candidate ORF (default: 6)")
  parser.add_argument("--genefilter", type=strlist, help="Only process given genes")
  parser.add_argument("--tpth", type=float, default=1, help="TIS p value threshold (default: 1)")
  parser.add_argument("--fpth", type=float, default=1, help="Frame p value threshold (default: 1)")
  parser.add_argument("--minpth", type=float, default=0.05, help="At least one of TIS or frame p value should be lower than this threshold (default: 0.05)")
  #parser.add_argument("--epth", type=float, default=1, help="Enrichment p value threshold")
  #parser.add_argument("--fspth", type=float, default=0.05, help="Fisher's p value threshold")
  #parser.add_argument("--qth", type=float, default=1, help="FDR q value threshold")
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")
  parser.add_argument("-v", "--verbose", action="count", help="Increase output verbosity")
  #parser.add_argument("--showtime", action="store_true", help="showtime")

  
def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')

use_message = '''
python tisRibo.py [options] -g gtfpath -t tisbampath -b ribobampath -f genomefapath -o output 
options:
-p : number of processes (1)
--parts parts : group transcript according to TIS reads density quantile (10) 
-a gtf_for_pred : GTF file for ORF prediction. Default: gtfpath
-e estpath : TIS background estimation result. (tisbampath+'.bgest.txt' if not exist)
-s inestpath : Input background estimation result file instead of instant estimation
--addchr : Auto add 'chr' for input GTF file
--tispara parapath : input parameter file for tisbam
--ribopara parapath : input parameter file for ribobam
--tis2ribo : Add TIS bam counts to ribo, if specified or -b not provided
--harr : The RPF data is treated with harringtonine (instead of LTM)
--harrwidth <int> : Flanking region for harr data TIS search, in codons (15) aa
--genefilter geneid : only process given genes
-v --verbose : show more details to STDOUT
--veryverbose : show extreme details to STDOUT
--tipth <float> : TIS p value threshold. Default: 0.05
--epth <float> : Enrichment p value threshold. Default: 1 #show all results and summarize in the next step
--fpth <float> : Frame p value threshold. Default: 1
--qth <float> : FDR q value threshold. Default: 1
'''

def run(args):
  '''Main function for ORF finding
  '''
  # prepare
  global tisbampaths, tisoffdict, ribobampaths, riboffdict, genomefapath, compatible
  global minaalen, enrichtest, slp, paras, verbose, alt, title, tis2ribo, gfilter
  global tpth, fpth, minpth #fspth
  #global showtime
  #showtime = args.showtime
  tisbampaths = args.tisbampaths
  ribobampaths = args.ribobampaths
  if len(tisbampaths) == 0 and len(ribobampaths) == 0 :
    print('No bam file input!')
    exit(1)
  genomefapath = args.genomefapath
  compatible = not args.nocompatible
  minaalen = args.minaalen
  enrichtest = args.enrichtest
  harrwidth = None
  if args.harrwidth is not None : harrwidth = args.harrwidth
  elif args.harr : harrwidth = 15
  verbose = args.verbose
  alt = args.alt
  if args.altcodons is not None : 
    alt = True
    orf.cstartlike = [c.upper() for c in args.altcodons]
  tpth, fpth, minpth = args.tpth, args.fpth, args.minpth # fspth
  tis2ribo = args.tis2ribo
  parts = [0.1 * (i+1) for i in range(args.nparts)]
  gfilter = None
  if args.genefilter is not None :
    gfilter = {}
    for gid in args.genefilter : gfilter[gid] = 1
  flank = 3 ##
  tisoffdict = find_offset(args.tisbampaths, args.tispara)
  riboffdict = find_offset(args.ribobampaths, args.ribopara)
  if len(args.ribobampaths) == 0 : 
    print('No ordinary riboseq data input. TIS data will also be used as ordinary riboseq data.')
    tis2ribo = True
  if len(args.tisbampaths) == 1 : 
    if args.inestpath is None : 
      path = args.tisbampaths[0] + '.bgest.txt'
      if isfile(path) : args.inestpath = path
      else : args.estpath = path
  if args.agenepath is None : args.agenepath = args.genepath
  #start = time.time()

  # load genome, fasta file indexing
  if args.verbose : print("Loading genome...")
  genome = fa.Fa(args.genomefapath)

  # TIS background estimation
  if len(args.tisbampaths) == 0 : 
    print('No input TIS data!')
    paras, slp = [(1,0.5)], [1] # No TIS input
  elif args.inestpath is None : #== '' :
    print ("Estimate TIS parameters...")
    if args.verbose : print("TIS estimation result will be saved to {}".format(args.estpath))
    from multiprocessing import Process
    import multiprocessing.pool
    class NoDaemonProcess(Process):
    # make 'daemon' attribute always return False
      def _get_daemon(self):
        return False
      def _set_daemon(self, value):
        pass
      daemon = property(_get_daemon, _set_daemon)
    class MyPool(multiprocessing.pool.Pool):
      Process = NoDaemonProcess

    pool = MyPool(1) # This is for memory efficiency
    paras, slp, data = pool.apply(ribo.estimateTISbg, args=(args.genepath, args.tisbampaths, args.genomefapath), kwds={'parts': parts, 'offdict': tisoffdict, 'numProc': args.numProc, 'verbose': args.verbose, 'geneformat': args.geneformat, 'harrwidth': harrwidth})
  #paras, slp, data = ribo.estimate_tis_bg_inframe(gtfpath, tisbampath, genomefapath, parts = parts, offdict = tisoffdict, addchr = addchr, numProc = numProc, verbose = verbose)
    estfile = open(args.estpath, 'w')
    for i in range(len(parts)):
      estfile.write("{}\t{}\t{}\t{}\t{}\n".format(paras[i][0], paras[i][1], parts[i], slp[i], data[i]))
    estfile.close()
    pool.close()
  else : 
    inestfile = open(args.inestpath, 'r')
    paras, slp = [], []
    for l in inestfile:
      lst = l.strip().split('\t')
      paras.append((float(lst[0]), float(lst[1])))
      slp.append(eval(lst[3]))
  
  inorf = None
  if args.input is not None :
    inorf = {}
    infile = open(args.input, 'r')
    for l in infile :
      lst = l.strip().split()
      tid, tis, stop = lst[0], int(lst[1]), int(lst[2])
      if gfilter is not None and tid not in gfilter : continue
      if tid not in inorf : inorf[tid] = []
      inorf[tid].append([tis, stop])
  print("Predicting...")
  profile = exp.Profile()
  #if enrichtest : title = ['TISGroup', 'TISCounts', 'TISPvalue', 'EnrichPvalue', 'EnrichPStatus']
  title = ['TISGroup', 'TISCounts', 'TISPvalue', 'RiboPvalue', 'RiboPStatus']
  j = 0
  #agenefile = open(args.agenepath,'r')
  trans_iter = io.transIter(args.agenepath, fileType = args.geneformat, chrs = genome.idx, verbose = args.verbose)
  para_iter = transPara(trans_iter, inorf)
  #para_iter = itertools.izip(gene_iter, itertools.repeat(paras), itertools.repeat(slp))
  if args.numProc <= 1 : pred_iter = itertools.imap(_pred_trans, para_iter)
  else : 
    from multiprocessing import Pool
    pool = Pool(processes = args.numProc - 1)
    pred_iter = pool.imap_unordered(_pred_trans, para_iter, chunksize = 5)
  for result in pred_iter:
    es, ji = result
    j += ji
    for e in es : 
      profile.add_exp(e)
      if verbose >= 2 : print e
  print("BH correcting...")
  profile.BHcorrection(2, total = j, append = True) # Calculate BH FDR of TIS p value
  profile.BHcorrection(3, total = j, append = True) # Frame p value
  #profile.BHcorrection(5, total = j, append = True) # Fisher's p value
  outfile = open(args.output,'w')
  s = "Gid\tTid\tSymbol\tGeneType\tGenomePos\tStartCodon\tStart\tStop\tTisType\t"
  s += '\t'.join(title)
  s += '\tTISQvalue\tRiboQvalue\tAALen\n'
  outfile.write(s)

  for e in profile:
    #e.data.append(e.q)
    #if e.data[3] > args.epth or e.data[4] > args.fpth or e.q > args.qth : continue
    outfile.write("%s\t%d\n" % (e, e.length)) #, e.sq))
  
  #end = time.time()
  #print('Time used: %s' % str(end - start))

def transPara(trans_iter, inorf):
  '''Generate parameters (trans, candidates/None) for function _pred_trans()
  '''
  if inorf is not None :
    for t in trans_iter:
      if t.id in inorf : yield t, inorf[t.id]
  else :
    #i = 0
    for t in trans_iter: 
      if gfilter is not None :
        if t.id not in gfilter or t.gid not in gfilter : continue
      yield t, None
      #i += 1
      #if i >= 10 : break
#offdict = None
def find_offset(bampaths, para):
  '''Get offset data for given bam data
  '''
  offlist = [None] * len(bampaths)
  for i, path in enumerate(bampaths):
    if para is None or i > len(para)-1 or para[i] == '' : path = path + '.para.py'
    else : path = para[i]
    if isfile(path) : 
      vessel = {}#'offdict': None}
      exec(open(path).read(), vessel) # execfile(path, vessel)
      offlist[i] = vessel['offdict']
  for i, od in enumerate(offlist):
    if od is None :
      print('No offset parameter file found for %s. Using default offset (12). ' %  bampaths[i])
  return offlist
    


def _pred_trans(ps): ### trans
  '''Main function of ORF prediction in given transcript
  '''
  #if showtime : timestart = time.time()
  t, candidates = ps
  es, j = [], 0
  tl = t.cdna_length()
  if tl < ribo.minTransLen : return es, j ##
  ttis = ribo.multiRibo(t, tisbampaths, offdict = tisoffdict, compatible = compatible)
  tribo = ribo.multiRibo(t, ribobampaths, offdict = riboffdict, compatible = compatible)
  #if showtime : time1 = time.time()
  score = ttis.abdscore()
  ip = ribo.pidx(score, slp) 
  genome = fa.Fa(genomefapath)
  if len(tisbampaths) > 0 and tis2ribo : tribo.merge(ttis) ##
  if verbose >= 2 : print g.id, t.id, ttis.total, tribo.total
  cds1 = t.cds_start(cdna = True) 
  cds2 = t.cds_stop(cdna = True) 
  tsq = genome.transSeq(t)
  # user provided candidates
  if candidates is not None : 
    for tis, stop in candidates:
      j += 1
      if len(tisbampaths) > 0 : tp = ttis.tis_test(tis, paras[ip][0], paras[ip][1])
      else : tp = None
      if enrichtest : rp = tribo.enrich_test(tis, stop)
      else : rp = tribo.frame_test(tis, stop)
      #fisher = stat.fisher_method([tp, rp])[0]
      if tp is not None and tp > tpth : continue 
      if rp > fpth : continue # or fisher > fspth 
      minp = rp
      if tp is not None and tp < minp : minp = tp
      if minp > minpth : continue
      tistype = tisType(tis, stop, cds1, cds2)
      orfstr = '{}\t{}\t{}'.format(tsq[tis:tis+3],tis,stop)
      tid = "%s\t%s\t%s\t%s\t%s:%d:%s\t%s\t%s" % (t.gid, t.id, t.symbol, t.genetype, t.chr, t.genome_pos(tis), t.strand, orfstr, tistype)
      values = [ip, ttis.cnts[tis], tp, rp, 'N']
      e = exp.Exp(tid, values)
      e.length = (stop - tis) / 3 - 1
      #e.sq = tsq[tis:stop]
      es.append(e)
    #if showtime : 
      #end = time.time()
      #print('%s\t%s\tTime_used:\t%s\t%s' % (t.symbol, t.id, str(time1 - timestart), str(end - time1)))
    return es, j
  # else : all possible ORFs
  orfs = orf.orflist(tsq, minaalen = minaalen)
  for o in orfs :
    starts = o.starts
    if alt : starts += o.altstarts
    starts.sort()
    ol = len(starts)
    tps = [None] * ol
    rps = [None] * ol
    #fishers = [None] * ol
    for i, tis in enumerate(starts) : 
      if len(tisbampaths) > 0 : tps[i] = ttis.tis_test(tis, paras[ip][0], paras[ip][1])
      if enrichtest : rps[i] = tribo.enrich_test(tis, o.stop)
      else : rps[i] = tribo.frame_test(tis, o.stop)
      #fishers[i] = stat.fisher_method([tps[i], rps[i]])[0]
    rst = pvalStatus(rps)
    for i, tis in enumerate(starts) : 
      if tps[i] is not None and tps[i] > tpth : continue
      if rps[i] > fpth : continue # or fishers[i] > fspth
      minp = rps[i]
      if tps[i] is not None and tps[i] < minp : minp = tps[i]
      if minp > minpth : continue
      tistype = tisType(tis, o.stop, cds1, cds2)
      orfstr = '{}\t{}\t{}'.format(tsq[tis:tis+3],tis,o.stop)
      tid = "%s\t%s\t%s\t%s\t%s:%d:%s\t%s\t%s" % (t.gid, t.id, t.symbol, t.genetype, t.chr, t.genome_pos(tis), t.strand, orfstr, tistype)
      values = [ip, ttis.cnts[tis], tps[i], rps[i], rst[i]] # , fishers[i]]
      e = exp.Exp(tid, values)
      e.length = (o.stop - tis) / 3 - 1
      #e.sq = tsq[tis:o.stop]
      es.append(e)
    j += ol
  #if showtime : 
      #end = time.time()
      #print('%s\t%s\tTime_used:\t%s\t%s' % (t.symbol, t.id, str(time1 - timestart), str(end - time1)))
  return es, j

def pvalStatus(ps) : 
  '''Find global best (T) & local best (L) TISs that correspond to same stop codon, suggesting potential TISs by ordinary riboseq data
  '''
  l = len(ps)
  st = ['N'] * l
  if l <= 0 : return st
  cp = ps[:] + [1]
  m, c = 1, 1
  for i, p in enumerate(cp): 
    if p is None : cp[i] = 1
    elif p < m : m, c = p, 1
    elif p == m : c += 1
  if m >= 1 : return st
  top = 'T{}'.format(c)
  last, c = 1, 1
  better = True
  for i, p in enumerate(cp): 
    #print (p,m)
    if p == m : st[i] = top
    elif p < last : 
      c = 1
      better = True
    elif p == last : c += 1
    else : 
      if better and last > m : 
        local = 'L{}'.format(c)
        for j in range(i-c, i) : st[j] = local
      better = False
      c = 1
    last = p
  return st
def tisType(start, stop, cds1, cds2):
  '''Type of the TIS relatve to annotated ORF
  '''
  if cds1 is None : return 'Novel'
  else :
    if start == cds1 : return 'Annotated'
    elif cds1 < stop <= cds2 and (stop - cds2) % 3 == 0 : # in frame
      if start > cds1 : return 'Truncated'
      else : return 'Extended'
    else : 
      if start < cds1 : return "5'UTR"
      elif start > cds2 : return "3'UTR"
      else : return "Inside"

if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())


    