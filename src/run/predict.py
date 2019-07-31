from ribotish.zbio import gtf, bam, ribo, stat, exp, tools, orf, fa, interval, io
import math, time, itertools
from os.path import isfile

def help():
  return "Predict ORFs from riboseq data"
def set_parser(parser):
  #### basic options ####
  parser.add_argument("-t", type=strlist, dest="tisbampaths", default=[], help="TIS enriched riboseq bam files, comma seperated")
  parser.add_argument("-b", type=strlist, dest="ribobampaths", default=[], help="Ordinary riboseq bam files, comma seperated")
  parser.add_argument("-g", type=str, dest="genepath", required=True, help='Gene annotation file for ORF prediction')
  parser.add_argument("-f", type=str, dest="genomefapath", required=True, help="Genome fasta file")
  parser.add_argument("-o", type=str, dest="output", default='predict.txt', required=True, help="Output result file")
  #### alt input options ####
  parser.add_argument("-i", type=str, dest="input", help="Only test input candidate ORFs, format: transID start stop (0 based, half open)")
  parser.add_argument("--geneformat", type=str, default='auto', help="Gene annotation file format (gtf, bed, gpd, gff, default: auto)")
  parser.add_argument("--tispara", type=strlist, help="Input offset parameter files for -t bam files")
  parser.add_argument("--ribopara", type=strlist, help="Input offset parameter files for -b bam files")
  parser.add_argument("--nparts", type=int, default=10, help="Group transcript according to TIS reads density quantile (default: 10)")
  parser.add_argument("-a", type=str, dest="agenepath", help="Gene file for known protein coding gene annotation and TIS background estimation instead of -g gene file")
  parser.add_argument("-e", type=str, dest="estpath", default='tisBackground.txt', help="Output TIS background estimation result (default: tisBackground.txt)")
  parser.add_argument("-s", type=str, dest="inestpath", help="Input background estimation result file instead of instant estimation")
  parser.add_argument("--transprofile", type=str, help="Output RPF P-site profile for each transcript")
  parser.add_argument("--inprofile", type=str, help="Input RPF P-site profile for each transcript, instead of reading bam reads, save time for re-running")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if annotation chr ids are not same as chr ids in bam/fasta files")
  parser.add_argument("--allresult", type=str, help="All result output without FDR q-value threshold (default: output + '_all.txt', 'off' to turn off)")
  # predict options
  parser.add_argument("--alt", action="store_true", help="Use alternative start codons (all codons with 1 base different from ATG)")
  parser.add_argument("--altcodons", type=strlist, help="Use provided alternative start codons, comma seperated, eg. CTG,GTG,ACG")
  parser.add_argument("--tis2ribo", action="store_true", help="Add TIS bam counts to ribo, if specified or -b not provided")
  parser.add_argument("--harr", action="store_true", help="The data is treated with harringtonine (instead of LTM)")
  parser.add_argument("--harrwidth", type=int, default=15, help="Flanking region for harr data, in codons (default: 15)")
  parser.add_argument("--enrichtest", action="store_true", help="Use enrich test instead of frame test")
  # output options
  parser.add_argument("--minaalen", type=int, default=6, help="Min amino acid length of candidate ORF (default: 6)")
  parser.add_argument("--genefilter", type=strlist, help="Only process given genes")
  parser.add_argument("--tpth", type=float, default=0.05, help="TIS p value threshold (default: 0.05)")
  parser.add_argument("--fpth", type=float, default=0.05, help="Frame p value threshold (default: 0.05)")
  parser.add_argument("--minpth", type=float, default=1, help="At least one of TIS or frame p value should be lower than this threshold (default: 0.05)")
  parser.add_argument("--fspth", type=float, default=0.05, help="Fisher's p value threshold")
  parser.add_argument("--fsqth", type=float, default=0.05, help="Fisher's FDR q value threshold")
  parser.add_argument("--framelocalbest", action="store_true", help="Only report local best frame test results")
  parser.add_argument("--framebest", action="store_true", help="Only report best frame test results")
  parser.add_argument("--longest", action="store_true", help="Only report longest possible ORF results")
  parser.add_argument("--seq", action="store_true", help="Report ORF sequences")
  parser.add_argument("--aaseq", action="store_true", help="Report amino acid sequences")
  parser.add_argument("--blocks", action="store_true", help="Report all exon block positions for predicted ORFs")
  # Reads filters
  parser.add_argument("--maxNH", type=int, default=5, help="Max NH value allowed for bam alignments (default: 5)")
  parser.add_argument("--minMapQ", type=float, default=1, help="Min MapQ value required for bam alignments (default: 1)")
  parser.add_argument("--secondary", action="store_true", help="Use bam secondary alignments")
  parser.add_argument("--paired", action="store_true", help="Reads are paired end")
  parser.add_argument("--nocompatible", action="store_true", help="Do not require reads compatible with transcript splice junctions")
  parser.add_argument("--compatiblemis", type=int, default=2, help="Missed bases allowed in reads compatibility check (default: 2)")
  # running options
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")
  parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase output verbosity")

def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')

try: imap = itertools.imap
except: imap = map

def run(args):
  '''Main function for ORF finding
  '''
  # prepare
  global tisbampaths, tisoffdict, ribobampaths, riboffdict, genomefapath, compatible, compatiblemis
  global minaalen, enrichtest, slp, paras, verbose, alt, title, tis2ribo, gfilter
  global tpth, fpth, minpth, fspth, framebest, framelocalbest, longest, transprofile, TIS_types #fspth
  global paired, seq, aaseq, blocks # showtime
  paired, seq, aaseq, blocks = args.paired, args.seq, args.aaseq, args.blocks
  ribo.maxNH, ribo.minMapQ, ribo.secondary = args.maxNH, args.minMapQ, args.secondary
  tisbampaths = args.tisbampaths
  ribobampaths = args.ribobampaths
  if len(tisbampaths) == 0 and len(ribobampaths) == 0 :
    print('No bam file input!')
    exit(1)
  genomefapath = args.genomefapath
  compatible = not args.nocompatible
  compatiblemis = args.compatiblemis
  minaalen = args.minaalen
  enrichtest = args.enrichtest
  transprofile = args.transprofile
  harrwidth = None
  TIS_types = ['Annotated', 'Truncated', 'Extended', "5'UTR", "3'UTR", 'Internal', 'Novel']
  if args.chrmap is not None :
    chrmap = {}
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]
    bam.chrmap = chrmap
    fa.chrmap = chrmap
  if args.harrwidth is not None : harrwidth = args.harrwidth
  elif args.harr : harrwidth = 15
  verbose = args.verbose
  alt = args.alt
  if args.altcodons is not None : 
    alt = True
    if args.altcodons[0].upper() == 'ALL' : orf.cstartlike = orf.allcodons
    else : orf.cstartlike = [c.upper() for c in args.altcodons]
  tpth, fpth, minpth, framebest, framelocalbest = args.tpth, args.fpth, args.minpth, args.framebest, args.framelocalbest # fspth
  fspth = args.fspth
  longest = args.longest
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
    print('No regular RiboSeq data input. TIS data will also be used as regular RiboSeq data.')
    tis2ribo = True
  if len(args.tisbampaths) == 1 : 
    if args.inestpath is None : 
      path = args.tisbampaths[0] + '.bgest.txt'
      if isfile(path) : args.inestpath = path
      else : args.estpath = path
  if args.agenepath is None : args.agenepath = args.genepath

  # load genome, fasta file indexing
  if args.verbose : print("{} Loading genome...".format(time.ctime()))
  genome = fa.Fa(args.genomefapath, verbose = args.verbose)

  # TIS background estimation
  if len(args.tisbampaths) == 0 : 
    print('No input TIS data!')
    paras, slp = [(1,0.5)], [1] # No TIS input
  elif args.inestpath is None : #== '' :
    print ("{} Estimating TIS background parameters...".format(time.ctime()))
    if args.verbose : print("TIS background estimation result will be saved to {}".format(args.estpath))
    if args.numProc > 1 : 
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
      paras, slp, data = pool.apply(ribo.estimateTISbg, args=(args.agenepath, args.tisbampaths, args.genomefapath), kwds={'parts': parts, 'offdict': tisoffdict, 'numProc': args.numProc, 'verbose': args.verbose, 'geneformat': args.geneformat, 'harrwidth': harrwidth, 'paired': paired})
      pool.close()
    else : 
      paras, slp, data = ribo.estimateTISbg(args.genepath, args.tisbampaths, args.genomefapath, parts = parts, offdict = tisoffdict, numProc = 1, verbose = verbose, geneformat = args.geneformat, harrwidth = harrwidth, paired = paired)
    estfile = open(args.estpath, 'w')
    for i in range(len(parts)):
      estfile.write("{}\t{}\t{}\t{}\t{}\n".format(paras[i][0], paras[i][1], parts[i], slp[i], data[i]))
    estfile.close()
    
  else : 
    inestfile = open(args.inestpath, 'r')
    paras, slp = [], []
    for l in inestfile:
      lst = l.strip().split('\t')
      paras.append((float(lst[0]), float(lst[1])))
      slp.append(eval(lst[3]))

  if args.inprofile is not None and not isfile(args.inprofile):
    print('inprofile {} not found!'.format(args.inprofile))
    if args.transprofile is None:
      transprofile = args.inprofile

  if args.numProc > 1 : 
    from multiprocessing import Pool
    pool = Pool(processes = args.numProc - 1) 

  cds_regions = {}
  known_tis = {}
  if args.agenepath != args.genepath :
    if verbose : print('Loading CDS annotation...')
    for g in io.geneIter(args.agenepath, fileType = args.geneformat, verbose = args.verbose):
      if genome.get_chrname(g.chr) is None: continue
      if g.chr not in cds_regions :
        cds_regions[g.chr] = {'+':[interval.Interval() for i in range(3)], '-':[interval.Interval() for i in range(3)]}
        known_tis[g.chr] = {'+':{}, '-':{}}
      for t in g.trans:
        cr = interval.cds_region_trans(t)
        for i in range(3) :
          cds_regions[t.chr][t.strand][i].lst += cr[i].lst
      #for t in g.trans :
        tis = t.cds_start(cdna = False)
        if tis is not None : known_tis[t.chr][t.strand][tis] = 1

  inorf = None
  if args.input is not None :
    if verbose : print('Loading candidates...')
    inorf = {}
    infile = open(args.input, 'r')
    for l in infile :
      lst = l.strip().split()
      tid, tis, stop = lst[0], int(lst[1]), int(lst[2])
      #if gfilter is not None and tid not in gfilter : continue
      if tid not in inorf : inorf[tid] = []
      inorf[tid].append([tis, stop])
  inprofile = None
  if args.inprofile is not None :
    if isfile(args.inprofile):
      if verbose : print('Loading transcript profile...')
      inprofile = {}
      for lst in io.splitIter(args.inprofile):
        try : gid, tid, tispf, ribopf = lst[0], lst[1], eval(lst[3]), eval(lst[4])
        except : continue
        if gid not in inprofile : inprofile[gid] = {}
        inprofile[gid][tid] = tispf, ribopf

  print("{} Predicting...".format(time.ctime()))
  profile = exp.Profile()
  title = ['TISGroup', 'TISCounts', 'TISPvalue', 'RiboPvalue', 'RiboPStatus']
  j = [0,0] # total number of ORF/TIS for BH correction
  gene_iter = io.geneIter(args.genepath, fileType = args.geneformat, verbose = args.verbose) # chrs = genome.idx
  para_iter = genePara(gene_iter, inorf, inprofile)
  if args.numProc <= 1 : pred_iter = imap(_pred_gene, para_iter)
  else : 
    #from multiprocessing import Pool
    #pool = Pool(processes = args.numProc - 1)
    pred_iter = pool.imap_unordered(_pred_gene, para_iter, chunksize = 5)
  if transprofile is not None : 
    tpfile = open(transprofile, 'w')
    tpfile.write('Gid\tTid\tSymbol\tTISProf\tRiboProf\n')

  for result in pred_iter:
    es, ji, tpfs, g = result
    j[0] += ji[0]
    j[1] += ji[1]
    for e in es : 
      profile.add_exp(e)
      if verbose >= 2 : print('{} {}'.format(time.ctime(), str(e)))
    if transprofile is not None : 
      for tid in tpfs :
        tpfile.write(io.tabjoin(tid, tpfs[tid])+'\n')
    if g.chr not in cds_regions : 
      cds_regions[g.chr] = {'+':[interval.Interval() for i in range(3)], '-':[interval.Interval() for i in range(3)]}
      known_tis[g.chr] = {'+':{}, '-':{}}
    for t in g.trans:
      cr = interval.cds_region_trans(t)
      for i in range(3) :
        cds_regions[t.chr][t.strand][i].lst += cr[i].lst
    #for t in g.trans : 
      tis = t.cds_start(cdna = False)
      if tis is not None : known_tis[t.chr][t.strand][tis] = 1

  for chr in cds_regions :
    for strand in cds_regions[chr] :
      for i in range(3) :
        cds_regions[chr][strand][i].check()
  print("{} Checking overlap with known CDS..".format(time.ctime()))
  for e in profile : 
    if e.tistype == 0 : continue
    elif e.gtis in known_tis[e.chr][e.strand] : e.id += ':Known'
    elif e.tistype > 1 : # ["5'UTR", "3'UTR", "Inside", "Novel", 'Extended'] :
      #coding_overlap = False
      for i in range(3) :
        its = cds_regions[e.chr][e.strand][i].intersect(e.cr[i]) # e.cr[i].intersect(cds_regions[e.chr][e.strand][i])
        if its.rlen() > 0 : 
          #coding_overlap = True
          e.id += ':CDSFrameOverlap'
          break

  print("{} BH correcting...".format(time.ctime()))
  profile.BHcorrection(2, total = j[1], append = True) # Calculate BH FDR of TIS p value
  profile.BHcorrection(3, total = j[0], append = True) # Frame p value
  i = 1
  if len(tisbampaths) == 0 : i = 0
  profile.BHcorrection(5, total = j[i], append = True) # Calculate BH FDR for Fisher's p value

  outfile = open(args.output,'w')
  s = "Gid\tTid\tSymbol\tGeneType\tGenomePos\tStartCodon\tStart\tStop\tTisType\t"
  s += '\t'.join(title)
  s += '\tFisherPvalue\tTISQvalue\tFrameQvalue\tFisherQvalue\tAALen'
  if seq : s += '\tSeq'
  if aaseq : s += '\tAASeq'
  if blocks: s += '\tBlocks'
  s += '\n'
  outfile.write(s)

  if args.allresult is not None and args.allresult.upper() == 'OFF' : allout = None
  elif args.fsqth == 1: allout = None
  else : 
    if args.allresult is None : 
      lst = args.output.split('.')
      if lst[-1] == 'txt' : args.allresult = args.output[:-4] + '_all.txt' 
      else : args.allresult = args.output + '_all.txt'
    allout = open(args.allresult, 'w')
    allout.write(s)

  for e in profile:
    #if e.q > args.fsqth : continue
    if len(tisbampaths) == 0 : e.data[5], e.data[8] = None, None # No Fisher's
    s = "%s\t%d" % (e, e.length)
    if seq : s += '\t' + e.sq
    if aaseq : s += '\t' + e.aa
    if blocks: s += '\t' + e.blocks
    s += '\n'
    if allout is not None : allout.write(s)
    if e.q <= args.fsqth : outfile.write(s) # "%s\t%d\n" % (e, e.length)) #, e.sq))
  
def check_overlap(e, known_tis, cds_regions) : 
  if e.tistype == 0 : return ''
  elif e.gtis in known_tis : return ':Known'
  elif e.tistype > 1 : # ["5'UTR", "3'UTR", "Inside", "Novel", 'Extended'] :
    for i in range(3) :
      its = e.cr[i].intersect(cds_regions[i])
      if its.rlen() > 0 : return ':CDSFrameOverlap'

def genePara(gene_iter, inorf, inprofile):
  '''Generate parameters (gene, candidates/None, profile) for function _pred_gene()
  '''
  if inorf is not None :
    for g in gene_iter:
      pf = {}
      if inprofile is not None and g.id in inprofile : pf = inprofile[g.id]
      cand = {}
      for t in g.trans:
        if t.id in inorf : cand[t.id] = inorf[t.id] #yield t, inorf[t.id]
      #if len(cand) > 0 : 
      yield g, cand, pf
  else :
    for g in gene_iter: 
      if gfilter is not None :
        if g.id not in gfilter : continue # or t.gid not in gfilter : continue
      pf = {}
      if inprofile is not None and g.id in inprofile : pf = inprofile[g.id]
      yield g, None, pf

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
    
def _pred_gene(ps): ### trans
  '''Main function of ORF prediction in given transcript
  '''
  g, candidates, pf = ps
  es, j = [], [0,0]
  tpfs = {} #trans profiles
  genome = fa.Fa(genomefapath)
  has_tis = len(tisbampaths) > 0
  if genome.get_chrname(g.chr) is None: return es, j, tpfs, g

  load = True
  if candidates is not None :
    if len(candidates) == 0 :  return es, j, tpfs, g
    if len(pf) >= len(candidates) : load = False
  if len(pf) >= len(g.trans) : load = False
  if load :
    tismbl = ribo.multiRiboGene(g, tisbampaths, offdict = tisoffdict, compatible = compatible, mis = compatiblemis, paired = paired)
    ribombl = ribo.multiRiboGene(g, ribobampaths, offdict = riboffdict, compatible = compatible, mis = compatiblemis, paired = paired)

  for t in g.trans:
    #if candidates is not None and t.id not in candidates : continue
    tl = t.cdna_length()
    if tl < ribo.minTransLen : continue # return es, j, tpfs, g ##
    #ttis = ribo.multiRibo(t, tisbampaths, offdict = tisoffdict, compatible = compatible)
    #tribo = ribo.multiRibo(t, ribobampaths, offdict = riboffdict, compatible = compatible)
    if not load :
      if t.id in pf : 
        tispf, ribopf = pf[t.id]
        ttis = ribo.Ribo(t)
        tribo = ribo.Ribo(t)
        if has_tis : ttis.dict2cnts(tispf)
        tribo.dict2cnts(ribopf)
      else : 
        print('Warning: transcript {} {} {} not in input trans profile! '.format(t.gid, t.id, t.symbol))
        continue
    else :
      ttis = ribo.Ribo(t, bamload = tismbl, compatible = compatible, mis = compatiblemis)
      tribo = ribo.Ribo(t, bamload = ribombl, compatible = compatible, mis = compatiblemis)

    score = ttis.abdscore()
    if score is not None: ip = ribo.pidx(score, slp)   
    else: ip = 0
    if verbose >= 2 : print(io.tabjoin(g.id, t.id, ttis.total, tribo.total))
    cds1 = t.cds_start(cdna = True) 
    cds2 = t.cds_stop(cdna = True) 
    tsq = genome.transSeq(t)

    if transprofile is not None:
      tid = '{}\t{}\t{}'.format(t.gid, t.id, t.symbol)
      tpfs[tid] = '{}\t{}'.format(ttis.cnts_dict_str(), tribo.cnts_dict_str())

    if has_tis and tis2ribo : tribo.merge(ttis) ##

  # user provided candidates
    if candidates is not None : 
      if t.id not in candidates : continue
      for tis, stop in candidates[t.id]:
        j[0] += 1
        j[1] += 1
        if has_tis : tp = ttis.tis_test(tis, paras[ip][0], paras[ip][1])
        else : tp = None
        if enrichtest : rp = tribo.enrich_test(tis, stop)
        else : rp = tribo.frame_test(tis, stop)
        if tp is not None and tp > tpth : continue 
        if rp > fpth : continue # or fisher > fspth 
        minp = rp
        if tp is not None and tp < minp : minp = tp
        if minp > minpth : continue
        fsp, fss = stat.fisher_method([tp, rp]) #
        if tp is not None and fsp > fspth : continue
        has_stop = tsq[stop-3:stop] in orf.cstop
        e = getResult(t, tis, stop, cds1, cds2, tsq, [ip, ttis.cnts[tis], tp, rp, 'N', fsp], has_stop)
        es.append(e)

    else : #all possible ORFs
      orfs = orf.orflist(tsq, minaalen = minaalen, tail = tl)
      for o in orfs :
        starts = o.starts
        if alt : starts += o.altstarts
        starts.sort()
        if longest and not has_tis : starts = starts[0:1]
        ol = len(starts)
        if ol == 0 : continue
        tps = [None] * ol
        rps = [1] * ol
        if has_tis : allz_tis = max(ttis.cnts[starts[0]:o.stop:3]) == 0 # all zeros
        else : allz_tis = True
        allz_ribo = max(tribo.cnts[starts[0]:o.stop:3]) == 0 # all zeros
        if allz_tis and allz_ribo : continue
        for i, tis in enumerate(starts) : 
          if has_tis : tps[i] = ttis.tis_test(tis, paras[ip][0], paras[ip][1])
          if not allz_ribo : allz_ribo = max(tribo.cnts[tis:o.stop:3]) == 0
          if not allz_ribo :
            if enrichtest : rps[i] = tribo.enrich_test(tis, o.stop)
            else : rps[i] = tribo.frame_test(tis, o.stop)
        rst = pvalStatus(rps)
        for i, tis in enumerate(starts) : 
          if tps[i] is not None and tps[i] > tpth : continue
          if rps[i] > fpth : continue # or fishers[i] > fspth
          minp = rps[i]
          if tps[i] is not None and tps[i] < minp : minp = tps[i]
          if minp > minpth : continue
          if tps[i] is None or tps[i] > minpth :
            if longest : 
              if i > 0 : continue
            else :
              if framelocalbest and rst[i] == 'N' : continue
              if framebest and rst[i][0] != 'T' : continue
          fsp, fss = stat.fisher_method([tps[i], rps[i]]) #
          if fsp > fspth : continue
          e = getResult(t, tis, o.stop, cds1, cds2, tsq, [ip, ttis.cnts[tis], tps[i], rps[i], rst[i], fsp], o.has_stop_codon)
          #tistype = tisType(tis, o.stop, cds1, cds2)
          #orfstr = '{}\t{}\t{}'.format(tsq[tis:tis+3],tis,o.stop)
          #tid = "%s\t%s\t%s\t%s\t%s:%d-%d:%s\t%s\t%s" % (t.gid, t.id, t.symbol, t.genetype, t.chr, t.genome_pos(tis), t.genome_pos(o.stop), t.strand, orfstr, tistype)
          #values = [ip, ttis.cnts[tis], tps[i], rps[i], rst[i]] # , fishers[i]]
          #e = exp.Exp(tid, values)
          #e.length = (o.stop - tis) / 3 - 1
          #e.sq = tsq[tis:o.stop]
          #e.chr, e.strand, e.tistype = t.chr, t.strand, tistype
          #if e.tistype == 'Extended' : e.cr = interval.cds_region_trans(t, tis, tis+3)
          #else : e.cr = interval.cds_region_trans(t, tis, o.stop)
          es.append(e)
        #if has_tis : 
        j[1] += ol
        j[0] += 1

  return es, j, tpfs, g

def getResult(t, tis, stop, cds1, cds2, tsq, values, has_stop = True):
  tistype = tisType(tis, stop, cds1, cds2)
  orfstr = '{}\t{}\t{}'.format(tsq[tis:tis+3], tis, stop)
  gtis, gstop = t.genome_pos(tis), t.genome_pos(stop)
  if t.strand != '-' : genomestr = '{}:{}-{}:{}'.format(t.chr, gtis, gstop, t.strand)
  else : genomestr = '{}:{}-{}:{}'.format(t.chr, gstop, gtis, t.strand)
  tid = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (t.gid, t.id, t.symbol, t.genetype, genomestr, orfstr, TIS_types[tistype])
  e = exp.Exp(tid, values)
  e.length = (stop - tis) / 3 - 1
  if not has_stop : e.length += 1
  if seq : e.sq = tsq[tis:stop]
  if aaseq : e.aa = orf.translate(e.sq)
  if blocks: 
    oi = interval.trans2interval(t, tis, stop)
    e.blocks = ','.join(['{}-{}'.format(r[0],r[1]) for r in oi.lst]) 
  e.chr, e.strand, e.tistype = t.chr, t.strand, tistype
  e.gtis, e.gstop = gtis, gstop
  if tistype == 2 : e.cr = interval.cds_region_trans(t, tis, tis+3) # Extended TIS region
  elif tistype > 2 : e.cr = interval.cds_region_trans(t, tis, stop) # new ORFs
  return e

def pvalStatus(ps) : 
  '''Find global best (T) & local best (L) TISs that correspond to same stop codon, suggesting potential TISs by regular RiboSeq data
  '''
  l = len(ps)
  st = ['N'] * l
  if l <= 0 : return st
  if longest and len(tisbampaths) == 0 : return st
  cp = ps[:] + [1]
  m, c = 1, 1
  for i, p in enumerate(cp): 
    if p is None : cp[i] = 1
    elif p < m : m, c = p, 1
    elif p == m : c += 1
  if m >= 1 : return st
  top = 'T'
  if c > 1 : top = 'T{}'.format(c)
  last, c = 1, 1
  better = True
  for i, p in enumerate(cp): 
    if p == m : st[i] = top
    elif p < last : 
      c = 1
      better = True
    elif p == last : c += 1
    else : 
      if better and last > m : 
        local = 'L'
        if c > 1 : local = 'L{}'.format(c)
        for j in range(i-c, i) : st[j] = local
      better = False
      c = 1
    last = p
  return st

def tisType(start, stop, cds1, cds2):
  '''Type of the TIS relatve to annotated ORF
  '''
  if cds1 is None : return 6 # 'Novel'
  else :
    if start == cds1 : return 0 # 'Annotated'
    elif cds1 < stop <= cds2 and (stop - cds2) % 3 == 0 : # in frame
      if start > cds1 : return 1 # 'Truncated'
      else : return 2 # 'Extended'
    else : 
      if start < cds1 : return 3 # "5'UTR"
      elif start > cds2 : return 4 # "3'UTR"
      else : return 5 # "Inside"

if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())


    
