'''
Main library for riboseq analysis in zbio/ribotish
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''
import itertools, math
from . import stat, bam, gtf, exp, orf, tools, io, fa, interval
from multiprocessing import Pool

codonSize = 3
maxNH = 5
secondary = False
minMapQ = 1 
minRlen = 6
minTransLen = 50
defOffset = 12
#nhead, ntail = defOffset, 30 - defOffset
nhead = ntail = 0

try: imap = itertools.imap
except: imap = map
try: izip = itertools.izip
except: izip = zip

def offset(read = 30, offdict = None):
  '''get P-site offset by read length, default: 12
  '''
  if offdict is None : return defOffset
  if type(read) == int : length = read # compatible with old versions!
  else : length = read.fragment_length()
  if 'm0' not in offdict or not read.is_m0() : od = offdict
  else : od = offdict['m0']
  if length in od : return od[length]
  else : return None # remove this length

def rstest_mw(x, y): 
  '''rank sum test p value of x > y, Not used because of some bug in mannwhitneyu
  '''
  from scipy.stats import ranksums, mannwhitneyu
  #n1, n2 = len(x), len(y)
  #mu = n1 * n2 / 2.
  st1, p1 = ranksums(x, y)
  try : st, p = mannwhitneyu(x, y)
  except : return 0.5
  if st1 > 0 : return p
  else: return 1 - p ###
def rstest(x, y, delta = 1e-4, show = False): 
  '''Exact rank sum test for intergers with ties
  '''
  rs =  stat.rankSumTiesExact(x, y, show = show)
  p = rs.test(delta=delta)
  return p
class Ribo: 
  '''riboseq profile for a transcript
  '''
  def __init__(self, trans, ribobam = None, bamload = None, offset = offset, offdict = None, compatible = True, mis = 2, downsample = 1.0, seed = 1, saverid = False, paired = False):
    self.length = trans.cdna_length()
    self.nhead, self.ntail = nhead, ntail
    self.trans = trans
    self.rids = {}
    #print('offdict',offdict)
    if bamload is not None : 
      self.cnts = bamload.transCounts(trans, compatible = compatible, mis = mis)
      self.total = sum(self.cnts)
      return
    else :
      self.cnts = [0] * self.length
      self.total = 0
    if ribobam is None : return
    if downsample < 1 : 
      import random
      random.seed(seed)
    for r in bam.transReadsIter(ribobam, trans, compatible = compatible, mis = mis, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, paired = paired) : #ribobam.fetch_reads(trans.chr, trans.start, trans.stop):
      l = r.fragment_length()
      off = offset(r, offdict)
      if off is None: continue
      i = trans.cdna_pos(r.genome_pos(off)) ## Default bias
      if i is None or i >= len(self.cnts) : continue
      if downsample < 1 and random.random() >= downsample : continue ## downsample
      self.cnts[i] += 1
      if saverid : self.rids[r.id] = i
    self.total = sum(self.cnts)
  def cnts_enum(self, start = nhead, stop = None):
    if stop is None : stop = self.length - self.ntail
    for i in range(start, stop):
      yield i, self.cnts[i]
  def merge(self, other):
    for i, c in other.cnts_enum():
      self.cnts[i] += c
    self.total += other.total
    return self.total
  def abdscore(self, start = None, stop = None, norm = 1000):
    '''abundant score, log reads per kilo nt
    '''
    if start is None and stop is None : total, length = self.total, self.length
    else :
      if start is None : start = nhead
      if stop is None : stop = self.length - ntail
      total = sum(self.cnts[start:stop])
      length = stop - start
    if total <= 0 : return None
    if length <= 0 : return None
    #import math
    return math.log(float(total) / length * norm)
  def frame_test(self, start, stop, delta=1e-4, glm = False): 
    '''basic frame test of in frame reads higher than out frame reads in the ORF
    start and stop should be in frame
    '''
    l = (stop - start) // 3
    inarr = [0] * l
    outarr = 2 * l * [0]
    for i in range(l):
      p = start + i * 3
      inarr[i] = self.cnts[p]
      outarr[i*2] = self.cnts[p + 1]
      outarr[i*2 + 1] = self.cnts[p + 2]
    if len(inarr) <= 0 or len(outarr) <= 0 : return 1
    if max(inarr) == 0 : return 1
    if glm : return stat.glmNBTest(inarr, outarr)
    return rstest(inarr, outarr, delta = delta)
  def enrich_test(self, start, stop, ts1 = None, ts2 = None, inframe = True, delta=1e-4, glm = False): 
    '''basic enrich test of start~stop reads higher than other reads in ts1~ts2
    if inframe, only test reads that in the ORF frame, and start, stop should be in frame
    '''
    if ts1 is None : ts1 = self.nhead
    if ts2 is None : ts2 = self.length-self.ntail
    if inframe : step = 3
    else : step = 1
    inarr = self.cnts[start:stop:step]
    outarr = self.cnts[start-step:ts1-1:-step] + self.cnts[stop:ts2:step]
    if len(outarr) <= 0: return None
    if max(inarr) == 0 : return 1
    if glm : return stat.glmNBTest(inarr, outarr)
    return rstest(inarr, outarr, delta = delta)
  def frame_test_region(self, region, frame, show = False, old = False, delta=1e-4, glm = False): 
    '''Frame test in a region
    '''
    inarr, outarr = [], []
    for i in region.num_iter():
      i = int(i)
      #print i
      if i % codonSize == frame : inarr.append(self.cnts[i])
      else : outarr.append(self.cnts[i])
    #print(self.trans.id, len(inarr), len(outarr))
    if len(inarr) <= 0 or len(outarr) <= 0 : return None
    if max(inarr) == 0 : return 1 ## 
    if old : return rstest_mw(inarr, outarr)
    elif glm : return stat.glmNBTest(inarr, outarr)
    return rstest(inarr, outarr, delta = delta, show = show)
  def enrich_test_region(self, r1, r2, frame, delta=1e-4, glm = False, show = False): 
    '''Enrich test between regions
    '''
    inarr, outarr = [], []
    for i in r1.num_iter():
      i = int(i)
      if i % codonSize == frame : inarr.append(self.cnts[i])
    for i in r2.num_iter():
      i = int(i)
      if i % codonSize == frame : outarr.append(self.cnts[i])
    if len(inarr) <= 0 or len(outarr) <= 0 : return None
    if max(inarr) == 0 : return 1
    if glm : return stat.glmNBTest(inarr, outarr)
    p = rstest(inarr, outarr, delta = delta)
    if show : 
      print('inarr: {}'.format(inarr))
      print('outarr: {}'.format(outarr))
    return p
  
  def multi_orf_test(self, orflist, glm = False): 
    '''multiple ORF detection for riboseq data, not used in final version
    '''
    from . import interval
    tid = self.trans.id
    blank = []
    for i in range(codonSize): ## Blank regions in different frames
      b = interval.Interval(start=self.nhead, stop=self.length-self.ntail, id=tid)
      blank.append(b)
    eps, fps = [1] * len(orflist), [1] * len(orflist)
    #indeps = [None] * len(orflist)
    gid = [None] * len(orflist) # Group ID
    grp = []
    grplst = [] # Group list
    # Group ORFs
    for i, o in enumerate(orflist): # o is orf.fixedORF object ##(start, stop)
      found = False
      for j, og in enumerate(grp):
        if o.stop == og.stop : # If belong to the group of og
          if og.start > o.start : grp[j] = o # Represented by the longest ORF
          gid[i] = j
          grplst[j].append(o)
          found = True
          break
      if not found : 
        grp.append(o)
        grplst.append([o])
        gid[i] = len(grp) - 1
    # Calculate ORF representative regions 
    for j, og in enumerate(grp):
      blank[og.frame()].sub_itv([og.start, og.stop]) # not blank in given frame
      grplst[j].sort()
      for k, o in enumerate(grplst[j]):
        if k == 0 : o.prev = None # To compare with upstream ORF in the same group
        else : o.prev = grplst[j][k-1]
        try : o.next = grplst[j][k+1] # Downstream ORF
        except : o.next = None
        #if o.next is None : stop = o.stop
        #else : stop = o.next.start  ## Inside group regions
        o.region = interval.Interval(start = o.start, stop = o.stop)
        o.indr = indep_region(o, blank)
    blankall = blank[0]
    for i in range(1, codonSize):
      blankall = blankall.intersect(blank[i])
    blank.append(blankall)
    show =False
    #if tid == 'ENST00000379198' : show = True
    for i, o in enumerate(orflist):
      #print self.trans.id, o, o.region, blank
      eps[i], fps[i] = self.efpvalues(o, blank, glm = glm, show = show) #calculate enrichment and frame p-values
    return eps, fps
    
  def efpvalues(self, orf, blank, glm = False, show = False): # For multi_orf_test
    if orf.indr.rlen() < minRlen : r1 = orf.region
    else : r1 = orf.indr
    #print(self.trans.id, r1)
    fp = self.frame_test_region(r1, orf.frame(), glm = glm)
    r2 = blank[codonSize]
    if r2.rlen() < minRlen or orf.indr.rlen() < minRlen : r2 = blank[orf.frame()]
    ep = self.enrich_test_region(r1, r2, orf.frame(), glm = glm, show = show)
    if show : 
      print(r1, r2, fp, ep, orf.frame())
      #print self.cnts[0:1000]
      #print self.cnts[1000:2000]
    ''' More complicated tests
    if orf.prev is None : 
      r2 = blank[codonSize]
      if r2.rlen() < minRlen or orf.indr.rlen() < minRlen : r2 = blank[orf.frame()]
      ep = self.enrich_test_region(r1, r2, orf.frame())
    else : 
      ep = None
      op = orf.prev
      while op is not None : # Should be larger than all upstream ORF regions
        if op.indr.rlen() < minRlen : r2 = op.region
        else : r2 = op.indr
        ep1 = self.enrich_test_region(r1, r2, orf.frame())
        if ep1 > ep : ep = ep1
        op = op.prev
    '''
    return ep, fp
  def max_orf(self, orf, blank): 
    '''for ribonly_orf_finder
    '''
    from . import interval
    indr = indep_region(orf, blank)
    mst, mp = None, 1
    for st in orf.sts :
      r = interval.Interval(start=st, stop=orf.st)
      r = r.intersect(indr)
      if r.rlen() < minRlen : continue
      fp = self.frame_test_region(r, orf.frame)
      if fp <= mp : mst, mp = st, fp
    return mst, mp

  def ribonly_orf_finder(self, orflist, alt = False, fpth = 0.05):
    '''find likely TIS by only riboseq data, not used, similar function in pvalStatus()
    However, this function seems better
    '''
    from . import interval
    blank = []
    for i in range(codonSize): ## Blank regions in different frames
      b = interval.Interval(start=self.nhead, stop=self.length-self.ntail)
      blank.append(b)
    os = []
    for o in orflist : # Initialize
      o.st = o.stop
      if alt : o.sts = o.starts + o.altstarts
      else : o.sts = o.starts[:]
      if len(o.sts) == 0 : continue
      o.sts.sort()
      if max(self.cnts[o.sts[0]:o.stop]) == 0 : continue
      o.region = interval.Interval(start = o.sts[0], stop = o.stop)
      o.frame = o.stop % 3
      o.fps = {}
      os.append(o)
    findnew = True
    lastaff = interval.Interval(start=self.nhead, stop=self.length-self.ntail)
    while findnew : # Search ORF
      findnew = False
      affected = interval.Interval()
      for o in os : 
        if len(o.sts) == 0 : continue
        if max(self.cnts[o.sts[0]:o.st]) == 0 : continue
        if lastaff.intersect(interval.Interval(o.sts[0], o.st)).rlen() == 0 : continue # Not affected
        lastst = o.st
        st, p = self.max_orf(o, blank)
        while p < fpth : # Max ORF length
          findnew = True
          o.fps[st] = p
          o.st = st
          o.sts = [s for s in o.sts if s < o.st]
          if len(o.sts) == 0 : break
          st, p = self.max_orf(o, blank)
        if lastst != o.st : 
          blank[o.frame].sub_itv([o.st, lastst])
          affected.add_itv([o.st, lastst])
      lastaff = affected
    results = []
    for o in orflist : # Summarize
      if o.st >= o.stop : continue
      o.region = interval.Interval(start=o.st, stop=o.stop)
      region = indep_region(o, blank)
      if region.rlen() < minRlen : region = o.region
      o.fp = min(o.fps.values()) ## min pvalue
      o.ep = self.enrich_test_region(region, blank[o.frame], o.frame)
      results.append(o)
    return results
  def tis_test(self, start, r, p, is_zt = False):
    '''negative binomial test for TIS site
    '''
    if is_zt : nb = stat.ZTNB(r, p)
    else : nb = stat.NegBinom(r, p)
    p = nb.pvalue(self.cnts[start])
    return p
  def is_summit(self, p, flank = 3):
    '''if the count is local max count. not used
    '''
    if p < nhead or p >= self.length - ntail: return False
    if self.cnts[p] == 0: return False
    for i in range(1, flank + 1):
      if p + i < self.length and self.cnts[p] < self.cnts[p + i]: return False
      if p - i >= 0 and self.cnts[p] < self.cnts[p - i]: return False
    return True
  def cnts_dict(self):
    '''only non-sero sites
    '''
    cd = {}
    for i in range(nhead, len(self.cnts) - ntail):
      if self.cnts[i] > 0 : cd[i] = self.cnts[i] 
    return cd
  def cnts_dict_str(self):
    '''string for cnts_dict
    '''
    cdstr = '{'
    for i in range(nhead, len(self.cnts) - ntail):
      if self.cnts[i] > 0 : cdstr += '{}:{}, '.format(i, self.cnts[i]) #cd[i] = self.cnts[i] 
    cdstr = cdstr.rstrip(', ') + '}'
    return cdstr
  def dict2cnts(self, d):
    '''load dict data into cnts
    '''
    for i in d :
      #if i < nhead : continue
      #if i > len(self.cnts) - ntail : continue
      self.cnts[i] = d[i]
  def top_summits_iter(self, start = None, stop = None, minratio = 0, flank = 3, is_zt = False, pth = 0.05, harr = False): 
    '''generate all possible sites, do not specify total number or test with r, p, not used in final version
    '''
    slist = []
    if start is None : start = nhead
    if stop is None : stop = self.length - ntail
    for i in range(start, stop) : #(len(self.cnts)):
      if self.cnts[i] == 0 : continue
      if not self.is_summit(i, flank = flank): continue
      slist.append((i, self.cnts[i]))
    l = len(slist)
    if l < 1 : return 
    slist.sort(key = lambda x: x[1], reverse = True)
    minc = slist[0][1] * minratio ## min reads cutoff, not used
    for i in range(l):
      if slist[i][1] >= minc : yield slist[i]
      #if self.tis_test(i, paras[ip][0], paras[ip][1]) :pass
    #if p < pth:
  def harrTest(self, start, stop, r, p, is_zt = False, harrwidth = 15):
    '''for Harritonine data, not used
    '''
    stop1 = start + harrwidth * 3
    if stop1 > stop : stop1 = stop
    ps = []
    for i in range(start, stop1, 3): ps.append(self.tis_test(i, r, p))
    #print ps
    return stat.fisher_method(ps)
  def harr_summits_iter(self, r, p, start = None, stop = None, pth = 0.05, harrwidth = 15): pass

#### For multiple ORF detection
def indep_region(orf, blank):
  ''' independent region for different overlapping ORFs, not used
  '''
  of = orf.frame
  r = orf.region
  for i in range(codonSize):
    if i != of : r = r.intersect(blank[i])
  return r

#### For TIS background estimation
def sl_idx(value, lst, parts): # group by number of genes
  l = len(lst) - 1
  for i in range(len(parts)):
    if lst[int(l * parts[i])] >= value : break
  return i
def pidx_uplim(i, lst, parts):
  l = len(lst) - 1
  return lst[int(l * parts[i])]
def pidx(value, lst):
  for i, s in enumerate(lst):
    if s >= value : break
  return i
def _merge_cnts(args): # for estimate_tis_bg_all
  g, bampath, offdict = args
  bamfile = bam.Bamfile(bampath, "rb")
  merge = g.merge_trans()
  ml = merge.cdna_length()
  if ml < minTransLen : return None ##
  mtis = Ribo(merge, bamfile, offdict = offdict, compatible = False)
  if mtis.total == 0 : return None
  score = mtis.abdscore()
  return merge, score, mtis.cnts_dict() #, mcds2

def _est_ztnb(args):
  da, maxqt = args
  maxcnt = da.quantile(maxqt) + 5 ##
  d = {}
  for i in range(1, maxcnt): ## max range 
    if i in da: d[i] = da[i]
  zt = stat.ztnb()
  zt.estimate(d, max_iter=10000)#, nlike=100)
  return zt.r, zt.p

def estimate_tis_bg_all(gtfpath, bampath, genomefapath, parts = [0.25, 0.5, 0.75], offdict = None, whole = True, maxqt = 0.95, skip_tis = True, alt_tis = False, tis_flank = 1, addchr = False, numProc = 1, verbose = False):
  '''estimate TIS background using counts from all regions of the transcript. not used
  '''
  parts.sort()
  if parts[-1] < 1 : parts.append(1)
  gtffile = open(gtfpath, 'r')
  genome = fa.Fa(genomefapath)
  fulldata, genes, sl = {}, [], []
  data = [exp.ReadDict() for i in range(len(parts))]
  
  gene_iter = gtf.gtfGene_iter(gtffile, addchr = addchr, chrs = genome.idx, verbose = verbose)
  para_iter = izip(gene_iter, itertools.repeat(bampath), itertools.repeat(offdict))
  if numProc <= 1 : merge_iter = imap(_merge_cnts, para_iter)
  else : 
    pool = Pool(processes = numProc - 1)
    merge_iter = pool.imap_unordered(_merge_cnts, para_iter, chunksize = 100)
    
  for result in merge_iter:
    if result is None : continue
    merge, score, cnts_dict = result
    sl.append(score)
    fulldata[merge.gid] = cnts_dict, score #, mcds1, mcds2
    genes.append(merge)
  sl.sort()
  if verbose : print('Group data...')
  s, ip = 0, 0
  for score, t in sl:
    s += t
    if s >= total * parts[ip] : # group by number of genes
      slp[ip] = score
      ip += 1
  if ip < len(parts) : slp[ip] = score
  for m in genes:
    cnts, score = fulldata[m.gid]
    ip = sl_idx(score, sl, parts)
    if whole : start = nhead ##
    elif len(mcds2) > 0: start = max(mcds2) ##The last stop codon
    else : continue ##
    ml = m.cdna_length()
    msq = tools.gtf2seq(genome, m)
    for i in range(start, ml - ntail):
      if skip_tis and orf.is_start(msq, i, alt = alt_tis, flank = tis_flank) : continue # i in mcds1: continue
      if i not in cnts: data[ip].record(0)
      else : data[ip].record(cnts[i])
      
  if verbose : print ('Estimate ZTNB parameters...')
  paras = []
  slp = [pidx_uplim(i, sl, parts) for i in range(len(parts))]
  if verbose : print (len(sl), slp)
  if numProc <= 1 : paras = list(map(_est_ztnb, [(da, maxqt) for da in data]))
  else : 
    paras = list(pool.map(_est_ztnb, [(da, maxqt) for da in data]))
    pool.close()
  return paras, slp, data

def multiRibo(trans, bampaths, offset = offset, offdict = None, compatible = True, mis = 2, paired = False): 
  '''Load multiple ribobam files to one object
  '''
  if type(bampaths) == list : 
    bamfiles = [bam.Bamfile(bp, "rb") for bp in bampaths]
    if offdict is None : offdict = [None] * len(bampaths)
  else : 
    bamfiles = [bam.Bamfile(bampaths, "rb")]
    offdict = [offdict]
  mribo = Ribo(trans) # Empty object
  for i in range(len(bamfiles)) :
    mribo.merge(Ribo(trans, bamfiles[i], offset = offset, offdict = offdict[i], compatible = compatible, mis = mis, paired = paired))
  return mribo
def multiRiboGene(gene, bampaths, offdict = None, compatible = True, mis = 2, paired = False): 
  '''Load multiple ribobam files to one object
  '''
  regions = interval.allTransRegions(gene.trans)
  if offdict is None : offdict = [None] * len(bampaths)
  mbl = bam.BamLoadChr(None, gene.chr, strand = gene.strand) # Empty object
  for i in range(len(bampaths)) :
    offunc = lambda r : r.genome_pos(offset(r, offdict[i]))
    #def offunc(r):
      #return offset(r, offdict[i])
    mbl.merge(bam.BamLoadChr(bampaths[i], chr = gene.chr, region = regions[gene.chr], strand = gene.strand, posFunc = offunc, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, paired = paired))
  return mbl

  if type(bampaths) == list : 
    bamfiles = [bam.Bamfile(bp, "rb") for bp in bampaths]
    if offdict is None : offdict = [None] * len(bampaths)
  else : 
    bamfiles = [bam.Bamfile(bampaths, "rb")]
    offdict = [offdict]
  mribo = Ribo(trans) # Empty object
  for i in range(len(bamfiles)) :
    mribo.merge(Ribo(trans, bamfiles[i], offdict = offdict[i], compatible = compatible, mis = mis))
  return mribo

def select(g): 
  '''select trans for TIS background estimation, replaced in gtf and bed
  '''
  maxlen, mt = 0, None
  if g.genetype == 'protein_coding' : 
    for t in g.trans:
      cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
      try : cdslen = cds2 - cds1 
      except : continue
      if cdslen % 3 != 0 : continue
      if cdslen > maxlen : maxlen, mt = cdslen, t
    return mt 
  else : return None ###
def _merge_cnts_cds(args): # updated for bampath list
  '''get TIS CDS background data for each transcript
  gene level score and transcript level inframe counts
  '''
  g, bampaths, offdict = args
  merge = g.merge_trans()
  ml = merge.cdna_length()
  if ml < minTransLen : return None ##
  mtis = multiRibo(merge, bampaths, offdict = offdict, compatible = False)
  if mtis.total == 0 : return None
  score = mtis.abdscore()
  t = select(g)
  if t is None : return None
  tl = t.cdna_length()
  cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
  tis = multiRibo(t, bampaths, offdict = offdict, compatible = False)
  return t, score, mtis.total, tis.cnts[cds1:cds2:3] #, mcds2
def _est_nb_lower(da): # not used
  nb = stat.NegBinom()
  return nb.fit_lower(da)
def _est_nb(da):
  nb = stat.NegBinom()
  nb.estimate(da)
  return nb.r, nb.p


def estimate_tis_bg_inframe(genepath, bampaths, genomefapath, parts = [0.25, 0.5, 0.75], offdict = None, skip_tis = True, alt_tis = True, addchr = False, numProc = 1, verbose = False, harrwidth = None):
  '''estimate TIS background using only ORF inframe reads, support only gtf, not used
  '''
  parts.sort()
  if parts[-1] < 1 : parts.append(1)
  gtffile = open(gtfpath, 'r')
  genome = fa.Fa(genomefapath)
  fulldata, genes, sl = {}, [], []
  data = [exp.ReadDict() for i in range(len(parts))]
  
  gene_iter = gtf.gtfgene_iter(gtffile, addchr = addchr, chrs = genome.idx, verbose = verbose)
  para_iter = izip(gene_iter, itertools.repeat(bampaths), itertools.repeat(offdict))
  if numProc <= 1 : merge_iter = imap(_merge_cnts_cds, para_iter)
  else : 
    pool = Pool(processes = numProc - 1)
    merge_iter = pool.imap_unordered(_merge_cnts_cds, para_iter, chunksize = 20)
  ctotal = 0
  for result in merge_iter:
    if result is None : continue
    t, score, total, cnts = result
    sl.append((score, total))
    ctotal += total
    fulldata[t.id] = cnts, score #, mcds1, mcds2
    genes.append(t)
  sl.sort()
  if verbose : print ('Group data...')
  slp = [None] * len(parts)
  s, ip = 0, 0
  for score, t in sl:
    s += t
    if s >= ctotal * parts[ip] : # group by sum of total reads
      slp[ip] = score
      ip += 1
  if ip < len(parts) : slp[ip] = score
  if verbose : print (slp)
  for t in genes:
    cnts, score = fulldata[t.id]
    ip = pidx(score, slp) #ip = sl_idx(score, sl, parts)
    start = t.cds_start(cdna = True)
    if start is None : print (t.gid, t.id)
    tl = t.cdna_length()
    tsq = tools.trans2seq(genome, t)
    for i, c in enumerate(cnts) : # range(start, tl - ntail):
      if harrwidth is not None and i < harrwidth : continue # Only skip known harr region
      if skip_tis and orf.is_start(tsq, start + 3*i, alt = alt_tis, flank = 0) : continue # i in mcds1: continue
      data[ip].record(c)
      
  if verbose : print ('Estimate NB parameters...')
  if numProc <= 1 : paras = list(map(_est_nb, data))
  else : 
    paras = list(pool.map(_est_nb, data)) #[(da, maxqt) for da in data])
    pool.close()
  return paras, slp, data

def _cdsCounts(args): # updated for bampath list
  '''get TIS CDS in frame background data for each transcript
  For estimateTISbg()
  '''
  t, args2 = args
  bampaths, offdict, genomefapath, harrwidth, skip_tis, alt_tis, paired = args2
  tl = t.cdna_length()
  if tl < minTransLen : return None
  cds1 = t.cds_start(cdna = True)
  if cds1 is None : return None
  cds2 = t.cds_stop(cdna = True)
  if cds2 is None or cds2 == cds1 : return None
  if cds2 < cds1 or (cds2 - cds1) % 3 > 0 : 
    print('Wrong CDS annotation: {} {} {}'.format(t.id, cds1, cds2))
    return None
  tis = multiRibo(t, bampaths, offdict = offdict, compatible = False, paired = paired) # not compatible
  score = tis.abdscore()
  if score is None: return None
  genome = fa.Fa(genomefapath)
  tsq = genome.transSeq(t) #tools.trans2seq(genome, t)
  tdata = exp.ReadDict()
  for i in range(cds1, cds2-3, 3) : # range(start, tl - ntail):
    if harrwidth is not None and i < harrwidth * 3 : continue # Only skip known harr region
    if skip_tis and orf.is_start(tsq, i, alt = alt_tis, flank = 0) : continue # i in mcds1: continue
    tdata.record(tis.cnts[i])
  return t, score, tis.total, tdata #tis.cnts[cds1:cds2-3:3] #, mcds2

def estimateTISbg(genepath, bampaths, genomefapath, parts = [0.25, 0.5, 0.75], offdict = None, skip_tis = True, alt_tis = True, addchr = False, numProc = 1, verbose = False, harrwidth = None, geneformat = 'auto', paired = False):
  '''estimate TIS background using only ORF inframe reads, corrent version
  '''
  parts.sort()
  if parts[-1] < 1 : parts.append(1)
  genome = fa.Fa(genomefapath)
  fulldata, genes, sl = {}, [], []
  l = len(parts)
  data = [exp.ReadDict() for i in range(l)]
  
  trans_iter = io.transSelectIter(genepath, fileType = geneformat, chrs = genome.idx, verbose = verbose)
  args2 = bampaths, offdict, genomefapath, harrwidth, skip_tis, alt_tis, paired
  para_iter = izip(trans_iter, itertools.repeat(args2)) # , itertools.repeat(offdict))
  if numProc <= 1 : merge_iter = imap(_cdsCounts, para_iter)
  else : 
    pool = Pool(processes = numProc - 1)
    merge_iter = pool.imap_unordered(_cdsCounts, para_iter, chunksize = 5)
  ctotal = 0
  for result in merge_iter:
    if result is None : continue
    t, score, total, tdata = result
    sl.append((score, total))
    ctotal += total
    fulldata[t.id] = tdata, score #, mcds1, mcds2
    genes.append(t)
  if len(sl) < len(parts):
    print('Too few TIS background samples! Known CDS annotation should be provided. See `-a` option')
    exit(1)
  sl.sort()
  if verbose : print ('Group data...')
  slp = [None] * l
  s, ip = 0, 0
  for score, t in sl:
    s += t
    if s >= ctotal * parts[ip] : 
      slp[ip] = score
      ip += 1
  if ip < l : slp[ip] = score
  if verbose : print (slp)
  for t in genes:
    tdata, score = fulldata[t.id]
    ip = pidx(score, slp) #ip = sl_idx(score, sl, parts)
    data[ip].merge(tdata)
  if verbose : print ('Estimate NB parameters...')
  for i, d in enumerate(data): # in case of empty parts...
    j = 1
    while d.size() == 0 : # len(d) == 0 :
      if i+j < l : d.merge(data[i+j])
      if i-j >= 0 : d.merge(data[i-j])
      j += 1
  if numProc <= 1 : paras = list(map(_est_nb, data))
  else : 
    paras = list(pool.map(_est_nb, data)) #[(da, maxqt) for da in data])
    pool.close()
  return paras, slp, data


def _selectMAXCDS(g): # old
  maxlen, mt = 0, None
  for t in g.trans:
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
    try : cdslen = cds2 - cds1 
    except : continue
    if cdslen % 3 != 0 : continue
    if cdslen > maxlen : maxlen, mt = cdslen, t
  return mt
def _max_cnts(args): # old
  g, bampath, offdict = args
  bamfile = bam.Bamfile(bampath, "rb")
  t = _selectMAXCDS(g)
  if t is None : return None
  tl = t.cdna_length()
  tis = ribo(t, bamfile, offdict = offdict, compatible = False)
  if tis.total == 0 : return None
  tis.score = tis.abdscore()
  return tis
def estimate_tis_bg(gtfpath, bampath, genomefapath, parts = [0.25, 0.5, 0.75], offdict = None, skip_tis = True, alt_tis = True, addchr = False, numProc = 1, verbose = False): # old
  parts.sort()
  if parts[-1] < 1 : parts.append(1)
  gtffile = open(gtfpath, 'r')
  genome = fa.Fa(genomefapath)
  fulldata, genes, sl = {}, [], []
  data = [exp.readdict() for i in range(len(parts))]
  
  gene_iter = gtf.gtfgene_iter(gtffile, addchr = addchr, chrs = genome.idx, verbose = verbose)
  para_iter = izip(gene_iter, itertools.repeat(bampath), itertools.repeat(offdict))
  if numProc <= 1 : merge_iter = imap(_max_cnts, para_iter)
  else : 
    pool = Pool(processes = numProc - 1)
    merge_iter = pool.imap_unordered(_max_cnts, para_iter, chunksize = 20)
  total = 0
  for result in merge_iter:
    if result is None : continue
    tis = result
    total += tis.total
    sl.append((tis.score, tis.total))
    fulldata[tis.trans.gid] = tis #.cnts, tis.score #, mcds1, mcds2
    genes.append(tis.trans)
  if verbose : print ('total = {}'.format(total))
  if verbose : print ('Sorting...')
  sl.sort()
  if verbose : print ('Group data...')
  s, ip, slp = 0, 0, [None] * len(parts)
  for score, t in sl:
    s += t
    if s >= total * parts[ip] : 
      slp[ip] = score
      ip += 1
  if ip < len(parts) : slp[ip] = score
  if verbose : print (slp)
  for t in genes:
    tis = fulldata[t.gid]
    ip = pidx(tis.score, slp)
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True)
    if cds1 is None : print ('No CDS', t.gid, t.id)
    tl = t.cdna_length()
    tsq = tools.gtf2seq(genome, t)
    for i in range(cds1+3, cds2, 3) :
      codon = tsq[i:i+3].upper()
      if codon in orf.cstart or codon in orf.cstartlike : continue
      data[ip].record(tis.cnts[i])
      
  if verbose : print ('Estimate NB parameters...')
  if numProc <= 1 : paras = list(map(_est_nb, data))
  else : 
    paras = list(pool.map(_est_nb, data)) #[(da, maxqt) for da in data])
    pool.close()
  return paras, slp, data

def frame_bias(arr): # ['0', '1', '2', '01', '02', '12', '012']
  m = max(arr[0:codonSize])
  s = ''
  if m == 0 : return s
  for i in range(codonSize):
    if arr[i] == m : s += str(i)
  return s

################### Quality Control ###########################

class lenDis:
  '''Calculate several distributions for quality control
  '''
  def __init__(self, lens, dis, tl = 0, cds1 = 0, cds2 = 0, offset = defOffset, flank = 0):
    d = dis[1] - dis[0]
    self.dis = dis
    self.lens = lens
    self.tl, self.flank = tl, flank
    self.tl2 = tl + flank * 2
    self.cds1, self.cds2 = cds1, cds2
    self.d1, self.d2 = {}, {} #sum of reads near start or stop condons
    self.d1d, self.d2d = {}, {} # distributions of d1, d2
    self.df = {} # distribution of frame
    self.dc = {}
    self.l = {}
    self.cnts = {}
    self.offset = offset
    for l in range(lens[0], lens[1]):
      self.l[l] = 0
      self.d1[l] = [0] * d
      self.d2[l] = [0] * d
      self.d1d[l] = [exp.ReadDict() for i in range(d)]
      self.d2d[l] = [exp.ReadDict() for i in range(d)]
      self.cnts[l] = [0] * self.tl2 # for transcript level, adding flanking positions
      self.df[l] = [exp.ReadDict() for i in range(codonSize)]
  def record(self, l, i, n = 1):
    '''record a read with length l and 5' end position i
    '''
    self.l[l] += n
    i2 = i + self.flank
    if 0 <= i2 < self.tl2: self.cnts[l][i2] += n
    ir = i - self.cds1
    if self.dis[0] <= ir < self.dis[1] : self.d1[l][ir - self.dis[0]] += n
    ir = i - self.cds2
    if self.dis[0] <= ir < self.dis[1] : self.d2[l][ir - self.dis[0]] += n
  def disFrame(self):
    '''calculate frame distribution
    '''
    for l in self.d1: 
      self.df[l] = [exp.ReadDict() for i in range(codonSize)] # reset
      for i in range(self.cds1, self.cds2 + codonSize, codonSize):
        io = i - self.offset - codonSize + self.flank # -15 -> stop
        if io < 0 : continue
        for i2 in range(codonSize):
          self.df[l][i2].record(self.cnts[l][io+i2]) #
    return self.df
  def disCDS(self, bins = 20):
    '''calculate CDS profile
    '''
    bins = int(bins)
    if bins < 2 : bins = 2
    self.bins = bins
    length = (self.cds2 - self.cds1) // codonSize + 1 # One more codon 
    binsize = float(length) / bins
    for l in self.d1 : 
      self.dc[l] = [[0.0] * bins for i in range(codonSize)]
    if binsize == 0 : return self.dc
    
    for l in self.d1 : 
      cfcnts = [[0] * length for j in range(codonSize)] # profile in 3 frames
      for si in range(length) :
        for j in range(codonSize) : 
          io = si * codonSize + j + self.cds1 - self.offset - codonSize + self.flank
          if io < 0 : continue
          cfcnts[j][si] = self.cnts[l][io]
      for i in range(bins):
        start = binsize * i
        stop = start + binsize
        si1, si2 = int(start), int(stop)
        if si1 == si2 : # within one base/codon
          for j in range(codonSize) : self.dc[l][j][i] = cfcnts[j][si1]
        else : # start, stop in different bases
          for j in range(codonSize) : 
            if si1 < start : self.dc[l][j][i] += cfcnts[j][si1] * (si1 + 1 - start) # start is between si1 and si1 + 1
            else : self.dc[l][j][i] += cfcnts[j][si1]
            for si in range(si1 + 1, si2) : self.dc[l][j][i] += cfcnts[j][si]
            if i < bins - 1 and si2 < stop : self.dc[l][j][i] += cfcnts[j][si2] * (stop - si2)
            if i == bins - 1 and si2 < length : self.dc[l][j][i] += cfcnts[j][si2]
            self.dc[l][j][i] /= binsize
    return self.dc
  def merge(self, other):
    '''merge all data from transcript together
    '''
    for l in other.d1: 
      self.l[l] += other.l[l]
      for di in range(other.dis[1] - other.dis[0]):
        self.d1d[l][di].record(other.d1[l][di])
        self.d2d[l][di].record(other.d2[l][di])
      if l not in self.df : self.df[l] = [exp.ReadDict() for i in range(codonSize)]
      for i2 in range(codonSize):
        self.df[l][i2].merge(other.df[l][i2]) #
      if l not in self.dc: self.dc[l] = [[0.0] * other.bins for i in range(codonSize)]
      for j in range(codonSize) :
        for i in range(other.bins):
          self.dc[l][j][i] += other.dc[l][j][i]
  def size(self):
    return sum(self.l.values())
  def write(self, outfile = None):
    md1, md2 = {}, {}
    mdf = {}
    for l in self.l: 
      md1[l] = [d.sum() for d in self.d1d[l]] #[0] * len(dis1[l])
      md2[l] = [d.sum() for d in self.d2d[l]] #[0] * len(dis2[l])
      mdf[l] = [d.sum() for d in self.df[l]]
    if outfile is not None : outfile.write('{}\n{}\n{}\n{}\n{}\n'.format(self.l, md1, md2, mdf, self.dc))
    return self.l, md1, md2, mdf, self.dc

def _lendis_gene(args):
  '''quality profile for each gene, not used
  '''
  g, bampath, lens, dis, ccds, minR, m0, cdsBins = args
  bamfile = bam.Bamfile(bampath, "rb")
  maxlen = 0
  mt = None
  for t in g.trans:
    if ccds and t.attr('ccds_id') == '' : continue
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
    try : cdslen = cds2 - cds1 
    except : continue
    if cdslen % 3 != 0 : continue
    if cdslen > maxlen : maxlen, mt = cdslen, t
  if mt is None : return None
  t = mt
  tl = t.cdna_length()
  cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
    
  td = lenDis(lens, dis, tl, cds1, cds2)
  if m0: tdm = lenDis(lens, dis, tl, cds1, cds2)
  
  tr = 0 # Total reads
  for r in bam.compatible_bam_iter(bamfile, t, mis = 2, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary):
    if m0 : ism0 = r.is_m0()
    l = r.fragment_length()
    if l < lens[0] or l >= lens[1]: continue # not in given length range
    i = t.cdna_pos(r.genome_pos(0)) # 5' end 
    if i is None : continue
    tr += 1
    if not m0 or not ism0 : 
      td.record(l, i)
    else : 
      tdm.record(l, i)
  if tr < minR : return None
  td.disFrame()
  td.disCDS(bins = cdsBins)
  td.cnts = {} # release memory
  result = [td]
  if m0 : 
    tdm.disFrame()
    tdm.disCDS(bins = cdsBins)
    tdm.cnts = {}
    result += [tdm]
  return result

def _lendis_trans(args):
  '''quality profile in each transcript
  '''
  t, bampath, lens, dis, ccds, minR, m0, cdsBins, paired = args
  #print(t)
  bamfile = bam.Bamfile(bampath, "rb")
  tl = t.cdna_length()
  cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
  flank = - dis[0] # in case transcript annotation has no 5'UTR (starts from CDS start)
  td = lenDis(lens, dis, tl, cds1, cds2, flank=flank)
  if m0: tdm = lenDis(lens, dis, tl, cds1, cds2, flank=flank)
  tr = 0 # Total reads
  for r in bam.transReadsIter(bamfile, t, compatible=False, maxNH=maxNH, minMapQ=minMapQ, secondary=secondary, paired=paired, flank=flank):
    if m0 : ism0 = r.is_m0()
    l = r.fragment_length()
    if l < lens[0] or l >= lens[1]: continue # not in given length range
    i = t.cdna_pos(r.genome_pos(0), flank=flank) # 5' end 
    if i is None : continue
    tr += 1
    if not m0 or not ism0 : 
      td.record(l, i)
    else : 
      tdm.record(l, i)
  if tr < minR : return None
  td.disFrame()
  td.disCDS(bins = cdsBins)
  td.cnts = {} # release memory
  result = [td]
  if m0 : 
    tdm.disFrame()
    tdm.disCDS(bins = cdsBins)
    tdm.cnts = {}
    result += [tdm]
  return result
def lendis(genepath, bampath, lens = [25,35], dis = [-40,20], ccds = False, minR = 1, m0 = True, cdsBins = 20, numProc = 1, addchr = False, verbose = False, geneformat = 'auto', paired = False):
  '''distributions of different reads length, for quality control
  '''
  #gtffile = open(gtfpath,'r')
  ad = lenDis(lens, dis)
  if m0: adm = lenDis(lens, dis)

  #gene_iter = gtf.gtfgene_iter(gtffile, addchr = addchr, verbose = verbose)
  trans_iter = io.transSelectIter(genepath, fileType = geneformat, verbose = verbose)
  rep = itertools.repeat
  para_iter = izip(trans_iter, rep(bampath), rep(lens), rep(dis), rep(ccds), rep(minR), rep(m0), rep(cdsBins), rep(paired))
  if numProc <= 1 : len_iter = imap(_lendis_trans, para_iter) #_lendis_gene
  else : 
    pool = Pool(processes = numProc - 1)
    len_iter = pool.imap_unordered(_lendis_trans, para_iter, chunksize = 20) #_lendis_gene
    
  for result in len_iter:
    if result is None : continue
    ad.merge(result[0])
    if m0 : adm.merge(result[1])
  if numProc > 1 : pool.close()
  results = [ad]
  if m0 : results += [adm]
  return results


def lendisM0(gtfpath, bampath, lens = [26,35], dis = [-40,20], minR = 1):
  '''old version, do not seperate mismatch at 0
  '''
  bamfile = bam.Bamfile(bampath, "rb")
  gtffile = open(gtfpath,'r')
  d = dis[1] - dis[0]
  dis1, dis2 = {}, {} #sum of reads near start or stop condons
  disf = {} # distribution of frame
  fbias = {} # frame bias
  fbkeys = ['0', '1', '2', '01', '02', '12', '012']
  lendis = {} # reads length distribution
  for l in range(lens[0], lens[1]):
    lendis[l] = 0
    dis1[l] = [None] * d
    dis2[l] = [None] * d
    disf[l] = [None] * codonSize
    fbias[l] = {}
    for di in range(d): 
      dis1[l][di] = exp.readdict()
      dis2[l][di] = exp.readdict()
    for i in range(codonSize):
      disf[l][i] = exp.readdict()
    for k in fbkeys:
      fbias[l][k] = 0
  for g in gtf.gtfgene_iter(gtffile):
    maxlen = 0
    mt = None
    for t in g.trans:
      cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
      try : cdslen = cds2 - cds1 
      except : continue
      if cdslen % 3 != 0 : continue
      if cdslen > maxlen : maxlen, mt = cdslen, t
    if mt is None : continue
    t = mt
    tl = t.cdna_length()
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
    tdis1, tdis2 = {}, {}
    cnts = {}
    for l in range(lens[0], lens[1]): 
      tdis1[l] = [0] * d
      tdis2[l] = [0] * d
      cnts[l] = [0] * tl
    tr = 0
    
    for r in bam.compatible_bam_iter(bamfile, t, mis = 2, maxNH = maxNH, minMapQ = minMapQ):
      if r.get_tag('MD')[0] != '0' : continue ## Only reads with mismatch at 1st nt
      l = r.fragment_length()
      if l not in tdis1: continue
      i = t.cdna_pos(r.genome_pos(0)) # 5' end 
      if i is None : continue
      tr += 1
      lendis[l] += 1
      cnts[l][i] += 1
      ir = i - cds1
      if dis[0] <= ir < dis[1] : tdis1[l][ir - dis[0]] += 1
      ir = i - cds2
      if dis[0] <= ir < dis[1] : tdis2[l][ir - dis[0]] += 1
    if tr < minR : continue
    for l in tdis1: 
      for di in range(d): 
        dis1[l][di].record(tdis1[l][di])
        dis2[l][di].record(tdis2[l][di])
      for i in range(cds1, cds2, codonSize):
        io = i-defOffset
        if io < 0 : continue
        for i2 in range(codonSize):
          disf[l][i2].record(cnts[l][io+i2]) #
        s = frame_bias(cnts[l][io:io+codonSize])
        if s != '' : fbias[l][s] += 1
  return lendis, dis1, dis2, disf, fbias

def quality(arr, threshold = 0.5, comment = False): 
  '''Quality estimation by RPF frame distribution
  '''
  good = 0.7
  m = max(arr)
  for i, a in enumerate(arr):
    if a == m : 
      frame = i
      break
  txt = "%.2f " % (m)
  use = True
  if m < threshold : 
    if comment : txt += 'fail'
    use = False
  elif comment : 
    if m < good : txt += 'pass'
    else : txt += 'good'
  return use, frame, txt

def getTIS(arr, dis = [-40,20], defOffset = defOffset) :
  '''get position of max count for TIS data
  '''
  m = max(arr)
  mis = [ i + dis[0] for i, n in enumerate(arr) if n == m]
  md, mp = 40, -defOffset
  for i in mis : 
    if abs(i + defOffset) < md : 
      md, mp = abs(i + defOffset), i
  return mp
def TISquality(arr, dis = [-40,20], defOffset = defOffset, flank = 3, threshold = 0.5, comment = False) :
  '''quality for TIS data
  '''
  good = 0.7
  mp = getTIS(arr, dis, defOffset)
  frame = mp % 3
  i = mp - dis[0]
  y0 = arr[i-1:i+2]
  ys = sum(y0)
  if ys > 0 : m = float(arr[i])/ys
  else : m = 0
  txt = " %.2f " % (m)
  use = True
  if abs(mp + defOffset) > flank : use = False
  elif m < threshold : use = False
  elif comment :
    if m < good : txt += 'pass'
    else : txt += 'good'
  return use, frame, txt, mp
# 
def get_offset(arr, dis = [-40,20], frame = 0, defOffset = defOffset, flank = 6, tis = False): 
  ''' Estimate RPF P site offset distance, return offset and threshold
  '''
  a0 = [x for i, x in enumerate(arr) if (i + dis[0]) % codonSize == frame and i + dis[0] <= - defOffset - flank]
  a1 = [x for i, x in enumerate(arr) if (i + dis[0]) % codonSize == frame and i + dis[0] >= - defOffset - flank]
  ai = [x for i, x in enumerate(arr) if (i + dis[0]) % codonSize == frame and i + dis[0] > - defOffset - flank and i + dis[0] < - defOffset + flank]
  #a0.sort()
  #a1.sort()
  #a0m = 1.0 * sum(a0) / len(a0)
  #a1m = 1.0 * sum(a1) / len(a1)
  a0m = max(a0) # / 2.0
  a1m = max(a1) # / 2.0
  aim = max(ai)
  #th = int(((a0m + 1) * (a1m + 1)) ** 0.5)
  th = a0m + int((a1m - a0m) / 6.0) ###
  for p in range(-defOffset - flank + 1, -defOffset + flank):
    if p % codonSize != frame : continue
    if tis :
      if arr[p - dis[0]] == aim : return -p, th
    elif arr[p - dis[0]] > th : return -p, th
  return None, th


def write_off_para(parafile, offdict):
  '''Generate python code of offset function for parameter file
  '''
  parafile.write('offdict = {}\n'.format(formatdict(offdict)))
def formatdict(d, tab=1):
  '''show dict in ordered manner
  '''
  s = '{'
  ks = sorted([k for k in d if type(k) == int]) + sorted([k for k in d if type(k) == str])
  for k in ks:
    if type(d[k]) == dict : s += '{}: {}, '.format(repr(k), formatdict(d[k], tab+1))
    else : s += '{}: {}, '.format(repr(k), repr(d[k]))
  s = s.strip(', ') + '}'
  return s

def TIStest_betaBinom(t1, t2, r1, r2, scale_t = 1, scale_r = 1, alt = 'two.tailed', prior = [1,1]):
  ''' Differental TIS test using Beta-Binomial model.
      t1 & t2 are TIS raw counts. r1 & r2 are RNASeq raw counts.
      scale_t & scale_r are scale factors of 2/1
  '''
  st2 = scale_t ** 0.5
  sr2 = scale_r ** 0.5
  alpha = (prior[0] + r1 * sr2) / st2
  beta = (prior[1] + r2 / sr2) * st2
  bb = stat.betaBinom(alpha, beta)
  #if alpha == beta == 0.0 : return 1
  return bb.pvalue(t1+t2, t1, alt = alt)
def TIStest_chi2(t1, t2, r1, r2, scale_t = 1, scale_r = 1, alt = 'two.tailed'):
  ''' Differental TIS test using Chi square test
      Parameters are same as TIStest_betaBinom
  '''
  st2, sr2 = scale_t ** 0.5, scale_r ** 0.5
  nt1, nt2 = t1 * st2, t2 / st2
  nr1, nr2 = r1 * sr2, r2 / sr2
  total = nt1 + nt2 + nr1 + nr2
  sum1 = nt1 + nr1 #, nt2 + nr2
  sumt = nt1 + nt2 #, nr1 + nr2
  p1, pt = sum1 / total, sumt / total
  p2, pr = 1 - p1, 1 - pt #sumt / total, sum2 / total
  obs = nt1, nt2, nr1, nr2
  exp = total * p1 * pt, total * p2 * pt, total * p1 * pr, total * p2 * pr
  chi2, pv = stat.chisquare(obs, exp, ddof=1)
  if math.isnan(pv) : 
    #print('nan for {} {} {} {}!'.format(t1, t2, r1, r2))
    return 1.0
  if alt == 'two.sided' : return pv
  if obs[0] >= exp[0] and alt in ('g', 'greater'): return pv / 2
  if obs[0] <= exp[0] and alt in ('l', 'less'): return pv / 2
  return 1 - pv / 2
def TIStest_FisherExact(t1, t2, r1, r2, scale_t = 1, scale_r = 1, alt = 'two.tailed') :
  ''' Differental TIS test using Fisher's exact test of rounded normalized counts
  '''
  st2, sr2 = scale_t ** 0.5, scale_r ** 0.5
  if t1 + t2 > r1 + r2 :
    nt1, nt2 = int(round(t1 * st2 / sr2)), int(round(t2 / st2 * sr2))
    nr1, nr2 = r1, r2
  else :
    nt1, nt2 = t1, t2
    nr1, nr2 = int(round(r1 * sr2 / st2)), int(round(r2 / sr2 * st2))
  total = nt1 + nt2 + nr1 + nr2
  sum1, sum2 = nt1 + nr1, nt2 + nr2
  sumt, sumr = nt1 + nt2, nr1 + nr2
  #p1, pt = sum1 / total, sumt / total
  #p2, pr = 1 - p1, 1 - pt 
  #exp = total * p1 * pt
  pv =  stat.hypergeo_test(total, sum1, sumt, nt1, alt)
  return pv
  #if alt == 'two.sided' : return pv
  #if nt1 >= exp and alt in ('g', 'greater'): return pv
  #if nt1 <= exp and alt in ('l', 'less'): return pv
  #return 1 - pv / 2

