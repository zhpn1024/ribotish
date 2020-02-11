'''
Bam file and reads processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''
import pysam
from . import interval

def changechr(chr):
  '''change between two chr versions
  '''
  if chr.isdigit() or chr in ('X','Y','M'): return 'chr' + chr
  elif chr == 'MT' : return 'chrM'
  elif chr == 'chrM' : return 'MT'
  elif chr[0:3] == 'chr' : return chr[3:]
  else : return chr

chrmap = {}
    
class Bamfile(pysam.Samfile):
  def __repr__(self):
    return 'pysam.Samfile '+self.filename
  def get_chrname(self, chr) : #, chrmap = chrmap):
    #sprint(chrmap)
    if chr in self.references : return chr
    if chr in chrmap :
      chr1 = chrmap[chr]
      if chr1 in self.references : return chr1
      else : 
        chr2 = changechr(chr1)
        if chr2 in self.references : return chr2
    chr1 = changechr(chr)
    if chr1 == chr : return None
    if chr1 in self.references : return chr1
    elif chr1 in chrmap : 
      chr2 = chrmap[chr1]
      if chr2 in self.references : return chr2
    return None
  def fetch_reads(self, chr, start, stop, maxNH = None, minMapQ = None, secondary = False, paired = False) : #, multiple_iterators=False):
    chr0, chr = chr, self.get_chrname(chr)
    #if chr not in self.references : 
      #chr = changechr(chr)
    if chr is None : 
      print("chr {} not found in Bamfile!".format(chr0))
      return #raise StopIteration
    rds = self.fetch(reference=chr, start=start, end=stop) #, multiple_iterators=multiple_iterators)
    r1, r2 = {}, {}
    for read in rds:
      try: 
        if maxNH is not None and read.get_tag('NH') > maxNH : continue
      except: pass
      if not secondary and read.is_secondary : continue
      if minMapQ is not None and read.mapping_quality < minMapQ : continue
      if paired and read.is_paired :
        id = read.query_name
        if id[-2:] in ('\1', '\2', '#1', '#2') : id = id[:-2] ##
        if read.is_read1 : 
          if id in r2 : 
            yield Bam(read, self, read2 = r2[id])
            del r2[id]
          else : r1[id] = read
        elif read.is_read2 : 
          if id in r1 : 
            yield Bam(r1[id], self, read2 = read)
            del r1[id]
          else : r2[id] = read
      else : yield Bam(read, self)
      #yield r
    if paired :
      for id in r1 : 
        if r1[id].mate_is_unmapped : continue
        if r1[id].next_reference_id != r1[id].reference_id : continue
        if r1[id].is_reverse != (not r1[id].mate_is_reverse) : continue
        if abs(r1[id].next_reference_start - r1[id].reference_start) > 1e6 : continue
        #print(id, 'r1')
        read = self.mate(r1[id])
        yield Bam(r1[id], self, read2 = read)
    #if len(r2) > 0 :
      for id in r2 : 
        if r2[id].mate_is_unmapped : continue
        if r2[id].next_reference_id != r2[id].reference_id : continue
        if r2[id].is_reverse != (not r2[id].mate_is_reverse) : continue
        if abs(r2[id].next_reference_start - r2[id].reference_start) > 1e6 : continue
        #print(id, 'r2')
        read = self.mate(r2[id])
        yield Bam(read, self, read2 = r2[id])

class Bam():#AlignedRead
  
  def __init__(self, read, bamfile, read2 = None):
    self.read = read
    self.ref = bamfile.references
    self.lens = bamfile.lengths
    self.read2 = read2
  #def __getattr__(self,attr):
    #a=getattr(self,attr)
    #return a
  
  @property
  def chr(self): 
    #print self.ref['chr1']
    return self.ref[self.read.tid]
  @property
  def start(self):
    return self.read.reference_start
  @property
  def fragment_start(self):
    if self.read2 is None : return self.read.reference_start
    else :
      if self.is_reverse() : return self.read2.reference_start
      else : return self.read.reference_start
    #except : return self.read.pos
  @property
  def stop(self):
    return self.read.reference_end
    #return int(self.read.reference_end)
  @property
  def fragment_stop(self):
    if self.read2 is None : return self.read.reference_end
    else :
      if self.is_reverse() : return self.read.reference_end
      else : return self.read2.reference_end
  @property
  def id(self): 
    return self.read.query_name
  @property
  def score(self): 
    try:
      return self.read.get_tag("AS")
    except:
      return 0.0
  @property
  def strand(self): 
    if self.read.is_reverse: return '-'
    else: return '+'
  @property
  def cigar(self): 
    try : return self.read.cigartuples
    except : return self.read.cigar
  @property
  def weight(self): 
    return 1
  def __len__(self): 
    return self.read.reference_length
  @property
  def length(self):
    return len(self)
  def fragment_length(self): 
    if self.read2 is None : return self.read.reference_length
    else : return abs(self.read.template_length)
  def __repr__(self):
    if self.read2 is None : return "Bam AlignedSegment {} {}".format(self.id, str(self))
    else : return "Bam paired AlignedSegment {} {}".format(self.id, str(self))
  def __str__(self):
    return "%s:%d-%d:%s" % (self.chr, self.fragment_start, self.fragment_stop, self.strand)

  def get_tag(self, k):
    return self.read.get_tag(k)
  def cdna_length(self): 
    return self.read.query_alignment_length
    #else : return self.fragment_length()

  def fragment_length(self):
    if self.read2 is None : return self.cdna_length()
    r1 = interval.Interval(itvs = self.read.get_blocks())
    r2 = interval.Interval(itvs = self.read2.get_blocks())
    ri = r2.intersect(r1)
    ril = ri.rlen()
    if ril > 0 : 
      if not self.is_reverse() :
        r1s = r1.intersect(interval.Interval(r1.start, ri.stop))
        r2s = r2.intersect(interval.Interval(ri.stop, r2.stop))
      else : 
        r1s = r1.intersect(interval.Interval(ri.start, r1.stop))
        r2s = r2.intersect(interval.Interval(r2.start, ri.start))
      return r1s.rlen() + r2s.rlen()
    else : # there may be problems here
      l = max(r1.start, r2.start) - min(r1.stop, r2.stop) # uncovered length, may be not accurate
      if l < 0 : l = 0 ###
      l += self.read.query_alignment_length + self.read2.query_alignment_length 
      
      return l
  def center(self): #Middle point, NEED REVISE!!
    return (self.start+self.stop)/2.0
  
  def is_reverse(self): 
    return self.read.is_reverse
  def is_paired(self):
    return self.read2 is not None
  def is_pair_gapped(self):
    if self.read2 is None : return False
    if not self.is_reverse : return self.read.reference_end < self.read2.reference_start
    else : return self.read2.reference_end < self.read.reference_start
  @property
  def end5(self): #5' end
    if self.read.is_reverse: return self.stop
    else : return self.start
  @property
  def end3(self): #3' end
    if self.read.is_reverse: return self.start
    else : return self.stop
  @property
  def fragment_end5(self): 
    if self.read2 is None : return self.end5
    if self.is_reverse: return self.fragment_stop
    else : return self.fragment_start
  @property
  def fragment_end3(self): 
    if self.read2 is None : return self.end3
    if self.is_reverse: return self.fragment_start
    else : return self.fragment_stop
  def __getattr__(self, name):
    if name == '__class__': return str
    else: return self.read.__getattribute__(name)
    #try : return getattr(self, name)
    #except AttributeError : return getattr(self.read, name)
  def cdna_pos(self, p, strict = False): 
    '''NOT READY because of insertions and deletions
    '''
    if p < self.start or p > self.stop:
      return None
    #p1 = p - self.start
    pos = 0
    for e in self.exons:
      if e.start <= p <= e.stop:
        if self.is_reverse():
          pos += e.stop - p
        else:
          pos += p - e.start
        return pos
      else:
        pos += len(e)
    return None
  
  def genome_pos(self, p, bias = 1):
    '''if bias is 1, the splice junction will be mapped to 5' end of downstream exon,
    if bias is 0, the splice junction will be mapped to 3' end of the upstream exon.
    Do not support paired end read2
    '''
    if p is None : return None
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    p1 = p
    if self.is_reverse():
      p1 = self.cdna_length() - p
      if bias == 0: bias = 1
      else: bias = 0
    pos = self.start
    for ctype, l in self.cigar:
      if ctype in [0,7,8]: # M: alignment match, =: sequence match, X: sequence mismatch
        if l - p1 >= bias:
          pos += p1
          return pos
        else:
          pos += l
          p1 -= l
      elif ctype in [1,]:#4]: # I: insertion to the reference, S: soft clipping
        if l - p1 >= bias: return pos
        else: p1 -= l
      elif ctype in [2,3]: # D: deletion from the reference, N: skipped region from the reference
        pos += l
    return pos
  @property
  def introns(self):
    from . import bed
    s = []
    p = 0
    pos = self.start
    i = 1
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        pos += l
        p += l
      elif ctype in [1]:#,4]:
        p += l
      elif ctype in [2]:
        pos += l
      elif ctype in [3]:
        b = bed.bed6([self.chr,pos,pos + l,self.id+"_intron"+str(i),p,self.strand])
        s.append(b)
        pos += l
        i += 1
    if self.is_reverse():
      return s[::-1]
    return s
  @property
  def exons(self): # NOT TESTED!! Insertions and deletions are skipped
    from . import bed
    s = []
    p = 0
    pos, d = self.start, 0
    i = 0
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        d += l
        p += l
      elif ctype in [1]:#,4]:
        p += l
      elif ctype in [2]:
        d += l
      elif ctype in [3]:
        b = bed.bed6([self.chr,pos,pos + d,self.id+"_exon"+str(i),p,self.strand])
        s.append(b)
        pos += d + l
        i += 1
    b = bed.bed6([self.chr,pos,pos + d,self.id+"_exon"+str(i),p,self.strand])
    s.append(b)
    if self.is_reverse():
      return s[::-1]
    return s
  def isCompatible(self, transitv, mis = 0):
    '''if compatible with given transcript, allow some (mis) incompatible bases
    '''
    if transitv.__class__.__name__ != 'Interval' : 
      transitv = interval.trans2interval(transitv)
    if self.read2 is None : return self.readCompatible(transitv, mis = mis)
    else : return self.readCompatible(transitv, mis = mis) and self.readCompatible(transitv, read = self.read2, mis = mis)
  def readCompatible(self, transitv, read = None, mis = 0):
    if read is None : read = self.read
    ri = interval.Interval(itvs = read.get_blocks())
    return ri.is_compatible(transitv, mis = mis)
    #ti = ri.intersect(transitv)
    #if ri.rlen() - ti.rlen() <= mis : return True
    #return False
    
  def isCompatible0(self, transitv, mis = 0):
    
    #from . import interval
    if transitv.__class__.__name__ != 'Interval' : 
      transitv = interval.trans2interval(transitv)
    r1 = self.start
    r2 = r1 + sum([l for ctype, l in self.cigar])
    ri = interval.Interval(r1, r2) # full reads region
    ti = ri.intersect(transitv) #interval.trans2interval(trans)
    start, stop = ti.start, ti.stop
    m = 0
    pos = self.start
    for ctype, l in self.cigar:
      if ctype in [0,7,8]: # M: alignment match, =: sequence match, X: sequence mismatch
        p1, p2 = max(pos, start), min(pos + l, stop)
        pos += l
        if p1 < p2 : 
          ci = interval.Interval(p1, p2)
          itv = ci.intersect(transitv)
          m += p2 - p1 - itv.rlen()
      elif ctype in [1]:#,4]: # I: insertion to the reference, S: soft clipping
        if start < pos < stop : m += l
      elif ctype in [2,3]: # D: deletion from the reference, N: skipped region from the reference
        p1, p2 = max(pos, start), min(pos + l, stop)
        pos += l
        if p1 < p2 : 
          ci = interval.Interval(p1, p2)
          itv = ci.intersect(transitv)
          m += itv.rlen()
      #print(ctype, l, m)
      if m > mis : return False
    return True
  def is_compatible(self, trans, mis = 0): # old version, slower
    m = 0
    last = -1
    for i in range(self.cdna_length()):
      pos = self.genome_pos(i, bias = 1)
      if pos < trans.start or pos > trans.stop : continue
      i2 = trans.cdna_pos(pos, strict = True)
      #if i2 is None : m += 1
      #else :
      if i2 is not None :
        if last >= 0 :
          if i2 != last + 1 :
            m += min(i, self.cdna_length() - i, abs(i2 - last - 1))
        else : m += min(i, i2)
        last = i2
      #print(i, i2, last, m)
      if m > mis : return False
    #print self, self.cigar, m
    if m > mis : return False
    else : return True
  def is_compatible0(self, trans = None, introns = [], mis = 0): # incorrect version...
    #if trans != None : exons = trans.exons
    if trans is not None : introns = trans.introns
    m = 0
    ris = self.introns
    if len(ris) == 0 : 
      for intr in introns:
        o = self.read.get_overlap(intr.start, intr.stop)
        if o > 0 : m += o
        if m > mis : return False
    else :
      for i in range(len(introns)):
        if not ris[0].is_upstream(introns[i]) : break
      c = True
      for ri in ris:
        if ri != introns[i] : 
          c = False
          o = self.read.get_overlap(introns[i].start, introns[i].stop)
          p = int(ri.score)
          m += max(o, min(p, self.cdna_length() - p))
          if m > mis : return False
        i += 1
    return True
  def is_inside(self, trans, mis = 0):
    '''only inside, not necessarily compatible
    '''
    l = self.cdna_length()
    al = 0
    for e in trans.exons:
      al += self.read.get_overlap(e.start, e.stop)
    if l - al > mis : return False ## 
    else: return True
  def is_m0(self): 
    '''if mismatch at 5' end
    '''
    if self.strand == '+' : 
      if self.read.get_tag('MD')[0] == '0' : return True # mismatch at 0
      #elif self.cigar[0][0] == 4 : return True
      #else : return False
    elif self.read2 is None :
      if self.read.get_tag('MD')[-1] == '0' : 
        if not self.get_tag('MD')[-2].isdigit() : return True
      #elif self.cigar[-1][0] == 4 : return True
      #else : return False
    else :
      if self.read2.get_tag('MD')[-1] == '0' : 
        if not self.read2.get_tag('MD')[-2].isdigit() : return True
    return False
  def mismatches(self):
    import re
    mm = {}
    md = self.read.get_tag('MD')
    ns = re.split('[ATGC^]+', md)
    ms = re.split('\d+', md)
    l = self.cdna_length()
    p = 0
    for i, n in enumerate(ns):
      n = int(n)
      if ms[i] != '':
        q = p
        p += len(ms[i])
        if ms[i].startswith('^'): p -= 1
        if self.is_reverse(): q = l - p
        mm[q] = ms[i]
      p += n
    return mm


def compatible_bam_iter(bamfile, trans, mis = 0, sense = True, maxNH = None, minMapQ = None, secondary = False, flank = 0): 
  '''compatible version of transReadsIter, slightly different
  '''
  chr = bamfile.get_chrname(trans.chr)
  if chr is None : return #raise StopIteration
    #print("chr {} not found in Bamfile!".format(chr0))
    #raise StopIteration
  #chr = trans.chr
  #if chr not in bamfile.references : 
      #chr = changechr(chr)
      #if chr not in bamfile.references : raise StopIteration
  if flank < 0: flank = 0
  rds = bamfile.fetch(reference=chr, start=trans.start-flank, end=trans.stop+flank)#, multiple_iterators=False)
  #introns = trans.introns
  for r in rds:
    read = Bam(r, bamfile)
    if sense and read.strand != trans.strand:
      #print read.id + " not sense"
      continue
    try: 
      if maxNH is not None and read.read.get_tag('NH') > maxNH : continue
    except: pass
    if not secondary and r.is_secondary : continue
    if minMapQ is not None and read.read.mapping_quality < minMapQ : continue
    o = read.read.get_overlap(trans.start-flank, trans.stop+flank)
    if o < read.cdna_length() - mis: 
      continue
    if read.is_compatible(trans = trans, mis = mis) :
      yield read
    #else:
      #print read.id, b, read.cdna_length()
def transReadsIter(bamfile, trans, compatible = True, mis = 0, sense = True, maxNH = None, minMapQ = None, secondary = False, paired = False, flank = 0):
  '''fetch reads from transcript exons
  '''
  chr = bamfile.get_chrname(trans.chr)
  if chr is None : 
    #print("chr {} not found in Bamfile!".format(trans.chr))
    return #raise StopIteration
  if compatible : 
    from . import interval
    transitv = interval.trans2interval(trans)
  if flank < 0: flank = 0
  used = {}
  for e in trans.exons : 
    start, stop = e.start, e.stop
    if flank > 0:
      if start == trans.start:
        start -= flank
        if start < 1: start = 1
      if stop == trans.stop: stop += flank
    rds = bamfile.fetch_reads(chr=chr, start=start, stop=stop, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, paired = paired)#, multiple_iterators=False)
    for read in rds: #yield read
      #read = Bam(r, bamfile)
      if (read.id, read.fragment_start) in used : continue
      if sense and read.strand != trans.strand : continue
      #try: 
        #if maxNH is not None and read.read.get_tag('NH') > maxNH : continue
      #except: pass
      #if not secondary and read.is_secondary : continue
      #if minMapQ is not None and read.read.mapping_quality < minMapQ : continue
      #o = read.read.get_overlap(trans.start, trans.stop)
      #if o < read.cdna_length() - mis: continue
      if read.fragment_start <= e.start or read.fragment_stop >= e.stop : 
        used[(read.id, read.fragment_start)] = 1 # for reads across exons
      if compatible and not read.isCompatible(transitv = transitv, mis = mis) : continue
      yield read

def end5(r):
  return r.end5

class BamLoad:
  def __init__(self, bampath, transRegions = None, posFunc = end5, numProc = 1, pool = None, maxNH = None, minMapQ = None, secondary = False, verbose = False):
    self.bamloads = {}
    self.bampath = bampath
    if bampath is None : return
    self.regions = transRegions
    #if verbose : print('Loading reads...')
    #import itertools
    if numProc >= 1 and pool is None :
      from multiprocessing import Pool
      pool = Pool(processes = numProc)
    if self.regions is not None :
      paras = [(bampath, chr, self.regions[chr], '.', posFunc, maxNH, minMapQ, secondary) for chr in self.regions]
      #for chr in self.regions:
        #self.bamloads[chr] = BamLoadChr(bampath, chr, self.regions[chr], '.', numProc, pool, posFunc, maxNH, minMapQ, secondary)
    else :
      bamfile = Bamfile(bampath)
      paras = [(bampath, chr, None, '.', posFunc, maxNH, minMapQ, secondary) for chr in bamfile.references]
    #if numProc <= 1 : load_iter = itertools.imap(getBamLoadChr, paras)
    #else : 
      #load_iter = pool.imap_unordered(getBamLoadChr, paras, chunksize = 1)
    loads = pool.map(getBamLoadChr, paras, chunksize = 1)
    #results = []
    #if self.regions is not None :
      #for chr in self.regions :
        #results.append(pool.apply_async(BamLoadChr, (bampath, chr, self.regions[chr], '.', posFunc, maxNH, minMapQ, secondary)))
    #else : 
      #bamfile = Bamfile(bampath)
      #for chr in bamfile.references :
        #results.append(pool.apply_async(BamLoadChr, (bampath, chr, None, '.', posFunc, maxNH, minMapQ, secondary)))
      #n = len(bamfile.references)
      #paras = [bampath]*n, bamfile.references, [None]*n, ['.']*n, [posFunc]*n, [maxNH]*n, [minMapQ]*n, [secondary]*n
      #paras = [(bampath, chr, None, '.', posFunc, maxNH, minMapQ, secondary) for chr in bamfile.references]
    #loads = map(BamLoadChr, paras)
    
    
    #loads = [pool.apply_async(BamLoadChr, (bampath, chr, self.regions[chr], '.', posFunc, maxNH, minMapQ, secondary)) for chr in self.regions]
    #loads = pool.imap_unordered(get_BamLoadChr, paras, chunksize = 1)
    #pool.close()
    #pool.join()
    for blc in loads:
      self.bamloads[blc.chr] = blc
    #pool.close()
    #pool.join()
    #for blc in loads : self.bamloads[blc.chr] = blc

  def merge(self, other):
    for chr in other.bamloads:
      if chr not in self.bamloads: self.bamloads[chr] = BamLoadChr(None, chr)#other.bamloads[chr]
      self.bamloads[chr].merge(other.bamloads[chr])
  def transCounts(self, trans, compatible = False, mis = 2):
    return self.bamloads[trans.chr].transCounts(trans, compatible = compatible, mis = mis)

def getBamLoadChr(args):
  bampath, chr, region, strand, posFunc, maxNH, minMapQ, secondary = args
  return BamLoadChr(bampath, chr, region, strand, posFunc, maxNH, minMapQ, secondary)

class BamLoadChr:
  def __init__(self, bampath, chr, region = None, strand = '.', posFunc = end5, maxNH = None, minMapQ = None, secondary = False, paired = False):
    self.chr = chr
    self.strand = strand
    self.bampath = bampath
    self.region = region
    if strand == '.' : self.data = {'+':{}, '-':{}}
    else : self.data = {strand:{}}
    if bampath is None : return
    bamfile = Bamfile(bampath)
    chr = bamfile.get_chrname(chr)
    if chr is None : #raise StopIteration
    #if chr not in bamfile.references : 
      #chr = changechr(chr)
      #if chr not in bamfile.references : 
      print("chr {} not found in Bamfile!".format(self.chr))
      return
    self.chr = chr
    if region is None : region = [[None, None]]
    for start, stop in region : 
      for r in bamfile.fetch_reads(chr, start, stop, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, paired = paired):
        if r.strand not in self.data : continue
        p = posFunc(r)
        #print p
        if p is None: continue
        if start is not None and p < start : continue
        if stop is not None and p > stop : continue
        if p not in self.data[r.strand] : self.data[r.strand][p] = {}
        if r.read2 is None : key =tuple(r.read.get_blocks()), #r.start, tuple(r.cigar)
        else : key = tuple(r.read.get_blocks()), tuple(r.read2.get_blocks())
        #print key
        if key not in self.data[r.strand][p] : self.data[r.strand][p][key] = 0
        self.data[r.strand][p][key] += 1
  def __repr__(self):
    return 'BamLoad {}:{}:{}'.format(self.bampath, self.chr, self.strand)
  def transCounts(self, trans, compatible = False, mis = 2):
    if self.strand != '.' and self.strand != trans.strand :
      print('Strand not match: {}, {}'.format(self, trans.id))
      return []
    d = self.data[trans.strand]
    if compatible : 
      from . import interval
      transitv = interval.trans2interval(trans)
    if trans.strand != '-' : step = 1
    else : step = -1
    cnts = [0] * trans.cdna_length()
    for e in trans.exons:
      p = e.end5
      i = trans.cdna_pos(p)
      for j in range(len(e)):
        if p in d :
          if compatible :
            for key in d[p]:
              if self._compatible(key, transitv, mis = mis) : cnts[i] += d[p][key]
          else : cnts[i] = sum(d[p].values())
        i += 1
        p += step
    return cnts
  def _compatible(self, key, transitv, mis = 2):
    '''if compatible with given transcript, allow some (mis) incompatible bases
    '''
    for blocks in key : 
      ri = interval.Interval(itvs = blocks)
      if not ri.is_compatible(transitv, mis = mis) : return False
    return True
  def _compatible0(self, key, transitv, mis = 2): # key structure changed
    
    from . import interval
    r1 = key[0]
    r2 = r1 + sum([l for ctype, l in key[1]])
    ri = interval.Interval(r1, r2) # full reads region
    ti = ri.intersect(transitv) #interval.trans2interval(trans)
    start, stop = ti.start, ti.stop
    m = 0
    pos = key[0] # self.start
    for ctype, l in key[1] : #self.cigar:
      if ctype in [0,7,8]: # M: alignment match, =: sequence match, X: sequence mismatch
        p1, p2 = max(pos, start), min(pos + l, stop)
        pos += l
        if p1 < p2 : 
          ci = interval.Interval(p1, p2)
          itv = ci.intersect(ti)
          m += p2 - p1 - itv.rlen()
      elif ctype in [1]:#,4]: # I: insertion to the reference, S: soft clipping
        if start < pos < stop : m += l
      elif ctype in [2,3]: # D: deletion from the reference, N: skipped region from the reference
        p1, p2 = max(pos, start), min(pos + l, stop)
        pos += l
        if p1 < p2 : 
          ci = interval.Interval(p1, p2)
          itv = ci.intersect(ti)
          m += itv.rlen()
      if m > mis : return False
    return True
  def merge(self, other):
    for s in other.data:
      if s not in self.data : self.data[s] = {}
      for p in other.data[s]:
        if p not in self.data[s] : self.data[s][p] = {}
        for key in other.data[s][p]:
          if key not in self.data[s][p] : self.data[s][p][key] = 0
          self.data[s][p][key] += other.data[s][p][key]
        
