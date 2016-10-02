import pysam
import bed

def changechr(chr):
  '''change between two chr versions
  '''
  if chr.isdigit() or chr in ('X','Y','M'): return 'chr' + chr
  elif chr == 'MT' : return 'chrM'
  elif chr == 'chrM' : return 'MT'
  elif chr[0:3] == 'chr' : return chr[3:]
  else : return chr

class Bamfile(pysam.Samfile):
  def __repr__(self):
    return 'pysam.Samfile '+self.filename
  def fetch_reads(self, chr, start, stop) : #, multiple_iterators=False):
    if chr not in self.references : 
      chr = changechr(chr)
      if chr not in self.references : raise StopIteration
    rds = self.fetch(reference=chr, start=start, end=stop) #, multiple_iterators=multiple_iterators)
    for read in rds:
      r = Bam(read, self)
      yield r

class Bam():#AlignedRead
  
  def __init__(self, read, bamfile):
    self.read = read
    self.ref = bamfile.references
    self.lens = bamfile.lengths
  #def __getattr__(self,attr):
    #a=getattr(self,attr)
    #return a
  
  @property
  def chr(self): 
    #print self.ref['chr1']
    return self.ref[self.read.tid]
  @property
  def start(self):
    return self.read.pos
  @property
  def stop(self):
    return int(self.read.aend)
  @property
  def id(self): 
    return self.read.qname
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
    return self.read.cigar
  
  def __len__(self): #All bed
    return self.read.alen
  @property
  def length(self):
    return len(self)
  
  def __repr__(self):
    return "Bam AlignedSegment " + self.id
  def __str__(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  
  def cdna_length(self): 
    return self.read.qlen
  
  def center(self): #Middle point, NEED REVISE!!
    return (self.start+self.stop)/2.0
  
  def is_reverse(self): 
    return self.read.is_reverse
  
  @property
  def end5(self): #5' end, all bed
    if self.read.is_reverse: return self.stop
    else : return self.start
  @property
  def end3(self): #3' end, all bed
    if self.read.is_reverse: return self.start
    else : return self.stop
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
    '''
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
      elif ctype in [1,4]: # I: insertion to the reference, S: soft clipping
        if l - p1 >= bias: return pos
        else: p1 -= l
      elif ctype in [2,3]: # D: deletion from the reference, N: skipped region from the reference
        pos += l
    return None
  @property
  def introns(self):
    s = []
    p = 0
    pos = self.start
    i = 1
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        pos += l
        p += l
      elif ctype in [1,4]:
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
    s = []
    p = 0
    pos, d = self.start, 0
    i = 0
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        d += l
        p += l
      elif ctype in [1,4]:
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
  def isCompatible(self, trans, mis = 0):
    '''if compatible with given transcript, allow some (mis) incompatible bases
    '''
    from . import interval
    transitv = interval.trans2interval(trans)
    m = 0
    pos = self.start
    for ctype, l in self.cigar:
      if ctype in [0,7,8]: # M: alignment match, =: sequence match, X: sequence mismatch
        p1, p2 = max(pos, trans.start), min(pos + l, trans.stop)
        pos += l
        if p1 < p2 : 
          ci = interval.Interval(p1, p2)
          itv = ci.intersect(transitv)
          m += p2 - p1 - itv.rlen()
      elif ctype in [1,4]: # I: insertion to the reference, S: soft clipping
        if trans.start < pos < trans.stop : m += l
      elif ctype in [2,3]: # D: deletion from the reference, N: skipped region from the reference
        p1, p2 = max(pos, trans.start), min(pos + l, trans.stop)
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
      if i2 is None : m += 1
      else :
        if last >= 0 :
          if i2 != last + 1 :
            m += min(i, self.cdna_length() - i, abs(i2 - last - 1))
        else : pass # m += min(i, i2)
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
      if self.get_tag('MD')[0] == '0' : return True # mismatch at 0
      else : return False
    elif self.get_tag('MD')[-1] == '0' : 
      if not self.get_tag('MD')[-2].isdigit() : return True
    return False
def compatible_bam_iter(bamfile, trans, mis = 0, sense = True, maxNH = None, minMapQ = None, secondary = False): 
  '''compatible version of transReadsIter, slightly different
  '''
  chr = trans.chr
  if chr not in bamfile.references : 
      chr = changechr(chr)
      if chr not in bamfile.references : raise StopIteration
  rds = bamfile.fetch(reference=chr, start=trans.start, end=trans.stop)#, multiple_iterators=False)
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
    o = read.read.get_overlap(trans.start, trans.stop)
    if o < read.cdna_length() - mis: 
      continue
    if read.is_compatible(trans = trans, mis = mis) :
      yield read
    #else:
      #print read.id, b, read.cdna_length()
def transReadsIter(bamfile, trans, compatible = True, mis = 0, sense = True, maxNH = None, minMapQ = None, secondary = False):
  '''fetch reads from transcript exons
  '''
  chr = trans.chr
  if chr not in bamfile.references : 
      chr = changechr(chr)
      if chr not in bamfile.references : raise StopIteration
  used = {}
  #trans.exons = trans.exons ##
  for e in trans.exons : 
    rds = bamfile.fetch_reads(chr=chr, start=e.start, stop=e.stop)#, multiple_iterators=False)
    for read in rds: #yield read
      #read = Bam(r, bamfile)
      if (read.id, read.start) in used : continue
      if sense and read.strand != trans.strand : continue
      try: 
        if maxNH is not None and read.read.get_tag('NH') > maxNH : continue
      except: pass
      if not secondary and read.is_secondary : continue
      if minMapQ is not None and read.read.mapping_quality < minMapQ : continue
      #o = read.read.get_overlap(trans.start, trans.stop)
      #if o < read.cdna_length() - mis: continue
      if read.start <= e.start or read.stop >= e.stop : used[(read.id, read.start)] = 1 # for reads across exons
      if compatible and not read.isCompatible(trans = trans, mis = mis) : continue
      yield read
