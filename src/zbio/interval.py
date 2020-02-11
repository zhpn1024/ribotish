'''
Intervals processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

def is_overlap(i1, i2):
  return i1[0] < i2[1] and i1[1] > i2[0]

class Interval:
  '''intervals in 1 dimentional axis
  all intervals are supposed to be [start, end) and start <= end
  '''
  def __init__(self, start = None, stop = None, id = '', itvs = [], adj_merge=True):
    self.id = id
    self.lst = [] # list of (start, stop) tuples
    self.adj_merge = adj_merge
    if itvs == []:
      if start is not None and start <= stop : self.lst.append((start, stop))
    else :
      self.lst += itvs
      #for itv in itvs:
        #self.lst.append(itv[:])
    self.check()
  def check(self):
    for itv in self.lst[:]:
      if itv[0] > itv[1] : self.lst.remove(itv)
    self.lst.sort()
    dl = [] # delete list
    j = 0
    for i in range(1, len(self.lst)):
      if self.lst[i][0] < self.lst[j][1] or (self.adj_merge and self.lst[i][0] == self.lst[j][1]) :
        if self.lst[j][1] < self.lst[i][1] : self.lst[j] = (self.lst[j][0], self.lst[i][1]) # self.lst[j][1] = self.lst[i][1]
        dl.append(i)
      else : j = i
    for i in dl[::-1]: del self.lst[i]
  def __len__(self):
    return len(self.lst)
  def rlen(self):
    l = 0
    for itv in self.lst: l += itv[1] - itv[0]
    return l
  def __repr__(self):
    s = ['interval '+ self.id + ':']
    if len(self.lst) == 0 : s.append('empty!')
    #for itv in self.lst:
    else: s += ['[{}, {})'.format(itv[0],itv[1]) for itv in self.lst]
    return ' '.join(s)
  def __str__(self):
    return ','.join(['{}-{}'.format(itv[0],itv[1]) for itv in self.lst])
  def __getitem__(self, i): 
    return self.lst[i]
  def __add__(self, other): # do not change self
    new = Interval(itvs = self, adj_merge=self.adj_merge)
    new.lst += list(other)
    #for itv in other:
      #new.lst.append(itv[:])
    new.check()
    return new
  def add(self, other):
    return self + other
  def add_itv(self, itv, check=True): # input single interval, change the object
    self.lst.append(tuple(itv))
    if check: self.check()
    return self
  def sub_itv(self, itv):
    if itv[0] > itv[1] : return self
    lst = []
    for i in self.lst:
      if is_overlap(i, itv) : 
        if i[0] < itv[0] : lst.append((i[0],itv[0]))
        if i[1] > itv[1] : lst.append((itv[1],i[1]))
      else : lst.append(i)
    self.lst = lst
    self.check()
    return self
  def __sub__(self, other): # to be optimized, not tested
    '''substract other intervals, self and other should both be sorted
    '''
    lst = []
    j = 0
    for itv in self.lst:
      while j < len(other) and other[j][1] <= itv[0] : j += 1
      j1 = j
      while j1 < len(other) and other[j1][0] < itv[1] : j1 += 1
      if j1 == j : lst.append(itv) # no overlap interval
      else : # overlap with j...j1-1
        if itv[0] < other[j][0] : lst.append((itv[0], other[j][0]))
        for i in range(j, j1-1) : lst.append((other[i][1], other[i+1][0]))
        if itv[1] > other[j1-1][1] : lst.append((other[j1-1][1], itv[1]))
    new = Interval(itvs = lst, adj_merge=self.adj_merge)
    return new
    #new = Interval(itvs = self)
    #for itv in other:
      #new.sub_itv(itv[:])
    #new.check()
    #return new
  def sub(self, other):
    return self - other
  def ints_itv(self, itv): # intersect 
    if itv[0] > itv[1] : 
      self.lst = []
      return self
    lst = []
    for i in self.lst:
      if is_overlap(i, itv): lst.append((max(i[0], itv[0]), min(i[1], itv[1])))
    self.lst = lst
    self.check()
    return self
  def intersect(self, other): # optimized
    '''intersect intervals, self and other should both be sorted
    '''
    if len(other) == 0 or len(self) == 0 : return Interval()
    os = other[0][0]
    osi, c = self.nearest(os)
    lst = []
    j = 0
    #for itv in self.lst:
    for i in range(osi, len(self.lst)) :
      itv = self.lst[i]
      while j < len(other) and other[j][1] <= itv[0] : j += 1
      j1 = j
      while j1 < len(other) and other[j1][0] < itv[1] :
        lst.append((max(itv[0], other[j1][0]), min(itv[1], other[j1][1])))
        j1 += 1
      if j >= len(other) : break
    new = Interval(itvs = lst, adj_merge=self.adj_merge)
    return new
    #return (self + other) - (self - other) - (other - self) # hehe
  def nearest(self, p): 
    ''' nearest interval index in self.lst, return index, cmp(p, lst[index])
    '''
    l = len(self.lst)
    if l == 0 : return -1, 1
    lasti = i = l // 2
    i0, i1 = 0, l
    while i0 < i1 : 
      if self.lst[i][0] == p : break
      if self.lst[i][0] < p : i0 = i
      else : i1 = i
      lasti = i
      i = (i0 + i1) // 2
      if i == lasti : break
    if self.lst[i][0] <= p <= self.lst[i][1] : return i, 0 # inside
    if self.lst[i][0] > p : return i, -1 # p is smaller than ith interval, only when p < self.start
    if self.lst[i][1] < p : return i, 1 # p is greater than ith interval
    
  def is_inside(self, p, left = True, right = False, strand = '+') : # is p inside interval
    if len(self.lst) == 0 : return False
    i, c = self.nearest(p)
    if c != 0 : return False
    if self.lst[i][0] < p < self.lst[i][1] : return True
    if strand == '-' : left, right = right, left
    if left and self.lst[i][0] == p : return True
    if right and self.lst[i][1] == p : return True
    return False
  def nearestPos(self, p, strand = '+', downstream = True):
    i, c = self.nearest(p)
    if c == 0 : return p
    larger = True
    if strand != '-' and downstream == False : larger = False
    if strand == '-' and downstream == True :  larger = False
    if larger : 
      if c > 0 : plus = 1
      else : plus = 0
      if i + plus < len(self.lst) : return self.lst[i+plus][0]
      else : return None
    else : 
      if c > 0 : minus = 0
      else : minus = 1
      if i - minus >= 0 : return self.lst[i-minus][1]
      else : return None
  def nearestStop(self, p, strand = '+', downstream = True):
    i, c = self.nearest(p)
    if c != 0 : return p
    larger = True
    if strand != '-' and downstream == False : larger = False
    if strand == '-' and downstream == True :  larger = False
    if larger : return self.lst[i][1]
    else : return self.lst[i][0]
  def num_iter(self, start = None, step = 1, strand = '+'):
    if self.rlen() <= 0 : return
    if start is None : start = self.end5(strand)
    i = start
    if strand == '-' : step = -step
    s1, s2 = self.start, self.stop
    while s1 <= i <= s2 : 
      if self.is_inside(i, strand=strand) : yield i
      i += step
    #for i in range(start, self.stop, step):
      #if self.is_inside(i) : yield i
  @property
  def start(self):
    if self.is_empty() : return None
    else : return self.lst[0][0]
  @property
  def stop(self):
    if self.is_empty() : return None
    else : return self.lst[-1][1]
  def end5(self, strand = '+'):
    if strand != '-' : return self.start
    else : return self.stop
  def end3(self, strand = '+'):
    if strand != '-' : return self.stop
    else : return self.start
  def is_empty(self):
    return len(self.lst) == 0
  def is_compatible(self, other, mis = 0) :
    p1 = max(self.start, other.start)
    p2 = min(self.stop, other.stop)
    if p2 - p1 <= mis : return True # not overlap
    si = self.intersect([(p1, p2)])
    oi = other.intersect([(p1, p2)])
    
    m = si.sub(oi).rlen()
    if m > mis : return False
    m += oi.sub(si).rlen()
    if m > mis : return False
    return True
  def contradict(self, other) :
    p1 = max(self.start, other.start)
    p2 = min(self.stop, other.stop)
    if p2 <= p1 : return None # not overlap, empty
    si = self.intersect([(p1, p2)])
    oi = other.intersect([(p1, p2)])
    return si.sub(oi), oi.sub(si)
  def genome_pos(self, p, strand = '+', bias = 1): # exon only
    if len(self) == 0: return None
    m = self.rlen() # cdna_length()
    if p < 0 or p > m: return None
    if strand == '-': 
      p = m - p
      bias = 1 - bias
    if p == 0 : return self.start
    if p == m: return self.stop
    p1 = p
    for itv in self.lst:
      l = itv[1] - itv[0]
      if l - p1 >= bias:
        return p1 + itv[0] # e.start
      else:
        p1 -= l
    return None
  def cdna_pos(self, p, strand = '+', strict = False):
    '''if strict is True, the 3' end of exon will be considered as not in the transcript,
    if strict is False, 3' end of exon will be considered as start of the next exon, 
    or transcript end (self.cdna_length()) if in the last exon.
    '''
    if len(self.lst) == 0 : return None
    i, c = self.nearest(p)
    if c != 0 : return None
    if strict:
      if strand != '-' and p == self.lst[i][1]: return None
      if strand == '-' and p == self.lst[i][0]: return None
    pos = sum([self.lst[j][1] - self.lst[j][0] for j in range(i)])
    pos += p - self.lst[i][0]
    if strand == '-': 
      pos = self.rlen() - pos
    return pos
 
def trans2interval(t, start = 0, stop = None):
  '''generate intervals for a transcript
  start, stop are cDNA positions
  '''
  i = Interval(id = t.chr, itvs = [(e.start, e.stop) for e in t.exons])
  if start == 0 and stop is None : return i
  s1 = t.end5
  if start > 0 : s1 = t.genome_pos(start)
  s2 = t.end3
  if stop is not None : s2 = t.genome_pos(stop)
  if s2 is None : s2 = t.end3 # = t.genome_pos(start), t.genome_pos(stop)
  itv = [s1, s2]
  itv.sort()
  i.ints_itv(itv)
  return i
  #itv = Interval(id = t.chr)
  #for e in t.exons:
    #itv.lst.append((e.start, e.stop))
  #itv.check()
  #return itv

def cds_region_trans(t, cds1 = None, cds2 = None):
  '''generate CDS regions in 3 frames of chromosome
  '''
  cr = [Interval() for i in range(3)]
  if cds1 is None : cds1 = t.cds_start(cdna = True) 
  if cds2 is None : cds2 = t.cds_stop(cdna = True)
  if cds1 is None : return cr
  if cds2 is None : 
    tl = t.cdna_length()
    cds2 = cds1 + (tl - cds1) // 3 * 3
  #print cds1, cds2
  if cds2 - cds1 == 0 : return cr
  wrong = False
  if (cds2 - cds1) % 3 > 0 : 
    print('Wrong CDS annotation: %s %s %d %d %d' % (t.gid, t.id, cds1, cds2, t.cdna_length()))
    wrong = True
    return cr ## Wrong CDS annotation
  thick = [t.genome_pos(cds1, 1), t.genome_pos(cds2, 0)]
  thick.sort()
  ts1, ts2 = thick #t.thick_start, t.thick_stop
  orf = t(start = ts1, stop = ts2)
  exons = t.exons[:]
  exons.sort()
  if wrong : 
    for e in t.cds : 
      df = int(e.frame) # Use frame annotation
      if t.strand == '+' : df = - df
      f = (e.end5 - df) % 3
      cr[f].lst.append((e.start, e.stop))
  else : 
    df = 0
    for e in exons:
      for ei in e.intersect(orf):
        f = (ei.start - df) % 3
        cr[f].lst.append((ei.start, ei.stop))
        df = (len(ei) + df) % 3
      #print f, df, ei
    if df != 0 : 
      print('Remain df is not 0 ! %s %s %d %d' % (t.gid, t.id, cds1, cds2))
      return [Interval() for i in range(3)]
  for r in cr : r.check()
  return cr
def cds_region_gene(g):
  cr = [Interval() for i in range(3)]
  for t in g.trans:
    tcr = cds_region_trans(t)
    for i in range(3):
      cr[i] += tcr[i]
  return cr

def allTransRegions(trans):
  regions = {}
  for t in trans :
    transitv = trans2interval(t)
    if t.chr not in regions : regions[t.chr] = Interval(id = t.chr)
    regions[t.chr].lst += transitv.lst
  for chr in regions: regions[t.chr].check()
  return regions

def str2interval(s, blocksep = ',', rangesep = '-', id = ''):
  itvs = [list(map(int, b.split(rangesep))) for b in s.split(blocksep)]
  return Interval(itvs=itvs, id=id)


