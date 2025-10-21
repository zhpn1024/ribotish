'''
Some useful tools
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

Nt=['A','C','G','T']

def rc(seq):
  '''reverse complement
  '''
  comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
       'B':"V", 'D':"H", 'H':"D", 'K':"M",
       'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
       'W':'W', 'N':'N', 'S':'S'}
  return ''.join([comps[x] for x in seq.upper()[::-1]])

def overlap(A, B):
  '''
  if A is overlapping with B.
  '''
  if chrcmp(A.chr, B.chr) != 0 : return False
  if (A.stop <= B.start) : return False
  if (B.stop <= A.start) : return False
  return True

def inside(A, B):
  '''
  if A is inside B.
  '''
  if(A.chr != B.chr) : return False
  if (B.stop < A.stop) : return False
  if (A.start < B.start) : return False
  return True
  
def distance(A,B):
  if A.chr!=B.chr: return None
  if overlap(A,B): return 0
  if A.start<B.start : return B.start-A.stop
  else: return A.start-B.stop
  
def centerin(A,B):
  '''
  If center of A is in B
  '''
  if overlap(A,B):
    return B.start<=A.center()<=B.stop
  else: return False

def end5(A):
  if A.strand=='-':
    return A.stop
  else:
    return A.start

def end3(A):
  if A.strand=='-':
    return A.start
  else:
    return A.stop

try: cmp
except:
  def cmp(a, b):
    return (a > b) - (a < b)

chr_order = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16,
             '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'M':25, 'MT':25,
             'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12,
             'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22,
             'chrX':23, 'chrY':24, 'chrM':25,
            }

def chrcmp(c1, c2):
  if c1 not in chr_order: return 1
  if c2 not in chr_order: return -1
  a, b = chr_order[c1], chr_order[c2]
  return (a > b) - (a < b)


def fa_iter(file): # also in fa.py
  id = sq = ""
  for l in file:
    l = l.strip()
    if l == '' : continue
    if l[0] == '>':
      if id != '':
        yield (id,sq)
      sq = ''
      id = l[1:]
    else:
      sq += l.replace('U','T').replace('u','t')
  yield (id,sq)

def cover_iter(bedIter, weight= lambda x: 1): 
  '''bed coverage
  '''
  current = []
  chr = ""
  start = 0
  clen = 0
  for bed in bedIter:
    bed.weight = weight(bed)
    #print bed
    if clen == 0 :
      chr = bed.chr
      start = bed.start
      current.append(bed)
      clen = 1
    elif chr == bed.chr:
      if start == bed.start: pass
      else:
        assert start < bed.start, "Records must be sorted!\n"
        while clen > 0 and current[0].stop <= bed.start:
          cover = 0
          for i in range(clen):
            cover += current[i].weight
          stop = current[0].stop
          yield (chr, start, stop, cover)
          #print clen
          start = stop
          i = 0
          while i < clen and current[i].stop <= stop:
            i += 1
          current[0:i] = []
          clen = len(current)
          #print "clean", i, clen
        if clen > 0 and current[0].stop > bed.start:
          cover = 0
          for i in range(clen):
            cover += current[i].weight
          stop = bed.start
          yield (chr, start, stop, cover)
          #print clen
          start = stop
      i = clen - 1
      while i >= 0 and bed.stop < current[i].stop:
        i -= 1
      current[i+1:i+1] = [bed]
      clen = len(current)
      if clen == 1: start = bed.start
    else:
      while clen > 0 :
        cover = 0
        for i in range(clen):
          cover += current[i].weight
        stop = current[0].stop
        yield (chr, start, stop, cover)
        start = stop
        i = 0
        while i < clen and current[i].stop <= stop:
          i += 1
        #for i in range(clen):
          #if current[i].stop > stop: break
        current[0:i] = []
        clen = len(current)
      chr = bed.chr
      start = bed.start
      current.append(bed)
      clen = 1
  
  while clen > 0 :
    cover = 0
    for i in range(clen):
      cover += current[i].weight
    stop = current[0].stop
    yield (chr, start, stop, cover)
    start = stop
    i = 0
    while i < clen and current[i].stop <= stop:
      i += 1
    #for i in range(clen):
      #if corrent[i].stop > stop: break
    current[0:i] = []
    clen = len(current)
    
def overlap_iter(bedIterA, bedIterB, func=overlap, ignoreStrand = True, chrcmp = chrcmp, counts = [0,0], alist = []): ### Need Revise!!
  '''overlap of two bed iteraters (files), both should be sorted.
  '''
  lst = alist # []
  ac = next(bedIterA) # .next() current bed A
  counts[0] += 1
  Aend = False
  for b in bedIterB:
    #print len(lst)
    counts[1] += 1
    j = -1
    cut = True
    for i in range(len(lst)):
      if chrcmp(lst[i].chr, b.chr) < 0: j = i
      elif chrcmp(lst[i].chr, b.chr) == 0:
        if cut and (lst[i].stop > b.start or func(lst[i], b)):
          j = i - 1
          cut = False          
        if func(lst[i], b):
          if ignoreStrand or lst[i].strand == b.strand : yield (lst[i], b)
    lst[0:(j+1)] = []
    if chrcmp(ac.chr, b.chr) > 0 or (chrcmp(ac.chr, b.chr) == 0 and ac.start > b.stop) : continue
    if func(ac, b):
      if ignoreStrand or ac.strand == b.strand : yield (ac, b)
    if Aend == False : lst.append(ac)
    c = 0
    for a in bedIterA:
      counts[0] += 1
      assert chrcmp(a.chr, ac.chr) > 0 or a.start >= ac.start, "Records must be sorted!\n"+str(a)
      ac = a
      c += 1
      if chrcmp(a.chr, b.chr) < 0: continue
      elif chrcmp(a.chr, b.chr) == 0:
        if a.stop < b.start: continue
        if func(a, b):
          if ignoreStrand or a.strand == b.strand : yield (a, b)
        elif a.start > b.stop: break
      else: break
      lst.append(a)
    if c == 0 : Aend = True
    if Aend and len(lst) == 0: break
    #print len(lst)

def rand_overlap_iter(bedIterA, bedListB, func=overlap, ignoreStrand = True): # ListB should have no overlap
  '''overlap of two bed iteraters (files), A is random, B is list, should be sorted and have no overlap.
  '''
  m = len(bedListB)
  for a in bedIterA:
    i0 = 0
    i1 = m-1
    #print a
    while i1 - i0 > 1 :
      i = (i1 + i0) // 2
      #print i, i0, i1
      if a > bedListB[i] : i0 = i
      else: i1 = i
    i = i0
    while func(a, bedListB[i]):
      if ignoreStrand or a.strand == bedListB[i].strand :
        yield (a, bedListB[i])
      i -= 1
      if i < 0: break
    i = i1
    while func(a, bedListB[i]):
      if ignoreStrand or a.strand == bedListB[i].strand :
        yield (a, bedListB[i])
      i += 1
      if i >= m: break
      #else : 
        #print 'i0', ignoreStrand, a.id, a.strand, bedListB[i0].id, bedListB[i0].strand
    #if func(a, bedListB[i1]):
      #if ignoreStrand or a.strand == bedListB[i1].strand :
        #yield (a, bedListB[i1])
      #else : 
        #print 'i1', ignoreStrand, a,id, a.strand, bedListB[i1].id, bedListB[i1].strand

def find_overlap(q, lst, func=overlap, ignoreStrand = True): 
  '''List should have no overlap
  '''
  m = len(lst)
  i0 = 0
  i1 = m-1
  while i1 - i0 > 1 :
    i = (i1 + i0) / 2
    if q > lst[i] : i0 = i
    else: i1 = i
  if func(q, lst[i0]):
    if ignoreStrand or q.strand == lst[i0].strand :
      return i0
  if func(q, lst[i1]):
    if ignoreStrand or q.strand == lst[i0].strand :
      return i1
  return -1

def bed2seq(seq, bed):
  if type(seq) == str : idx = False
  else : idx = True
  s = ''
  for e in bed.exons:
    if idx : es = seq.fetch(e.chr, start = e.start, stop = e.stop)
    else : es = seq[e.start:e.stop]
    if e.strand == '-':
      es = rc(es)
    s += es
  return s

def gtf2seq(seq, gtftrans):
  if type(seq) == str : idx = False
  else : idx = True
  s = ''
  for e in gtftrans.exons:
    if idx : es = seq.fetch(e.chr, start = e.start, stop = e.stop)
    else : es = seq[e.start:e.stop]
    if e.strand == '-':
      es = rc(es)
    s += es
  return s

def trans2seq(seq, trans):
  if type(seq) == str : idx = False
  else : idx = True
  s = ''
  for e in trans.exons:
    if idx : es = seq.fetch(e.chr, start = e.start, stop = e.stop)
    else : es = seq[e.start:e.stop]
    if e.strand == '-':
      es = rc(es)
    s += es
  return s


def exon2seq(seq, e):
  if type(seq) == str : s = seq[e.start:e.stop]
  else : s = seq.fetch(e.chr, start = e.start, stop = e.stop)
  if e.strand == '-':
    s = rc(s)
  return s

def cdna_pos(trans, p):
  exons = trans.exons
  exons.sort(reverse = trans.is_reverse())
  if p < trans.start or p > trans.stop:
    return None
  #p1 = p - trans.start
  pos = 0
  for e in exons:
    if e.start <= p <= e.stop:
      if trans.is_reverse():
        pos += e.stop - p
      else:
        pos += p - e.start
      return pos
    else:
      pos += len(e)
  return None

def genome_pos(trans, p, bias = 1):
  if p < 0 or p > trans.cdna_length():
    return None
  p1 = p
  pos = trans.start
  #l=range(len())
  for e in trans.exons:
    if len(e) - p1 >= bias:
      if trans.is_reverse():
        pos = e.stop - p1
      else:
        pos = p1 + e.start
      return pos
    else:
      p1 -= len(e)
  return None

def non_overlap_lists(lst):
  nolists=[]
  for b in lst:
    f = False
    for l in nolists:
      if overlap(b, l[-1]):
        l.append(b)
        f = True
        break
    if f : continue
    nolists.append([b])
  return nolists

def range_to_bins(start, end):
  binOffsets = [512+64+8+1, 64+8+1, 8+1, 1, 0]
  binFirstShift = 17    # How much to shift to get to finest bin.
  binNextShift = 3     # How much to shift to get to next larger bin.
  startBin = start
  endBin = end - 1
  startBin >>= binFirstShift
  endBin >>= binFirstShift
  for bo in binOffsets:
    if startBin == endBin: return bo + startBin
    startBin >>= binNextShift
    endBin >>= binNextShift

def bins_overlap_range(start, end):
  binOffsets = [512+64+8+1, 64+8+1, 8+1, 1, 0]
  binFirstShift = 17    # How much to shift to get to finest bin. 
  binNextShift = 3     # How much to shift to get to next larger bin. 
  bins = []
  startBin = start
  endBin = end - 1
  startBin >>= binFirstShift
  endBin >>= binFirstShift
  for bo in binOffsets:
    for j in range(startBin, endBin+1):
      bins.append(bo + j);
    startBin >>= binNextShift
    endBin >>= binNextShift
  return bins;


def numround(f, r1 = 4, r2 = 2, exact = False):
  f = float(f)
  if f == 0: return '0.0'
  if abs(f) < 0.001:
    s = ('%.{}e'.format(r2) % (f)) .replace('e-0', 'e-')
    if not exact:
      s = s.replace('0e', 'e')
      s = s.replace('.e', '.0e')
  else:
    s = '%.{}f'.format(r1) % (f)
    if not exact:
      s = s.rstrip('0')
      if s.endswith('.'): s += '0'
  return s

def downsample_even(arr, n):
  l = len(arr)
  if l <= n: return arr
  out = []
  m = l * n * 2
  #print(m)
  d = l * 2 # float(l) / n
  s = l # d / 2
  for i in range(n):
    j = s // n // 2 # int(s)
    #print(s, j)
    out.append(arr[j])
    s += d
  return out
