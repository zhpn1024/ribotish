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
  if(A.chr != B.chr) : return False
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
    
def overlap_iter(bedIterA, bedIterB, func=overlap, ignoreStrand = True, chrcmp = cmp, counts = [0,0]): ### Need Revise!!
  '''overlap of two bed iteraters (files), both should be sorted.
  '''
  lst = []
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





