'''
Open reading frame (ORF) processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

from . import fa

codonSize = 3
cstart = ['ATG']
cstartlike = ['TTG', 'CTG', 'GTG', 'AAG', 'AGG', 'ACG', 'ATT', 'ATA', 'ATC'] #AAG & AGG may be removed
#caltstart = cstartlike
cstop = ['TGA', 'TAA', 'TAG']

senseframes = [1, 2, 3]
antiframes = [-1, -2, -3]
frames = [1, 2, 3, -1, -2, -3]

class Orf:
  '''multiple start with same stop (in the same ORF)
  '''
  def __init__(self, lst = [], frame = 0, stop = -1):
    if len(lst) != 0:
      #lst = l.strip().split('\t')
      self.frame = int(lst[0])
      if '' == lst[1] : 
        self.starts = []
      else : self.starts = list(map(int, lst[1].split(',')))
      if lst[2] == '' : self.altstarts = []
      else : self.altstarts = list(map(int, lst[2].split(',')))
      self.stop = int(lst[3])
    else:
      self.frame = frame
      self.starts = []
      self.altstarts = []
      self.stop = stop
      self.has_stop_codon = True # For Non-stop end
      #self.end = self.stop - 3 
  def __str__(self):
    return "%d\t%s\t%s\t%d" % (self.frame, ','.join(map(str, self.starts)), ','.join(map(str, self.altstarts)), self.stop)
  def __repr__(self):
    return "Orf object: " + str(self)
  def has_start(self):
    return len(self.starts) + len(self.altstarts) > 0
  def has_stop(self):
    return self.stop >= 0
  def __len__(self):
    if not self.is_complete() : return 0
    if self.stop > 0 : return self.stop - self.start()
    #elif self.end > 0 : return self.end - self.start()
    else : return 0
  def __cmp__(self, other):
    return cmp(len(self), len(other)) or cmp(self.start(),other.start)
  def length(self):
    return len(self)
  def aa_len(self):
    l = len(self)
    if l == 0 : return 0
    if self.has_stop_codon : l -= 3
    return l // codonSize
  #@property
  def start(self, alt = True):
    if not self.has_start() : return None
    if alt : return min(self.starts + self.altstarts)
    else :
      if len(self.starts) > 0 : return min(self.starts)
      else : return None
  def has_strictstart(self):
    return len(self.starts) > 0
  def is_complete(self):
    return self.has_start() and self.has_stop()
  def orf_by_ribo(self, frame, start, stop, srange = 4): #look for TIS near ribo suggested ORF
    if self.frame != frame: return None
    if self.stop < start: return None
    if min(self.starts + self.altstarts) > stop : return None
    dmin = srange * codonSize
    sm = -1
    for s in self.starts :
      d = abs(s - start)
      if d < dm :
        dm = d
        sm = s
    if sm > 0 : return (s, self.stop)
    for s in self.altstarts :
      d = abs(s - start)
      if d < dm :
        dm = d
        sm = s
    if sm > 0 : return (s, self.stop)
    return None
  def filtByLen(self, minaalen, tail = -1):
    end = self.stop - 3
    if self.stop < 0 : 
      end = tail
      start = self.start()
      if start is None : return
      end -= (end - start) % 3 # keep in frame
      self.stop = end
      self.has_stop_codon = False
    #if stop < 0 : return
    th = minaalen * 3
    rm = False
    for i, s in enumerate(self.starts):
      if end - s < th : 
        rm = True
        break
    if rm : self.starts[i:] = []
    rm = False
    for i, s in enumerate(self.altstarts):
      if end - s < th : 
        rm = True
        break
    if rm : self.altstarts[i:] = []

def allorf(seq, strand = '+', minaalen = 0, tail = -1) :
  '''find all possible ORFs in the sequence
  '''
  seq = seq.upper().replace('U','T')
  if strand == '+' : fr = senseframes
  elif strand == '-' : 
    fr = antiframes
    antiseq = fa.rc(seq)
  else: 
    fr = frames
    antiseq = fa.rc(seq)
  
  length = len(seq)
  for f in fr:
    if f > 0: s = seq
      #fa = f
    else: s = antiseq
    fa = abs(f)
    o = Orf(frame = f)
    for i in range(fa-1, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cstart: o.starts.append(i)
      elif codon in cstartlike: o.altstarts.append(i)
      if codon in cstop:
        #o.end = i
        o.stop = i + codonSize
        o.filtByLen(minaalen = minaalen, tail = tail)
        if o.has_start():
          yield o
        o = Orf(frame = f)
    o.filtByLen(minaalen = minaalen, tail = tail)
    if o.has_start() : yield o
def orflist(seq, strand = '+', sort = True, minaalen = 0, tail = -1):
  ol = list(allorf(seq, strand = strand, minaalen = minaalen, tail = tail))
  #for o in allorf(seq, strand = strand, minaalen = minaalen, tail = tail):
    #ol.append(o)
  if sort : ol.sort(reverse = True)
  return ol

class FixedOrf():
  '''only one start and one stop
  '''
  def __init__(self, start, stop = -1, has_stop_codon = True):
    self.start = start
    self.stop = stop
    self.has_stop_codon = has_stop_codon
  def __len__(self):
    if self.stop > 0 : return self.stop - self.start
    #elif self.end > 0 : return self.end - self.start
    else : return 0
  def length(self):
    return len(self)
  def __repr__(self):
    return "Fixed ORF object: %d - %d" % (self.start, self.stop)
  def frame(self): # 0 based
    return self.start % 3
  def __cmp__(self, other):
    return cmp(self.start, other.start) or cmp(self.stop, other.stop)
  def orfstr(self, seq):
    return "%s|%d-%d" % (seq[self.start:self.start+3], self.start, self.stop)
def orf_by_pos(seq, pos): ### Unkown start codon, only to find stop codon
  for i in range(pos, len(seq), codonSize):
    try: codon = seq[i:i+codonSize]
    except: FixedOrf(start = pos, stop = i, has_stop_codon = False)
    if codon in cstop: return FixedOrf(start = pos, stop = i + 3)
      #i += codonSize
      #break
  #o = FixedOrf(start = pos, stop = i)
  #return o
          
def orfDict(orflist, alt = True):
  '''dict of start -> stop
  '''
  od = {}
  for o in orflist:
    if o.stop < 0 : continue
    for s in o.starts: od[s] = o.stop
    if alt:
      for s in o.altstarts: od[s] = o.stop
  return od
def nearest_start(s, od, flank = 1):
  '''nearest start in ORF dict from s
  '''
  if s in od: return s
  for i in range(1, flank + 1):
    if s + i in od : return s + i ###
    if s - i in od : return s - i
  return None

codonTable = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L', 'TCT':'S','TCC':'S','TCA':'S','TCG':'S', 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*', 'TGT':'C','TGC':'C','TGA':'*','TGG':'W', 
               'CTT':'L','CTC':'L','CTA':'L','CTG':'L', 'CCT':'P','CCC':'P','CCA':'P','CCG':'P', 'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q', 'CGT':'R','CGC':'R','CGA':'R','CGG':'R', 
               'ATT':'I','ATC':'I','ATA':'I','ATG':'M', 'ACT':'T','ACC':'T','ACA':'T','ACG':'T', 'AAT':'N','AAC':'N','AAA':'K','AAG':'K', 'AGT':'S','AGC':'S','AGA':'R','AGG':'R', 
               'GTT':'V','GTC':'V','GTA':'V','GTG':'V', 'GCT':'A','GCC':'A','GCA':'A','GCG':'A', 'GAT':'D','GAC':'D','GAA':'E','GAG':'E', 'GGT':'G','GGC':'G','GGA':'G','GGG':'G', 
              }
AACodon = {'F':['TTC','TTT'], 'L':['CTG','CTC','CTT','TTG','TTA','CTA'], 'M':['ATG'], 'I':['ATC','ATT','ATA'], 'V':['GTG','GTC','GTT','GTA'],
           'S':['AGC','TCC','TCT','TCA','AGT','TCG'], 'P':['CCC','CCT','CCA','CCG'], 'T':['ACC','ACA','ACT','ACG'], 'A':['GCC','GCT','GCA','GCG'],
           'Y':['TAC','TAT'], '*':['TAA','TAG','TGA'], 'H':['CAC','CAT'], 'Q':['CAG','CAA'], 'N':['AAC','AAT'], 'K':['AAG','AAA'], 'D':['GAC','GAT'], 'E':['GAG','GAA'],
           'C':['TGC','TGT'], 'W':['TGG'], 'R':['AGA','AGG','CGG','CGC','CGA','CGT'], 'G':['GGC','GGA','GGG','GGT']
           } # ordered by human codon usage

allcodons = list(codonTable.keys())

def translate(seq):
  aa = ""
  for i in range(0, len(seq), codonSize):
    try : a = codonTable[seq[i:i+codonSize].upper().replace('U','T')]
    except : a = 'X'
    aa += a
  return aa

def senseMut(seq, mutkey, cds1 = 0, cds2 = None) : # pass
  if cds2 is None : cds2 = len(seq)
  if type(mutkey) == str : keys = [mutkey]
  else : keys = mutkey
  seq = seq.upper()
  change = {}
  has_key = True
  while has_key :
    has_key = False
    for key in keys : 
      p = seq.find(key)
      while p >= 0 : 
        has_key = True
        if p not in change : change[p] = 0
        mut = sense_mut1(seq, p, p+len(key), cds1, cds2, change[p])
        change[p] += 1
        if mut is None : 
          print('Cannot find sense mutation for {} at {} of {}!'.format(key, p, seq))
          return None
        newseq, i, c = mut
        print('{} {} changed to {}, {}'.format(i, seq[i:i+len(c)], c, newseq))
        seq = newseq #, i, c = mut
        p = seq.find(key)
  return seq

bases = ['T', 'C', 'A', 'G']
def sense_mut1(seq, start, stop, cds1, cds2, n = 0) : # mutate any one base in a region. n is time of retries.
  #incds = True
  if start >= cds2 or stop <= cds1 : # incds = False
    return base_mut(seq, start, stop, n)
  else : 
    n1 = (cds1 - start) * 4
    if n1 < 0 : n1 = 0
    if n < n1 : return base_mut(seq, start, cds1, n) # upstream mutation
    n2 = (stop - cds2) * 4
    if n2 < 0 : n2 = 0
    if n < n1 + n2 : return base_mut(seq, cds2, stop, n-n1) # downstream mutation
    return codon_mut(seq, start, stop, cds1, cds2, n-n1-n2) # codon mutation
########
def base_mut(seq, start, stop, n = 0) : 
  if n >= (stop - start) * 4 : return None # exceed limit
  i = int((start + stop) / 2)
  shift = n / 4 # 4 bases
  if shift % 2 > 0 : i -= int(shift / 2) + 1
  else : i += int(shift / 2)
  bi = n % 4
  newseq = seq[0:i] + bases[bi] + seq[i+1:]
  return newseq, i, bases[bi]

def codon_mut(seq, start, stop, cds1, cds2, n = 0) : # pass
  ni, ns = n, 0
  for i in range(cds1, cds2, 3) : 
    if i + 3 <= start : continue
    if i >= stop : return None
    codon = seq[i:i+3]
    aa = codonTable[codon]
    aac = AACodon[aa]
    if ni >= len(aac) : 
      ni -= len(aac)
      continue
    c2 = aac[ni]
    #for j in range(max(start-i,0), min(i+3-stop,3)) : 
    newseq = seq[0:i] + c2 + seq[i+3:]
    return newseq, i, c2

def is_start(seq, pos, alt = False, flank = 0):
  '''if start / alt start codon is nearby
  '''
  for i in range(pos - flank, pos + flank + 1):
    try : c = seq[i:i+codonSize]
    except : continue
    if c in cstart : return True
    elif alt and c in cstartlike : return True
  return False
def is_startlike(seq, pos, flank = 0):
  return is_start(seq, pos, alt = True, flank = flank)


