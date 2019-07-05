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

def cmp3(a, b):
  ''' cmp for python3'''
  return (a > b) - (a < b)

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
    return cmp(len(self), len(other)) or cmp(self.start(),other.start())
  def cmp(self, other):
    return cmp3(len(self), len(other)) or cmp3(self.start(),other.start())
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0

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
    af = abs(f)
    o = Orf(frame = f)
    for i in range(af-1, length, codonSize):
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
  def cmp(self, other):
    return cmp3(self.start, other.start) or cmp3(self.stop, other.stop)
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0

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

codonFreqHuman = {'TTT':17.6, 'TTC':20.3, 'TTA':7.7,  'TTG':12.9, 'TCT':15.2, 'TCC':17.7, 'TCA':12.2, 'TCG':4.4,
                  'TAT':12.2, 'TAC':15.3, 'TAA':1.0,  'TAG':0.8,  'TGT':10.6, 'TGC':12.6, 'TGA':1.6,  'TGG':13.2,
                  'CTT':13.2, 'CTC':19.6, 'CTA':7.2,  'CTG':39.6, 'CCT':17.5, 'CCC':19.8, 'CCA':16.9, 'CCG':6.9,
                  'CAT':10.9, 'CAC':15.1, 'CAA':12.3, 'CAG':34.2, 'CGT':4.5,  'CGC':10.4, 'CGA':6.2,  'CGG':11.4,
                  'ATT':16.0, 'ATC':20.8, 'ATA':7.5,  'ATG':22.0, 'ACT':13.1, 'ACC':18.9, 'ACA':15.1, 'ACG':6.1,
                  'AAT':17.0, 'AAC':19.1, 'AAA':24.4, 'AAG':31.9, 'AGT':12.1, 'AGC':19.5, 'AGA':12.2, 'AGG':12.0,
                  'GTT':11.0, 'GTC':14.5, 'GTA':7.1,  'GTG':28.1, 'GCT':18.4, 'GCC':27.7, 'GCA':15.8, 'GCG':7.4,
                  'GAT':21.8, 'GAC':25.1, 'GAA':29.0, 'GAG':39.6, 'GGT':10.8, 'GGC':22.2, 'GGA':16.5, 'GGG':16.5,
              }

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
  shift = int(n / 4) # 4 bases
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

def senseMut2(seq, mutkey, cds1 = 0, cds2 = None) : # pass
  if cds2 is None : cds2 = len(seq)
  if type(mutkey) == str : keys = [mutkey]
  else : keys = mutkey
  seq = seq.upper()
  change = {}
  cands = {}
  has_key = True
  while has_key :
    has_key = False
    for key in keys :
      p = seq.find(key)
      while p >= 0 :
        has_key = True
        if p not in change : 
          change[p] = 0
          cands[p] = []
        newseq = sense_mut2(seq, p, p+len(key), cds1, cds2, cands[p], change[p])
        if newseq is None :
          print('Cannot find sense mutation for {} at {} of {}!'.format(key, p, seq))
          return None
        i, c = cands[p][change[p]][1:3]
        print('{} {} changed to {}, {}'.format(i, seq[i:i+len(c)], c, newseq))
        seq = newseq #, i, c = mut
        change[p] += 1
        p = seq.find(key)
  return seq

def mutated_seq(seq, cand):
  p, i, s = cand
  l = len(s)
  return seq[0:i] + s + seq[i+l:]

def sense_mut2(seq, start, stop, cds1, cds2, candsarr, n = 0) : # mutate any one base in a region. n is time of retries.
  #incds = True
  if n > 0: 
    if n < len(candsarr): return mutated_seq(seq, candsarr[n])
    else: return None
  if start >= cds2 or stop <= cds1 : # incds = False
    return base_mut2(seq, start, stop, candsarr, 0)
  else :
    if cds1 > start:
      base_mut2(seq, start, cds1, candsarr, 0) # upstream mutation
    if stop > cds2:
      base_mut2(seq, cds2, stop, candsarr, 0) # downstream mutation
    codon_mut2(seq, start, stop, cds1, cds2, candsarr, 0) # codon mutation
    return mutated_seq(seq, candsarr[0])

def base_mut2(seq, start, stop, candsarr, n = 0, priority = -100):
  if n >= (stop - start) * 4 : return None # exceed limit
  if n > 0: return mutated_seq(seq, candsarr[n])
  i0 = int((start + stop) / 2)
  for i in range(start, stop):
    dist = i - i0
    if dist <= 0:
      p = priority - dist*2
    else:
      p = priority + dist*2 + 1
    for j in range(4):
      candsarr.append((p+0.1*j, i, bases[j]))
  candsarr.sort()
  return mutated_seq(seq, candsarr[0])

def codon_mut2(seq, start, stop, cds1, cds2, candsarr, n = 0, priority = 0):
  if n > 0:
    if n < len(candsarr): return mutated_seq(seq, candsarr[n])
    else: return None
  for i in range(cds1, cds2, 3) :
    if i + 3 <= start : continue
    if i >= stop : break
    codon = seq[i:i+3]
    if codon not in codonTable: continue ## not a codon?
    aa = codonTable[codon]
    aac = AACodon[aa]
    for c in aac:
      p = priority - codonFreqHuman[c]
      candsarr.append((p, i, c))
  candsarr.sort()
  return mutated_seq(seq, candsarr[0])

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


