'''
Fasta format file processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

from os.path import isfile

def faIter(file, protein = False):
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
      if not protein : l = l.replace('U','T').replace('u','t')
      sq += l
  yield (id, sq)

def rc(seq):
  '''reverse complement
  '''
  if seq is None: return None
  comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
       'B':"V", 'D':"H", 'H':"D", 'K':"M",
       'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
       'W':'W', 'N':'N', 'S':'S'}
  return ''.join([comps[x] for x in seq.upper()[::-1]])

def changechr(chr):
  '''change between two chr versions
  '''
  if chr[0].isdigit() or chr in ('X','Y','M'): return 'chr' + chr
  elif chr == 'MT' : return 'chrM'
  elif chr == 'chrM' : return 'MT'
  elif chr[0:3] == 'chr' : return chr[3:]
  else : return chr

chrmap = {}

class Faidx:
  '''index for genome fasta file
  '''
  def __init__(self, fid, length = 0, pos = 0, ls = 0, ll = 0):
    '''
    ls : sequence length in one line
    ll : line length
    '''
    self.id = fid
    self.length = length
    self.pos = pos
    self.ls = ls
    self.ll = ll
    self.is_endls = False
    self.is_endll = False
  def __str__(self):
    return "%s\t%d\t%d\t%d\t%d" % (self.id, self.length, self.pos, self.ls, self.ll)
  
  def updatels(self, ls):
    if ls <= 0 : raise Exception("Empty sequence!")
    if self.ls == 0 : self.ls = ls
    elif self.is_endls : raise Exception("Sequence has ended: %d" % (ls))
    elif ls > self.ls :
      raise Exception("Length error: %d -> %d" % (self.ls, ls))
    elif ls < self.ls : self.is_endls = True
  def updatell(self, ll):
    #print ll
    if ll <= 0 : raise Exception("Empty line!")
    if self.ll == 0 : self.ll = ll
    elif self.is_endll : raise Exception("Sequence has ended: %d" % (ll))
    elif ll > self.ll :
      raise Exception("Length error: %d -> %d" % (self.ll, ll))
    elif ll < self.ll : self.is_endll = True
      #print 'End: ', self.ll, ll

class Fa:
  '''indexed genome fasta file
  '''
  def __init__(self, fapath, verbose = False):
    self.file = open(fapath, 'r')
    idxpath = fapath + '.fai'
    if isfile(idxpath) : self.load_idx(idxpath)
    else : 
      if verbose : print('Making fasta index for {}...'.format(fapath))
      self.make_idx(idxpath)
  def make_idx(self, idxpath):
    self.file.seek(0)
    self.idxarr = []
    pos, flen = 0, 0
    for l in self.file:
      pos += len(l)
      if l[0] == '>' : 
        if flen > 0 : self.idxarr[-1].length = flen
        fid = l.strip().split()[0][1:]
        #print fid
        self.idxarr.append(Faidx(fid, pos = pos))
        flen = 0
      else : 
        ll, ls = len(l), len(l.strip())
        flen += ls
        self.idxarr[-1].updatels(ls)
        self.idxarr[-1].updatell(ll)
        #if ll > self.idxarr[-1][4] : self.idxarr[-1][4] = ll
    self.idxarr[-1].length = flen
    outfile = open(idxpath, 'w')
    self.idx = {}
    for i, idx in enumerate(self.idxarr):
      outfile.write('{}\n'.format(idx)) # print >>outfile, idx
      self.idx[idx.id] = idx
    outfile.close()
  def load_idx(self, idxpath):
    self.idx = {}
    for l in open(idxpath):
      lst = l.strip().split('\t')
      if len(lst) < 5 : continue
      idx = Faidx(lst[0])
      idx.length, idx.pos, idx.ls, idx.ll = list(map(int, lst[1:5]))
      self.idx[idx.id] = idx
  def get_chrname(self, chr) : #, chrmap = chrmap):
    #sprint(chrmap)
    if chr in self.idx : return chr
    mt = {'M':'MT', 'MT':'M'}
    if chr in mt and mt[chr] in self.idx: return mt[chr]
    if chr in chrmap :
      chr1 = chrmap[chr]
      if chr1 in self.idx : return chr1
      else : 
        chr2 = changechr(chr1)
        if chr2 in self.idx : return chr2
    chr1 = changechr(chr)
    if chr1 == chr : return None
    if chr1 in self.idx : return chr1
    elif chr1 in chrmap : 
      chr2 = chrmap[chr1]
      if chr2 in self.idx : return chr2
    return None
  def seekpos(self, fid, p):
    idx = self.idx[fid]
    pos = idx.pos # chr start
    pos += int(p / idx.ls) * idx.ll # line start
    pos += p % idx.ls # base start
    return pos

  def get_pos(self, fid, start = 0, stop = -1, length = -1):
  #def fetch(self, fid, start = 0, stop = -1, length = -1):
    chr = self.get_chrname(fid)
    if chr is None : 
      print("fa id {} not found in fasta file!".format(fid))
      return None, None, None, None
    fid = chr
    if start < 0 : start = start % self.idx[fid].length
    if length < 0 : 
      if stop < 0 : stop = stop % self.idx[fid].length
      length = stop - start
    return fid, start, stop, length

  def fetch(self, fid, start = 0, stop = -1, length = -1):
    fid, start, stop, length = self.get_pos(fid, start, stop, length)
    if fid is None: return None
    if length <= 0 : return ''
    self.file.seek(self.seekpos(fid, start))
    seq = ''
    idx = self.idx[fid]
    for l in self.file:
      l = l.strip()
      if length > len(l) :
        seq += l
        length -= len(l)
      else : 
        seq += l[0:length]
        break
    return seq
  def __iter__(self):
    for fid in self.idx:
      yield fid
    #return self.idx.itervalues()
  def exonSeq(self, e):
    s = self.fetch(e.chr, start = e.start, stop = e.stop)
    if e.strand == '-': s = rc(s)
    return s
  def transSeq(self, trans):
    s = []
    for e in trans.exons: s.append(self.exonSeq(e))
    return ''.join(s)

  def intervalSeq(self, fid, regions, strand='+'):
    chr = self.get_chrname(fid)
    if chr is None :
      print("fa id {} not found in fasta file!".format(fid))
      return None
    s = []
    for itv in regions:
      s.append(self.fetch(chr, itv[0], itv[1]))
    sq = ''.join(s)
    if strand == '-': sq = rc(sq)
    return sq


