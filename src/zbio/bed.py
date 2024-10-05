'''
Bed format annotation processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''
def cmp3(a, b):
  ''' cmp for python3'''
  return (a > b) - (a < b)

class Bed3:
  ''' Bed3 class for bed3 annotation format, father of Bed6
  '''
  Header=('chr','start','stop')
  Format=(str,int,int)
  n=3
  def __init__(self, x, bin = False): #Fit all bed
    '''init for all bed format
    x: can be a bed format string, Bed object, dict or list
    '''
    if type(x)==str:
      l = x.strip().split('\t')
      if bin: l[0:1] = []
    elif type(x) == list or type(x) == tuple: l = x
    #elif type(x) == type(self) : l = [getattr(x, i) for i in self.Header]
    elif type(x)==dict: l=[x[i] for i in self.Header]
    else: l = [getattr(x, i) for i in self.Header]
    lst=[]
    for i in range(self.n):
      try: lst.append(self.Format[i](l[i]))
      except: 
        if i in (6, 7) : lst.append(lst[1]) # blockstart/stop
        elif i in (8,) : lst.append('0') # rgb
        elif i == 9 : lst.append(1)
        elif i == 10 : lst.append((lst[2]-lst[1],))
        elif i == 11 : lst.append((0,))
        elif i in (3,4,5) : lst.append('.')
        else : raise ValueError
    self.items = tuple(lst)
    self.symbol = self.gid = self.tid = self.id ###
    self.genetype = ''
    if type(x) == dict:
      for i in x :
        if i not in self.Header : setattr(self, i, x[i])
  def __str__(self): #Bed3 and Bed6
    return '\t'.join(map(str,self.items))
  def __repr__(self): #All bed
    return 'Bed'+str(self.n)+' object:\n'+str(self)+'\n'
  def short_str(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  @property
  def chr(self): #All bed
    return self.items[0]
  @property
  def start(self): #All bed
    return self.items[1]
  @property
  def stop(self): #All bed
    return self.items[2]
  @property
  def id(self): #Bed3 only
    return "noname"
  #@property
  #def tid(self):
    #return self.id
  #@property
  #def gid(self):
    #return self.tid
  #@property
  #def symbol(self):
    #return self.gid
  #@property
  #def genetype(self):
    #return ""
  @property
  def score(self): #Bed3 only
    return 0.0
  @property
  def strand(self): #Bed3 only
    return "."
  def __len__(self): #All bed
    return self.stop - self.start
  @property
  def length(self):
    return len(self)
  @property
  def end(self): #All bed
    return self.stop
  @property
  def dict(self): #All bed
    return dict([(self.Header[i],self.items[i]) for i in range(self.n)])
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start

  def is_reverse(self): #All bed
    return self.strand=='-'
  def bed(self): #Bed copy, Bed3 only
    '''make a copy of self
    '''
    return Bed3(self)
  def __getitem__(self,i): #All bed
    return self.items[i]
  def __cmp__(self, other): #All bed
    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
  def cmp(self, other):
    return cmp3(self.chr,other.chr) or cmp3(self.start,other.start) or cmp3(self.stop,other.stop)
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0
  def __sub__(self, other):
    return sub(self, other)
  def intersect(self, other):
    return intersect(self, other)
  def union(self, other):
    return union(self, other)
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return sep.join(self.Header)
  def head(self,n=10): #Bed3 only
    if n>len(self):
      n=len(self)
    return Bed3([self.chr,self.start,self.start+n])
  def tail(self,n=10): #Bed3 only
    if n>len(self):
      n=len(self)
    return Bed3([self.chr,self.stop-n,self.stop])
  def __call__(self,**args): 
    '''New bed that modified from self. For all bed
    '''
    cp=self.bed()
    if args=={}: return cp
    d=self.dict
    for k in args: d[k] = args[k]
      #else: assert 0, k+' not defined!\n'
    cp.__init__(d)
    return cp
  def center(self): #Middle point, Bed3 and Bed6
    return (self.start+self.stop)/2.0
  def cdna_length(self): #Bed3 and Bed6
    return len(self)
  @property
  def exons(self): #Bed3 and Bed6
    l=[]
    e = self(id=self.id+"_Exon_1")
    e.symbol, e.gid, e.tid = self.symbol, self.gid, self.tid
    l.append(e)
    return l
  @property
  def introns(self): # No intron, empty
    return []
  @property
  def trans(self):
    return [self]
  @property
  def thick_start(self):
    return self.start
  @property
  def thick_stop(self):
    return self.start
  @property
  def itemRgb(self):
    return 0
  @property
  def blockCount(self):
    return 1
  @property
  def blockSizes(self):
    return (len(self),)
  @property
  def blockStarts(self):
    return (0,)

  def merge_trans(self):
    return self
  def is_contain(self, p, strict = False): # all bed
    if self.start <= p <= self.stop :
      if not strict : return True
      else : return p != self.end3
    else : return False
  def is_upstream(self, p, strict = False) : # all bed
    if self.strand != '-' : return p < self.start
    else : return p > self.stop
  def is_downstream(self, p, strict = False) : # all bed
    if self.strand != '-' : 
      if strict : return p >= self.stop
      else : return p > self.stop
    else : 
      if strict : return p <= self.start
      else : return p < self.start
  def is_exon(self, p, strict = False): #Bed3 and Bed6
    return self.is_contain(p)
  def is_intron(self, p, strict = False): #False, Bed3 and Bed6
    return False
  def is_sense(self,other): #In same strand? All bed
    return self.strand==other.strand

  def flank_pos(self, p, flank = 100):
    if flank < 0: return None
    if self.start - flank <= p <= self.start:
      if self.strand != '-': return p - self.start # upstream
      else: return self.start - p + self.cdna_length()
    if self.stop <= p <= self.stop + flank:
      if self.strand != '-': return p - self.stop + self.cdna_length()
      else: return self.stop - p # upstream
    return None

  def cdna_pos(self, p, strict = False, flank = 0): #Bed3 and Bed6
    if self.is_contain(p, strict): return abs(p-self.end5)
    elif flank > 0: return self.flank_pos(p, flank)
    else: return None
  def genome_pos(self, p, bias = 1): #Bed3 and Bed6
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    if not self.is_reverse() : return self.start + p
    else : return self.stop - p
  def is_overlap(self, other):
    if(self.chr != other.chr) : return False
    if (self.stop <= other.start) : return False
    if (other.stop <= self.start) : return False
    return True

  def is_compatible(self, other, ignoreStrand = False): # test needed
    if not ignoreStrand :
      if self.strand != other.strand : return False
    if not self.is_overlap(other) : return True
    start = max(self.start, other.start)
    stop = min(self.stop, other.stop)
    int1 = []
    for i in self.introns:
      if i.start <= start : return False
      if i.stop >= stop : return False
      int1.append(i)
    int2 = []
    for i in other.introns:
      if i.stop <= start : continue
      if i.start <= start : return False
      if i.stop >= stop : return False
      int2.append(i)
    if len(int1) != len(int2) : return False
    int1.sort()
    int2.sort()
    for i in range(len(int1)):
      if int1[i] != int2[i] : return False
    return True

class Bed6(Bed3):
  Header=('chr','start','stop','id','score','strand')
  Format=(str,int,int,str,str,str)
  n=6
  #def __init__(self,x):
  #def __str__(self):
  @property
  def id(self): #Bed6 and Bed12
    return self.items[3]
  @property
  def score(self): #Bed6 and Bed12
    return self.items[4]
  @property
  def strand(self): #Bed6 and Bed12
    return self.items[5]
  @property
  def anti_strand(self):
    if self.strand == '+': return '-'
    elif self.strand == '-': return '+'
    else : return '.'
  def bed(self): #Only Bed6, make a copy
    return Bed6(self)
  
  def head(self, n=10): #Bed6 and Bed12
    if n>len(self): n=len(self)
    if self.strand != '-':
      start=self.start
      stop=start+n
    else:
      stop=self.stop
      start=stop-n
      if start<0: start=0
    return Bed6([self.chr,start,stop,self.id+'_Head_'+str(n),self.score,self.strand])
  
  def tail(self,n=10): #Bed6 and Bed12
    if n>len(self): n=len(self)
    if self.strand == '-':
      start=self.start
      stop=start+n
    else:
      stop=self.stop
      start=stop-n
      if start<0: start=0
    return Bed6([self.chr,start,stop,self.id+'_Tail_'+str(n),self.score,self.strand])
  
  #def is_upstream(self, other): #other is upstream of self
    #if self.strand == '-' : return other.start >= self.stop
    #else : return other.stop <= self.start

def com2tup(s): #comma string to tuple
  if type(s)==tuple: return s
  if type(s)==list: return tuple(s)
  return tuple(map(int,s.strip().strip(',').split(',')))
def tup2com(t): #tuple to comma string
  if type(t)==str: return t
  return ','.join(map(str,t))+','

class Bed12(Bed6):
  
  Header=('chr','start','stop','id','score','strand',"thick_start","thick_stop","itemRgb","blockCount","blockSizes","blockStarts")
  Format=(str,int,int,str,str,str,int,int,com2tup,int,com2tup,com2tup)
  n=12
  @property
  def thick_start(self): #Bed12
    return self.items[6]
  @property
  def thick_stop(self):
    return self.items[7]
  #@property
  def cds_start(self, cdna = False): #Bed12
    if self.items[6] == self.items[7] : return None ##
    if self.strand != '-' : s = self.items[6]
    else : s = self.items[7]
    if not cdna : return s
    else : return self.cdna_pos(s)
  #@property
  def cds_stop(self, cdna = False):
    if self.items[6] == self.items[7] : return None ##
    if self.strand != '-' : s = self.items[7]
    else : s = self.items[6]
    if not cdna : return s
    else : return self.cdna_pos(s)
  @property
  def itemRgb(self):
    return self.items[8]
  @property
  def blockCount(self):
    return self.items[9]
  @property
  def blockSizes(self):
    return self.items[10]
  @property
  def blockStarts(self):
    return self.items[11]

  @property
  def blockStops(self):
    return tuple(map(lambda x,y:x+y,self.blockStarts,self.blockSizes))
  def blockStop(self,i): 
    return self.blockStarts[i]+self.blockSizes[i]
  
  def __str__(self): #Bed string, Bed12 only
    l=list(self.items)
    l[8]=tup2com(l[8])
    l[10]=tup2com(l[10])
    l[11]=tup2com(l[11])
    return '\t'.join(map(str,l))

  def genePredStr(self, geneName = True, extended = False):
    lst = []
    if geneName : lst += [self.gid]
    lst += [self.id, self.chr, self.strand, self.start, self.stop, self.thick_start, self.thick_stop, self.blockCount]
    exonStarts = [x + self.start for x in self.blockStarts]
    exonStops = [x + self.start for x in self.blockStops]
    lst += [tup2com(exonStarts), tup2com(exonStops)]
    if extened : 
      lst += [self.name2, tup2com(self.exonFrames), self.cdsStartStat, self.cdsEndStat]
    return '\t'.join(map(str, lst))
  def bed(self): #Bed copy, Bed12 only
    return Bed12(self)
  
  def cdna_length(self): #Bed12
    #l = 0
    #for i in self.blockSizes: l += i
    return sum(self.blockSizes)
  
  def center(self): #Middle point, Bed12 only
    len=self.cdna_length()/2.0
    for i in range(self.blockCount):
      len-=self.blockSizes[i]
      if len < 0 : break
    return self.start+self.blockStarts[i]+self.blockSizes[i]+len
  @property
  def exons(self): #Bed12
    if hasattr(self, '_exons') : return self._exons
    a=[]
    if self.strand=="-": step, j = -1, self.blockCount
    else: step, j = 1, 1
    for i in range(self.blockCount):
      start = self.start + self.blockStarts[i]
      end = start + self.blockSizes[i]
      id=self.id+"_Exon_"+str(j)
      j+=step
      e = Bed6((self.chr,start,end,id,self.score,self.strand))
      e.symbol, e.gid, e.tid = self.symbol, self.gid, self.tid ##
      a.append(e)
    if self.strand=="-": self._exons = a[::-1]
    else: self._exons = a
    return self._exons
  @property
  def introns(self):
    if hasattr(self, '_introns') : return self._introns
    a=[]
    if self.strand=="-": step, j = -1, self.blockCount - 1
    else : step, j = 1, 1
    for i in range(self.blockCount-1):
      start=self.start+self.blockStarts[i]+self.blockSizes[i]
      end=self.start+self.blockStarts[i+1]
      id=self.id+"_Intron_"+str(j)
      j+=step
      e = Bed6((self.chr,start,end,id,self.score,self.strand))
      e.symbol, e.gid, e.tid = self.symbol, self.gid, self.tid ##
      a.append(e)
    if self.strand=="-": self._introns = a[::-1]
    else : self._introns = a
    return self._introns

  def abs2relative(self):
    '''change absolute blockStart values to relative values
    '''
    blockstarts = [st - self.start for st in self.blockStarts]
    return self(blockStarts = blockstarts)

  def cdna_pos(self, p, strict = False, flank = 0):
    '''if strict is True, the 3' end of exon will be considered as not in the transcript,
    if strict is False, 3' end of exon will be considered as start of the next exon, 
    or transcript end (self.cdna_length()) if in the last exon.
    '''
    if not self.is_contain(p, strict) : # return None
      if flank <= 0: return None
      else: return self.flank_pos(p, flank)
    p1 = p - self.start # genome relative posistion
    pos = 0 # cDNA pos
    l=list(range(self.blockCount))
    if self.is_reverse(): l = l[::-1]
    for i in l:
      if self.blockStarts[i]<=p1<=self.blockStop(i):
        if strict :
          if self.strand != '-' and p1 == self.blockStop(i): return None
          if self.strand == '-' and p1 == self.blockStarts[i]: return None
        if self.is_reverse(): pos+=self.blockStop(i)-p1
        else: pos+=p1-self.blockStarts[i]
        return pos
      else: pos+=self.blockSizes[i]
      #print(pos)
    else: return None

  def genome_pos(self, p, bias=1):
    '''if bias is 1, the splice junction will be mapped to 5' end of downstream exon,
    if bias is 0, the splice junction will be mapped to 3' end of the upstream exon.
    '''
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    p1=p
    pos=self.start
    l=list(range(self.blockCount))
    if self.is_reverse(): l=l[::-1]
    for i in l:
      if self.blockSizes[i]-p1>=bias:
        if self.is_reverse(): pos+=self.blockStop(i)-p1
        else: pos+=p1+self.blockStarts[i]
        return pos
      else: p1-=self.blockSizes[i]
    else: return None 
  
  #def cdna_type(p):
    
  def type(self, p, genome = True):
    ''' 5UTR, CDS or 3UTR
    '''
    if genome == False:
      pg = self.genome_pos(p)
      return self.type(pg)
    if type(p) != int : return None
    if p > self.stop or p < self.start : return None
    if p < self.thick_start:
      if self.is_reverse(): return "3UTR"
      else: return "5UTR"
    elif  p > self.thick_stop:
      if self.is_reverse(): return "5UTR"
      else: return "3UTR"
    else: return "CDS"

  def is_exon(self, p, strict = False): #If in exon, Bed12
    if not self.is_contain(p, strict = strict): return False
    l=list(range(self.blockCount))
    if self.is_reverse() : l=l[::-1]
    p1=p-self.start
    for i in l:
      if self.blockStarts[i]<=p1<=self.blockStop(i):
        if strict :
          if self.strand != '-' and p1 == self.blockStop(i): return False
          if self.strand == '-' and p1 == self.blockStarts[i]: return False
        return True
    else: return False
    #return self.is_contain(p)
  def is_intron(self, p, strict = False): #If in intron
    if not self.is_contain(p, strict = strict): return False
    l=list(range(1,self.blockCount))
    if self.is_reverse() : l=l[::-1]
    p1=p-self.start
    for i in l:
      if self.blockStarts[i] >= p1 >= self.blockStop(i-1):
        if strict :
          if self.strand != '-' and p1 == self.blockStarts[i]: return False
          if self.strand == '-' and p1 == self.blockStop(i-1): return False
        return True
    else: return False

def bed3_iter(file, bin = False):
  for l in file:
    try: yield Bed3(l, bin)
    except ValueError: pass


def bed6_iter(file, bin = False):
  for l in file:
    try: yield Bed6(l, bin)
    except ValueError: pass

def bed12_iter(file, bin = False, filt = [], chrs = None, verbose = False):
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ''
  for l in file:
    try: b = Bed12(l, bin)
    except ValueError: continue
    if b.chr != chr :
      if chrs is not None and b.chr not in chrs : continue
      if verbose : print(b.chr)
      chr = b.chr
    if len(filt) > 0 and b.id not in dfilt: continue
    yield b
    
def bed12SelectIter(file, select = None, **kwargs):
  for b in bed12_iter(file, **kwargs):
    if select is None : yield b
    else : 
      if select(b) is not None : yield b
        
def bed12_fetch(file, id, bin = False, verbose = False):
  for b in bed12_iter(file, bin = bin, verbose = verbose) :
    if b.id == id : return b


class refGene(Bed12): # start with bins
  def __init__(self, x):
    if type(x) == list: l = x
    else: l = x.strip().split('\t')
    self.exonStarts = com2tup(l[9])
    self.exonEnds = com2tup(l[10])
    start = int(l[4])
    blockSizes = [self.exonEnds[i] - self.exonStarts[i] for i in range(len(self.exonStarts))]
    blockStarts = [self.exonStarts[i] - start for i in range(len(self.exonStarts))]
    l2 = [l[2],l[4],l[5],l[1],l[11],l[3],l[6],l[7],'0',l[8],tup2com(blockSizes),tup2com(blockStarts)]
    Bed12.__init__(self, l2)
    self.bin = int(l[0])
    self.gid = self.name2 = l[12]
    self.exonFrames = com2tup(l[15])
    self.cdsStartStat = l[13]
    self.cdsEndStat = l[14]
    
  def __repr__(self): 
    return 'refGene object:\n'+str(self)+'\n'
    
def refGene_iter(file):
  for l in file:
    try: yield refGene(l)
    except ValueError: pass
    
def gpd_iter(file, filt = [], chrs = None, verbose = False): #refFlat_iter(file):
  '''GenePred with Gene Names
  '''
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ''
  for l in file:
    lst = l.strip().split()
    starts = com2tup(lst[9])
    stops = com2tup(lst[10])
    st1 = tuple(map(lambda x: x - int(lst[4]), starts))
    sizes = tuple(map(lambda x, y: y - x, starts, stops))
    lstb = [lst[2],lst[4],lst[5],lst[1],lst[0],lst[3],lst[6],lst[7],"0,0,0",lst[8],sizes,st1]
    bed = Bed12(lstb)
    bed.symbol = bed.gid = bed.score ##
    if bed.chr != chr : 
      if chrs is not None and bed.chr not in chrs : continue
      if verbose : print(bed.chr)
      chr = bed.chr
    if len(filt) > 0 and bed.id not in dfilt: continue
    yield bed
refFlat_iter = gpd_iter

class gpdGene:
  '''genepred gene
  '''
  def __init__(self, id = '', b = None):
    if b is None :
      self.trans = []
      self.chr, self.start, self.stop, self.strand = '', None, None, ''
      self.id = id
    else : 
      self.trans = [b]
      self.chr, self.start, self.stop, self.strand = b.chr, b.start, b.stop, b.strand
      self.id = b.gid
  def add_trans(self, tr):
    if self.strand in ('+', '-') and tr.strand != self.strand :
      print('Inconsistent trans strand: {} {} {} {}'.format(tr.gid, tr.tid, tr.strand, self.strand))
      self.strand = '.'
    elif self.strand == '': self.strand = tr.strand
    self.trans.append(tr)
  def __repr__(self):
    s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + Exon.__repr__(self)
    for t in self.trans:
      s += '\n\t' + t.__repr__()
    return s
  def check(self):
    if len(self.trans) > 0 :
      if self.chr == '' : self.chr = self.trans[0].chr
      if self.trand == '' : self.strand = self.trans[0].strand
      if self.start is None : self.start, self.stop = self.trans[0].start, self.trans[0].stop
    for t in self.trans:
      #t.check()
      if t.start < self.start: self.start = t.start
      if t.stop > self.stop: self.stop = t.stop
  def merge_trans(self):
    '''generate a new transcript that merge all transcript exons in the gene
    '''
    es = []
    merge = []
    for t in self.trans:
      es += t.exons
    es.sort()
    me = es[0]
    for i in range(1, len(es)):
      if me.stop >= es[i].start:
        me = union(me, es[i])[0]
      else: 
        merge.append(me)
        me = es[i]
    merge.append(me)
    #merge.check()
    return exons2bed12(merge)
def selectMaxCDS(trans): 
  '''select transcript with max CDS length
  '''
  maxlen, mt = 0, None
  for t in trans:
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
    try : cdslen = cds2 - cds1 
    except : continue
    if cdslen % 3 != 0 : continue
    if cdslen > maxlen : maxlen, mt = cdslen, t
  return mt

def gpdGeneIter(file, **kwargs):
  chr, gt = '', {}
  for b in gpd_iter(file, **kwargs):
    if b.chr != chr :
      for gid in gt : 
        yield gt[gid]
        #t = selectMaxCDS(gt[gid])
        #if t is not None : yield t
      chr, gt = '', {}
    if b.gid not in gt : gt[b.gid] = gpdGene(b = b)
    else : gt[b.gid].add_trans(b)
  for gid in gt : yield gt[gid]
    #t = selectMaxCDS(gt[gid])
    #if t is not None : yield t


def gpdSelectIter(file, select = selectMaxCDS, **kwargs):
  chr, gt = '', {}
  for b in gpd_iter(file, **kwargs):
    if b.chr != chr :
      for gid in gt : 
        t = selectMaxCDS(gt[gid])
        if t is not None : yield t
      chr, gt = '', {}
    if b.gid not in gt : gt[b.gid] = []
    gt[b.gid].append(b)
  for gid in gt : 
    t = selectMaxCDS(gt[gid])
    if t is not None : yield t
      
def gpd_fetch(file, id):
  for b in gpd_iter(file):
    if b.id == id : return b

def short_bed(s, name = ''): #like chr1:1-200:+
  lst = s.strip().split(":")
  l1 = lst[1].split("-")
  if name == "": name = s.strip()
  if len(lst) >= 3: strand = lst[2]
  else: stand = '.'
  return Bed6([lst[0], l1[0], l1[1], name, 0, strand])

def exons2bed12(exons, id = 'exons_merge', thick_start = None, thick_stop = None) : #From exons to bed12
  exons.sort()
  chr, start, stop = exons[0].chr, exons[0].start, exons[-1].stop
  strand, score = exons[0].strand, exons[0].score
  if thick_start is None or thick_stop is None :
    thick_start = thick_stop = start
  itemRgb = 0
  blockCount = len(exons)
  blockSizes, blockStarts = [], []
  for e in exons:
    blockSizes.append(len(e))
    blockStarts.append(e.start - start)
  return Bed12([chr, start, stop, id, score, strand, thick_start, thick_stop, itemRgb, blockCount, tuple(blockSizes), tuple(blockStarts)])

def sub(a, b): # a - b
  if a.chr != b.chr : return [a]
  if a.start >= b.stop or a.stop <= b.start : return [a]
  out = []
  if a.start < b.start : out.append(a(stop = b.start))
  if b.stop < a.stop: out.append(a(start = b.stop))
  return out

def intersect(a, b): # a & b
  out = []
  if a.chr != b.chr : return out
  if a.start >= b.stop or a.stop <= b.start : return out
  out.append(a(start = max(a.start, b.start), stop = min(a.stop, b.stop)))
  return out

def union(a, b): # a U b
  if a.chr != b.chr or a.start > b.stop or a.stop < b.start :
    return [a, b]
  return [a(start = min(a.start, b.start), stop = max(a.stop, b.stop))]

def test():
  return 1


