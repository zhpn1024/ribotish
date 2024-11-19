'''
Gtf format annotation processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

def add_chr(chr):
  if len(chr) <= 3 and (chr.isdigit() or chr in ('X','Y',)): chr = 'chr' + chr
  elif chr in ('M','MT','MtDNA','mitochondrion_genome','Mito') : chr = 'chrM' 
  elif len(chr) == 2 and chr[0].isdigit(): chr = 'chr' + chr
  elif chr in ('I','II','III','IV','V','VI','VII','VIII','IX','XI','XII','XIII','XIV','XV','XVI'): chr = 'chr' + chr
  return chr
def rm_chr(chr):
  if chr[3:].isdigit() or chr in ('chrX','chrY'): chr = chr[3:]
  elif chr == 'chrM' : chr = 'MT'
  return chr
def changechr(chr):
  '''change between two chr versions
  '''
  if chr.isdigit() or chr in ('X','Y','M'): return 'chr' + chr
  elif chr == 'MT' : return 'chrM'
  elif chr == 'chrM' : return 'MT'
  elif chr[0:3] == 'chr' : return chr[3:]
  else : return chr
def cmp3(a, b):
  ''' cmp for python3'''
  return (a > b) - (a < b)

class Exon:
  '''exon in gtf format, 0 based, similar as bed
  '''
  def __init__(self, lst, gff = False, addchr = False):
    self.chr, self.strand, self.start, self.stop = lst[0], lst[6], int(lst[3]) - 1, int(lst[4])
    if addchr and self.chr[0:3] != 'chr' :
      if self.chr.isdigit() or self.chr in ('X','Y','M',): self.chr = 'chr' + self.chr
      elif self.chr == 'MT' : self.chr = 'chrM' 
    self.type, self.score, self.attrstr = lst[2], lst[5], lst[8]
    self.gff, self.addchr, self.frame = gff, addchr, lst[7]
    self.lst = lst
    self.exons = [self]
  def __repr__(self):
    return self.chr + ':'+str(self.start)+'-'+str(self.stop)+':'+self.strand
  def __str__(self, sep = '\t', new = False):
    if not new : return sep.join(self.lst)
    lst = [self.chr, self.lst[1], self.type, str(self.start+1), str(self.stop), str(self.score), self.strand, str(self.frame), self.attrstr]
    return sep.join(lst)
  def string(self, sep = '\t', new = False):
    return self.__str__(sep = sep, new = new)
  def short_str(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  def __cmp__(self, other):
    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
  def cmp(self, other):
    return cmp3(self.chr,other.chr) or cmp3(self.start,other.start) or cmp3(self.stop,other.stop)
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0
  def __len__(self): 
    return self.stop - self.start
  def attr(self, key):
    '''get attributes in the last column
    '''
    s = self.attrstr
    if self.gff : 
      p1 = s.find(key+'=')
      if p1 < 0 : p1 = s.find(key+':')
      if p1 < 0 : return ''
      else : p1 += len(key) #+ 1
      p2 = s.find(';', p1)
      if p2 == -1 : p2 = None
      if s[p1] == '=' : return s[p1+1:p2]
      elif s[p1] == ':' : return s[p1+1:p2].split(',')[0]
      else : return ''
    p1 = s.find(key + ' ')
    if p1 < 0 : return ''
    else : p1 += len(key) + 1
    p2 = s.find(';', p1)
    if p2 == -1 : p2 = None
    try: return eval(s[p1:p2])
    except: return s[p1:p2]
  @property
  def length(self):
    return len(self)
  def cdna_length(self): 
    return len(self)
  def is_reverse(self): 
    return self.strand=='-'
  @property
  def anti_strand(self):
    if self.strand == '+': return '-'
    elif self.strand == '-': return '+'
    else : return '.'
  @property
  def gid(self):
    try : return self.gid_c # gid cache
    except : 
      if self.gff and self.attr('GeneID') != '' : self.gid_c = self.attr('GeneID')
      else : self.gid_c = self.attr('gene_id')
    return self.gid_c
  @property
  def tid(self):
    try : return self.tid_c
    except : 
      self.tid_c = self.attr('transcript_id')
      if self.gff: # and self.tid_c == '':
        if self.lst[2] in ('transcript', 'mRNA', 'tRNA', 'rRNA'):
          self.tid_c = self.id
        else:
          p = self.attr('Parent')
          if p.startswith('rna-'): self.tid_c = p[4:]
          elif p.startswith('gene-'): self.tid_c = p[5:]
          elif p.startswith('rna'): self.tid_c = p
          elif p.startswith('gene'): self.tid_c = p
      return self.tid_c
  @property
  def symbol(self):
    try : return self.sym_c
    except : 
      self.sym_c = self.attr('gene_name')
      if self.gff and self.sym_c == '': self.sym_c = self.attr('gene')
    return self.sym_c
  @property
  def id(self):
    try: return self._id
    except:
      self._id = self.attr('exon_id')
      if self._id == '': self._id = self.attr('ID')
      return self._id
  @id.setter
  def id(self, value):
    self._id = value
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start
  @property
  def genetype(self):
    a = self.attr('gene_biotype')
    if a != '' : return a
    a = self.attr('gene_type')
    if a != '' : return a
    if self.tid[0:3] == 'NM_' : return "protein_coding"
    return a
  def short_str(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  def __call__(self, **args):
    '''New exon that modified from self.
    '''
    lst = self.lst[:]
    new = Exon(lst, self.gff, self.addchr)
    for k in args:
      new.__dict__[k]=args[k]
    return new
  def __sub__(self, other):
    return sub(self, other)
  def intersect(self, other):
    return intersect(self, other)
  def union(self, other):
    return union(self, other)
  def is_contain(self, p, strict = False): #if i in exon
    if self.start <= p <= self.stop :
      if not strict : return True
      else : return p != self.end3
    else : return False
  def is_exon(self, p, strict = False): 
    return self.is_contain(p, strict = strict)
  def is_intron(self, p, strict = False): # False
    return False
  def is_upstream(self, p, strict = False) :
    if self.strand != '-' : return p < self.start
    else : return p > self.stop
  def is_downstream(self, p, strict = False) :
    if self.strand != '-' : 
      if strict : return p >= self.stop
      else : return p > self.stop
    else : 
      if strict : return p <= self.start
      else : return p < self.start
  def flank_pos(self, p, flank = 100):
    if flank < 0: return None
    if self.start - flank <= p <= self.start:
      if self.strand != '-': return p - self.start # upstream
      else: return self.start - p + self.cdna_length()
    if self.stop <= p <= self.stop + flank:
      if self.strand != '-': return p - self.stop + self.cdna_length()
      else: return self.stop - p # upstream
    return None

  def cdna_pos(self, p, strict = False, flank = 0): # exon only
    if self.is_contain(p, strict): return abs(p-self.end5)
    elif flank > 0: return self.flank_pos(p, flank)
    else: return None
  def genome_pos(self, p, bias = 1): # exon only
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    if not self.is_reverse() : return self.start + p
    else : return self.stop - p
  def copy(self, other):
    '''copy everything from other exon
    '''
    self.chr, self.strand, self.start, self.stop = other.chr, other.strand, other.start, other.stop
    self.type, self.score, self.attrstr = other.type, other.score, other.attrstr
    #self.attr = other.attr
    self.gff, self.addchr, self.frame = other.gff, other.addchr, other.frame
    self.lst = other.lst
  #@property
  #def exons(self):
    #return [self]
  def cds_start(self, cdna = False): return None
  def cds_stop(self, cdna = False): return None
    
def attr(s, gff = False):
  '''get attr dict from attr string
  '''
  if type(s) == dict: return s
  if gff : 
    l = s.strip(';').split(';')
    a = {}
    for att in l:
      l2 = att.split('=')
    #print att, l2
      a[l2[0]] = l2[1]
      if l2[0] == "Dbxref" :
        l3 = l2[1].split(',')
        for att2 in l3:
          l4 = att2.split(':')
          a[l4[0]] = l4[1]
    return a
  l = s.strip(';').split('; ')
  a = {}
  for att in l:
    att = att.strip()
    p = att.find(' ')
    a[att[0:p]] = eval(att[p+1:])
  return a
def attrstr(a) :  # Not tested!
  s = ''
  for k in a :
    s += k + ' '
    if type(a[k]) == str : s += '"' + a[k] + '"'
    else : s += str(a[k])
    s += '; '
  return s
def tup2com(t): #tuple to comma string
  if type(t)==str: return t
  return ','.join(map(str,t))+','

class gtfTrans(Exon):
  '''gtf transcript
  '''
  def __init__(self, lst, gff = False, addchr = False): 
    Exon.__init__(self, lst, gff, addchr)
    if self.tid != '': self.id = self.tid
    elif self.id != '': self.tid_c = self.id
    #self.type = lst[1]
    self.exons = []
    self.cds = []
    self.utr = []
    self.start_codon = None
    self.stop_codon = None
    self.other = []
  def add_exon(self, e): 
    if e.strand != self.strand : 
      print('Wrong exon strand: {} {} {} {}'.format(e.gid, e.tid, e.strand, self.strand))
      e.strand = self.strand
    if e.type == 'exon': self.exons.append(e)
    elif e.type == 'CDS': self.cds.append(e)
    elif e.type == 'UTR': self.utr.append(e)
    elif e.type == 'start_codon': self.start_codon = e
    elif e.type == 'stop_codon': self.stop_codon = e
    else: self.other.append(e)
    #self.check()
  def __repr__(self):
    s = self.symbol + " transcript_id " + self.id + ', ' + str(len(self.exons)) + " exons, " + Exon.__repr__(self)
    return s
  def __str__(self, sep = '\t', new = False):
    arr = self.exons + self.cds + self.utr + self.other
    if self.start_codon is not None : arr.append(self.start_codon)
    if self.stop_codon is not None : arr.append(self.stop_codon)
    arr.sort(reverse = self.is_reverse())
    s = ''
    if self.type == 'transcript' : s += Exon.__str__(self, new = new) + '\n'
    for e in arr : 
      s += e.__str__(sep = sep, new = new) + '\n'
    return s.strip()
  def string(self, sep = '\t', new = False):
    return self.__str__(sep = sep, new = new)
  def allRecords(self):
    arr = self.exons + self.cds + self.utr + self.other
    if self.start_codon is not None : arr.append(self.start_codon)
    if self.stop_codon is not None : arr.append(self.stop_codon)
    arr.sort(reverse = self.is_reverse())
    if self.type == 'transcript' : arr[0:0] = [self]
    return arr
  def genePredStr(self, geneName = True, extended = False, slice = None):
    '''generate genepred format
    '''
    lst = []
    if geneName : lst += [self.gid]
    start, stop, thickstart, thickstop = self.start, self.stop, self.thick_start, self.thick_stop
    if slice is None :
      exons = self.exons
    else :
      if start < slice[0] : start = slice[0]
      if thickstart < slice[0] : thickstart = slice[0]
      if stop > slice[1] : stop = slice[1]
      if thickstop > slice[1] : thickstop = slice[1]
      exons = []
      for e in self.exons:
        if e.stop <= slice[0] : continue
        if e.start >= slice[1] : continue
        exons.append(e(start = max(e.start, slice[0]), stop = min(e.stop, slice[1])))
    lst += [self.id, self.chr, self.strand, start, stop, thickstart, thickstop, len(exons)]
    exonStarts = [x.start for x in exons]
    exonStops = [x.stop for x in exons]
    exonStarts.sort()
    exonStops.sort()
    lst += [tup2com(exonStarts), tup2com(exonStops)]
    return '\t'.join(map(str, lst))
  def check(self, final = True):
    if final and len(self.exons) == 0:
      if len(self.cds) > 0: self.exons += self.cds
    for e in self.exons:
      if e.start < self.start: self.start = e.start
      if e.stop > self.stop: self.stop = e.stop
    self.exons.sort(reverse = self.is_reverse())
    self.cds.sort(reverse = self.is_reverse())
  #@property
  def cds_start(self, cdna = False): 
    if self.start_codon is None and len(self.cds) == 0 : return None
    if self.start_codon is None :
      self.cds.sort(reverse = self.is_reverse())
      gs = self.cds[0].end5
      if self.stop_codon is not None : 
        cds2 = self.cds_stop(cdna = True)
        cds1 = self.cdna_pos(gs)
        cs = cds2 - (cds2 - cds1) // 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      try : frame = int(self.cds[0].frame)
      except : frame = 0
      if frame > 0 : 
        cs = self.cdna_pos(gs) + frame
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      if cdna : return self.cdna_pos(gs)
      else : return gs
    sc = self.start_codon
    if len(sc) != 3 : 
      for e in self.exons:
        if sc.end5 == e.end5 : 
          cs = self.cdna_pos(sc.end3) - 3
          if cdna : return cs
          else : return self.genome_pos(cs, 1)
        if sc.end3 == e.end3 :
          break
      else : 
        print('Start codon error : %s %s %s' % (self.gid, self.id, sc.short_str()))
    gs = sc.end5
    if cdna : return self.cdna_pos(gs)
    else : return gs
  def cds_stop(self, cdna = False):
    if self.stop_codon is None and len(self.cds) == 0 : return None
    if self.stop_codon is None :
      self.cds.sort(reverse = self.is_reverse())
      gs = self.cds[-1].end3
      if self.start_codon is not None : # determine by paired start codon
        cds1 = self.cds_start(cdna = True)
        cds2 = self.cdna_pos(gs)
        cs = cds1 + (cds2 - cds1) // 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 0)
      try : frame = int(self.cds[-1].frame) # determine by frame
      except : frame = None
      if frame is not None : 
        cds1 = self.cdna_pos(self.cds[-1].end5) + frame
        cds2 = self.cdna_pos(gs)
        cs = cds1 + (cds2 - cds1) // 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      if cdna : return self.cdna_pos(gs) #cs
      else : return gs # self.genome_pos(cs, 0)
    sc = self.stop_codon # Determine by stop codon
    if len(sc) != 3 : 
      for e in self.exons:
        if sc.end5 == e.end5 : break
        elif sc.end3 == e.end3 : # bug elif 
          cs = self.cdna_pos(sc.end5) + 3
          if cdna : return cs
          else : return self.genome_pos(cs, 0)
      else : 
        print('Stop codon error : %s %s %s' % (self.gid, self.id, sc.short_str()))
    gs = sc.end3
    if cdna : return self.cdna_pos(gs)
    else : return gs
  @property
  def thick_start(self):
    if self.is_reverse(): s = self.cds_stop(cdna = False)
    else : s = self.cds_start(cdna = False)
    if s is None : return self.start
    else : return s
  @property
  def thick_stop(self):
    if self.is_reverse(): s = self.cds_start(cdna = False)
    else : s = self.cds_stop(cdna = False)
    if s is None : return self.start
    else : return s
  @property
  def itemRgb(self):
    return '0,0,0'
  @property
  def blockCount(self):
    return len(self.exons)
  @property
  def blockSizes(self):
    es = self.exons[:]
    if self.is_reverse() : es.sort()
    return ','.join([str(len(e)) for e in es])
  @property
  def blockStarts(self):
    es = self.exons[:]
    if self.is_reverse() : es.sort()
    return ','.join([str(e.start - self.start) for e in es])
  def cdna_length(self): 
    l = 0
    for e in self.exons:
      l += len(e)
    return l
  def cds_length(self): 
    try : return self.cds_stop(cdna = True) - self.cds_start(cdna = True) ##
    except: return 0
  @property
  def introns(self):
    if hasattr(self, '_introns') : return self._introns
    introns = []
    last = -1
    for i, e in enumerate(self.exons):
      if last >= 0 : 
        lst = [last, e.end5]
        lst.sort()
        it = Exon([self.chr,'','intron',lst[0]+1,lst[1],'',self.strand,'',self.attrstr], gff=self.gff)
        it.id += '_intron{}'.format(i)
        introns.append(it)
      last = e.end3
    self._introns = introns
    return introns
  def cdna_pos(self, p, strict = False, flank = 0):
    '''if strict is True, the 3' end of exon will be considered as not in the transcript,
    if strict is False, 3' end of exon will be considered as start of the next exon, 
    or transcript end (self.cdna_length()) if in the last exon.
    if flank is not 0, TSS and TTS flanking regions are returned with negative distance and
    cDNA length + distance
    '''
    if not self.is_contain(p, strict) :
      if flank <= 0: return None
      else: return self.flank_pos(p, flank)
    pos = 0
    for e in self.exons:
      if e.is_upstream(p, strict) : return None
      if e.is_contain(p, strict) : return pos + abs(p - e.end5) # e.start <= p <= e.stop:
      pos += len(e)
    return None
  def is_exon(self, p, strict = False):
    if not self.is_contain(p, strict) : return False
    for e in self.exons:
      if e.is_upstream(p, strict) : return False
      if e.is_contain(p, strict) : return True # e.start <= p <= e.stop:
    return False
  def genome_pos(self, p, bias = 1):
    '''if bias is 1, the splice junction will be mapped to 5' end of downstream exon,
    if bias is 0, the splice junction will be mapped to 3' end of the upstream exon.
    '''
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    p1 = p
    pos = self.start
    for e in self.exons:
      if len(e) - p1 >= bias:
        if self.is_reverse():
          pos = e.stop - p1
        else:
          pos = p1 + e.start
        return pos
      else:
        p1 -= len(e)
    return None
  def copytrans(self, o):
    self.exons, self.cds, self.utr, self.start_codon, self.stop_codon, self.others = o.exons, o.cds, o.utr, o.start_codon, o.stop_codon, o.others
    
class gtfGene(Exon):
  '''gtf gene
  '''
  def __init__(self, lst, gff = False, addchr = False):
    Exon.__init__(self, lst, gff, addchr)
    self.id = self.gid
    self.trans = []
    #self.type = lst[1]
  def add_trans(self, tr):
    if self.strand in ('+', '-') and tr.strand != self.strand :
      print('Inconsistent trans strand: {} {} {} {}'.format(tr.gid, tr.tid, tr.strand, self.strand))
      #tr.strand = self.strand
      self.strand = '.'
    self.trans.append(tr)
    #self.check()
  def __repr__(self):
    s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + Exon.__repr__(self)
    for t in self.trans:
      s += '\n\t' + t.__repr__()
    return s
  def check(self):
    for t in self.trans:
      t.check()
      if t.start < self.start: self.start = t.start
      if t.stop > self.stop: self.stop = t.stop
  def merge_trans(self):
    '''generate a new transcript that merge all transcript exons in the gene
    '''
    es = []
    merge = gtfTrans(self.lst)
    for t in self.trans:
      es += t.exons
    es.sort()
    me = es[0]
    for i in range(1, len(es)):
      if me.stop >= es[i].start:
        me = union(me, es[i])[0]
      else: 
        merge.add_exon(me)
        me = es[i]
    merge.add_exon(me)
    merge.check()
    return merge
  def max_cds(self):
    return selectMaxCDS(self)
    
def load_gtf(fin, filt = [], gff = False, addchr = False, verbose = False):
  genes = {}
  trans = {}
  symbol = {}
  dfilt = {}
  if len(filt) > 0: 
    for gid in filt : dfilt[gid] = True
  chr = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if len(lst) < 2 : continue
    if lst[0] != chr :
      ks = list(dfilt.keys())
      for k in ks : 
        if k in genes or k in trans or k in symbol : del(dfilt[k])
        #if genes[gid].symbol in dfilt : del(dfilt[genes[gid].symbol])
      chr = lst[0]
      if verbose : print(chr)
    #for k in dfilt:
      #if not dfilt[gid] : continue
      #if lst[8].find(k) >= 0 : break
    #else : continue
    if lst[2] == 'region' : continue
    e = Exon(lst, gff, addchr)
    if len(filt) > 0 and e.gid not in dfilt and e.tid not in dfilt and e.symbol not in dfilt: 
      continue
    if lst[2] == 'gene':
      g = gtfGene(lst, gff, addchr)
      if g.id in genes: g.trans = genes[g.id].trans
      genes[g.id] = symbol[g.symbol] = g
    elif lst[2] == 'transcript':
      t = gtfTrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = symbol[g.symbol] = g
      if t.id in trans: trans[t.id].copy(t)
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
      e = Exon(lst, gff, addchr)
      if e.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = symbol[g.symbol] = g
      if e.tid not in trans:
        t = gtfTrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
  for g in genes : genes[g].check()
  return genes, trans

def fetch_gtf(fin, gid = '', tid = '', gff = False, addchr = False):
  genes = {}
  trans = {}
  if gid == '' and tid == '' : return genes, trans
  ch = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if len(lst) < 2 : continue
    if ch != '' and ch != lst[0] : break
    if gid != '' and lst[8].find(gid) < 0 : continue
    if tid != '' and lst[8].find(tid) < 0 : continue
    if lst[2] == 'region' : continue
    e = Exon(lst, gff, addchr)
    if e.tid != tid  and e.gid != gid : continue
    if e.gid == '' : continue
    if lst[2] == 'gene':
      g = gtfGene(lst, gff, addchr)
      if g.id in genes: g.trans = genes[g.id].trans
      if g.id == gid : genes[g.id] = g
    elif lst[2] == 'transcript':
      t = gtfTrans(lst, gff, addchr)
      if t.id != tid  and t.gid != gid : continue 
      if t.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = g
      if t.id in trans: trans[t.id].copy(t)
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
      if e.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = g
      if e.tid not in trans:
        t = gtfTrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
    ch = g.chr
  for g in genes : genes[g].check()
  return genes, trans

def gtfgene_iter(fin, filt = [], gff = False, addchr = False, chrs = None, verbose = False):
  genes, genes2, trans, gidlist = {}, {}, {}, [] # genes2: genes with headline, for better output
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if len(lst) < 2 : continue
    if lst[2] == 'region' : continue
    if lst[0] != chr:
      if chrs is not None and lst[0] not in chrs and changechr(lst[0]) not in chrs: 
        #print('{}\t{}\t{}'.format(lst[0], changechr(lst[0]), chrs))
        continue
      if verbose : print(lst[0])
      for gid in gidlist:
        if gid in genes : 
          genes[gid].check()
          yield genes[gid]
      genes, genes2, trans, gidlist = {}, {}, {}, []
      chr = lst[0]
    e = Exon(lst, gff, addchr)
    if len(filt) > 0 and e.gid not in dfilt: 
        continue
    if e.gid not in genes : # check gene output
      dgs = []
      for gid in genes2 : 
        if e.start > genes2[gid].stop : 
          genes2[gid].check()
          yield genes2[gid]
          dgs.append(gid)
      for gid in dgs:
        del genes2[gid]
        del genes[gid]
    if lst[2] == 'gene':
      g = gtfGene(lst, gff, addchr)
      if g.id in genes: g.trans = genes[g.id].trans # Already exists due to disorder
      else : gidlist.append(g.id)
      genes[g.id] = g
      genes2[g.id] = g
      #if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] in ('transcript', 'mRNA', 'tRNA', 'rRNA'):
      t = gtfTrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      if t.id in trans: trans[t.id].copy(t) # Already exists due to disorder
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
      if e.gid == '' or e.tid == '' : continue
      if e.gid not in genes:
        g = gtfGene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      if e.tid not in trans:
        t = gtfTrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
  for gid in gidlist:
    if gid not in genes : continue
    genes[gid].check()
    yield genes[gid]

def selectMaxCDS(g): 
  '''select transcript with max CDS length
  '''
  maxlen, mt = 0, None
  for t in g.trans:
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
    try : cdslen = cds2 - cds1 
    except : continue
    if cdslen % 3 != 0 : continue
    if cdslen > maxlen : maxlen, mt = cdslen, t
  return mt
def gtftransSelectIter(fin, select = selectMaxCDS, **kwargs):
  for g in gtfgene_iter(fin, **kwargs):
    t = select(g)
    if t is not None : yield t
  
def gtftrans_iter(fin, filt = [], gff = False, addchr = False, chrs = None, verbose = False): # To be update
  genes, trans, trans2, tidlist = {}, {}, {}, []
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if len(lst) < 2 : continue
    if lst[2] == 'region' : continue
    if lst[0] != chr:
      if chrs is not None and lst[0] not in chrs and changechr(lst[0]) not in chrs: continue
      if verbose : print(lst[0])
      for tid in tidlist:
        if tid in trans : 
          trans[tid].check()
          yield trans[tid]
      genes, trans, trans2, tidlist = {}, {}, {}, []
      chr = lst[0]
    if lst[2] == 'gene': continue ## not gene structure
    e = Exon(lst, gff, addchr)
    if len(filt) > 0 and e.tid not in dfilt: 
        continue
    if e.tid not in trans : # check gene output
      dgs = [] # genes to delete
      for tid in trans2 : 
        if e.start > trans2[tid].stop : 
          trans2[tid].check()
          yield trans2[tid]
          dgs.append(tid)
      for tid in dgs:
        del trans2[tid]
        del trans[tid]
    if lst[2] == 'transcript':
      t = gtfTrans(lst, gff, addchr)
      if len(filt) > 0 and t.tid not in dfilt: 
        continue
      if t.id in trans: trans[t.id].copy(t)
      else : 
        #genes[t.gid].add_trans(t)
        trans[t.id] = t
        trans2[t.id] = t
        tidlist.append(t.id)
    else:
      #e = exon(lst, gff, addchr)
      if e.tid not in trans:
        t = gtfTrans(lst, gff, addchr)
        #genes[t.gid].add_trans(t)
        trans[t.id] = t
        tidlist.append(t.id)
        #print e.attr['transcript_id'], t.id
      trans[e.tid].add_exon(e)
  for tid in tidlist:
    if tid not in trans : continue
    trans[tid].check()
    yield trans[tid]
    #for t in genes[gid].trans:
      #yield t
  #return genes, trans

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

