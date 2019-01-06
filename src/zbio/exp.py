'''
Expression profile processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''
def cmp3(a, b):
  ''' cmp for python3'''
  if a is None:
    if b is None: return 0
    else: return -1
  return bool(a > b) - bool(a < b)

import math
class Exp(): #values for one gene/trans/probe
  '''expression values for one gene/transcript
  '''
  def __init__(self, id, data, anno = ''): #sample and expression list
    self.id = id
    #self.sample = sample
    self.data = data
    self.anno = anno
    self.value = [] #sort index
  def __str__(self, showanno = False, sep = '\t'):
    s = self.id + sep
    if showanno : s += self.anno + sep
    return s + sep.join(map(str, self.data))
  def string(self, showanno = False, sep = '\t'):
    return self.__str__(showanno, sep)
  def __repr__(self):
    return self.headerline() + "\n" + str(self)
  def __len__(self):
    return len(self.data)
  def __cmp__(self, other):
    c = 0
    for i in range(len(self.value)):
      c = c or cmp(self.value[i], other.value[i])
      if c != 0: break
    return c
  def cmp(self, other):
    c = 0
    for i in range(len(self.value)):
      c = c or cmp3(self.value[i], other.value[i])
      if c != 0: break
    return c
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0

  def headerline(self, showanno = False, sep='\t'):# Header string
    sep=str(sep)
    s = 'id' + sep
    if showanno : s += 'anno' # + sep
    #s += sep.join(map(str, self.sample))
    return s
  
class Gene:
  def __init__(self, gid):
    self.id = gid
    #exp.__init__(self, gid, sample, data)
    self.trans = []
  @property
  def gid(self):
    return self.id
  def add_trans(self, t):
    self.trans.append(t)
  def headerline(self,sep='\t'):# not finished!
    sep=str(sep)
    if len(self.trans) == 0:
      return 'gid' + sep + sep.join(map(str, self.sample))
    else:
      return 'gid' + sep + 'tid' + sep + sep.join(map(str, self.sample))
  
def subarr(lst, ids = [], head = 0): # selected items in the list
  if len(ids) == 0 : return lst[head:]
  return [lst[i] for i in ids]
  #a = []
  #for i in ids:
        #a.append(lst[i])
  #return a
def gtexp_iter(expfile, gi = 0, ti = 1, sep = '\t', ei = [], skip = 0):
  for i in range(skip) : l = next(expfile)
  #if len(sample) == 0 :
    #l = next(expfile) #.next()
    #lst = l.strip('\n').split(sep)
    #sample = subarr(lst, ei, ti+1)
  for l in expfile:
    lst = l.strip('\n').split(sep)
    n = len(lst)
    gid = lst[gi]
    tid = lst[ti]
    data = list(map(float, subarr(lst, ei, ti+1)))
    t = Exp(tid, data, anno = gid)
    t.gid = gid
    yield t
def gtexp_load(expfile, gi = 0, ti = 1, sep = '\t', ei = [], skip = 0):
  gs = {}
  for t in gtexp_iter(expfile, gi, ti, sep, ei, skip):
    if t.gid not in gs:
      gs[t.gid] = gene(t.gid)
    gs[t.gid].add_trans(t)
  return gs

def exp_iter(expfile, ii = 0, dsi = -1, annoi = -1, sep = '\t', ei = [], skip = 0, innerskip = []):
  if dsi < 0 : dsi = max(ii, annoi) + 1 # supposed data start id
  for i in range(skip) : l = next(expfile) #.expfile.next()
  if False and header :
    l = next(expfile) #.expfile.next()
    i += 1
    lst = l.strip('\n').split(sep)
    sample = subarr(lst, ei, dsi)
  for l in expfile:
    i += 1
    if i in innerskip : continue
    lst = l.strip('\n').split(sep)
    n = len(lst)
    id = lst[ii]
    anno = ''
    if annoi >= 0 : anno = lst[annoi]
    data = list(map(float, subarr(lst, ei, dsi)))
    e = exp(id, data, anno)
    yield e

class Profile():
  '''expression profile of a list of genes
  '''
  def __init__(self):
    self.exps = {}
    self.ids = []
  def add_exp(self, e):
    self.exps[e.id] = e
    self.ids.append(e.id)
  def __len__(self):
    return len(self.exps)
  def __iter__(self):
    for eid in self.ids:
      yield self.exps[eid]
  def __getitem__(self, i):
    return self.exps[self.ids[i]]
  def BHcorrection(self, pid = -1, total = -1, append = False):
    lst = list(self.exps.values())
    n = len(lst)
    if total < 0: total = n
    for e in lst:
      #if len(e.value) < 1: e.value.append(1)
      e.value = [e.data[pid]]
    lst.sort()
    qc = 1
    for i in range(n-1, -1, -1):
      if lst[i].value[0] is None : q = None
      else : 
        q = float(lst[i].value[0]) * total / (i+1)
        if q > qc : q = qc
      lst[i].q = q
      if append : lst[i].data.append(q)
      qc = q
    return lst
  def write(self, outfile, showanno = False, sep = '\t'):
    for eid in self.ids:
      outfile.write(self.exps[eid].string(showanno, sep) + '\n')
      #showanno = False
  def TMM(self, i1 = 0, i2 = 1, mtrim = 0.3, atrim = 0.1): # The Trimmed Mean of M-values by edgeR, return log2 scale factor
    exps = [] # list(self.exps.values())
    #n = len(exps)
    #nmt = int(mtrim * n) + 1 # m trim 0.3
    #nat = int(atrim * n) + 1 # a trim 0.1
    N1, N2 = 0, 0
    #s = 0
    for eid in self.exps:
      e = self.exps[eid]
      N1 += e.data[i1]
      N2 += e.data[i2]
      if e.data[i1] <= 0 or e.data[i2] <= 0 : continue
      e.M = math.log(1.0 * e.data[i1] / e.data[i2], 2)
      e.A = 0.5 * math.log(e.data[i1] * e.data[i2], 2)
      #e.V = 1.0 / e.data[i1] + 1.0 / e.data[i2]
      #s += e.M
      #N1 += e.data[i1]
      #N2 += e.data[i2]
      #e.data += [m, a, v]
      exps.append(e)
      e.select = True
      e.value[0:1] = [e.M] ### sort1 = m
    #s /= len(exps)
    #print('mean of M: {}'.format(s))
    exps.sort()
    n = len(exps)
    nmt = int(mtrim * n) + 1 # m trim 0.3
    nat = int(atrim * n) + 1 # a trim 0.05

    #print('median of M: {}'.format(exps[int(len(exps))/2].M))
    for i in range(nmt): exps[i].select = False
    for i in range(n-nmt, n) : exps[i].select = False
  
    for e in exps: e.value[0:1] = [e.A] # sort2 = a
    exps.sort()
    for i in range(nat): exps[i].select = False
    for i in range(n-nat, n) : exps[i].select = False

    s = w = 0
    for e in exps:
      if not e.select : continue
      e.V = 1.0 / e.data[i1] + 1.0 / e.data[i2] - 1.0 / N1 - 1.0 / N2 # weight
      s += e.M / e.V
      w += 1 / e.V
    f = s / w
    #print f, s, w, n
    return f

class ReadDict(dict):
  '''dict for read counts
  read count -> times of the count appear
  '''
  def __init__(self, d = {}):
    dict.__init__(self)
    for i in d: self[i] = d[i]
  def sum(self):
    s = 0
    for i in self: s += i * self[i]
    return s
  def size(self):
    #s = 0 
    #for i in self: s += self[i]
    return sum(self.values())
  def value(self, n, default = 0):
    '''allow non-exist keys
    '''
    if n in self : return self[n]
    else : return default
  def mean(self):
    return 1.0 * self.sum() / self.size()
  def quantile(self, r = 0.5):
    size = self.size()
    size *= r
    #ks = list(self.keys())
    #ks.sort()
    ks = sorted(self)
    maxi = len(ks) - 1
    for i, k in enumerate(ks):
      if i == maxi : return k
      size -= self[k]
      if size < 0 : return k
      elif size == 0 : return (k + ks[i+1])/2.0
    return ks[-1]
  def geomean(self, add = 1):
    s = 0
    for i in self:
      s += math.log(i+add, 2) * self[i]
    s /= self.size()
    return 2 ** s - add
  def median(self):
    return self.quantile(0.5)
  def record(self, read, n = 1):
    if read not in self: self[read] = 0
    self[read] += n
  def merge(self, other):
    for read in other : self.record(read, other[read])
  def string(self):
    ks = sorted(self)
    s = '{'
    for k in ks:
      s += "%d:%d, " % (k, self[k])
    s = s.rstrip(', ') + '}'
    return s

        
