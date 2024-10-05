# -*- coding: utf-8 -*-
'''
File processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

'''
读文件迭代器，自动识别gz，自动strip和split，支持跳过前N行
example: 
for lst in io.splitIter('infile.txt.gz', skip=1):
  print(lst[0])
'''
def splitIter(filePath, sep = '\t', gz = False, skip = 0, title = None, encoding = None, maxsplit = -1):
  if type(filePath) is str :
    filePath = filePath.strip()
    if filePath.split('.')[-1].lower() == 'gz' : gz = True
    if gz :
      import gzip
      try: infile = gzip.open(filePath, 'rt', encoding=encoding)
      except: infile = gzip.open(filePath)
    else :
      try: infile = open(filePath, 'r', encoding=encoding)
      except: infile = open(filePath)
  else : infile = filePath
  for i in range(skip):
    l = next(infile)
  if title is not None : 
    l = next(infile)
    #if gz: l = l.decode()
    lst = l.rstrip('\n').split(sep, maxsplit)
    title[:] = lst
  for l in infile : 
    #if gz: l = l.decode()
    lst = l.rstrip('\n').split(sep, maxsplit)
    yield lst

'''
类似splitIter，但不split
'''
def lineIter(filePath, gz = False, skip = 0, encoding = None, strip = False):
  if type(filePath) is str :
    filePath = filePath.strip()
    if filePath.split('.')[-1].lower() == 'gz' : gz = True
    if gz :
      import gzip
      try: infile = gzip.open(filePath, 'rt', encoding=encoding)
      except: infile = gzip.open(filePath)
    else :
      try: infile = open(filePath, 'r', encoding=encoding)
      except: infile = open(filePath)
  else : infile = filePath
  for i in range(skip):
    l = next(infile)
  for l in infile :
    #if gz: l = l.decode()
    if strip: l = l.rstrip('\n')
    yield l

'''
转录本注释迭代器，返回每个转录本对象，参考bed, gtf
example:
for trans in io.transIter('anno.gtf.gz'):
  print(trans.id, trans.gid, trans.cdna_length())
'''
def transIter(filePath, fileType = 'auto', gz = False, **kwargs):
  '''yield all transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rt')
  else : infile = open(filePath, 'r')
  if fileType == 'auto' : fileType = suffixType(filePath, gz)
  if fileType == 'bed' :
    from . import bed
    return bed.bed12_iter(infile, **kwargs)
  elif fileType == 'gtf' :
    from . import gtf
    return gtf.gtftrans_iter(infile, **kwargs)
  elif fileType == 'gff' :
    from . import gtf
    return gtf.gtftrans_iter(infile, gff = True, **kwargs)
  elif fileType == 'gpd' :
    from . import bed
    return bed.gpd_iter(infile, **kwargs)
  else : raise IOError('Unknown trans file type: {}'.format(fileType))

'''
基因注释迭代器，返回每个基因对象，参考bed, gtf
example:
for gene in io.geneIter('anno.gtf.gz'):
  print(gene.id, gene.chr, gene.start, gene.strand, gene.genetype)
'''

def geneIter(filePath, fileType = 'auto', gz = False, **kwargs):
  '''yield all transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rt')
  else : infile = open(filePath, 'r')
  if fileType == 'auto' : fileType = suffixType(filePath, gz)
  if fileType == 'bed' :
    from . import bed
    return bed.bed12_iter(infile, **kwargs)
  elif fileType == 'gtf' :
    from . import gtf
    return gtf.gtfgene_iter(infile, **kwargs)
  elif fileType == 'gff' :
    from . import gtf
    return gtf.gtfgene_iter(infile, gff = True, **kwargs)
  elif fileType == 'gpd' :
    from . import bed
    return bed.gpdGeneIter(infile, **kwargs)
  else : raise IOError('Unknown trans file type: {}'.format(fileType))

'''
从基因注释中提取指定转录本，不适合大批量操作
example:
trans = io.transFetch('anno.gtf.gz', 'ENST000111')
'''

def transFetch(filePath, tid, fileType = 'auto', gz = False, **kwargs):
  '''fetch given transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rt')
  else : infile = open(filePath, 'r')
  if fileType == 'auto' : fileType = suffixType(filePath, gz)
  if fileType == 'bed' :
    from . import bed
    return bed.bed12_fetch(infile, id = tid, **kwargs)
  elif fileType == 'gtf' :
    from . import gtf
    gene, trans = gtf.fetch_gtf(infile, tid = tid, **kwargs)
    return trans[tid]
  elif fileType == 'gff' :
    from . import gtf
    gene, trans =  gtf.fetch_gtf(infile, tid = tid, gff = True, **kwargs)
    return trans[tid]
  elif fileType == 'gpd' :
    from . import bed
    return bed.gpd_fetch(infile, id = tid, **kwargs)
  else : raise IOError('Unknown trans file type: {}'.format(fileType))
  
'''
从基因注释中迭代返回每个基因的代表转录本
example:
for trans in io.transSelectIter('anno.gtf.gz'):
  print(trans.id, trans.gid)
'''

def transSelectIter(filePath, fileType = 'auto', gz = False,  **kwargs):
  '''yield selected transcript for each gene in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rt')
  else : infile = open(filePath, 'r')
  if fileType == 'auto' : fileType = suffixType(filePath, gz)
  if fileType == 'bed' :
    from . import bed
    return bed.bed12SelectIter(infile, **kwargs) # no gene name/id annotation
  elif fileType == 'gtf' :
    from . import gtf
    return gtf.gtftransSelectIter(infile, **kwargs)
  elif fileType == 'gff' :
    from . import gtf
    return gtf.gtftransSelectIter(infile, gff = True, **kwargs)
  elif fileType == 'gpd' :
    from . import bed
    return bed.gpdSelectIter(infile, **kwargs)
  else : raise IOError('Unknown trans file type: {}'.format(fileType))

'''
从文件名获取后缀文件类型，包括bed, gtf, gff, gpd，支持带gz
'''
def suffixType(filePath, gz = False):
  i = -1
  if gz : i = -2
  suffix = filePath.split('.')[i].lower()
  if suffix in ('bed', 'bed12') : fileType = 'bed'
  elif suffix in ('gtf', ) : fileType = 'gtf'
  elif suffix in ('gff', 'gff3') : fileType = 'gff'
  elif suffix in ('gpd', 'genepred') : fileType = 'gpd'
  else : raise IOError('Unknown trans file format: {}'.format(filePath))
  return fileType

'''
自动tab连接变量和列表，最后元素为'\n'可合并换行
example:
outfile.write(io.tabjoin(arr1, n1, name2, dict2, '\n'))
'''
def tabjoin(a, *args): #, sep = '\t'):
  sep = '\t'
  if type(a) == str: arr = [a]
  elif hasattr(a, '__iter__') : #s = sep.join(map(str, a))
    arr = list(map(str, a))
  else : #s = str(a)
    arr = [str(a)]
  for a in args :
    if type(a) == str: arr.append(a)
    elif hasattr(a, '__iter__') : arr += map(str, a)
    else : arr.append(str(a))
  if len(arr) > 1 and arr[-1] == '\n':
    arr[-2] = '{}{}'.format(arr[-2], arr[-1])
    arr = arr[:-1]
  return sep.join(arr)

'''
并行处理输入的迭代器，p < 1 时为单线程，结果为skip时跳过
需自定义process函数process(i, args)，适用于process耗时的场景

'''
def multiProcIter(input, process, args, p = 2, skip = None):
  p = int(p)
  if p < 1:
    for i in input:
      r = process(i, args)
      if r is not skip: yield r
  else:
    from multiprocessing import Pool
    pool = Pool(p)
    jobs, res = [], []
    for i in input: # range(n):
      if len(jobs) == p:
        if len(res) != '':
          for r in res:
            if r is not skip: yield r
        res = [j.get() for j in jobs]
        jobs = []
      jobs.append(pool.apply_async(process, (i, args, )))
    if len(res) != '':
      for r in res:
        if r is not skip: yield r
    res = [j.get() for j in jobs]
    for r in res:
      if r is not skip: yield r
    pool.close()
    pool.join()

class OrderedList:

  def __init__(self, data = None):
    if data is None: data = []
    self.data = data
    self.data.sort()
    self.i = 0
    self.l = len(self.data)

  def __len__(self):
    return self.l

  def check(self):
    self.data.sort()
    self.i = 0
    self.l = len(self.data)

  def find(self, a):
    while self.i < self.l and self.data[self.i] < a:
      self.i += 1
    if self.i == self.l: return False
    if self.data[self.i] == a: return True
    else: return False


