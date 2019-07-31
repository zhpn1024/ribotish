'''
File processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

def splitIter(filePath, sep = '\t', gz = False, skip = 0, title = None):
  if type(filePath) is str :
    filePath = filePath.strip()
    if filePath.split('.')[-1].lower() == 'gz' : gz = True
    if gz :
      import gzip
      infile = gzip.open(filePath, 'rb')
    else : infile = open(filePath, 'r')
  else : infile = filePath
  for i in range(skip):
    l = next(infile)
  if title is not None : 
    l = next(infile)
    lst = l.rstrip('\n').split(sep)
    title[:] = lst
  for l in infile : 
    lst = l.rstrip('\n').split(sep)
    yield lst

def transIter(filePath, fileType = 'auto', gz = False, **kwargs):
  '''yield all transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rb')
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

def geneIter(filePath, fileType = 'auto', gz = False, **kwargs):
  '''yield all transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rb')
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


def transFetch(filePath, tid, fileType = 'auto', gz = False, **kwargs):
  '''fetch given transcript in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rb')
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
  
def transSelectIter(filePath, fileType = 'auto', gz = False,  **kwargs):
  '''yield selected transcript for each gene in gene annotation file
  '''
  if filePath.split('.')[-1].lower() == 'gz' : gz = True
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rb')
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

def tabjoin(a, *args) : #, sep = '\t'):
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
  return sep.join(arr)
