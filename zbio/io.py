def splitIter(filePath, sep = '\t', gz = False):
  if gz :
    import gzip
    infile = gzip.open(filePath, 'rb')
  else : infile = open(filePath, 'r')
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
    import bed
    return bed.bed12_iter(infile, **kwargs)
  elif fileType == 'gtf' :
    import gtf
    return gtf.gtftrans_iter(infile, **kwargs)
  elif fileType == 'gff' :
    import gtf
    return gtf.gtftrans_iter(infile, gff = True, **kwargs)
  elif fileType == 'gpd' :
    import bed
    return bed.gpd_iter(infile, **kwargs)
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
    import bed
    return bed.bed12_fetch(infile, id = tid, **kwargs)
  elif fileType == 'gtf' :
    import gtf
    gene, trans = gtf.fetch_gtf(infile, tid = tid, **kwargs)
    return trans[tid]
  elif fileType == 'gff' :
    import gtf
    gene, trans =  gtf.fetch_gtf(infile, tid = tid, gff = True, **kwargs)
    return trans[tid]
  elif fileType == 'gpd' :
    import bed
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
    import bed
    return bed.bed12SelectIter(infile, **kwargs) # no gene name/id annotation
  elif fileType == 'gtf' :
    import gtf
    return gtf.gtftransSelectIter(infile, **kwargs)
  elif fileType == 'gff' :
    import gtf
    return gtf.gtftransSelectIter(infile, gff = True, **kwargs)
  elif fileType == 'gpd' :
    import bed
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
