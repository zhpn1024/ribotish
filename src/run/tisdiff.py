from zbio import gtf, bam, ribo, stat, exp, tools, orf, fa, interval, io
import math, time, itertools
from os.path import isfile

def help():
  return "Differential TIS calling"
def set_parser(parser):
  #### basic input options ####
  parser.add_argument("-1", type=str, dest="tis1path", required=True, help="Predict result of group 1 TIS data")
  parser.add_argument("-2", type=str, dest="tis2path", required=True, help="Predict result of group 2 TIS data")
  parser.add_argument("-a", type=strlist, dest="tis1bampaths", required=True, help="Group 1 TIS enriched riboseq bam files, comma seperated")
  parser.add_argument("-b", type=strlist, dest="tis2bampaths", required=True, help="Group 2 TIS enriched riboseq bam files, comma seperated")
  parser.add_argument("-g", type=str, dest="genepath", required=True, help='Gene annotation file')
  #parser.add_argument("-f", type=str, dest="genomefapath", required=True, help="genome fasta file")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output result file")
  #### alt input options ####
  parser.add_argument("--tis1para", type=strlist, help="Input offset parameter files for group 1 bam files")
  parser.add_argument("--tis2para", type=strlist, help="Input offset parameter files for group 2 bam files")
  parser.add_argument("--geneformat", type=str, default='auto', help="Gene annotation file format (gtf, bed, gpd, gff, default: auto)")
  parser.add_argument("--maxNH", type=int, default=5, help="Max NH value allowed for bam alignments (default: 5)")
  parser.add_argument("--minMapQ", type=float, default=1, help="Min MapQ value required for bam alignments (default: 1)")
  parser.add_argument("--secondary", action="store_true", help="Use bam secondary alignments")
  parser.add_argument("--nocompatible", action="store_true", help="Do not require reads compatible with transcript splice junctions")
  parser.add_argument("--compatiblemis", type=int, default=2, help="Missed bases allowed in reads compatibility check")

  # scatter plot result output
  parser.add_argument("--plotout", type=str, help="Scatter plot output pdf file")
  parser.add_argument("--figsize", type=int2, default=(8,8), help="Scatter plot figure size (default: 8,8)")
  # other options
  parser.add_argument("--qi", type=int, default=15, help="Index of TIS q value (default: 15)")
  parser.add_argument("-f", type=float, dest="foldchange", default=1.5, help="Minimum fold change threshold (default: 1.5)")
  parser.add_argument("--ipth", type=float, default=0.05, help="Input TIS p value threshold (default: 0.05)")
  parser.add_argument("--iqth", type=float, default=0.05, help="Input TIS q value threshold (default: 0.05)")
  parser.add_argument("--opth", type=float, default=0.05, help="Output TIS diff p value threshold (default: 0.05)")
  parser.add_argument("--oqth", type=float, default=0.05, help="Output TIS diff q value threshold (default: 0.05)")
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")
  parser.add_argument("-v", "--verbose", action="count", help="Increase output verbosity")
  
def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')
def int2(s):
  '''Convert comma seperated string to tuple
  '''
  lst = eval('['+s+']')
  return tuple(map(int, lst))


#profile = exp.profile()
title = ['TIS1count', 'TIS1sig', 'TIS2count', 'TIS2sig']

def sig(data):
  '''if data fit thresholds
  '''
  return data[1] < pth and data[2] < qth
def get_tis(gp):
  '''get TIS genome position
  '''
  lst = gp.split(':')
  lst2 = lst[1].split('-')
  if lst[2] != '-' : lst[1] = lst2[0]
  else : lst[1] = lst2[-1]
  return ':'.join(lst)
def run(args):
  '''Main function for differential TIS
  '''
  global pth, qth, tis1bampaths, tis2bampaths, tis1offdict, tis2offdict, compatible, compatiblemis
  pth, qth = args.pth, args.qth
  tis1bampaths = args.tis1bampaths
  tis2bampaths = args.tis2bampaths
  ribo.maxNH, ribo.minMapQ, ribo.secondary = args.maxNH, args.minMapQ, args.secondary
  compatible = not args.nocompatible
  compatiblemis = args.compatiblemis

  if len(tis1bampaths) == 0 or len(tis2bampaths) == 0 :
    print('Missing bam file input!')
    exit(1)
  tis1offdict = find_offset(args.tis1bampaths, args.tis1para)
  tis2offdict = find_offset(args.tis2bampaths, args.tis2para)
  tis1file = open(args.tis1path, 'r')
  tis2file = open(args.tis2path, 'r')
  compatible = not args.nocompatible
  if args.verbose : print("Loading data...")
  trans = {}
  gname = {}
  t1, t2 = {}, {}
  n = 0
  for l in tis1file :
    lst = l.strip().split()
    try : tis = (lst[1], int(lst[6]))
    except : continue
    cnt, pval, qval =int(lst[10]), float(lst[11]), float(lst[args.qi])
    #tis, ip, cnt, pval, qval = autodetect(lst)
    #if pval >= pth or qval >= qth : break ## Should be ordered
    t1[tis] = cnt, pval, qval
    #if lst[0] not in g1 : g1[
    #trans[tis] = lst[0]
    if sig(t1[tis]) : n += 1
    lst[4] = get_tis(lst[4])
    gname[tis] =  '\t'.join(lst[:9]) # information for the TIS
  if args.verbose : print("{} TISs in 1.".format(n))
  n = 0
  for l in tis2file :
    lst = l.strip().split()
    try : tis = (lst[1], int(lst[6]))
    except : continue
    cnt, pval, qval =int(lst[10]), float(lst[11]), float(lst[args.qi])
    #if pval >= pth or qval >= qth : break
    t2[tis] = cnt, pval, qval
    #trans[tis] = lst[0]
    if sig(t2[tis]) : n += 1
    lst[4] = get_tis(lst[4])
    gname[tis] = '\t'.join(lst[:9])
  if args.verbose : print("{} TISs in 2.".format(n))
  profile = exp.Profile()
  profile2 = exp.Profile() # for TMM
  g1, g2 = {}, {} # TIS genes need to be analyzed
  cm = [[],[]] # Common TIS for plot
  for tis in t1 :
    if tis in t2 :
      cm[0].append(t1[tis][0])
      cm[1].append(t2[tis][0])
      values = [t1[tis][0], sig(t1[tis]), t2[tis][0], sig(t2[tis])]
      if values[1] or values[3] :
        e = exp.Exp(gname[tis], values)
        profile.add_exp(e)
        if values[1] and values[3] : profile2.add_exp(e)
    #a1.append(math.log(t1[tis]+1, 2))
    #a2.append(math.log(t2[tis]+1, 2))
    elif sig(t1[tis]) : 
      if tis[0] not in g2 : g2[tis[0]] = []
      g2[tis[0]].append(tis[1])
    #print tis, gene[tis]

  for tis in t2 :
    if tis in t1 : pass
    elif sig(t2[tis]) : 
      if tis[0] not in g1 : g1[tis[0]] = []
      g1[tis[0]].append(tis[1])
  if args.verbose : print("{} TISs in both.".format(len(profile2)))
  if args.verbose : print ('Estimate scale factor...')
  f = profile2.TMM(i1 = 0, i2 = 2)
  if args.verbose : print ('TMM f = {}'.format(f))

  if args.verbose : print ("Reading bams...")
  u1, u2 = [[],[]], [[],[]] # Unique TIS lists for plot
  trans_iter = io.transIter(args.genepath, fileType = args.geneformat, verbose = args.verbose)
  para_iter = transPara(trans_iter, g1, g2)
  if args.numProc <= 1 : pred_iter = itertools.imap(_get_tis, para_iter)
  else : 
    from multiprocessing import Pool
    pool = Pool(processes = args.numProc - 1)
    pred_iter = pool.imap_unordered(_get_tis, para_iter, chunksize = 20)
  for result in pred_iter: 
    tid, r1, r2 = result
    for i, c in r1 : 
      tis = tid, i
      values = [c, False, t2[tis][0], True]
      e = exp.Exp(gname[tis], values)
      profile.add_exp(e)
    for i, c in r2 : 
      tis = tid, i
      values = [t1[tis][0], True, c, False]
      e = exp.Exp(gname[tis], values)
      profile.add_exp(e)
  if args.verbose : print("{} TISs in total.".format(len(profile)))
  if args.verbose : print ('Diff test...')
  exps = profile.exps.values()
  r = 2 ** (-f)
  p = 1 / (r + 1)
  for e in exps:
    x, y = e.data[0], e.data[2]
    n = x + y
    if x == 0 : fc, alt = 'INF', 'less'
    elif y == 0 : fc, alt = 0, 'greater'
    else : 
      fc = y / (1.0 * x * r) # / y
      if r * x < y : alt = 'less'
      elif r * x > y : alt = 'greater'
    #e.data.append(stat.ACtest(e.data[0], e.data[2], r, alt = alt))
    #e.data.append(stat.FCtest(e.data[0], e.data[2], r, FC, alt = alt))
    e.data.append(fc)
    e.data.append(stat.binom_test(x, n, p = p, alt = alt))#stat.FCtest(e.data[0], e.data[2], r, FC, alt = alt))

  result = profile.BHcorrection(5) # for binom test
  for e in exps:
    e.data.append(e.q)

  if args.verbose : print ('Output...')
  outfile = open(args.output, 'w')
  s = "Gid\tTid\tSymbol\tGeneType\tGenomePos\tStartCodon\tStart\tStop\tTisType\t"
  s += '\t'.join(title)
  s += '\tFoldChange\tDiffPvalue\tDiffQvalue\n'
  outfile.write(s)
  for e in profile :
    fc = e.data[4]
    e.is_q = e.is_fc = True
    if fc != 'INF' and fc != 0 and max(fc, 1/fc) < args.foldchange : 
      e.is_fc = False
      #continue
    if e.data[5] > args.opth or e.data[6] > args.oqth : 
      e.is_q = False
      #continue
    if e.is_q and e.is_fc : 
      outfile.write(str(e)+'\n')
  #profile.write(outfile)

  # Plot
  if args.plotout != '' :
    if args.verbose : print ("Ploting...")
    from zbio import plot
    plot.figure(figsize = args.figsize)
    qd1 = [math.log(e.data[0]+1,2) for e in exps if e.is_q and e.is_fc]
    qd2 = [math.log(e.data[2]+1,2) for e in exps if e.is_q and e.is_fc]
    pd1 = [math.log(e.data[0]+1,2) for e in exps if e.is_q and not e.is_fc]
    pd2 = [math.log(e.data[2]+1,2) for e in exps if e.is_q and not e.is_fc]
    nd1 = [math.log(e.data[0]+1,2) for e in exps if not e.is_q]
    nd2 = [math.log(e.data[2]+1,2) for e in exps if not e.is_q]
    plot.scatter(qd1, qd2, alpha=0.2, edgecolors='none', color='r', label='q < {} & FC > {}'.format(args.oqth, args.foldchange))
    plot.scatter(pd1, pd2, alpha=0.2, edgecolors='none', color='y', label='q < {} & FC <= {}'.format(args.oqth, args.foldchange))
    plot.scatter(nd1, nd2, alpha=0.2, edgecolors='none', color='g', label='q >= {}'.format(args.oqth))
    plot.legend(loc='upper left', frameon=False)
    d = - f / 2
    m1 = max(qd1+pd1+nd1)
    m2 = max(qd2+pd2+nd2)
    m = min(m1, m2)
    plot.plot([0-d,m-d],[0+d,m+d], color='k', linestyle = ':')
    d = - (f + math.log(args.foldchange, 2)) / 2
    plot.plot([0-d,m-d],[0+d,m+d], color='r', linestyle = ':')
    d = - (f - math.log(args.foldchange, 2)) / 2
    plot.plot([0-d,m-d],[0+d,m+d], color='r', linestyle = ':')
    plot.save(args.plotout)
    '''
    exps = profile2.exps.values()
    plot.figure(figsize = args.figsize)
    plot.axhline(f)
    ms = [e.M for e in exps if e.select]
    aa = [e.A for e in exps if e.select]
    plot.scatter(aa, ms, alpha=0.1, edgecolors='none',color='r')
    ms = [e.M for e in exps if not e.select]
    aa = [e.A for e in exps if not e.select]
    plot.scatter(aa, ms, alpha=0.1, edgecolors='none',color='b')
    plot.save('matmm.pdf')
'''

    
def transPara(trans_iter, g1, g2):
  '''Generate parameters (trans, g1t, g2t) for function _get_tis()
  '''
  for t in trans_iter :
    g1t = g2t = None
    if t.id in g1 : g1t = g1[t.id]
    if t.id in g2 : g2t = g2[t.id]
    if g1t is None and g2t is None : continue
    yield t, g1t, g2t

def _get_tis(ps) : 
  '''get non-significant TIS counts
  '''
  t, g1t, g2t = ps
  r1, r2 = [], []
  if g1t is not None :
    t1tis = ribo.multiRibo(t, tis1bampaths, offdict = tis1offdict, compatible = compatible, mis = compatiblemis)
    #ttis = ribo.Ribo(merge, tisbam1file, offdict = tis1offdict, compatible = False)
    for i in g1t : r1.append((i, t1tis.cnts[i]))
  if g2t is not None :
    t2tis = ribo.multiRibo(t, tis2bampaths, offdict = tis2offdict, compatible = compatible, mis = compatiblemis)
    for i in g2t : r2.append((i, t2tis.cnts[i]))
  return t.id, r1, r2

def find_offset(bampaths, para):
  '''Get offset data for given bam data
  '''
  offlist = [None] * len(bampaths)
  for i, path in enumerate(bampaths):
    if para is None or i > len(para)-1 or para[i] == '' : path = path + '.para.py'
    else : path = para[i]
    if isfile(path) : 
      vessel = {}#'offdict': None}
      exec(open(path).read(), vessel) # execfile(path, vessel)
      offlist[i] = vessel['offdict']
  for i, od in enumerate(offlist):
    if od is None :
      print('No offset parameter file found for %s. Using default offset (12). ' %  bampaths[i])
  return offlist
    
