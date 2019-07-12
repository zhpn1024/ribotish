from ribotish.zbio import gtf, bam, ribo, stat, exp, tools, orf, fa, interval, io
import math, time, itertools, sys
from os.path import isfile

def help():
  return "Differential TIS calling"
def set_parser(parser):
  #### basic input options ####
  parser.add_argument("-1", type=strlist, dest="tis1paths", required=True, help="Prediction results of group 1 TIS data")
  parser.add_argument("-2", type=strlist, dest="tis2paths", required=True, help="Prediction results of group 2 TIS data")
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
  parser.add_argument("--paired", action="store_true", help="Reads are paired end")

  parser.add_argument("--l1", type=strlist, dest="tis1labels", default=[], help="Labels for group 1 TIS data")
  parser.add_argument("--l2", type=strlist, dest="tis2labels", default=[], help="Labels for group 2 TIS data")

  parser.add_argument("--nocompatible", action="store_true", help="Do not require reads compatible with transcript splice junctions")
  parser.add_argument("--compatiblemis", type=int, default=2, help="Missed bases allowed in reads compatibility check")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if annotation chr ids are not same as chr ids in bam/fasta files")
  parser.add_argument("--normcomm", action="store_true", help="Use common TISs instead of union TISs for normalization")
  parser.add_argument("--normanno", action="store_true", help="Use only annotated TISs for normalization")

  parser.add_argument("--rnaseq", type=str, help="RNASeq count input")
  parser.add_argument("--scalefactor", type=float, help="Input TIS scale factor of 2/1 (default: auto)")
  parser.add_argument("--rnascale", type=float, help="Input RNASeq scale factor of 2/1 (default: auto)")
  parser.add_argument("--chi2", action="store_true", help="Use chisquare test instead of Fisher's exact test for mRNA referenced test")
  parser.add_argument("--betabinom", action="store_true", help="Use beta-binom test instead of Fisher's exact test for mRNA referenced test")

  parser.add_argument("--export", type=str, help="Export count table for differential analysis with other tools")

  # scatter plot result output
  parser.add_argument("--plotout", type=str, help="Scatter plot output pdf file")
  parser.add_argument("--figsize", type=int2, default=(8,8), help="Scatter plot figure size (default: 8,8)")
  parser.add_argument("--plotma", type=str, help="TIS normalization MA plot output pdf file")

  # other options
  parser.add_argument("--qi", type=int, default=15, help="Index of TIS q value, 0 based (default: 15)")
  parser.add_argument("-f", type=float, dest="foldchange", default=1.5, help="Minimum fold change threshold (default: 1.5)")
  parser.add_argument("--ipth", type=float, default=0.05, help="Input TIS p value threshold (default: 0.05)")
  parser.add_argument("--iqth", type=float, default=0.05, help="Input TIS q value threshold (default: 0.05)")
  parser.add_argument("--opth", type=float, default=0.05, help="Output TIS diff p value threshold (default: 0.05)")
  parser.add_argument("--oqth", type=float, default=0.05, help="Output TIS diff q value threshold (default: 0.05)")
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")
  parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase output verbosity")
  
def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')

def int2(s):
  '''Convert comma seperated string to tuple
  '''
  lst = eval('['+s+']')
  return tuple(map(int, lst))

try: imap = itertools.imap
except: imap = map

def sig(data):
  '''if data fit thresholds
  '''
  return data[1] < ipth and data[2] < iqth
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
  global ipth, iqth, tis1bampaths, tis2bampaths, tis1offdict, tis2offdict, compatible, compatiblemis, paired
  ipth, iqth = args.ipth, args.iqth
  tis1bampaths = args.tis1bampaths
  tis2bampaths = args.tis2bampaths
  ribo.maxNH, ribo.minMapQ, ribo.secondary = args.maxNH, args.minMapQ, args.secondary
  compatible = not args.nocompatible
  compatiblemis = args.compatiblemis
  paired = args.paired
  if len(tis1bampaths) < len(args.tis1paths) or len(tis2bampaths) < len(args.tis2paths) : # == 0 :
    print('Missing bam file input!')
    exit(1)
  if args.chrmap is not None :
    chrmap = {}
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]
    bam.chrmap = chrmap
    fa.chrmap = chrmap

  global tis1bampathslist, tis2bampathslist, tis1offdictlist, tis2offdictlist
  tis1bampathslist = [s.split(';') for s in args.tis1bampaths]
  tis2bampathslist = [s.split(';') for s in args.tis2bampaths]
  if args.tis1para is None : tis1paralist = [None] * len(args.tis1paths)
  else : tis1paralist = [s.split(';') for s in args.tis1para]
  if args.tis2para is None : tis2paralist = [None] * len(args.tis2paths)
  else : tis2paralist = [s.split(';') for s in args.tis2para]

  tis1offdictlist = [find_offset(bampaths, para) for bampaths, para in zip(tis1bampathslist, tis1paralist)]
  tis2offdictlist = [find_offset(bampaths, para) for bampaths, para in zip(tis2bampathslist, tis2paralist)]

  if len(args.tis1labels) < len(args.tis1paths) : 
    for i in range(len(args.tis1labels), len(args.tis1paths)):
      args.tis1labels.append(args.tis1paths[i])
  if len(args.tis2labels) < len(args.tis2paths) :
    for i in range(len(args.tis2labels), len(args.tis2paths)):
      args.tis2labels.append(args.tis2paths[i])
  title = args.tis1labels + args.tis2labels
  tis_title = ['TIS_'+lab for lab in title]
  rna_title = ['RNA_'+lab for lab in title]
  l = len(title)

  if args.rnaseq is not None :
    if args.verbose : print("Loading RNASeq data...")
    rna_profile = exp.Profile()
    for lst in io.splitIter(args.rnaseq):
      try : values = list(map(int, lst[1:]))
      except : 
        try : 
          values = list(map(float, lst[1:]))
          print('Error: RNASeq data should be integers {}.'.format(lst))
          sys.exit(1)
        except ValueError: pass
        continue
      m = min(values)
      if m < 0 : 
        print('Error: RNASeq data should be non-negative integers {}.'.format(lst))
        sys.exit(1)
      e = exp.Exp(lst[0], values)
      rna_profile.add_exp(e)
    if args.verbose : print("{} genes.".format(len(rna_profile)))

  if args.verbose : print("Loading {} TIS data...".format(l))
  gname, gpos, gsig = {}, {}, {}
  tall = [] # {}, {}
  gid = {}
  anno = {} # annotated TIS
  for i, fname in enumerate(args.tis1paths + args.tis2paths) :
    n = 0
    tdata = {}
    for lst in io.splitIter(fname) :
      try : tis = (lst[1], int(lst[6]))
      except : continue
      gid[lst[1]] = lst[0]
      cnt, pval, qval =int(lst[10]), float(lst[11]), float(lst[args.qi])
      tdata[tis] = cnt, pval, qval
      lst[4] = get_tis(lst[4])
      if lst[8] == 'Annotated' : anno[lst[4]] = 1 # genome position
      if sig(tdata[tis]) : 
        n += 1
        #lst[4] = get_tis(lst[4])
        gname[tis] = '\t'.join(lst[:9]) # information for the TIS
        gpos[tis] = lst[4]
        if tis not in gsig : gsig[tis] = [0] * len(tis_title)
        gsig[tis][i] = 1
    if args.verbose : print("{} TISs in {}.".format(n, fname))
    tall.append(tdata)

  profile = exp.Profile()
  profile2 = exp.Profile() # uniq TISs for TMM
  trans_for_bam = {} # TIS genes need to be analyzed
  es = {}
  uniq_gpos = {}
  for tis in gname : # t1 :
    values = []
    for i, tdata in enumerate(tall) : 
      if tis not in tdata : 
        if tis[0] not in trans_for_bam : trans_for_bam[tis[0]] = [{} for j in range(l)]
        trans_for_bam[tis[0]][i][tis[1]] = None 
        values.append(0) #(None)
      else : values.append(tdata[tis][0])
    if args.rnaseq is not None : 
      if tis[0] in rna_profile.exps : values += rna_profile.exps[tis[0]].data # trans level
      elif gid[tis[0]] in rna_profile.exps : values += rna_profile.exps[gid[tis[0]]].data # gene level
      else :
        print('Warning: transcript {} {} is not found in RNA file!'.format(gid[tis[0]], tis[0]))
        values += [0] * l # len(title)

    e = exp.Exp(gname[tis], values)
    e.tis = tis
    es[tis] = e
    profile.add_exp(e)
    if gpos[tis] not in uniq_gpos :
      if args.normanno and gpos[tis] not in anno : continue
      if args.normcomm : 
        for tdata in tall :
          if tis not in tdata or not sig(tdata[tis]) : break
        else : profile2.add_exp(e)
      #elif args.normanno : 
        #if gpos[tis] in anno : profile2.add_exp(e)
      else : 
        profile2.add_exp(e)
      uniq_gpos[gpos[tis]] = 1

  if args.verbose : print("{} TISs in total.".format(len(profile)))
  if not args.normcomm :
    uns = [0] * len(tis_title)
    for i in range(len(tis_title)) :
      uns[i] = len([i for e in profile2 if gsig[e.tis][i] == 1])
    m = min(uns)
    elst = list(profile2.exps.values())
    profile3 = exp.Profile()
    for i in range(len(tis_title)) :
      for e in elst : e.value[0:2] = [gsig[e.tis][i], e.data[i]]
      elst.sort(reverse = True)
      for j in range(m) :
        if elst[j].id not in profile3.exps : profile3.add_exp(elst[j])
    profile2 = profile3

  if args.verbose : print ("Reading bams...")
  trans_iter = io.transIter(args.genepath, fileType = args.geneformat, verbose = args.verbose, filt = trans_for_bam)
  para_iter = transPara(trans_iter, trans_for_bam)
  if args.numProc <= 1 : pred_iter = imap(_get_tis, para_iter)
  else : 
    from multiprocessing import Pool
    pool = Pool(processes = args.numProc - 1)
    pred_iter = pool.imap_unordered(_get_tis, para_iter, chunksize = 5)
  for result in pred_iter: 
    tid, pos_cnt = result # r1, r2 = result
    for i, pc in enumerate(pos_cnt):
      for pos in pc :
        tis = tid, pos
        es[tis].data[i] = pc[pos]

  if len(args.tis1paths) > 1 or len(args.tis2paths) > 1 or args.export is not None : 
    if args.export is None : args.export = 'tisdiff_export.txt'
    if args.verbose : print('Export TIS counts table to {}.'.format(args.export))
    exfile = open(args.export, 'w')
    if args.rnaseq is None : exfile.write(io.tabjoin('TIS', tis_title)+'\n') # args.tis1labels, args.tis2labels)+'\n')
    else :  exfile.write(io.tabjoin('TIS', tis_title, rna_title)+'\n')
    for tis in gname:
      s = '{}_{}_{}\t'.format(tis[0], tis[1], gpos[tis])
      s += io.tabjoin(es[tis].data)
      exfile.write(s+'\n')
    return
  if args.scalefactor is not None : scale = args.scalefactor
  else :
    if args.verbose : print ('Estimate scale factor...')
    f = profile2.TMM(i1 = 0, i2 = 1)
    if args.verbose : print ('TIS TMM log2 f = {}'.format(f))
    scale = 2 ** (-f)
  if args.rnaseq is not None : 
    if args.rnascale is not None : scale_r = args.rnascale
    else :
      fr = rna_profile.TMM(i1 = 0, i2 = 1) # for only one replicate
      if args.verbose : print ('RNASeq TMM log2 f = {}'.format(fr))
      scale_r = 2 ** (-fr)

  if args.verbose : print ('Diff test...')

  exps = profile.exps.values()
  for e in exps:
    if args.rnaseq is None :
      p = 1 / (scale + 1)
      x, y = e.data[0], e.data[1] # [2]
      n = x + y
      if x == 0 : fc, alt = 'INF', 'less'
      elif y == 0 : fc, alt = 0, 'greater'
      else : 
        fc = y / (1.0 * x * scale) # / y
        if scale * x <= y : alt = 'less'
        else : alt = 'greater'
      pv = stat.binom_test(n, x, p = p, alt = alt)
    else : 
      x, y, r1, r2 = e.data[0:4]
      if x == 0 : fc, alt = 'INF', 'less'
      elif y == 0 or r1 == 0: fc, alt = 0, 'greater'
      elif r2 == 0 : fc, alt = 'INF', 'less'
      else : 
        fc = y / (1.0 * x * scale) / (r2 / (1.0 * r1 * scale_r))
        if fc >= 1 : alt = 'less' # test x
        else : alt = 'greater'
      if args.chi2 : pv = ribo.TIStest_chi2(x, y, r1, r2, scale, scale_r, alt = alt)
      elif args.betabinom : pv = ribo.TIStest_betaBinom(x, y, r1, r2, scale, scale_r, alt = alt)
      else : pv = ribo.TIStest_FisherExact(x, y, r1, r2, scale, scale_r, alt = alt)
    pv *= 2 # two tailed
    if pv > 1 : pv = 1
    e.data.append(fc)
    e.data.append(pv) 

  result = profile.BHcorrection(-1, append = True) # (5) 

  if args.verbose : print ('Output...')
  outfile = open(args.output, 'w')
  s = "Gid\tTid\tSymbol\tGeneType\tGenomePos\tStartCodon\tStart\tStop\tTisType\t"
  s += '\t'.join(tis_title)
  if args.rnaseq is not None : s += '\t' + '\t'.join(rna_title)
  s += '\tFoldChange\tDiffPvalue\tDiffQvalue\n'
  outfile.write(s)
  for e in profile :
    fc = e.data[-3]
    e.is_q = e.is_fc = True
    if fc != 'INF' and fc != 0 and max(fc, 1/fc) < args.foldchange : 
      e.is_fc = False
    if e.data[-2] > args.opth or e.data[-1] > args.oqth : 
      e.is_q = False
    if e.is_q and e.is_fc : 
      outfile.write(str(e)+'\n')

  # Plot
  if args.plotout is not None :
    if args.verbose : print ("Ploting...")
    from zbio import plot
    plot.figure(figsize = args.figsize)
    if args.rnaseq is not None : 
      qd1 = [math.log(e.data[0]+1,2) - math.log(e.data[2]+1,2) for e in exps if e.is_q and e.is_fc]
      qd2 = [math.log(e.data[1]+1,2) - math.log(e.data[3]+1,2) for e in exps if e.is_q and e.is_fc]
      pd1 = [math.log(e.data[0]+1,2) - math.log(e.data[2]+1,2) for e in exps if e.is_q and not e.is_fc]
      pd2 = [math.log(e.data[1]+1,2) - math.log(e.data[3]+1,2) for e in exps if e.is_q and not e.is_fc]
      nd1 = [math.log(e.data[0]+1,2) - math.log(e.data[2]+1,2) for e in exps if not e.is_q]
      nd2 = [math.log(e.data[1]+1,2) - math.log(e.data[3]+1,2) for e in exps if not e.is_q]
    else :
      qd1 = [math.log(e.data[0]+1,2) for e in exps if e.is_q and e.is_fc]
      qd2 = [math.log(e.data[1]+1,2) for e in exps if e.is_q and e.is_fc]
      pd1 = [math.log(e.data[0]+1,2) for e in exps if e.is_q and not e.is_fc]
      pd2 = [math.log(e.data[1]+1,2) for e in exps if e.is_q and not e.is_fc]
      nd1 = [math.log(e.data[0]+1,2) for e in exps if not e.is_q]
      nd2 = [math.log(e.data[1]+1,2) for e in exps if not e.is_q]
    plot.scatter(qd1, qd2, alpha=0.1, edgecolors='none', color='r', label='q < {} & FC > {}'.format(args.oqth, args.foldchange))
    plot.scatter(pd1, pd2, alpha=0.1, edgecolors='none', color='y', label='q < {} & FC <= {}'.format(args.oqth, args.foldchange))
    plot.scatter(nd1, nd2, alpha=0.1, edgecolors='none', color='g', label='q >= {}'.format(args.oqth))
    plot.legend(loc='upper left', frameon=False)
    plot.xlabel(title[0])
    plot.ylabel(title[1])
    if args.rnaseq is not None : d = (fr - f) / 2 
    else : d = - f / 2
    m1 = max(min(qd1+pd1+nd1), min(qd2+pd2+nd2))
    m2 = min(max(qd1+pd1+nd1), max(qd2+pd2+nd2))
    plot.plot([m1-d,m2-d],[m1+d,m2+d], color='k', linestyle = ':')
    d2 = d - math.log(args.foldchange, 2) / 2
    plot.plot([m1-d2,m2-d2],[m1+d2,m2+d2], color='r', linestyle = ':')
    d2 = d + math.log(args.foldchange, 2) / 2
    plot.plot([m1-d2,m2-d2],[m1+d2,m2+d2], color='r', linestyle = ':')
    plot.save(args.plotout)

    if args.plotma is not None :
      exps = profile2.exps.values()
      plot.figure(figsize = args.figsize)
      plot.axhline(f)
      ms = [e.M for e in exps if hasattr(e, 'select') and e.select]
      aa = [e.A for e in exps if hasattr(e, 'select') and e.select]
      plot.scatter(aa, ms, alpha=0.1, edgecolors='none',color='r')
      ms = [e.M for e in exps if hasattr(e, 'select') and not e.select]
      aa = [e.A for e in exps if hasattr(e, 'select') and not e.select]
      plot.scatter(aa, ms, alpha=0.1, edgecolors='none',color='b')
      plot.save(args.plotma)


    
def transPara(trans_iter, trans_for_bam) : #g1, g2):
  '''Generate parameters (trans, g1t, g2t) for function _get_tis()
  '''
  for t in trans_iter :
    if t.id in trans_for_bam : yield t, trans_for_bam[t.id]

def _get_tis(ps) : 
  '''get non-significant TIS counts
  '''
  t, pos_cnt = ps # g1t, g2t = ps
  bamlist = tis1bampathslist + tis2bampathslist
  offlist = tis1offdictlist + tis2offdictlist
  for i, pc in enumerate(pos_cnt) :
    if len(pc) == 0 : continue
    ttis = ribo.multiRibo(t, bamlist[i], offdict=offlist[i], compatible=compatible, mis=compatiblemis, paired=paired)
    for pos in pc : pc[pos] = ttis.cnts[pos]
  return t.id, pos_cnt

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
    
