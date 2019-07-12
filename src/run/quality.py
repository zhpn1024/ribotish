#import sys, getopt
from ribotish.zbio import ribo, plot, io, bam

def help():
  return "Quality control for riboseq data"
def set_parser(parser):
  # basic input options
  parser.add_argument("-b", type=str, dest="ribobampath", required=True, help="Riboseq bam file")
  parser.add_argument("-g", type=str, dest="genepath", help='Gene annotation file')
  parser.add_argument("-o", type=str, dest="output", help="Output data file (default: ribobampath[:-4]+ '_qual.txt')")
  # alt input options
  parser.add_argument('-t', "--tis", action="store_true", help="The data is TIS enriched (for LTM & Harritonine)")
  parser.add_argument("-i", type=str, dest="input", help="Input previous output file, do not read gene file and bam file again")
  parser.add_argument("--geneformat", type=str, default='auto', help="Gene annotation file format (gtf, bed, gpd, gff, default: auto)")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if annotation chr ids are not same as chr ids in bam/fasta files")

  # quality result output
  parser.add_argument("-f", type=str, dest="figpdfpath", help="Output pdf figure file (default: ribobampath[:-4]+ '_qual.pdf')")
  parser.add_argument("-r", type=str, dest="parapath", help="Output offset parameter file (default: ribobampath[:-4]+ '.para.py')")
  # other options
  parser.add_argument("-l", type=int2, dest="lens", default=(25,35), help="Range of tag length (default: 25,35)")
  parser.add_argument("-d", type=int2, dest="dis", default=(-40,20), help="Position range near start codon or stop codon (default: -40,20)")
  parser.add_argument("--bins", type=int, default=20, help="Bins for cds profile (default: 20)")
  parser.add_argument("--nom0", action="store_true", help="Do not consider reads with mismatch at position 0 as a new group")
  parser.add_argument("--th", type=float, default=0.5, help="Threshold for quality (default: 0.5)")
  # Reads filters
  parser.add_argument("--maxNH", type=int, default=1, help="Max NH value allowed for bam alignments (default: 1)")
  parser.add_argument("--minMapQ", type=float, default=1, help="Min MapQ value required for bam alignments (default: 1)")
  parser.add_argument("--secondary", action="store_true", help="Use bam secondary alignments")
  parser.add_argument("--paired", action="store_true", help="Reads are paired end")
  
  parser.add_argument("--colorblind", action="store_true", help="Use a color style readable for color blind people ('#F00078,#00F000,#0078F0')")
  parser.add_argument("--colors", type=strlist, help="User specified Matplotlib accepted color codes for three frames (default: 'r,g,b')")

  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes (default: 1)")
  parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase output verbosity")

def int2(s):
  '''Convert comma seperated string to tuple
  '''
  lst = eval('['+s+']')
  return tuple(map(int, lst))

def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')

#m0 = True
def run(args):
  '''Main function for quality control
  '''
  global m0
  m0 = True
  if args.nom0 : m0 = False
  ribo.maxNH, ribo.minMapQ, ribo.secondary = args.maxNH, args.minMapQ, args.secondary
  global fbcols
  fbcols = ['r','g','b','r','g','b']
  if args.colorblind:
    fbcols = ['#F00078', '#00F000', '#0078F0'] # ['#FA007D', '#00E800', '#007DFA']
  if args.colors is not None:
    fbcols = args.colors

  if args.input is None :
    minR = 1
    if args.genepath is None : 
      print('Error: missing -g gene annotation file!')
      exit(1)
    if args.ribobampath is None : 
      print('Error: missing -b riboseq bam file!')
      exit(1)
    if args.output is None : 
      args.output = args.ribobampath[:-4] + '_qual.txt'
    if args.chrmap is not None :
      chrmap = {}
      for lst in io.splitIter(args.chrmap, sep=None):
        if len(lst) < 2: continue
        chrmap[lst[0]] = lst[1]
        chrmap[lst[1]] = lst[0]
      bam.chrmap = chrmap
    #fa.chrmap = chrmap

    # read data
    results = ribo.lendis(args.genepath, args.ribobampath, lens = args.lens, dis = args.dis, minR = minR, m0 = m0, paired = args.paired,
                          cdsBins = args.bins, numProc = args.numProc, verbose = args.verbose, geneformat = args.geneformat)
    # write quality data
    outfile = open(args.output, 'w')
    s = 0
    for ad in results :
      s += ad.size()
      ad.write(outfile)
    print('Counted reads: {}'.format(s))
    outfile.close()
    if s == 0:
      print('Error: no reads found! Check read length or protein coding annotation.')
      exit(1)
    args.input = args.output
  if args.figpdfpath is None : 
    args.figpdfpath = args.ribobampath[:-4] + '_qual.pdf'
  if args.parapath is None : 
    args.parapath = args.ribobampath + '.para.py'
  # get quality plot and quality parameter file
  qualityPlot(args)

#fbcols = ['r','g','b','r','g','b']
#if args.colorblind:
  #fbcols = ['#FF008C', '#00FF00', '008CFF']

def qualityPlot(args):
  ''' Quality plot and offset parameters
  '''
  global m0
  if args.figpdfpath is None : args.figpdfpath = args.ribobampath[:-4] + '_qual.pdf'
    #print('Need to provide -f output pdf figure file name!')
    #exit(1)
  elif args.verbose : print('Quality pdf figure will be saved to {}'.format(args.figpdfpath))
  if args.parapath is None : args.parapath = args.ribobampath + '.para.py'
    #print('Need to provide -r output offset parameter file name!')
    #exit(1)
  elif args.verbose : print('Offset parameter will be saved to {}'.format(args.parapath))
  offdict = {}
  infile = open(args.input,'r')

  lines = [l.rstrip() for l in infile]
  lendis = eval(lines[0])
  dis1 = eval(lines[1])
  dis2 = eval(lines[2])
  disf = eval(lines[3])
  disc = eval(lines[4])
  if m0 and len(lines) >= 10 : 
    lendism0 = eval(lines[5])
    dis1m0 = eval(lines[6])
    dis2m0 = eval(lines[7])
    disfm0 = eval(lines[8])
    discm0 = eval(lines[9])
    m0 = True
    offdict['m0'] = {}
  else : m0 = False

  fig = plot.figure(figsize = (16, 16))
  k = list(lendis.keys())
  k.sort()
  #x = list(range(*args.dis))
  if m0 : gs0 = plot.matplotlib.gridspec.GridSpec(1, 2, top = 0.92, bottom = 0.86)
  else : gs0 = plot.matplotlib.gridspec.GridSpec(1, 1, top = 0.92, bottom = 0.86)

  ax = plot.subplot(gs0[0, 0])
  x1 = [l - 0.4 for l in lendis.keys()]
  ax.bar(x1, lendis.values(), alpha=0.6, align='edge')
  formatax(ax)
  ty = nearest(max(max(lendis.values()),1), up=False)
  ax.set_yticks([0, ty])
  ax.set_yticklabels([0,numk(ty)])
  if m0 : ax.set_title("5' end match RPFs\n\nRPF length distribution")
  else : ax.set_title("RPFs\n\nRPF length distribution")
  if m0 : 
    ax = plot.subplot(gs0[0, 1])
    x1 = [l - 0.4 for l in lendism0.keys()]
    ax.bar(x1, lendism0.values(), alpha=0.6, align='edge')
    formatax(ax)
    ty = nearest(max(max(lendism0.values()),1), up=False)
    ax.set_yticks([0, ty])
    ax.set_yticklabels([0,numk(ty)])
    ax.set_title("5' end mismatch RPFs\n\nRPF length distribution")

  if m0 : gs = plot.matplotlib.gridspec.GridSpec(len(k), 12, top = 0.80, bottom = 0.05, left = 0.07, right = 0.98, wspace=0.7, hspace=0.3)
  else : gs = plot.matplotlib.gridspec.GridSpec(len(k), 6, top = 0.80, bottom = 0.05, left = 0.1, right = 0.95, wspace=0.35, hspace=0.3)

  for i, l in enumerate(k):
    plot4(gs, i, l, disf, dis1, dis2, disc, offdict, start = 0, args=args)
    if m0 : plot4(gs, i, l, disfm0, dis1m0, dis2m0, discm0, offdict['m0'], start = 6, args=args)
  plot.save(args.figpdfpath) 

  parafile = open(args.parapath, 'w') 
  ribo.write_off_para(parafile, offdict)
  total = sum(lendis.values())
  if m0: total += sum(lendism0.values())
  n = nt = 0
  for l in offdict :
    if l in disf : 
      n += sum(disf[l])
      if args.tis : nt += dis1[l][-args.dis[0] - offdict[l]]
  if m0 :
    for l in offdict['m0'] :
      n += sum(disfm0[l])
      if args.tis : nt += dis1m0[l][-args.dis[0] - offdict['m0'][l]]
  if args.verbose : 
    print('Total RPF counts in mRNA: {}'.format(total))
    print('Effective RPF counts: {}'.format(n))
    if args.tis : print('Effective TIS counts: {}'.format(nt))

def nearest(n, up = True, step = 10, levels = range(1, 10)):
  lasti = i = e = 1
  while n >= i :
    for l in levels:
      i = e * l
      if n <= i : break
      lasti = i
    if n <= i : break
    e *= step
  if up : return i
  else : return lasti
def numk(n):
  step = 1000
  s = ['k', 'M', 'G', 'T']
  i = -1
  while n >= step :
    n /= step
    i += 1
    if i >= len(s) - 1 : break
  if i >= 0 : return str(n) + s[i]
  else : return str(n)
def frange(start=0,stop=1,step=1):
  lst = []
  i = start
  while i < stop:
    lst.append(i)
    i += step
  return lst
def stdisplot(ax, x, y, lab = '', hali = 'left', vali = 'top', frame = 3, f0 = 0, color = ['r','g','b'], size='medium', tcol = 'k'):
  for i in range(frame):
    i0 = (f0 + i) % frame
    ax.bar(x[i0::frame], y[i0::frame], width = 0.6, color = color[i], edgecolor = color[i],alpha=0.6, align='edge')
  if hali == 'left' : tx = min(x)
  else : tx = max(x)
  ty = nearest(max(y+[1]), up=False) ##
  ax.text(tx ,max(y+[1])*0.7 , lab,horizontalalignment=hali,verticalalignment=vali, size=size, color = tcol)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(False)
  ax.tick_params(direction='in', top = False, right = False)
  ax.set_yticks([0, ty])
  ax.set_yticklabels([0,numk(ty)])
def cdsplot(ax, disc):
  m = 1
  for i in range(3):
    ax.plot(disc[i], color = fbcols[i])
    mi = max(disc[i])
    if m < mi : m = mi
  formatax(ax)
  ty = nearest(m, up=False)
  ax.set_yticks([0, ty])
  ax.set_yticklabels([0,numk(ty)])
  ax.set_xticks([])
def formatax(ax):
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(False)
  ax.tick_params(direction='in', top = False, right = False)


def plot4(gs, i, l, disf, dis1, dis2, disc, offdict, args, start = 0):
  y0 = disf[l]
  ys = sum(y0)
  if ys > 0 : y = [float(n)/ys for n in y0]
  else : y = y0
  ax0 = plot.subplot(gs[i, start])
  ax0.bar(frange(0.2, len(y)), y, width = 0.6, color = fbcols, edgecolor = fbcols, alpha=0.6, align='edge')
  use, frame, txt = ribo.quality(y, threshold = args.th)
  ax0.text(3.1 ,0.9 ,txt, horizontalalignment= 'right',verticalalignment='top', color='k')
  formatax(ax0)
  #ax0.set_xticks(range(4))
  ax0.set_xticks([])
  ax0.set_yticks(frange(stop=1.1, step=0.5))
  if start == 0 : 
    ax0.set_ylabel(str(l)+'nt',rotation=0,size='x-large')
    if m0 : ax0.yaxis.set_label_coords(-0.9, 0.4)
    else : ax0.yaxis.set_label_coords(-0.4, 0.4)

  if args.tis : 
    use, tisframe, tistxt, mp = ribo.TISquality(dis1[l], dis = args.dis, threshold = args.th)
    if tisframe != frame : use = False
    #print use, frame, tistxt
  txt = ''
  if i == 0 : txt = "offset:"
  if use : 
    offset, th = ribo.get_offset(dis1[l], dis = args.dis, frame = frame, tis = args.tis)
    if offset is not None :
      txt += str(offset)
      offdict[l] = offset
    #else : txt = ''
  else : txt += "NA"
  tcol = 'k'
  if use : tcol = fbcols[frame]
  ax1 = plot.subplot(gs[i, start+1:start+3])
  x = list(range(*args.dis))
  stdisplot(ax1, x, dis1[l], lab = txt, hali = 'left', f0 = -args.dis[0], color=fbcols, tcol = tcol)
  
  if use and not args.tis : 
    plot.hlines(th, -18, -6, color=tcol, linestyles ='dashed')
  if args.tis : 
    #mp = ribo.getTIS(dis1[l], [d1,d2])
    ax1.text(mp+1 ,max(dis1[l]+[1])*0.7 , tistxt, horizontalalignment='left',verticalalignment='top', size='medium', color = tcol)
  ax2 = plot.subplot(gs[i, start+3:start+5])
  
  stdisplot(ax2, x, dis2[l], lab = '', hali = 'right', f0 = -args.dis[0], color=fbcols, tcol = tcol)
  ax3 = plot.subplot(gs[i, start+5])
  cdsplot(ax3, disc[l])
  if args.tis : 
    bins = len(disc[l][tisframe])
    mean = sum(disc[l][tisframe]) / bins
    if mean == 0 : enrich = 0
    else : enrich = dis1[l][mp - args.dis[0]] / mean
    if i == 0 : enrtxt = 'TIS\nenrich\n%.1f' % enrich
    else : enrtxt = '%.1f' % enrich
    ax3.text(bins ,ax3.get_ylim()[1] , enrtxt, horizontalalignment='right',verticalalignment='top', size='medium', color = tcol)
  if i == 0 : ax0.set_title('Frame distr.\n')
  if i == 0 : ax1.set_title('RPF near start codon\n')
  if i == 0 : ax2.set_title('RPF near stop codon\n')
  if i == 0 : ax3.set_title('CDS Profile\n')
  if i == gs._nrows - 1 : ax0.set_xlabel('Frame in codon', size='medium')
  if i == gs._nrows - 1 : ax1.set_xlabel('Dist. from start codon (nt)', size='medium')
  if i == gs._nrows - 1 : ax2.set_xlabel('Dist. from stop codon (nt)', size='medium')
  #if i == 0 and m0 : ax1.text(0, 0, "5' end match RPFs", ha='left', va='bottom')

if __name__=="__main__":
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())
