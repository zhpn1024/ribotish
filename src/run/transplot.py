from ribotish.zbio import ribo, plot, bam, io, orf, tools, fa
from os.path import isfile

def help():
  return "Plot riboseq profile for a given transcript"
def set_parser(parser):
  # basic input options
  parser.add_argument("-g", type=str, dest="genepath", required=True, help='Gene annotation file')
  parser.add_argument("-t", type=str, dest="tid", required=True, help='Transcript id')
  parser.add_argument("-b", type=strlist, dest="bampaths", required=True, help="List of riboseq bam files, comma seperated")
  parser.add_argument("-f", type=str, dest="genomefapath", required=True, help="Genome fasta file")
  parser.add_argument("-o", type=str, dest="output", help="Output pdf figure file")
  # alt input options
  parser.add_argument("-l", type=strlist, dest="labels", help="Labels for riboseq bam files, comma seperated")
  parser.add_argument("--para", type=strlist, help="Input offset parameter files for -b bam files")
  parser.add_argument("--scale", type=strlist, help="Input scale parameters for bam files, comma seperated")
  parser.add_argument("--rna", type=strlist, help="Which bam files are RNASeq instead of RiboSeq, 0 based, comma seperated")
  parser.add_argument("--rnaoffset", type=int, help="Offset value for RNASeq reads (default: 12)")
  parser.add_argument("--rnacol", type=strlist, help="Color for RNASeq tracks, comma seperated, matching --rna option (default: black)")
  parser.add_argument("--ymax", type=strlist, help="Max y scale for tracks, comma seperated (default: auto)")

  parser.add_argument("-s", type=int2, dest="size", help="Figure size (default: 12,auto)")
  parser.add_argument("-r", type=int2, dest="range", help="Range shown on the transcript, format: start,stop (default: full transcript)")
  parser.add_argument("--geneformat", type=str, default='auto', help="Gene annotation file format (gtf, bed, gpd, gff, default: auto)")
  parser.add_argument("--alt", action="store_true", help="Use alternative start codons (all codons with 1 base different from ATG)")
  parser.add_argument("--altcodons", type=strlist, help="Use provided alternative start codons, comma seperated, eg. CTG,GTG,ACG")
  parser.add_argument("--morecds", type=str, help="More cds regions plot in ORF track, format: start1-stop1,start2-stop2...")
  parser.add_argument("--morecdsgp", type=str, help="More cds regions in genome position, format: start1-stop1,start2-stop2...")
  parser.add_argument("--morecdsbox", action="store_true", help="Add box to more CDS regions")
  parser.add_argument("--morecdslabel", type=str, help="Label for more cds regions plot in ORF track, comma seperated")
  parser.add_argument("--markpept", type=str, help="Mark peptide position in morecds, relative to morecds start format: start1-stop1,start2-stop2...")

  parser.add_argument("--colorblind", action="store_true", help="Use a color style readable for color blind people ('#F00078,#00F000,#0078F0')")
  parser.add_argument("--colors", type=strlist, help="User specified Matplotlib accepted color codes for three frames (default: 'r,g,b')")

  parser.add_argument("--paired", action="store_true", help="Reads are paired end")
  parser.add_argument("--nocompatible", action="store_true", help="Do not require reads compatible with transcript splice junctions")
  parser.add_argument("--compatiblemis", type=int, default=2, help="Missed bases allowed in reads compatibility check (default: 2)")
    
def int2(s):
  '''Convert comma seperated string to tuple
  '''
  lst = eval('['+s+']')
  return tuple(map(int, lst))

def strlist(s):
  '''Convert comma seperated file name string to list
  '''
  return s.split(',')

def rnaoffset(r, off):
  return off

def run(args):
  '''Main function for trans plot
  '''
  # load data
  t = io.transFetch(args.genepath, tid = args.tid, fileType = args.geneformat)
  l = len(args.bampaths)
  if args.para is None: args.para = [None] * l
  elif len(args.para) < l: args.para += [None] * (l - len(args.para))
  #riboffdict = find_offset(args.bampaths, args.para)
  global rnaoff
  rnaoff = 12
  if args.rnaoffset is not None:
    rnaoff = args.rnaoffset
  ribos = []
  if args.rna is None: args.rna = []
  else: 
    args.rna = list(map(int, args.rna))
    rnacol = {}
    if args.rnacol is not None:
      for i, c in enumerate(args.rnacol):
        rnacol[args.rna[i]] = c
  compatible = not args.nocompatible
  ymax = [None] * l
  if args.ymax is not None:
    for i, ym in enumerate(args.ymax):
      try: ymax[i] = float(ym)
      except: pass
  morecds = None
  morecdslabel = None
  if args.morecds is not None or args.morecdsgp is not None:
    morecds = []
    if args.morecds is not None:
      for s in args.morecds.split(','):
        morecds.append(list(map(int, s.split('-'))))
    if args.morecdsgp is not None:
      for s in args.morecdsgp.split(','):
        gp = list(map(int, s.split('-')))
        cp = [t.cdna_pos(gp[0]), t.cdna_pos(gp[1])]
        if None in cp: continue
        cp.sort()
        morecds.append(cp)
    if args.morecdslabel is not None:
      morecdslabel = args.morecdslabel.split(',')
  markpept = None
  if morecds is not None and args.markpept is not None:
    markpept = []
    for i, s in enumerate(args.markpept.split(',')):
      mp0, mp1 = list(map(int, s.split('-')))
      markpept.append([morecds[i][0]+mp0, morecds[i][0]+mp1])
  for i, s in enumerate(args.bampaths):
    bampathslist = s.split(';')
    if i not in args.rna: 
      paralist = args.para[i]
      if paralist is not None: 
        paralist = paralist.split(';')
      offdictlist = find_offset(bampathslist, paralist)
    #ribobamfile = bam.Bamfile(ribopath, "rb")
    if i in args.rna:
      tribo = ribo.multiRibo(t, bampathslist, offset=rnaoffset, offdict = [rnaoff]*len(bampathslist), compatible = compatible, mis = args.compatiblemis, paired = args.paired)
      #tribo = ribo.Ribo(t, ribobamfile, offset=rnaoffset, offdict = rnaoff)
    else:
      tribo = ribo.multiRibo(t, bampathslist, offdict = offdictlist, compatible = compatible, mis = args.compatiblemis, paired = args.paired)
      #tribo = ribo.Ribo(t, ribobamfile, offdict = riboffdict[i])#riboffdict
    ribos.append(tribo)
  # prepare parameters
  scales = [1] * l
  if args.scale is not None:
    for i, s in enumerate(args.scale):
      try:
        s = float(s)
        scales[i] = s
      except: pass
  tl = t.cdna_length()
  rst, ren = 0, tl # start, end
  if args.range is not None :
    if rst < args.range[0] : rst = args.range[0]
    if len(args.range) < 2: ren = tl
    else: ren = min(args.range[1], tl)
  rlen = ren - rst #rst, rlen = 0, 400#tl
  cds1, cds2 = t.cds_start(True), t.cds_stop(True)
  nsub = l + 1
  if args.size is None : size = 12, 2 * nsub + 1
  else : size = args.size
  anno = t.symbol + ' ' + t.id + ' ' + t.genetype
  rids = args.bampaths
  if args.labels is not None : 
    for i, lab in enumerate(args.labels):
      if lab is not None and len(lab) > 0 : rids[i] = lab
  if args.altcodons is not None : 
    args.alt = True
    orf.cstartlike = [c.upper() for c in args.altcodons]
  col = ['r', 'g', 'b']
  if args.colorblind:
    col = ['#F00078', '#00F000', '#0078F0'] # ['#FA007D', '#00E800', '#007DFA']
  if args.colors is not None:
    col = args.colors
  blk = ['k', 'k', 'k']
  # plot reads
  fp, axarr = plot.subplots(nsub, sharex=True, figsize= size)
  showlegend = True
  for i in range(len(ribos)):
    if i == 0 : title = anno+'\n'+rids[i]
    else : title = rids[i]
    if i in args.rna:
      rcol = blk
      if i in rnacol: rcol = [rnacol[i]] * 3
      plot.riboShow(axarr[i], t, ribos[i].cnts, start = rst, stop = ren, title = title, col=rcol, scale=scales[i], showlegend = False, showframe=False, ymax=ymax[i])
    else:
      plot.riboShow(axarr[i], t, ribos[i].cnts, start = rst, stop = ren, title = title, col=col, scale=scales[i], showlegend=showlegend, ymax=ymax[i])
      showlegend = False
  # plot ORF
  genome = fa.Fa(args.genomefapath)
  tsq = tools.gtf2seq(genome, t)
  orfs = orf.orflist(tsq, tail=tl)
  plot.orfShow(axarr[-1], orfs, start = rst, stop = ren, cds = [cds1,cds2], alt = args.alt, col=col, morecds=morecds, morecdsbox=args.morecdsbox, morecdslabel=morecdslabel, markpept=markpept)
  # plot splice site
  x, l = [], 0
  for e in t.exons[:-1]:
    l += len(e)
    li = l - rst
    if 0 <= li <= rlen : x.append(li)
  axarr[-1].plot(x,[0.12]*len(x),'m^', alpha = 1)#,color='m')
  axarr[-1].spines['top'].set_visible(False)
  axarr[-1].spines['left'].set_visible(False)
  axarr[-1].spines['right'].set_visible(False)
  axarr[-1].get_yaxis().set_visible(False)
  axarr[-1].set_xticklabels(list(map(int,axarr[0].get_xticks() + rst)))
  # save figure
  output = args.output
  if output is None : output = '%s.pdf' % t.symbol
  plot.savefig(output)

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

if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())

