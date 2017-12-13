'''
Some plot functions
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from matplotlib.pylab import *

def plotTrans(t, ypos = 0, intv = None, r = [0.1, 0.3], color = 'blue',rid = -0.5):
  '''plot transcript
  '''
  plot([t.start,t.stop],[ypos,ypos],color=color)
  # arrows
  x, y = [], []
  if intv is None : intv = len(t) / 20
  for i in range(t.start+intv//2, t.stop-intv//3, intv):
    x.append(i)
    y.append(ypos)
  arr = '>'
  if t.is_reverse() : arr = "<"
  plot(x,y,'w'+arr)
  # blocks
  x = [[],[]]
  y = [[],[]]
  orf = t(start = t.thick_start, stop = t.thick_stop)
  for e in t.exons:
    for es in e - orf : # UTR 
      x[0].append(es.start)
      y[0].append(len(es))
    for ei in e.intersect(orf): # CDS
      x[1].append(ei.start)
      y[1].append(len(ei))
  bar(x[0],[r[0]*2]*len(x[0]),width=y[0],bottom=ypos-r[0],edgecolor=color,color=color)
  bar(x[1],[r[1]*2]*len(x[1]),width=y[1],bottom=ypos-r[1],edgecolor=color,color=color)
  text((t.start+t.stop)//2, ypos+rid, t.id)

def save(filename):
  savefig(filename, transparent=True)

def riboShow(ax, trans, cnts, start = 0, stop = -1, ymax = None, scale = 1, col = ['r','g','b'], title = '', showlegend = False, showframe = True, bottom = 0.8, height = 0.1):
  '''plot riboseq profile
  '''
  if stop < start : stop = trans.cdna_length()
  rlen = stop - start
  lx = [[], [], []]
  ly = [[], [], []]
  m = 1
  for j in range(rlen):
    p = j + start
    if cnts[j + start] > 0 : 
      i = p % 3
      lx[i].append(j)
      y = cnts[p] * scale
      if ymax is not None and ymax > 0 and y > ymax: y = ymax
      ly[i].append(y)
      if m < y : m = y
  ylab = 'Count'
  if scale != 1 : ylab = 'Scaled count'
  if ymax is None or ymax < 0 : ymax = m
  for i in range(3):
    ax.bar(lx[i], ly[i], color=col[i], width=1, edgecolor=col[i], log=False, alpha=0.4, label='Frame '+str(i+1))
  
  [ax.spines[side].set_visible(False) for side in ('right','top','bottom')]
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  ax.set_xlim((0, rlen))
  ax.set_ylim((0, ymax))
  ax.set_ylabel(ylab)
  ax.set_title(title)
  if showlegend :
    try : ax.legend(loc='best', frameon=False)
    except : pass
  if not showframe : return
  fx = [[],[],[]]
  fy = [[],[],[]]
  fw = [[],[],[]]
  for i in range(3):
    last = False
    for j in range(i, stop, 3) :
      if j < start : continue
      top = True
      if cnts[j] <= 0 : top = False
      if j-1 >= 0 and cnts[j] <= cnts[j-1] : top = False
      if j+1 < stop and cnts[j] <= cnts[j+1] : top = False
      if top : ##
        if last : tw += 3
        else : tx, tw = j-start, 1
        last = True
      else : 
        if last and tw > 3 :
          fx[i].append(tx)
          fy[i].append(ymax * height)
          fw[i].append(tw)
        last = False
    if last and tw > 3 :
      fx[i].append(tx)
      fy[i].append(ymax * 0.1)
      fw[i].append(tw)
  for i in range(3):
    ax.bar(fx[i], fy[i], color=col[i], bottom = ymax * bottom,width=fw[i], alpha=0.2, linewidth = 0)
  
def orfShow(ax, orfs, start = 0, stop = -1, col = ['r','g','b'], cds = [None, None], title = 'Potential ORFs in 3 reading frames', alt = True):
  '''plot possible ORFs
  '''
  if stop < start : stop = trans.cdna_length()
  rlen = stop - start
  # ORF in 3 frames
  lx = [[],[],[]]
  ly = [[],[],[]]
  orf_s = []
  for o in orfs:
    if not alt and len(o.starts) == 0 : continue
    if 0 <= o.stop <= start or o.start(alt=alt) >= stop : continue
    orf_s.append(o)
    o1 = o.start(alt=alt) - start
    if o1 < 0: o1 = 0
    lx[o.frame-1].append(o1) # (o.start(alt=alt) - start)
    if o.has_stop():
      o2 = min(stop, o.stop) - max(o.start(alt=alt), start)
    else:
      o2 = stop - max(o.start(), start)
    ly[o.frame-1].append(o2)
    #if o.has_stop(): ly[o.frame-1].append(o.stop-o.start(alt=alt))
    #else : ly[o.frame-1].append(rlen)
  for i in range(3):
    ax.bar(lx[i], [0.2]*len(lx[i]), color=col[i], bottom=2-i+0.4, width=ly[i], alpha=0.4, linewidth=0)
  # annotated ORF
  if cds[0] is not None and not (cds[0]>stop or cds[1]< start): 
    i = cds[0] % 3
    newcds = [c - start for c in cds]
    if newcds[0] < 0: newcds[0] = 0
    if cds[1] > stop: newcds[1] = rlen
    ax.text(max(newcds[0],0), 2-i+0.8, 'Annotated ORF', color=col[i])
    ax.bar(newcds[0], [0.4] ,color=col[i], bottom=2-i+0.3, width=newcds[1]-newcds[0], alpha=0.3, edgecolor=col[i], linewidth=2)
  # start & stop codons
  lx = [[],[],[]]
  ly = [[],[],[]]
  lz = [[],[],[]]
  for o in orf_s:
    lx[o.frame-1] += [s - start for s in o.starts if start<=s<stop]
    if alt : ly[o.frame-1] += [s - start for s in o.altstarts if start<=s<stop] # orf.altstarts
    if o.has_stop() and start<=o.stop-3<stop: lz[o.frame-1] += [o.stop - 3 - start]
  for i in range(3):
    ax.bar(ly[i], [0.4]*len(ly[i]), color='yellow', bottom=2-i+0.3, width=3, alpha=0.6, edgecolor='yellow')
    ax.bar(lx[i], [0.4]*len(lx[i]), color='lime', bottom=2-i+0.3, width=3, alpha=1, edgecolor='lime')
    ax.bar(lx[i], [0.04]*len(lx[i]), color='w', bottom=2-i+0.48, width=3, edgecolor='w')
    ax.bar(lz[i], [0.4]*len(lz[i]), color='red', bottom=2-i+0.3, width=3, alpha=1, edgecolor='red')
    ax.bar(lz[i], [0.04]*len(lz[i]), color='k', bottom=2-i+0.48, width=3, edgecolor='k')
  ax.set_title(title)
  
