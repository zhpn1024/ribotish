'''
Statistics tools
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

import math
from scipy.stats import nbinom, chisquare # chisqprob
from scipy.stats import chi2
from scipy.special import betaln, betainc, gammaln
#logarr = [None] # log(N)
#logsumarr = [0] # log(N!)
def logsum(n): # log(N!)
  return gammaln(n+1)
def logarr_ext(n) : #, logarr = logarr, logsumarr = logsumarr): # prepare log values
  global logarr, logsumarr
  l = len(logarr)
  if l < n + 1 : 
    logarr += [None] * (n + 1 - l)
    logsumarr += [None] * (n + 1 - l)
    for i in range(l, n + 1):
      logarr[i] = math.log(i)
      #print i, logarr[i]
      logsumarr[i] = logsumarr[i-1] + logarr[i]
  return logarr, logsumarr
def betainc1(a, b, x):
  if a == 0: return 1
  return betainc(a, b, x)
def data_count(data):
  total, cnt = 0, 0
  for k in data:
    total += k * data[k]
    cnt += data[k]
  return total, cnt
def mean_var(data):
  total, tsq, cnt = 0, 0, 0
  for k in data:
    total += k * data[k]
    cnt += data[k]
    tsq += k * k * data[k]
  mean = float(total) / cnt
  var = float(tsq) / cnt - mean ** 2
  return mean, var
def load_data(arr):
  data = {}
  for k in arr:
    if k not in data : data[k] = 0
    data[k] += 1
  return data
def fisher_method(ps):
  n, fs = 0, 0
  for p in ps:
    if p is None : continue
    if p == 0 : return 0.0, -1 ###
    fs += - 2 * math.log(p)
    n += 1
  fp = chi2.sf(fs, 2 * n) # chisqprob(fs, 2 * n)
  return fp, fs

def combination_log(n, k, show = False): 
  '''N choose K, log combination number, NATURAL LOG
  '''
  n, k = int(n), int(k)
  if n < 0 : return None
  if k > n or k < 0: return None #None is log(0)
  #logarr_ext(n)
  lpr = logsum(n) - logsum(k) - logsum(n-k) # logsumarr[n] - logsumarr[k] - logsumarr[n-k]
  return lpr

def ACprob(x, y, r = 1): 
  '''p(y|x) = r^y C(x+y,y) / (1+r)^(x+y+1) , r = N2/N1 by Audic and Claverie 
  '''
  lp = combination_log(x+y, x)
  lp += math.log(r) * y
  lp -= math.log(1 + r) * (1 + x + y)
  return math.exp(lp)
def ACtest(x, y, r = 1, alt = 'auto', double = True): # Diff expression test by Audic and Claverie 
  if alt in ('g', 'greater') : n1, n2, nr = x, y, r  # if x > y ?
  elif alt in ('l', 'less') : n1, n2, nr = y, x, 1.0/r
  elif x * r < y : n1, n2, nr = y, x, 1.0/r # n2 is smaller than n1
  else : n1, n2, nr = x, y, r
  pv = 0
  for i in range(n2 + 1):
    pv += ACprob(n1, i, nr)
    if double and pv >= 0.5 : return 1
    #print pv
  if double : pv *= 2 # Two tailed
  if pv > 1 : pv = 1
  return pv

def FCtest(x, y, r = 1, fc = 1.5, alt = 'auto', double = True): 
  '''Test whether exp diff > fold change (binom test). For data with no replicate|dispersion 
  '''
  n = x + y
  if y > 0 : fcr = 1.0 * x * r / y
  else : fcr = fc + 1 ###
  if alt in ('g', 'greater') and fcr <= fc : return 1
  elif alt in ('l', 'less') and fcr >= 1.0 / fc : return 1
  if 1.0 / fc <= fcr <= fc : return 1
  if fcr > fc : 
    p = fc / (fc + r)
    pv = binom_test(x, n, p = p, alt = "g")
  elif fcr < 1.0 / fc : 
    p = 1.0 / (1 + r * fc)
    pv = binom_test(x, n, p = p, alt = "l")
  if double : pv *= 2  # Doubling the smaller tail
  if pv > 1 : pv = 1
  return pv
    
def hypergeo_log(N, K, n, k): 
  return combination_log(K, k) + combination_log(N-K, n-k) - combination_log(N, n)

def hypergeo(N, K, n, k):
  if k < 0 or k < K+n-N: return 0
  if k > n or k > K: return 0
  lp = hypergeo_log(N, K, n, k)
  if lp is None : return 0
  return math.exp(lp)
  
def hypergeo_test(N, K, n, k, alt = 'two.tailed') :
  p = hypergeo(N, K, n, k)
  if alt in ('g', 'greater') :
    for i in range(k, min(K, n)) :
      p += hypergeo(N, K, n, i+1)
      if p > 1 : return 1
    return p
  elif alt in ('l', 'less') :
    for i in range(max(0, n + K - N), k) :
      p += hypergeo(N, K, n, i)
      if p > 1 : return 1
    return p
  p *= 2
  if p > 1 : p = 1
  return p

def binomial(n, k, p = 0.5, show=False):
  if k > n or k < 0: return 0
  if p < 0 : p = 0
  if p > 1 : p = 1
  if k * 2 > n:
    k = n - k
    p = 1 - p
  q = 1 - p
  pr = 1.0
  nk = n - k
  if k == 0 : return q ** n
  t = int(nk / k)
  if t * k < nk: t += 1
  qi = 0
  for i in range(k):
    pi = 1.0
    pi *= n - i
    pi /= k - i
    pi *= p
    for ti in range(t):
      if qi >= nk: break
      pi *= q
      qi += 1
    pr *= pi
    if show: print (pr)
  return pr
def binom_log(n, k, p = 0.5, show = False): #log probability value
  if n < 0 : return None
  if k > n or k < 0: return None #None is log(0)
  if p <= 0 :
    if k == 0 : return 0 # log(1)
    else : return None
  elif p >= 1 : 
    if k == n : return 0
    else : return None
  #q = 1 - p
  lpr = math.log(p) * k + math.log(1-p) * (n-k)
  lpr += combination_log(n, k)
  return lpr
def binom_test0(n, k, p = 0.5, alt = "g", log = True, show=False): 
  '''binomial test, no two sided yet!
  '''
  if not log : return binomTest0(n, k, p, alt, show) # if log, p are calculated with log values
  lpk = binom_log(n, k, p)
  if show : print (lpk)
  if lpk is None : 
    if alt[0] == 'g' and p >= 1: return 1
    if alt[0] != 'g' and p <= 0: return 1
    return 0
  elif lpk == 0 : return 1
  pv = math.exp(lpk) ### log reverse
  q = 1 - p
  lp = math.log(p)
  lq = math.log(q)
  if alt[0] == 'g':
    for i in range(k, n):
      lpk += lp + math.log(n - i) - lq - math.log(i + 1) #r = p * (n - i) / q / (i + 1)
      pv += math.exp(lpk)
  else:
    for i in range(k, 0, -1):
      lpk += lq + math.log(i) - lp - math.log(n - i + 1) #r = q * i / p / (n - i + 1)
      pv += math.exp(lpk)
  return pv
def binom_test(n, k, p = 0.5, alt = "g"):
  if alt[0] == 'l':
    return betainc1(n-k, k+1, 1-p)
  elif alt[0] == 'g':
    return betainc1(k, n-k+1, p)
  else:
    if k <= n * p:
      pv = betainc1(n-k, k+1, 1-p)
    else:
      pv = betainc1(k, n-k+1, p)
    return min(pv*2, 1)

def binomTest0(n, k, p = 0.5, alt = "g", show=False): # No two sided yet!
  '''binomial test no log version
  '''
  pk = binomial(n, k, p)
  if show : print (pk)
  if pk == 1 : return 1
  if pk == 0 : 
    if alt[0] == 'g' and p >= 1: return 1
    if alt[0] != 'g' and p <= 0: return 1
    return 0
  pv = pk
  q = 1 - p
  if alt[0] == 'g': 
    for i in range(k, n):
      #r = 1.0
      r = p * (n - i) / q / (i + 1)
      pk *= r
      pv += pk
  else:
    for i in range(k, 0, -1):
      r = q * i / p / (n - i + 1)
      pk *= r
      pv += pk
  return pv
class NegBinom: 
  '''Negative binomial distribution
  Number of 'failures' before 'r' 'successes' with success probability 'p'
  Note that when p is higher, the expect is lower, consistent with nbinom in scipy.stats.
  '''
  rMax = 1e8
  rMin = 1e-8
  Delta = 1e-8
  p_record = {}
  def __init__(self, r = 1.0, p = 0.5):
    self.p = p
    self.r = r
  @property
  def q(self):
    return 1 - self.p
  def expect(self):
    return self.q * self.r / self.p
  def variance(self):
    return self.expect() / self.p
  def logpmf(self, k = 0):
    return nbinom.logpmf(k, self.r, self.p)
  def pmf(self, k = 0):
    return nbinom.pmf(k, self.r, self.p)
  def cdf(self, k = 0):
    return nbinom.cdf(k, self.r, self.p)
  def pvalue(self, k = 0, record = True):
    if record : 
      key = (self.r, self.p)
      if key not in self.p_record : self.p_record[key] = {}
      elif k in self.p_record[key] : return self.p_record[key][k]
    if k == 0: p = 1
    else: p = betainc(k, self.r, self.q)
    '''
    p = 1 - nbinom.cdf(k-1, self.r, self.p)
    if p > self.Delta : 
      if record : self.p_record[key][k] = p
      return p
    p = nbinom.pmf(k, self.r, self.p)
    if p == 0 : 
      if record : self.p_record[key][k] = p
      return p
    ka = k + 1
    pa = nbinom.pmf(ka, self.r, self.p)
    p += pa
    while pa / p > self.Delta:
      #p += pa
      ka += 1
      pa = nbinom.pmf(ka, self.r, self.p)
      p += pa
    #p += pa'''
    if record : self.p_record[key][k] = p
    return p
      
  def estimate(self, data): #data dict value:counts
    '''estimate parameters from input data (dict type, value -> counts)
    '''
    total, cnt = data_count(data)
    rmax, rmin = self.rMax, self.rMin
    rmid = math.sqrt(rmax * rmin)
    while (rmax - rmin) / rmid >= self.Delta:
      #print rmax, rmin
      score = self.r_log_like_score(data, rmid)
      #print score, cnt, rmid
      if score > 0 : rmin = rmid
      elif score < 0 : rmax = rmid
      else : break
      #rmid = (rmax + rmin) / 2
      rmid = math.sqrt(rmax * rmin)
    self.r = rmid
    #self.p = total / (self.r * cnt + total)
    self.p = self.r / (self.r + 1.0*total/cnt)
    return self.r, self.p

  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, r = -1):
    '''log likelyhood score as function of r
    '''
    if r < 0 : r = self.r
    total, cnt = data_count(data)
    s1, d = 0, 0
    for i in range(max(data.keys()) + 1):
      if i in data : s1 += d * data[i]
      d += 1.0 / (r + i)
    score = s1 / cnt + math.log(r / (r + 1.0*total/cnt))
    return score
  def expected(self, k, size):
    ''' expected number under given sample size
    '''
    p = self.pmf(k)
    return size * p
  def estimate_truncated(self, data, size, max_iter = 1e4, nlike = 10, report = False): 
    '''estimate with truncated data. Not finished!
    '''
    total, cnt = data_count(data)
    lastllh = 0
    i = 0
    km = max(data)
    for j in range(int(max_iter)) :
      ps = 0
      for k in data : ps += self.pmf(k)
      #size = round(total / ps)
      k = 0
      # generate expected full data
      d = {}
      while k <= km or d[k-1] >= 1 : 
        if k in data : d[k] = data[k]
        else : d[k] = round(self.expected(k, size)) ## round?
        k += 1
      self.estimate(d)
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(d)
      diff = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if diff < self.Delta : break
      if report : 
        s = '{} {} {} {} {} {}'.format(diff, llh, self.r, self.p, k, size)
        for i in range(15) : s += ' %d:%d' % (i, d[i])
        s += '...'
        print(s)
      lastllh = llh
    return self.r, self.p
  def estimate_by_012(self, n0, n1, n2, start = 0) : 
    '''estimate by first 3 numbers
    '''
    r1 = (start + 1) * float(n1) / n0 
    r2 =  (start + 2) * float(n2) / n1
    p1 = r2 - r1
    self.p = 1 - p1 # 
    self.r = (r1 - start * p1) / p1
    return self.r, self.p
  def fit_lower(self, data, nmax = 40, pth = 0.01, start = 0) : # in test, not used
    r, p = self.estimate_by_012(data[start], data[start+1], data[start+2], start = start)#, nlike=100)
    if r > 0 and 0 < p < 0.9 : 
      #self.r, self.p = r, p
      return r, p
    else : 
      import numpy as np
      from scipy.optimize import leastsq
      def res(p, x, y) : 
        return p[0] * x + p[1] - y
      m = max(data)
      if m > nmax : m = nmax # max
      s = data[start] + data[start+1] + data[start+2]
      y = [(start + 1) * float(data[start+1]) / data[start], (start + 2) * float(data[start+2]) / data[start+1]]
      w = [(data[start+1] * data[start])**0.5, (data[start+1] * data[start+2])**0.5]
      for i in range(start+3, m): 
        y.append(i * float(data.value(i)) / data.value(i-1))
        w.append((data[i] * data[i-1])**0.5)
        wm = min([wi for wi in w if wi > 0])
        w2 = [int(round(wi / wm)) for wi in w]
        s += data.value(i)
        xl = sum([[j] * w2[j] for j in range(i)], [])
        yl = sum([[y[j]] * w2[j] for j in range(i)], [])
        xa = np.array(xl, dtype = float)
        ya = np.array(yl, dtype = float)
        rst = leastsq(res, [1,1], args=(xa, ya))
        p = 1 - rst[0][0]
        r = rst[0][1] / rst[0][0]
        print (i, (r, p), len(xl), w)#, y
        if r <= 0 or p >= 1 or p <= 0 : continue
        if 0.2 < p < 1 : break
        self.r, self.p = r, p
        #f1 = self.pmf(i) / self.cdf(i)
        f1 = self.pmf(i+1) / (self.cdf(i+1) - self.cdf(start-1))
        f2 = self.pmf(i+2) / (self.cdf(i+2) - self.cdf(start-1))
        #pv1 = binom_test(data.value(i), s, p = f1, alt = 'g')
        pv1 = binom_test(data.value(i+1), s+data.value(i+1), p = f1, alt = 'g')
        pv2 = binom_test(data.value(i+2), s+data.value(i+1)+data.value(i+2), p = f2, alt = 'g')
        #for j in range(i) : print j, self.pmf(j),
        print (data.value(i+1), s+data.value(i+1), f1, 'pv1 =', pv1)
        print (data.value(i+2), s+data.value(i+1)+data.value(i+2), f2, 'pv2 =', pv2)
        #fp, fs = fisher_method([pv1, pv2])
        fp = max(pv1, pv2)
        #print fp
        if fp < pth : break #and pv2 < pth : break
      self.r, self.p = r, p
      return r, p
  def fit_linear(self, data, poisson = False, maxi = 40, minc = 3, start = 0, total = None, show = False) : # in test, not used
    if poisson : N2 = 0
    else : 
      if total is None : total = sum(data.values())
      N2 = 2.0 / total
    rate, rvar = {}, {}
    for i in data : 
      if i - 1 < start : continue
      if i-1 not in data or data[i-1] <= 0 : continue
      if i not in data or data[i] <= 0 : continue
      if data[i] < minc : break
      #if data.value(i-1) == 0 or data.value(i) == 0 : continue
      if i >= maxi : break
      rate[i] = float(i) * data[i] / data[i-1]
      rvar[i] = rate[i] ** 2 * (1.0 / data[i] + 1.0 / data[i] - N2)
    #print 'rate =', rate
    #print 'rvar =', rvar
    c1, c2, c3, c4, c5 = 0,0,0,0,0
    a, b = {}, {}
    r, p = {}, {}
    #ia = rate.keys()
    #ia.sort()
    ia = sorted(rate)
    for n, i in enumerate(ia) : 
      c1 += 1 / rvar[i]
      c2 += i / rvar[i]
      c3 += i * i / rvar[i]
      c4 += rate[i] / rvar[i]
      c5 += i * rate[i] / rvar[i]
      if n > 0 : 
        a[i] = (c3 * c4 - c2 * c5) / (c1 * c3 - c2 * c2)
        b[i] = (c2 * c4 - c1 * c5) / (c2 * c2 - c1 * c3)
        p[i] = 1 - b[i]
        r[i] = a[i] / b[i] + 1
        if show : print (i, a[i], b[i], r[i], p[i])
    self.r, self.p = r[i], p[i]
    return a, b, r, p, rate, rvar
    r, p = self.estimate_by_012(data[start], data[start+1], data[start+2], start = start)#, nlike=100)
    if r > 0 and 0 < p < 0.9 : 
      #self.r, self.p = r, p
      return r, p
    else : 
      import numpy as np
      from scipy.optimize import leastsq
      def res(p, x, y) : 
        return p[0] * x + p[1] - y
      m = max(data)
      if m > nmax : m = nmax # max
      s = data[start] + data[start+1] + data[start+2]
      y = [(start + 1) * float(data[start+1]) / data[start], (start + 2) * float(data[start+2]) / data[start+1]]
      w = [(data[start+1] * data[start])**0.5, (data[start+1] * data[start+2])**0.5]
      for i in range(start+3, m): 
        y.append(i * float(data.value(i)) / data.value(i-1))
        w.append((data[i] * data[i-1])**0.5)
        wm = min([wi for wi in w if wi > 0])
        w2 = [int(round(wi / wm)) for wi in w]
        s += data.value(i)
        xl = sum([[j] * w2[j] for j in range(i)], [])
        yl = sum([[y[j]] * w2[j] for j in range(i)], [])
        xa = np.array(xl, dtype = float)
        ya = np.array(yl, dtype = float)
        rst = leastsq(res, [1,1], args=(xa, ya))
        p = 1 - rst[0][0]
        r = rst[0][1] / rst[0][0]
        print (i, (r, p), len(xl), w) #, y
        if r <= 0 or p >= 1 or p <= 0 : continue
        if 0.2 < p < 1 : break
        self.r, self.p = r, p
        #f1 = self.pmf(i) / self.cdf(i)
        f1 = self.pmf(i+1) / (self.cdf(i+1) - self.cdf(start-1))
        f2 = self.pmf(i+2) / (self.cdf(i+2) - self.cdf(start-1))
        #pv1 = binom_test(data.value(i), s, p = f1, alt = 'g')
        pv1 = binom_test(data.value(i+1), s+data.value(i+1), p = f1, alt = 'g')
        pv2 = binom_test(data.value(i+2), s+data.value(i+1)+data.value(i+2), p = f2, alt = 'g')
        #for j in range(i) : print j, self.pmf(j),
        print (data.value(i+1), s+data.value(i+1), f1, 'pv1 =', pv1)
        print (data.value(i+2), s+data.value(i+1)+data.value(i+2), f2, 'pv2 =', pv2)
        #fp, fs = fisher_method([pv1, pv2])
        fp = max(pv1, pv2)
        #print fp
        if fp < pth : break #and pv2 < pth : break
      self.r, self.p = r, p
      return r, p
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    #obs.append(ob)
    #exs.append(ex)
    obs[-1] += ob
    exs[-1] += ex
    print (obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs))
    return chisquare(obs, exs)
  
class ZTNB(NegBinom):
  '''Zero truncated negative binomial
  '''
  def logpmf(self, k = 1):
    if k < 1 : return NegBinom.pmf(self, -1)
    p0 = NegBinom.pmf(self, 0)
    lp = NegBinom.logpmf(self, k)
    return lp - math.log(1 - p0)
  def pmf(self, k = 1):
    if k < 1 : return 0
    p0 = NegBinom.pmf(self, 0)
    p = NegBinom.pmf(self, k)
    return p / (1- p0)
  def expected_zeros(self, size):
    p0 = NegBinom.pmf(self, 0)
    return size * p0 / (1 - p0)
  def estimate(self, data, max_iter = 1e4, nlike = 10, report = False):
    '''estimate ZTNB parameters from input data (dict type, value -> counts)
    '''
    d1 = {}
    for i in data : 
      if i > 0 : d1[i] = data[i]
    total, cnt = data_count(d1)
    lastllh = 0
    i = 0
    for j in range(int(max_iter)) :
      ez = self.expected_zeros(cnt)
      d1[0] = ez
      NegBinom.estimate(self, d1)
      d1[0] = 0
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(d1)
      d = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if d < self.Delta : break
      if report : print (d, llh, self.r, self.p, ez)
      lastllh = llh
    return self.r, self.p
  def pvalue(self, k = 1):
    if k <= 1 : return 1
    p0 = NegBinom.pmf(self, 0)
    p = NegBinom.pvalue(self, k)
    return p / (1- p0)
  def expect(self):
    p0 = NegBinom.pmf(self, 0)
    return NegBinom.expect(self) / (1 - p0)
  def variance(self):
    nb = NegBinom(self.r, self.p)
    p0 = nb.pmf(0)
    return (nb.variance() + nb.expect() ** 2) / (1 - p0) - self.expect() ** 2
    #return self.expect / self.p

class Poisson: 
  '''Poisson distribution
  '''
  lMax = 1e8
  lMin = 1e-8
  Delta = 1e-8
  def __init__(self, l = 1.0):
    self.l = l # lambda
  def expect(self):
    return self.l
  def variance(self):
    return self.l
  def logpmf(self, k = 0):
    logarr_ext(k)
    lpr = k * math.log(self.l) - self.l
    lpr -= logsumarr[k]
    return lpr
  def pmf(self, k = 0):
    return math.exp(self.logpmf(k))
  def cdf(self, k = 0): #, logarr = logarr):
    if k < 0 : return 0
    #logarr_ext(k)
    lpr = self.logpmf(0)
    logl = math.log(self.l)
    cdf = math.exp(lpr)
    for i in range(1, k + 1):
      lpr += logl - math.log(i) # logarr[i]
      cdf += math.exp(lpr)
    return cdf
  def pvalue(self, k = 0):
    return 1 - self.cdf(k-1)
  def estimate(self, data): #data dict value:counts
    '''estimate Poisson parameter from input data (dict type, value -> counts)
    '''
    total, cnt = data_count(data)
    self.l = float(total) / cnt
    return self.l
  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, l = -1):
    if l < 0 : l = self.l
    total, cnt = data_count(data)
    return float(total) / l / cnt -1
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    obs[-1] += ob
    exs[-1] += ex
    print (obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs))
    return chisquare(obs, exs)

class ZTPoisson(Poisson): 
  '''Zero truncated poisson distribution, NOT finished!
  '''
  def expect(self):
    p0 = Poisson.pmf(self, 0)
    return self.l / (1 - p0)
  def variance(self):
    p0 = Poisson.pmf(self, 0)
    return (self.l + self.l ** 2) / (1 - p0) - self.expect() ** 2
  def logpmf(self, k = 0):
    if k < 1 : return Poisson.pmf(self, -1)
    p0 = math.exp(Poisson.logpmf(self, 0))
    lp = Poisson.logpmf(self, k)
    return lp - math.log(1 - p0)
  #def pmf(self, k = 0):
    #return math.exp(self.logpmf(k))
  def cdf(self, k = 1): #, logarr = logarr):
    if k <= 0 : return 0
    #logarr_ext(k)
    lpr = self.logpmf(1)
    logl = math.log(self.l)
    cdf = math.exp(lpr)
    for i in range(2, k + 1):
      lpr += logl - math.log(i) # logarr[i]
      cdf += math.exp(lpr)
    return cdf
  #def pvalue(self, k = 0):
    #return 1 - self.cdf(k - 1)
  def estimate(self, data, max_iter = 1e4, nlike = 10): ###### NOT FINISHED!!
    total, cnt = data_count(data)
    lastllh = 0
    i = 0
    for j in range(max_iter) :
      ez = self.expected_zeros(cnt)
      data[0] = ez
      negbinom.estimate(self, data)
      data[0] = 0
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(data)
      d = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if d < self.Delta : break
      print (d, llh, self.r, self.p, ez)
      lastllh = llh
    return self.r, self.p
  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, r = -1):
    if r < 0 : r = self.r
    total, cnt = data_count(data)
    #dr = scipy.special.digamma(r)
    #score = math.log(self.q) - dr
    s1, d = 0, 0
    for i in range(max(data.keys()) + 1):
      #d += 1.0 / (r + i)
      if i in data : s1 += d * data[i]
      d += 1.0 / (r + i)
    score = s1 / cnt + math.log(r / (r + 1.0*total/cnt))
    return score
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    #obs.append(ob)
    #exs.append(ex)
    obs[-1] += ob
    exs[-1] += ex
    print (obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs))
    return chisquare(obs, exs)

class rankSumTiesExact: 
  '''exact rank sum test of x < y (one-sided) for ribo
  '''
  def __init__(self, x, y, show = False):
    '''x, y are list of samples
    '''
    self.a = list(x) + list(y)
    self.x = list(x)
    self.N, self.n = len(self.a), len(x)
    self.cd = countDict(self.a) # count dict
    self.xcd = countDict(x)
    self.rd = rankDict(self.cd)
    self.xrs =  self.rankSum(self.xcd)
    self.ks = sorted(self.cd)
    self.l = len(self.ks)
    self.vs = [self.count(i) for i in range(self.l)]
    if show : 
      print (self.cd, self.xcd, self.complexity()) #print (self.tieRatio()) print (self.complexity()) #print self.numStats()
  def copy(self, other):
    self.a = other.a
    self.x = other.x
    self.N, self.n = other.N, other.n
    self.cd = other.cd # count dict
    self.xcd = other.xcd
    self.rd = other.rd
    self.xrs =  other.xrs
    self.ks = other.ks
    self.l = other.l
    self.vs = other.vs
  def complexity(self):
    '''sum of log(n+1) for each ties
    '''
    #logarr_ext(max(self.cd.values()) + 1)
    complog = sum([math.log(self.count(i)+1) for i in range(self.l)])
    return complog
  def shuffleTest(self, n = 1000, show = False):
    import random
    c = 0
    a = self.a[:]
    for i in range(n):
      random.shuffle(a)
      rcd = countDict(a[0:self.n])
      rs =  self.rankSum(rcd)
      if rs >= self.xrs: c += 1
      if show : print (a, rs, xrs, c)
      #if c >= 100 : return float(c) / (i+1) # 2 effective numbers
    #if i == 19 and c > 10 : return float(c) / 20
    #if i == 99 and c > 20 : return float(c) / 100
    return float(c) / n
  def count(self, i):
    return self.cd[self.ks[i]]
  def rank(self, i):
    return self.rd[self.ks[i]]
  def numStats(self):
    #logarr_ext(max(self.cd.values()) + 1)
    return self._numStats(self.N, self.n, 0)
  def _numStats(self, N, n, i): # NOT CORRECT!!
    n1 = N - n
    if n < 0 or n1 < 0 : return 0
    if n == 0 or n1 == 0 : return 1
    if i == self.l - 2 : return min(n, n1) + 1
    if max(self.vs[i:]) == 1 : return math.exp(combination_log(N,n))
    t1, t2 = self.count(i), N - self.count(i)
    if t1 >= n >= t2 or t1 >= n1 >= t2 : 
      d = math.exp(sum([math.log(self.count(k)+1) for k in range(i+1, self.l)]))
      return d
      #print 'multi all', d, N, n, i
    s = 0
    for ni in range(max(0, n-t1),min(self.count(i), n) + 1) :
      d = self._numStats(N - self.count(i), n - ni, i+1)
      s += d
      #print d,s,N - ni, n - ni, i+1
    return s
  def numStatsRaw(self): # use complexity() instead
    s = 1
    for i in range(self.l) :
      s *= self.count(i) + 1
    return s
  def tieRatio(self):
    '''tie correction ratio
    '''
    s = sum([v**3 - v for v in self.cd.values() if v > 1])
    return s / float(self.N ** 3 - self.N)
  def test(self, th = 10, delta = 1e-4):
    '''if either size <= th, use fastTest
    elif complexity <= th, use fastTest
    else use normal mwtest
    '''
    if self.n <= th or self.N - self.n <= th : return self.fastTest(delta = delta) #, True
    if self.complexity() <= th: return self.fastTest(delta = delta) # * logarr[2]
    return self.mwtest() #, False
  def mwtest(self, use_continuity = True, show = False):
    '''normal mwtest, alt = greater, one-tailed
    '''
    from scipy.stats import norm
    n1, n2, n = self.n, self.N - self.n, self.N
    mu = n1 * (n + 1) / 2.0
    s = sum([v**3 - v for v in self.cd.values() if v > 1])
    if show : print ('effect size {}'.format(n - s / float(n*(n-1))))
    var = n1 * n2 * (n ** 3 - n - s) / 12.0 / n / (n - 1)
    if use_continuity : z = (self.xrs - mu - 0.5) / var ** 0.5
    else : z = (self.xrs - mu) / var ** 0.5
    p = norm.sf(z)
    if show : print (self.xrs, mu, self.xrs-mu+n1*(n1+1)/2.0, var, z, p)
    return p
    #else : return 1 - p # one sided
  def isExtreme(self, rs, twotailed = False, alt = 'g', delta = 1e-5): # not used. if the given rank sum is farther than xrs
    if twotailed : 
      if not hasattr(self, 'mu'):
        n1, n2, n = self.n, self.N - self.n, self.N
        self.mu = n1 * (n + 1) / 2.0
        self.th = abs(self.xrs - self.mu)
      return abs(prs - self.mu) >= self.th - delta
    elif alt == 'g' : return rs >= self.xrs - delta
    else : return rs <= self.xrs + delta
  def exactTest(self, twotailed = False, show = False): 
    '''get exact p-value for ranksum permutation test, time consuming.
    '''
    pval = 0
    if twotailed : 
      n1, n2, n = self.n, self.N - self.n, self.N
      mu = n1 * (n + 1) / 2.0
      th = abs(self.xrs - mu)
    for pcd in self.multiHypergeoIter(self.n, (), 0, self.N):
      prs =  self.rankSum(pcd)
      #if show : print prs, pcd
      if twotailed : 
        if abs(prs - mu) < th - 0.0001 : continue
      elif prs < self.xrs - 0.0001 : continue
      p = multiHypergeoProb(pcd, self.cd)
      pval += p
      if show : print (prs, pcd, p)
    if pval > 1.0 : pval = 1.0
    return pval
  def fastTest(self, show = False, delta = 1e-4):
    '''fast version of exactTest, reduced time at some cost of accuracy.
    '''
    self.pval = 0
    self.lp0 = - combination_log(self.N, self.n)
    for pcd in self.multiHypergeoFastIter(self.n, (), 0, self.N, 0, delta=delta):
      p = multiHypergeoMergeProb(pcd, self.cd, self.N, self.n, lp0 = self.lp0) # merged probability
      self.pval += p
      if show : print (p, pcd)
    if self.pval > 1.0 : self.pval = 1.0
    return self.pval
  def multiHypergeoIter(self, n, vs, i, N):
    if n == 0 : yield zipDict(self.ks, vs) # no more, return
    elif n > N : return # impossible
    elif i + 1 == len(self.ks) : # the last key
      if n <= self.cd[self.ks[i]] : yield zipDict(self.ks, vs + (n,))
      return
    else :
      N1 = N - self.count(i)
      for j in range(max(0, n-N1), min(n, self.count(i)) + 1):
        for pcd in self.multiHypergeoIter(n-j, vs+(j,), i+1, N1):
          yield pcd
  def multiHypergeoFastIter(self, n, vs, i, N, rs, delta=1e-4): 
    '''fast iter, only select rank sum higher than x rank sum conditions
    '''
    if n > N : return # impossible
    elif n == 0 : yield zipDict(self.ks, vs) # no more, return
    elif i + 1 == len(self.ks) : # the last key
      if n <= self.count(i) : yield zipDict(self.ks, vs + (n,))
      return
    else :
      N1 = N - self.count(i)
      jd, ju = max(0, n-N1), min(n, self.count(i)) + 1
      rsu, rsd = {}, {}
      j1, j2 = jd, ju
      while j1 < j2 : # looking for threshold of all lower conditions
        j = (j1 + j2) // 2 # + 1 for upper int
        rsu[j] = self.rankSumUpLimit(vs+(j,), rs+j*self.rank(i))
        if rsu[j] < self.xrs - 0.0001 : j2 = j # only keep possible conditions
        else : j1 = j + 1
      ju = j1
      j1, j2 = jd, ju
      while j1 < j2 : # looking for threshold of all higher conditions
        j = (j1 + j2) // 2 # + 1 for upper int
        rsd[j] = self.rankSumDownLimit(vs+(j,), rs+j*self.rank(i))
        if rsd[j] >= self.xrs - 0.0001 : j1 = j + 1 # only keep partial possible conditions
        else : j2 = j
      for j in range(jd, j1) : yield zipDict(self.ks, vs+(j,))
      ps = [multiHypergeoMergeProb(zipDict(self.ks, vs+(j,)), self.cd, self.N, self.n, lp0=self.lp0) for j in range(j1, ju)]
      for jo in orderIter(ps, reverse = True): #range(j1, ju):
        j = j1 + jo
        vsj, rsj = vs+(j,), rs+j*self.rank(i)
        if self.pval > 0 and ps[jo] / self.pval < delta : # whether prob. high enough to calculate in detail
          self.pval += ps[jo] / 2 # 0 order approximation; ps[jo] * ratio # 1st order approximation
        else :
          for pcd in self.multiHypergeoFastIter(n-j, vsj, i+1, N1, rsj, delta=delta):
            yield pcd
  def rankSum(self, cd):
    return sum([self.rd[k] * cd[k] for k in cd])
  def rankSumUpLimit(self, vs, rs): # rank sum up & down limit
    n, i = self.n - sum(vs), len(self.ks) - 1
    while n >= 1 :
      if n > self.count(i) : d = self.count(i)
      else : d = n
      rs += d * self.rank(i)
      n -= d
      i -= 1
    return rs
  def rankSumDownLimit(self, vs, rs): # rank sum up & down limit
    n, i = self.n - sum(vs), len(vs)
    while n >= 1 : 
      #print n, i
      if n > self.count(i) : d = self.count(i)
      else : d = n
      rs += d * self.rank(i)
      n -= d
      i += 1
    return rs
def countDict(arr):
  '''count array to dict
  '''
  data = {}
  for d in arr:
    if d not in data : data[d] = 0
    data[d] += 1
  return data
def rankDict(data):
  '''count array to rank dict
  '''
  ks = sorted(data)
  total = 0
  rd = {}
  for k in ks:
    total += data[k]
    rd[k] = total - (data[k] - 1) / 2.0
  return rd

def zipDict(ks, vs):
  data = {}
  for i, v in enumerate(vs) : data[ks[i]] = v
  return data

def orderIter(arr, reverse = False):
  '''generate order of given array. e.g. the index of smallest item come first
  '''
  ad = {}
  for i, v in enumerate(arr): 
    if v not in ad : ad[v] = []
    ad[v].append(i)
  a = sorted(ad, reverse = reverse)
  for v in a : 
    for i in ad[v] : yield i

def order(arr, reverse = False):
  return [i for i in orderIter(arr, reverse = reverse)]

def multiHypergeoProb(pcd, cd, show = False):
  '''pmf of multi-hypergeomic distribution, cd provides total numbers, pcd provides selected numbers
  '''
  lp = - combination_log(sum(cd.values()), sum(pcd.values()))
  if show : print (sum(cd.values()), sum(pcd.values()), lp)
  for k in pcd :
    lp += combination_log(cd[k], pcd[k])
    if show : print (cd[k], pcd[k], lp)
  return math.exp(lp)
def multiHypergeoMergeProb(pcd, cd, N, n, lp0 = None, show = False):
  '''keys not in pcd/cd are not considered. total numbers are provided by N & n
  '''
  if lp0 is None : lp = - combination_log(N, n)
  else : lp = lp0
  if show : print (lp)
  for k in pcd :
    lp += combination_log(cd[k], pcd[k])
    n -= pcd[k]
    N -= cd[k]
    if show : print (cd[k], pcd[k], N, n, lp)
  lp += combination_log(N, n)
  return math.exp(lp)

def glmNBTest(x, y):
  '''test difference of two array of counts by generized linear model with negative binomial family
  require statsmodels module
  '''
  import statsmodels.formula.api as smf
  import statsmodels.api as sm
  if max(x) == 0 : return 1.
  data = {}
  data['data'] = list(x) + list(y)
  data['grp'] = [1]*len(x) + [0]*len(y)
  model=smf.glm("data ~ grp", data=data, family=sm.families.NegativeBinomial()).fit()
  p = model.pvalues[1] / 2
  if model.tvalues[1] >= 0 : return p
  else : return 1-p

class betaBinom:
  ''' Beta-Binomial distribution model
  '''
  def __init__(self, a, b): # alpha & beta
    self.a = a
    self.b = b
  def logpmf(self, n, k):
    return combination_log(n, k) + betaln(k + self.a, n - k + self.b) - betaln(self.a, self.b)
  def pmf(self, n, k):
    return math.exp(self.logpmf(n, k))
  def cdf(self, n, k):
    p = 0
    for i in range(k+1):
      p += self.pmf(n, i)
    if p > 1 : p = 1
    return p
  def pvalue(self, n, k, alt = 'two.sided'):
    cdf = self.cdf(n, k)
    if alt in ('l', 'less') : return cdf
    else :
      p = 1 - cdf + self.pmf(n, k)
      if alt in ('g', 'greater') : return p
      else :
        m = min(cdf, p) * 2
        if m > 1 : m = 1
        return m
    #if alt not in ('g', 'greater') : return cdf
    #else :
      #p = self.pmf(n, k)
      #return 1 - cdf + p

