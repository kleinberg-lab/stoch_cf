# Causal Explanation in Stochastic Models using Counterfactuals
# Algorithm Module

import math
import random
import numpy as np
from scipy import stats

# Assumptions and Simplifications
# 1) discrete time 

class Algorithm1(object):
  def __init__(self,mu,w_0,tau,n,epsilon,Z,H):
    self.mu = mu
    self.w_0 = w_0
    self.tau = tau
    self.n = n
    self.epsilon = epsilon
    self.Z = Z
    self.H = H
    self.t_now = len(w_0)

  def discount(self,t_ref):
    t_diff = t_ref - self.tau
    return t_diff

  def characteristic_samples(self,A,B):
    # for actual events A
    assert(A(self.w_0))
    W_cf = []
    W_def = []
    k = 0
    kC = 0
    kD = 0
    while k < self.n:
      while True:
        t = self.discount(self.Z(A,self.w_0))
        w_cf = self.mu(self.w_0,t)
        if not A(w_cf):
          break
      while True:
        w_def = self.mu(self.w_0,t)
        if A(w_def):
          break
      print int(t)
      W_cf.append(w_cf)
      W_def.append(w_def)
      k = k + 1
      if B(w_cf):
        kC += 1
      if B(w_def):
        kD += 1
      print "Generated %d samples so far... (B true in %d/%d)"%(k,kC,kD)
    return (W_cf, W_def)

  def explain(self,B,A):
    # CF of B on A
    (W_cf, W_def) = self.characteristic_samples(A,B)
    P_cf = len(filter(B,W_cf))/float(self.n)
    P_def = len(filter(B,W_def))/float(self.n)
    P_delta = P_def - P_cf
    p_P = lltest(P_cf*self.n,P_def*self.n,self.n)
    print "P_cf = %f"%P_cf
    print "P_def = %f"%P_def
    print "P_delta = %s"%str(P_delta)
    print "P_cf != P_def with p = %s"%str(p_P)
    def mean_time(E,distr):
      times = [self.Z(E,w) for w in distr if E(w)]
      mu = np.mean(times)
      sigma = np.std(times)
      return mu, sigma, times
    t_mu_cf, t_sigma_cf, times_cf = mean_time(B,W_cf)
    t_mu_def, t_sigma_def, times_def = mean_time(B,W_def)
    t_delta = t_mu_def - t_mu_cf
    p_t = stats.ttest_ind(times_cf, times_def, equal_var = False)[1]
    print "t_cf = %f (%f)"%(t_mu_cf,t_sigma_cf)
    print "t_def = %f (%f)"%(t_mu_def,t_sigma_def)
    print "t_delta = %s"%str(t_delta)
    print "t_cf != t_def with p = %s"%str(p_t)
    def mean_intensity(E,distr):
      intensities = [self.H(E,w) for w in distr if E(w)]
      mu = np.mean(intensities)
      sigma = np.std(intensities)
      return mu, sigma, intensities
    h_mu_cf, h_sigma_cf, intensities_cf = mean_intensity(B,W_cf)
    h_mu_def, h_sigma_def, intensities_def = mean_intensity(B,W_def)
    h_delta = h_mu_def - h_mu_cf
    p_h = stats.ttest_ind(intensities_cf, intensities_def, equal_var = False)[1]
    print "h_cf = %f (%f)"%(h_mu_cf,h_sigma_cf)
    print "h_def = %f (%f)"%(h_mu_def,h_sigma_def)
    print "h_delta = %s"%str(h_delta)
    print "h_cf != h_def with p = %s"%str(p_h)
    return P_delta, p_P, t_delta, p_t, times_cf, times_def, W_cf, W_def

def lltest(k1h,k2h,N):
  # returns p-value for p1 != p2
  #  using log-likelihood ratio statistic with chi2 test
  N = float(N)
  k1t = N - k1h
  k2t = N - k2h
  try:
    return stats.chi2_contingency([[k1h,k2h],[k1t,k2t]],lambda_="log-likelihood")[1]
  except:
    return 1.0


