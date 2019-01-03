#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('../clint.mpl')

import waveLibs as wl

def main():
    """ trying different distributions for the m2s238 efficiency functions """

    # asymmetric cdf of xgauss
    # dt = 1
    # ts = np.arange(0,1000,dt)
    # mu, sig, tau = 500, 50, -100
    # wf = wl.evalXGaus(ts,mu,sig,tau)
    # dw = [0]
    # for i in range(1,len(ts)):
    #     dw.append((wf[i]-wf[i-1])/dt)
    # iw = [wf[0]]
    # for i in range(1,len(ts)):
    #     iw.append(wf[i]+iw[i-1])
    # plt.plot(ts, wf/np.sum(wf))
    # plt.plot(ts, dw/np.sum(dw))
    # plt.plot(ts, iw)
    # plt.show()

    # from scipy.stats import skewnorm
    # a=3
    # x = np.linspace(skewnorm.ppf(0.01, a), skewnorm.ppf(0.99, a), 100)
    # plt.plot(x, skewnorm.cdf(x, a),'r-', label='skewnorm cdf')
    # plt.show()

    # compare logistic to genlogistic (fits a little better)
    # from scipy.stats import logistic
    # yL = amp * logistic.cdf(xE, mu, sig) # ok this is identical to the manual one
    # plt.plot(xE, yL, "-k", lw=5, label='regular')
    # from scipy.stats import genlogistic
    # skews = np.arange(0.9,2.0,0.1)
    # cmap = plt.cm.get_cmap('hsv',len(skews)+1)
    # for i,skew in enumerate(skews):
    #     yG = amp * genlogistic.cdf(xE, skew, mu, sig)
    #     plt.plot(xE, yG, "-", c=cmap(i), label='skew=%.1f' % skew)
    # plt.xlabel("Energy")
    # plt.ylim(0,1)
    # plt.legend(loc=4)

    # best-fit to all det's with logistic
    xE = np.arange(0, 30, 0.1)
    mu, sig, amp = 0.6, 5.54, 0.91
    from scipy.stats import logistic
    yE = amp * logistic.cdf(xE,mu,sig)
    plt.plot(xE, yE, "-b", lw=5, label='regular')

    # weibull is pretty good
    # from scipy.stats import weibull_min as wb
    # c is like sharpness (sig)
    # loc is the last zero (not mu)
    # c, loc, scale = 1, 1, 5
    # plt.plot(xE, 0.9 * wb.cdf(xE,c,loc,scale),'-',c='r')
    # plt.legend(loc=4)

    # xgauss cdf
    from scipy.stats import exponnorm as xg
    k, loc, scale, amp = 10, -1, 0.5, 1
    yG = amp * xg.cdf(xE, k, loc, scale)
    plt.plot(xE, yG, "-", c='r')
    plt.axhline(0.9,c='g',alpha=0.5)
    plt.ylim(0,1)
    plt.show()






if __name__=="__main__":
    main()