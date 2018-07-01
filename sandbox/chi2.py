import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt

df = 1

x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)

y = chi2.cdf(x, df)

# ysum = np.cumsum(y)/np.sum(y)

ciy = y[np.where(y>=0.9)][0]
cix = x[np.where(y>=0.9)][0]

plt.plot(x, y, '-b')
plt.axhline(ciy, c='g', label="%.2f" % ciy)
plt.axvline(cix, c='g', label="%.2f" % cix)
plt.legend()

print(cix/2)


# ax.plot(x, chi2.cdf(x, df), '-b', lw=3,  label='chi2')

# ymax = chi2.cdf(0.9, df)
# ax.axvline(ymax, c='g')

# print(chi2.pdf(0.90, 1))
# vals = chi2.ppf([0.001, 0.5, 0.999], df)
# print(vals)


# Double_t Yat_Xmax = 0.5*ROOT::Math::chisquared_quantile(fInterval->ConfidenceLevel(),1);

plt.show()