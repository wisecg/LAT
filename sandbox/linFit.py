

    # fs10 = 1.74 kev
    # fs12 = 1.84 kev
    # fs24 = 2.31 kev
    # fs28 = 2.46 kev

import plotly.plotly as py
import plotly.graph_objs as go

# MatPlotlib
import matplotlib.pyplot as plt
from matplotlib import pylab

# Scientific libraries
from numpy import arange,array,ones
from scipy import stats


y = [1.74, 1.84, 2.31, 2.46]
x = [10., 12., 24., 28.]

# Generated linear fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

print slope

# line = slope*x+intercept
#
# plt.plot(x,y,'o', x, line)
# pylab.title('Linear Fit with Matplotlib')
# ax = plt.gca()
# ax.set_axis_bgcolor((0.898, 0.898, 0.898))
# fig = plt.gcf()
# py.plot_mpl(fig, filename='linear-Fit-with-matplotlib')