#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def main():

    fig = plt.figure(figsize=(8,7),facecolor='w')

    xData = np.arange(0,20000,10)

    A, t0, P, tau = 10, 10000, 150, 1000
    yData = 0.5 * np.tanh((xData-t0)/tau) + P

    plt.plot(xData,yData)
    plt.show()


if __name__=="__main__":
    main()