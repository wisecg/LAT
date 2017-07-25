#!/usr/local/bin/python
import pywt
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import waveLibs as wl

""" NOTE: These pywavelets pages were super useful.
    http://pywavelets.readthedocs.io/en/latest/ref/wavelet-packets.html
    http://pywavelets.readthedocs.io/en/latest/regression/wp.html

    Best basis/cost functions:
    http://www.bearcave.com/misl/misl_tech/wavelets/packet/index.html
"""

def main():

    reconstruct()

def reconstruct():
    """ Take only the lowest-frequency component of the WPT. """

    npzfile = np.load("./data/rawInputs.npz")
    raw = npzfile['arr_0']
    wave, waveTS, dataE, dataRT = raw[0],raw[1],raw[2],raw[3]

    wp = pywt.WaveletPacket(data=wave, wavelet='haar', mode='symmetric',maxlevel=3)
    new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
    new_wp['aaa'] = wp['aaa'].data
    newWave = new_wp.reconstruct(update=False)

    # plots
    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_title("Energy %.1f keV  RiseTime %.1f ns" % (dataE, dataRT))
    a1.set_ylabel("ADC [A.U.]",y=0.95, ha='right')
    a1.set_xlabel("Time (ns)",x=0.95, ha='right')
    a1.plot(waveTS,wave,color='blue',alpha=0.8,label='Original Waveform')
    a1.plot(np.arange(0,len(newWave)*10,10),newWave,color='red',label='Lvl3 Haar Denoised WF')
    a1.legend(loc=4)
    plt.show(block=False)
    plt.savefig("./plots/denoised-wf.pdf")

def bestBasis():
    npzfile = np.load("./data/rawInputs.npz")
    raw = npzfile['arr_0']
    wave, waveTS, dataE, dataRT = raw[0],raw[1],raw[2],raw[3]
    level = 5

    # create wp, loop over nodes and calculate cost
    wp = pywt.WaveletPacket(wave, 'db2', 'symmetric',maxlevel=level)
    costDic = {}
    for l in range(1, level+1):
        nodes = wp.get_level(l, order='freq')

        for n in wp.get_leaf_nodes():
            cost = entropyCost(n)
            # print n.path, "len: ", len(n.data), "cost",cost

            avg = thresholdCost(n)
            print n.path, "len: ", len(n.data), "avg %.3f" %avg

            costDic[n.path] = cost

    # empty wavelet packet to fill
    # new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')

    # IDEA: try doing best basis selection for the non-low frequency parts
    # (the ones that have at least one 'd')

    # Now compare costs
    # print "now doing second loop"
    # bestTree = {}
    # # Have to re-declare to make the loop work again, for some reason
    # wp = pywt.WaveletPacket(wave, 'db2', 'symmetric', maxlevel=level)
    # for l in range(1, level):
    #     print "now scanning level",l
    #     nodes = wp.get_level(l, order='freq')
    #
    #     for n in wp.get_leaf_nodes():
    #         subnode1 = n.path+'a'
    #         subnode2 = n.path+'d'
    #
    #         print "cost"


def thresholdCost(node,thresh=1):

    # cost=0
    avg=0
    for j in range(len(node.data)):
        avg += node.data[j]
    avg = avg / len(node.data)
    return avg



def entropyCost(node):
    """ Calculate cost with shannon entropy. """

    pTot = np.sum( np.power( np.abs(node.data) ,2) )
    cost = 0
    for j in range(len(node.data)):
        p_i = np.power( np.abs( node.data[j] ), 2 ) / pTot

        # skip zeros and nan's
        if p_i > 0 and not np.isnan(p_i):
            cost += p_i * np.log( p_i )

    return -1.0 * cost




def reconstructExample():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    # Now create a new Wavelet Packet and set its nodes with some data.
    new_wp = pywt.WaveletPacket(data=None, wavelet='db1', mode='symmetric')
    new_wp['aa'] = wp['aa'].data
    new_wp['ad'] = [-2., -2.]

    # For convenience, Node.data gets automatically extracted from the Node object:
    new_wp['d'] = wp['d']

    # And reconstruct the data from the aa, ad and d packets.
    print(new_wp.reconstruct(update=False))
    # [ 1.  2.  3.  4.  5.  6.  7.  8.]

    # If the update param in the reconstruct method is set to False, the node's data will not be updated.
    print(new_wp.data)
    # None

    # Otherwise, the data attribute will be set to the reconstructed value.
    print(new_wp.reconstruct(update=True))
    # [ 1.  2.  3.  4.  5.  6.  7.  8.]

    print(new_wp.data)
    # [ 1.  2.  3.  4.  5.  6.  7.  8.]

    print([n.path for n in new_wp.get_leaf_nodes(False)])
    # ['aa', 'ad', 'd']

    print([n.path for n in new_wp.get_leaf_nodes(True)])
    # ['aaa', 'aad', 'ada', 'add', 'daa', 'dad', 'dda', 'ddd']

def manipulateNodes():
    """ http://pywavelets.readthedocs.io/en/latest/regression/wp.html """

    # Create sample data
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    ts = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    # First, start with a tree decomposition at level 2. Leaf nodes in the tree are:
    dummy = wp.get_level(2)
    for n in wp.get_leaf_nodes(False):
        print(n.path, format_array(n.data))

    # access a single node
    node = wp['ad']
    print(node)

    # To remove a node from the WP tree, use Python's del obj[x] (Node.__delitem__):
    del wp['ad']

    # The leaf nodes that left in the tree are:
    for n in wp.get_leaf_nodes():
        print(n.path, format_array(n.data))

    # And the reconstruction is:
    # print wp.reconstruct()
    recon = wp.reconstruct()

    # Now restore the deleted node value.
    # wp['ad'].data = node.data

    # Printing leaf nodes and tree reconstruction confirms the original state of the tree:

    # for n in wp.get_leaf_nodes(False):
        # print(n.path, format_array(n.data))

    # aa [  5.  13.]
    # ad [-2. -2.]
    # da [-1. -1.]
    # dd [ 0.  0.]

    # Print data reconstruction
    # print(wp.reconstruct())
    # [ 1.  2.  3.  4.  5.  6.  7.  8.]

    fig = plt.figure(figsize=(8,8), facecolor='w')

    a1 = plt.subplot(211)
    a1.plot(ts,x,color='blue')

    a2 = plt.subplot(212)
    a2.plot(ts,recon,color='red')

    plt.show()

def format_array(a):
    a = np.where(np.abs(a) < 1e-5, 0, a)
    return np.array2string(a, precision=5, separator=' ', suppress_small=True)

def mjdWaveform():

    # Load a data and a temp waveform
    npzfile = np.load("./data/rawInputs.npz")
    raw = npzfile['arr_0']
    wave, waveTS, dataE, dataRT = raw[0],raw[1],raw[2],raw[3]

    wp = pywt.WaveletPacket(wave, 'db2', 'symmetric', maxlevel=4)
    nodes = wp.get_level(4, order='freq', decompose=True)
    yWT = np.array([n.data for n in nodes], 'd')
    yWT = abs(yWT)
    waveletYTrans = yWT  # final coeff's (used in process-waveforms)
    rec = wp.reconstruct(update=True)

    # plots
    fig = plt.figure(figsize=(8,8), facecolor='w')

    a1 = plt.subplot(311)
    a1.plot(waveTS,wave,color='blue')

    a2 = plt.subplot(312)
    a2.imshow(waveletYTrans, interpolation='nearest',
        aspect="auto", origin="lower",extent=[0, 1, 0, len(waveletYTrans)],cmap='jet')

    a3 = plt.subplot(313)
    a3.plot(waveTS,rec,color='red')

    plt.show()

def waveletTransform(signalRaw, level=4, wavelet='db2', order='freq'):
    """ Use PyWavelets to do a wavelet transform. """
    wp = pywt.WaveletPacket(signalRaw, wavelet, 'symmetric', maxlevel=level)
    nodes = wp.get_level(level, order=order)
    yWT = np.array([n.data for n in nodes], 'd')
    yWT = abs(yWT)
    return wp.data, yWT

def blogPost():
    """ this doesn't use wavelet packets, it just uses regular wavelets.  Adapted from
    http://connor-johnson.com/2016/01/24/using-pywavelets-to-remove-high-frequency-noise/
    """
    level=4

    npzfile = np.load("./data/rawInputs.npz")
    raw = npzfile['arr_0']
    wave, waveTS, dataE, dataRT = raw[0],raw[1],raw[2],raw[3]

    wavelet = 'db8'

    # calculate the wavelet coefficients
    coeff = pywt.wavedec( wave, wavelet, mode="per" )

    # calculate a threshold
    from statsmodels.robust import mad
    sigma = mad( coeff[-level] )

    # changing this threshold also changes the behavior,
    # but I have not played with this very much
    uthresh = sigma * np.sqrt( 2*np.log( len( wave ) ) )
    coeff[1:] = ( pywt.threshold( i, value=uthresh, mode="soft" ) for i in coeff[1:] )

    # reconstruct the signal using the thresholded coefficients
    recon = pywt.waverec( coeff, wavelet, mode="per" )

    # plot
    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.plot(waveTS,wave,color='blue', alpha=0.5)
    a1.plot(waveTS,recon,color='red', alpha=0.5)
    plt.show()


if __name__ == "__main__":
    main()
