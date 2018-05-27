#!/usr/bin/env python3
"""
    % Email chain: "Threshold info for skim files"
    % File: ~/Documents/mjd/MJD\ Reports/david_threshold_erfns.pdf

    % Also David's main talk: http://mjcalendar.npl.washington.edu/indico/event/2351/material/slides/0.pdf
    # using positive and negative trigger data?
    2016-9-15-P3LQG_Run60001538 C2P1D3 P427488

    Questions:
    - Have you worked out how to take into account the effect of the trap-max nature of the on-board threshold?

    - How different are the threshold functions from the standard erf shape?

      -- david:  We could evaluate that by fitting the trap-max E=0 noise peak with a Gaussian and looking at the residuals.
                 I think it's pretty darn close; but how much discrepancy would we want to allow?
      -- jason:  It's an issue of where the analysis threshold will be.  To first order, the lower it is, the more sensitive we
                 are to non-Gaussian.  It also matters if the non-gaussianity is similar or different on all the channels.
                 If it's different, a lot of the differences will average away in a joint fit over many channels.
                 Another option besides fitting gaussians to the E=0 peaks is to fit your threshold curves to erfs.
                 It should be mathematically identical but the latter puts the residuals on the scale that matters and might
                 make it easier to evaluate.
      -- david: (does study)
                The residual on the threshold curve is less than about 1%.  Do we think that's good enough?
      -- jason: Nice work.  You can try adding a skew gaussian, but i would think that 1\% uncertainty is good enough for now,
                especially since it wiggles about the efficiency curve, so when you sum up over many slight left-right shifts
                in the erf position they will largely cancel.

    David's study:
        Take one detector, two different trap filters (trap max and fixed-time pickoff) w/ same rise and flat times.
        Fit with Gaussians
        Integrate to make erfn-type threshold curves
        Look at residuals
        Erfn residuals are less than 1\%.
        Call it good.

    Clint's take:
        For my thesis we are staying 3 sigma above the trigger threshold (about 99%.)
        So we don't need to adjust the erf shape,
        and we've learned from David's study that the onboard trap-max trigger doesn't bias the results
        (we get the same result as if we used the fixed-time pickoff.)
        Also, we have great stats on this measurement for every background run, so we know we're staying at 99% efficiency.
        Ok, done.  Now let's make some plots.
"""
import numpy as np
from scipy.special import erf

def erf(x, mu, sig):
    return (1/2) * (1 + erf( (x-mu) / (np.sqrt(2) * sig) ) )


def main():

    print("hi")
    mu = 1.448
    sig = 0.305



if __name__=="__main__":
    main()
