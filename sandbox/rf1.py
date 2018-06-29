#!/usr/bin/env python3
import ROOT

def main():

    # Build Gaussian PDF
    x = ROOT.RooRealVar( 'x', 'x', -10, 10 )
    mean = ROOT.RooRealVar( 'mean', 'mean of gaussian', -1 )
    sigma = ROOT.RooRealVar( 'sigma', 'width of gaussian', 3 )
    gauss = ROOT.RooGaussian( 'gauss', 'gaussian PDF', x, mean, sigma )

    # Plot PDF
    c = ROOT.TCanvas("c","c",800,600)
    xframe = x.frame()
    gauss.plotOn( xframe, ROOT.Components(10))
    xframe.Draw()
    c.Print("../plots/rf-test.pdf")

    # # Generate a toy MC set
    # data = gauss.generate( RooArgSet(x), 10000 )
    # # Plot PDF and toy data overlaid
    # xframe2 = x.frame()
    # data.plotOn( xframe2, RooLinkedList() )
    # gauss.plotOn( xframe2 )
    # xframe2.Draw()
    # # Fit PDF to toy
    # mean.setConstant( ROOT.kFALSE )
    # sigma.setConstant( ROOT.kFALSE )
    # gauss.fitTo( data, 'mh' )
    # # Print final value of parameters
    # mean.Print()

if __name__=="__main__":
    main()