"""Exam assignment.
"""

import argparse
import  ROOT
from . import *

#Include the header file
ROOT.gInterpreter.Declare('#include "include/exam_assignment.h"')


# Enable multi-threading
ROOT.gInterpreter.ProcessLine("ROOT::EnableImplicitMT()")



timer = ROOT.TStopwatch()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Program that find the Higgs boson \
    in the decay channel H->ZZ->4l with CMS Open data.')
    parser.add_argument('-nofast', action='store_const', default=True, const=False, \
    help="No-fast mode take data from raw file and it does the data selection. \
    This saves also data selected to data path.")
    parser.add_argument('-local', action='store_const', default=True, const=False, \
    help="Local mode take data from raw file in local and it does the data selection.")
    parser.add_argument('-time', action='store_const', default=True, const=False, \
    help="Print the execution time.")
    args = parser.parse_args()

    if(args.time==False):    timer.Start()

    raw_to_data_selected(args.nofast,args.local)

    data_rdf = ROOT.RDataFrame("Events", (f for f in ["data/2e2muToZZ.root",\
        "data/4muToZZ.root", "data/4eToZZ.root"]))

    Spectrum = data_rdf.Histo1D(("Spectrum", "", 36, 70, 180), "InvariantMass")
    h = Spectrum.GetValue()
    h.SetMarkerStyle(8)
    h.SetLineColor(1)
    h.GetYaxis().SetTitle("N_{Events}")
    h.GetYaxis().SetTitleSize(0.07)
    h.SetStats(0)

    c1 = ROOT.TCanvas()
    c1.cd()

    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
    upper_pad.Draw()
    lower_pad.Draw()
    upper_pad.cd()

    h.DrawCopy("PE1")

    montecarlo_selection(args.nofast,args.local)

    montecarlo_rdf = ROOT.RDataFrame("Events", (f for f in ["montecarlo/montecarlo_2e2muToZZ.root",\
        "montecarlo/montecarlo_4muToZZ.root", "montecarlo/montecarlo_4eToZZ.root"]))

    montecarlo_spectrum = montecarlo_rdf.Histo1D(("Spectrum", "", 36, 70, 180), "InvariantMass", \
    "weight")
    h1 = montecarlo_spectrum.GetValue()
    h1.SetLineColor(1)
    h1.SetFillColor(7)
    h1.SetFillStyle(3003)

    h1.DrawCopy("HISTO SAME")

    # Add legend
    legend = ROOT.TLegend(0.72, 0.70, 0.92, 0.88)
    legend.SetFillColor(0)
    legend.SetBorderSize(1)
    legend.SetTextSize(0.03)
    legend.AddEntry(h, "Data", "PE1")
    legend.AddEntry(h1, "ZZ montecarlo", "f")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.06)
    text.DrawLatex(0.38, 0.84, "CMS Open Data")
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.DrawLatex(0.38, 0.78, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}")

    lower_pad.cd()
    ratioplot = h.Clone()
    ratioplot.Add(h1, -1)
    for i in range(1, h.GetNbinsX()):
        ratioplot.SetBinError(i, h.GetBinError(i))
    
    #Create observables
    x = ROOT.RooRealVar("x", "x", 70, 180)
    x.setRange("signal",120,130)
    #Create Gaussian
    mean = ROOT.RooRealVar("mean", "mean of gaussian", 1, 120, 130)
    sigma = ROOT.RooRealVar("sigma", "width of gaussian", 1, 0.01, 10)
    gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)
    #Create a binned dataset that imports contents of ratioplot hist and associates
    #its contents to observable 'x'
    dh = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(x), ROOT.RooFit.Import(ratioplot))
    # Fit a Gaussian p.d.f to the data
    gauss.fitTo(dh,  ROOT.RooFit.Save(True), ROOT.RooFit.Range("signal"))
    #Now plot the data and the fit result
    xframe = x.frame()
    dh.plotOn(xframe)
    gauss.plotOn(xframe, ROOT.RooFit.Range(70, 180), ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(9))
    gauss.plotOn(xframe)
    #Optional
    gauss.paramOn(xframe, ROOT.RooFit.Layout(0.8,0.9), ROOT.RooFit.Format("NEU",1))
    xframe.SetTitle("")
    xframe.GetXaxis().SetTitle("m_{4l} (GeV)")
    xframe.GetXaxis().SetTitleSize(0.12)
    xframe.GetXaxis().SetLabelSize(0.08)
    xframe.GetYaxis().SetLabelSize(0.07)
    xframe.GetYaxis().SetTitle("Data - Bkg.")
    xframe.GetYaxis().SetTitleSize(0.09)
    xframe.GetYaxis().SetTitleOffset(0.4)
    xframe.Draw()
    #Print values of mean and sigma (that now reflect fitted values and
    #errors)
    mean.Print()
    sigma.Print()

    if(args.time==False):
        timer.Stop()
        print("Execution time {:0.2f} m".format(timer.RealTime()/60))

