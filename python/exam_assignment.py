"""Exam assignment.
"""

import argparse


import  ROOT

#Include the header file
ROOT.gInterpreter.Declare('#include "include/exam_assignment.h"')


# Enable multi-threading
ROOT.gInterpreter.ProcessLine("ROOT::EnableImplicitMT()")


#Python functions

def filter_z_mass(rdf):
    '''This function select good Z mass events'''

    rdf_z_a_mass = rdf.Filter("Z_mass[0] > 40 && Z_mass[0] < 120",\
    "Mass of first Z candidate in [40, 120]")
    rdf_z_b_mass = rdf_z_a_mass.Filter("Z_mass[1] > 12 && Z_mass[1] < 120",\
    "Mass of second Z candidate in [12, 120]")


    return rdf_z_b_mass




def selection_2e2mu(rdf):
    '''This function select events with e+ e- mu+ mu-'''

    rdf_2e2mu = rdf.Filter("nElectron >= 2 && nMuon >= 2",\
    "First selection with at least 2 electrons and 2 muons")

    rdf_pt = rdf_2e2mu.Filter("All(Muon_pt>5) && All(Electron_pt>7)", "Pt cuts")

    rdf_eta = rdf_pt.Filter("All(abs(Electron_eta)<2.5) && All(abs(Muon_eta)<2.4)",\
    "Direction cuts (geometrical acceptance)")

    rdf_iso = rdf_eta.Filter("All(abs(Electron_pfRelIso03_all)<0.40) && \
    All(abs(Muon_pfRelIso04_all)<0.40)", "Require good isolation")

    #Definition of 3d impact parameter significance for Electrons
    rdf_el_3dip = rdf_iso.Define("Electron_3DIP", \
    "sqrt(Electron_dxy*Electron_dxy + Electron_dz*Electron_dz)")
    rdf_el_3dips = rdf_el_3dip.Define("Electron_3DIPS",\
    "Electron_3DIP/sqrt(Electron_dxyErr*Electron_dxyErr + Electron_dzErr*Electron_dzErr)")

    rdf_el_track = rdf_el_3dips.Filter("All(Electron_3DIPS<4) && All(abs(Electron_dxy)<0.5)\
    && All(abs(Electron_dz)<1.0)", "Muons and electrons come from the same vertex")

    #Definition of 3d impact parameter significance for Muons
    rdf_mu_3dip = rdf_el_track.Define("Muon_3DIP", "sqrt(Muon_dxy*Muon_dxy + Muon_dz*Muon_dz)")
    rdf_mu_3dips = rdf_mu_3dip.Define("Muon_3DIPS",\
    "Muon_3DIP/sqrt(Muon_dxyErr*Muon_dxyErr + Muon_dzErr*Muon_dzErr)")

    rdf_mu_track = rdf_mu_3dips.Filter("All(Muon_3DIPS<4) && All(abs(Muon_dxy)<0.5) && \
    All(abs(Muon_dz)<1.0)", "Muons and electrons come from the same vertex")

    rdf_charge = rdf_mu_track.Filter("Sum(Electron_charge) \
     == 0 && Sum(Muon_charge) == 0", "Selection for total 0 charge and 0 leptons flavour")

    #Definition of Z masses for filtering
    rdf_z_mass = rdf_charge.Define("Z_mass", "calculation_Z_mass_2el2mu(Electron_pt, \
    Electron_eta, Electron_phi, Electron_mass, Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
    rdf_cut = filter_z_mass(rdf_z_mass)


    return rdf_cut

def selection_4mu(rdf):
    '''This function select events with 2mu+ 2mu-'''

    rdf_4mu = rdf.Filter("nMuon >= 4", "First selection with at least 4 muons")

    rdf_kin = rdf_4mu.Filter("All(Muon_pt>5) && All(abs(Muon_eta)<2.4)", "Kinematics cuts")

    rdf_iso = rdf_kin.Filter("All(abs(Muon_pfRelIso04_all)<0.4)", "Require good isolation")

    #Definition of 3d impact parameter significance
    rdf_3dip = rdf_iso.Define("Muon_3DIP", "sqrt(Muon_dxy*Muon_dxy + Muon_dz*Muon_dz)")
    rdf_3dips = rdf_3dip.Define("Muon_3DIPS", \
    "Muon_3DIP/sqrt(Muon_dxyErr*Muon_dxyErr + Muon_dzErr*Muon_dzErr)")

    rdf_pv = rdf_3dips.Filter("All(Muon_3DIPS<4) && \
    All(abs(Muon_dxy)<0.5) && All(abs(Muon_dz)<1.0)", "Muons come from the same vertex")

    rdf_charge = rdf_pv.Filter("nMuon==4 && Sum(Muon_charge == 1) == 2 && \
    Sum(Muon_charge == -1)==2", "Selection for total 0 charge(two positive and two negative muons)")

    #Definition of Z masses for filtering
    rdf_z_mass = rdf_charge.Define("Z_mass", "calculation_Z_mass_4l\
    (Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge)")
    rdf_cut = filter_z_mass(rdf_z_mass)


    return rdf_cut

def selection_4e(rdf):
    '''This function select events with 2e+ 2e-'''

    rdf_4e = rdf.Filter("nElectron >= 4", "First selection with at least 4 electrons")

    rdf_kin = rdf_4e.Filter("All(Electron_pt>7) && All(abs(Electron_eta)<2.5)", "Kinematics cuts")

    rdf_iso = rdf_kin.Filter("All(abs(Electron_pfRelIso03_all)<0.40)", "Require good isolation")

    #Definition of 3d impact parameter significance for Electrons
    rdf_3dip = rdf_iso.Define("Electron_3DIP", \
    "sqrt(Electron_dxy*Electron_dxy + Electron_dz*Electron_dz)")
    rdf_3dips = rdf_3dip.Define("Electron_3DIPS", \
    "Electron_3DIP/sqrt(Electron_dxyErr*Electron_dxyErr + Electron_dzErr*Electron_dzErr)")

    rdf_pv = rdf_3dips.Filter("All(Electron_3DIPS<4) && \
    All(abs(Electron_dxy)<0.5) && All(abs(Electron_dz)<1.0)", "Electrons come from the same vertex")

    rdf_charge = rdf_pv.Filter("nElectron==4 && Sum(Electron_charge == 1) == 2 && \
    Sum(Electron_charge == -1) == 2", \
    "Selection for total 0 charge (two positive and two negative electrons)")

    #Definition of Z masses for filtering
    rdf_z_mass = rdf_charge.Define("Z_mass", "calculation_Z_mass_4l \
    (Electron_pt, Electron_eta, Electron_phi, Electron_mass, Electron_charge)")
    rdf_cut = filter_z_mass(rdf_z_mass)


    return rdf_cut


def raw_to_data_selected(fast_active=True , local_active=True):
    '''This is the selection function'''

    if fast_active == False:

        if local_active == False:
            rdf = ROOT.RDataFrame("Events", "dati.root")
        else:
            data_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
            rdf = ROOT.RDataFrame("Events", (data_path + f for f in ["Run2012B_DoubleMuParked.root",\
            "Run2012C_DoubleMuParked.root", "Run2012B_DoubleElectron.root", \
            "Run2012C_DoubleElectron.root"]))
        
        

        two_electrons_two_muons_rdf = selection_2e2mu(rdf)
        four_electrons_rdf = selection_4e(rdf)
        four_muons_rdf = selection_4mu(rdf)

        mass_2e2mu = two_electrons_two_muons_rdf.Define("InvariantMass", "invariant_mass_2el2mu \
        (Electron_pt, Electron_eta, Electron_phi, Electron_mass, Muon_pt, Muon_eta, Muon_phi, \
        Muon_mass)").Filter("InvariantMass > 70", "Mass of four leptons greater than 70 GeV")

        mass_4mu = four_muons_rdf.Define("InvariantMass", "invariant_mass_4l(Muon_pt, Muon_eta, \
        Muon_phi, Muon_mass)").Filter("InvariantMass > 70", "M of four leptons greater than 70 GeV")

        mass_4e = four_electrons_rdf.Define("InvariantMass", "invariant_mass_4l(Electron_pt,\
        Electron_eta, Electron_phi, Electron_mass)").Filter("InvariantMass > 70",\
        "Mass of four leptons greater than 70 GeV")

        mass_2e2mu.Snapshot("Events", "data/2e2muToZZ.root")
        mass_4mu.Snapshot("Events", "data/4muToZZ.root")
        mass_4e.Snapshot("Events", "data/4eToZZ.root")

def montecarlo_selection(fast_active=True , local_active=True):
    '''This is the selection function'''

    if fast_active == False:

        #Background
        if local_active == False:
            montecarlo_2e2mu_rdf = ROOT.RDataFrame("Events", "monte_2e2mu.root")
            montecarlo_4mu_rdf = ROOT.RDataFrame("Events", "monte_4mu.root")
            montecarlo_4e_rdf = ROOT.RDataFrame("Events", "monte_4e.root")
        else:
            montecarlo_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
            montecarlo_2e2mu_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo2e2mu.root")
            montecarlo_4mu_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo4mu.root")
            montecarlo_4e_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo4e.root")


        # Weights
        luminosity = 11580.0  # Integrated luminosity of the data samples
        scale_zz_4l = 1.386  #Scale factor for ZZ to four leptons

        cross_sec_zz_2e2mu = 0.18  #Standard Model cross-section for ZZ to 2e2mu
        number_event_zz_2e2mu = montecarlo_2e2mu_rdf.Count().GetValue()#Num of simul evs ZZ to 2e2mu
        weight_2e2mu = luminosity * cross_sec_zz_2e2mu * scale_zz_4l / number_event_zz_2e2mu

        cross_sec_zz_4mu = 0.077  #Standard Model cross-section for ZZ to 4mu
        number_event_zz_4mu = montecarlo_4mu_rdf.Count().GetValue() #Num of simul evs ZZ to 4mu
        weight_4mu = luminosity * cross_sec_zz_4mu * scale_zz_4l / number_event_zz_4mu

        cross_sec_zz_4e = 0.077  #Standard Model cross-section for ZZ to 4e
        number_event_zz_4e = montecarlo_4e_rdf.Count().GetValue() #Num of simul evs ZZ to 4e
        weight_4e = luminosity * cross_sec_zz_4e * scale_zz_4l / number_event_zz_4e

        montecarlo_2e2mu_rdf = selection_2e2mu(montecarlo_2e2mu_rdf)
        montecarlo_4e_rdf = selection_4e(montecarlo_4e_rdf)
        montecarlo_4mu_rdf = selection_4mu(montecarlo_4mu_rdf)

        montecarlo_mass_2e2mu = montecarlo_2e2mu_rdf.Define("InvariantMass", "invariant_mass_2el2mu\
        (Electron_pt, Electron_eta, Electron_phi,Electron_mass, Muon_pt, Muon_eta, Muon_phi, \
        Muon_mass)").Filter("InvariantMass > 70", "Mass of four leptons greater than 70 GeV")\
        .Define("weight", "{}".format(weight_2e2mu))

        montecarlo_mass_4mu = montecarlo_4mu_rdf.Define("InvariantMass", "invariant_mass_4l(Muon_pt,\
        Muon_eta, Muon_phi, Muon_mass)").Filter("InvariantMass > 70", "M of four leptons greater \
        than 70 GeV").Define("weight", "{}".format(weight_4mu))

        montecarlo_mass_4e = montecarlo_4e_rdf.Define("InvariantMass", "invariant_mass_4l(\
        Electron_pt, Electron_eta, Electron_phi,Electron_mass)").Filter("InvariantMass > 70",\
        "Mass of four leptons greater than 70 GeV").Define("weight", "{}".format(weight_4e))

        montecarlo_mass_2e2mu.Snapshot("Events", "montecarlo/montecarlo_2e2muToZZ.root")
        montecarlo_mass_4mu.Snapshot("Events", "montecarlo/montecarlo_4muToZZ.root")
        montecarlo_mass_4e.Snapshot("Events", "montecarlo/montecarlo_4eToZZ.root")

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

