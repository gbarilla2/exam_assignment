import argparse
import os

import  ROOT



ROOT.gInterpreter.Declare('#include "include/prova.h"')


# Enable multi-threading
ROOT.gInterpreter.ProcessLine("ROOT::EnableImplicitMT()")


#Python functions

def filter_Z_mass(rdf):
    '''This function select good Z mass events'''

    rdf_Z_a_mass = rdf.Filter("Z_mass[0] > 40 && Z_mass[0] < 120", "Mass of first Z candidate in [40, 120]")
    rdf_Z_b_mass = rdf_Z_a_mass.Filter("Z_mass[1] > 12 && Z_mass[1] < 120", "Mass of second Z candidate in [12, 120]")


    return rdf_Z_b_mass




def selection_2e2mu(rdf):
    '''This function select events with e+ e- mu+ mu-'''
    
    rdf_2e2mu = rdf.Filter("nElectron >= 2 && nMuon >= 2", "First selection with at least 2 electrons and 2 muons")
    
    rdf_pt = rdf_2e2mu.Filter("All(Muon_pt>5) && All(Electron_pt>7)", "Pt cuts")
    
    rdf_eta = rdf_pt.Filter("All(abs(Electron_eta)<2.5) && All(abs(Muon_eta)<2.4)", "Direction cuts (geometrical acceptance)")
    
    rdf_iso = rdf_eta.Filter("All(abs(Electron_pfRelIso03_all)<0.40) && All(abs(Muon_pfRelIso04_all)<0.40)", "Require good isolation")
    
    #Definition of 3d impact parameter significance for Electrons
    rdf_el_3dip = rdf_iso.Define("Electron_3DIP", "sqrt(Electron_dxy*Electron_dxy + Electron_dz*Electron_dz)")    
    rdf_el_3dips = rdf_el_3dip.Define("Electron_3DIPS", "Electron_3DIP/sqrt(Electron_dxyErr*Electron_dxyErr + Electron_dzErr*Electron_dzErr)")
    
    rdf_el_track = rdf_el_3dips.Filter("All(Electron_3DIPS<4) && All(abs(Electron_dxy)<0.5) && All(abs(Electron_dz)<1.0)", "Muons and electrons come from the same vertex")
    
    #Definition of 3d impact parameter significance for Muons                                
    rdf_mu_3dip = rdf_el_track.Define("Muon_3DIP", "sqrt(Muon_dxy*Muon_dxy + Muon_dz*Muon_dz)")
    rdf_mu_3dips = rdf_mu_3dip.Define("Muon_3DIPS", "Muon_3DIP/sqrt(Muon_dxyErr*Muon_dxyErr + Muon_dzErr*Muon_dzErr)")
    
    rdf_mu_track = rdf_mu_3dips.Filter("All(Muon_3DIPS<4) && All(abs(Muon_dxy)<0.5) && All(abs(Muon_dz)<1.0)", "Muons and electrons come from the same vertex")

    rdf_charge = rdf_mu_track.Filter(" Sum(Electron_charge) == 0 && Sum(Muon_charge) == 0", "Selection for total 0 charge and 0 leptons flavour")
    
    #Definition of Z masses for filtering
    rdf_Z_mass = rdf_charge.Define("Z_mass","calculation_Z_mass_2el2mu(Electron_pt, Electron_eta, Electron_phi,"
                                " Electron_mass, Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
    rdf_cut = filter_Z_mass(rdf_Z_mass)
    
    
    return rdf_cut

def selection_4mu(rdf):
    '''This function select events with 2mu+ 2mu-'''
    
    rdf_4mu = rdf.Filter("nMuon >= 4", "First selection with at least 4 muons")
    
    rdf_kin = rdf_4mu.Filter("All(Muon_pt>5) && All(abs(Muon_eta)<2.4)", "Kinematics cuts")
    
    rdf_iso = rdf_kin.Filter("All(abs(Muon_pfRelIso04_all)<0.4)", "Require good isolation")
    
    #Definition of 3d impact parameter significance 
    rdf_3dip = rdf_iso.Define("Muon_3DIP", "sqrt(Muon_dxy*Muon_dxy + Muon_dz*Muon_dz)")
    rdf_3dips = rdf_3dip.Define("Muon_3DIPS", "Muon_3DIP/sqrt(Muon_dxyErr*Muon_dxyErr + Muon_dzErr*Muon_dzErr)")
    
    rdf_pv = rdf_3dips.Filter("All(Muon_3DIPS<4) && All(abs(Muon_dxy)<0.5) && All(abs(Muon_dz)<1.0)", "Muons come from the same vertex")
    
    rdf_charge = rdf_pv.Filter("nMuon==4 && Sum(Muon_charge == 1) == 2 && Sum(Muon_charge == -1) == 2", "Selection for total 0 charge (two positive and two negative muons)")
    
    #Definition of Z masses for filtering
    rdf_Z_mass = rdf_charge.Define("Z_mass","calculation_Z_mass_4l(Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge)")
    rdf_cut = filter_Z_mass(rdf_Z_mass)
    
    
    return rdf_cut

def selection_4e(rdf):
    '''This function select events with 2e+ 2e-'''
    
    rdf_4e = rdf.Filter("nElectron >= 4", "First selection with at least 4 electrons")
    
    rdf_kin = rdf_4e.Filter("All(Electron_pt>7) && All(abs(Electron_eta)<2.5)", "Kinematics cuts")
    
    rdf_iso = rdf_kin.Filter("All(abs(Electron_pfRelIso03_all)<0.40)", "Require good isolation")
    
    #Definition of 3d impact parameter significance for Electrons
    rdf_3dip = rdf_iso.Define("Electron_3DIP", "sqrt(Electron_dxy*Electron_dxy + Electron_dz*Electron_dz)")    
    rdf_3dips = rdf_3dip.Define("Electron_3DIPS", "Electron_3DIP/sqrt(Electron_dxyErr*Electron_dxyErr + Electron_dzErr*Electron_dzErr)")
    
    rdf_pv = rdf_3dips.Filter("All(Electron_3DIPS<4) && All(abs(Electron_dxy)<0.5) && All(abs(Electron_dz)<1.0)", "Electrons come from the same vertex")
    
    rdf_charge = rdf_pv.Filter("nElectron==4 && Sum(Electron_charge == 1) == 2 && Sum(Electron_charge == -1) == 2", "Selection for total 0 charge (two positive and two negative electrons)")
    
    #Definition of Z masses for filtering
    rdf_Z_mass = rdf_charge.Define("Z_mass","calculation_Z_mass_4l(Electron_pt, Electron_eta, Electron_phi, Electron_mass, Electron_charge)")
    rdf_cut = filter_Z_mass(rdf_Z_mass)
    
    
    return rdf_cut
    

def raw_to_data_selected( fast_active = True):
    if ( fast_active == False):
        
        data_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
        rdf = ROOT.RDataFrame("Events",(data_path + f for f in ["Run2012B_DoubleMuParked.root", "Run2012C_DoubleMuParked.root","Run2012B_DoubleElectron.root", "Run2012C_DoubleElectron.root"]))


        TwoElectronsTwoMuons_rdf = selection_2e2mu(rdf)
        FourElectrons_rdf = selection_4e(rdf)
        FourMuons_rdf = selection_4mu(rdf)

        mass_2e2mu = TwoElectronsTwoMuons_rdf.Define("InvariantMass","invariant_mass_2el2mu(Electron_pt, Electron_eta, Electron_phi,"
                                        " Electron_mass, Muon_pt, Muon_eta, Muon_phi, Muon_mass)").Filter("InvariantMass > 70", "Mass of four leptons greater than 70 GeV")
                                        
        mass_4mu = FourMuons_rdf.Define("InvariantMass","invariant_mass_4l(Muon_pt, Muon_eta, Muon_phi, Muon_mass)").Filter("InvariantMass > 70", "Mass of four leptons greater than 70 GeV")

        mass_4e = FourElectrons_rdf.Define("InvariantMass","invariant_mass_4l(Electron_pt, Electron_eta, Electron_phi,Electron_mass)").Filter("InvariantMass > 70", "Mass of four leptons greater than 70 GeV")

        mass_2e2mu.Snapshot("Events","data/2e2muToZZ.root")
        mass_4mu.Snapshot("Events","data/4muToZZ.root")
        mass_4e.Snapshot("Events","data/4eToZZ.root")
            



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Program that find the Higgs boson in the decay channel H->ZZ->4l with CMS Open data.')
    parser.add_argument('--nofast', action='store_const', default=True, const=False, help="No-fast mode take data from raw file and it does the data selection. This saves also data selected in 3 files to data path.")
    args = parser.parse_args()
    
    
    raw_to_data_selected(args.nofast)




'''Spectrum= mass_2e2mu.Histo1D(("Spectrum","",36,70,180),"InvariantMass")
Spectrum4mu= mass_4mu.Histo1D(("Spectrum","",36,70,180),"InvariantMass")
Spectrum4e= mass_4e.Histo1D(("Spectrum","",36,70,180),"InvariantMass")

h = Spectrum.GetValue()
h.Add(Spectrum4mu.GetValue())
h.Add(Spectrum4e.GetValue())

h.DrawCopy("PE1")'''



