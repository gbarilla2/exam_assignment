import  ROOT

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()


#Python functions


def selection_2e2mu(rdf):
    '''This function select events with e+ e- mu+ mu-'''
    
    rdf_2e2mu = rdf.Filter("nElectron >= 2 && nMuon >= 2", "First selection with at least 2 electrons and 2 muons")
    
    rdf_pt = rdf_2e2mu.Filter("Muon_pt > 5 && Electron_pt > 7", "Pt cuts")
    
    rdf_eta = rdf_pt.Filter("abs(Muon_eta) < 2.4 && abs(Electron_eta) < 2.5", "Direction cuts (geometrical acceptance)")
    
    rdf_charge = rdf_eta.Filter("Sum(Electron_charge) == 0 && Sum(Muon_charge) == 0", "Selection for total 0 charge")
    
    return rdf_charge

def selection_4mu(rdf):
    '''This function select events with 2mu+ 2mu-'''
    
    rdf_4mu = rdf.Filter("nMuon >= 4", "First selection with at least 4 muons")
    
    rdf_pt = rdf_4mu.Filter("Muon_pt > 5", "Pt cuts")
    
    rdf_eta = rdf_pt.Filter("abs(Muon_eta) < 2.4", "Direction cuts (geometrical acceptance)")
    
    rdf_charge = rdf_eta.Filter("Sum(Muon_charge == 1) == 2 && Sum(Muon_charge == -1) == 2", "Selection for total 0 charge (two positive and two negative muons)")
    
    return rdf_charge

def selection_4e(rdf):
    '''This function select events with 2e+ 2e-'''
    
    rdf_4e = rdf.Filter("nElectron >= 4", "First selection with at least 4 electrons")
    
    rdf_pt = rdf_4e.Filter("Electron_pt > 7", "Pt cuts")
    
    rdf_eta = rdf_pt.Filter("abs(Electron_eta) < 2.5", "Direction cuts (geometrical acceptance)")
    
    rdf_charge = rdf_eta.Filter("Sum(Electron_charge == 1) == 2 && Sum(Electron_charge == -1) == 2", "Selection for total 0 charge (two positive and two negative electrons)")
    
    return rdf_charge
    
    
    
    
f= ROOT.TFile.Open("root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root")


rdf = ROOT.RDataFrame(f.Events).Range(10000)



cppcode='''
        float fourmass(float pt, float eta, float phi, float mass\
                ,float pt2, float eta2, float phi2, float mass2\
                ,float pt3, float eta3, float phi3, float mass3\
                ,float pt4, float eta4, float phi4, float mass4){
        TLorentzVector p0,p1,p2,p3;
            p0.SetPtEtaPhiM(pt,eta,phi,mass);
            p1.SetPtEtaPhiM(pt2,eta2,phi2,mass2);
            p2.SetPtEtaPhiM(pt3,eta3,phi3,mass3);
            p3.SetPtEtaPhiM(pt4,eta4,phi4,mass4);
            return (p0+p1+p2+p3).M();
         }

'''

ROOT.gInterpreter.ProcessLine(cppcode)

mass=twoOCMuonstwoOCElectrons.Define("DiMuonDiElectronMass","fourmass(Muon_pt[mu0],Muon_eta[mu0],Muon_phi[mu0],0.106,Muon_pt[mu1],Muon_eta[mu1],Muon_phi[mu1],0.106,\
                                                                      Electron_pt[e0],Electron_eta[e0],Electron_phi[e0],0.0005,Electron_pt[e1],Electron_eta[e1],Electron_phi[e1],0.0005)")

Spectrum= mass.Histo1D(("Spectrum","The dimuondielectron mass spectrum in CMS",400,0,150),"DiMuonDiElectronMass")
Spectrum.Draw()



