
import  ROOT

__all__ = ["filter_z_mass", "selection_2e2mu", "selection_4mu",
           "selection_4e", "raw_to_data_selected","montecarlo_selection"]


def filter_z_mass(rdf):
    '''This function select good Z mass events
    :param rdf: The ``RDataFrame`` from which to take the Events `TTree`.
    :type input_file: ROOT.RDataFrame
    :rtype: ROOT.RDataFrame
    Makes a ``ROOT.RDataFrame`` appling readme cuts to all analysis 
    by calling ``Filter`` and returns the result. 
    The cuts used are: (present on the readme file)
    
    *  Mass of the lepton pair that closest to Z mass 40 < m < 120 GeV;
    *  Mass of the other lepton pair 12 < m < 120 GeV;
    '''

    rdf_z_a_mass = rdf.Filter("Z_mass[0] > 40 && Z_mass[0] < 120",\
    "Mass of first Z candidate in [40, 120]")
    rdf_z_b_mass = rdf_z_a_mass.Filter("Z_mass[1] > 12 && Z_mass[1] < 120",\
    "Mass of second Z candidate in [12, 120]")


    return rdf_z_b_mass




def selection_2e2mu(rdf):
    '''This function select events with e+ e- mu+ mu-.
    
    :param rdf: The ``RDataFrame`` from which to take the Events `TTree`.
    :type input_file: ROOT.RDataFrame
    :rtype: ROOT.RDataFrame
    Makes a ``ROOT.RDataFrame`` appling readme cuts to all analysis 
    by calling ``Filter`` and returns the result. 
    Nothing will be run as only lazy action are requested.
    The cuts used are: (present on the readme file)
    
    *  Selection of only events with exactly two opposite charge muons and two opposite 
       charge electrons;
    *  Transverse momentum of muons > 5 GeV;
    *  Transverse momentum of electrons > 7 GeV;
    *  Pseudorapidity of muons < 2.4;
    *  Pseudorapidity of electrons < 2.5;
    * 3D impact parameter significance < 4;  
    * Longitudinal impact parameter w.r.t. primary vertex < 0.5 cm;
    * Transverse impact parameter w.r.t. primary vertex < 1 cm;
    
    To filter will be created ``Electron_3DIP``, ``Electron_3DIPS``,
    ``Muon_3DIP``, ``Muon_3DIPS`` and ``Z_mass`` columns
    '''

    rdf_2e2mu = rdf.Filter("nElectron == 2 && nMuon == 2",\
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
    '''This function select events with 2mu+ 2mu-
    
    :param rdf: The ``RDataFrame`` from which to take the Events `TTree`.
    :type input_file: ROOT.RDataFrame
    :rtype: ROOT.RDataFrame
    Makes a ``ROOT.RDataFrame`` appling readme cuts to all analysis 
    by calling ``Filter`` and returns the result. 
    Nothing will be run as only lazy action are requested.
    The cuts used are: (present on the readme file)
    
    *  Selection of only events with exactly four opposite charge muons;
    *  Transverse momentum of muons > 5 GeV;
    *  Pseudorapidity of muons < 2.4;
    * 3D impact parameter significance < 4;  
    * Longitudinal impact parameter w.r.t. primary vertex < 0.5 cm;
    * Transverse impact parameter w.r.t. primary vertex < 1 cm;
    
    To filter will be created ``Muon_3DIP``, ``Muon_3DIPS`` and ``Z_mass`` columns
    '''

    rdf_4mu = rdf.Filter("nMuon == 4", "First selection with 4 muons")

    rdf_kin = rdf_4mu.Filter("All(Muon_pt>5) && All(abs(Muon_eta)<2.4)", "Kinematics cuts")

    rdf_iso = rdf_kin.Filter("All(abs(Muon_pfRelIso04_all)<0.4)", "Require good isolation")

    #Definition of 3d impact parameter significance
    rdf_3dip = rdf_iso.Define("Muon_3DIP", "sqrt(Muon_dxy*Muon_dxy + Muon_dz*Muon_dz)")
    rdf_3dips = rdf_3dip.Define("Muon_3DIPS", \
    "Muon_3DIP/sqrt(Muon_dxyErr*Muon_dxyErr + Muon_dzErr*Muon_dzErr)")

    rdf_pv = rdf_3dips.Filter("All(Muon_3DIPS<4) && \
    All(abs(Muon_dxy)<0.5) && All(abs(Muon_dz)<1.0)", "Muons come from the same vertex")

    rdf_charge = rdf_pv.Filter("Sum(Muon_charge == 1) == 2 && \
    Sum(Muon_charge == -1)==2", "Selection for total 0 charge(two positive and two negative muons)")

    #Definition of Z masses for filtering
    rdf_z_mass = rdf_charge.Define("Z_mass", "calculation_Z_mass_4l\
    (Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge)")
    rdf_cut = filter_z_mass(rdf_z_mass)


    return rdf_cut

def selection_4e(rdf):
    '''This function select events with 2e+ 2e-
    
    :param rdf: The ``RDataFrame`` from which to take the Events `TTree`.
    :type input_file: ROOT.RDataFrame
    :rtype: ROOT.RDataFrame
    Makes a ``ROOT.RDataFrame`` appling readme cuts to all analysis 
    by calling ``Filter`` and returns the result. 
    Nothing will be run as only lazy action are requested.
    The cuts used are: (present on the readme file)
    
    *  Selection of only events with exactly four opposite charge electrons;
    *  Transverse momentum of electrons > 7 GeV;
    *  Pseudorapidity of electrons < 2.5;
    * 3D impact parameter significance < 4;  
    * Longitudinal impact parameter w.r.t. primary vertex < 0.5 cm;
    * Transverse impact parameter w.r.t. primary vertex < 1 cm;
    
    To filter will be created ``Electron_3DIP``, ``Electron_3DIPS``,
    and ``Z_mass`` columns
    '''
    
    rdf_4e = rdf.Filter("nElectron == 4", "First selection with 4 electrons")

    rdf_kin = rdf_4e.Filter("All(Electron_pt>7) && All(abs(Electron_eta)<2.5)", "Kinematics cuts")

    rdf_iso = rdf_kin.Filter("All(abs(Electron_pfRelIso03_all)<0.40)", "Require good isolation")

    #Definition of 3d impact parameter significance for Electrons
    rdf_3dip = rdf_iso.Define("Electron_3DIP", \
    "sqrt(Electron_dxy*Electron_dxy + Electron_dz*Electron_dz)")
    rdf_3dips = rdf_3dip.Define("Electron_3DIPS", \
    "Electron_3DIP/sqrt(Electron_dxyErr*Electron_dxyErr + Electron_dzErr*Electron_dzErr)")

    rdf_pv = rdf_3dips.Filter("All(Electron_3DIPS<4) && \
    All(abs(Electron_dxy)<0.5) && All(abs(Electron_dz)<1.0)", "Electrons come from the same vertex")

    rdf_charge = rdf_pv.Filter("Sum(Electron_charge == 1) == 2 && \
    Sum(Electron_charge == -1) == 2", \
    "Selection for total 0 charge (two positive and two negative electrons)")

    #Definition of Z masses for filtering
    rdf_z_mass = rdf_charge.Define("Z_mass", "calculation_Z_mass_4l \
    (Electron_pt, Electron_eta, Electron_phi, Electron_mass, Electron_charge)")
    rdf_cut = filter_z_mass(rdf_z_mass)


    return rdf_cut


def raw_to_data_selected(fast_active=True , local_active=True):
    '''This function read the data from a .root input file
    and take the the Events `TTree`.
    
    :param fast_active: Optional argument parser
    :type fast_active: bool
    :param local_active: Optional argument parser
    :type local_active: bool
    
    .. role:: bluetext
    
    If fast_active is :bluetext:`False` (disabled) the function redo the
    selection on the raw data instead if it is :bluetext:`True` skip this 
    function directly without doing anything.
    If local_active is :bluetext:`False` the function takes data from local
    instead if is True it takes directly from CMS open data.
    
    If you want to run it with local_active activated you must 
    have saved the data in the root path with the following 
    command-line in a python3 env (ROOT required):
    
    .. code-block:: python
        
        import ROOT
        
        data_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
        rdf = ROOT.RDataFrame("Events", (data_path + f for f in ["Run2012B_DoubleMuParked.root",\
        "Run2012C_DoubleMuParked.root", "Run2012B_DoubleElectron.root", \
        "Run2012C_DoubleElectron.root"]))
        rdf.Snapshot("Events", "dati.root")
        
    During this will be created "InvariantMass" column to filter 
    the invariant mass of the system by requiring it greater 
    than 70 GeV.
    
    '''

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
    '''This function read the data from a .root input file
    and take the the Events `TTree`.
    
    :param fast_active: Optional argument parser
    :type fast_active: bool
    :param local_active: Optional argument parser
    :type local_active: bool
    
    .. role:: bluetext
    
    If fast_active is :bluetext:`False` (disabled) the function redo the
    selection on the raw data instead if it is :bluetext:`True` skip this 
    function directly without doing anything.
    If local_active is :bluetext:`False` the function takes data from local
    instead if is True it takes directly from CMS open data.
    
    If you want to run it with local_active activated you must 
    have saved the data in the root path with the following 
    command-line in a python3 env (ROOT required):
    
    .. code-block:: python
        
        import ROOT
        
        montecarlo_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
        montecarlo_2e2mu_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo2e2mu.root")
        montecarlo_4mu_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo4mu.root")
        montecarlo_4e_rdf = ROOT.RDataFrame("Events", montecarlo_path + "ZZTo4e.root")
        montecarlo_2e2mu_rdf.Snapshot("Events", "monte_2e2mu.root")
        montecarlo_4mu_rdf.Snapshot("Events", "monte_4mu.root")
        montecarlo_4e_rdf.Snapshot("Events", "monte_4e.root")
        
    During this will be created "weight" and "InvariantMass" 
    columns to filter the invariant mass of the system by
    requiring it greater than 70 GeV.    
    '''

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