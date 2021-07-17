"""Svago.
"""



import  ROOT
import numpy


#Include the header file
ROOT.gInterpreter.Declare('#include "include/prova.h"')


# Enable multi-threading
ROOT.gInterpreter.ProcessLine("ROOT::EnableImplicitMT()")


#Python functions



timer = ROOT.TStopwatch()

if __name__ == "__main__":

    data_path = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/"
    rdf = ROOT.RDataFrame("Events", data_path + "Run2012C_DoubleMuParked.root")

    timer.Start()
    rdf_4e = rdf.Filter("nElectron >= 4", "First selection with at least 4 electrons")\
                .Histo1D(("Prova", "", 500, 0, 500), "Electron_pt")            
    h = rdf_4e.GetValue()
    h.Draw()
    timer.Stop()
    print(timer.RealTime())

    timer.Start()
    temp = rdf.AsNumpy(columns=["nElectron", "Electron_pt"])
    timer.Stop()
    timer.Continue()
    print(timer.RealTime())
    mask = temp["nElectron"] >= 4
    el_pt = temp["Electron_pt"][mask]
    timer.Stop()
    timer.Continue()
    print(timer.RealTime())
    h1 = ROOT.TH1D("h1","h1 title",500, 0, 500)
    for i in range (0, len(el_pt)):
        for j in range (0,len(el_pt[i])):
            h1.Fill(el_pt[i][j])
    h1.SetLineColor(2)
    h1.Draw("SAME")
    timer.Stop()
    print(timer.RealTime())
