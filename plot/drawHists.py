import math
import gc
import sys
import ROOT
import numpy as np
import copy
import os
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1;")
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.gStyle.SetOptStat(0)
from array import array
from ROOT import TColor
from ROOT import TGaxis
from ROOT import THStack
import gc
import argparse
TGaxis.SetMaxDigits(2)


def stackPlots(hists, SignalHists, Fnames, ch = "channel", reg = "region", year='2017', var="sample", varname="v"):
    ff='FakeFactorStudy'
    mr='MR'
    if not os.path.exists(ff):
       os.makedirs(ff)
    if not os.path.exists(ff + '/' + mr):
       os.makedirs(ff + '/' + mr)
    if not os.path.exists(ff + '/' + mr + '/' + ch):
       os.makedirs(ff + '/' + mr + '/' + ch)
    if not os.path.exists(ff + '/' + mr + '/' + ch +'/'+reg):
       os.makedirs(ff + '/' + mr + '/' + ch +'/'+reg)
    hs = ROOT.THStack("hs","")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    showData = True

    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.55,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

    pad1=ROOT.TPad("pad1", "pad1", 0, 0.315, 1, 0.99 , 0)#used for the hist plot
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.305 , 0)#used for the ratio plot
    pad1.Draw()
    pad2.Draw() 
    pad2.SetGridy()
    pad2.SetTickx()
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad2.SetTopMargin(0.1)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetFillStyle(0)
    pad1.SetFillStyle(0)
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad2.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)

    y_min=0
    y_max=1.6*hs.GetStack().Last().GetMaximum()
#y_max=1.6*dummy.GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.2)
    dummy.SetTitle("")
    dummy.GetYaxis().SetTitle('Events')
    dummy.GetXaxis().SetLabelSize(0)
    dummy.GetYaxis().SetTitleOffset(0.8)
    dummy.GetYaxis().SetTitleSize(0.07)
    dummy.GetYaxis().SetLabelSize(0.04)
    dummy.GetYaxis().SetRangeUser(y_min,y_max)
    dummy.GetYaxis().SetNoExponent()
    dummy.Draw("e")
    hs.Draw("histSAME")
    for H in SignalHists:
        H.SetLineWidth(2)
        H.SetFillColor(0)
        H.SetLineStyle(9)
        H.Draw("histSAME")
    dummy.Draw("eSAME")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")

    Lumi = '137.42'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
    label_cms="CMS Preliminary"
    Label_cms = ROOT.TLatex(0.15,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+ch+" ("+reg+")")
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")

    if showData:
       legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        legend.AddEntry(hists[num],Fnames[num],'F')
    for H in range(len(SignalHists)):
        legend.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
    legend.Draw("same")

    if (showData) and (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.2,0.75,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
        Label_DM.Draw("same")

    pad1.Update()

    pad2.cd()
    SumofMC = hs.GetStack().Last()
    dummy_ratio = hists[0].Clone()
    dummy_ratio.SetTitle("")
    dummy_ratio.SetMarkerStyle(20)
    dummy_ratio.SetMarkerSize(1.2)
    dummy_ratio.GetXaxis().SetTitle(varname)
#    dummy_ratio.GetXaxis().CenterTitle()
    dummy_ratio.GetYaxis().CenterTitle()
    dummy_ratio.GetXaxis().SetMoreLogLabels()
    dummy_ratio.GetXaxis().SetNoExponent()
    dummy_ratio.GetYaxis().SetNoExponent()
    dummy_ratio.GetXaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetYaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetXaxis().SetTitleFont(42)
    dummy_ratio.GetYaxis().SetTitleFont(42)
    dummy_ratio.GetXaxis().SetTickLength(0.05)
    dummy_ratio.GetYaxis().SetTickLength(0.05)
    dummy_ratio.GetXaxis().SetLabelSize(0.115)
    dummy_ratio.GetYaxis().SetLabelSize(0.089)
    dummy_ratio.GetXaxis().SetLabelOffset(0.02)
    dummy_ratio.GetYaxis().SetLabelOffset(0.01)
    dummy_ratio.GetYaxis().SetTitleOffset(0.42)
    dummy_ratio.GetXaxis().SetTitleOffset(1.1)
    dummy_ratio.GetYaxis().SetNdivisions(504)    
    dummy_ratio.GetYaxis().SetRangeUser(0,2)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    if not showData:
        dummy_ratio.Reset("ICE")
    dummy_ratio.Draw()
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    canvas.Print(ff + '/' + mr + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()


def stackPlotsFF(hists, SignalHists, Fnames, f="FFregion", ch = "channel", reg = "region", year='2017', var="sample", varname="v"):
    ff='FakeFactorStudy'
    if not os.path.exists(ff):
       os.makedirs(ff)
    if not os.path.exists(ff + '/' + f):
       os.makedirs(ff + '/' + f)
    if not os.path.exists(ff + '/' + f + '/' + ch):
       os.makedirs(ff + '/' + f + '/' + ch)
    if not os.path.exists(ff + '/' + f + '/' + ch + '/' + reg):
       os.makedirs(ff + '/' + f + '/' + ch + '/' + reg)
    hs = ROOT.THStack("hs","")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    showData = True

    
    canvas = ROOT.TCanvas(year+f+ch+reg+var,year+f+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.7,0.55,0.9,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

    pad1=ROOT.TPad("pad1", "pad1", 0, 0.315, 1, 0.99 , 0)#used for the hist plot
    pad2=ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.305 , 0)#used for the ratio plot
    pad1.Draw()
    pad2.Draw() 
    pad2.SetGridy()
    pad2.SetTickx()
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad2.SetTopMargin(0.1)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetFillStyle(0)
    pad1.SetFillStyle(0)
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad2.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)

    y_min=0
    y_max=1.6*hists[0].GetMaximum()
    dummy.SetMarkerStyle(20)
    dummy.SetMarkerSize(1.2)
    dummy.SetTitle("")
    dummy.GetYaxis().SetTitle('Events')
    dummy.GetXaxis().SetLabelSize(0)
    dummy.GetYaxis().SetTitleOffset(0.8)
    dummy.GetYaxis().SetTitleSize(0.07)
    dummy.GetYaxis().SetLabelSize(0.04)
    dummy.GetYaxis().SetRangeUser(y_min,y_max)
    dummy.GetYaxis().SetNoExponent()
    dummy.Draw("e")
    hs.Draw("histSAME")
    for H in SignalHists:
        H.SetLineWidth(2)
        H.SetFillColor(0)
        H.SetLineStyle(9)
        H.Draw("histSAME")
    dummy.Draw("eSAME")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")

    Lumi = '137.42'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
    label_cms="CMS Preliminary"
    Label_cms = ROOT.TLatex(0.15,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+f+" / "+ch+" ("+reg+")")
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")

    if showData:
       legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        legend.AddEntry(hists[num],Fnames[num],'F')
    for H in range(len(SignalHists)):
        legend.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
    legend.Draw("same")

    if (showData) and (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.2,0.75,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
        Label_DM.SetNDC()
        Label_DM.SetTextFont(42)
        Label_DM.Draw("same")

    pad1.Update()

    pad2.cd()
    SumofMC = hs.GetStack().Last()
    dummy_ratio = hists[0].Clone()
    dummy_ratio.SetTitle("")
    dummy_ratio.SetMarkerStyle(20)
    dummy_ratio.SetMarkerSize(1.2)
    dummy_ratio.GetXaxis().SetTitle(varname)
#    dummy_ratio.GetXaxis().CenterTitle()
    dummy_ratio.GetYaxis().CenterTitle()
    dummy_ratio.GetXaxis().SetMoreLogLabels()
    dummy_ratio.GetXaxis().SetNoExponent()
    dummy_ratio.GetYaxis().SetNoExponent()
    dummy_ratio.GetXaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetYaxis().SetTitleSize(0.04/0.3)
    dummy_ratio.GetXaxis().SetTitleFont(42)
    dummy_ratio.GetYaxis().SetTitleFont(42)
    dummy_ratio.GetXaxis().SetTickLength(0.05)
    dummy_ratio.GetYaxis().SetTickLength(0.05)
    dummy_ratio.GetXaxis().SetLabelSize(0.115)
    dummy_ratio.GetYaxis().SetLabelSize(0.089)
    dummy_ratio.GetXaxis().SetLabelOffset(0.02)
    dummy_ratio.GetYaxis().SetLabelOffset(0.01)
    dummy_ratio.GetYaxis().SetTitleOffset(0.42)
    dummy_ratio.GetXaxis().SetTitleOffset(1.1)
    dummy_ratio.GetYaxis().SetNdivisions(504)    
    dummy_ratio.GetYaxis().SetRangeUser(0,2)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    if not showData:
        dummy_ratio.Reset("ICE")
    dummy_ratio.Draw()
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    canvas.Print(ff + '/' + f + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()

year=['2017']
regions=["ll","llMetl30","llMetl30Jetg0","llJetg1Bl2"]
channels=["MR_e","MR_mu"];
variables=["FlepPt","FlepEta","FlepPhi","TlepPt","TlepEta","TlepPhi"]
variablesName=["Fakeable lepton p_{T} [GeV]", "Fakeable lepton #eta", "Fakeable lepton #phi", "Tight lepton p_{T} [GeV]", "Tight lepton #eta", "Tight lepton #phi"]


FF=["CR","AR1","AR2"]
channelsFF=["eee","mumumu"];
regionsFF=["lll","lllMetl30","lllMetg30"]
variablesFF=["lep1Pt","lep1Eta","jet1Pt","jet1Eta","njet","nbjet","Met","nVtx"]
variablesNameFF=["Leading lepton p_{T} [GeV]", "Leading lepton #eta", "Leading jet p_{T} [GeV]", "Leading jet #eta", "njet", "nbjet", "MET [GeV]", "nVtx"]


# set up an argument parser
parser = argparse.ArgumentParser()

parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--l', dest = 'LOCATION', default= '/afs/cern.ch/user/a/asparker/public/LFVTopCode_MyFork/TopLFV/hists/')

ARGS = parser.parse_args()

verbose = ARGS.VERBOSE
HistAddress = ARGS.LOCATION

Samples = ['data.root','TTV.root','WZ.root', 'ZZ.root', 'TTbar.root', 'others.root', 'SMEFTfr_ST_vector_emutu.root', 'SMEFTfr_TT_vector_emutu.root']
SamplesName = ['data','TTV','WZ', 'ZZ', 'TTbar', 'others', 'ST_vector_emutu', 'TT_vector_emutu']
SamplesNameLatex = ['data', 'TTV', 'WZ', 'ZZ', 'TTbar', 'others', 'ST\_vector\_emutu', 'TT\_vector\_emutu']

SamplesNameFF = ['data','TTV','WZ', 'ZZ', 'Nonprompt Lepton', 'ST_vector_emutu', 'TT_vector_emutu']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kPink, ROOT.kCyan-6]

## MR
Hists = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channels):
            l2=[]
            for numreg, namereg in enumerate(regions):
                l3=[]
                for numvar, namevar in enumerate(variables):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h.SetFillColor(colors[f])
                    if 'SMEFTfr' not in Samples[f]:
                        h.SetLineColor(colors[0])
                    else:
                        h.SetLineColor(colors[f])
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Hists.append(l0)       

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numvar, namevar in enumerate(variables):
                HH=[]
                HHsignal=[]
                for f in range(len(Samples)):
                    if 'SMEFTfr' in Samples[f]:
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        HH.append(Hists[numyear][f][numch][numreg][numvar])
                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])

## CR, AR1 and AR2
HistsFF = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        print HistAddress + nameyear+ '_' + Samples[f]
        for numf, namef in enumerate(FF):
            l2=[]
            for numch, namech in enumerate(channelsFF):
                l3=[]
                for numreg, namereg in enumerate(regionsFF):
                    l4=[]
                    for numvar, namevar in enumerate(variablesFF):
                        h= Files[f].Get(namef + '_' + namech + '_' + namereg + '_' + namevar)
                        h.SetFillColor(colors[f])
                        if 'SMEFTfr' not in Samples[f]:
                            h.SetLineColor(colors[0])
                        else:
                            h.SetLineColor(colors[f])
                        l4.append(h)
                    l3.append(l4)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    HistsFF.append(l0)

for numyear, nameyear in enumerate(year):
    for numf, namef in enumerate(FF):
        for numch, namech in enumerate(channelsFF):
            for numreg, namereg in enumerate(regionsFF):
                for numvar, namevar in enumerate(variablesFF):
                    HHFF=[]
                    HHsignalFF=[]
                    for f in range(len(Samples)):
                        if 'SMEFTfr' in Samples[f]:
                            HHsignalFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])
                        else:
                            HHFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])

                    stackPlotsFF(HHFF, HHsignalFF, SamplesName, namef, namech, namereg, nameyear,namevar,variablesNameFF[numvar])

ROOT.gStyle.SetPaintTextFormat("7.3f")
## FF calculation
FF1D = []
for numreg, namereg in enumerate(regions):
    l0=[]
    for numch, namech in enumerate(channels):
        x=""
        if (numreg == 0) and (numch == 0):
           x="("
        if (numreg == (len(regions)-1)) and (numch == (len(channels)-1)):
           x=")"
        c = ROOT.TCanvas("c","c",600,600)
        lt = ROOT.TLatex()
        c.Divide(2,2)
        HHFakeablePt = Hists[0][0][numch][numreg][0].Clone()
        HHFakeableEta = Hists[0][0][numch][numreg][1].Clone()
        HHTightPt = Hists[0][0][numch][numreg][3].Clone()
        HHTightEta = Hists[0][0][numch][numreg][4].Clone()
        for f in range(1,len(Samples)-4): # subtracted by promt contribution
            HHFakeablePt.Add(Hists[0][f][numch][numreg][0],-1)
            HHFakeableEta.Add(Hists[0][f][numch][numreg][1],-1)
            HHTightPt.Add(Hists[0][f][numch][numreg][3],-1)
            HHTightEta.Add(Hists[0][f][numch][numreg][4],-1)
        FR_Pt = HHFakeablePt.Clone()
        FR_Pt.SetTitle("")
        FR_Pt.GetXaxis().SetNoExponent()
        FR_Pt.SetName("FR_Pt_"+namech+namereg)
        FR_Pt.GetYaxis().SetTitle("Fake Rate")
        FR_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
        FR_Pt.Divide(HHTightPt, HHFakeablePt, 1.0, 1.0, "B") # FR=Tight/Fakeable
        c.cd(1)
        FR_Pt.Draw("text0")
        FR_Pt.SetAxisRange(0,2.5,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" / "+namereg)
        FR_Eta = HHFakeableEta.Clone()
        FR_Eta.SetTitle("")
        FR_Eta.GetXaxis().SetNoExponent()
        FR_Eta.SetName("FR_Eta_"+namech+namereg)
        FR_Eta.GetYaxis().SetTitle("Fake Rate")
        FR_Eta.GetXaxis().SetTitle("Lepton #eta")
        FR_Eta.Divide(HHTightEta, HHFakeableEta, 1.0, 1.0, "B")
        c.cd(2)
        FR_Eta.Draw("text0")
        FR_Eta.SetAxisRange(0,2.5,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" / "+namereg)
        HHFakeablePt.Add(HHTightPt,-1)  # FF=Tight/(Fakeable-Tight)
        FF_Pt = FR_Pt.Clone()
        FF_Pt.SetTitle("")
        FF_Pt.SetName("FF_Pt_"+namech+namereg)
        FF_Pt.GetYaxis().SetTitle("Fake Factor")
        FF_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
        FF_Pt.Divide(HHTightPt, HHFakeablePt, 1.0, 1.0, "B")
        c.cd(3)
        FF_Pt.Draw("text0")
        FF_Pt.SetAxisRange(0,2.5,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" / "+namereg)
        HHFakeableEta.Add(HHTightEta,-1)
        FF_Eta = FR_Eta.Clone()
        FF_Eta.SetTitle("")
        FF_Eta.SetName("FF_Eta_"+namech+namereg)
        FF_Eta.GetYaxis().SetTitle("Fake Factor")
        FF_Eta.GetXaxis().SetTitle("Lepton #eta")
        FF_Eta.Divide(HHTightEta, HHFakeableEta, 1.0, 1.0, "B")
        c.cd(4)
        FF_Eta.Draw("text0")
        FF_Eta.SetAxisRange(0,2.5,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" / "+namereg)
        c.Print("FakeRateAndFakeFactor.pdf"+x,"pdf")
        l0.append(FF_Pt)
    FF1D.append(l0)

d = ROOT.TCanvas("d","d",800,600)
ROOT.gStyle.SetPaintTextFormat("7.0f")
# Prepare 2D and 3D AR histograms
HistsAR1 = []
HistsAR2 = []
for numyear, nameyear in enumerate(year):
    l0AR1=[]
    l0AR2=[]
    Files = []
    for f in range(len(Samples)):
        l1AR1=[]
        l1AR2=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channelsFF):
            l2AR1=[]
            l2AR2=[]
            for numreg, namereg in enumerate(regionsFF):
                l3AR1=[]
                l3AR2=[]
                for numvar, namevar in enumerate(variablesFF):
                    hAR1= Files[f].Get('AR1_2D_' + namech + '_' + namereg + '_' + namevar)
                    hAR2= Files[f].Get('AR2_3D_' + namech + '_' + namereg + '_' + namevar)
                    l3AR1.append(hAR1)
                    l3AR2.append(hAR2)
                l2AR1.append(l3AR1)
                l2AR2.append(l3AR2)
            l1AR1.append(l2AR1)
            l1AR2.append(l2AR2)
        l0AR1.append(l1AR1)
        l0AR2.append(l1AR2)
    HistsAR1.append(l0AR1)
    HistsAR2.append(l0AR2)

HistsAR1[0][0][1][1][0].Draw("text0")
d.Print("test.pdf","pdf")

#Prepare 2D and 3D FF
FF2D=[]
FF3D=[]
for numreg, namereg in enumerate(regions):
    l02D=[]
    l03D=[]
    for numch, namech in enumerate(channels):
        l12D=[]
        l13D=[]
        tmp_FF=FF1D[numreg][numch].Clone()
        for numvar, namevar in enumerate(variablesFF):
            tmp_FF2D=HistsAR1[0][0][0][0][numvar].Clone()
            tmp_FF2D.Reset("ICE")
            tmp_FF3D=HistsAR2[0][0][0][0][numvar].Clone()
            tmp_FF3D.Reset("ICE")
            for ithbin in range (1,tmp_FF2D.GetNbinsX()+1):
                for jthbin in range (1,tmp_FF2D.GetNbinsY()+1):
                    tmp_FF2D.SetBinContent(ithbin,jthbin,tmp_FF.GetBinContent(jthbin))
                    tmp_FF2D.SetBinError(ithbin,jthbin,tmp_FF.GetBinError(jthbin))
                    for kthbin in range (1,tmp_FF3D.GetNbinsZ()+1):
                        tmp_FF3D.SetBinContent(ithbin,jthbin,kthbin,tmp_FF.GetBinContent(jthbin)*tmp_FF.GetBinContent(kthbin))
                        tmp_FF3D.SetBinError(ithbin,jthbin,kthbin,math.sqrt(tmp_FF.GetBinError(jthbin)**2+tmp_FF.GetBinError(kthbin)**2))
            l12D.append(tmp_FF2D)
            l13D.append(tmp_FF3D)
        l02D.append(l12D)
        l03D.append(l13D)
    FF2D.append(l02D)
    FF3D.append(l03D)

# Start making stack plots
MRreg = 3 # Pick your favorite FF MR reg
for numch, namech in enumerate(channelsFF):
    for numreg, namereg in enumerate(regionsFF):
        for numvar, namevar in enumerate(variablesFF):
            tmp_AR1=HistsAR1[0][0][numch][numreg][numvar].Clone()
            tmp_AR2=HistsAR2[0][0][numch][numreg][numvar].Clone()
            for f in range(1,len(Samples)-4): # subtracted by promt contribution
                tmp_AR1.Add(HistsAR1[0][f][numch][numreg][numvar],-1)
                tmp_AR2.Add(HistsAR2[0][f][numch][numreg][numvar],-1)
            tmp_CR1=tmp_AR1.Clone()
            tmp_CR1.Multiply(tmp_AR1, FF2D[MRreg][numch][numvar], 1.0, 1.0, "B")
            tmp_CR2=tmp_AR2.Clone()
            tmp_CR2.Multiply(tmp_AR2, FF3D[MRreg][numch][numvar], 1.0, 1.0, "B")
            tmp_CR=HistsFF[0][0][0][0][0][numvar].Clone()
            tmp_CR.Reset("ICE")
            for ithbin in range (1,tmp_CR1.GetNbinsX()+1):
                EventCounts=0#Nonprompt Contribution in CR
                Error=0
                for jthbin in range (1,tmp_CR1.GetNbinsY()+1):
                    EventCounts=EventCounts+tmp_CR1.GetBinContent(ithbin,jthbin)
                    Error=Error+(tmp_CR1.GetBinError(ithbin,jthbin))**2
                    for kthbin in range (1,tmp_CR2.GetNbinsZ()+1):
                        EventCounts=EventCounts-tmp_CR2.GetBinContent(ithbin,jthbin,kthbin)
                        Error=Error+(tmp_CR2.GetBinError(ithbin,jthbin,kthbin))**2
                tmp_CR.SetBinContent(ithbin,EventCounts)
                tmp_CR.SetBinError(ithbin,math.sqrt(Error))
            tmp_CR.SetFillColor(colors[5])
            tmp_CR.SetLineColor(colors[0])
            HHCorr=[]
            HHsignalCorr=[]
            for f in range(len(Samples)):
                if ('TTbar' not in Samples[f]) and ('others' not in Samples[f]):
                   if 'SMEFTfr' in Samples[f]:
                       HHsignalCorr.append(HistsFF[0][f][0][numch][numreg][numvar])
                   else:
                       HHCorr.append(HistsFF[0][f][0][numch][numreg][numvar])
            HHCorr.append(tmp_CR)
            namef='CR_DataDriven'
            stackPlotsFF(HHCorr, HHsignalCorr, SamplesNameFF, namef, namech, namereg, year[0],namevar,variablesNameFF[numvar])




