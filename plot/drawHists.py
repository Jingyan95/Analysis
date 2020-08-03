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

def cutFlowTable(hists, samples, regions, ch, year,caption='2017', nsig=6):
    mcSum = list(0 for i in xrange(0,len(regions)))
    if not ch==1:
        showData = True
    else:
        showData = False
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table*}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{' + caption +"}\n"
    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{|l|' + ''.join([('' if ('Metl' in c) or ('geqX' in c) or ('leq' in c) or ('lllB' in c) else 'l|') for c in regions]).strip() +'}' + "\n"
    table += '\\hline' + "\n"
    table += 'Samples ' + ''.join([('' if ('Metl' in c) or ('geqX' in c) or ('leq' in c) or ('lllB' in c) else ' & '+c) for c in regions]).strip() + '\\\\' + "\n"
    table += '\\hline' + "\n"

    for ids, s in enumerate(samples):
        if ids==0:
            continue
        for idr, r in enumerate(regions):
            if ids<nsig:
               mcSum[idr] += hists[year][ids][ch][idr][2].Integral()

    for ids, s in enumerate(samples):
        if ids==0:
            continue
        table += s
        for idr, r in enumerate(regions):
            if ('Metl' not in r) and ('geqX' not in r) and ('leq' not in r) and ('lllB' not in r):
               if ids<nsig:
                  if mcSum[idr]==0:
                     mcSum[idr]=1
                  table += (' & ' + str(round(hists[year][ids][ch][idr][2].Integral(),2)) + '[' + str(round((100*hists[year][ids][ch][idr][2].Integral()/mcSum[idr]),2)) + '\%]')
               else:
                  table += (' & ' + str(round(hists[year][ids][ch][idr][2].Integral(),2)) )
        table += '\\\\' + "\n"    
    table += '\\hline' + "\n"
    table += 'Prediction '
    for idr, r in enumerate(mcSum):
        if ('Metl' not in regions[idr]) and ('geqX' not in regions[idr]) and ('leq' not in regions[idr]) and ('lllB' not in regions[idr]):
           table += (' & ' + str(round(r,2)))
    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    if showData:
      table += 'Data '
      for idr, r in enumerate(regions):
          if ('Metl' not in regions[idr]) and ('geqX' not in regions[idr]) and ('leq' not in regions[idr]) and ('lllB' not in regions[idr]):
             table += (' & ' + str(hists[year][0][ch][idr][2].Integral()))
      table += '\\\\' + "\n"
      table += '\\hline' + "\n"
      table += 'Data$/$Pred. '
      for idr, r in enumerate(mcSum):
          if r==0:
             r=0.1
          if ('Metl' not in regions[idr]) and ('geqX' not in regions[idr]) and ('leq' not in regions[idr]) and ('lllB' not in regions[idr]):
             table += (' & ' + str(round(hists[year][0][ch][idr][2].Integral()/r,2)))
      table += '\\\\' + "\n"
      table += '\\hline' + "\n"
    table += '\\end{tabular}}' + "\n"
    table += '\\end{table*}' + "\n"
#    table += '\\end{sidewaystable*}' + "\n"
    print table

def stackPlots(hists, SignalHists, Fnames, ch = "channel", reg = "region", year='2017', var="sample", varname="v"):
    if not os.path.exists(year):
       os.makedirs(year)
    if not os.path.exists(year + '/' + ch):
       os.makedirs(year + '/' + ch)
    if not os.path.exists(year + '/' + ch +'/'+reg):
       os.makedirs(year + '/' + ch +'/'+reg)
    hs = ROOT.THStack("hs","")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    if ('emul' not in ch) or ('B' not in reg):
        showData = True
    else:
        showData = False
        dummy.Reset("ICE")

    
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
    Label_cms = ROOT.TLatex(0.2,0.92,label_cms)
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
    dummy_ratio.GetYaxis().SetRangeUser(0.8,1.2)
    dummy_ratio.Divide(SumofMC)
    dummy_ratio.SetStats(ROOT.kFALSE)
    dummy_ratio.GetYaxis().SetTitle('Data/Pred.')
    if not showData:
        dummy_ratio.Reset("ICE")
    dummy_ratio.Draw()
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    canvas.Print(year + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()


def stackPlotsFF(hists, SignalHists, Fnames, f="FFregion", ch = "channel", reg = "region", year='2017', var="sample", varname="v"):
    ff='FakeFactor'
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
    if ('emul' not in ch) or ('B' not in reg):
        showData = True
    else:
        showData = False
#dummy.Reset("ICE")

    
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
    Label_cms = ROOT.TLatex(0.2,0.92,label_cms)
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
    dummy_ratio.GetYaxis().SetRangeUser(0.8,1.2)
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


def compareHists(hists,Fnames, ch = "channel", reg = "region", var="sample", varname="v"):
    for num in range(len(hists)):
        if (hists[num].Integral() <= 0):
            return  
    Fol = 'compareHists'
    if not os.path.exists(Fol):
       os.makedirs(Fol)
    if not os.path.exists(Fol + '/' + ch):
       os.makedirs(Fol + '/' + ch)
    if not os.path.exists(Fol + '/' + ch +'/'+reg):
       os.makedirs(Fol + '/' + ch +'/'+reg)
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
        hists[num].Scale(1/hists[num].Integral())

    canvas = ROOT.TCanvas(ch+reg+var,ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    legend = ROOT.TLegend(0.6,0.7,0.85,0.88)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)

    pad1=ROOT.TPad("pad1", "pad1", 0.05, 0.05, 1, 0.99 , 0)#used for the hist plot
    pad1.Draw()
    pad1.cd()
    pad1.SetLogx(ROOT.kFALSE)
    pad1.SetLogy(ROOT.kFALSE)

    y_min=0
    y_max=1.2* max(hists[0].GetMaximum(), hists[1].GetMaximum(), hists[2].GetMaximum())
    hists[0].SetTitle("")
    hists[0].GetYaxis().SetTitle('Fraction')
    hists[0].GetXaxis().SetLabelSize(0.03)
    hists[0].GetXaxis().SetNoExponent()
    hists[0].GetYaxis().SetTitleOffset(0.8)
    hists[0].GetYaxis().SetTitleSize(0.05)
    hists[0].GetYaxis().SetLabelSize(0.04)
    hists[0].GetYaxis().SetRangeUser(y_min,y_max)
    hists[0].GetXaxis().SetTitle(varname)
    hists[0].Draw("Hist")
    hists[0].SetLineWidth(2)
    hists[0].SetFillColor(0)
    for H in range(1,len(hists)):
        hists[H].SetLineWidth(2)
        hists[H].SetFillColor(0)
        hists[H].Draw("histSAME")
    hists[0].Draw("AXISSAMEY+")
    hists[0].Draw("AXISSAMEX+")

    for num in range(0,len(hists)):
        legend.AddEntry(hists[num],Fnames[num],'L')
    legend.Draw("same")

    pad1.Update()
    canvas.Print(Fol + '/' + ch +'/'+reg+'/'+var + ".png")
    del canvas
    gc.collect()

#year=['2016','2017','2018','All']
year=['2017']
regions=["lll","lllOnZ","lllOffZ","lllB0","lllB1", "lllBgeq2", "lllMetl20", "lllMetg20", "lllMetl20Jet1B1", "lllMetl20Jet2B1", "lllMetg20Jetgeq1B0", "lllMetg20Jet1B1", "lllMetg20Jet2B1", "lllMetg20Jetgeq3B1", "lllMetg20Jetgeq2B2"]
channels=["eee", "emul", "mumumu"];
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","lep3Pt","lep3Eta","lep3Phi",
           "LFVePt","LFVeEta","LFVePhi","LFVmuPt","LFVmuEta","LFVmuPhi","balPt","balEta","balPhi","Topmass",
           "llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx","llMZw","LFVTopmass"]
#variables=["lep1Pt"]
variablesName=["Leading lepton p_{T} [GeV]","Leading lepton #eta","Leading lepton #varphi","2nd-Leading lepton p_{T} [GeV]","2nd-Leading lepton #eta","2nd-Leading lepton #varphi","3rd-Leading lepton p_{T} [GeV]","3rd-Leading lepton #eta","3rd-Leading lepton #varphi","cLFV electron p_{T} [GeV]","cLFV electron #eta","cLFV electron #varphi","cLFV muon p_{T} [GeV]","cLFV muon #eta","cLFV muon #varphi","Bachelor lepton p_{T} [GeV]","Bachelor lepton #eta","Bachelor lepton #varphi","Standard top mass [GeV]","M(ll) [GeV]","p_{T}(ll) [GeV]","#Delta R(ll)","#Delta #varphi(ll)","Leading jet p_{T} [GeV]","Leading jet #eta","Leading jet #varphi","Number of jets","Number of b-tagged jets","MET [GeV]","#varphi(MET)","Number of vertices", "M(ll) (OSSF) [GeV]", "cLFV top mass [GeV]"]

FFregions=["lll", "lllMetl20", "lllMetg20","lllOnZ","lllOffZ","lllMetg20B1"]
FFvariables=["FakeableEPt","FakeableMuPt","TightEPt","TightMuPt","AntiEPt","AntiMuPt"]
FFvariablesName=["Fakeable Electron p_{T} [GeV]","Fakeable Muon p_{T} [GeV]","Tight Electron p_{T} [GeV]","Tight Muon p_{T} [GeV]","anti-Electron p_{T} [GeV]","anti-Muon p_{T} [GeV]"]
FF=["MR", "AR"]



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

SamplesNameFF = ['data','TTV','WZ', 'ZZ', 'Nonprompt Electron', 'Nonprompt Muon', 'ST_vector_emutu', 'TT_vector_emutu']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6]

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

#                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])

#HH1=[]
#HHsignal1=[]
#
#HHsignal1.append(Hists[0][6][2][2][0])
#HHsignal1.append(Hists[0][7][2][2][0])
#
#Nonpromp = []
#Nonpromp.append(ROOT.TFile.Open('nonprompt.root'))
#
#f_e = 'fake_e'
#f_mu = 'fake_mu'
#
#h_fake_e=Nonpromp[0].Get(f_e)
#h_fake_mu=Nonpromp[0].Get(f_mu)
#
#h_fake_e.SetFillColor(colors[4])
#h_fake_e.SetLineColor(colors[4])
#
#h_fake_mu.SetFillColor(colors[5])
#h_fake_mu.SetLineColor(colors[5])
#
#HH1.append(Hists[0][0][2][0][0])
#HH1.append(Hists[0][1][2][0][0])
#HH1.append(Hists[0][2][2][0][0])
#HH1.append(Hists[0][3][2][0][0])
#HH1.append(h_fake_e)
#HH1.append(h_fake_mu)
#
#stackPlots(HH1, HHsignal1, SamplesNameFF, channels[2], regions[0], year[0],variables[0],variablesName[0])



### Fake Factor study
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
            for numch, namech in enumerate(channels):
                l3=[]
                for numreg, namereg in enumerate(FFregions):
                    l4=[]
                    for numvar, namevar in enumerate(FFvariables):
                        h= Files[f].Get(namef + '_' + namech + '_' + namereg + '_' + namevar)
                        h.SetFillColor(colors[f])
                        h.SetLineColor(colors[f])
                        l4.append(h)
                    l3.append(l4)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    HistsFF.append(l0)

for numyear, nameyear in enumerate(year):
    for numf, namef in enumerate(FF):
        for numch, namech in enumerate(channels):
            for numreg, namereg in enumerate(FFregions):
                for numvar, namevar in enumerate(FFvariables):
                    HHFF=[]
                    HHsignalFF=[]
                    for f in range(len(Samples)):
                        if 'SMEFTfr' in Samples[f]:
                            HHsignalFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])
                        else:
                            HHFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])
                
                    stackPlotsFF(HHFF, HHsignalFF, SamplesName, namef, namech, namereg, nameyear,namevar,FFvariablesName[numvar])

for numreg, namereg in enumerate(FFregions):
    for numch, namech in enumerate(channels):
        x=""
        if (numreg == 0) and (numch == 0):
           x="("
        if (numreg == (len(FFregions)-1)) and (numch == (len(channels)-1)):
           x=")"
        ROOT.gStyle.SetPaintTextFormat("7.3f")
        c = ROOT.TCanvas("c","c",600,600)
        lt = ROOT.TLatex()
        c.Divide(2,2)
        HHFakeableE = HistsFF[0][0][0][numch][numreg][0].Clone()
        HHFakeableMu = HistsFF[0][0][0][numch][numreg][1].Clone()
        HHTightE = HistsFF[0][0][0][numch][numreg][2].Clone()
        HHTightMu = HistsFF[0][0][0][numch][numreg][3].Clone()
        for f in range(1,len(Samples)-4): # subtracted by nonpromt contribution
            HHFakeableE.Add(HistsFF[0][f][0][numch][numreg][0],-1)
            HHFakeableMu.Add(HistsFF[0][f][0][numch][numreg][1],-1)
            HHTightE.Add(HistsFF[0][f][0][numch][numreg][2],-1)
            HHTightMu.Add(HistsFF[0][f][0][numch][numreg][3],-1)
        FR_E = HHTightE.Clone()
        FR_E.SetTitle("")
        FR_E.GetXaxis().SetNoExponent()
        FR_E.SetName("FR_E_"+namech+namereg)
        FR_E.GetYaxis().SetTitle("Fake Rate")
        FR_E.GetXaxis().SetTitle("Electron p_{T} [GeV]")
        FR_E.Divide(HHTightE, HHFakeableE, 1.0, 1.0, "B") # FR=Tight/Fakeable
        c.cd(1)
        FR_E.Draw("text0")
        FR_E.SetAxisRange(0,1,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" channel "+namereg)
        FR_Mu = HHTightMu.Clone()
        FR_Mu.SetTitle("")
        FR_Mu.GetXaxis().SetNoExponent()
        FR_Mu.SetName("FR_Mu_"+namech+namereg)
        FR_Mu.GetYaxis().SetTitle("Fake Rate")
        FR_Mu.GetXaxis().SetTitle("Muon p_{T} [GeV]")
        FR_Mu.Divide(HHTightMu, HHFakeableMu, 1.0, 1.0, "B")
        c.cd(2)
        FR_Mu.Draw("text0")
        FR_Mu.SetAxisRange(0,1,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" channel "+namereg)
        HHFakeableE.Add(HHTightE,-1)  # FF=Tight/(Fakeable-Tight)
        FF_E = FR_E.Clone()
        FF_E.SetTitle("")
        FF_E.SetName("FF_E_"+namech+namereg)
        FF_E.GetYaxis().SetTitle("Fake Factor")
        FF_E.GetXaxis().SetTitle("Electron p_{T} [GeV]")
        FF_E.Divide(HHTightE, HHFakeableE, 1.0, 1.0, "B")
        c.cd(3)
        FF_E.Draw("text0")
        FF_E.SetAxisRange(0,1,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" channel "+namereg)
        HHFakeableMu.Add(HHTightMu,-1)
        FF_Mu = FR_Mu.Clone()
        FF_Mu.SetTitle("")
        FF_Mu.SetName("FF_Mu_"+namech+namereg)
        FF_Mu.GetYaxis().SetTitle("Fake Factor")
        FF_Mu.GetXaxis().SetTitle("Muon p_{T} [GeV]")
        FF_Mu.Divide(HHTightMu, HHFakeableMu, 1.0, 1.0, "B")
        c.cd(4)
        FF_Mu.Draw("text0")
        FF_Mu.SetAxisRange(0,1,"Y")
        lt.DrawLatexNDC(0.2,0.8,namech+" channel "+namereg)
        c.Print("FakeRateAndFakeFactor.pdf"+x,"pdf")


le = '\\documentclass{article}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\begin{document}' + "\n"

print le
for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        cutFlowTable(Hists, SamplesNameLatex, regions, numch, numyear, nameyear + ' ' + namech, len(Samples)-2 )
print '\\end{document}' + "\n"


for numreg, namereg in enumerate(regions):
    for numvar, namevar in enumerate(variables):
        HH=[]
        HHname=[]
        for f in range(len(Samples)):
            if 'SMEFTfr' in Samples[f] or 'TTbar' in Samples[f]:
                HH.append(Hists[0][f][1][numreg][numvar])
                HHname.append(SamplesName[f])
#        compareHists(HH,HHname, 'emul', namereg,namevar,variablesName[numvar])
