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


def stackPlots(hists, SignalHists, Fnames, ch = "channel", reg = "region", eta = "eta", year='2017', var="sample", varname="v"):
    ff='FakeFactorStudy_emul'
    mr='MR'
    if not os.path.exists(ff):
       os.makedirs(ff)
    if not os.path.exists(ff + '/' + mr):
       os.makedirs(ff + '/' + mr)
    if not os.path.exists(ff + '/' + mr + '/' + ch):
       os.makedirs(ff + '/' + mr + '/' + ch)
    if not os.path.exists(ff + '/' + mr + '/' + ch +'/'+reg):
       os.makedirs(ff + '/' + mr + '/' + ch +'/'+reg)
    if not os.path.exists(ff + '/' + mr + '/' + ch +'/'+reg +'/'+eta):
       os.makedirs(ff + '/' + mr + '/' + ch +'/'+reg +'/'+eta)
    hs = ROOT.THStack("hs","")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) + hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) + SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for num in range(1,len(hists)):
        hs.Add(hists[num])
    
    dummy = hists[0].Clone()
    if ('emul' in ch) and ('OnZ' not in reg):
        showData = False
        dummy.Reset("ICE")
    else:
        showData = True

    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    if showData:
        legend = ROOT.TLegend(0.66,0.74,0.78,0.88)
    else:
        legend = ROOT.TLegend(0.66,0.74,0.78,0.833)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0);
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend2 = ROOT.TLegend(0.78,0.74,0.9,0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0);
    legend2.SetTextFont(42)
    legend2.SetTextSize(0.04)
    legend3 = ROOT.TLegend(0.66,0.64,0.81,0.74)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0);
    legend3.SetTextFont(42)
    legend3.SetTextSize(0.04)

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
    if showData:
       y_max=1.6*hists[0].GetMaximum()
    else:
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
    for H in range(len(SignalHists)):
        SignalHists[H].SetLineWidth(2)
        SignalHists[H].SetFillColor(0)
        SignalHists[H].SetLineStyle(9)
        if H>0:
            SignalHists[H].Scale(10)
        SignalHists[H].Draw("histSAME")
    dummy.Draw("eSAME")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")

    Lumi = '137'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
    label_cms="CMS"
    Label_cms = ROOT.TLatex(0.15,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    label_cms1="Work in Progress"
    Label_cms1 = ROOT.TLatex(0.22,0.92,label_cms1)
    Label_cms1.SetNDC()
    Label_cms1.SetTextSize(0.042);
    Label_cms1.SetTextFont(52)
    Label_cms1.Draw()
    Label_lumi = ROOT.TLatex(0.72,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+ch)
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")
    Label_region = ROOT.TLatex(0.2,0.75,reg+" / "+eta)
    Label_region.SetNDC()
    Label_region.SetTextFont(42)
    Label_region.Draw("same")

    if showData:
       legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        if num<(len(hists)-3):
           legend.AddEntry(hists[num],Fnames[num],'F')
        else:
           legend2.AddEntry(hists[num],Fnames[num],'F')
    for H in range(len(SignalHists)):
        if H==0:
            legend3.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
        else:
            legend3.AddEntry(SignalHists[H], Fnames[len(hists)+H]+" (x10)",'L')
    legend.Draw("same")
    legend2.Draw("same")
    legend3.Draw("same")

    if (showData) and (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.2,0.7,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
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
    canvas.Print(ff + '/' + mr + '/' + ch +'/'+ reg +'/'+ eta +'/'+var + ".png")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) - hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) - SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for H in range(len(SignalHists)):
        if H>0:
            SignalHists[H].Scale(0.1)
    del canvas
    gc.collect()


def stackPlotsFF(hists, SignalHists, error, errorRatio, Fnames, f="FFregion", ch = "channel", reg = "region", year='2017', var="sample", varname="v"):
    ff='FakeFactorStudy_emul'
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
    if ('emul' in ch) and ('OnZ' not in reg):
        showData = False
        dummy.Reset("ICE")
    else:
        showData = True

    
    canvas = ROOT.TCanvas(year+f+ch+reg+var,year+f+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    if showData:
        legend = ROOT.TLegend(0.64,0.74,0.76,0.88)
    else:
        legend = ROOT.TLegend(0.64,0.74,0.76,0.833)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0);
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend2 = ROOT.TLegend(0.76,0.74,0.88,0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0);
    legend2.SetTextFont(42)
    legend2.SetTextSize(0.04)
    legend3 = ROOT.TLegend(0.64,0.64,0.79,0.74)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0);
    legend3.SetTextFont(42)
    legend3.SetTextSize(0.04)

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
    if showData:
       y_max=1.6*hists[0].GetMaximum()
    else:
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
    dummy.Draw("ex0")
    hs.Draw("histSAME")
    for H in range(len(SignalHists)):
        SignalHists[H].SetLineWidth(2)
        SignalHists[H].SetFillColor(0)
        SignalHists[H].SetLineStyle(9)
        if H>0:
           SignalHists[H].Scale(10)
        SignalHists[H].Draw("histSAME")
    dummy.Draw("ex0SAME")
    dummy.Draw("AXISSAMEY+")
    dummy.Draw("AXISSAMEX+")

    error.SetFillColor(13)
    error.SetLineColor(13)
    error.SetFillStyle(3004)
    error.Draw("2")

    Lumi = '137'
    if (year == '2016'):
        Lumi = '35.92'
    if (year == '2017'):
        Lumi = '41.53'
    if (year == '2018'):
        Lumi = '59.97'
    label_cms="CMS"
    Label_cms = ROOT.TLatex(0.15,0.92,label_cms)
    Label_cms.SetNDC()
    Label_cms.SetTextFont(61)
    Label_cms.Draw()
    label_cms1="Work in Progress"
    Label_cms1 = ROOT.TLatex(0.22,0.92,label_cms1)
    Label_cms1.SetNDC()
    Label_cms1.SetTextSize(0.042);
    Label_cms1.SetTextFont(52)
    Label_cms1.Draw()
    Label_lumi = ROOT.TLatex(0.72,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+f+" / "+ch)
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")
    Label_region = ROOT.TLatex(0.2,0.75,reg)
    Label_region.SetNDC()
    Label_region.SetTextFont(42)
    Label_region.Draw("same")

    x=3

    if showData:
       legend.AddEntry(dummy,Fnames[0],'ep')
    for num in range(1,len(hists)):
        if num<(len(hists)-x):
           legend.AddEntry(hists[num],Fnames[num],'F')
        else:
           legend2.AddEntry(hists[num],Fnames[num],'F')
    for H in range(len(SignalHists)):
        if H==0:
            legend3.AddEntry(SignalHists[H], Fnames[len(hists)+H],'L')
        else:
            legend3.AddEntry(SignalHists[H], Fnames[len(hists)+H]+" (x10)",'L')
    legend.Draw("same")
    legend2.Draw("same")
    legend3.Draw("same")

    if (showData) and (hs.GetStack().Last().Integral()>0):
        Label_DM = ROOT.TLatex(0.2,0.7,"Data/MC = " + str(round(hists[0].Integral()/hs.GetStack().Last().Integral(),2)))
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
    dummy_ratio.Draw('ex0')
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    errorRatio.SetFillColor(13)
    errorRatio.SetLineColor(13)
    errorRatio.SetFillStyle(3004)
    errorRatio.Draw("2")
    canvas.Print(ff + '/' + f + '/' + ch +'/'+reg+'/'+var + ".png")
    for num in range(len(hists)):
        hists[num].SetBinContent(hists[num].GetXaxis().GetNbins(), hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()) - hists[num].GetBinContent(hists[num].GetXaxis().GetNbins()+1))
    for num in range(len(SignalHists)):
        SignalHists[num].SetBinContent(SignalHists[num].GetXaxis().GetNbins(),SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()) - SignalHists[num].GetBinContent(SignalHists[num].GetXaxis().GetNbins()+1))
    for H in range(len(SignalHists)):
        if H>0:
            SignalHists[H].Scale(0.1)
    del canvas
    gc.collect()

year=['2017']
regions=["ll","llMetg20","llMetg20Jetgeq1Bleq1","llJetgeq2Bleq1","llMetg20Bgeq1","llJetgeq1B0","llMetl20Jetgeq1B0","llMetg20B2"]
channels=["MR_e","MR_mu"];
etaregs=["barrel","transition","endcap"]
variables=["FlepPt","FlepEta","FlepPhi","TlepPt","TlepEta","TlepPhi"]
variablesName=["Fakeable lepton p_{T} [GeV]", "Fakeable lepton #eta", "Fakeable lepton #phi", "Tight lepton p_{T} [GeV]", "Tight lepton #eta", "Tight lepton #phi"]


FF=["SR","AR1","AR2"]
channelsFF=["e","mu","emu"];
regionsFF=["lll","lllOnZ"]
variablesFF=["lep1Pt","lep1Eta","jet1Pt","jet1Eta","njet","nbjet","Met","nVtx"]
variablesNameFF=["Leading lepton p_{T} [GeV]", "Leading lepton #eta", "Leading jet p_{T} [GeV]", "Leading jet #eta", "njet", "nbjet", "MET [GeV]", "nVtx"]
channelsZZ=["eeee","mumumumu"];
regionsZZ=["lll","lllMetl20","lllMetg20","lllOnZ","lllOffZ","lllOffZMetg20Jetgeq1Bleq1"]


# set up an argument parser
parser = argparse.ArgumentParser()

parser.add_argument('--v', dest='VERBOSE', default=True)

ARGS = parser.parse_args()

verbose = ARGS.VERBOSE
HistAddress = '/eos/cms/store/user/jingyan/TEST/'
HistAddress_emul = '/eos/cms/store/user/jingyan/'

Samples = ['data.root','TTV.root','WZ.root', 'ZZ.root', 'TTbar.root', 'others.root', 'LFVStVecU.root', 'LFVTtVecU.root']
SamplesName = ['data','TTV','WZ', 'ZZ', 'TTbar', 'others', 'ST_vector_emutu', 'TT_vector_emutu']
SamplesNameLatex = ['data', 'TTV', 'WZ', 'ZZ', 'TTbar', 'others', 'ST\_vector\_emutu', 'TT\_vector\_emutu']

SamplesNameFF = ['data','TTV','WZ', 'ZZ', 'Nonprompt-e', 'Nonprompt-mu', 'ST_vector_emutu', 'TT_vector_emutu']

colors =  [ROOT.kBlack,ROOT.kYellow,ROOT.kGreen,ROOT.kBlue-3,ROOT.kRed-4,ROOT.kOrange-3, ROOT.kOrange-6, ROOT.kCyan-6, ROOT.kGray, ROOT.kViolet-9, ROOT.kViolet]

MakePlots = False
ZZcorr=False
FFsys=True

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
                for numeta, nameeta in enumerate(etaregs):
                    l4=[]
                    for numvar, namevar in enumerate(variables):
                        h= Files[f].Get(namech + '_' + namereg + '_' + nameeta + '_' + namevar)
                        h.SetFillColor(colors[f])
                        if 'LFV' not in Samples[f]:
                            h.SetLineColor(colors[0])
                        else:
                            h.SetLineColor(colors[f])
                        l4.append(h)
                    l3.append(l4)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    Hists.append(l0)       

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        for numreg, namereg in enumerate(regions):
            for numeta, nameeta in enumerate(etaregs):
                for numvar, namevar in enumerate(variables):
                    HH=[]
                    HHsignal=[]
                    for f in range(len(Samples)):
                        if 'LFV' in Samples[f]:
                            HHsignal.append(Hists[numyear][f][numch][numreg][numeta][numvar])
                        else:
                            HH.append(Hists[numyear][f][numch][numreg][numeta][numvar])
                    if MakePlots:
                       stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameeta, nameyear,namevar,variablesName[numvar])

## SR
HistsFF = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress_emul + nameyear+ '_' + Samples[f]))
        print HistAddress_emul + nameyear+ '_' + Samples[f]
        for numf, namef in enumerate(FF):
            l2=[]
            if numf>0:
               continue
            for numch, namech in enumerate(channelsFF):
                l3=[]
                if numch>0:
                   continue
                for numreg, namereg in enumerate(regionsFF):
                    l4=[]
                    for numvar, namevar in enumerate(variablesFF):
                        h= Files[f].Get('emul_' + namereg + '_' + namevar)
                        h.SetFillColor(colors[f])
                        if 'LFV' not in Samples[f]:
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
        if numf>0:
            continue
        for numch, namech in enumerate(channelsFF):
            if numch==0:
                namech = 'emul'
            else:
                continue
            for numreg, namereg in enumerate(regionsFF):
                for numvar, namevar in enumerate(variablesFF):
                    HHFF=[]
                    HHsignalFF=[]
                    for f in range(len(Samples)):
                        if 'LFV' in Samples[f]:
                            HHsignalFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])
                        else:
                            HHFF.append(HistsFF[numyear][f][numf][numch][numreg][numvar])
                    Error=ROOT.TGraphAsymmErrors()
                    ErrorRatio=ROOT.TGraphAsymmErrors()
                    if MakePlots:
                       stackPlotsFF(HHFF, HHsignalFF, Error, ErrorRatio, SamplesName, namef, namech, namereg, nameyear,namevar,variablesNameFF[numvar])

## ZZ->4l CR
HistsZZ = []
for numyear, nameyear in enumerate(year):
    l0=[]
    Files = []
    for f in range(len(Samples)):
        l1=[]
        Files.append(ROOT.TFile.Open(HistAddress + nameyear+ '_' + Samples[f]))
        print HistAddress + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channelsZZ):
            l2=[]
            for numreg, namereg in enumerate(regionsZZ):
                l3=[]
                for numvar, namevar in enumerate(variablesFF):
                    h= Files[f].Get(namech + '_' + namereg + '_' + namevar)
                    h.SetFillColor(colors[f])
                    if 'LFV' not in Samples[f]:
                        h.SetLineColor(colors[0])
                    else:
                        h.SetLineColor(colors[f])
                    l3.append(h)
                l2.append(l3)
            l1.append(l2)
        l0.append(l1)
    HistsZZ.append(l0)

for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channelsZZ):
        for numreg, namereg in enumerate(regionsZZ):
            for numvar, namevar in enumerate(variablesFF):
                HHZZ=[]
                HHsignalZZ=[]
                for f in range(len(Samples)):
                    if 'LFV' in Samples[f]:
                        HHsignalZZ.append(HistsZZ[numyear][f][numch][numreg][numvar])
                    else:
                        HHZZ.append(HistsZZ[numyear][f][numch][numreg][numvar])
                Error=ROOT.TGraphAsymmErrors()
                ErrorRatio=ROOT.TGraphAsymmErrors()
                namef='CR_ZZTo4l'
                if MakePlots:
                   stackPlotsFF(HHZZ, HHsignalZZ, Error, ErrorRatio, SamplesName, namef, namech, namereg, nameyear,namevar,variablesNameFF[numvar])

ROOT.gStyle.SetPaintTextFormat("15.3f")
## TF calculation
TF = []
for numch, namech in enumerate(channelsZZ):
    l0=[]
    for numreg, namereg in enumerate(regionsFF):
        l1=[]
        for numvar, namevar in enumerate(variablesFF):
            HHZZData = HistsZZ[0][0][numch][numreg][numvar].Clone()
            HHZZMC = HistsZZ[0][3][numch][numreg][numvar].Clone()
            HHTF = HHZZData.Clone()
            HHTF.Divide(HHZZData, HHZZMC, 1.0, 1.0, "B")
            l1.append(HHTF)
        l0.append(l1)
    TF.append(l0)

## FF calculation
FF1D = []
FF1Dup = []
FF1Ddown = []
Trash = [] #trash bin to store Histos
LTrash = [] #trash bin to store Legends
ff='FakeFactorStudy_emul'
c = ROOT.TCanvas("c","c",600,800)
c.Divide(2,3)
c2 = ROOT.TCanvas("c2","c2",600,800)
c2.Divide(2,3)
c3 = ROOT.TCanvas("c3","c3",600,800)
c3.Divide(2,3)
for numreg, namereg in enumerate(regions):
    l0=[]
    l0up=[]
    l0down=[]
    x=""
    if (numreg == 0):
        x="("
    if (numreg == (len(regions)-1)):
        x=")"
    for numch, namech in enumerate(channels):
        l1=[]
        l1up=[]
        l1down=[]
        for numeta, nameeta in enumerate(etaregs):
            lt = ROOT.TLatex()
            legend = ROOT.TLegend(0.11,0.65,0.4,0.78)
            legend.SetTextFont(42)
            legend.SetTextSize(0.025)
            HHFakeablePt = Hists[0][0][numch][numreg][numeta][0].Clone()
            HHTightPt = Hists[0][0][numch][numreg][numeta][3].Clone()
            TTbarFakeablePt = Hists[0][4][numch][numreg][numeta][0].Clone()
            TTbarTightPt = Hists[0][4][numch][numreg][numeta][3].Clone()
            OthersFakeablePt = Hists[0][5][numch][numreg][numeta][0].Clone()
            OthersTightPt = Hists[0][5][numch][numreg][numeta][3].Clone()
            DataFakeablePt = Hists[0][0][numch][numreg][numeta][0].Clone()
            DataTightPt = Hists[0][0][numch][numreg][numeta][3].Clone()
            HHFakeablePtUp = Hists[0][0][numch][numreg][numeta][0].Clone()
            HHTightPtUp = Hists[0][0][numch][numreg][numeta][3].Clone()
            HHFakeablePtDown = Hists[0][0][numch][numreg][numeta][0].Clone()
            HHTightPtDown = Hists[0][0][numch][numreg][numeta][3].Clone()
            for f in range(1,len(Samples)-4): # subtracted by promt contribution
                HHFakeablePt-=Hists[0][f][numch][numreg][numeta][0]
                HHTightPt-=Hists[0][f][numch][numreg][numeta][3]
                HHFakeablePtUp-=1.2*Hists[0][f][numch][numreg][numeta][0]
                HHTightPtUp-=1.2*Hists[0][f][numch][numreg][numeta][3]
                HHFakeablePtDown-=0.8*Hists[0][f][numch][numreg][numeta][0]
                HHTightPtDown-=0.8*Hists[0][f][numch][numreg][numeta][3]
            
            FR_Pt = HHFakeablePt.Clone()
            FR_Pt.SetTitle("")
            FR_Pt.SetName("FR_Pt_"+namech+namereg+nameeta)
            FR_Pt.GetXaxis().SetNoExponent()
            FR_Pt.GetYaxis().SetNoExponent()
            FR_Pt.GetYaxis().SetTitle("Fake Rate")
            FR_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
            FR_Pt.Divide(HHTightPt, HHFakeablePt, 1.0, 1.0, "B")
            FR_Pt.SetAxisRange(0,1.1,"Y")
            
            c3.cd(2*numeta+numch+1)
            FR_Pt.Draw("lp")
            lt.DrawLatexNDC(0.12,0.8,namech+" / "+namereg+" / "+nameeta)
        
            FF_Pt = HHFakeablePt.Clone()
            FF_Pt.SetTitle("")
            FF_Pt.SetName("FF_Pt_"+namech+namereg+nameeta)
            FF_Pt.GetXaxis().SetNoExponent()
            FF_Pt.GetYaxis().SetNoExponent()
            FF_Pt.GetYaxis().SetTitle("Fake Factor")
            FF_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
            FF_Pt.Divide(HHTightPt, HHFakeablePt-HHTightPt, 1.0, 1.0, "B")
            
            FF_TTbar_Pt = TTbarFakeablePt.Clone()
            FF_TTbar_Pt.SetMarkerColor(2)
            FF_TTbar_Pt.SetLineColor(2)
            FF_TTbar_Pt.SetTitle("")
            FF_TTbar_Pt.SetName("FF_TTbar_Pt_"+namech+namereg+nameeta)
            FF_TTbar_Pt.GetXaxis().SetNoExponent()
            FF_TTbar_Pt.GetYaxis().SetNoExponent()
            FF_TTbar_Pt.GetYaxis().SetTitle("Fake Factor")
            FF_TTbar_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
            FF_TTbar_Pt.Divide(TTbarTightPt, TTbarFakeablePt-TTbarTightPt, 1.0, 1.0, "B")
            
            FF_Others_Pt = OthersFakeablePt.Clone()
            FF_Others_Pt.SetMarkerColor(3)
            FF_Others_Pt.SetLineColor(3)
            FF_Others_Pt.SetTitle("")
            FF_Others_Pt.SetName("FF_Others_Pt_"+namech+namereg+nameeta)
            FF_Others_Pt.GetXaxis().SetNoExponent()
            FF_Others_Pt.GetYaxis().SetNoExponent()
            FF_Others_Pt.GetYaxis().SetTitle("Fake Factor")
            FF_Others_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
            FF_Others_Pt.Divide(OthersTightPt, OthersFakeablePt-OthersTightPt, 1.0, 1.0, "B")

            FF_Data_Pt = DataFakeablePt.Clone()
            FF_Data_Pt.SetTitle("")
            FF_Data_Pt.SetName("FF_Data_Pt_"+namech+namereg+nameeta)
            FF_Data_Pt.GetXaxis().SetNoExponent()
            FF_Data_Pt.GetYaxis().SetNoExponent()
            FF_Data_Pt.GetYaxis().SetTitle("Fake Factor")
            FF_Data_Pt.GetXaxis().SetTitle("Lepton p_{T} [GeV]")
            FF_Data_Pt.Divide(DataTightPt, DataFakeablePt-DataTightPt, 1.0, 1.0, "B")
            
            FF_Pt_Up = HHFakeablePtUp.Clone()
            FF_Pt_Up.Divide(HHTightPtUp, HHFakeablePtUp-HHTightPtUp, 1.0, 1.0, "B")
            
            FF_Pt_Down = HHFakeablePtDown.Clone()
            FF_Pt_Down.Divide(HHTightPtDown, HHFakeablePtDown-HHTightPtDown, 1.0, 1.0, "B")
            
            FF_Diff_Pt = FF_TTbar_Pt.Clone()
            FF_Diff_Pt.GetYaxis().SetTitle("Fractional Difference")
            FF_Diff_Pt.Divide(FF_TTbar_Pt-FF_Pt, FF_Pt, 1.0, 1.0, "B")
            FF_Diff_Pt.SetAxisRange(-1.5,1.5,"Y")
            
            c2.cd(2*numeta+numch+1)
            FF_Diff_Pt.Draw("lp")
            lt.DrawLatexNDC(0.12,0.8,namech+" / "+namereg+" / "+nameeta)
            
            binwidth= array( 'd' )
            bincenter= array( 'd' )
            yvalue= array( 'd' )
            yerrup= array( 'd' )
            yerrdown= array( 'd' )
            
            FF_PtUp = FF_Pt.Clone()
            FF_PtDown = FF_Pt.Clone()
            
            FF_sys_sample = FF_Pt.Clone()
            FF_sys_sample.Reset("ICE")
            FF_sys_ek = FF_sys_sample.Clone()
            FF_sys_stat = FF_sys_sample.Clone()
       
            for b in range(FF_Pt_Up.GetNbinsX()+1):
                StatError = FF_Pt.GetBinError(b+1)
                if (numeta==0):
                    SampleError = 0.4*abs(FF_Pt.GetBinContent(b+1))
                else:
                    SampleError = 0.2*abs(FF_Pt.GetBinContent(b+1))
                binwidth.append(FF_Pt_Up.GetBinWidth(b+1)/2)
                bincenter.append(FF_Pt_Up.GetBinCenter(b+1))
                if FF_Pt_Up.GetBinContent(b+1)>0:
                   error1 = abs(FF_Pt_Up.GetBinContent(b+1)-FF_Pt.GetBinContent(b+1))
                else:
                   error1 = 0.0000001
                if FF_Pt_Down.GetBinContent(b+1)>0:
                   error2 = abs(FF_Pt_Down.GetBinContent(b+1)-FF_Pt.GetBinContent(b+1))
                else:
                   error2 = 0.0000001
                if error1>error2:
                   EwError = error2
                else:
                   EwError = error1
                FF_PtUp.SetBinContent(b+1,FF_Pt.GetBinContent(b+1)+math.sqrt(EwError**2+SampleError**2+StatError**2))
                FF_PtDown.SetBinContent(b+1,FF_Pt.GetBinContent(b+1)-math.sqrt(EwError**2+SampleError**2+StatError**2))
#                FF_sys_sample.SetBinContent(b+1,SampleError/FF_Pt.GetBinContent(b+1))
#                FF_sys_ek.SetBinContent(b+1,EwError/FF_Pt.GetBinContent(b+1))
#                FF_sys_stat.SetBinContent(b+1,StatError/FF_Pt.GetBinContent(b+1))
            
#            if (numreg==4) and (numch==1) and (numeta==0):
#               x = ROOT.TCanvas("x","x",800,600)
#               x.cd()
#               FF_sys_sample.SetLineColor(1)
#               FF_sys_sample.GetYaxis().SetTitle('Uncertainty (%)')
#               FF_sys_sample.SetFillColor(0)
#               FF_sys_ek.SetLineColor(2)
#               FF_sys_ek.SetFillColor(0)
#               FF_sys_stat.SetLineColor(3)
#               FF_sys_sample.Scale(100)
#               FF_sys_ek.Scale(100)
#               FF_sys_stat.Scale(100)
#               FF_sys_stat.SetFillColor(0)
#               FF_sys_sample.SetAxisRange(0,110,"X")
#               FF_sys_sample.Draw("hist")
#               FF_sys_ek.Draw("histsame")
#               FF_sys_stat.Draw("histsame")
#               legendx = ROOT.TLegend(0.11,0.65,0.5,0.78)
#               legendx.AddEntry(FF_sys_sample,"Sample Dependence",'l')
#               legendx.AddEntry(FF_sys_ek,"Prompt Background Crossection",'l')
#               legendx.AddEntry(FF_sys_stat,"Statistical Uncertainty",'l')
#               legendx.Draw("same")
#               x.Print("NonpromptUncertainty.pdf","pdf")
            
            c.cd(2*numeta+numch+1)
            FF_Pt.SetAxisRange(0,2,"Y")
            FF_Pt.Draw("lp")
            FF_TTbar_Pt.Draw("lp same")
            FF_Others_Pt.Draw("lp same")
            legend.AddEntry(FF_Pt,"Data",'lp')
            legend.AddEntry(FF_TTbar_Pt,"TTbar MC",'lp')
            legend.AddEntry(FF_Others_Pt,"Z+jets",'lp')
            legend.Draw("same")
            lt.DrawLatexNDC(0.12,0.8,namech+" / "+namereg+" / "+nameeta)
            l1.append(FF_Pt)
            l1up.append(FF_PtUp)
            l1down.append(FF_PtDown)
            Trash += [FF_TTbar_Pt,FF_Others_Pt,FF_Data_Pt,FF_Diff_Pt,FR_Pt]
            LTrash += [legend]
        l0.append(l1)
        l0up.append(l1up)
        l0down.append(l1down)
    if not os.path.exists(ff):
           os.makedirs(ff)
    c.Print(ff+"/FakeFactor.pdf"+x,"pdf")
    c2.Print(ff+"/FractionalDifference.pdf"+x,"pdf")
    c3.Print(ff+"/FakeRate.pdf"+x,"pdf")
    FF1D.append(l0)
    FF1Dup.append(l0up)
    FF1Ddown.append(l0down)

#ROOT.gStyle.SetPaintTextFormat("7.0f")
#d = ROOT.TCanvas("d","d",800,600)
#FF1D[3][0][0].SetAxisRange(0,2,"Y")
#FF1D[3][0][0].SetMarkerStyle(8)
#FF1D[3][0][1].SetMarkerStyle(8)
#FF1D[3][0][2].SetMarkerStyle(8)
#FF1D[3][0][0].SetMarkerColor(1)
#FF1D[3][0][1].SetMarkerColor(2)
#FF1D[3][0][2].SetMarkerColor(3)
#FF1D[3][0][0].SetLineColor(1)
#FF1D[3][0][1].SetLineColor(2)
#FF1D[3][0][2].SetLineColor(3)
#FF1D[3][0][0].Draw("lp")
#FF1D[3][0][1].Draw("lp same")
#FF1D[3][0][2].Draw("lp same")
#l = ROOT.TLegend(0.20,0.65,0.50,0.75)
#l.SetFillColor(0)
#l.SetLineColor(0)
#l.SetTextSize(0.04)
#l.AddEntry(FF1D[4][1][0],"Barrel","lep")
#l.AddEntry(FF1D[4][1][1],"Transition","lep")
#l.AddEntry(FF1D[4][1][2],"Endcap","lep")
#l.SetTextFont(42)
#l.Draw()
#d.Print("FF_Electron.pdf","pdf")

d = ROOT.TCanvas("d","d",800,600)
d.Divide(2,3)
for numch, namech in enumerate(channels):
    for numeta, nameeta in enumerate(etaregs):
        d.cd(2*numeta+numch+1)
        lt = ROOT.TLatex()
        l = ROOT.TLegend(0.20,0.65,0.50,0.75)
        l.SetFillColor(0)
        l.SetLineColor(0)
        l.SetTextSize(0.04)
        for numreg, namereg in enumerate(regions):
            if ((numch==0) and ((numreg==3) or (numreg>5))) or ((numch==1) and (numreg>4)):
                FF1D[numreg][numch][numeta].SetAxisRange(0,10,"Y")
                FF1D[numreg][numch][numeta].SetMarkerStyle(8)
                FF1D[numreg][numch][numeta].SetMarkerSize(0.3)
                if (numreg<6):
                   FF1D[numreg][numch][numeta].SetMarkerColor(1)
                   FF1D[numreg][numch][numeta].SetLineColor(1)
#                   for b in range(FF1D[numreg][numch][numeta].GetNbinsX()):
#                       content=FF1D[numreg][numch][numeta].GetBinContent(b+1)
#                       Error1=abs(FF1D[6][numch][numeta].GetBinContent(b+1)-content)
#                       Error2=abs(FF1D[7][numch][numeta].GetBinContent(b+1)-content)
#                       if Error1>Error2:
#                          FF1Dup[numreg][numch][numeta].SetBinContent(b+1,content+Error1)
#                          FF1Ddown[numreg][numch][numeta].SetBinContent(b+1,content-Error1)
#                       else:
#                          FF1Dup[numreg][numch][numeta].SetBinContent(b+1,content+Error2)
#                          FF1Ddown[numreg][numch][numeta].SetBinContent(b+1,content-Error2)
                elif (numreg==6):
                    FF1D[numreg][numch][numeta].SetMarkerColor(2)
                    FF1D[numreg][numch][numeta].SetLineColor(2)
                elif (numreg==7):
                    FF1D[numreg][numch][numeta].SetMarkerColor(3)
                    FF1D[numreg][numch][numeta].SetLineColor(3)
                if (numreg<6):
                    FF1D[numreg][numch][numeta].Draw("lp")
                else:
                    FF1D[numreg][numch][numeta].Draw("lp same")
                if (numreg<6):
                    l.AddEntry(FF1D[numreg][numch][numeta],"FF MR","lep")
                elif (numreg==6):
                    l.AddEntry(FF1D[numreg][numch][numeta],"Z+ Jets CR","lep")
                elif (numreg==7):
                    l.AddEntry(FF1D[numreg][numch][numeta],"TTbar CR","lep")
        l.SetTextFont(42)
        l.Draw("same")
        LTrash += [l]
        lt.DrawLatexNDC(0.12,0.8,namech+" / "+nameeta)
d.Print("SampleDependence.pdf","pdf")



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
        Files.append(ROOT.TFile.Open(HistAddress_emul + nameyear+ '_' + Samples[f]))
        print HistAddress_emul + nameyear+ '_' + Samples[f]
        for numch, namech in enumerate(channelsFF):
            l2AR1=[]
            l2AR2=[]
            for numreg, namereg in enumerate(regionsFF):
                l3AR1=[]
                l3AR2=[]
                for numvar, namevar in enumerate(variablesFF):
                    l4AR1=[]
                    l4AR2=[]
                    for numeta1, nameeta1 in enumerate(etaregs):
                        hAR1= Files[f].Get('AR1_2D_' + namech + '_' + namereg + '_' + nameeta1 + '_' + namevar)
                        l4AR1.append(hAR1)
                        l5AR2=[]
                        for numeta2, nameeta2 in enumerate(etaregs):
                            hAR2= Files[f].Get('AR2_3D_' + namech + '_' + namereg + '_' + nameeta1 + '_' + nameeta2 + '_' + namevar)
                            l5AR2.append(hAR2)
                        l4AR2.append(l5AR2)
                    l3AR1.append(l4AR1)
                    l3AR2.append(l4AR2)
                l2AR1.append(l3AR1)
                l2AR2.append(l3AR2)
            l1AR1.append(l2AR1)
            l1AR2.append(l2AR2)
        l0AR1.append(l1AR1)
        l0AR2.append(l1AR2)
    HistsAR1.append(l0AR1)
    HistsAR2.append(l0AR2)

#Prepare 2D and 3D FF
FF2D=[]
FF3D=[]
FF2Dup=[]
FF3Dup=[]
FF2Ddown=[]
FF3Ddown=[]
for numreg, namereg in enumerate(regions):
    l02D=[]
    l03D=[]
    l02Dup=[]
    l03Dup=[]
    l02Ddown=[]
    l03Ddown=[]
    for numch, namech in enumerate(channelsFF):
        l12D=[]
        l13D=[]
        l12Dup=[]
        l13Dup=[]
        l12Ddown=[]
        l13Ddown=[]
        for numvar, namevar in enumerate(variablesFF):
            l22D=[]
            l23D=[]
            l22Dup=[]
            l23Dup=[]
            l22Ddown=[]
            l23Ddown=[]
            for numeta1, nameeta1 in enumerate(etaregs):
                if numch>1:
                   tmp_FF=FF1D[numreg][0][numeta1].Clone()
                   tmp_FFup=FF1Dup[numreg][0][numeta1].Clone()
                   tmp_FFdown=FF1Ddown[numreg][0][numeta1].Clone()
                else:
                   tmp_FF=FF1D[numreg][numch][numeta1].Clone()
                   tmp_FFup=FF1Dup[numreg][numch][numeta1].Clone()
                   tmp_FFdown=FF1Ddown[numreg][numch][numeta1].Clone()
                tmp_FF2D=HistsAR1[0][0][0][0][numvar][0].Clone()
                tmp_FF2D.Reset("ICE")
                tmp_FF3D=HistsAR2[0][0][0][0][numvar][0][0].Clone()
                tmp_FF3D.Reset("ICE")
                tmp_FF2Dup=HistsAR1[0][0][0][0][numvar][0].Clone()
                tmp_FF2Dup.Reset("ICE")
                tmp_FF3Dup=HistsAR2[0][0][0][0][numvar][0][0].Clone()
                tmp_FF3Dup.Reset("ICE")
                tmp_FF2Ddown=HistsAR1[0][0][0][0][numvar][0].Clone()
                tmp_FF2Ddown.Reset("ICE")
                tmp_FF3Ddown=HistsAR2[0][0][0][0][numvar][0][0].Clone()
                tmp_FF3Ddown.Reset("ICE")
                for ithbin in range (1,tmp_FF2D.GetNbinsX()+2):
                    for jthbin in range (1,tmp_FF2D.GetNbinsY()+2):
                        tmp_FF2D.SetBinContent(ithbin,jthbin,tmp_FF.GetBinContent(jthbin))
                        tmp_FF2D.SetBinError(ithbin,jthbin,tmp_FF.GetBinError(jthbin))
                        tmp_FF2Dup.SetBinContent(ithbin,jthbin,tmp_FFup.GetBinContent(jthbin))
                        tmp_FF2Ddown.SetBinContent(ithbin,jthbin,tmp_FFdown.GetBinContent(jthbin))
                l22D.append(tmp_FF2D)
                l22Dup.append(tmp_FF2Dup)
                l22Ddown.append(tmp_FF2Ddown)
                l33D=[]
                l33Dup=[]
                l33Ddown=[]
                for numeta2, nameeta2 in enumerate(etaregs):
                    if numch>1:
                       tmp_FF_copy=FF1D[numreg][1][numeta2].Clone()
                       tmp_FF_copyup=FF1Dup[numreg][1][numeta2].Clone()
                       tmp_FF_copydown=FF1Ddown[numreg][1][numeta2].Clone()
                    else:
                       tmp_FF_copy=FF1D[numreg][numch][numeta2].Clone()
                       tmp_FF_copyup=FF1Dup[numreg][numch][numeta2].Clone()
                       tmp_FF_copydown=FF1Ddown[numreg][numch][numeta2].Clone()
                    for ithbin in range (1,tmp_FF3D.GetNbinsX()+2):
                        for jthbin in range (1,tmp_FF3D.GetNbinsY()+2):
                            for kthbin in range (1,tmp_FF3D.GetNbinsZ()+2):
                                tmp_FF3D.SetBinContent(ithbin,jthbin,kthbin,tmp_FF.GetBinContent(jthbin)*tmp_FF_copy.GetBinContent(kthbin))
                                tmp_FF3D.SetBinError(ithbin,jthbin,kthbin,math.sqrt(tmp_FF.GetBinError(jthbin)**2+tmp_FF_copy.GetBinError(kthbin)**2))
                                tmp_FF3Dup.SetBinContent(ithbin,jthbin,kthbin,tmp_FFup.GetBinContent(jthbin)*tmp_FF_copyup.GetBinContent(kthbin))
                                tmp_FF3Ddown.SetBinContent(ithbin,jthbin,kthbin,tmp_FFdown.GetBinContent(jthbin)*tmp_FF_copydown.GetBinContent(kthbin))
                    l33D.append(tmp_FF3D)
                    l33Dup.append(tmp_FF3Dup)
                    l33Ddown.append(tmp_FF3Ddown)
                l23D.append(l33D)
                l23Dup.append(l33Dup)
                l23Ddown.append(l33Ddown)
            l12D.append(l22D)
            l13D.append(l23D)
            l12Dup.append(l22Dup)
            l13Dup.append(l23Dup)
            l12Ddown.append(l22Ddown)
            l13Ddown.append(l23Ddown)
        l02D.append(l12D)
        l03D.append(l13D)
        l02Dup.append(l12Dup)
        l03Dup.append(l13Dup)
        l02Ddown.append(l12Ddown)
        l03Ddown.append(l13Ddown)
    FF2D.append(l02D)
    FF3D.append(l03D)
    FF2Dup.append(l02Dup)
    FF3Dup.append(l03Dup)
    FF2Ddown.append(l02Ddown)
    FF3Ddown.append(l03Ddown)

fout = ROOT.TFile(year[0]+"_Nonprompt.root","recreate")

# Start making stack plots
for numreg, namereg in enumerate(regionsFF):
    for numvar, namevar in enumerate(variablesFF):
        tmp_VR=HistsFF[0][0][0][0][0][numvar].Clone()
        tmp_VR.Reset("ICE")
        tmp_VRup=HistsFF[0][0][0][0][0][numvar].Clone()
        tmp_VRup.Reset("ICE")
        tmp_VRdown=HistsFF[0][0][0][0][0][numvar].Clone()
        tmp_VRdown.Reset("ICE")
        tmp_VR_e=HistsFF[0][0][0][0][0][numvar].Clone()
        tmp_VR_e.Reset("ICE")
        tmp_VR_mu=HistsFF[0][0][0][0][0][numvar].Clone()
        tmp_VR_mu.Reset("ICE")
        for numch, namech in enumerate(channelsFF):
            if (numch==0):
                MRreg = 3
            else:
                MRreg = 3
            for numeta1, nameeta1 in enumerate(etaregs):
                tmp_AR1=HistsAR1[0][0][numch][numreg][numvar][numeta1].Clone()
                for f in range(1,len(Samples)-4): # subtracted by promt contribution
                    tmp_AR1.Add(HistsAR1[0][f][numch][numreg][numvar][numeta1],-1)
                tmp_VR1=tmp_AR1.Clone()
                tmp_VR1up=tmp_AR1.Clone()
                tmp_VR1down=tmp_AR1.Clone()
                tmp_VR1.Multiply(tmp_AR1, FF2D[MRreg][numch][numvar][numeta1], 1.0, 1.0, "B")
                tmp_VR1up.Multiply(tmp_AR1, FF2Dup[MRreg][numch][numvar][numeta1], 1.0, 1.0, "B")
                tmp_VR1down.Multiply(tmp_AR1, FF2Ddown[MRreg][numch][numvar][numeta1], 1.0, 1.0, "B")
                for ithbin in range (1,tmp_VR1.GetNbinsX()+2):
                    for jthbin in range (1,tmp_VR1.GetNbinsY()+2):
                        if tmp_VR1.GetBinContent(ithbin,jthbin)>0:#Throw away negative numbers
                           tmp_VR.SetBinContent(ithbin,tmp_VR.GetBinContent(ithbin)+tmp_VR1.GetBinContent(ithbin,jthbin))
                           if numch==0:
                              tmp_VR_e.SetBinContent(ithbin,tmp_VR_e.GetBinContent(ithbin)+tmp_VR1.GetBinContent(ithbin,jthbin))
                           else:
                              tmp_VR_mu.SetBinContent(ithbin,tmp_VR_mu.GetBinContent(ithbin)+tmp_VR1.GetBinContent(ithbin,jthbin))
                           tmp_VR.SetBinError(ithbin,math.sqrt(tmp_VR.GetBinError(ithbin)**2+tmp_VR1.GetBinContent(ithbin,jthbin)**2))
                        if tmp_VR1up.GetBinContent(ithbin,jthbin)>0:#Throw away negative numbers
                           tmp_VRup.SetBinContent(ithbin,tmp_VRup.GetBinContent(ithbin)+tmp_VR1up.GetBinContent(ithbin,jthbin))
                        if tmp_VR1down.GetBinContent(ithbin,jthbin)>0:#Throw away negative numbers
                           tmp_VRdown.SetBinContent(ithbin,tmp_VRdown.GetBinContent(ithbin)+tmp_VR1down.GetBinContent(ithbin,jthbin))
                for numeta2, nameeta2 in enumerate(etaregs):
                    tmp_AR2=HistsAR2[0][0][numch][numreg][numvar][numeta1][numeta2].Clone()
                    for f in range(1,len(Samples)-4): # subtracted by promt contribution
                        tmp_AR2.Add(HistsAR2[0][f][numch][numreg][numvar][numeta1][numeta2],-1)
                    tmp_VR2=tmp_AR2.Clone()
                    tmp_VR2up=tmp_AR2.Clone()
                    tmp_VR2down=tmp_AR2.Clone()
                    tmp_VR2.Multiply(tmp_AR2, FF3D[MRreg][numch][numvar][numeta1][numeta2], 1.0, 1.0, "B")
                    tmp_VR2up.Multiply(tmp_AR2, FF3Dup[MRreg][numch][numvar][numeta1][numeta2], 1.0, 1.0, "B")
                    tmp_VR2down.Multiply(tmp_AR2, FF3Ddown[MRreg][numch][numvar][numeta1][numeta2], 1.0, 1.0, "B")
                    for ithbin in range (1,tmp_VR2.GetNbinsX()+2):
                        for jthbin in range (1,tmp_VR2.GetNbinsY()+2):
                            for kthbin in range (1,tmp_VR2.GetNbinsZ()+2):
                                if (tmp_VR2.GetBinContent(ithbin,jthbin,kthbin)>0) and (tmp_VR.GetBinContent(ithbin)-tmp_VR2.GetBinContent(ithbin,jthbin,kthbin)>0):
                                   tmp_VR.SetBinContent(ithbin,tmp_VR.GetBinContent(ithbin)-tmp_VR2.GetBinContent(ithbin,jthbin,kthbin))
                                   if numch==0:
                                      tmp_VR_e.SetBinContent(ithbin,tmp_VR_e.GetBinContent(ithbin)-tmp_VR2.GetBinContent(ithbin,jthbin,kthbin))
                                   else:
                                      tmp_VR_mu.SetBinContent(ithbin,tmp_VR_mu.GetBinContent(ithbin)-tmp_VR2.GetBinContent(ithbin,jthbin,kthbin))
                                   tmp_VR.SetBinError(ithbin,math.sqrt(tmp_VR.GetBinError(ithbin)**2+tmp_VR1.GetBinContent(ithbin,jthbin,kthbin)**2))
                                if (tmp_VR2up.GetBinContent(ithbin,jthbin,kthbin)>0) and (tmp_VRup.GetBinContent(ithbin)-tmp_VR2up.GetBinContent(ithbin,jthbin,kthbin)>0):
                                    tmp_VRup.SetBinContent(ithbin,tmp_VRup.GetBinContent(ithbin)-tmp_VR2up.GetBinContent(ithbin,jthbin,kthbin))
                                if (tmp_VR2down.GetBinContent(ithbin,jthbin,kthbin)>0) and (tmp_VRdown.GetBinContent(ithbin)-tmp_VR2down.GetBinContent(ithbin,jthbin,kthbin)>0):
                                    tmp_VRdown.SetBinContent(ithbin,tmp_VRdown.GetBinContent(ithbin)-tmp_VR2down.GetBinContent(ithbin,jthbin,kthbin))
        tmp_VR.SetFillColor(colors[8])
        tmp_VR.SetLineColor(colors[0])
        tmp_VR.SetName('emul_'+ namereg + '_' + namevar)
        tmp_VR_e.SetFillColor(colors[9])
        tmp_VR_e.SetLineColor(colors[0])
        tmp_VR_mu.SetFillColor(colors[10])
        tmp_VR_mu.SetLineColor(colors[0])
        tmp_VRup.SetName('emul_'+ namereg + '_' + namevar+ '_nonprompt_Up')
        tmp_VRdown.SetName('emul_'+ namereg + '_' + namevar+ '_nonprompt_Down')
        tmp_VR.Write()
        tmp_VRup.Write()
        tmp_VRdown.Write()
        HHCorr=[]
        HHsignalCorr=[]
        HHprit=HistsFF[0][0][0][0][numreg][numvar].Clone()
        HHprit.Reset("ICE")
        for f in range(len(Samples)):
            if ('TTbar' not in Samples[f]) and ('others' not in Samples[f]):
               if 'LFV' in Samples[f]:
                   HHsignalCorr.append(HistsFF[0][f][0][0][numreg][numvar])
               elif 'ZZ' in Samples[f]:
                   if ZZcorr:
                      HistsFF[0][f][0][0][numreg][numvar].Multiply(HistsFF[0][f][0][0][numreg][numvar], TF[0][0][numvar], 1.0, 1.0, "B")
                   HHCorr.append(HistsFF[0][f][0][0][numreg][numvar])
               else:
                   HHCorr.append(HistsFF[0][f][0][0][numreg][numvar])
            if ('ZZ' in Samples[f]) or ('TTV' in Samples[f]) or ('WZ' in Samples[f]):
                HHprit.Add(HistsFF[0][f][0][0][numreg][numvar])
        HHCorr.append(tmp_VR_e)
        HHCorr.append(tmp_VR_mu)
        HHpritup=HHprit.Clone()
        HHpritdown=HHprit.Clone()
        HHprit.Add(tmp_VR)
        HHpritup.Add(tmp_VRup)
        HHpritdown.Add(tmp_VRdown)
        HHprit.SetBinContent(HHprit.GetXaxis().GetNbins(),HHprit.GetBinContent(HHprit.GetXaxis().GetNbins()) + HHprit.GetBinContent(HHprit.GetXaxis().GetNbins()+1))
        HHpritup.SetBinContent(HHpritup.GetXaxis().GetNbins(),HHpritup.GetBinContent(HHpritup.GetXaxis().GetNbins()) + HHpritup.GetBinContent(HHpritup.GetXaxis().GetNbins()+1))
        HHpritdown.SetBinContent(HHpritdown.GetXaxis().GetNbins(),HHpritdown.GetBinContent(HHpritdown.GetXaxis().GetNbins()) + HHpritdown.GetBinContent(HHpritdown.GetXaxis().GetNbins()+1))
        binwidth= array( 'd' )
        bincenter= array( 'd' )
        yvalue= array( 'd' )
        yerrup= array( 'd' )
        yerrdown= array( 'd' )
        yvalueRatio= array( 'd' )
        yerrupRatio= array( 'd' )
        yerrdownRatio= array( 'd' )
        for b in range(HHprit.GetNbinsX()):
            if HHprit.GetBinContent(b+1)==0:
               HHprit.SetBinContent(b+1,0.0000001)
            if HHpritup.GetBinContent(b+1)==0:
               HHpritup.SetBinContent(b+1,0.0000001)
            if HHpritdown.GetBinContent(b+1)==0:
               HHpritdown.SetBinContent(b+1,0.0000001)
            content=HHprit.GetBinContent(b+1)
            contentup=HHpritup.GetBinContent(b+1)
            contentdown=HHpritdown.GetBinContent(b+1)
            binwidth.append(HHprit.GetBinWidth(b+1)/2)
            bincenter.append(HHprit.GetBinCenter(b+1))
            yvalue.append(content)
            ymax=content
            ymin=content
            if contentup>ymax:
               ymax=contentup
            if contentdown>ymax:
               ymax=contentdown
            if contentup<ymin:
               ymin=contentup
            if contentdown<ymin:
               ymin=contentdown
            yerrup.append(ymax-content)
            yerrupRatio.append((ymax-content)/content)
            yerrdown.append(content-ymin)
            yerrdownRatio.append((content-ymin)/content)
            yvalueRatio.append(1)
        if FFsys:
           Error=ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalue,binwidth,binwidth,yerrdown,yerrup)
           ErrorRatio=ROOT.TGraphAsymmErrors(len(bincenter),bincenter,yvalueRatio,binwidth,binwidth,yerrdownRatio,yerrupRatio)
        else:
           Error=ROOT.TGraphAsymmErrors()
           ErrorRatio=ROOT.TGraphAsymmErrors()
        namef='SR_DataDriven'
            #if MakePlots:
        stackPlotsFF(HHCorr, HHsignalCorr, Error, ErrorRatio, SamplesNameFF, namef, 'emul', namereg, year[0],namevar,variablesNameFF[numvar])

fout.Close()
