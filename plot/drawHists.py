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
import ctypes
import gc
import argparse
TGaxis.SetMaxDigits(2)

def cutFlowTable(hists, samples, regions, ch, year,caption='2017', nsig=6):
    mcSum = list(0 for i in xrange(0,len(regions)))
    mcErr = list(0 for i in xrange(0,len(regions)))
    if not ch==1:
        showData = True
    else:
        showData = False
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table*}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{' + caption +"}\n"
    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{|l|' + ''.join([('' if (c!='lll') and (c!='lllOffZ') and (c!='lllOffZMetg20') and (c!='lllOffZMetg20Jetgeq1Bleq1') else 'l|') for c in regions]).strip() +'}' + "\n"
    table += '\\hline' + "\n"
    table += 'Samples ' + ''.join([('' if (c!='lll') and (c!='lllOffZ') and (c!='lllOffZMetg20') and (c!='lllOffZMetg20Jetgeq1Bleq1') else ' & '+c) for c in regions]).strip() + '\\\\' + "\n"
    table += '\\hline' + "\n"

    for ids, s in enumerate(samples):
        error_M=ctypes.c_double(0.)
        if ids==0:
            continue
        for idr, r in enumerate(regions):
            if ids<nsig:
               mcSum[idr] += hists[year][ids][ch][idr][2].IntegralAndError(0, (hists[year][ids][ch][idr][2].GetNbinsX() + 1), error_M, "")
               mcErr[idr] += (error_M.value)**2

    for ids, s in enumerate(samples):
        if ids==0:
            continue
        table += s
        for idr, r in enumerate(regions):
            if (regions[idr]=='lll') or (regions[idr]=='lllOffZ') or (regions[idr]=='lllOffZMetg20') or (regions[idr]=='lllOffZMetg20Jetgeq1Bleq1'):
               error_A=ctypes.c_double(0.)
               if ids<nsig:
                  if mcSum[idr]==0:
                     mcSum[idr]=1
                  table += (' & ' + str(int(hists[year][ids][ch][idr][2].IntegralAndError(0, (hists[year][ids][ch][idr][2].GetNbinsX() + 1), error_A, "")))+'$\pm$'+str(int(error_A.value))+ '[' + str(round((100*hists[year][ids][ch][idr][2].Integral()/mcSum[idr]),2)) + '\%]')
               else:
                  table += (' & ' + str(int(hists[year][ids][ch][idr][2].IntegralAndError(0, (hists[year][ids][ch][idr][2].GetNbinsX() + 1), error_A, "")))+'$\pm$'+str(int(error_A.value)) )
        table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += 'Prediction '
    for idr, r in enumerate(mcSum):
        if (regions[idr]=='lll') or (regions[idr]=='lllOffZ') or (regions[idr]=='lllOffZMetg20') or (regions[idr]=='lllOffZMetg20Jetgeq1Bleq1'):
           table += (' & ' + str(int(r)) + '$\pm$' + str(int(math.sqrt(mcErr[idr]))) )
    table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    if showData:
      table += 'Data '
      for idr, r in enumerate(regions):
          error_D=ctypes.c_double(0.)
          if (regions[idr]=='lll') or (regions[idr]=='lllOffZ') or (regions[idr]=='lllOffZMetg20') or (regions[idr]=='lllOffZMetg20Jetgeq1Bleq1'):
             table += (' & ' + str(int(hists[year][0][ch][idr][2].IntegralAndError(0, (hists[year][0][ch][idr][2].GetNbinsX() + 1), error_D, "")))+'$\pm$'+str(int(error_D.value)) )
      table += '\\\\' + "\n"
      table += '\\hline' + "\n"
      table += 'Data$/$Pred. '
      for idr, r in enumerate(mcSum):
          if r==0:
             r=0.1
          if (regions[idr]=='lll') or (regions[idr]=='lllOffZ') or (regions[idr]=='lllOffZMetg20') or (regions[idr]=='lllOffZMetg20Jetgeq1Bleq1'):
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
        if num>0:
            SignalHists[num].Scale(10)
    for num in range(1,len(hists)):
        hs.Add(hists[num])

    dummy = hists[0].Clone()
    if ('emul' not in ch) or ('OnZ' in reg):
        showData = True
    else:
        showData = False
        dummy.Reset("ICE")

    
    canvas = ROOT.TCanvas(year+ch+reg+var,year+ch+reg+var,50,50,865,780)
    canvas.SetGrid();
    canvas.SetBottomMargin(0.17)
    canvas.cd()

    if showData:
       legend = ROOT.TLegend(0.68,0.74,0.8,0.88)
    else:
       legend = ROOT.TLegend(0.68,0.74,0.8,0.833)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0);
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend2 = ROOT.TLegend(0.8,0.74,0.92,0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0);
    legend2.SetTextFont(42)
    legend2.SetTextSize(0.04)
    legend3 = ROOT.TLegend(0.68,0.64,0.83,0.74)
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
    for H in SignalHists:
        H.SetLineWidth(2)
        H.SetFillColor(0)
        H.SetLineStyle(9)
        H.Draw("histSAME")
    dummy.Draw("ex0SAME")
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
    Label_lumi = ROOT.TLatex(0.71,0.92,Lumi+" fb^{-1} (13 TeV)")
    Label_lumi.SetNDC()
    Label_lumi.SetTextFont(42)
    Label_lumi.Draw("same")
    Label_channel = ROOT.TLatex(0.2,0.8,year +" / "+ch)
    Label_channel.SetNDC()
    Label_channel.SetTextFont(42)
    Label_channel.Draw("same")
    Label_region = ROOT.TLatex(0.2,0.75,reg)
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
    dummy_ratio.Draw("ex0")
    dummy_ratio.Draw("AXISSAMEY+")
    dummy_ratio.Draw("AXISSAMEX+")
    canvas.Print(year + '/' + ch +'/'+reg+'/'+var + ".png")
    for num in range(len(SignalHists)):
        if num>0:
           SignalHists[num].Scale(0.1)
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
regions=["lll","lllOnZ","lllOffZ","lllOffZB0","lllOffZB1", "lllOffZBgeq2", "lllOffZMetl20", "lllOffZMetg20", "lllOffZMetg20Jetgeq1", "lllOffZMetg20Jet1B1", "lllOffZMetg20Jet2B1","lllOffZMetg20Jetgeq1Bleq1","lllOnZMetg20Jetgeq1Bleq1"]
channels=["eee", "emul", "mumumu"];
variables=["lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","lep3Pt","lep3Eta","lep3Phi",
           "LFVePt","LFVeEta","LFVePhi","LFVmuPt","LFVmuEta","LFVmuPhi","balPt","balEta","balPhi","Topmass",
           "lllM","lllPt","lllHt","lllMt","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx",
           "llZM","llZPt","llZDr","llZDphi","LFVTopmass","llM","llPt","llDr","llDphi",
           "Ht","Ms","ZlDr","ZlDphi","JeDr","JmuDr","tM"]
#variables=["lep1Pt"]
variablesName=["Leading lepton p_{T} [GeV]","Leading lepton #eta","Leading lepton #varphi","2nd-Leading lepton p_{T} [GeV]","2nd-Leading lepton #eta","2nd-Leading lepton #varphi","3rd-Leading lepton p_{T} [GeV]","3rd-Leading lepton #eta","3rd-Leading lepton #varphi","cLFV electron p_{T} [GeV]","cLFV electron #eta","cLFV electron #varphi","cLFV muon p_{T} [GeV]","cLFV muon #eta","cLFV muon #varphi","Bachelor lepton p_{T} [GeV]","Bachelor lepton #eta","Bachelor lepton #varphi","Standard top mass [GeV]","M(lll) [GeV]","p_{T}(lll) [GeV]","H_{t}(lll) [GeV]","M_{t}(lll) [GeV]","Leading jet p_{T} [GeV]","Leading jet #eta","Leading jet #varphi","Number of jets","Number of b-tagged jets","MET [GeV]","#varphi(MET)","Number of vertices", "M(ll) (OSSF) [GeV]", "p_{T}(ll) (OSSF) [GeV]", "Dr(ll) (OSSF)", "Dphi(ll) (OSSF)", "cLFV top mass [GeV]", "M(ll) (cLFV) [GeV]", "p_{T}(ll) (cLFV) [GeV]", "Dr(ll) (cLFV)", "Dphi(ll) (cLFV)", "H_{t} [GeV]", "M_{s} [GeV]", "Dr(Zl)", "Dphi(Zl)", "Dr(Je)", "Dr(Jmu)", "M_{T} (W candidate)"]


# set up an argument parser
parser = argparse.ArgumentParser()

parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--l', dest = 'LOCATION', default= '/eos/cms/store/user/jingyan/')

ARGS = parser.parse_args()

verbose = ARGS.VERBOSE
HistAddress = ARGS.LOCATION

Samples = ['data.root','TTV.root','WZ.root', 'ZZ.root', 'TTbar.root', 'others.root', 'LFVStVecU.root', 'LFVTtVecU.root']
SamplesName = ['data','TTV','WZ', 'ZZ', 'TTbar', 'others', 'ST_vector_emutu', 'TT_vector_emutu']
SamplesNameLatex = ['data', 'TTV', 'WZ', 'ZZ', 'TTbar', 'others', 'ST\_vector\_emutu', 'TT\_vector\_emutu']
#Samples = ['data.root','DYg50.root','DYl50.root', 'Fake.root', 'LFVStVecU.root', 'LFVTtVecU.root']
#SamplesName = ['data','DY50','DY10to50', 'Others', 'ST_vector_emutu', 'TT_vector_emutu']
#SamplesNameLatex = ['data', 'DY50', 'DY10to50', 'Others', 'ST\_vector\_emutu', 'TT\_vector\_emutu']

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
                    if 'LFV' not in Samples[f]:
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
                    if 'LFV' in Samples[f]:
                        HHsignal.append(Hists[numyear][f][numch][numreg][numvar])
                    else:
                        HH.append(Hists[numyear][f][numch][numreg][numvar])

                stackPlots(HH, HHsignal, SamplesName, namech, namereg, nameyear,namevar,variablesName[numvar])

le = '\\documentclass{article}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\begin{document}' + "\n"

print le
for numyear, nameyear in enumerate(year):
    for numch, namech in enumerate(channels):
        cutFlowTable(Hists, SamplesNameLatex, regions, numch, numyear, nameyear + ' ' + namech, len(Samples)-2 )
print '\\end{document}' + "\n"

