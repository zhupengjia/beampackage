#!/usr/bin/env python
import glob,os,re,urllib
from numpy import *
from array import *
from odsread import odsread
try:import ROOT
except:print "Error!! pyroot didn't compiled! please recompile your root!"

#find harp peak from raw harp scan data file and generate rootfile for all raw files
def genharppeak(datapath,rootfilepath):
    filename=['IHA1H04.03052012_20:31:53', 'IHA1H04.03052012_20:44:24', 'IHA1H04.03052012_21:11:54', 'IHA1H04.03052012_21:50:51',\
      'IHA1H04.03052012_22:26:25', 'IHA1H04.03052012_22:48:37', 'IHA1H04.03052012_23:32:13', 'IHA1H04.03052012_23:59:13',\
      'IHA1H04.03062012_00:07:05', 'IHA1H04.03062012_00:24:34', 'IHA1H04.03062012_00:29:08', 'IHA1H04.03062012_00:49:45',\
      'IHA1H04.03062012_00:29:08', 'IHA1H04.03062012_00:49:45', 'IHA1H04.03062012_02:38:29', 'IHA1H04.03062012_03:15:49',\
      'IHA1H04.03062012_03:36:30', 'IHA1H04.03062012_03:50:50', 'IHA1H04.03062012_04:41:34', 'IHA1H04.03062012_04:42:31',\
      'IHA1H04.03062012_04:44:19', 'IHA1H04.03062012_04:48:46', 'IHA1H04.03062012_04:51:22', 'IHA1H04.03132012_00:44:04',\
      'IHA1H04.03132012_01:17:54', 'IHA1H04.03132012_01:46:04', 'IHA1H04.03132012_02:09:43', 'IHA1H04.03132012_02:27:57',\
      'IHA1H04.03142012_03:10:48', 'IHA1H04.03142012_03:55:28', 'IHA1H04.03142012_04:28:20', 'IHA1H04.03142012_04:42:02',\
      'IHA1H04.03282012_07:35:23', 'IHA1H04.03282012_10:17:06', 'IHA1H04.03282012_10:41:04', 'IHA1H04.03282012_11:18:24',\
      'IHA1H04.04112012_03:19:33', 'IHA1H04.04112012_04:47:52', 'IHA1H04.04112012_05:00:09', 'IHA1H04.04112012_05:41:31',\
      'IHA1H04.04112012_06:06:24', 'IHA1H04.04262012_16:26:15', 'IHA1H04.05032012_01:09:53', 'IHA1H04.05032012_02:11:06',\
      'IHA1H04.05032012_02:21:20', 'IHA1H04.05032012_03:21:32', 'IHA1H04.05032012_04:11:08', 'IHA1H04.05032012_04:17:41',\
      'IHA1H04.05032012_05:10:32', 'IHA1H04.05032012_06:32:55', 'IHA1H04.05032012_07:12:46', 'IHA1H04.05092012_23:57:52',\
      'IHA1H04.05102012_00:14:19', 'IHA1H04.05102012_10:30:13', 'IHA1H04.05102012_11:27:33', 'IHA1H04.05102012_11:48:42',\
      'IHA1H05A.03052012_20:36:58', 'IHA1H05A.03052012_20:40:02', 'IHA1H05A.03052012_21:13:27', 'IHA1H05A.03052012_21:49:38',\
      'IHA1H05A.03052012_22:37:59', 'IHA1H05A.03052012_23:17:06', 'IHA1H05A.03052012_23:58:06', 'IHA1H05A.03062012_00:10:00',\
      'IHA1H05A.03062012_00:31:07', 'IHA1H05A.03062012_00:50:47', 'IHA1H05A.03062012_02:39:58', 'IHA1H05A.03062012_03:17:57',\
      'IHA1H05A.03062012_03:37:23', 'IHA1H05A.03062012_03:52:43', 'IHA1H05A.03062012_04:46:48', 'IHA1H05A.03132012_00:43:13',\
      'IHA1H05A.03082012_01:20:01', 'IHA1H05A.03082012_02:04:54', 'IHA1H05A.03082012_02:14:09', 'IHA1H05A.03082012_03:01:07',\
      'IHA1H05A.03082012_03:18:50', 'IHA1H05A.03082012_01:29:30', 'IHA1H05A.03082012_01:51:34', 'IHA1H05A.03082012_02:19:25',\
      'IHA1H05A.03082012_02:34:07', 'IHA1H05A.03082012_02:35:37', 'IHA1H05A.03082012_03:23:00',\
      'IHA1H05A.03132012_01:14:22', 'IHA1H05A.03132012_01:44:34', 'IHA1H05A.03132012_02:11:55', 'IHA1H05A.03132012_02:28:39',\
      'IHA1H05A.03142012_03:12:09', 'IHA1H05A.03142012_03:56:07', 'IHA1H05A.03142012_04:29:13', 'IHA1H05A.03142012_04:42:46',\
      'IHA1H05A.03162012_17:50:55', 'IHA1H05A.03162012_17:51:29', 'IHA1H05A.03162012_20:08:50', 'IHA1H05A.03162012_21:07:49',\
      'IHA1H05A.03162012_21:09:38', 'IHA1H05A.03162012_21:37:59', 'IHA1H05A.03282012_07:36:26', 'IHA1H05A.03282012_10:17:42',\
      'IHA1H05A.03282012_10:35:45', 'IHA1H05A.03282012_10:41:41', 'IHA1H05A.03282012_10:55:02', 'IHA1H05A.03282012_11:19:00',\
      'IHA1H05A.04112012_03:17:56', 'IHA1H05A.04112012_04:20:39', 'IHA1H05A.04112012_04:46:13', 'IHA1H05A.04112012_04:58:34',\
      'IHA1H05A.04112012_05:31:07', 'IHA1H05A.04112012_05:39:56', 'IHA1H05A.04112012_06:00:25', 'IHA1H05A.04112012_06:04:52',\
      'IHA1H05A.04262012_16:24:16', 'IHA1H05A.04262012_16:29:15', 'IHA1H05A.04262012_18:53:58', 'IHA1H05A.04262012_20:07:42',\
      'IHA1H05A.04262012_21:02:54', 'IHA1H05A.04262012_21:12:26', 'IHA1H05A.04262012_21:18:01', 'IHA1H05A.05032012_01:08:51',\
      'IHA1H05A.05032012_02:10:24', 'IHA1H05A.05032012_02:25:54', 'IHA1H05A.05032012_03:14:07', 'IHA1H05A.05032012_03:20:53',\
      'IHA1H05A.05032012_04:12:56', 'IHA1H05A.05032012_04:17:01', 'IHA1H05A.05032012_05:12:10', 'IHA1H05A.05032012_06:34:29',\
      'IHA1H05A.05032012_07:12:06', 'IHA1H05A.05092012_23:53:59', 'IHA1H05A.05102012_00:12:22', 'IHA1H05A.05102012_10:28:33',\
      'IHA1H05A.05102012_10:51:17', 'IHA1H05A.05102012_11:09:53', 'IHA1H05A.05102012_11:26:01', 'IHA1H05A.05102012_11:46:32']
    filename.sort()
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch(True)
    rootfile = ROOT.TFile(os.path.join(rootfilepath,"harpdata.root"),"RECREATE")
    s=ROOT.TSpectrum()
    leaves,tree,branches,Vleaves,Vharp=[],[],[],[],[]

    for i in range(len(filename)):
      fakename=filename[i].replace(":","").replace(".","_")
          #make tree
      tree.append(ROOT.TTree("harp%i"%i,filename[i]))
      leaves.append(["index","pos","sig"])
      Vleaves.append("index/F:pos/F:sig/F")
      Vharp.append(array("f",[0.0,0.0,0.0]))
      branches.append(tree[i].Branch(fakename,Vharp[i],Vleaves[i]))
      #fill tree
      rawname=filename[i].replace(":","_")
      hfile=os.path.join(datapath,rawname)
      print "filling %s"%hfile
      for line in open(hfile,'r'):
            if line.startswith("#"): continue
            field=line.split()
            Vharp[i][0]=float(field[0])
            Vharp[i][1]=float(field[1])
            Vharp[i][2]=float(field[2])
            #print i,Vharp[i]
            tree[i].Fill()
      tree[i].Write("",ROOT.TObject.kOverwrite)
    #draw pic
    c,h=[],[]
    t=ROOT.TLatex()
    t.SetTextAlign(10)
    t.SetTextSize(0.03)
    for j in range(len(filename)):
      fakename="G"+filename[j].replace(":","").replace(".","_")
      c.append(ROOT.TCanvas(fakename,"harp signal %s"%filename[j],1280,800))
      c[j].Divide(2,1)
      c[j].cd(1)
      fitpol=ROOT.TF1("fitpol","pol1",0,40000)
      tree[j].Fit("fitpol","pos:index","pos>0","QR")
      p0=fitpol.GetParameter(0)
      p1=fitpol.GetParameter(1)
      tree[j].Draw("pos:index","pos>0")
      c[j].cd(2)
      tree[j].Draw("sig:%f+index*%f"%(p0,p1))
      graph=ROOT.gPad.GetPrimitive("Graph")
      mean=graph.GetMean(2)
      h.append(ROOT.TProfile("h%i"%j,filename[j],300,5,90))
      tree[j].Draw("sig-%f:%f+index*%f>>h%i"%(mean,p0,p1,j))
      nfound,threshold=0,0.5
      while nfound<3:
          nfound=s.Search(h[j],2,"same",threshold)
          threshold-=0.05
      peak=[]
      for l in range(nfound):peak.append(s.GetPositionX()[l])
      #h[j].Write()
      peak.sort()
      yaxismax=h[j].GetMaximum()
      yaxismin=h[j].GetMinimum()
      tree[j].Draw("sig-%f:pos"%mean,"","sameL")
      for k in range(nfound):
          print peak[k],
          t.DrawLatex(7,yaxismax-0.035*k*(yaxismax-yaxismin),"%i, %f"%(k+1,peak[k]))
      print "\n___________%i/%i %s___________________\n"%((j+1),len(filename),filename[j])
      c[j].Write()
    rootfile.Close()






