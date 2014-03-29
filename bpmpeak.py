#!/usr/bin/env python
import os,re,pickle,sys
from ROOT import gPad,gROOT,TSpectrum,TF1,TCanvas,TFile,TProfile,TH1F,TTree,TObject
from array import array
from bcmconst import *
from signalfilter import *

def checkrunavail(path,run,checkorder=["bpm","ring","g2p"]):
    rootfiles=filter(lambda p:"%i"%run in p,os.listdir(path))
    filelist=[]
    for order in checkorder:
	pattern="%s_?[LR]?_%i.root"%(order,run)
	t=[]
	for f in rootfiles:
	    if re.match(pattern,f):t.append(f)
	t.sort()
	filelist+=t
    if len(filelist)<1:
	print "sorry no file found for run %s in %s"%(run,path)
	return False
    return os.path.join(path,filelist[0])

def find_peak(tree,leaf,cut=""):
    peaky=0
    #fitgaus=TF1("fitgaus","gaus")
    try:
	tree.Draw(leaf+">>h%s(300)"%leaf[-1],cut)
    except:
	print "error! can't find leaf %s in tree, please check"%(leaf)
	sys.exit()
    #find peak
    h1=gPad.GetPrimitive("h%s"%leaf[-1])
    h1.SetTitle(leaf)
    hmean=h1.GetMean()
    hsigma=h1.GetRMS()
    cut2="%s>%i&&%s<%i"%(leaf,hmean-4*hsigma,leaf,hmean+4*hsigma)
    if len(cut)>1:cut2="%s&&%s"%(cut,cut2)
    tree.Draw(leaf+">>hh%s(300)"%leaf[-1],cut2)
    h1=gPad.GetPrimitive("hh%s"%leaf[-1])
    hmean=h1.GetMean()
    hsigma=h1.GetRMS()
    s=TSpectrum()
    nfound=s.Search(h1,1,"nobackground")
    for j in range(nfound):
	peakx_new=s.GetPositionX()[j]
	peaky_new=s.GetPositionY()[j]
	if peaky_new>peaky:
	    peakx=peakx_new
	    peaky=peaky_new
    #fit
    if hsigma>1:fitcut="%s>%i&&%s<%i"%(leaf,peakx-hsigma,leaf,peakx+hsigma)
    else:fitcut=""
    if len(cut)>1:fitcut="%s&&%s"%(cut,fitcut)
    fitgaus=TF1("fitgaus","gaus",peakx-hsigma,peakx+hsigma)
    tree.Fit("fitgaus",leaf,fitcut,"QRsame")
    mean=fitgaus.GetParameter("Mean")
    sigma=fitgaus.GetParameter("Sigma")
    return (mean,sigma)

#if ped is list, threshold for 8 bpm antenna,else it is used for checking if it is pedestal run(=1 or 0)
def bpm_peak(run,ped=0,rootpath=os.getenv("REPLAY_OUT_PATH"),treename="T"):
    pkldir=os.path.join(rootpath,"pkl")
    pklpath=os.path.join(pkldir,"bpmpeak_%i_%s_%i.pkl"%(ped,treename,run))
    if not os.path.exists(pkldir):os.mkdir(pkldir)
    elif os.path.exists(pklpath):
	result=pickle.load(open(pklpath,"rb"))
	return result
    runpath=checkrunavail(rootpath,run)
    if not runpath:return False
    if not decodefromrootfile(runpath,treename):sys.exit()
    hbpm=[0]*8
    #get average curr
    nocurr=0.01 #below this current will deal as no signal
    currshift=500
    rawcurr=pickle.load(open(os.path.join(pkldir,"bpmraw_curr_%i.pkl"%run),"rb"))
    curr=filter(lambda x:x>nocurr,rawcurr)
    avecurr=sum(curr)/len(curr) if len(curr)>0 else 0
    bpmmean,bpmsigma=[0]*8,[0]*8
    #ped or signal cut
    if ped>0:
	curr1=[1 if c<nocurr else 0 for c in rawcurr]
	curr2=[0]*currshift+curr1[:-currshift]
	bpmavail=map(lambda c1,c2:c1*c2,curr1,curr2)
    else:
	bpmavail=pickle.load(open(os.path.join(pkldir,"bpmraw_bpmavail_%i.pkl"%run),"rb"))
    #fill tree
    #bpmrootfile=TFile(os.path.join(rootpath,"bpmpeak_%i.root"%run),"RECREATE")
    bpmtree=TTree("bpm","bpm")
    currleaf="curr"
    bpmleaves=[currleaf]
    Vbpmleaves=currleaf
    Vbpm=array("f",[0.0]*9)
    for i in range(8):
	bpmleaves.append("bpmraw%i"%i)
	Vbpmleaves+=":%s/F"%("bpmraw%i"%i)
    bpmbranch=bpmtree.Branch("bpmraw",Vbpm,Vbpmleaves)
    print "filling bpm raw tree for getting peak for run %s"%run   
    numdata=len(rawcurr)
    rawbpm=pickle.load(open(os.path.join(pkldir,"bpmraw_sbpm_%i.pkl"%run),"rb"))
    #only use the most stable period event
    if ped<1:
	#split beam
	posrms={}
	rms_pure=array("f",[])
	aveevents,aveeventsbg=1000,300
	Vbpmave,Vbpmavebg=Cavestack(aveevents),Cavestack(aveeventsbg)
	for i in range(numdata):
	    if i%10000==0:print "generated %i rms events, %i events left"%(i,numdata-i)
	    x=(rawbpm[0][i]-rawbpm[2][i])/(rawbpm[0][i]+rawbpm[2][i])
	    if bpmavail[i]<0.5:
		rms_pure.append(-0.03)
		continue
	    Vbpmave.push(x)
	    Vbpmavebg.push(x)
	    rms1,rms2=Vbpmave.rms(0),Vbpmave.rms(1)
	    rmsbg1,rmsbg2=Vbpmavebg.rms(0),Vbpmavebg.rms(1)
	    rms_pure.append(abs(rms1-rmsbg1)+abs(rms2-rmsbg2))
	del Vbpmave,Vbpmavebg
	moveentries=getposrms(rms_pure,"pic/rms_%i.png"%(run))
	#get most stable period
	periodid,lenperiod=0,0
	for i in range(len(moveentries)):
	    if moveentries[i][1]-moveentries[i][0]>lenperiod:
		periodid=i
		lenperiod=moveentries[i][1]-moveentries[i][0]
	entryrange=moveentries[periodid]
    else:entryrange=[0,numdata]
    #fill tree
    for i in range(entryrange[0],entryrange[1]):
	if bpmavail[i]:
	    if i%10000==0:print "filling %i events, %i left"%(i,numdata-i)
	    Vbpm[0]=float(rawcurr[i])
	    for j in range(8):
		Vbpm[j+1]=float(rawbpm[j][i])
	    bpmtree.Fill()
    #bpmtree.Write("",TObject.kOverwrite)
    
    #fit
    gROOT.SetBatch(True)
    c1=TCanvas("c1","BPM peak",1920,1080)
    c1.Divide(4,2)
    for i in range(8):
	c1.cd(i+1)
	print "getting bpm peak for raw %i"%(i+1)
	try:
	    (bpmmean[i],bpmsigma[i])=find_peak(bpmtree,bpmleaves[i+1])
	except:
	    bpmrootfile.Close()
    if not os.path.exists("pic"):os.makedirs("pic")
    c1.Print("pic/peak_%i_%i.png"%(run,ped),"png")
    #c1.WaitPrimitive()
    #bpmrootfile.Close()
    result=(bpmmean,bpmsigma)
    pickle.dump(result,open(pklpath,"wb",-1),-1)
    del rawcurr,rawbpm,bpmavail
    return result
	