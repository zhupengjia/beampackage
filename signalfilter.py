#!/usr/bin/env python
import re,os,glob,sys,gc,ctypes
import numpy as np
try:import cPickle as pickle
except:import pickle
from ROOT import TFile,TTree,TObject,gROOT,TSpectrum,TCanvas,gPad,TProfile
from array import array
#from pylab import plot,show,subplot
from bcmconst import *
from runinfo import pkgsetting
try:from scipy.signal import butter,freqz,lfilter
except Exception as err:
    print err
    print "sorry no scipy module found from your computer,it is needed to filter the bpm raw data infomation, please install it first"

def lowfilter(raw,cutfreq):
    n=4 # filter order
    fs=960.015 #sample rate,here is helicity rate
    fc=2*cutfreq/fs #Normalize LPF cutoff frequency to Nyquist frequency
    b,a=butter(n,fc)
    #w,h=freqz(b,a,n)
    sf=lfilter(b,a,raw)
    return np.float32(sf)

def signalave(raw,avefreq):
    fs=960.015 #trigger rate,here is helicity rate
    aveevents=int(fs/avefreq)
    Vave=Cavestack(aveevents)
    rawlen=len(raw)
    averaw=np.zeros(rawlen,dtype=np.float32)
    for i in range(rawlen):
	Vave.push(raw[i])
	averaw[i]=Vave.ave()
    del Vave
    return averaw

def getmemory():
    try:
	for line in open("/proc/meminfo","r"):
	    if "MemTotal" in line:
		return int(re.split("[:kB]","".join(re.split("\s",line)))[1])*1000
    except:return 2054132000

class Cavestack:
    def __init__(self,size=7500,ransize=0):
	try:
	    self._a=ctypes.CDLL(os.path.join(os.path.split(__file__)[0],"Cavestack.so"))
	except Exception as e:
	    print e
	    print "attention!!!!!!!!!!please run \"make\" in beampackage folder!!!!!!!!!!"
	    sys.exit()
	if ransize<=0:
	    self._a_new=self._a.avestack_new(ctypes.c_int(size))
	    self.avestack_del=self._a.avestack_del
	    self.avestack_isfull=self._a.avestack_isfull
	    self.avestack_isempty=self._a.avestack_isempty
	    self.avestack_getsize=self._a.avestack_getsize
	    self.avestack_empty=self._a.avestack_empty
	    self.avestack_push=self._a.avestack_push
	    self.avestack_ave=self._a.avestack_ave
	    self.avestack_rms=self._a.avestack_rms
	else:
	    self._a_new=self._a.ravestack_new(ctypes.c_int(size),ctypes.c_int(ransize))
	    self.avestack_del=self._a.ravestack_del
	    self.avestack_isfull=self._a.ravestack_isfull
	    self.avestack_isempty=self._a.ravestack_isempty
	    self.avestack_getsize=self._a.ravestack_getsize
	    self.avestack_empty=self._a.ravestack_empty
	    self.avestack_push=self._a.ravestack_push
	    self.avestack_ave=self._a.ravestack_ave
	    self.avestack_rms=self._a.ravestack_rms
	self.avestack_isfull.restype=ctypes.c_bool
	self.avestack_isempty.restype=ctypes.c_bool
	self.avestack_getsize.restype=ctypes.c_int
	self.avestack_ave.restype=ctypes.c_float
	self.avestack_rms.restype=ctypes.c_float
    def __del__(self):
        self.avestack_del(self._a_new)
    def isfull(self):
        return self.avestack_isfull(self._a_new)
    def isempty(self):
        return self.avestack_isempty(self._a_new)
    def getsize(self):
        return self.avestack_getsize(self._a_new)
    def empty(self):
        self.avestack_empty(self._a_new)
    def push(self,x,y=0):
        self.avestack_push(self._a_new,ctypes.c_float(x),ctypes.c_float(y))
    def ave(self,xy=0):
        return self.avestack_ave(self._a_new,ctypes.c_int(xy))
    def rms(self,xy=0):
        return self.avestack_rms(self._a_new,ctypes.c_int(xy))

#get raw data from rootfile,save it to pkl file and return as dict type,bpm signal dealt with filter
#filter1 used to get average pos,filter2 used to get raw pos that can see slow raster signal
def decodefromrootfile(runpath,treename="T",firstevent=-1,lastevent=-1,forceredecode=0,buildtree=0):
    rootfilepath,runfilename=os.path.split(runpath)
    run=int(re.split("[_]",re.split("[.]",runfilename)[-2])[-1])
    if run<100:run=int(re.split("[_]",re.split("[.]",os.path.split(runpath)[1])[-2])[-2])
    arm="L" if run<20000 else "R"
    #get pkl file
    pkldir=os.path.join(rootfilepath,"pkl")
    if not os.path.exists(pkldir):os.makedirs(pkldir)
    pklpathn=[["rbpm","sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail"],["raster","clock","event"]]
    pklpath=[{},{}]
    pklexists=[{},{}] #file exists?
    pklredecode=[False,False]
    for i in range(2):
	for j in range(len(pklpathn[i])):
	    n=pklpathn[i][j]
	    pklpath[i][n]=os.path.join(pkldir,"bpmraw_%s_%i.pkl"%(n,run))
	    if os.path.exists(pklpath[i][n]):pklexists[i][n]=True
	    else:pklexists[i][n]=False
	pklredecode[i]=forceredecode or not all(pklexists[i].values())
    #check event
    eventtolerate=1000
    if not pklredecode[0]:
	hapevent=pickle.load(open(pklpath[0]["hapevent"],"rb"))
	if not((firstevent<0 or hapevent[0]-firstevent<eventtolerate) and (lastevent<0 or lastevent-hapevent[-1]<eventtolerate)):
	    pklredecode[0]=True
	del hapevent
    if not pklredecode[1]:
	event=pickle.load(open(pklpath[1]["event"],"rb"))
	if not((firstevent<0 or event[0]-firstevent<eventtolerate) and (lastevent<0 or lastevent-event[-1]<eventtolerate)):
	    pklredecode[1]=True
	del event
    redecode=any(pklredecode)
    pkgpara=pkgsetting()
    if run in pkgpara.hapunavail:fastbus=True
    else:fastbus=False
    if redecode:
	#get redundant rootfiles
	rootfiles=glob.glob(os.path.join(rootfilepath,runfilename.replace(".root","_*.root")))
	rootfiles.append(runpath)
	rootfiles.sort()
	#raw leaves
	eventleaf="fEvtHdr.fEvtNum" #event number
	if pklredecode[0]:
	    thappexprefix="happex.%s."%arm
	    numringleaf="%snumring"%thappexprefix
	    bpmrawleaf,hapbcmrawleaf,scalerbcmrawleaf=[],[0,0],[0,0]
	    for i in range(8): #for bpm
		whichbpm="A" if i<4 else "B"
		if fastbus:bpmrawleaf.append("%surb.BPM%s.rawcur.%i"%(arm,whichbpm,i%4+1))
		else:bpmrawleaf.append("%sBPM%s.rawcur.%i"%(thappexprefix,whichbpm,i%4+1))
	    hapbcmrawleaf[0]="%sbcm_up"%thappexprefix #for bcm
	    hapbcmrawleaf[1]="%sbcm_down"%thappexprefix
	    scalerbcmrawleaf[0]="evleft_bcm_upr" if run<20000 else "evright_bcm_upr"
	    scalerbcmrawleaf[1]="evleft_bcm_downr" if run<20000 else "evright_bcm_downr"
	    #bcm const
	    const=bcmconst()
	    hapbcmconst=[const.getconst(run,"happex",a) for a in ["up","down"]]
	    hapbcmavail=[True if a else False for a in hapbcmconst]
	    hapbcmavailall=hapbcmavail[0] or hapbcmavail[1]
	    scalerbcmconst=[const.getconst(run,"sis3800",a,"slow") for a in ["up","down"]]
	    scalerbcmavail=[True if a else False for a in scalerbcmconst]
	if pklredecode[1]:
	    rasterrawleaf=["%srb.Raster.rawcur.x"%arm,"%srb.Raster.rawcur.y"%arm] #for raster
	    rasterrawleaf+=["%srb.Raster.rawcurSL.x"%arm,"%srb.Raster.rawcurSL.y"%arm]
	    clkrawleaf=arm+"clk.fastclk" #for clock
	#get raw data
	try:
	    if pklredecode[0]:
		bpmraw=[array("i",[]),array("i",[]),array("i",[]),array("i",[]),\
		    array("i",[]),array("i",[]),array("i",[]),array("i",[])]
		hapevent=array("I",[])
		curr=array("f",[])
	    if pklredecode[1]:
		rasterraw=[array("i",[]),array("i",[]),array("i",[]),array("i",[])]
		event=array("I",[])
		clkraw=array("I",[])
		rasteravail=True
		clkavail=True
	    for runpath in rootfiles:
		print "decoding raw data from rootfile %s..."%runpath
		rootfile=TFile(runpath,"READ")
		if rootfile.IsZombie():
		    print "Error! no file %s exists!"%runpath
		    if runpath==rootfiles[0]:sys.exit()
		    else:continue
		tree=rootfile.Get(treename)
		events=tree.GetEntries()
		#get leaf and branch
		try:leventleaf=tree.GetLeaf(eventleaf)
		except:
		    print "Error!! no leaf %s in your rootfile!!"%eventleaf
		    return False
		if pklredecode[0]:
		    try:lnumringleaf=tree.GetLeaf(numringleaf)
		    except:
			print "Error!! no leaf %s in your rootfile!!"%numringleaf
			return False
		    try:
			if fastbus:
			    lbpmrawleaf=[tree.GetLeaf(l) for l in bpmrawleaf]
			else:    
			    bpmrawbranch=[tree.GetBranch(l) for l in bpmrawleaf]
			    lbpmrawleaf=[b.GetLeaf("data") for b in bpmrawbranch]
		    except:
			print "Error!! no leaf %s in your rootfile!!"%bpmrawleaf[0]
			return False
		    try:
			hapbcmrawbranch=[tree.GetBranch(l) for l in hapbcmrawleaf]
			lhapbcmrawleaf=[b.GetLeaf("data") for b in hapbcmrawbranch]
		    except:
			print "Error!! no leaf %s in your rootfile!!"%hapbcmrawleaf[0]
			print "will try to use scaler bcm info instead since you didn't replay happex bcm info"
			hapbcmavailall=False
		    if not hapbcmavailall:
			try:lscalerbcmrawleaf=[tree.GetLeaf(scalerbcmrawleaf[i]) for i in range(2)]
			except:
			    print "Error!! no leaf %s in your rootfile!!"%scalerbcmrawleaf[0]
			    return False
		if pklredecode[1]:	
		    try:lclkrawleaf=tree.GetLeaf(clkrawleaf)
		    except:
			print "Error!! no leaf %s in your rootfile!will leave it as empty!"%clkrawleaf
			clkavail=False
		    try:lrasterrawleaf=[tree.GetLeaf(rasterrawleaf[i]) for i in range(4)]
		    except:
			print "Error!! no leaf %s in your rootfile!will leave it as empty!"%rasterrawleaf[0]
			rasteravail=False
		    
		#decode from rootfile
		for e in range(events):
		    if e%1000==0:
			print "decoding %i events, %i left"%(e,events-e)
		    tree.GetEntry(e)
		    try:ee=int(leventleaf.GetValue())
		    except Exception as err:
			print err
			print "Error!! no leaf %s in your rootfile!!"%eventleaf
			return False
		    if pklredecode[1]:
			event.append(ee)
			if clkavail:
			    try:clkraw.append(int(lclkrawleaf.GetValue()))
			    except:
				clkavail=False
				print "Error!! no leaf %s in your rootfile!will leave it as empty!"%clkrawleaf
				if not pklredecode[0] and pklexists[1]["raster"] and pklexists[1]["event"]:
				    pklredecode[1]=False
				    break
			if rasteravail:
			    try:
				for i in range(4):
				    rasterraw[i].append(int(lrasterrawleaf[i].GetValue()))
			    except:
				rasteravail=False
				print "Error!! no leaf %s in your rootfile!will leave it as empty!"%rasterrawleaf[0]
				if not pklredecode[0] and pklexists[1]["clock"] and pklexists[1]["event"]:
				    pklredecode[1]=False
				    break
		    if pklredecode[0]:
			bcmraw=[False,False]
			if not hapbcmavailall:
			    for i in range(2):
				if scalerbcmavail[i]:
				    try:
					bcmraw[i]=lscalerbcmrawleaf[i].GetValue()
				    except Exception as err:
					print err
					print "Error!! no leaf %s in your rootfile!!"%scalerbcmrawleaf[i]
					return False
				    bcmraw[i]=getcurr(bcmraw[i],scalerbcmconst[i],"sis3800","slow")
			    bcmraw=[a for a in bcmraw if a]
			if fastbus:
			    bpmrawfastbus=[int(lbpmrawleaf[i].GetValue()) for i in range(8)]
			#for happex data
			numring=int(lnumringleaf.GetValue())
			if numring<1:continue
			for i in range(numring):
			    if fastbus:
				for j in range(8):
				    bpmraw[j].append(bpmrawfastbus[j])
			    else:
				for j in range(8):
				    bpmraw[j].append(int(lbpmrawleaf[j].GetValue(i)))
			    if hapbcmavailall:
				for j in range(2):
				    if hapbcmavail[j]:
					bcmraw[j]=lhapbcmrawleaf[j].GetValue(i)
					bcmraw[j]=getcurr(bcmraw[j],hapbcmconst[j],"happex")
				bcmraw=[a for a in bcmraw if a]
			    curr.append(sum(bcmraw)/len(bcmraw))
			    hapevent.append(ee)
		try:
		    del tree,rootfile,leventleaf,
		    if pklredecode[0]:
			del lbpmrawleaf,hapbcmrawbranch,lhapbcmrawleaf,lnumringleaf
			if not fastbus:bpmrawbranch
		    if pklredecode[1]:
			del lclkrawleaf,lrasterrawleaf
		except:continue
	except Exception as err:
	    print "Error!! something wrong with your rootfiles!! Please check the error message below and check your rootfile!!"
	    print err
	    return False
	if pklredecode[0]:
	    for i in range(8):
		bpmraw[i]=np.asarray(bpmraw[i],dtype=np.int32)
	    try:pickle.dump(bpmraw,open(pklpath[0]["rbpm"],"wb",-1),-1)
	    except:sys.exit("\n\n\n\nError!failed to dump for raw bpm data,do you have write permission in dir %s?"%pkldir)
	    #deal signal
	    fbpmraw_sslow,fbpmraw_slow,fbpmraw_fast=[0]*8,[0]*8,[0]*8
	    for i in range(7,-1,-1):
		print "filtering for bpm channel %i"%(i+1)
		if pkgpara.filter1<500:
		    if "ave" in pkgpara.filter1type or fastbus:
			fbpmraw_slow[i]=signalave(bpmraw[-1],pkgpara.filter1)
		    else:
			fbpmraw_slow[i]=lowfilter(bpmraw[-1],pkgpara.filter1)
		else:
		    fbpmraw_slow[i]=bpmraw[-1]
		if pkgpara.filter2<500:
		    if "ave" in pkgpara.filter2type or fastbus:
			fbpmraw_fast[i]=signalave(bpmraw[-1],pkgpara.filter2)
		    else:
			fbpmraw_fast[i]=lowfilter(bpmraw[-1],pkgpara.filter2)
		else:
		    fbpmraw_fast[i]=bpmraw[-1]
		if pkgpara.filter3<500:
		    if "ave" in pkgpara.filter3type or fastbus:
			fbpmraw_sslow[i]=signalave(bpmraw[-1],pkgpara.filter3)
		    else:
			fbpmraw_sslow[i]=lowfilter(bpmraw[-1],pkgpara.filter3)
		else:
		    fbpmraw_sslow[i]=bpmraw[-1]
		del bpmraw[-1]
		#gc.collect() #recycle memory
	    #build bpm avail variable
	    curravail=0.02
	    scurrshift=int(1000/pkgpara.filter1+0.9)
	    fcurrshift=int(1000/pkgpara.filter2+0.9)
	    sscurrshift=int(1000/pkgpara.filter3+0.9)
	    #at least cut 2000 events(2s)
	    minshift=2000
	    scurrshift=minshift if scurrshift<minshift else scurrshift
	    fcurrshift=minshift if fcurrshift<minshift else fcurrshift
	    sscurrshift=minshift if sscurrshift<minshift else sscurrshift
	    curr1=[1 if c>curravail else 0 for c in curr]
	    scurr2=[0]*scurrshift+curr1[:-scurrshift]
	    fcurr2=[0]*fcurrshift+curr1[:-fcurrshift]
	    sscurr2=[0]*sscurrshift+curr1[:-sscurrshift]
	    sbpmavail,fbpmavail,ssbpmavail=array("i",[]),array("i",[]),array("i",[])
	    for i in range(len(curr1)):
		sbpmavail.append(curr1[i]*scurr2[i])
		fbpmavail.append(curr1[i]*fcurr2[i])
		ssbpmavail.append(curr1[i]*sscurr2[i])
	    #dump
	    try:
		pickle.dump(fbpmraw_slow,open(pklpath[0]["sbpm"],"wb",-1),-1)
		pickle.dump(fbpmraw_sslow,open(pklpath[0]["ssbpm"],"wb",-1),-1)
		pickle.dump(fbpmraw_fast,open(pklpath[0]["fbpm"],"wb",-1),-1)
		pickle.dump(np.asarray(curr,dtype=np.float32),open(pklpath[0]["curr"],"wb",-1),-1)
		pickle.dump(np.asarray(sbpmavail,dtype=np.bool),open(pklpath[0]["bpmavail"],"wb",-1),-1)
		pickle.dump(np.asarray(fbpmavail,dtype=np.bool),open(pklpath[0]["fbpmavail"],"wb",-1),-1)
		pickle.dump(np.asarray(ssbpmavail,dtype=np.bool),open(pklpath[0]["sbpmavail"],"wb",-1),-1)
		pickle.dump(np.asarray(hapevent,dtype=np.uint32),open(pklpath[0]["hapevent"],"wb",-1),-1)
	    except:sys.exit("\n\n\n\nError!failed to dump for filtered bpm data,do you have write permission in dir %s?"%pkldir)
	    del fbpmraw_slow,fbpmraw_sslow,fbpmraw_fast,curr,hapevent,sbpmavail,fbpmavail,ssbpmavail
	if pklredecode[1]:
	    if clkavail:
		try:pickle.dump(np.asarray(clkraw,dtype=np.uint32),open(pklpath[1]["clock"],"wb",-1),-1)
		except:sys.exit("\n\n\n\nError!failed to dump for clock data,do you have write permission in dir %s?"%pkldir)
	    if rasteravail:
		for i in range(4):
		    rasterraw[i]=np.asarray(rasterraw[i],dtype=np.int32)
		try:pickle.dump(rasterraw,open(pklpath[1]["raster"],"wb",-1),-1)
		except:sys.exit("\n\n\n\nError!failed to dump for rasterdata,do you have write permission in dir %s?"%pkldir)
	    try:pickle.dump(np.asarray(event,dtype=np.uint32),open(pklpath[1]["event"],"wb",-1),-1)
	    except:sys.exit("\n\n\n\nError!failed to dump for event data,do you have write permission in dir %s?"%pkldir)
	    del rasterraw,clkraw,event
	gc.collect() #recycle memory
    if buildtree:fillbpmrawtree(run,rootfilepath)
    return True

def fillbpmrawtree(run,rootpath,fileprefix="bpmraw"):
    print "filling bpm raw trees for run %i"%run   
    ##pkl file list
    pkldir=os.path.join(rootpath,"pkl")
    pklfilen=[["rbpm","sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail"],["raster","clock","event"]]
    datatypes=["bpm","raster"]
    for p in range(len(pklfilen)):
	for f in range(len(pklfilen[p])):pklfilen[p][f]="raw_"+pklfilen[p][f]
    for a in ["a","b"]:
	pklfilen[0].append("pos_sbpm%shall"%(a))
	pklfilen[0].append("pos_ssbpm%sbpm"%(a))
	for b in ["bpm","rot"]:
	    for c in ["s","f"]:
		pklfilen[0].append("pos_%sbpm%s%s"%(c,a,b))
    tgtpklfiles=glob.glob(os.path.join(pkldir,"bpmpos_tgt*_%i.pkl"%run))
    for p in tgtpklfiles:
	if os.path.getsize(p)>1000:
	    fn=os.path.split(p)[1]
	    pklfilen[0].append("pos_"+fn.split("_")[1])
    pklfilen[0].append("pos_rms")
    #get pkl file size and group them as total size 
    insertgroup=[]
    reverseram=3e8 #reverse memory,300mb
    maxsetram=2e9 #2gb memory maximum
    maxram=getmemory()
    maxram=maxram if maxram<maxsetram else maxsetram  
    maxram=maxram-reverseram
    for p in range(len(pklfilen)):
	insertgroup.append([[]])
	totalsize=0
	for f in range(len(pklfilen[p])):
	    fn=os.path.join(rootpath,"pkl/bpm%s_%i.pkl"%(pklfilen[p][f],run))
	    if not os.path.exists(fn):continue
	    totalsize+=os.path.getsize(fn)
	    if totalsize>=maxram:
		insertgroup[p].append([])
		totalsize=0
	    insertgroup[p][-1].append(pklfilen[p][f])
    #insert
    firstbranch=[True,True]
    firstopen=True
    for p in range(len(insertgroup)):
	for f in range(len(insertgroup[p])):
	    Nfile=len(insertgroup[p][f]) 
	    if Nfile<1:continue
	    if firstopen:
		bpmrootfile=TFile(os.path.join(rootpath,"%s_%i.root"%(fileprefix,run)),"RECREATE")
	    else:
		bpmrootfile=TFile(os.path.join(rootpath,"%s_%i.root"%(fileprefix,run)),"UPDATE")
	    print "filling %i-%i group data for run %i"%(p,f,run)
	    #create tree
	    if firstbranch[p]:
		tree=TTree(datatypes[p],datatypes[p])
	    else:
		tree=bpmrootfile.Get(datatypes[p])
	    #create branches
	    data,branch,Vdata=[0]*Nfile,[0]*Nfile,[0]*Nfile
	    numentry,numvar=[],[]
	    for pf in range(Nfile):
		fn=os.path.join(rootpath,"pkl/bpm%s_%i.pkl"%(insertgroup[p][f][pf],run))
		data[pf]=pickle.load(open(fn,"rb"))
		dataleaves=insertgroup[p][f][pf]
		if "-" in dataleaves:dataleaves=dataleaves.replace("-","m")
		branchname=dataleaves
		numdata=len(data[pf])
		if numdata>10:nvalues=1 #check if have more than 1 variables in a pkl file
		else:
		    nvalues=numdata
		    numdata=len(data[pf][0])
		numentry.append(numdata)
		numvar.append(nvalues)
		if nvalues<1:continue
		elif nvalues==1:dataleavessub=[""]
		elif nvalues==2:dataleavessub=["x","y"]
		elif nvalues==3:dataleavessub=["x","y","z"]
		elif nvalues==5:dataleavessub=["x","y","z","theta","phi"]
		else:dataleavessub=[str(i+1) for i in range(nvalues)]
		dataleaves=[dataleaves+a for a in dataleavessub]
		Vdata[pf]=array("f",[0.0]*nvalues)
		Vdataleaves=""
		for i in range(len(dataleaves)):
		    if i>0:Vdataleaves+=":%s/F"%dataleaves[i]
		    else:Vdataleaves+="%s/F"%dataleaves[i]
		branch[pf]=tree.Branch(branchname,Vdata[pf],Vdataleaves)
	    #fill
	    numdata=min(numentry)
	    for i in range(numdata):
		if not firstbranch[p]:tree.GetEntry(i)
		if i%10000==0:
		    print "filling %i-%i %i events,%i left"%(p,f,i,numdata-i)
		for pf in range(Nfile):
		    if numvar[pf]<2:
			Vdata[pf][0]=float(data[pf][i])
		    else:
			for j in range(numvar[pf]):
			    Vdata[pf][j]=float(data[pf][j][i])
		if firstbranch[p]:tree.Fill()
		else:
		    for pf in range(Nfile):branch[pf].Fill()
	    tree.Write("",TObject.kOverwrite)
	    bpmrootfile.Close()
	    if firstbranch[p]:firstbranch[p]=False
	    if firstopen:firstopen=False
	    del data,Vdata,tree,branch
	    gc.collect() #recycle memory

def getposrms(rms,rmspic="rms.png"):
    #from pylab import plot,show
    #plot(rms)
    #show()
    mineventsplit=1000
    gROOT.SetBatch(True)
    #fill rms tree
    rmstree=TTree("bpmrms","bpmrms")
    leaves=["rms"]
    Vleaves="rms/F"
    Vrms=array("f",[0.0])
    branch=rmstree.Branch("rms",Vrms,Vleaves)
    entries=len(rms)
    for i in range(entries):
	Vrms[0]=rms[i]
	rmstree.Fill()
    #find rms peaks
    s=TSpectrum()
    c1=TCanvas("c1","BPM rms",1024,768)
    rmstree.Draw("%s:Entry$"%leaves[0],"%s>0"%(leaves[0]))
    graph=gPad.GetPrimitive("Graph")
    meanrms=graph.GetRMS(2)
    arms=TProfile("arms","arms",3000,0,entries)
    rmstree.Draw("%s:Entry$>>arms"%leaves[0],"%s>0.08"%(leaves[0]),"same")
    nfound=s.Search(arms,10,"same",0.1)
    peakx=[]
    for j in range(nfound):
	peakx.append(int(s.GetPositionX()[j]))
    c1.Print(rmspic,"png")
    #find beam trip
    trippeaks,tripbackpeaks=[],[]
    oldtrip,oldtripback=0,0
    for i in range(len(rms)):
	if rms[i]<0:
	    if i-oldtrip>mineventsplit:trippeaks.append(i)
	    oldtrip=i
	else:
	    if i-oldtripback>mineventsplit:tripbackpeaks.append(i)
	    oldtripback=i
    #get rid of close trip
    trips=peakx+trippeaks+tripbackpeaks
    trips.append(0)
    trips.append(entries-1)
    trips.sort()
    finaltrips=[]
    for i in range(len(trips)):
	if i>0 and trips[i]-trips[i-1]<mineventsplit:continue
	else:finaltrips.append(trips[i])
    #split beam move
    splitentries=[]
    for i in range(len(finaltrips)-1):
	if i==0 and finaltrips[i+1] in tripbackpeaks:continue
	elif i==len(finaltrips)-2 and finaltrips[i] in trippeaks:break
	elif finaltrips[i] in trippeaks and finaltrips[i+1] in tripbackpeaks:
	    continue
	splitentries.append(finaltrips[i:i+2])
    return splitentries

def getposslowmove(pos,esplit):
    splitentries=[]
    mineventsmove=5000
    aveevents=1000
    tolepos=0.3 #min tolerant position
    for s in esplit:
	if s[1]-s[0]<mineventsmove*2:
	    splitentries.append(s)
	    continue
	avepos=np.mean(pos[s[0]:s[0]+mineventsmove])
	bpmave=Cavestack(aveevents)
	splitpoint=s[0]
	for e in range(s[0],s[1]):
	    if e>=len(pos):break
	    bpmave.push(pos[e])
	    if e<splitpoint+mineventsmove:continue
	    if abs(bpmave.ave()-avepos)>tolepos:
		avepos=bpmave.ave()
		splitentries.append([splitpoint,int(e-aveevents/2.)])
		splitpoint=int(e-aveevents/2.)
	splitentries.append([splitpoint,s[1]])
    return splitentries    
    
    