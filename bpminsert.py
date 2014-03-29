#!/usr/bin/env python
import os,re,shutil,numpy,sys,time,ctypes
try:import cPickle as pickle
except:import pickle
from urllib import urlretrieve
from ROOT import gROOT,gPad,gDirectory,TCanvas,TFile,TObject,TTree
from array import array
from bpmcalib import bpmconstread
from bcmconst import *
from harppos import bpmpos
from runinfo import runinfo
from signalfilter import *

def removebranch(tree,branchname):
    removeb=tree.GetBranch(branchname)
    removel=tree.GetLeaf(branchname)
    if removeb:
	tree.GetListOfBranches().Remove(removeb)
	tree.GetListOfBranches().Compress()
    if removel:
	tree.GetListOfLeaves().Remove(removel)
	tree.GetListOfLeaves().Compress()

def getcurrfromraw(runpath,treename="T"):
    if not decodefromrootfile(runpath,treename):return False
    nocurr=0.002 #below this current will deal as no current
    rootpath,runfilename=os.path.split(runpath)
    run=int(re.split("[_]",re.split("[.]",runfilename)[-2])[-1])
    if run<100:run=int(re.split("[_]",re.split("[.]",os.path.split(runpath)[1])[-2])[-2])
    pklpath=os.path.join(rootpath,"pkl/bpmraw_curr_%i.pkl"%run)
    rawdata=pickle.load(open(pklpath,"rb"))
    curr=filter(lambda x:x>nocurr,rawdata)
    if len(curr)<len(rawdata)/50.:
	return 0 #if no current at 98% of time, will treat as no current
    avecurr=sum(curr)/len(curr)
    del rawdata
    return avecurr*1000
    
def bpmbeamdriftprep(runorbit,tgtz):
    if runorbit<=0:return 0
    FourFloat=ctypes.c_float*4
    path=os.path.join(os.getenv("BEAMDBPATH"),"pyDB")
    posmap,map_x,map_y,map_theta,map_phi={},{},{},{},{}
    for z in tgtz:
	intz=int(z)
	try:
	    posmap[intz]=ctypes.CDLL(os.path.join(path,"bpm_%i_%i.so"%(runorbit,intz)))
	except Exception as e:
	    print e
	    print "sorry no position mapping funtion for orbit %i and target z %f,does %s file exists?\
	    if yes, please compile it; if no, please contact pengjia"%(runorbit,z,os.path.join(path,"bpm_%i_%i.c"%(runorbit,z)))
	    sys.exit()
	    return False
	map_x[intz]=posmap[intz].target_x
	map_y[intz]=posmap[intz].target_y
	map_theta[intz]=posmap[intz].target_theta
	map_phi[intz]=posmap[intz].target_phi
	map_x[intz].restype=ctypes.c_float
	map_y[intz].restype=ctypes.c_float
	map_theta[intz].restype=ctypes.c_float
	map_phi[intz].restype=ctypes.c_float
    return map_x,map_y,map_theta,map_phi,FourFloat

def bpminsertparaprep(runpath,treename="T"):
    rootfilepath,runfilename=os.path.split(runpath)
    run=int(re.split("[_]",re.split("[.]",runfilename)[-2])[-1])
    if run<100:run=int(re.split("[_]",re.split("[.]",os.path.split(runpath)[1])[-2])[-2])
    #get run info
    period=runinfo()
    runorbit=period.orbit(run)
    bpmasurvey=period.bpmasurvey(run)
    bpmbsurvey=period.bpmbsurvey(run)
    survey=(bpmasurvey,bpmbsurvey)
    #get run current
    currepics=period.current(run)
    if not currepics:
	print "sorry no current info for run %i in database,will try to find it from rawdata first"%run
	currepics=getcurrfromraw(runpath,treename)
	if currepics:
	    dbdir=os.getenv("BEAMDBPATH")
	    if dbdir==None:
		print "please define BEAMDBPATH in your env"
		return False
	    pydb=os.path.join(dbdir,"pyDB")
	    currdb=os.path.join(pydb,"runcurr.pkl")
	    currinfo=pickle.load(open(currdb,"rb"))
	    currinfo[run]=currepics
	    print "checked run:%i, current:%f nA"%(run,currepics)
	    try:
		pickle.dump(currinfo,open(currdb,"wb",-1))
		print "updated currinfo database,please share with other people for this file %s or send it to pengjia so that other people don't need to run it again."%currdb
	    except:
		print "sorry can not update currinfo database, please check if you have permission to write in %s."%pydb
	else:
	    from select import select
	    print "sorry can not get curr info for the run %i,please input current for this run by hand. this window will disappear and will exit program after 20 secs.please input(nA):"%run
	    i,o,e =select([sys.stdin],[],[],20)
	    if i:
		try:currepics=float(sys.stdin.readline().strip())
		except:
		    sys.exit("oh no you put incorrect value,please input number next time,will exit now.")
	    else:
		sys.exit("don't have curr info for run %i,will exit insert"%run)
    print "average current is %inA"%currepics
    #get bpm calibration constant
    tmp=bpmconstread(run,currepics)
    if not tmp:return False
    pedestal=tmp["ped"]
    tgtz=tmp["tgtz"]
    calconsta=tmp["bpmaconst"]
    calconstb=tmp["bpmbconst"]
    #for non straight through,import mapping function
    posmapfun=bpmbeamdriftprep(runorbit,tgtz)
    para={"run":run,"ped":pedestal,"survey":survey,"tgtz":tgtz,"orbit":runorbit,\
	"calconsta":calconsta,"calconstb":calconstb,"posmapfun":posmapfun}
    return para

def bpminsert(runpath,eventsinsert=0,treename="T",separatefile=0,forceredecode=0):
    para=bpminsertparaprep(runpath,treename)
    if not para:return False
    run,tgtz=para["run"],para["tgtz"]
    #get extra file
    rootfilepath,runfilename=os.path.split(runpath)
    arm="L" if run<20000 else "R"
    rootfiles=glob.glob(os.path.join(rootfilepath,runfilename.replace(".root","_*.root")))
    rootfiles.append(runpath)
    rootfiles.sort()
    if len(rootfiles)<1:
	print "can not find file %s,will try to insert as separate file"%runpath
	separatefile=1
    eventleaf="fEvtHdr.fEvtNum" #event number
    #get total rootfile events
    firstevent,lastevent=-1,-1
    if not separatefile:
	print "getting total events for rootfiles"
	#first event
	for i in range(len(rootfiles)):
	    rootfile=TFile(rootfiles[i],"READ")
	    if rootfile.IsZombie():continue
	    try:
		tree=rootfile.Get(treename)
		tree.GetEntry(0)
		firstevent=tree.GetLeaf(eventleaf).GetValue()
		break
	    except:continue
	#last event
	for i in range(len(rootfiles)-1,-1,-1):
	    rootfile=TFile(rootfiles[i],"READ")
	    if rootfile.IsZombie():continue
	    try:
		tree=rootfile.Get(treename)
		tree.GetEntry(tree.GetEntries()-1)
		lastevent=tree.GetLeaf(eventleaf).GetValue()
		break
	    except:continue
	if lastevent-firstevent<1000:
	    print "less than 1000 events found in your rootfile:",rootfiles,"will try to insert as separate"
	    separatefile=1
    #get raw data
    getposfromraw(para,runpath,treename,firstevent,lastevent,forceredecode) 
    if separatefile:
	fillbpmrawtree(run,rootfilepath,"bpmpos")
	return
    #pkl file list
    pkldir=os.path.join(rootfilepath,"pkl")
    if not os.path.exists(pkldir):os.makedirs(pkldir)
    pklfilen,hapdata=["bpmavail","curr"],[True,True]
    for p in range(len(pklfilen)):
	pklfilen[p]="raw_"+pklfilen[p]
    for a in ["a","b"]:
	pklfilen.append("pos_sbpm%srot"%a)
	hapdata.append(True)  # it is from happex daq
    for z in tgtz:
	pklfilen.append("pos_tgt%i"%(int(z)))
	hapdata.append(True)
    pklfilen.append("raw_raster")
    hapdata.append(False)
    #get pkl file size and group them as total size 
    reverseram=3e8 #reverse memory,300mb
    maxsetram=2e9 #2gb memory maximum
    maxram=getmemory()
    maxram=maxram if maxram<maxsetram else maxsetram
    maxram=maxram-reverseram
    insertgroup=[[[],[]]]
    totalsize=0
    for p in range(len(pklfilen)):
	fn=os.path.join(rootfilepath,"pkl/bpm%s_%i.pkl"%(pklfilen[p],run))
	if not os.path.exists(fn):continue
	totalsize+=os.path.getsize(fn)
	if totalsize>=maxram:
	    insertgroup.append([[],[]])
	    totalsize=0
	#group as happex
	if hapdata[p]:insertgroup[-1][0].append(pklfilen[p])
	else:insertgroup[-1][1].append(pklfilen[p])
    #get raw event data
    pklhapevent=os.path.join(pkldir,"bpmraw_hapevent_%i.pkl"%run)
    pklevent=os.path.join(pkldir,"bpmraw_event_%i.pkl"%run)
    pklsplit=os.path.join(pkldir,"bpmpos_split_%i.pkl"%run)
    hapevent=pickle.load(open(pklhapevent,"rb"))
    eevent=pickle.load(open(pklevent,"rb"))
    rawhapentries=len(hapevent)
    rawentries=len(eevent)
    #for average position
    tgtleaves=[]
    tgtz.sort()
    tgtz=[int(z) for z in tgtz]
    for z in tgtz:
	    for x in ["x","y","z","theta","phi"]:
		tgtleaves.append(arm+"rb.tgtave_%s_%s"%(str(z).replace("-","m"),x))
    try:
	beamsplit=pickle.load(open(pklsplit,"rb"))
	aveinsert=True
    except:aveinsert=False
    #backup
    for runpath in rootfiles:
	print "back up rootfiles..."
	try:
	    if not os.path.exists(runpath+"_bak"):
		shutil.copy(runpath,runpath+"_bak")
	    else:shutil.copy(runpath+"_bak",runpath)
	except:
	    print "no file %s exists"%runpath
	    return False
    #insert 
    for p in range(len(insertgroup)):
	Nfile=[len(insertgroup[p][0]),len(insertgroup[p][1])]
	if Nfile<1:continue
	#load data
	data,dataleaves=[[0]*Nfile[0],[0]*Nfile[1]],[[0]*Nfile[0],[0]*Nfile[1]]
	numentry,numvar=[[],[]],[[],[]]
	for h in range(2):
	    for f in range(Nfile[h]):
		dataleaves[h][f]=insertgroup[p][h][f]
		fn=os.path.join(rootfilepath,"pkl/bpm%s_%i.pkl"%(dataleaves[h][f],run))
		data[h][f]=pickle.load(open(fn,"rb"))
		if "_" in dataleaves[h][f]:
		    dataleaves[h][f]=reduce(lambda x,y:x+"_"+y,dataleaves[h][f].split("_")[1:])
		if "tgt" in dataleaves[h][f]:
		    dataleaves[h][f]=dataleaves[h][f].replace("tgt","tgt_")
		    dataleaves[h][f]=dataleaves[h][f].replace("-","m")
		dataleaves[h][f]=arm+"rb."+dataleaves[h][f]
		numdata=len(data[h][f])
		if numdata>10:nvalues=1 #check if have more than 1 variables in a pkl file
		else:
		    nvalues=numdata
		    numdata=len(data[h][f][0])
		numentry[h].append(numdata)
		numvar[h].append(nvalues)
		if nvalues<1:continue
		elif nvalues==1:dataleavessub=[""]
		elif nvalues==2:dataleavessub=["x","y"]
		elif nvalues==3:dataleavessub=["x","y","z"]
		elif nvalues==5:dataleavessub=["x","y","z","theta","phi"]
		else:dataleavessub=[str(i+1) for i in range(nvalues)]
		if nvalues>1:dataleaves[h][f]=[dataleaves[h][f]+"_"+a for a in dataleavessub]
		else:dataleaves[h][f]=[dataleaves[h][f]+a for a in dataleavessub]
	#insert for each rootfiles
	eoffset=0
	for runpath in rootfiles:
	    print "insert group %i data for file %s..."%(p,runpath)
	    rootfile=TFile(runpath,"UPDATE")
	    if rootfile.IsZombie():
		print "Error!! file %s is abnormal, please check if it is a valid rootfile! pass it..."%runpath
		continue
	    try:
		tree=rootfile.Get(treename)
		events=tree.GetEntries()
		leventleaf=tree.GetLeaf(eventleaf)
	    except:
		print "Error!! file %s is abnormal, please check if it is a valid rootfile! pass it..."%runpath
		continue
	    #build branch
	    Vdata,branch=[[],[]],[[],[]]
	    for h in range(2):
		for f in range(Nfile[h]):
		    Vdata[h].append([])
		    branch[h].append([])
		    for j in range(numvar[h][f]):
			Vdata[h][f].append(array("f",[0.0]))
			branch[h][f].append(tree.Branch(dataleaves[h][f][j],Vdata[h][f][j],dataleaves[h][f][j]+"/F"))
	    tt1=time.time()
	    #for average position
	    if p==0 and aveinsert:
		aveVdata,avebranch=[],[]
		splitpos=0
		avenvalues=len(tgtleaves)
		nsplits=len(beamsplit["splitevents"])
		for i in range(avenvalues):
		    aveVdata.append(array("f",[1.0]))
		    avebranch.append(tree.Branch(tgtleaves[i],aveVdata[i],tgtleaves[i]+"/F"))
	    else:aveinsert=False
	    #fill
	    for e in range(events):
		tree.GetEntry(e)
		if e%1000==0:
		    print "filling group %i %i events,%i left"%(p,e,events-i)
		if e>eventsinsert and eventsinsert>0:break
		#get event number for this entry
		try:rootevent=leventleaf.GetValue() 
		except:
		    print "Error!! no leaf %s in your rootfile!!"%eventleaf
		    return False
		#compare event
		#for data in happex
		#get first matched event
		if e+eoffset<rawhapentries and Nfile[0]>0:
		    rawhapevent=hapevent[e+eoffset]
		    if rawhapevent>=rootevent:
			for i in range(-1,-50,-1):
			    if e+eoffset+i<0:break
			    rawhapevent=hapevent[e+eoffset+i]
			    if rootevent>rawhapevent:
				eoffset=eoffset+i+1
				break
		    else:
			for i in range(1,50):
			    if e+eoffset+i>=rawhapentries:break
			    rawhapevent=hapevent[e+eoffset+i]
			    if rootevent<=rawhapevent:
				eoffset=eoffset+i
				break
		    #get last matched event
		    if hapevent[e+eoffset]==rootevent:
			numring=1
			for i in range(1,50):
			    if e+eoffset+i>=rawhapentries:break
			    rawhapevent=hapevent[e+eoffset+i]
			    if rawhapevent==rootevent:numring+=1
			    else:break
		    else:numring=0
		    #insert happex data
		    if numring>0:
			for f in range(Nfile[0]):
			    if numvar[0][f]<2:
				if "bpmavail" in dataleaves[0][f]:
				    Vdata[0][f][0][0]=all(data[0][f][e+eoffset:e+eoffset+numring])
				else:
				    Vdata[0][f][0][0]=float(numpy.mean(data[0][f][e+eoffset:e+eoffset+numring]))
			    else:
				for i in range(numvar[0][f]):
				    Vdata[0][f][i][0]=float(numpy.mean(data[0][f][i][e+eoffset:e+eoffset+numring]))
		    eoffset=eoffset+numring
		#fill branch
		for h in range(2):
		    for f in range(Nfile[h]):
			for i in range(numvar[h][f]):
			    branch[h][f][i].Fill()
		#average position insert
		if aveinsert:
		    if splitpos>=nsplits:
			for i in range(avenvalues):
			    aveVdata[i][0]=-1000
			    avebranch[i].Fill()
			continue
		    if rootevent>=beamsplit["splitevents"][splitpos][1]:
			splitpos+=1
			if splitpos>=nsplits:
			    for i in range(nvalues):
				aveVdata[i][0]=-1000
				avebranch[i].Fill()
			    continue
		    if rootevent<beamsplit["splitevents"][splitpos][0]:
			for i in range(nvalues):
			    aveVdata[i][0]=-1000
			    avebranch[i].Fill()
			continue
		    k=0
		    for z in tgtz:
			for i in range(5):
			    aveVdata[k][0]=beamsplit["tgt%i"%z][splitpos][i]
			    avebranch[k].Fill()
			    k+=1
	    tree.Write("",TObject.kOverwrite)
	    tt2=time.time()
	    print "total insert time: %i secs, %i events per sec"%(tt2-tt1,events/(tt2-tt1))
	    rootfile.Close()
	    del data,Vdata,dataleaves,tree,branch
	    gc.collect() #recycle memory
		
def getposfromraw(para,runpath,treename,firstevent=-1,lastevent=-1,forceredecode=0):
    ##get run info
    calconst={}
    run,pedestal,survey,tgtz,runorbit,calconst["a"],calconst["b"],posmapfun=para["run"],para["ped"],\
	para["survey"],para["tgtz"],para["orbit"],para["calconsta"],para["calconstb"],\
	    para["posmapfun"]
    if not calconst["a"]:
	print "sorry no constant found for BPM,please check"
	return False
    rootfilepath,runfilename=os.path.split(runpath)
    ##pkl file list
    pkldir=os.path.join(rootfilepath,"pkl")
    if not os.path.exists(pkldir):
	try:os.makedirs(pkldir)
	except:sys.exit("Error!! failed to make dir in %s,please make sure you have write permission in there"%rootfilepath)
    #for raw data
    pklpathn_raw=["sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail","raster","clock","event"]
    pklpath_raw,pklpath_bpm,pklpath_tgt={},{},{}
    for n in pklpathn_raw:pklpath_raw[n]=os.path.join(pkldir,"bpmraw_%s_%i.pkl"%(n,run))
    #for pos at bpm
    pklpathn_bpm,pklpathn_tgt=[],[]
    for a in ["a","b"]:
	pklpathn_bpm.append("sbpm%shall"%(a))
	pklpath_bpm[pklpathn_bpm[-1]]=os.path.join(pkldir,"bpmpos_%s_%i.pkl"%(pklpathn_bpm[-1],run))
	pklpathn_bpm.append("ssbpm%sbpm"%(a))
	pklpath_bpm[pklpathn_bpm[-1]]=os.path.join(pkldir,"bpmpos_%s_%i.pkl"%(pklpathn_bpm[-1],run))
	for b in ["bpm","rot"]:
	    for c in ["s","f"]:
		pklpathn_bpm.append("%sbpm%s%s"%(c,a,b))
		pklpath_bpm[pklpathn_bpm[-1]]=os.path.join(pkldir,"bpmpos_%s_%i.pkl"%(pklpathn_bpm[-1],run))
    #for pos at tgt
    for z in tgtz:
	pklpathn_tgt.append("tgt%i"%(int(z)))
	pklpath_tgt[pklpathn_tgt[-1]]=os.path.join(pkldir,"bpmpos_%s_%i.pkl"%(pklpathn_tgt[-1],run))
    #for beam split
    pklpath_rms=os.path.join(pkldir,"bpmpos_rms_%i.pkl"%run)
    pklpath_split=os.path.join(pkldir,"bpmpos_split_%i.pkl"%run)
    
    ##get position at bpm
    bpmposexists=False
    #get raw data
    if not decodefromrootfile(runpath,treename,firstevent,lastevent,forceredecode):sys.exit()
    hapevent=pickle.load(open(pklpath_raw["hapevent"],"rb"))
    rawdataentries=len(hapevent)
    #del hapevent
    #gc.collect() #recycle memory
    #check if events same as raw
    if os.path.exists(pklpath_bpm[pklpathn_bpm[0]]) and not forceredecode:
	data_bpm=pickle.load(open(pklpath_bpm[pklpathn_bpm[0]],"rb"))
	if len(data_bpm[0])==rawdataentries:
	    bpmposexists=True
	del data_bpm
	#gc.collect() #recycle memory
    if not bpmposexists:
	print "calculating position at bpm for run %i....................."%run
	#bpm initial
	xypmorder={"a":[2,3,0,1],"b":[6,7,4,5]}
	bpminfo={"a":bpmpos(survey[0]),"b":bpmpos(survey[1])}
	for c in ["s","f","ss"]:
	    #get pure raw data
	    rawbpm=pickle.load(open(pklpath_raw["%sbpm"%c],"rb"))
	    for i in range(8):rawbpm[i]=rawbpm[i]-pedestal[i]
	    for a in ["a","b"]:
		#pos at bpm in bpm coordinate
		tt1=time.time()
		print "getting bpm %s position in bpm coordinate with %s filter"%(a,c)
		posbpm_bpm=getrealpos([[rawbpm[xypmorder[a][0]],rawbpm[xypmorder[a][1]]],\
		    [rawbpm[xypmorder[a][2]],rawbpm[xypmorder[a][3]]]],calconst[a])
		try:pickle.dump(posbpm_bpm,open(pklpath_bpm["%sbpm%sbpm"%(c,a)],"wb",-1),-1)
		except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%pkldir)
		tt2=time.time()
		print tt2-tt1,"secs used"
		#pos at bpm in rot coordinate
		if c!="ss":
		    print "getting bpm %s position in rot coordinate with %s filter"%(a,c)
		    posbpm_rot=bpminfo[a].posbpmrotate(posbpm_bpm)
		    try:pickle.dump(posbpm_rot,open(pklpath_bpm["%sbpm%srot"%(c,a)],"wb",-1),-1)
		    except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%pkldir)
		    tt3=time.time()
		    print tt3-tt2,"secs used"
		del posbpm_bpm
		#pos at bpm in hall coordinate
		if c=="s":
		    print "getting bpm %s position in hall coordinate with %s filter"%(a,c)
		    posbpm_hall=bpminfo[a].posbpmrot2hall(posbpm_rot)
		    try:pickle.dump(posbpm_hall,open(pklpath_bpm["%sbpm%shall"%(c,a)],"wb",-1),-1)
		    except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%pkldir)
		    del posbpm_hall
		    tt4=time.time()
		    print tt4-tt3,"secs used"
		if c!="ss":del posbpm_rot
		#gc.collect() #recycle memory
    #split beam move
    skipcal=False
    if os.path.exists(pklpath_rms) and not forceredecode:
	rms_pure=pickle.load(open(pklpath_rms,"rb"))
	if len(rms_pure)==rawdataentries and os.path.exists(pklpath_split):
	    skipcal=True
	    beamsplit=pickle.load(open(pklpath_split,"rb"))
	del rms_pure
	gc.collect()
    if not skipcal:
	posrms={}
	print "getting beam sharp move info for run %i........."%run
	tt10=time.time()
	aveevents,aveeventsbg=1000,300
	sbpmabpm=pickle.load(open(pklpath_bpm["sbpmabpm"],"rb"))
	bpmavail=pickle.load(open(pklpath_raw["bpmavail"],"rb"))
	rawevents=len(sbpmabpm[0])
	rms_pure=numpy.zeros(rawevents,dtype=numpy.float32)
	Vbpmave,Vbpmavebg=Cavestack(aveevents),Cavestack(aveeventsbg)
	for i in range(rawevents):
	    if i%10000==0:print "generated %i rms events, %i events left"%(i,rawevents-i)
	    if not bpmavail[i]:
		rms_pure[i]=-0.03
		continue
	    Vbpmave.push(sbpmabpm[0][i],sbpmabpm[1][i])
	    Vbpmavebg.push(sbpmabpm[0][i],sbpmabpm[1][i])
	    rms1,rms2=Vbpmave.rms(0),Vbpmave.rms(1)
	    rmsbg1,rmsbg2=Vbpmavebg.rms(0),Vbpmavebg.rms(1)
	    rms_pure[i]=abs(rms1-rmsbg1)+abs(rms2-rmsbg2)
	del sbpmabpm,bpmavail,Vbpmave,Vbpmavebg
	try:pickle.dump(rms_pure,open(pklpath_rms,"wb",-1),-1)
	except:sys.exit("\n\n\n\nError!failed to dump,do you have write permission in dir %s?"%pkldir)
	gc.collect() #recycle memory
	moveentries=getposrms(rms_pure)
	del rms_pure
	#gc.collect() #recycle memory
	#totposmove
	print "getting beam slow move info for run %i........."%run
	ssbpmabpm=pickle.load(open(pklpath_bpm["ssbpmabpm"],"rb"))
	totpos=numpy.sqrt(ssbpmabpm[0]**2+ssbpmabpm[1]**2)
	moveentries=getposslowmove(totpos,moveentries)
	del ssbpmabpm,totpos
	gc.collect() #recycle memory
	#get corresponded events
	#hapevent=pickle.load(open(pklpath_raw["hapevent"],"rb"))
	moveevents=[map(lambda x:hapevent[x],e) for e in moveentries]
	#del hapevent
	#gc.collect() #recycle memory
	beamsplit={"splitevents":moveevents,"splithapentries":moveentries}
	pickle.dump(beamsplit,open(pklpath_split,"wb",-1))
	tt11=time.time()
	print tt11-tt10,"secs used"
    ##get position at target
    #check which z is not inserted
    print "checking if tgt position calculated......"
    notinsertedz=[]
    for z in tgtz:
	intz=int(z)
	fn="tgt%i"%(intz)
	if os.path.exists(pklpath_tgt[fn]) and not forceredecode:
	    data_tgt=pickle.load(open(pklpath_tgt[fn],"rb"))
	    lendata=len(data_tgt[0])
	    del data_tgt
	    #gc.collect() #recycle memory
	    if lendata==rawdataentries:continue
	notinsertedz.append(z)
    lennotz=len(notinsertedz)
    if lennotz>0:
	#calculate target position
	print "calculating position at target for run %i........."%run
	tt8=time.time()
	tgtpos={}
	tgtfn,intz=[],[]
	for z in notinsertedz:
	    intz.append(int(z))
	    tgtfn.append("tgt%i"%(intz[-1]))
	if runorbit<=0:
	    sbpmahall=pickle.load(open(pklpath_bpm["sbpmahall"],"rb"))
	    sbpmbhall=pickle.load(open(pklpath_bpm["sbpmbhall"],"rb"))
	    p=tgtpos_straight(sbpmahall,sbpmbhall,notinsertedz)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmahall,sbpmbhall,p
	    #gc.collect() #recycle memory
	else:
	    sbpmarot=pickle.load(open(pklpath_bpm["sbpmarot"],"rb"))
	    sbpmbrot=pickle.load(open(pklpath_bpm["sbpmbrot"],"rb"))
	    p=tgtpos_field(sbpmarot,sbpmbrot,notinsertedz,posmapfun)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmarot,sbpmbrot,p
	    #gc.collect() #recycle memory
	#dump target position data
	try:
	    for i in range(lennotz):
		pickle.dump(tgtpos[intz[i]],open(pklpath_tgt[tgtfn[i]],"wb",-1),-1)
	except:sys.exit("\n\n\n\nError!failed to dump,do you have write permission in dir %s?"%pkldir)
	tt9=time.time()
	print tt9-tt8,"secs used"
	#get average target position
	print "getting average positions..."
	for z in notinsertedz:
	    intz=int(z)
	    fn="tgt%i"%(intz)
	    beamsplit[fn]=[]
	    for i in range(len(beamsplit["splithapentries"])):
		entryrange=beamsplit["splithapentries"][i]
		avepos=[0,0,0,0,0]
		for k in range(5):
		    avepos[k]=numpy.mean(tgtpos[intz][k][entryrange[0]:entryrange[1]])
		beamsplit[fn].append(avepos)
	pickle.dump(beamsplit,open(pklpath_split,"wb",-1))
	tt10=time.time()
	print tt10-tt9,"secs used"
	del tgtpos
	gc.collect() #recycle memory

def tgtpos_straight(bpma,bpmb,tgtz):
    #calculated by using bpm hall coordinate position
    ab=map(lambda a,b:b-a,bpma,bpmb)
    theta=numpy.float32(numpy.arctan2(ab[1],ab[2]))
    phi=numpy.float32(numpy.arctan2(ab[0],ab[2]))
    #theta_err=sqrt(bpma_err[1]**2+bpmb_err[1]**2)/(ab[2]+ab[1]**2/ab[2])
    #phi_err=sqrt(bpma_err[0]**2+bpmb_err[0]**2)/(ab[2]+ab[0]**2/ab[2])
    pos={}
    for z in tgtz:
	intz=int(z)
	pos[intz]=[]
	#pos_err[intz]=[]
	for j in range(2):pos[intz].append(bpma[j]+(z-bpma[2])/ab[2]*(bpmb[j]-bpma[j]))
	#for j in range(2):pos_err[intz].append(sqrt(bpma_err[j]**2+((z-bpma[2])/ab[2])**2*(bpmb_err[j]**2+bpma_err[j]**2)))
	pos[intz].append(z) if isinstance(theta,numpy.float32) else pos[intz].append(numpy.asarray([z]*len(theta),dtype=numpy.float32))
	pos[intz].append(theta)
	pos[intz].append(phi)
	#pos_err[intz].append(0)
	#pos_err[intz].append(theta_err)
	#pos_err[intz].append(phi_err)
 #   print pos,"\n",pos_err,"\n\n\n"
    return pos

def tgtpos_field(bpma,bpmb,tgtz,posmapfun):
    #calculated by using bpm rot coordinate position
    #posmapfun:map_x,map_y,map_theta,map_phi,FourFloat
    if isinstance(bpma[0],numpy.float32):   
	pos,pos_err={},{}
	x=posmapfun[4](bpma[0],bpma[1],bpmb[0],bpmb[1])
	#dx=posmapfun[4](bpma_err[0],bpma_err[1],bpmb_err[0],bpmb_err[1])
	dx=posmapfun[4](0,0,0,0)
	dxsave=dx
	for z in tgtz:
	    intz=int(z)
	    #pos_err[intz]=[]
	    ###x###
	    tgtx=numpy.float32(posmapfun[0][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    dx=dxsave
	    ###y###survey
	    tgty=numpy.float32(posmapfun[1][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    dx=dxsave
	    ###theta###
	    tgttheta=numpy.float32(posmapfun[2][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    dx=dxsave
	    ###phi###
	    tgtphi=numpy.float32(posmapfun[3][intz](x,dx,4))
	    pos[intz]=[tgtx,tgty,z,tgttheta,tgtphi] if isinstance(tgtx,numpy.float32) else [tgtx,tgty,numpy.asarray([z]*len(theta),dtype=numpy.float32),tgttheta,tgtphi]
	    #pos_err[intz].append(dx[0])
	    dx=dxsave
    else:
	lendata=len(bpma[0])
	tgtz=[int(z) for z in tgtz]
	pos={}
	posz=numpy.asarray([z]*lendata,dtype=numpy.float32)
	for z in tgtz:
	    pos[z]= [numpy.zeros(lendata,dtype=np.float32), numpy.zeros(lendata,dtype=np.float32), posz, numpy.zeros(lendata,dtype=np.float32), numpy.zeros(lendata,dtype=np.float32)]
	dx=posmapfun[4](0,0,0,0)
	for i in range(lendata):
	    if i%10000==0:print "calculating position at target %i/%i"%(i,lendata)
	    x=posmapfun[4](bpma[0][i],bpma[1][i],bpmb[0][i],bpmb[1][i])
	    for z in tgtz:
		pos[z][0][i]=posmapfun[0][z](x,dx,4)
		pos[z][1][i]=posmapfun[1][z](x,dx,4)
		pos[z][3][i]=posmapfun[2][z](x,dx,4)
		pos[z][4][i]=posmapfun[3][z](x,dx,4)
    return pos

def getrealpos(raw,calconst):
    ar,gx,gy=calconst[0]
    px=numpy.matrix(calconst[1]).T
    py=numpy.matrix(calconst[2]).T
    xdiff_sum=(raw[0][0]-gx*raw[0][1])/(raw[0][0]+gx*raw[0][1])
    ydiff_sum=(raw[1][0]-gy*raw[1][1])/(raw[1][0]+gy*raw[1][1])
    xbyb2=xdiff_sum**2+ydiff_sum**2
    xb2x=1./xbyb2-1./numpy.sqrt(xbyb2)*numpy.sqrt(1./xbyb2-1)
    xdata=ar*xdiff_sum*xb2x
    ydata=ar*ydiff_sum*xb2x
    M=numpy.matrix([xdata,ydata,1]) if calconst[3]<5 else numpy.matrix([xdata,ydata,1,xdata*xdata,ydata*ydata])
    x=numpy.array(M*px)[0].tolist()[0]
    y=numpy.array(M*py)[0].tolist()[0]
    #theta=atan2(ydiff_sum,xdiff_sum)
    #rho=ar*(1/numpy.sqrt(xbyb2)-numpy.sqrt(1/xbyb2-1))
    #xn=rho*cos(theta)
    #yn=rho*sin(theta)
    return [x,y]


