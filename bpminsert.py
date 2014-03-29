#!/usr/bin/env python
import os,re,shutil,numpy,sys,time,ctypes,bisect
try:import cPickle as pickle
except:import pickle
from urllib import urlretrieve
try:import ROOT
except:print "Error!! pyroot didn't compiled! please recompile your root!"
from array import array
from bcmconst import *
from harppos import bpmpos
from runinfo import runinfo,getpklpath,rasterconstread
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
    
def bpmbeamdriftprep(runorbit,tgtz):
    if runorbit<=0:return 0
    FourFloat=ctypes.c_float*4
    path=os.path.join(os.getenv("BEAMDBPATH"),"pyDB")
    posmap,map_x,map_y,map_theta,map_phi={},{},{},{},{}
    for z in tgtz:
	intz=int(z)
	try:
	    posmap[intz]=ctypes.CDLL(os.path.join(path,"bpm_%i_%i.so"%(runorbit,intz)))
	    map_x[intz]=posmap[intz].target_x
	    map_y[intz]=posmap[intz].target_y
	    map_theta[intz]=posmap[intz].target_theta
	    map_phi[intz]=posmap[intz].target_phi
	    map_x[intz].restype=ctypes.c_float
	    map_y[intz].restype=ctypes.c_float
	    map_theta[intz].restype=ctypes.c_float
	    map_phi[intz].restype=ctypes.c_float
	except Exception as e:
	    print e
	    print "sorry no position mapping funtion for orbit %i and target z %f,does %s file exists?\
	    if yes, please compile it; if no, please contact pengjia.will not calculate this target position now..."%(runorbit,z,os.path.join(path,"bpm_%i_%i.c"%(runorbit,z)))
	    posmap[intz]=False
	    map_x[intz],map_y[intz],map_theta[intz],map_phi[intz]=False,False,False,False
    return map_x,map_y,map_theta,map_phi,FourFloat

def bpminsertparaprep(runpath,treename="T",forcefastbus=False):
    rootfilepath,runfilename=os.path.split(runpath)
    run=int(re.split("[_]",re.split("[.]",runfilename)[-2])[-1])
    if run<100:run=int(re.split("[_]",re.split("[.]",os.path.split(runpath)[1])[-2])[-2])
    #get run info
    period=runinfo()
    runorbit=period.orbit(run)
    #get bpm calibration constant
    tmp=period.bpmconstread(run,forcefastbus)
    if not tmp["a"] or not tmp["b"]:
	print "sorry no constant found for run %i"%run
	return False
    pedestal=tmp["a"]["ped"]+tmp["b"]["ped"]
    bpmoffset=tmp["a"]["offset"]+tmp["b"]["offset"]
    tgtz=tmp["a"]["tgtz"]
    calconst={"a":tmp["a"]["const"],"b":tmp["b"]["const"]}
    fitorder={"a":tmp["a"]["fitorder"],"b":tmp["b"]["fitorder"]}
    #get raster calibration constant
    rasterconst=rasterconstread(run)
    #for non straight through,import mapping function
    posmapfun=bpmbeamdriftprep(runorbit,tgtz)
    para={"run":run,"ped":pedestal,"offset":bpmoffset,"tgtz":tgtz,"orbit":runorbit,\
	"calconst":calconst,"fitorder":fitorder,"posmapfun":posmapfun,"rasterconst":rasterconst}
    return para

def bpminsert(runpath,eventsinsert=0,treename="T",separatefile=0,forceredecode=0,forcefastbus=False):
    para=bpminsertparaprep(runpath,treename,forcefastbus)
    if not para:return False
    run,tgtz=para["run"],para["tgtz"]
    tgtz=[int(z) for z in tgtz]
    #get extra file
    rootfilepath,runfilename=os.path.split(runpath)
    pp=getpklpath(rootfilepath)
    info=runinfo()
    arm="L" if run<20000 else "R"
    rootfiles=glob.glob(os.path.join(rootfilepath,runfilename.replace(".root","_*.root")))
    rootfiles.append(runpath)
    rootfiles.sort()
    if len(rootfiles)<1:
	print "can not find file %s,will try to insert as separate file"%runpath
	separatefile=1
    #backup rootfile
    for runpath in rootfiles:
	print "back up rootfiles..."
	try:
	    if not os.path.exists(runpath+"_bak"):
		shutil.copy(runpath,runpath+"_bak")
	    else:shutil.copy(runpath+"_bak",runpath)
	except:
	    print "no file %s exists"%runpath
	    return False
    #scan rootfile events
    eventleaf="fEvtHdr.fEvtNum" #event number
    rootfileinfo={}
    firstevent,lastevent=-1,-1
    if separatefile==0:
	t1=time.time()
	for rf in rootfiles:
	    print "scan total events for rootfile %s"%rf
	    rootfile=ROOT.TFile(rf,"READ")
	    rootfileinfo[rf]={}
	    if rootfile.IsZombie():
		print "Error!! file %s is abnormal, please check if it is a valid rootfile! pass it..."%rf
		rootfileinfo[rf]["avail"]=False
		continue
	    try:
		tree=rootfile.Get(treename)
		events=tree.GetEntries()
		leventleaf=tree.GetLeaf(eventleaf)
	    except:
		print "Error!! file %s is abnormal, please check if it is a valid rootfile! pass it..."%rf
		rootfileinfo[rf]["avail"]=False
		continue
	    rootfileinfo[rf]["avail"]=True
	    rootfileinfo[rf]["events"]=events
	    evnums=numpy.zeros(events,dtype=numpy.int32)
	    for e in range(events):
		if e%10000==0:
		    print "scanning events %i %i"%(e,events-e)
		tree.GetEntry(e)
		evnums[e]=leventleaf.GetValue()
	    rootfileinfo[rf]["evnum"]=evnums
	    if firstevent<0 or evnums[0]<firstevent:
		firstevent=evnums[0]
	    if lastevent<0 or evnums[-1]>lastevent:
		lastevent=evnums[-1]
	    del evnums
	    rootfile.Close()
	if lastevent-firstevent<1000:
	    print "less than 1000 events found in your rootfile:",rootfiles,"will try to insert as separate"
	    separatefile=1
	gc.collect()
	t2=time.time()
	print "use time:",t2-t1
    #get raw data
    getposfromraw(para,runpath,treename,firstevent,lastevent,forceredecode=forceredecode,forcefastbus=forcefastbus) 
    if separatefile>0:
	fillbpmrawtree(run,rootfilepath,"bpmpos")
	return
    elif separatefile<0:return
    #sort events
    hapevent=pickle.load(open(pp.getpath("raw","hapevent",run),"rb"))
    eevent=pickle.load(open(pp.getpath("raw","event",run),"rb"))
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    rootfileinfo[rf]["eevoff"]=bisect.bisect_left(eevent,rootfileinfo[rf]["evnum"][0])
	    rootfileinfo[rf]["lhapev"]=numpy.searchsorted(hapevent,rootfileinfo[rf]["evnum"])
	    #rootfileinfo[rf]["rhapev"]=numpy.searchsorted(hapevent,rootfileinfo[rf]["evnum"],side="right")
	    rootfileinfo[rf]["lhapev"]=numpy.where(rootfileinfo[rf]["lhapev"]<len(hapevent),rootfileinfo[rf]["lhapev"],len(hapevent)-1)
	    #rootfileinfo[rf]["rhapev"]=numpy.where(rootfileinfo[rf]["rhapev"]<len(hapevent),rootfileinfo[rf]["rhapev"],len(hapevent)-1)
    #check raster on/off
    Iraster=[info.getsqlinfo(run,s) for s in ["FastRasterIx","FastRasterIy","SlowRasterIx","SlowRasterIy"]] 
    Iraster=[False if i==None or i==0 else True for i in Iraster]
    Iraster=[any(Iraster[:2]),any(Iraster[2:])]
    #raster const
    if any(Iraster):
	rasterconst=para["rasterconst"]
	rasterslope,rasterped=[{},{},{},{}],[{},{},{},{}]
	if rasterconst:
	    rasterconst["tgtz"]=[int(z) for z in rasterconst["tgtz"]]
	    for i in range(len(rasterconst["tgtz"])):
		z=rasterconst["tgtz"][i]
		rasterslope[0][z]=rasterconst["fstxslope"][i]
		rasterslope[1][z]=rasterconst["fstyslope"][i]
		rasterslope[2][z]=rasterconst["slxslope"][i]
		rasterslope[3][z]=rasterconst["slyslope"][i]
		rasterped[0][z]=rasterconst["fstxped"][i]
		rasterped[1][z]=rasterconst["fstyped"][i]
		rasterped[2][z]=rasterconst["slxped"][i]
		rasterped[3][z]=rasterconst["slyped"][i]
	else:
	    Iraster=[False,False]
    else:rasterconst=False
    #load pickle file
    bpmavail=pickle.load(open(pp.getpath("raw","bpmavail",run),"rb"))
    curr=pickle.load(open(pp.getpath("raw","curr",run),"rb"))
    if rasterconst:
	rasterraw=pickle.load(open(pp.getpath("raw","raster",run),"rb"))
	rastercenter=[(r.max()+r.min())/2. for r in rasterraw]
    tgtpos={}
    for z in tgtz:
	tgtpos[z]=pickle.load(open(pp.getpath("pos","tgt%i"%z,run),"rb"))
    #get pos for each rootfile
    print "get positions for each event for rootfile"
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    rootfileinfo[rf]["curr"]=curr[rootfileinfo[rf]["lhapev"]]
	    rootfileinfo[rf]["bpmavail"]=bpmavail[rootfileinfo[rf]["lhapev"]]
	    for z in tgtz:
		rootfileinfo[rf][z]=[0,0,0,0]
		rootfileinfo[rf][z][0]=tgtpos[z][0][rootfileinfo[rf]["lhapev"]]
		rootfileinfo[rf][z][1]=tgtpos[z][1][rootfileinfo[rf]["lhapev"]]
		rootfileinfo[rf][z][2]=tgtpos[z][3][rootfileinfo[rf]["lhapev"]]
		rootfileinfo[rf][z][3]=tgtpos[z][4][rootfileinfo[rf]["lhapev"]]
		if Iraster[0] and z in rasterconst["tgtz"]:
		    #for fast raster, x and y need to exchange
		    rootfileinfo[rf][z][0]=rootfileinfo[rf][z][0]+\
			(rasterraw[1][rootfileinfo[rf]["eevoff"]:rootfileinfo[rf]["eevoff"]+rootfileinfo[rf]["events"]]\
			    -rastercenter[1])*rasterslope[0][z]
		    rootfileinfo[rf][z][1]=rootfileinfo[rf][z][1]+\
			(rasterraw[0][rootfileinfo[rf]["eevoff"]:rootfileinfo[rf]["eevoff"]+rootfileinfo[rf]["events"]]\
			    -rastercenter[0])*rasterslope[1][z]
		if Iraster[1] and z in rasterconst["tgtz"]:
		    #for slow raster, x and y need to exchange
		    rootfileinfo[rf][z][0]=rootfileinfo[rf][z][0]+\
			(rasterraw[3][rootfileinfo[rf]["eevoff"]:rootfileinfo[rf]["eevoff"]+rootfileinfo[rf]["events"]]\
			    -rastercenter[3])*rasterslope[2][z]
		    rootfileinfo[rf][z][1]=rootfileinfo[rf][z][1]+\
			(rasterraw[2][rootfileinfo[rf]["eevoff"]:rootfileinfo[rf]["eevoff"]+rootfileinfo[rf]["events"]]\
			    -rastercenter[2])*rasterslope[3][z]
	    del rootfileinfo[rf]["lhapev"]#,rootfileinfo[rf]["rhapev"]
    del curr,bpmavail,tgtpos
    if rasterconst:del rasterraw
    gc.collect()
    #tgt leaves
    ltgt={}
    for z in tgtz:
	key=str(z).replace("-","m")
	ltgt[z]=[]
	for x in ["x","y","theta","phi"]:
	    ltgt[z].append(arm+"rb.tgt_%s_%s"%(key,x))
    #insert
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    rootfile=ROOT.TFile(runpath,"UPDATE")
	    tree=rootfile.Get(treename)
	    #branch
	    Vdata=[array("i",[0]),array("f",[0])]
	    branch=[tree.Branch(arm+"rb.bpmavail",Vdata[0],arm+"rb.bpmavail/I"),\
		tree.Branch(arm+"rb.curr",Vdata[1],arm+"rb.curr/F")]
	    nbranch=2
	    for z in tgtz:
		for x in range(4):
		    nbranch+=1
		    Vdata.append(array("f",[0]))
		    branch.append(tree.Branch(ltgt[z][x],Vdata[nbranch-1],ltgt[z][x]+"/F"))
	    print "inserting bpm info to file %s"%rf
	    tt1=time.time()
	    for e in range(rootfileinfo[rf]["events"]):
		if e%1000==0:print "filling %i events,%i left"%(e,rootfileinfo[rf]["events"]-e)
		Vdata[0][0]=rootfileinfo[rf]["bpmavail"][e]
		Vdata[1][0]=rootfileinfo[rf]["curr"][e]
		vd=1
		for z in tgtz:
		    for x in range(4):
			vd+=1
			Vdata[vd][0]=rootfileinfo[rf][z][x][e]
		for i in range(nbranch):
		    branch[i].Fill()
	    tree.Write("",ROOT.TObject.kOverwrite)
	    tt2=time.time()
	    print "use time:",tt2-tt1
	    rootfile.Close()		
		
def getposfromraw(para,runpath,treename,firstevent=-1,lastevent=-1,forceredecode=0,forcefastbus=False):
    ##get run info
    calconst={}
    run,pedestal,bpmoffset,tgtz,runorbit,calconst,fitorder,posmapfun=\
	para["run"],para["ped"],para["offset"],para["tgtz"],para["orbit"],para["calconst"],para["fitorder"],para["posmapfun"]
    if not calconst["a"]:
	print "sorry no constant found for BPM,please check"
	return False
    period=runinfo()
    rootfilepath,runfilename=os.path.split(runpath)
    if not period.ifhapavail(run) or forcefastbus:fastbus=True
    else:fastbus=False
    #pedestal for run, read from pedestal.pkl from database
    tmp=period.pedestal(run,fastbus)
    if tmp["a"]:
	for i in range(4):pedestal[i]=tmp["a"][i]
    if tmp["b"]:
	for i in range(4):pedestal[i+4]=tmp["b"][i]
    ##pkl file list
    pp=getpklpath(rootfilepath)
    rawprefix="raw"
    posprefix="pos"
    #for raw data
    pklpathn_raw=["sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail","raster","clock","event"]
    if fastbus:pklpathn_raw.append("fbbpm")
    #for pos at bpm
    pklpathn_bpm=[]
    for a in ["a","b"]:
	pklpathn_bpm.append("sbpm%shall"%(a))
	pklpathn_bpm.append("fbpm%shall"%(a))
	pklpathn_bpm.append("ssbpm%sbpm"%(a))
	for b in ["bpm","rot"]:
	    if fastbus:
		pklpathn_bpm.append("fbbpm%s%s"%(a,b))
	    for c in ["s","f"]:
		pklpathn_bpm.append("%sbpm%s%s"%(c,a,b))
    #for beam split
    pklpathn_rms="rms"
    pklpathn_split="split"
    ##get position at bpm
    if forcefastbus:
	if not os.path.exists(pp.getpath(rawprefix,"fbbpm",run)):    
	    forceredecode=1
    bpmposexists=False
    d=decode(runpath,treename,firstevent,lastevent,forceredecode,0,forcefastbus)
    d.autodecode()
    hapevent=pickle.load(open(pp.getpath(rawprefix,"hapevent",run),"rb"))
    rawdataentries=len(hapevent)
    #del hapevent
    #gc.collect() #recycle memory
    #check if events same as raw
    if os.path.exists(pp.getpath(posprefix,pklpathn_bpm[0],run)) and not forceredecode:
	data_bpm=pickle.load(open(pp.getpath(posprefix,pklpathn_bpm[0],run),"rb"))
	if len(data_bpm[0])==rawdataentries:
	    bpmposexists=True
	del data_bpm
	#gc.collect() #recycle memory
    if not bpmposexists:
	print "calculating position at bpm for run %i....................."%run
	#bpm initial
	bpminfo={"a":bpmpos(run,"a"),"b":bpmpos(run,"b")}
	for c in ["s","f","ss"]:
	    #get pure raw data
	    rawbpm=pickle.load(open(pp.getpath(rawprefix,"%sbpm"%c,run),"rb"))
	    for i in range(8):rawbpm[i]=rawbpm[i]-pedestal[i]-bpmoffset[i]
	    for a in ["a","b"]:
		ab=0 if a=="a" else 4
		#pos at bpm in bpm coordinate
		tt1=time.time()
		print "getting bpm %s position in bpm coordinate with %s filter"%(a,c)
		posbpm_bpm=getrealpos(rawbpm[ab:ab+4],calconst[a],fitorder[a])
		try:pickle.dump(posbpm_bpm,open(pp.getpath(posprefix,"%sbpm%sbpm"%(c,a),run,1),"wb",-1),-1)
		except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
		tt2=time.time()
		print tt2-tt1,"secs used"
		#pos at bpm in rot coordinate
		if c!="ss":
		    print "getting bpm %s position in rot coordinate with %s filter"%(a,c)
		    posbpm_rot=bpminfo[a].posbpmrotate(posbpm_bpm)
		    try:pickle.dump(posbpm_rot,open(pp.getpath(posprefix,"%sbpm%srot"%(c,a),run,1),"wb",-1),-1)
		    except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
		    tt3=time.time()
		    print tt3-tt2,"secs used"
		del posbpm_bpm
		#pos at bpm in hall coordinate
		if c!="ss":
		    print "getting bpm %s position in hall coordinate with %s filter"%(a,c)
		    posbpm_hall=bpminfo[a].posbpmrot2hall(posbpm_rot)
		    try:pickle.dump(posbpm_hall,open(pp.getpath(posprefix,"%sbpm%shall"%(c,a),run,1),"wb",-1),-1)
		    except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
		    del posbpm_hall
		    tt4=time.time()
		    print tt4-tt3,"secs used"
		if c!="ss":del posbpm_rot
		#gc.collect() #recycle memory
	#for fastbus bpm
	if fastbus:
	    print "calculating fastbus bpm position"
	    tt1=time.time()
	    fbrawbpm=pickle.load(open(pp.getpath(rawprefix,"fbbpm",run),"rb"))
	    for i in range(8):fbrawbpm[i]=fbrawbpm[i]-pedestal[i]-bpmoffset[i]
	    for a in ["a","b"]:
		ab=0 if a=="a" else 4
		posbpm_bpm=getrealpos(fbrawbpm[ab:ab+4],calconst[a],fitorder[a])
		fbposbpm_rot=bpminfo[a].posbpmrotate(fbposbpm_bpm)
		try:
		    pickle.dump(fbposbpm_bpm,open(pp.getpath(posprefix,"fbbpm%sbpm"%a,run,1),"wb",-1),-1)
		    pickle.dump(fbposbpm_rot,open(pp.getpath(posprefix,"fbbpm%srot"%a,run,1),"wb",-1),-1)
		except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
	    del fbrawbpm,fbposbpm_bpm,fbposbpm_rot
	    tt2=time.time()
	    print tt2-tt1,"secs used"
    #split beam move
    skipcal=False
    if os.path.exists(pp.getpath(posprefix,pklpathn_rms,run)) and not forceredecode:
	rms_pure=pickle.load(open(pp.getpath(posprefix,pklpathn_rms,run),"rb"))
	if len(rms_pure)==rawdataentries and os.path.exists(pp.getpath(posprefix,pklpathn_split,run)):
	    skipcal=True
	    beamsplit=pickle.load(open(pp.getpath(posprefix,pklpathn_split,run),"rb"))
	del rms_pure
	gc.collect()
    if not skipcal:
	posrms={}
	print "getting beam sharp move info for run %i........."%run
	tt10=time.time()
	aveevents,aveeventsbg=1000,300
	sbpmabpm=pickle.load(open(pp.getpath(posprefix,"sbpmabpm",run),"rb"))
	bpmavail=pickle.load(open(pp.getpath(rawprefix,"bpmavail",run),"rb"))
	
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
	try:pickle.dump(rms_pure,open(pp.getpath(posprefix,pklpathn_rms,run,1),"wb",-1),-1)
	except:sys.exit("\n\n\n\nError!failed to dump,do you have write permission in dir %s?"%rootfilepath)
	gc.collect() #recycle memory
	moveentries=getposrms(rms_pure)
	del rms_pure
	#gc.collect() #recycle memory
	#totposmove
	print "getting beam slow move info for run %i........."%run
	ssbpmabpm=pickle.load(open(pp.getpath(posprefix,"ssbpmabpm",run),"rb"))
	totpos=numpy.sqrt(ssbpmabpm[0]**2+ssbpmabpm[1]**2)
	moveentries=getposslowmove(totpos,moveentries)
	del ssbpmabpm,totpos
	gc.collect() #recycle memory
	#get corresponded events
	#hapevent=pickle.load(open(pklpath_raw["hapevent"],"rb"))
	moveevents=[map(lambda x:hapevent[x],e) for e in moveentries]
	#del hapevent
	#gc.collect() #recycle memory
	beamsplit={"splitevents":moveevents,"splitentries":moveentries}
	pickle.dump(beamsplit,open(pp.getpath(posprefix,pklpathn_split,run,1),"wb",-1))
	tt11=time.time()
	print tt11-tt10,"secs used"
    ##get position at target
    #check which z is not inserted
    print "checking if tgt position calculated......"
    notinsertedz=[]
    for z in tgtz:
	intz=int(z)
	fn="tgt%i"%(intz)
	if os.path.exists(pp.getpath(posprefix,fn,run)) and not forceredecode:
	    data_tgt=pickle.load(open(pp.getpath(posprefix,fn,run),"rb"))
	    lendata=len(data_tgt[0])
	    del data_tgt
	    #gc.collect() #recycle memory
	    if lendata==rawdataentries:continue
	if runorbit<=0 or posmapfun[0][intz]:
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
	    sbpmahall=pickle.load(open(pp.getpath(posprefix,"sbpmahall",run),"rb"))
	    sbpmbhall=pickle.load(open(pp.getpath(posprefix,"sbpmbhall",run),"rb"))
	    p=tgtpos_straight(sbpmahall,sbpmbhall,notinsertedz)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmahall,sbpmbhall,p
	    #gc.collect() #recycle memory
	else:
	    sbpmarot=pickle.load(open(pp.getpath(posprefix,"sbpmarot",run),"rb"))
	    sbpmbrot=pickle.load(open(pp.getpath(posprefix,"sbpmbrot",run),"rb"))
	    p=tgtpos_field(sbpmarot,sbpmbrot,notinsertedz,posmapfun)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmarot,sbpmbrot,p
	    #gc.collect() #recycle memory
	#dump target position data
	try:
	    for i in range(lennotz):
		pickle.dump(tgtpos[intz[i]],open(pp.getpath(posprefix,tgtfn[i],run,1),"wb",-1),-1)
	except:sys.exit("\n\n\n\nError!failed to dump,do you have write permission in dir %s?"%rootfilepath)
	tt9=time.time()
	print tt9-tt8,"secs used"
	#get average target position
	print "getting average positions..."
	for z in notinsertedz:
	    intz=int(z)
	    fn="tgt%i"%(intz)
	    beamsplit[fn]=[]
	    for i in range(len(beamsplit["splitentries"])):
		entryrange=beamsplit["splitentries"][i]
		avepos=[0,0,0,0,0]
		for k in range(5):
		    avepos[k]=numpy.mean(tgtpos[intz][k][entryrange[0]:entryrange[1]])
		beamsplit[fn].append(avepos)
	pickle.dump(beamsplit,open(pp.getpath(posprefix,pklpathn_split,run,1),"wb",-1))
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
	if isinstance(theta,numpy.float32):pos[intz].append(z)
	else:pos[intz].append(numpy.asarray([z]*len(theta),dtype=numpy.float32))
	pos[intz].append(theta)
	pos[intz].append(phi)
	#pos_err[intz].append(0)
	#pos_err[intz].append(theta_err)
	#pos_err[intz].append(phi_err)
    return pos

def tgtpos_field(bpma,bpmb,tgtz,posmapfun):
    #calculated by using bpm rot coordinate position
    #posmapfun:map_x,map_y,map_theta,map_phi,FourFloat
    try:lendata=len(bpma[0])
    except:lendata=0
    if lendata<1:
	pos={}
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
	    ###y###
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



