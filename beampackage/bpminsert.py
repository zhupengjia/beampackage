#!/usr/bin/env python
import os,re,shutil,numpy,sys,time,ctypes,bisect
from urllib import urlretrieve
try:import ROOT
except:print "Error!! pyroot didn't compiled! please recompile your root!"
from array import array
from bcmconst import *
from harppos import bpmpos
from runinfo import runinfo,getpklpath,rasterconstread,zload,zdump
from signalfilter import *

#remove the existed branch
def removebranch(tree,branchname):
    removeb=tree.GetBranch(branchname)
    removel=tree.GetLeaf(branchname)
    if removeb:
	tree.GetListOfBranches().Remove(removeb)
	tree.GetListOfBranches().Compress()
    if removel:
	tree.GetListOfLeaves().Remove(removel)
	tree.GetListOfLeaves().Compress()
    
#prepare the transported function from bpms to target, combine with *.c code
def bpmbeamdriftprep(runorbit,tgtz,path=os.path.join(os.getenv("BEAMDBPATH"),"pyDB")):
    if runorbit<=0:return 0
    FourFloat=ctypes.c_float*4
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

#prepare the parameters before insert
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
	pedestal,bpmoffset,tgtzc,calconst,fitorder=False,False,[0],{"a":False,"b":False},False
    else:
	pedestal=tmp["a"]["ped"]+tmp["b"]["ped"]
	bpmoffset=tmp["a"]["offset"]+tmp["b"]["offset"]
	tgtzc=tmp["a"]["tgtz"]
	fitorder={"a":tmp["a"]["fitorder"],"b":tmp["b"]["fitorder"]}
	calconst={"a":tmp["a"]["const"],"b":tmp["b"]["const"]}
    #check which target 
    tgtz=[]
    for z in period.gettgtz(run):
	if z in tgtzc:
	   tgtz.append(z)
    #get raster calibration constant
    rasterconst=rasterconstread(run)
    #for non straight through,import mapping function
    posmapfun=bpmbeamdriftprep(runorbit,tgtz)
    para={"run":run,"ped":pedestal,"offset":bpmoffset,"tgtz":tgtz,"orbit":runorbit,\
	"calconst":calconst,"fitorder":fitorder,"posmapfun":posmapfun,"rasterconst":rasterconst}
    return para

#insert the bpm info to existed rootfile
#additional should be a dict, like {leafname:pklfilename}, pklfile should be a numpy array and aligned with beam events, or {leafname:[varpklfile,evnumpklfile]}, varpklfile is a numpy array and aligned with events array in evnumpklfile.pklfile can also add id/key like pklfilename~evnum, means pklfilename['evnum'] will be read 
def bpminsert(runpath,eventsinsert=0,treename="T",separatefile=0,forceredecode=0,forcefastbus=False,backup=True,additional={}):
    constavail=True
    para=bpminsertparaprep(runpath,treename,forcefastbus)
    if not para["calconst"]["a"]:
	constavail=False
    run,tgtz=para["run"],para["tgtz"]
    tgtz=[int(z) for z in tgtz]
    #get extra file
    rootfilepath,runfilename=os.path.split(runpath)
    pp=getpklpath(rootfilepath)
    info=runinfo()
    current=info.current(run)
    targetoffset=info.getsqlinfo(run,"TargetOffsetmm")
    arm="L" if run<20000 else "R"
    rootfiles=glob.glob(os.path.join(rootfilepath,runfilename.replace(".root","_*.root")))
    rootfiles.append(runpath)
    rootfiles.sort()
    if len(rootfiles)<1:
	print "can not find file %s,will try to insert as separate file"%runpath
	separatefile=1
    #backup rootfile
    if separatefile==0:
	for runpath in rootfiles:
	    try:
		if not os.path.exists(runpath+"_bak"):
		    if backup:
			print "back up rootfile %s"%runpath
			shutil.copy(runpath,runpath+"_bak")
		else:
		    print "restore rootfile %s"%runpath
		    shutil.copy(runpath+"_bak",runpath)
	    except Exception as err:
		raise Exception("file %s abnormal"%runpath)
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
    if constavail:
	getposfromraw(para,runpath,treename,firstevent,lastevent,forceredecode=forceredecode,forcefastbus=forcefastbus) 
	if separatefile>0:
	    fillbpmrawtree(run,rootfilepath,"bpmpos")
	    return
	elif separatefile<0:return
    else:
	if separatefile!=0:return
    #load additional pkls
    def recpkl(v,keys): #recursion to get final value from var
	if len(keys)<1:return v
	elif len(keys)<2:
	    try:return v[keys[0]]
	    except:return None
	else:
	    return recpkl(v[keys[0]],keys[1:])
    def loadpkl(v):#load value in pkl's deep structure
	v=v.split("~")
	if len(v)<2:keys=[]
	else:keys=v[1:]
	if not os.path.exists(v[0]):return None
	return recpkl(zload(v[0]),keys)
    for k in additional.keys():
	if isinstance(additional[k],str):
	    additional[k]=loadpkl(additional[k])
	    if additional[k]==None:del additional[k]
	elif len(additional[k])<2:
	    additional[k]=loadpkl(additional[k][0])
	    if additional[k]==None:del additional[k]
	else:
	    additional[k]=[loadpkl(additional[k][i]) for i in range(2)]
	    if None in additional[k]:del additional[k]
    additionalkeys=additional.keys() 
    #sort events
    hapevent=zload(pp.getpath("raw","hapevent",run))
    eevent=zload(pp.getpath("raw","event",run))
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    #rootfileinfo[rf]["eevoff"]=bisect.bisect_left(eevent,rootfileinfo[rf]["evnum"][0]) event offset
	    rootfileinfo[rf]["leev"]=numpy.searchsorted(eevent,rootfileinfo[rf]["evnum"]) #sort event
	    #numpy.where use rootfileinfo if condition else len(eevent)
	    rootfileinfo[rf]["leev"]=numpy.where(rootfileinfo[rf]["leev"]<len(eevent),rootfileinfo[rf]["leev"],len(eevent)-1) 
	    rootfileinfo[rf]["lhapev"]=numpy.searchsorted(hapevent,rootfileinfo[rf]["evnum"]) #sort event
	    rootfileinfo[rf]["lhapev"]=numpy.where(rootfileinfo[rf]["lhapev"]<len(hapevent),rootfileinfo[rf]["lhapev"],len(hapevent)-1) #use rootfileinfo if condition else len(hapevent)
	    for k in additionalkeys:
		if len(additional[k])<5:#for the additional info that dodn't aligned with beam event
		    rootfileinfo[rf]["ev"+k]=numpy.searchsorted(additional[k][1],rootfileinfo[rf]["evnum"])
		    rootfileinfo[rf]["ev"+k]=numpy.where(rootfileinfo[rf]["ev"+k]<len(additional[k][1]),rootfileinfo[rf]["ev"+k],len(additional[k][1])-1)
    #check raster on/off
    Iraster=[info.getsqlinfo(run,s) for s in ["FastRasterIx","FastRasterIy","SlowRasterIx","SlowRasterIy"]]
    Iraster=[False if i==None or i==0 else True for i in Iraster]
    Iraster=[any(Iraster[:2]),any(Iraster[2:])]
    #raster const
    if any(Iraster):
	rasterconst=para["rasterconst"]
	rasterslope=[{},{},{},{}]
	if rasterconst:
	    rasterconst["tgtz"]=[int(z) for z in rasterconst["tgtz"]]
	    for i in range(len(rasterconst["tgtz"])):
		z=rasterconst["tgtz"][i]
		rasterslope[0][z]=rasterconst["fstxslope"][i]
		rasterslope[1][z]=rasterconst["fstyslope"][i]
		rasterslope[2][z]=rasterconst["slxslope"][i]
		rasterslope[3][z]=rasterconst["slyslope"][i]
	    #additional z
	    for z in tgtz:
		if not z in rasterconst["tgtz"]:
		    rasterconst["tgtz"].append(z)
		    for i in range(4):
			rasterslope[i][z]=rasterslope[i][rasterconst["tgtz"][0]]
	else:
	    Iraster=[False,False]
    else:rasterconst=False
    #load pickle file
    if constavail:
	bpmavail=zload(pp.getpath("raw","bpmavail",run))
	curr=zload(pp.getpath("raw","curr",run))
	if rasterconst:
	    rasterraw=zload(pp.getpath("raw","raster",run))
	    rastercenter=[(r.max()+r.min())/2. for r in rasterraw]
	tgtpos={}
	for z in tgtz:
	    tgtpos[z]=zload(pp.getpath("pos","tgt%i"%z,run))
    #get pos for each rootfile
    print "get positions for each event for rootfile"
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    if constavail:
		rootfileinfo[rf]["curr"]=curr[rootfileinfo[rf]["lhapev"]].astype(numpy.float32)
		rootfileinfo[rf]["bpmavail"]=bpmavail[rootfileinfo[rf]["lhapev"]].astype(numpy.int32)
	    else:
		rootfileinfo[rf]["curr"]=numpy.zeros(len(rootfileinfo[rf]["evnum"]),dtype=numpy.float32)+current
		rootfileinfo[rf]["bpmavail"]=numpy.zeros(len(rootfileinfo[rf]["evnum"]),dtype=numpy.int32)
	    for z in tgtz:
		rootfileinfo[rf][z]=[0,0,0,0,0,0]
		if constavail:
		    rootfileinfo[rf][z][0]=tgtpos[z][0][rootfileinfo[rf]["lhapev"]]
		    rootfileinfo[rf][z][1]=tgtpos[z][1][rootfileinfo[rf]["lhapev"]]
		    rootfileinfo[rf][z][2]=tgtpos[z][3][rootfileinfo[rf]["lhapev"]]
		    rootfileinfo[rf][z][3]=tgtpos[z][4][rootfileinfo[rf]["lhapev"]]
		    #for 2Hz pos
		    rootfileinfo[rf][z][4]=rootfileinfo[rf][z][0]
		    rootfileinfo[rf][z][5]=rootfileinfo[rf][z][1]
		    if Iraster[0] and z in rasterconst["tgtz"]:
			#for fast raster, x and y need to exchange
			rootfileinfo[rf][z][0]=rootfileinfo[rf][z][0]+(rasterraw[1][rootfileinfo[rf]["leev"]]-rastercenter[1])*rasterslope[0][z]
			rootfileinfo[rf][z][1]=rootfileinfo[rf][z][1]+(rasterraw[0][rootfileinfo[rf]["leev"]]-rastercenter[0])*rasterslope[1][z]
		    if Iraster[1] and z in rasterconst["tgtz"]:
			#for slow raster, x and y need to exchange
			rootfileinfo[rf][z][0]=rootfileinfo[rf][z][0]+(rasterraw[3][rootfileinfo[rf]["leev"]]-rastercenter[3])*rasterslope[2][z]
			rootfileinfo[rf][z][1]=rootfileinfo[rf][z][1]+(rasterraw[2][rootfileinfo[rf]["leev"]]-rastercenter[2])*rasterslope[3][z]
		    for i in range(6):
			rootfileinfo[rf][z][i]=rootfileinfo[rf][z][i].astype(numpy.float32)
		else:
		    for i in range(6):
			rootfileinfo[rf][z][i]=numpy.zeros(len(rootfileinfo[rf]["evnum"]),dtype=numpy.float32)
		    rootfileinfo[rf][z][1]+=targetoffset
		    rootfileinfo[rf][z][5]+=targetoffset
	    for k in additionalkeys:
		if len(additional[k])<5:
		    #rootfileinfo[rf][k]=additional[k][rootfileinfo[rf]["eevoff"]:rootfileinfo[rf]["eevoff"]+rootfileinfo[rf]["events"]]
		    rootfileinfo[rf][k]=additional[k][0][rootfileinfo[rf]["ev"+k]]
		else:
		    rootfileinfo[rf][k]=additional[k][rootfileinfo[rf]["leev"]]
		rootfileinfo[rf][k]=rootfileinfo[rf][k].astype(numpy.float32)
	    if constavail:del rootfileinfo[rf]["lhapev"],rootfileinfo[rf]["leev"]
    del additional
    if constavail:del curr,bpmavail,tgtpos
    if rasterconst:del rasterraw
    gc.collect()
    #tgt leaves
    ltgt={}
    for z in tgtz:
	key=str(z).replace("-","m")
	ltgt[z]=[]
	for x in ["x","y","theta","phi"]:
	    ltgt[z].append(arm+"rb.tgt_%s_%s"%(key,x))
	#for 2Hz pos
	for x in ["x","y"]:
	    ltgt[z].append(arm+"rb.tgtave_%s_%s"%(key,x))
    #insert
    for rf in rootfileinfo.keys():
	if rootfileinfo[rf]["avail"]:
	    rootfile=ROOT.TFile(rf,"UPDATE")
	    tree=rootfile.Get(treename)
	    #remove old branch
	    print "remove old branch"
	    removebranch(tree,arm+"rb.bpmavail")
	    removebranch(tree,arm+"rb.curr")
	    for z in tgtz:
		for x in range(6):#for 2Hz pos
		    removebranch(tree,ltgt[z][x])
	    for k in additionalkeys:
		removebranch(tree,k)
	    #branch
	    Vdata=[array("i",[0]),array("f",[0])]
	    branch=[tree.Branch(arm+"rb.bpmavail",Vdata[0],arm+"rb.bpmavail/I"),\
		tree.Branch(arm+"rb.curr",Vdata[1],arm+"rb.curr/F")]
	    nbranch=2
	    for z in tgtz:
		for x in range(6):#for 2Hz pos
		    nbranch+=1
		    Vdata.append(array("f",[0]))
		    branch.append(tree.Branch(ltgt[z][x],Vdata[nbranch-1],ltgt[z][x]+"/F"))
	    for k in additionalkeys:
		nbranch+=1
		Vdata.append(array("d",[0]))
		branch.append(tree.Branch(k,Vdata[nbranch-1],k+"/D"))
	    print "inserting bpm info to file %s"%rf
	    tt1=time.time()
	    for e in range(rootfileinfo[rf]["events"]):
		if e%1000==0:print "filling %i events,%i left"%(e,rootfileinfo[rf]["events"]-e)
		Vdata[0][0]=rootfileinfo[rf]["bpmavail"][e]
		Vdata[1][0]=rootfileinfo[rf]["curr"][e]
		vd=1
		for z in tgtz:
		    for x in range(6):#for 2Hz pos
			vd+=1
			Vdata[vd][0]=rootfileinfo[rf][z][x][e]
		for k in additionalkeys:
		    vd+=1
		    Vdata[vd][0]=rootfileinfo[rf][k][e]
		for i in range(nbranch):
		    branch[i].Fill()
	    tree.Write("",ROOT.TObject.kOverwrite)
	    tt2=time.time()
	    print "done for inserting file %s"%rf
	    print "use time:",tt2-tt1
	    rootfile.Close()	
	    del tree,branch
	
#calculate the positions at bpms and target
def getposfromraw(para,runpath,treename="T",firstevent=-1,lastevent=-1,forceredecode=0,forcefastbus=False):
    ##get run info
    calconst={}
    run,pedestal,bpmoffset,tgtz,runorbit,calconst,fitorder,posmapfun=\
	para["run"],para["ped"],para["offset"],para["tgtz"],para["orbit"],para["calconst"],para["fitorder"],para["posmapfun"]
    if not calconst["a"]:
	print "Error!!! sorry no constant found for BPM,please check"
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
    hapevent=zload(pp.getpath(rawprefix,"hapevent",run))
    rawdataentries=len(hapevent)
    #del hapevent
    #gc.collect() #recycle memory
    #check if events same as raw
    if os.path.exists(pp.getpath(posprefix,pklpathn_bpm[0],run)) and not forceredecode:
	data_bpm=zload(pp.getpath(posprefix,pklpathn_bpm[0],run))
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
	    rawbpm=zload(pp.getpath(rawprefix,"%sbpm"%c,run))
	    for i in range(8):rawbpm[i]=rawbpm[i]-pedestal[i]-bpmoffset[i]
	    for a in ["a","b"]:
		ab=0 if a=="a" else 4
		#pos at bpm in bpm coordinate
		tt1=time.time()
		print "getting bpm %s position in bpm coordinate with %s filter"%(a,c)
		posbpm_bpm=getrealpos(rawbpm[ab:ab+4],calconst[a],fitorder[a])
		try:zdump(posbpm_bpm,pp.getpath(posprefix,"%sbpm%sbpm"%(c,a),run,1))
		except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
		tt2=time.time()
		print tt2-tt1,"secs used"
		#pos at bpm in rot coordinate
		if c!="ss":
		    print "getting bpm %s position in rot coordinate with %s filter"%(a,c)
		    posbpm_rot=bpminfo[a].posbpmrotate(posbpm_bpm)
		    try:zdump(posbpm_rot,pp.getpath(posprefix,"%sbpm%srot"%(c,a),run,1))
		    except:sys.exit("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
		    tt3=time.time()
		    print tt3-tt2,"secs used"
		del posbpm_bpm
		#pos at bpm in hall coordinate
		if c!="ss":
		    print "getting bpm %s position in hall coordinate with %s filter"%(a,c)
		    posbpm_hall=bpminfo[a].posbpmrot2hall(posbpm_rot)
		    try:zdump(posbpm_hall,pp.getpath(posprefix,"%sbpm%shall"%(c,a),run,1))
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
	    fbrawbpm=zload(pp.getpath(rawprefix,"fbbpm",run))
	    for i in range(8):fbrawbpm[i]=fbrawbpm[i]-pedestal[i]-bpmoffset[i]
	    for a in ["a","b"]:
		ab=0 if a=="a" else 4
		fbposbpm_bpm=getrealpos(fbrawbpm[ab:ab+4],calconst[a],fitorder[a])
		fbposbpm_rot=bpminfo[a].posbpmrotate(fbposbpm_bpm)
		try:
		    zdump(fbposbpm_bpm,pp.getpath(posprefix,"fbbpm%sbpm"%a,run,1))
		    zdump(fbposbpm_rot,pp.getpath(posprefix,"fbbpm%srot"%a,run,1))
		except:raise Exception("Error!failed to dump,do you have write permission in dir %s?"%rootfilepath)
	    del fbrawbpm,fbposbpm_bpm,fbposbpm_rot
	    tt2=time.time()
	    print tt2-tt1,"secs used"
    #split beam move
    skipcal=False
    if os.path.exists(pp.getpath(posprefix,pklpathn_rms,run)) and not forceredecode:
	rms_pure=zload(pp.getpath(posprefix,pklpathn_rms,run))
	if len(rms_pure)==rawdataentries and os.path.exists(pp.getpath(posprefix,pklpathn_split,run)):
	    skipcal=True
	    beamsplit=zload(pp.getpath(posprefix,pklpathn_split,run))
	del rms_pure
	gc.collect()
    if not skipcal:
	posrms={}
	print "getting beam sharp move info for run %i........."%run
	tt10=time.time()
	aveevents,aveeventsbg=1000,300
	sbpmabpm=zload(pp.getpath(posprefix,"sbpmabpm",run))
	bpmavail=zload(pp.getpath(rawprefix,"bpmavail",run))
	
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
	try:zdump(rms_pure,pp.getpath(posprefix,pklpathn_rms,run,1))
	except:sys.exit("\n\n\n\nError!failed to dump,do you have write permission in dir %s?"%rootfilepath)
	gc.collect() #recycle memory
	moveentries=getposrms(rms_pure)
	del rms_pure
	#gc.collect() #recycle memory
	#totposmove
	print "getting beam slow move info for run %i........."%run
	ssbpmabpm=zload(pp.getpath(posprefix,"ssbpmabpm",run))
	totpos=numpy.sqrt(ssbpmabpm[0]**2+ssbpmabpm[1]**2)
	moveentries=getposslowmove(totpos,moveentries)
	del ssbpmabpm,totpos
	gc.collect() #recycle memory
	#get corresponded events
	#hapevent=zload(pklpath_raw["hapevent"])
	moveevents=[map(lambda x:hapevent[x],e) for e in moveentries]
	#del hapevent
	#gc.collect() #recycle memory
	beamsplit={"splitevents":moveevents,"splitentries":moveentries}
	zdump(beamsplit,pp.getpath(posprefix,pklpathn_split,run,1))
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
	    data_tgt=zload(pp.getpath(posprefix,fn,run))
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
	    sbpmahall=zload(pp.getpath(posprefix,"sbpmahall",run))
	    sbpmbhall=zload(pp.getpath(posprefix,"sbpmbhall",run))
	    p=tgtpos_straight(sbpmahall,sbpmbhall,notinsertedz)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmahall,sbpmbhall,p
	    #gc.collect() #recycle memory
	else:
	    sbpmarot=zload(pp.getpath(posprefix,"sbpmarot",run))
	    sbpmbrot=zload(pp.getpath(posprefix,"sbpmbrot",run))
	    p=tgtpos_field(sbpmarot,sbpmbrot,notinsertedz,posmapfun)
	    for z in intz:tgtpos[z]=p[z]
	    del sbpmarot,sbpmbrot,p
	    #gc.collect() #recycle memory
	#dump target position data
	try:
	    for i in range(lennotz):
		zdump(tgtpos[intz[i]],pp.getpath(posprefix,tgtfn[i],run,1))
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
	zdump(beamsplit,pp.getpath(posprefix,pklpathn_split,run,1))
	tt10=time.time()
	print tt10-tt9,"secs used"
	del tgtpos
	gc.collect() #recycle memory

#get the target position from two bpms with straight through setting
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

#get the target position from two bpms with target field setting
def tgtpos_field(bpma,bpmb,tgtz,posmapfun):
    #calculated by using bpm rot coordinate position
    #posmapfun:map_x,map_y,map_theta,map_phi,FourFloat
    try:lendata=len(bpma[0])
    except:lendata=0
    if lendata<1:
	pos={}
	x=posmapfun[4](bpma[0],bpma[1],bpmb[0],bpmb[1])
	#dx=posmapfun[4](bpma_err[0],bpma_err[1],bpmb_err[0],bpmb_err[1])
	for z in tgtz:
	    intz=int(z)
	    #pos_err[intz]=[]
	    ###x###
	    dx=posmapfun[4](0,0,0,0)
	    tgtx=numpy.float32(posmapfun[0][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    ###y###
	    dx=posmapfun[4](0,0,0,0)
	    tgty=numpy.float32(posmapfun[1][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    ###theta###
	    dx=posmapfun[4](0,0,0,0)
	    tgttheta=numpy.float32(posmapfun[2][intz](x,dx,4))
	    #pos_err[intz].append(dx[0])
	    ###phi###
	    dx=posmapfun[4](0,0,0,0)
	    tgtphi=numpy.float32(posmapfun[3][intz](x,dx,4))
	    pos[intz]=[tgtx,tgty,z,tgttheta,tgtphi] if isinstance(tgtx,numpy.float32) else [tgtx,tgty,numpy.asarray([z]*len(theta),dtype=numpy.float32),tgttheta,tgtphi]
	    #pos_err[intz].append(dx[0])
    else:
	lendata=len(bpma[0])
	tgtz=[int(z) for z in tgtz]
	pos={}
	posz=numpy.asarray([z]*lendata,dtype=numpy.float32)
	for z in tgtz:
	    pos[z]= [numpy.zeros(lendata,dtype=np.float32), numpy.zeros(lendata,dtype=np.float32), posz, numpy.zeros(lendata,dtype=np.float32), numpy.zeros(lendata,dtype=np.float32)]
	for i in range(lendata):
	    if i%10000==0:print "calculating position at target %i/%i"%(i,lendata)
	    x=posmapfun[4](bpma[0][i],bpma[1][i],bpmb[0][i],bpmb[1][i])
	    for z in tgtz:
		dx=posmapfun[4](0,0,0,0)
		pos[z][0][i]=posmapfun[0][z](x,dx,4)
		dx=posmapfun[4](0,0,0,0)
		pos[z][1][i]=posmapfun[1][z](x,dx,4)
		dx=posmapfun[4](0,0,0,0)
		pos[z][3][i]=posmapfun[2][z](x,dx,4)
		dx=posmapfun[4](0,0,0,0)
		pos[z][4][i]=posmapfun[3][z](x,dx,4)
    return pos



