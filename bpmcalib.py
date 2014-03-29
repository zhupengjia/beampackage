#!/usr/bin/env python
import os,re,pickle,shutil,glob
from math import isnan
from pacupdate import *
from harppos import *
from bpmfit import *
from bpmpeak import bpm_peak
from runinfo import runinfo,harpinfo

def getmidpoint(pos):
    r=map(lambda p:p[0]**2+p[1]**2,pos)
    return r.index(min(r))

def filternan(a):
    b=[x for x in a if not math.isnan(x)]
    b,nanindex=[],[]
    for i in range(len(a)):
	if math.isnan(a[i]):nanindex.append(i)
	else:b.append(a[i])
    return b,nanindex

def reducelist(a,index):
    b=[]
    for i in range(len(a)):
	if i in index:continue
	else:b.append(a[i])
    return b

def calibrateonebpm(gxy,ar,pos,peak,err,ped):
    #purepeak:1st level:x,y;2nd level:n pos;3rd level:x+,x-
    #pos:1st level:n pos;2nd level:x,y
    purepeak=map(lambda p1,p2:map(lambda p3:map(lambda p4,p5:float(p4-p5),p3,p2),p1),peak,ped)
    xdiff_sum=map(lambda p:(p[0]-gxy[0]*p[1])/(p[0]+gxy[0]*p[1]),purepeak[0])
    ydiff_sum=map(lambda p:(p[0]-gxy[1]*p[1])/(p[0]+gxy[1]*p[1]),purepeak[1])
    xbyb2=map(lambda p1,p2:p1**2+p2**2,xdiff_sum,ydiff_sum)
    xb2x,nanindex=filternan(map(lambda p:1/p-1/sqrt(p)*sqrt(1/p-1),xbyb2))
    xdiff_sum=reducelist(xdiff_sum,nanindex)
    ydiff_sum=reducelist(ydiff_sum,nanindex)
    xdata=map(lambda p1,p2:ar*p1*p2,xdiff_sum,xb2x)
    ydata=map(lambda p1,p2:ar*p1*p2,ydiff_sum,xb2x)
    xharp=reducelist(map(lambda p:p[0],pos),nanindex)
    yharp=reducelist(map(lambda p:p[1],pos),nanindex)
    #err
    purepeak[0]=reducelist(purepeak[0],nanindex)
    purepeak[1]=reducelist(purepeak[1],nanindex)
    err=list(err)
    err[0]=reducelist(err[0],nanindex)
    err[1]=reducelist(err[1],nanindex)
    xbyb2=reducelist(xbyb2,nanindex)
    dxdiff_sum_dp0=map(lambda p:2*gxy[0]*p[1]/((1+gxy[0]*p[1]/p[0])*p[0])**2,purepeak[0])
    dxdiff_sum_dp1=map(lambda p:2*gxy[0]/p[0]/((1+gxy[0]*p[1]/p[0]))**2,purepeak[0])
    dydiff_sum_dp0=map(lambda p:2*gxy[1]*p[1]/((1+gxy[1]*p[1]/p[0])*p[0])**2,purepeak[1])
    dydiff_sum_dp1=map(lambda p:2*gxy[1]/p[0]/((1+gxy[1]*p[1]/p[0]))**2,purepeak[1])
    dxdiff_sum=map(lambda p0,p1,p2:sqrt(p0[0]*p0[0]*p1*p1+p0[1]*p0[1]*p2*p2),err[0],dxdiff_sum_dp0,dxdiff_sum_dp1)
    dydiff_sum=map(lambda p0,p1,p2:sqrt(p0[0]*p0[0]*p1*p1+p0[1]*p0[1]*p2*p2),err[1],dydiff_sum_dp0,dydiff_sum_dp1)
    dxbyb2=map(lambda p1,p2,p3,p4:sqrt(4*p1*p1*p2*p2+4*p3*p3*p4*p4),xdiff_sum,dxdiff_sum,ydiff_sum,dydiff_sum)
    #print xbyb2,dxbyb2
    dxb2x=map(lambda p1,p2:(-1./p1**2+0.5*sqrt(1./p1**4-1./p1**3)+0.5/sqrt(p1**4-p1**5))*p2,xbyb2,dxbyb2)
    dxdata=map(lambda p1,p2,p3,p4:sqrt((ar*p3*p2)**2+(ar*p1*p4)**2),xdiff_sum,dxdiff_sum,xb2x,dxb2x)
    dydata=map(lambda p1,p2,p3,p4:sqrt((ar*p3*p2)**2+(ar*p1*p4)**2),ydiff_sum,dydiff_sum,xb2x,dxb2x)
    #fit
    xfit=bpmfit(xharp,(xdata,ydata),(dxdata,dydata))
    px,pxerr=xfit.fit()
    yfit=bpmfit(yharp,(xdata,ydata),(dxdata,dydata))
    py,pyerr=yfit.fit()
    return px,py

def calibrateonebpm_ds(pos,peak,ped):
    index=getmidpoint(pos)
    ar=34.925
    purepeak=map(lambda p1,p2:map(lambda p3:map(lambda p4,p5:float(p4-p5),p3,p2),p1),peak,ped)
    gx=map(lambda p1,p2:p1[0]*(1-2/ar*p2[0])/(p1[1]*(1+2/ar*p2[0])),purepeak[0],pos)
    gy=map(lambda p1,p2:p1[0]*(1-2/ar*p2[1])/(p1[1]*(1+2/ar*p2[1])),purepeak[1],pos)
    return gx[index],gy[index]

def getar(gxy,pos,peak,ped):
    #get radius of bpm
    index=getmidpoint(pos)
    purepeak=map(lambda p1,p2:map(lambda p3:map(lambda p4,p5:float(p4-p5),p3,p2),p1),peak,ped)
    arx=map(lambda p1,p2:2*p2[0]*(p1[0]+gxy[0]*p1[1])/(p1[0]-gxy[0]*p1[1]),purepeak[0],pos)
    ary=map(lambda p1,p2:2*p2[1]*(p1[0]+gxy[1]*p1[1])/(p1[0]-gxy[1]*p1[1]),purepeak[1],pos)
    print index,arx,ary
    return (arx[index]+ary[index])/2.

def calibrateone(gxy,bpmapos,bpmbpos,peaks,ped):
    ar=34.925
    #bpm A,note bpm coordinate x is hall coordinate y
    peakx,peaky,pedx,pedy,errx,erry=[],[],[],[],[],[]
    for i in range(len(bpmapos)):
	peakx.append(peaks[i][0][2:4])
	errx.append(peaks[i][1][2:4])
	peaky.append(peaks[i][0][:2])
	erry.append(peaks[i][1][:2])
    pedx=ped[0][2:4]
    pedy=ped[0][0:2]
    #gxy[0],gxy[1]=calibrateonebpm_ds(bpmapos,(peakx,peaky),(pedx,pedy))
    #ar=getar(gxy[:2],bpmapos,(peakx,peaky),(pedx,pedy))
    (cax,cay)=calibrateonebpm(gxy[:2],ar,bpmapos,(peakx,peaky),(errx,erry),(pedx,pedy))
    #bpm B
    peakx,peaky,pedx,pedy,errx,erry=[],[],[],[],[],[]
    for i in range(len(bpmbpos)):
	peakx.append(peaks[i][0][6:])
	errx.append(peaks[i][1][6:])
	peaky.append(peaks[i][0][4:6])
	erry.append(peaks[i][1][4:6])
    pedx=ped[0][6:]
    pedy=ped[0][4:6]
    #gxy[2],gxy[3]=calibrateonebpm_ds(bpmbpos,(peakx,peaky),(pedx,pedy))
    #ar=getar(gxy[2:],bpmbpos,(peakx,peaky),(pedx,pedy))
    (cbx,cby)=calibrateonebpm(gxy[2:],ar,bpmbpos,(peakx,peaky),(errx,erry),(pedx,pedy))
    return [ar,gxy,cax,cay,cbx,cby]

#keywords used for seaching the right harp data and calibration run,including period, arm, current,...
def calibrate(keywords,rootpath=os.getenv("REPLAY_OUT_PATH"),treename="T"):
    debug=0
    #get harp data
    hd=harpinfo()
    tmp=hd.getdata(keywords)
    if not tmp:return False
    peak_04=[tmp[i]["harp04"] for i in range(tmp["ndata"])]
    peak_05=[tmp[i]["harp05"] for i in range(tmp["ndata"])]
    run=[tmp[i]["run"][0] for i in range(tmp["ndata"])]
    pedrun=tmp["pedrun"][0]
    #currrun=tmp["currrun"][0]
    print "calibrate bpm with run",run,"and pedestal run %i,"%pedrun,"keywords: ",keywords
    #get survey data
    period=runinfo()
    harp04survey=period.harp04survey(run[0])
    harp05survey=period.harp05survey(run[0])
    bpmasurvey=period.bpmasurvey(run[0])
    bpmbsurvey=period.bpmbsurvey(run[0])
    harpsurvey=[harp04survey,harp05survey]
    #get position in bpm for harp data
    harp=harppos(harpsurvey)
    bpma=bpmpos(bpmasurvey)
    bpmb=bpmpos(bpmbsurvey)
    bpmapos_bpm,bpmbpos_bpm,pos_harp04,pos_harp05=[],[],[],[]
    for i in range(len(peak_04)):
	pos_harp04.append(harp.getpos_04(peak_04[i]))
	pos_harp05.append(harp.getpos_05(peak_05[i]))
	bpmapos_bpm.append(bpma.getpos_bpm(pos_harp04[i],pos_harp05[i]))
	bpmbpos_bpm.append(bpmb.getpos_bpm(pos_harp04[i],pos_harp05[i]))
    #get peak data
    peaks=[]
    for i in range(len(run)):
	tmp=bpm_peak(run[i],0,rootpath,treename)
	if not tmp:return False
	peaks.append(tmp)
    #pedpeaks=bpm_peak(pedrun,1,rootpath,treename)
    pedpeaks=[[2125.10698642,3231.40291772,8823.49384253,5821.47710681,1381.83193059,6919.63705969,10105.2690181,2495.87669436 ],[0]*8]
    gxy=[0.8871,0.9949,0.9267,1.1283]
    #gxy=[1,1,1,1]
    
    if not pedpeaks:return False
    #get calibration data
    #calibrate gx,gy for 2 bpms
    #gxyfit=gxgyfit(currrun,rootpath,bpmapos_bpm,bpmbpos_bpm,peaks,pedpeaks)
    #gxy=gxyfit.fit()
    #calibrate a,b,c
    out=calibrateone(gxy,bpmapos_bpm,bpmbpos_bpm,peaks,pedpeaks)
    if debug:print out
    else:bpmconstsave(run,out,pedpeaks,keywords)
    return out

def bpmconstsave(run,const,pedpeaks,keywords):
    dbdir=os.getenv("BEAMDBPATH")
    if dbdir==None:
	print "please define BEAMDBPATH in your env"
	return False
    pydb=os.path.join(dbdir,"pyDB")
    if not os.path.exists(pydb):os.makedirs(pydb)
    run=sorted(run)
    #save const
    filename=os.path.join(pydb,"bpm_%i.dat"%(run[0]))
    #print "please input available run period for calibration run %i ..."%sorted(run)[0]
    #runperiod=raw_input("run period(for example 5490-5495,5496 5700,5754):")
    runperiod=""
    for r in run:runperiod+="%i,"%r
    runperiod=runperiod[:-1]
    #curravail=raw_input("avail for curr(nA):")
    #period=runinfo()
    #currepics=period.current(run[0])
    #def curr100(curr):
	#return 100 if curr>100 else curr
    #curravail=str(int((100 if currepics>100 else currepics)/25)*25)

    datfile=open(filename,"w")
    datfile.write("All of the survey info directly come from survey data,please read survey report to get the detail info about the coordinate\n")
    datfile.write("Please contact pengjia immediately if you have any questions(email,gtalk,phone...)\n")
    datfile.write("keywords: ")
    for keyword in keywords:
	datfile.write("%s "%keyword)
	if "nA" in keyword:curravail=keyword[:-2]
    datfile.write("\n\n")
    datfile.write("------------------------------------------------\n\n")
    datfile.write("avail run period:%s\n"%runperiod)
    datfile.write("avail curr(nA):%s\n"%curravail)
    #datfile.write("avail curr(nA):100 50 25\n")
    datfile.write("target z position(mm,support multi):-14.135 0 14.135 -10.81 -13.6271\n")
    datfile.write("pedestal peak:%f %f %f %f %f %f %f %f\n"%tuple(pedpeaks[0]))
    datfile.write("bpma ar,gx,gy:%.15f %.15f %.15f\n"%(const[0],const[1][0],const[1][1]))
    datfile.write("bpmb ar,gx,gy:%.15f %.15f %.15f\n"%(const[0],const[1][2],const[1][3]))
    for i in range(4):
	ab="a" if i<2 else "b"
	xy="x" if i%2==0 else "y"
	datfile.write("bpm%s %s a,b,c:"%(ab,xy))
	for j in range(len(const[i+2])):
	    datfile.write("%.15f "%const[i+2][j])
	datfile.write("\n")

def sortrun(frun,runs):
    diffrun=map(lambda r:(abs(frun-r),r),runs)
    runs=sorted(diffrun,key=lambda r:r[0])
    return [run[1] for run in runs if run[0]<15000]

def bpmconstread(run,runcurr=50):
    if os.getenv("BEAMUPDATECHECK")=="1":updatepackage()
    dbdir=os.getenv("BEAMDBPATH")
    if dbdir==None:
	print "please define BEAMDBPATH in your env"
	return False
    pydb=os.path.join(dbdir,"pyDB")
    if not os.path.exists(pydb):os.makedirs(pydb)
    datfiles=glob.glob(os.path.join(pydb,"bpm_*.dat"))
    calruns=[]
    for datfile in datfiles:
	calruns.append(int(re.split("[_.]",os.path.split(datfile)[1])[1]))
    calruns=sortrun(run,calruns)
    def constread(filename,dbread=False):
	datfile=open(filename,"r")
	if not isinstance(dbread,dict):
	    dbread={"period":False,"curr":False,"tgtz":False,"ped":False,"bpmaconst":False,"bpmbconst":False}
	importfile=False
	for line in datfile:
	    reline=re.split("[:\n]",line)
	    data=re.split(" ",reline[1].strip())
	    if "run period" in reline[0] and not dbread["period"]:
		#runperiod=[int(x) for x in data]
		runperiods=re.split(",",reline[1].strip())
		runperiod=[]
		for p in runperiods:
		    if p.isdigit():runperiod.append(int(p))
		    else:
			ptmp=[int(x) for x in re.split("[ ~-]",p)]
			runperiod+=range(ptmp[0],ptmp[1]+1)
		dbread["period"]=runperiod
	    #elif "run orbit" in reline[0] and not dbread["orbit"]:
		#dbread["orbit"]=int(data[0])
	    elif "curr" in reline[0] and not dbread["curr"]:
		dbread["curr"]=[float(x) for x in data]
	    elif "target z" in reline[0] and not dbread["tgtz"]:
		dbread["tgtz"]=[float(x) for x in data]
	    elif "pedestal" in reline[0] and not dbread["ped"]:
		dbread["ped"]=[float(x) for x in data]
	    elif "bpma" in reline[0]:
		#if not dbread["bpmasurvey"]:
		 #   if all(x in reline[0] for x in ["survey","position"]):a_survey_pos=[float(x) for x in data]
		 #   elif all(x in reline[0] for x in ["survey","angle"]):a_survey_ang=[float(x) for x in data]
		if not dbread["bpmaconst"]:
		    if all(x in reline[0] for x in ["ar","gx","gy"]):a_rgxy=[float(x) for x in data]
		    elif all(x in reline[0] for x in ["x","a","b","c"]):ax=[float(x) for x in data]
		    elif all(x in reline[0] for x in ["y","a","b","c"]):ay=[float(x) for x in data]
	    elif "bpmb" in reline[0]:
		#if not dbread["bpmbsurvey"]:
		 #   if all(x in reline[0] for x in ["survey","position"]):b_survey_pos=[float(x) for x in data]
		 #   elif all(x in reline[0] for x in ["survey","angle"]):b_survey_ang=[float(x) for x in data]
		if not dbread["bpmbconst"]:
		    if all(x in reline[0] for x in ["ar","gx","gy"]):b_rgxy=[float(x) for x in data]
		    elif all(x in reline[0] for x in ["x","a","b","c"]):bx=[float(x) for x in data]
		    elif all(x in reline[0] for x in ["y","a","b","c"]):by=[float(x) for x in data]
	    if "import" in reline[0]:
		importfile=os.path.join(os.path.split(filename)[0],re.split(" ",reline[0].strip())[1])
		break
	#if locals().has_key("a_survey_pos"):dbread["bpmasurvey"]=[a_survey_pos,a_survey_ang]
	#if locals().has_key("b_survey_pos"):dbread["bpmbsurvey"]=[b_survey_pos,b_survey_ang]
	if locals().has_key("ay"):dbread["bpmaconst"]=[tuple(a_rgxy),ax,ay,len(ax)]
	if locals().has_key("by"):dbread["bpmbconst"]=[tuple(b_rgxy),bx,by,len(bx)]
	if importfile:return constread(importfile,dbread)
	else:return dbread
#	return (runperiod,curr,ped,\
#	    ([a_survey_pos,a_survey_ang],[b_survey_pos,b_survey_ang],(tgtz,runorbits)),\
#	    [tuple(a_rgxy),ax,ay,tuple(b_rgxy),bx,by])
    for calrun in calruns:
	const=constread(os.path.join(pydb,"bpm_%i.dat"%calrun))
	if run in const["period"] and any(abs(x-(100 if runcurr>100 else runcurr))<25 for x in const["curr"]):
	    del const["period"]
	    del const["curr"]
	    return const
    print "sorry no bpm calibration constant available for run %i,please contact pengjia"%run
    return False

def checkcalibrunavail(keywords,rootpath=os.getenv("REPLAY_OUT_PATH"),treename="T"):
    hd=harpinfo()
    tmp=hd.getdata(keywords)
    if not tmp:return False
    runs=[tmp[i]["run"][0] for i in range(tmp["ndata"])]
    notreplayed=[]
    for run in runs:
	path=os.path.join(os.getenv("REPLAY_OUT_PATH"),"bpm_%i.root"%run)
	bakpath=os.path.join(os.getenv("REPLAY_OUT_PATH"),"bpm_%i.root_bak"%run)
	if not os.path.exists(path):
	    if os.path.exists(bakpath):shutil.copy(bakpath,path)
	    else:notreplayed.append(run)
	else:
	    if os.path.getsize(path)<1e6:notreplayed.append(run)
    if len(notreplayed)<1:return True
    else:
	for r in notreplayed:print r,
	print "need to replay first"
	return False