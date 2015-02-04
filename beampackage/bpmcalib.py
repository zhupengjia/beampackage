#!/usr/bin/env python
import os,re,numpy
from harppos import *
from bpmfit import *
from signalfilter import decode,getrealpos
from runinfo import *

#class to calibrate bpm
class bpmcalib:
    def __init__(self,keywords=False,treename="T",rootpath=os.getenv("REPLAY_OUT_PATH"), onlyb=False,forcefastbus=False,forceredecode=False,ab=False):
	self.keywords=keywords
	self.treename=treename
	self.rootpath=rootpath
	self.pklpath=getpklpath(rootpath)
	self.onlyb=onlyb
	self.forcefastbus=forcefastbus
	self.forceredecode=forceredecode
	self.period=runinfo()
	if not ab:
	    self.ab=[1] if self.onlyb else [0,1]
	else:self.ab=ab
	self.chan=[2,3,0,1]
	self.bpmraw,self.constfile={},{}
	
	if self.keywords:
	    self.filtertype="s"
	    self.getharpinfo()
	    self.getsurvey()
	    self.getcalibconf()
	    self.gethardpos()
	
    #get harp information
    def getharpinfo(self):
	harpdata=harpinfo()
	tmp=harpdata.getdata(self.keywords)
	if not tmp:raise Exception("no harp data found")
	if not self.onlyb:self.peak_04=[tmp[i]["harp04"] for i in range(tmp["ndata"])]
	self.peak_05=[tmp[i]["harp05"] for i in range(tmp["ndata"])]
	self.run=[tmp[i]["run"][0] for i in range(tmp["ndata"])]
	self.pedrun=tmp["pedrun"][0]
	try:self.availruns={"a":tmp["availa"],"b":tmp["availb"]}
	except:self.availruns=False
	print "calibrate bpm with run",self.run,"and pedestal run %i,"%self.pedrun,"keywords: ",self.keywords
	
    #get survey information
    def getsurvey(self,run=False):
	if not run:run=self.run[0]
	#get position in bpm for harp data
	self.harp=harppos(run)
	self.bpm={"a":bpmpos(run,"a"),"b":bpmpos(run,"b")}
	self.currepics=self.period.current(run)
	
    #get calibration configure
    def getcalibconf(self):
	self.calibconf=calibsetting(self.keywords)
	self.datanum=self.calibconf.datanum
	
    #get position at bpm from harp data
    def gethardpos(self):
	pos1,pos2,posbpma,posbpmb=[],[],[],[]
	for i in range(len(self.peak_05)):
	    if self.onlyb:
		bpmaposhall=self.bpm["a"].posbpm2hall(self.calposfromraw(self.run[i],0))
		pos1.append([numpy.mean(x) for x in bpmaposhall])
	    else:pos1.append(self.harp.getpos_04(self.peak_04[i]))
	    pos2.append(self.harp.getpos_05(self.peak_05[i]))
	    if not self.onlyb:
		posbpma.append(self.bpm["a"].getpos_bpm(pos1[i],pos2[i]))
		posbpmb.append(self.bpm["b"].getpos_bpm(pos1[i],pos2[i]))
	    else:
		#print self.bpm["b"].
		posbpmb.append(self.harp.getpos_05_local(self.peak_05[i]))
	hardpos=[posbpma,posbpmb]
	r=map(lambda p:p[0]**2+p[1]**2,posbpmb)
	self.centerid=r.index(min(r))
	self.hardpos=[]
	for i in self.ab:
	    self.hardpos.append(hardpos[i])
	#print out
	print "hard position is:"
	for i in range(len(self.hardpos[0])):
	    print self.run[i],
	    for j in range(len(self.hardpos)):
		for p in self.hardpos[j][i]:
		    print "%1.2f,"%p,
		print "\t",
	    print
	    
    #get bpm calibration constant, used for calibrating bpm B only with A and harp info
    def getcalibconst(self,ab,run=False):
	if not run:run=self.run[0]
	tmp=self.period.bpmconstread(run,self.forcefastbus)[ab]
	if not tmp:
	    print "can not find const for run %i"%run
	    return False
	#pedestal for run, read from pedestal.pkl from database
	pedtmp=False
	if not self.period.ifautogain(run):
	    pedtmp=self.period.pedestal(run,self.forcefastbus)[ab]	
	if not pedtmp:pedtmp=tmp["ped"]
	pedestal=map(lambda a,b:a+b,pedtmp,tmp["offset"])
	calconst=tmp["const"]
	fitorder=tmp["fitorder"]
	self.constfile[ab]=tmp["constfile"]
	return pedestal,calconst,fitorder

    #calculate position from raw data,ab=0 for bpm a ,1 for b
    def calposfromraw(self,run,ab,rotate=False):
	self.getrawdata(run)
	ab="a" if ab==0 else "b"
	tmp=self.getcalibconst(ab)
	if not tmp:return False,False
	ped,const,fitorder=tmp
	raw=self.bpmraw[run][:4] if ab=="a" else self.bpmraw[run][4:]
	raw=[raw[c]-ped[c] for c in range(4)]
	x,y=getrealpos(raw,const,fitorder)
	x,y=x[x>-100],y[y>-100]
	x,y=x[x<100],y[y<100]
	if rotate:x,y,z=self.bpm[ab].posbpmrotate([x,y])
	return x,y
	
    #get raw data
    def getrawdata(self,run,ped=False,eventcut=False):
	if run in self.bpmraw.keys():return
	bpmrawpkl=self.pklpath.getpath("raw","%sbpm"%self.filtertype,run)
	currpkl=self.pklpath.getpath("raw","curr",run)
	availpkl=self.pklpath.getpath("raw","bpmavail",run)
	if not os.path.exists(bpmrawpkl):
	    runpath=checkrunavail(self.rootpath,run)
	    if not runpath:raise Exception("no file found for run %i"%run)
	    d=decode(runpath,self.treename,forcefastbus=self.forcefastbus,forceredecode=self.forceredecode)
	    d.autodecode()
	raw=zload(bpmrawpkl)
	#ped or signal cut
	if ped:
	    curr=zload(currpkl)
	    #get average curr
	    nocurr=0.01 #below this current will deal as no signal
	    currshift=500
	    curr=curr<nocurr
	    curr1=numpy.concatenate((numpy.zeros(currshift),curr[:-currshift]))
	    bpmavail=curr*curr1
	else:
	    bpmavail=zload(availpkl)
	#event cut
	if not eventcut:
	    ecut=getbpmeventcut()
	    eventcut=ecut.getcut(run,self.forcefastbus)
	#filter the unwanted event
	if eventcut:
	    if (len(bpmavail)-eventcut[1])>=0:
		cut=numpy.asarray([0]*eventcut[0]+[1]*(eventcut[1]-eventcut[0])\
		    +[0]*(len(bpmavail)-eventcut[1]),dtype=numpy.int32)
	    else:
		cut=numpy.asarray([0]*eventcut[0]+[1]*(len(bpmavail)-eventcut[0])\
		    ,dtype=numpy.int32)
	    raw=[x+bpmavail*1e6+cut*1e6-2e6 for x in raw]
	else:raw=[x+bpmavail*1e6-1e6 for x in raw]
	raw=[x[x>-1e4] for x in raw]
	self.bpmraw[run]=raw
	
    #get center point
    def getcenterpoint(self,pos):
	r=map(lambda p:p[0]**2+p[1]**2,pos)
	return r.index(min(r))
	
    #get bpm raw beak and calibration configure
    def bpmpeak(self):
	for r in self.run:self.getrawdata(r)
	#ped peaks and offset
	if self.calibconf.pedpeaks:
	    self.pedpeaks=self.calibconf.pedpeaks
	else:
	    pedtmp=self.period.pedestal(self.run[0],self.forcefastbus)
	    if pedtmp["a"] and pedtmp["b"]:
		self.pedpeaks=pedtmp["a"]+pedtmp["b"]
	    else:
		self.getrawdata(self.pedrun,True)
		self.pedpeaks=[numpy.mean(r) for r in self.bpmraw[self.pedrun]]
	if self.calibconf.offset:self.offset=self.calibconf.offset
	else:self.offset=[0]*8
	self.peaks=map(lambda r:[numpy.asarray([numpy.mean(x)]*self.datanum,dtype=numpy.float32) for x in self.bpmraw[r]],self.run)
	if self.calibconf.gxy:self.gxy=self.calibconf.gxy
	else:self.gxy=[False,False,False,False]
	
    #calibrate gx and gy
    def calibrategxgy(self,pos,peak,ped,offset):
	index=self.getcenterpoint(pos)
	ar=self.ar
	purepeak=map(lambda p1,p2,p0:map(lambda p3:map(lambda p4,p5,p6:p4-p5-p6,p3,p2,p0),p1),peak,ped,offset)
	gx=map(lambda p1,p2:p1[0]*(1-2/ar*p2[0])/(p1[1]*(1+2/ar*p2[0])),purepeak[0],pos)
	gy=map(lambda p1,p2:p1[0]*(1-2/ar*p2[1])/(p1[1]*(1+2/ar*p2[1])),purepeak[1],pos)
	return gx[index],gy[index]
	
    #calibrate one bpm
    def calibrateone(self,gxy,pos,peak,ped,offset):
	#purepeak:1st level:x,y;2nd level:n pos;3rd level:x+,x-
	#pos:1st level:n pos;2nd level:x,y
	purepeak=map(lambda p1,p2,p0:map(lambda p3:map(lambda p4,p5,p6:p4-p5-p6,p3,p2,p0),p1),peak,ped,offset)
	xdiff_sum=map(lambda p:(p[0]-gxy[0]*p[1])/(p[0]+gxy[0]*p[1]),purepeak[0])
	ydiff_sum=map(lambda p:(p[0]-gxy[1]*p[1])/(p[0]+gxy[1]*p[1]),purepeak[1])
	xbyb2=map(lambda p1,p2:p1**2+p2**2,xdiff_sum,ydiff_sum)
	xb2x=map(lambda p:1/p-1/numpy.sqrt(p)*numpy.sqrt(1/p-1),xbyb2)
	xdata=map(lambda p1,p2:self.ar*p1*p2,xdiff_sum,xb2x)
	ydata=map(lambda p1,p2:self.ar*p1*p2,ydiff_sum,xb2x)
	xharp=map(lambda p:p[0],pos)
	yharp=map(lambda p:p[1],pos)
	#filternan
	nanxdata=[all(x) for x in numpy.isnan(xdata)]
	nanydata=[all(x) for x in numpy.isnan(ydata)]
	nanxy=nanxdata or nanydata
	for i in range(len(nanxy)-1,-1,-1):
	    if nanxy[i]:
		del xdata[i],ydata[i],xharp[i],yharp[i]	
	#fit
	centerid=self.getcenterpoint(pos)
	xfit=bpmfit(self.keywords,0,xharp,(xdata,ydata),centerid)
	px,pxerr,pxval=xfit.fit()
	yfit=bpmfit(self.keywords,1,yharp,(ydata,xdata),centerid)
	py,pyerr,pyval=yfit.fit()
	
	return px,py,pxerr,pyerr,pxval,pyval	
	
    #calibrate
    def calibrate(self,rawconst=False):
	self.ar=34.925
	self.cx,self.cy,self.ex,self.ey,self.px,self.py=[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]
	self.bpmpeak()
	#read const input
	if rawconst:
	    if "offset" in rawconst.keys():
		for i in range(len(self.ab)):
		    for j in range(4):
			self.pedpeaks[self.ab[i]*4+j]=rawconst["ped"][i*4+j]
			self.offset[self.ab[i]*4+j]=rawconst["offset"][i*4+j]
	    if "gxgy" in rawconst.keys():
		for i in range(len(self.ab)):
		    for j in range(2):
			self.gxy[self.ab[i]*2+j]=rawconst["gxgy"][i*2+j]
	#calibrate	
	for i in range(len(self.ab)):
	    peak=([],[])
	    xchan=self.chan[0]+4*self.ab[i]
	    ychan=self.chan[2]+4*self.ab[i]
	    for j in range(len(self.hardpos[i])):
		peak[0].append(self.peaks[j][xchan:xchan+2])
		peak[1].append(self.peaks[j][ychan:ychan+2])
	    ped=(self.pedpeaks[xchan:xchan+2],self.pedpeaks[ychan:ychan+2])
	    print "-----------------------------------------------------------------------------------------"
	    offset=(self.offset[xchan:xchan+2],self.offset[ychan:ychan+2])
	    #get gxy
	    gxychan=self.ab[i]*2
	    if not self.gxy[gxychan]:
		self.gxy[gxychan],self.gxy[gxychan+1]=\
		    self.calibrategxgy(self.hardpos[i],peak,ped,offset)
	    #calibrate a,b,c
	    self.cx[self.ab[i]],self.cy[self.ab[i]],self.ex[self.ab[i]],self.ey[self.ab[i]],self.px[self.ab[i]],self.py[self.ab[i]]=self.calibrateone(self.gxy[gxychan:gxychan+2],self.hardpos[i],peak,ped,offset)
	    
    #save const a or b,used for constsave function
    def __constsaveone(self,ab):
	dbdir=os.getenv("BEAMDBPATH")
	if dbdir==None:
	    print "please define BEAMDBPATH in your env"
	    return False
	pydb=os.path.join(dbdir,"pyDB")
	if not os.path.exists(pydb):os.makedirs(pydb)
	run=sorted(self.run)
	if not self.period.ifhapavail(run[0]) or self.forcefastbus:fastbus=True
	else:fastbus=False
	#save const
	if fastbus:
	    filename=os.path.join(pydb,"bpm%sfb_%i.dat"%(ab,run[0]))
	else:
	    filename=os.path.join(pydb,"bpm%s_%i.dat"%(ab,run[0]))
	if self.availruns:runperiod=self.availruns[ab]
	else:
	    runperiod=""
	    for r in run:runperiod+="%i,"%r
	    runperiod=runperiod[:-1]
	fitorder=self.calibconf.fitorder
	datfile=open(filename,"w")
	datfile.write("All of the survey info directly come from survey data,please read survey report to get the detail info about the coordinate\n")
	datfile.write("Please contact pengjia immediately if you have any questions(email,gtalk,phone...)\n")
	datfile.write("keywords: ")
	for keyword in self.keywords:
	    datfile.write("%s "%keyword)
	    if "nA" in keyword:curravail=keyword[:-2]
	datfile.write("\n\n")
	datfile.write("------------------------------------------------\n\n")
	datfile.write("avail run period:%s\n"%runperiod)
	try:datfile.write("avail curr(nA):%i\n"%(int(self.currepics)))
	except:datfile.write("avail curr(nA):%s\n"%(curravail))
	datfile.write("target z position(mm,support multi):-14.135 0 14.135 -10.81 -13.6271\n")
	ped=self.pedpeaks[:4] if ab=="a" else self.pedpeaks[4:]
	offset=self.offset[:4] if ab=="a" else self.offset[4:]
	datfile.write("pedestal peak:%f %f %f %f\n"%tuple(ped))
	datfile.write("offset:%f %f %f %f\n"%tuple(offset))
	abnum=0 if ab=="a" else 1
	datfile.write("bpm%s ar,gx,gy:%.15f %.15f %.15f\n"%(ab,self.ar,self.gxy[abnum*2],self.gxy[abnum*2+1]))
	datfile.write("fitorder:%i %i\n"%(fitorder[0],fitorder[1]))
	cxy=[self.cx[abnum],self.cy[abnum]]
	exy=[self.ex[abnum],self.ey[abnum]]
	for i in range(2):
	    xy="x" if i==0 else "y"
	    datfile.write("bpm%s %s a,b,c:"%(ab,xy))
	    for j in range(len(cxy[i])):
		datfile.write("%.15f "%cxy[i][j])
	    datfile.write("\n")
	#for i in range(2):
	    #xy="x" if i==0 else "y"
	    #datfile.write("bpm%s %s para error:"%(ab,xy))
	    #for j in range(len(exy[i])):
		#datfile.write("%.15f "%exy[i][j])
	    #datfile.write("\n")
	datfile.write("fval:%.7f %.7f"%(self.px[abnum],self.py[abnum]))
	datfile.write("\n")
	datfile.close()
	#print constant
	print "\n\n-----------------------------------------"
	for line in open(filename,"r"):print line.strip()
	print "-----------------------------------------\n\n"
	    
    #save constant
    def constsave(self):
	dbdir=os.getenv("BEAMDBPATH")
	if not self.onlyb:self.__constsaveone("a")
	self.__constsaveone("b")
	
    #check calibration constant
    def calibcheck(self):
	try:
	    from pylab import savefig,figure
	    from matplotlib.colors import LogNorm
	    from matplotlib import pyplot
	    from matplotlib.ticker import MultipleLocator, FormatStrFormatter  
	except:
	    print "sorry the matplotlib package is needed for plotting!"
	    return 
	fig =figure(figsize=(5.0*len(self.ab), 5.0), dpi=100)
	axes=[]
	majorLocator= MultipleLocator(1)
	minorLocator= MultipleLocator(0.2)	
	for i in range(len(self.ab)):
	    axes.append(fig.add_subplot(1,len(self.ab),i+1))
	    axes[i].clear()
	    xall,yall=numpy.zeros(0,dtype=numpy.float32),numpy.zeros(0,dtype=numpy.float32)
	    for r in self.run:
		x,y=self.calposfromraw(r,self.ab[i])
		xall=numpy.concatenate((xall,x))
		yall=numpy.concatenate((yall,y))
	    xymax=max([abs(min([xall.min(),yall.min()])),abs(max([xall.max(),yall.max()]))])*1.2
	    histrange=[[-xymax,xymax],[-xymax,xymax]]
	    axes[i].hist2d(xall,yall,bins=300,range=histrange,norm=LogNorm())
	    #harp pos
	    hardpos=[[x[0] for x in self.hardpos[i]],[x[1] for x in self.hardpos[i]]]
	    axes[i].plot(hardpos[0],hardpos[1],"+",markersize=50.,fillstyle="none")
	    axes[i].xaxis.set_major_locator(majorLocator)
	    axes[i].yaxis.set_major_locator(majorLocator)
	    axes[i].xaxis.set_minor_locator(minorLocator)
	    axes[i].yaxis.set_minor_locator(minorLocator)
	try:
	    fig.suptitle("%inA,using %s"%(self.curr,self.constfile["b"]))
	    construn=re.split("[_.]",self.constfile["b"])[1]
	    savefig("pic/points%i_%inA_%s.png"%(sorted(self.run)[0],self.curr,construn))
	except:
	    savefig("pic/points%i.png"%(sorted(self.run)[0]))
	
    def ovalfun(self,x,a,b,c):
	#par:a center,b radius,c entries radius
	return c*numpy.sqrt(1-(x-a)*(x-a)/(b*b))
	
    #check calibration constant by using slow raster
    def calibcheck_raster(self,run):
	try:
	    from pylab import savefig,figure
	    from matplotlib.colors import LogNorm
	    from matplotlib import pyplot,mlab
	    from matplotlib.ticker import MultipleLocator, FormatStrFormatter  
	    from matplotlib.patches import Ellipse
	    from scipy.optimize import curve_fit
	except:
	    print "sorry the matplotlib package is needed for plotting!"
	    return 
	#if not self.keywords:
	self.run=[run]
	self.filtertype="f"
	self.getsurvey(run)
	fig =figure(figsize=(5.0*len(self.ab), 5.0), dpi=100)
	axes=[]
	for i in range(2):
	    tmp=self.calposfromraw(run,i,rotate=True)
	    if not tmp:continue
	    x,y=tmp
	    axes.append(fig.add_subplot(121+i))
	    axes[i].clear()
	    center=[numpy.mean(x),numpy.mean(y)]
	    xyrange=[x.max()-x.min(),y.max()-y.min()]
	    xymin=min([x.min(),y.min()])
	    xymax=max([x.max(),y.max()])
	    histrange=[[xymin,xymax],[xymin,xymax]]
	    axislable=int(xymax/10.)*2
	    if axislable<1:axislable=1 
	    majorLocator= MultipleLocator(axislable)
	    minorLocator= MultipleLocator(axislable/5.)
	    axes[i].hist2d(x,y,bins=300,range=histrange,norm=LogNorm())
	    axes[i].xaxis.set_major_locator(majorLocator)
	    axes[i].yaxis.set_major_locator(majorLocator)
	    axes[i].xaxis.set_minor_locator(minorLocator)
	    axes[i].yaxis.set_minor_locator(minorLocator)
	try:    
	    fig.suptitle("%inA,using %s"%(self.curr,self.constfile["b"]))
	    construn=re.split("[_.]",self.constfile["b"])[1]
	    savefig("pic/calibcheck%i_%inA_%s.png"%(run,self.curr,construn))
	except:
	    savefig("pic/calibcheck%i.png"%(run))