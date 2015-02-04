#!/usr/bin/env python
import re,os,time,glob,sys,numpy
try:import ROOT
except:print "Error!! pyroot didn't compiled! please recompile your root!"
from array import array
from harppos import bpmpos
from runinfo import runinfo,getpklpath,sortrun,zload,zdump
from bpminsert import tgtpos_straight,tgtpos_field,bpminsertparaprep
from signalfilter import decode

class rasterrecon:
    def __init__(self):
	ROOT.gROOT.SetBatch(True)
	self.fastclkrate=103910
	#self.fastclkrate=103858
	self.shape=ROOT.TF1("pol1","pol1")
	self.rastertype=""
	self.rasterchan=[0,1]#chan in pkl file
	
    def runinit(self,runpath=False,forcefastbus=False):
	rootpath,runfilename=os.path.split(runpath)
	self.run=int(re.split("[_]",re.split("[.]",runfilename)[-2])[-1])
	self.pklpath=getpklpath(rootpath)
	self.clkpkl=self.pklpath.getpath("raw","clock",self.run)
	if not os.path.exists(self.clkpkl):
	    d=decode(runpath,forcefastbus=forcefastbus)
	    d.autodecode()
	if not os.path.exists(self.clkpkl):
	    print self.clkpkl
	    print "sorry no clock info in rootfile,please replay again"
	    return False
	self.clkraw=array("d",zload(self.clkpkl))
	self.events=len(self.clkraw)
	return True
	
    #general get size function
    def gengetsize(self,multifits,chan,eventrange=[0,1000000000]):
	c1=ROOT.TCanvas("c1","%s raster size"%self.rastertype,800,600)
	if chan<4:
	    h1=ROOT.TH1F("h1","raster size",100,min(self.rasterraw[chan]),max(self.rasterraw[chan]))
	    for i in range(eventrange[0],eventrange[1]):
		h1.Fill(self.rasterraw[chan][i])
	else:
	    ab=int((chan-4)%2)
	    xy=int((chan-4)/2)
	    pos=self.bpmpos[ab][xy][eventrange[0]:eventrange[1]]+(self.bpmavail[eventrange[0]:eventrange[1]]*1e8-1e8)
	    pos=pos[pos>-100]
	    posmean=numpy.mean(pos)
	    posrms=numpy.sqrt(numpy.mean(pow(pos-posmean,2)))
	    h1=ROOT.TH1F("h1","raster size at bpm",300,posmean-3*posrms,posmean+3*posrms)
	    for p in pos:h1.Fill(p)
	#par:c center,r radius,a entries amplitude
	def shapefit(h1,c0,c1,r0,r1,a0,a1,e1,e2):
	    self.shape.SetParLimits(0,c0,c1)
	    self.shape.SetParLimits(1,r0,r1)
	    self.shape.SetParLimits(2,a0,a1)
	    h1.Fit("shape","Qsame","",e1,e2)
	    c=self.shape.GetParameter(0)
	    r=self.shape.GetParameter(1)
	    a=self.shape.GetParameter(2)
	    return c,r,a
  
	c=h1.GetMean()
	hmax=h1.GetMaximum()
	a=hmax
	xmax=h1.GetXaxis().GetXmax()
	xmin=h1.GetXaxis().GetXmin()
	r=(xmax-xmin)/2.
	#multitimes fit:
	c,r,a=multifits(h1,shapefit,c,r,a,xmin,xmax,hmax)
	#if not os.path.exists("pic"):os.makedirs("pic")
	if chan<4:
	    c1.Print("pic/%srastersize_%i_%i_%i_%i.png"%(self.rastertype,chan,self.run,eventrange[0],eventrange[1]),"png")
	else:
	    ab="a" if ab==0 else "b"
	    xy="x" if xy==0 else "y"
	    c1.Print("pic/%sbpm%s%ssize_%i_%i_%i.png"%(self.rastertype,ab,xy,self.run,eventrange[0],eventrange[1]),"png")
	#c1.WaitPrimitive()
	return {"radius":r,"center":c}
    
    def getsize(self,eventrange=[0,1000000000]):
	self.rasterpkl=self.pklpath.getpath("raw","raster",self.run)
	if not os.path.exists(self.rasterpkl):
	    print "sorry no raster info in rootfile,please replay again"
	    return False
	self.rasterraw=[array("d",a) for a in zload(self.rasterpkl)]
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	self.size=[]
	for i in range(2):
	    self.size.append(self.gengetsize(self.rasterchan[i],eventrange))
	return self.size

    def getsizeinbpm(self,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	bpms=["a","b"]
	self.bpmpkl,self.bpmpos=[],[]
	for ab in bpms:
	    self.bpmpkl.append(self.pklpath.getpath("pos","fbpm%srot"%ab,self.run))
	    self.bpmpos.append(zload(self.bpmpkl[-1]))
	    self.bpmavailpkl=self.pklpath.getpath("raw","fbpmavail",self.run)
	    self.bpmavail=zload(self.bpmavailpkl)
	self.sizeinbpm=[]
	for i in range(2):#x,y
	    self.sizeinbpm.append([])
	    for j in range(2):#bpma,b
		self.sizeinbpm[i].append(self.gengetsize(4+i*2+j,eventrange))
	return self.sizeinbpm
    
class slrecon(rasterrecon):
    def __init__(self):
	rasterrecon.__init__(self)
	self.amperiod=1/60.
	self.period=1/99.412
	self.amphase=[0,0]
	self.phase=[0,0]
	self.fitclockrate=[0,0] #fast clock rate calculated from fit
	#par:0 center,1 radius,2 entries radius
	self.shape=ROOT.TF1("shape","[2]*sqrt(1-1/([1]*[1])*(x-[0])*(x-[0]))") #oval
	self.rastertype="slow"
	self.rasterchan=[2,3]
	
    #general get size function
    def gengetsize(self,chan,eventrange=[0,1000000000]):
	#multitimes fit:
	def multifits(hist,fitfun,c,r,a,xmin,xmax,hmax):
	    #three times fit
	    c,r,a=fitfun(hist,c-r*0.2,c+r*0.2, r*0.5,r, a*0.1,a, xmin,xmax)
	    c,r,a=fitfun(hist,c-r*0.2,c+r*0.2, r-r*0.2,r+r*0.2, a*0.5,hmax, c-r,c+r)
	    c,r,a=fitfun(hist,c-r*0.2,c+r*0.2, r-r*0.2,r+r*0.2, a*0.5,hmax, c-r,c+r)
	    return c,r,a
	return rasterrecon.gengetsize(self,multifits,chan,eventrange)

    def __slrasterfunc(self,x,par):
	#par:am phase,phase,am period,period,max am,center am,fastclkrate
	if not(par[2] and par[3] and par[6]):return -10000000
	time=x[0]/par[6]
	a=par[4]/numpy.sqrt(par[2]) #normalize
	#time=x[0]
	dp=(time+par[0])%(par[2]*4)
	if dp<=par[2]:
	    am=numpy.sqrt(dp)
	elif dp>par[2] and dp<=par[2]*2:
	    
	    am=numpy.sqrt(abs(2*par[2]-dp))
	elif dp>par[2]*2 and dp<=par[2]*3:
	    am=-numpy.sqrt(abs(-2*par[2]+dp))
	else:
	    am=-numpy.sqrt(abs(4*par[2]-dp))
	return a*am*numpy.sin(2*numpy.pi*time/par[3]+par[1])+par[5]
    	
    def getphase(self,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	clkrange=[self.clkraw[eventrange[0]],self.clkraw[eventrange[1]-1]]
	clkoffset=self.clkraw[eventrange[0]]
	#clkoffset=0
	clkrawoff=array("d",[a-clkoffset for a in self.clkraw[eventrange[0]:eventrange[1]]])
	wave=ROOT.TF1("wave",self.__slrasterfunc,clkrange[0],clkrange[1],7)
	for i in range(0,2):
	    h1=ROOT.TGraph(eventrange[1]-eventrange[0],clkrawoff,self.rasterraw[self.rasterchan[i]][eventrange[0]:eventrange[1]])
	    wave.SetParLimits(0,-2*self.amperiod,2*self.amperiod)
	    wave.SetParLimits(1,-numpy.pi,numpy.pi)
	    wave.FixParameter(2,self.amperiod)
	    wave.FixParameter(3,self.period)
	    wave.FixParameter(4,self.size[i]["radius"])
	    wave.FixParameter(5,self.size[i]["center"])
	    wave.FixParameter(6,self.fastclkrate)
	    h1.Fit("wave","Q")
	    self.amphase[i]=wave.GetParameter(0)
	    self.phase[i]=wave.GetParameter(1)
	    amphaseoffset=(clkoffset/self.fastclkrate)%(self.amperiod*4)
	    phaseoffset=2*numpy.pi*((clkoffset/self.fastclkrate)%self.period)/self.period
	    self.amphase[i]=self.amphase[i]-amphaseoffset
	    self.phase[i]=self.phase[i]-phaseoffset
	    #c1.Print("pic/%srastershape_%i_%i_%i_%i.png"%(self.rastertype,self.rasterchan[i],self.run,eventrange[0],eventrange[1]),"png")
	return self.amphase,self.phase
    
    def getphase2(self,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	eventrange0=[eventrange[0],eventrange[0]+2000]
	clkrange0=[self.clkraw[eventrange0[0]],self.clkraw[eventrange0[1]]]
	wave0=ROOT.TF1("wave0",self.__slrasterfunc,clkrange0[0],clkrange0[1],7)
	for i in range(0,2):
	    c1=ROOT.TCanvas("c1","%s raster phase"%self.rastertype,800,600)
	    #first fit, for phase
	    print "xy=%i,fit first time,for phase"%(i)
	    clkoffset=self.clkraw[eventrange0[0]]
	    clkrawoff=array("d",[a-clkoffset for a in self.clkraw[eventrange0[0]:eventrange0[1]]])
	    h0=ROOT.TGraph(eventrange0[1]-eventrange0[0],clkrawoff,self.rasterraw[self.rasterchan[i]][eventrange0[0]:eventrange0[1]])
	    #h0.Draw("AL")
	    wave0.SetParLimits(0,-2*self.amperiod,2*self.amperiod)
	    wave0.SetParLimits(1,-numpy.pi,numpy.pi)
	    wave0.FixParameter(2,self.amperiod)
	    wave0.FixParameter(3,self.period)
	    wave0.FixParameter(4,self.size[i]["radius"])
	    wave0.FixParameter(5,self.size[i]["center"])
	    wave0.FixParameter(6,self.fastclkrate)
	    h0.Fit("wave0")
	    self.amphase[i]=wave0.GetParameter(0)
	    self.phase[i]=wave0.GetParameter(1)
	    amphaseoffset=(clkoffset/self.fastclkrate)%(self.amperiod*4)
	    phaseoffset=2*numpy.pi*((clkoffset/self.fastclkrate)%self.period)/self.period
	    self.amphase[i]=self.amphase[i]-amphaseoffset
	    self.phase[i]=self.phase[i]-phaseoffset
	    #second fit,for fast clock rate
	    fitevents=1000
	    eventadd=10000
	    self.fitclockrate[i]=self.fastclkrate
	    eventrange1=[eventrange0[0]+eventadd ,eventrange0[0]+eventadd+fitevents]
	    if eventrange[1]<eventrange1[1]:eventrange1=[eventrange[1]-fitevents,eventrange[1]]
	    fittimes=0
	    raterange=100
	    trytimes=0
	    while(eventrange[1]-eventrange1[1]>2000 and eventrange1[1]<eventrange[0]+100000):
		print "xy=%i,fit %i time,for fast clock rate,eventrange:%i~%i,raterange:%f+-%f"%(i,fittimes+1,eventrange1[0],eventrange1[1],self.fitclockrate[i],raterange)
		clkrange1=[self.clkraw[eventrange1[0]],self.clkraw[eventrange1[1]]]
		wave1=ROOT.TF1("wave1",self.__slrasterfunc,clkrange1[0],clkrange1[1],7)
		clkoffset=self.clkraw[eventrange1[0]]
		clkrawoff=array("d",[a-clkoffset for a in self.clkraw[eventrange1[0]:eventrange1[1]]])
		h1=ROOT.TGraph(eventrange1[1]-eventrange1[0],clkrawoff,self.rasterraw[self.rasterchan[i]][eventrange1[0]:eventrange1[1]])
		#h1.Draw("AL")
		amphaseoffset=(clkoffset/self.fastclkrate)%(self.amperiod*4)
		phaseoffset=2*numpy.pi*((clkoffset/self.fastclkrate)%self.period)/self.period
		amphaseoff=self.amphase[i]+amphaseoffset
		phaseoff=self.phase[i]+phaseoffset
		wave1.FixParameter(0,amphaseoff)
		wave1.FixParameter(1,phaseoff)
		wave1.FixParameter(2,self.amperiod)
		wave1.FixParameter(3,self.period)
		wave1.FixParameter(4,self.size[i]["radius"])
		wave1.FixParameter(5,self.size[i]["center"])
		wave1.SetParLimits(6,self.fitclockrate[i]-raterange,self.fitclockrate[i]+raterange)
		h1.Fit("wave1")
		tmprate=wave1.GetParameter(6)
		tmperr=wave1.GetParError(6)
		if tmperr<0.001:
		    if trytimes>5:
			trytimes=0
			pass
		    else:
			raterange=(raterange/1.1 if raterange/1.1>1 else 1)
			trytimes+=1
			continue
		else:
		    self.fitclockrate[i]=tmprate
		raterange=(raterange/2.5 if raterange/2.5>1 else 1)
		fittimes+=1
		eventadd=eventadd*2
		eventrange1=[eventrange1[0]+eventadd,eventrange1[0]+eventadd+fitevents]
		if eventrange[1]<eventrange1[1]:eventrange1=[eventrange[1]-fitevents,eventrange[1]]
	#c1.Print("pic/%sraster_%i.png"%(self.rastertype,self.run),"png")
	#c1.WaitPrimitive()
	return self.amphase,self.phase,self.fitclockrate
    
    def getphaseall(self):
	if os.path.exists("aa.pkl"):
	    self.phaseall=zload("aa.pkl")
	    return
	step=2000
	self.eventpkl=self.pklpath.getpath("raw","event",self.run)
	self.event=zload(self.eventpkl)
	splitentries,splitevents,amphases,phases=[],[],[],[]
	for i in range(0,self.events,step):
	    istep=i+step if i+step<self.events else self.events-1
	    print "fit for slow raster shape %i-%i,totally %i events"%(i,istep,self.events)
	    clkrange=[self.clkraw[i],self.clkraw[istep-1]]
	    clkoffset=self.clkraw[i]
	    clkrawoff=array("d",[a-clkoffset for a in self.clkraw[i:istep]])
	    wave=ROOT.TF1("wave",self.__slrasterfunc,clkrange[0],clkrange[1],7)
	    amphase,phase=[0,0],[0,0]
	    for j in range(0,2):
		h1=ROOT.TGraph(step,clkrawoff,self.rasterraw[self.rasterchan[j]][i:istep])
		wave.SetParLimits(0,-2*self.amperiod,2*self.amperiod)
		wave.SetParLimits(1,-numpy.pi,numpy.pi)
		wave.FixParameter(2,self.amperiod)
		wave.FixParameter(3,self.period)
		wave.FixParameter(4,self.size[j]["radius"])
		wave.FixParameter(5,self.size[j]["center"])
		wave.FixParameter(6,self.fastclkrate)
		fitresult=h1.Fit("wave","Q")
		if int(fitresult)>0:
		    amphase[j],phase[j]=False,False
		    continue
		amphase[j]=wave.GetParameter(0)
		phase[j]=wave.GetParameter(1)
		#print i,amphase[j],phase[j],int(fitresult)
		amphaseoffset=(clkoffset/self.fastclkrate)%(self.amperiod*4)
		phaseoffset=2*numpy.pi*((clkoffset/self.fastclkrate)%self.period)/self.period
		amphase[j]-=amphaseoffset
		phase[j]-=phaseoffset
		if amphase[j]<-2*self.amperiod:amphase[j]+=4*self.amperiod
		if phase[j]<-numpy.pi:phase[j]+=2*numpy.pi
	    splitentries.append([i,istep])
	    splitevents.append([self.event[i],self.event[istep]])
	    amphases.append(amphase)
	    phases.append(phase)
	self.phaseall={"splitentries":splitentries,"splitevents":splitevents,"amphase":amphases,"phase":phases}
	zdump(self.phaseall,"aa.pkl")
	    
    def getrebuiltvalue(self,xy,clkvalue):
	#get rebuilt adc raster value for each event
	return self.__slrasterfunc([clkvalue],[self.amphase[xy],self.phase[xy],self.amperiod,self.period,self.size[xy]["radius"],self.size[xy]["center"],self.fitclockrate[xy]])

    def checkphase(self,xy,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	clkrange=[self.clkraw[eventrange[0]],self.clkraw[eventrange[1]-1]]
	c1=ROOT.TCanvas("c1","%s phase check"%self.rastertype,800,600)
	g1=ROOT.TGraph(eventrange[1]-eventrange[0],self.clkraw[eventrange[0]:eventrange[1]],self.rasterraw[xy+self.rasterchan[0]][eventrange[0]:eventrange[1]])
	g1.Draw("AL")
	for i in range(len(self.phaseall["splitentries"])):
	    if eventrange[0] in range(self.phaseall["splitentries"][i][0],self.phaseall["splitentries"][i][1]):
		amph=self.phaseall["amphase"][i]
		ph=self.phaseall["phase"][i]
		print i,self.phaseall["splitentries"][i],amph[xy],ph[xy]
		break
	if not amph[xy]:return
	#amph=self.amphase
	#ph=self.phase
	wave=ROOT.TF1("wave",self.__slrasterfunc,clkrange[0],clkrange[1],7)
	wave.SetParameters(amph[xy],ph[xy],self.amperiod,self.period,self.size[xy]["radius"],self.size[xy]["center"],float(self.fastclkrate))
	wave.Draw("A*same")
	#c1.Print("pic/%srasterphasecheck_%i_%i.png"%(self.rastertype,self.run,i),"png")
	c1.Print("pic/%srasterphasecheck_%i.png"%(self.rastertype,self.run),"png")
	del c1,wave
	#c1.WaitPrimitive()
	
class fastrecon:
    def __init__(self):
	rasterrecon.__init__(self)
	self.period=1/2870.
	self.phase=[0,0]
	self.shape=ROOT.TF1("shape",self.__rectangle,0,9000,3) #rectangle
	self.rastertype="fast"
	self.rasterchan=[0,1]
	
    def __rectangle(self,x,par):
	#par:center,radius,amplitude
	ledge=par[0]-par[1]
	redge=par[0]+par[1]
	if x[0]<ledge or x[0]>redge:return 0
	elif x[0]>ledge and x[0]<redge:return par[2]
	elif x[0]==ledge or x[0]==redge:return par[2]/2.
	else:return 0
	
    #general get size function
    def gengetsize(self,chan,eventrange=[0,1000000000]):
	#multitimes fit:
	def multifits(hist,fitfun,c,r,a,xmin,xmax,hmax):
	    #two times fit
	    c,r,a=fitfun(hist,c-r*0.05,c+r*0.05,r*0.5,r,a*0.5,a,xmin,xmax)
	    c,r,a=fitfun(hist,c-r*0.05,c+r*0.05,r*0.8,r*2,a*0.8,a*1.2,xmin,xmax)
	    return c,r,a
	return rasterrecon.gengetsize(self,multifits,chan,eventrange)
	
    def __fastrasterfunc(self,x,par):
	#par:phase,period,downratio/upratio,max amp,center,fastclkrate
	if not (par[1] and par[2] and par[5]):return -10000000
	time=x[0]/par[5]
	dp=(time+par[0])%par[1]
	splittimeup=par[1]/(1.+1./par[2])
	splittimedown=splittimeup+(par[1]-splittimeup)/2.
	if dp<=splittimeup:tria=dp-splittimeup/2.
	else:tria=par[2]*(splittimedown-dp)
	return par[3]/par[1]*tria+par[4]
    
    def getphase(self,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	clkrange=[self.clkraw[eventrange[0]],self.clkraw[eventrange[1]-1]]
	clkoffset=self.clkraw[eventrange[0]]
	clkrawoff=array("d",[a-clkoffset for a in self.clkraw[eventrange[0]:eventrange[1]]])
	eventcut="Entry$>=%i&&Entry$<=%i"%tuple(eventrange)
	wave=ROOT.TF1("wave",self.__fastrasterfunc,clkrange[0],clkrange[1],6)
	for i in range(0,1):
	    c1=ROOT.TCanvas("c1","%s raster phase"%self.rastertype,800,600)
	    h1=ROOT.TGraph(eventrange[1]-eventrange[0],clkrawoff,self.rasterraw[self.rasterchan[i]][eventrange[0]:eventrange[1]])
	    wave.SetParLimits(0,-self.period/2.,self.period/2.)
	    wave.SetParLimits(1,self.period*0.9,self.period*1.3)
	    #wave.SetParLimits(2,0.8,3)
	    wave.FixParameter(2,1)
	    #wave.SetParLimits(3,self.size[i]["radius"]*0.8,self.size[i]["radius"]*2)
	    wave.FixParameter(3,self.size[i]["radius"])
	    wave.FixParameter(4,self.size[i]["center"])
	    wave.FixParameter(5,self.fastclkrate)
	    h1.Fit("wave")
	    self.phase[i]=wave.GetParameter(0)
	    self.period=wave.GetParameter(1)
	    phaseoffset=(clkoffset/self.fastclkrate)%(self.period)
	    c1.Print("pic/%sraster_%i.png"%(self.rastertype,self.run),"png")
	#c1.WaitPrimitive()
	
    def checkphase(self,xy,eventrange=[0,1000000000]):
	ROOT.gROOT.SetBatch(True)
	eventrange[1]=eventrange[1] if eventrange[1]<self.events else self.events
	clkrange=[self.clkraw[eventrange[0]],self.clkraw[eventrange[1]-1]]
	c1=ROOT.TCanvas("c1","%s phase check"%self.rastertype,800,600)
	g1=ROOT.TGraph(eventrange[1]-eventrange[0],self.clkraw[eventrange[0]:eventrange[1]],self.rasterraw[xy+self.rasterchan[0]][eventrange[0]:eventrange[1]])
	g1.Draw("AL")
	wave=ROOT.TF1("wave",self.__fastrasterfunc,clkrange[0],clkrange[1],6)
	wave.SetParameters(self.phase[xy],self.period,1,self.size[xy]["radius"],self.size[xy]["center"],self.fastclkrate)
	wave.Draw("A*same")
	c1.Print("pic/%sraster_%i_%i_%i.png"%(self.rastertype,self.run,eventrange[0],eventrange[1]),"png")
	del c1,wave
	#c1.WaitPrimitive()
	
class rastercalib:
    #general raster size,parameter: raster class, raster x,y magnet z pos,...
    def __rastersizecalib(self,raster,rasterz,rastertype,runpaths,eventsplits):
	ROOT.gROOT.SetBatch(True)
	#get run info
	para=bpminsertparaprep(runpaths[0])
	survey=para["survey"]
	self.run=para["run"]
	if not hasattr(self,'tgtz'):
	    self.tgtz=para["tgtz"]
	#get data
	size,sizeinbpm,g1=[],[],[]
	datanum=0
	rasters=[]
	for i in range(len(runpaths)):
	    if i>=len(eventsplits):break
	    rasters.append(raster)
	    rasters[-1].runinit(runpaths[i])
	    for s in eventsplits[i]:
		s=[int(a) for a in s]
		size.append(rasters[-1].getsize(s))
		sizeinbpm.append(rasters[-1].getsizeinbpm(s))
		datanum+=1
	rx,ry,bax,bay,bbx,bby=array('d'),array('d'),array('d'),array('d'),array('d'),array('d')
	rxc,ryc,baxc,bayc,bbxc,bbyc=array('d'),array('d'),array('d'),array('d'),array('d'),array('d')
	for i in range(datanum):
	    rxtmp,rytmp=size[i][0]["radius"],size[i][1]["radius"]
	    baxtmp,baytmp=sizeinbpm[i][0][0]["radius"],sizeinbpm[i][0][1]["radius"]
	    bbxtmp,bbytmp=sizeinbpm[i][1][0]["radius"],sizeinbpm[i][1][1]["radius"]
	    rxc.append(size[i][0]["center"])
	    ryc.append(size[i][1]["center"])
	    if any([abs(baxtmp)>50,abs(baytmp)>50,abs(bbxtmp)>50,abs(bbytmp)>50]):
		continue
	    rx.append(rxtmp)
	    ry.append(rytmp)
	    bax.append(baxtmp)
	    bay.append(baytmp)
	    bbx.append(bbxtmp)
	    bby.append(bbytmp)
	    #baxc.append(sizeinbpm[i][0][0]["center"])
	    #bayc.append(sizeinbpm[i][0][1]["center"])
	    #bbxc.append(sizeinbpm[i][1][0]["center"])
	    #bbyc.append(sizeinbpm[i][1][1]["center"])
	rastercenter=[sum(rxc)/len(rxc),sum(ryc)/len(ryc)]
	# get raster size in target
	#bpma=bpmpos(survey[0])
	#bpmb=bpmpos(survey[1])
	#for i in range(len(eventsplits)):
	    #if runorbit<=0:
		#rasteredge_axp=bpma.posbpmhall2hall([baxc[i]+bax[i],bayc[i]])
		#rasteredge_axm=bpma.posbpmhall2hall([baxc[i]-bax[i],bayc[i]])
		#rasteredge_ayp=bpma.posbpmhall2hall([baxc[i],bayc[i]+bay[i]])
		#rasteredge_aym=bpma.posbpmhall2hall([baxc[i],bayc[i]-bay[i]])
		#rasteredge_bxp=bpmb.posbpmhall2hall([bbxc[i]+bbx[i],bbyc[i]])
		#rasteredge_bxm=bpmb.posbpmhall2hall([bbxc[i]-bbx[i],bbyc[i]])
		#rasteredge_byp=bpmb.posbpmhall2hall([bbxc[i],bbyc[i]+bby[i]])
		#rasteredge_bym=bpmb.posbpmhall2hall([bbxc[i],bbyc[i]-bby[i]])
		#tgtpos_xp=tgtpos_straight(rasteredge_axp,rasteredge_bxp,[0]*3,[0]*3,self.tgtz)
		#tgtpos_xm=tgtpos_straight(rasteredge_axm,rasteredge_bxm,[0]*3,[0]*3,self.tgtz)
		#tgtpos_yp=tgtpos_straight(rasteredge_ayp,rasteredge_byp,[0]*3,[0]*3,self.tgtz)
		#tgtpos_ym=tgtpos_straight(rasteredge_aym,rasteredge_bym,[0]*3,[0]*3,self.tgtz)
	    #else:
		#rasteredge_axp=bpma.posbpmhall2rot([baxc[i]+bax[i],bayc[i]])
		#rasteredge_axm=bpma.posbpmhall2rot([baxc[i]-bax[i],bayc[i]])
		#rasteredge_ayp=bpma.posbpmhall2rot([baxc[i],bayc[i]+bay[i]])
		#rasteredge_aym=bpma.posbpmhall2rot([baxc[i],bayc[i]-bay[i]])
		#rasteredge_bxp=bpmb.posbpmhall2rot([bbxc[i]+bbx[i],bbyc[i]])
		#rasteredge_bxm=bpmb.posbpmhall2rot([bbxc[i]-bbx[i],bbyc[i]])
		#rasteredge_byp=bpmb.posbpmhall2rot([bbxc[i],bbyc[i]+bby[i]])
		#rasteredge_bym=bpmb.posbpmhall2rot([bbxc[i],bbyc[i]-bby[i]])
		#tgtpos_xp=tgtpos_field(rasteredge_axp,rasteredge_bxp,[0]*3,[0]*3,self.tgtz,posmapfun)
		#tgtpos_xm=tgtpos_field(rasteredge_axm,rasteredge_bxm,[0]*3,[0]*3,self.tgtz,posmapfun)
		#tgtpos_yp=tgtpos_field(rasteredge_ayp,rasteredge_byp,[0]*3,[0]*3,self.tgtz,posmapfun)
		#tgtpos_ym=tgtpos_field(rasteredge_aym,rasteredge_bym,[0]*3,[0]*3,self.tgtz,posmapfun)
	    #print rasteredge_axp,rasteredge_bxp,tgtpos_xp[0][1]
	    #print rasteredge_axm,rasteredge_bxm,tgtpos_xm[0][1]
	    #print rasteredge_ayp,rasteredge_byp,tgtpos_yp[0][1]
	    #print rasteredge_aym,rasteredge_bym,tgtpos_ym[0][1]
	    #print rasteredge_axp[0]-rasteredge_axm[0],rasteredge_bxp[0]-rasteredge_bxm[0],tgtpos_xp[0][1][0]-tgtpos_xm[0][1][0]
	    #print rasteredge_ayp[1]-rasteredge_aym[1],rasteredge_byp[1]-rasteredge_bym[1],tgtpos_yp[0][1][1]-tgtpos_ym[0][1][1]
	#raster x y vs bpm x y
	g1.append(ROOT.TGraph(len(rx),rx,bax))
	g1.append(ROOT.TGraph(len(ry),ry,bay))
	g1.append(ROOT.TGraph(len(rx),rx,bbx))
	g1.append(ROOT.TGraph(len(ry),ry,bby))
	c1=ROOT.TCanvas("c1","%s raster size calibration"%rastertype,1280,800)
	c1.Divide(2,2)
	pol1=ROOT.TF1("pol1","pol1")
	slope,ped,rbpmvstgt,sizeconst=[],[],[],[]
	for i in range(4):
	    c1.cd(i+1)
	    g1[i].Draw("A*")
	    g1[i].Fit("pol1")
	    ped.append(pol1.GetParameter(0))
	    slope.append(pol1.GetParameter(1))
	    #ratio of raster size in bpm vs in target
	    rbpmvstgt.append([])
	    sizeconst.append([])
	    for j in range(len(self.tgtz)):
		rbpmvstgt[-1].append((rasterz[i%2]+self.tgtz[j])/(rasterz[i%2]-survey[int(i/2)][0][2]))
		sizeconst[-1].append([rbpmvstgt[-1][-1]*slope[i],rbpmvstgt[-1][-1]*ped[i],rastercenter[i%2]])
	c1.Print("pic/%scalib_%i.png"%(rastertype,self.run),"png")
	#only use bpm A right now,will edit after improve bpm resolution
	return sizeconst[:2]
	
    def slowrastersizecalib(self,runpaths,eventsplits):
	rasterz=[20888.667,20059.467] #slow raster x,y magnet z position(g2p hall coordinate,from yves' orbit file)
	self.slowsizeconst=self.__rastersizecalib(slrecon(),rasterz,"slow",runpaths,eventsplits)
	
    def fastrastersizecalib(self,runpaths,eventsplits):
	rasterz=[22122.967,21622.967] #fast raster x,y magnet z position(g2p hall coordinate,from yves' orbit file)
	self.fastsizeconst=self.__rastersizecalib(fastrecon(),rasterz,"fast",runpaths,eventsplits)
	
    def rasterconstsave(self):
	dbdir=os.getenv("BEAMDBPATH")
	if dbdir==None:
	    print "please define BEAMDBPATH in your env"
	    return False
	pydb=os.path.join(dbdir,"pyDB")
	if not os.path.exists(pydb):os.makedirs(pydb)
	filename=os.path.join(pydb,"raster_%i.dat"%(self.run))
	period=runinfo()
	runrange=period.runrange(self.run)
	datfile=open(filename,"w")
	datfile.write("raster calibration constant database, including size, fastclock\n")
	datfile.write("\n")
	datfile.write("------------------------------------------------\n\n")
	datfile.write("avail run period:%i-%i\n"%tuple(runrange))
	datfile.write("fastclock rate:103920\n")
	datfile.write("#set to 1 if want to rebuild raster function instead of using ADC value directly, first is fast raster,second is slow raster\n")
	datfile.write("function rebuild:0 0\n")
	datfile.write("target z position(mm,dont change it):")
	for z in self.tgtz: datfile.write("%f "%z)
	datfile.write("\n")
	if hasattr(self,"slowsizeconst"):
	    datfile.write("slow raster x slope:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.slowsizeconst[0][i][0])
	    datfile.write("\n")
	    datfile.write("slow raster x ped:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.slowsizeconst[0][i][1])
	    datfile.write("\n")
	    datfile.write("slow raster y slope:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.slowsizeconst[1][i][0])
	    datfile.write("\n")
	    datfile.write("slow raster y ped:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.slowsizeconst[1][i][1])
	    datfile.write("\n")
	    datfile.write("slow raster x center:%i\n"%self.slowsizeconst[0][0][2])
	    datfile.write("slow raster y center:%i\n"%self.slowsizeconst[1][0][2])
	if hasattr(self,"fastsizeconst"):
	    datfile.write("fast raster x slope:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.fastsizeconst[0][i][0])
	    datfile.write("\n")
	    datfile.write("fast raster x ped:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.fastsizeconst[0][i][1])
	    datfile.write("\n")
	    datfile.write("fast raster y slope:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.fastsizeconst[1][i][0])
	    datfile.write("\n")
	    datfile.write("fast raster y ped:")
	    for i in range(len(self.tgtz)):datfile.write("%f "%self.fastsizeconst[1][i][1])
	    datfile.write("\n")
	    datfile.write("fast raster x center:%f\n"%self.fastsizeconst[0][0][2])
	    datfile.write("fast raster y center:%f\n"%self.fastsizeconst[1][0][2])
    
if __name__ == '__main__':
    sl=slrecon()
    runpath=os.path.join(os.getenv("REPLAY_OUT_PATH"),"bpm_5555.root")
    sl.runinit(runpath)
    sl.getsize()
    sl.getphase([1000,3000])
    sl.checkphase(0,[7000,7100])