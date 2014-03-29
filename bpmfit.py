#!/usr/bin/env python
import os,sys
try:import cPickle as pickle
except:import pickle
from ROOT import TMinuit,TVirtualFitter,Long,Double,TTree,TCanvas,gROOT,TF1,TBranch,TGraph,gPad,kFALSE,NULL
from numpy import sqrt,mean,arange,matrix,log,asarray
from array import array
from pylab import plot,subplot,hist,savefig,figure,ylim

def dominuit(fcn,Npar,startvalues,step,parlimitslo,parlimitshi,parswitch):
    pMinuit=TMinuit(Npar)
    pMinuit.SetFCN(fcn)
    arglist=array('d',10*[0.])
    ierflg=Long(0)
    #Set print level
    #SET PRIntout  <level>
    #Sets the print level, determining how much output will be
    #produced. Allowed values and their meanings are displayed
    #after a SHOw PRInt command, and are currently <level>=:
    #[-1]  no output except from SHOW commands
    #[0]  minimum output
    #[1]  default value, normal output
    #[2]  additional output giving intermediate results.
    #[3]  maximum output, showing progress of minimizations.
    arglist[0]=-1
    pMinuit.mnexcm("SET PRIntout",arglist,1,ierflg)
    #Set warnings
    #-SET NOWarnings-
    #Supresses Minuit warning messages.
    #-SET WARnings-
    #Instructs Minuit to output warning messages when suspicious
    #conditions arise which may indicate unreliable results.
    #This is the default.
    arglist[0]=0
    pMinuit.mnexcm("SET NOWarnings",arglist,0,ierflg)
    #Set starting values and step sizes for parameters
    for i in range(Npar):
	pMinuit.mnparm(i,"Par%i"%(i+1),Double(startvalues[i]),Double(step[i]),0,0,ierflg)
    arglist[0]=1
    pMinuit.mnexcm("SET ERR",arglist,1,ierflg)
    #Set limits for parameters
    for i in range(Npar):
	if (parswitch>>i)&0x1:
	    arglist[0]=i+1
	    arglist[1]=parlimitslo[i]
	    arglist[2]=parlimitshi[i]
	    pMinuit.mnexcm("SET LIMits",arglist,3,ierflg)
	else:
	    arglist[0]=i+1
	    pMinuit.mnexcm("FIX",arglist,1,ierflg)
    #Now ready for minimization step
    arglist[0] = 500# Number of calls to FCN before giving up
    arglist[1] = 1. #tolerance
    pMinuit.mnexcm("MIGRAD",arglist,2,ierflg)
    #Store the fitting result
    dumx,dumxerr=Double(0),Double(0)
    fitparams,fiterrors=[],[]
    for i in range(Npar):
	pMinuit.GetParameter(i,dumx, dumxerr)
	fitparams.append(float(dumx))
	fiterrors.append(float(dumxerr))
    #Print results
    amin,edm,errdef=Double(0.18),Double(0.19),Double(0.20)
    nvpar,nparx,icstat=Long(1),Long(2),Long(3)
    pMinuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat)
  #  pMinuit.mnprin(3,amin)
    del pMinuit
    return fitparams,fiterrors

def dominuit2(fcn,Npar,startvalues,step,parlimitslo,parlimitshi,parswitch):
    TVirtualFitter.SetDefaultFitter("Minuit")
    fitter=TVirtualFitter.Fitter(NULL,Npar)
    fitter.SetFCN(fcn)
    #starting values
    for i in range(Npar):
	if (parswitch>>i)&0x1:
	    fitter.SetParameter(i,"Par%i"%(i+1),Double(startvalues[i]),Double(step[i]),Double(parlimitslo[i]),Double(parlimitshi[i]))
	else:
	    fitter.FixParameter(i)
    arglist=array('d',10*[0.])
    arglist[0]=1
    fitter.ExecuteCommand("SET PRIntout",arglist,1)
    arglist[0] = 500# Number of calls to FCN before giving up
    arglist[1] = 0.01 #tolerance
    fitter.ExecuteCommand("MINIMIZE", arglist, 2)
    fitparams,fiterrors=[],[]
    for i in range(Npar):
	fitparams.append(float(fitter.GetParameter(i)))
	fiterrors.append(float(fitter.GetParError(i)))
    return fitparams,fiterrors
    

class bpmfit:
    def __init__(self,real,data,ddata):
	self.realpos=real
	self.datax,self.datay=data[0],data[1]
	self.ddatax,self.ddatay=ddata[0],ddata[1]
    def minuitfun(self,npar,gin,f,par,iflag):
	chisq,delta=0,0
	for i in range(len(self.datax)):
	    x,y=self.datax[i],self.datay[i]
	    dx,dy=self.ddatax[i],self.ddatay[i]
	    realp=self.realpos[i]
	    calpos=par[0]*x+par[1]*y+par[2]
	    dcalpos=sqrt(dx*dx+dy*dy)
	    delta=(realp-calpos)/dcalpos
	    chisq+=delta**2
	f[0]=chisq
    def minuitfun2(self,npar,gin,f,par,iflag):
	chisq,delta=0,0
	for i in range(len(self.datax)):
	    x,y=self.datax[i],self.datay[i]
	    dx,dy=self.ddatax[i],self.ddatay[i]
	    realp=self.realpos[i]
	    calpos=par[0]*x+par[1]*y+par[2]+par[3]*x*x+par[4]*y*y
	    dcalpos=sqrt(dx*dx+dy*dy)
	    delta=(realp-calpos)
	    chisq+=delta**2
	f[0]=chisq    
    def fit(self):
	#if len(self.datax)<5:
	return dominuit2(self.minuitfun,3,[-2,-2,-3],[0.0001]*3,[-2,-2,-3],[2,2,3],7)
	#else:
	 #   return dominuit(self.minuitfun2,5,[-2,-2,-3,-1,-1],[0.0001]*5,[-2,-2,-5,-1,-1],[2,2,5,1,1],31)
	
    
class gxgyfit:
    def __init__(self,run,rootpath,bpmapos_bpm,bpmbpos_bpm,peaks,pedpeaks):
	gROOT.SetBatch(True)
	self.bpmpos=map(lambda p1,p2:p1[:2]+p2[:2],bpmapos_bpm,bpmbpos_bpm)
	self.peaks=map(lambda p0:map(lambda p1,p2:p1-p2,p0[0],pedpeaks[0]),peaks)
	pkldir=os.path.join(rootpath,"pkl")
	#get raw data for bpm
	self.rawbpm=pickle.load(open(os.path.join(pkldir,"bpmraw_sbpm_%i.pkl"%run),"rb"))
	self.bpmavail=pickle.load(open(os.path.join(pkldir,"bpmraw_bpmavail_%i.pkl"%run),"rb"))
	self.numdata=len(self.rawbpm[0])
	self.run=run
	for i in range(8):
	    self.rawbpm[i]=(self.rawbpm[i]-pedpeaks[0][i])*self.bpmavail
	self.rawbpm=[a[a>0] for a in self.rawbpm]
	    
    def minuitfun(self,npar,gin,f,par,iflag):
	rms=sqrt(mean((self.rawbpm[self.chan[0]]-par[0]*self.rawbpm[self.chan[1]])/(self.rawbpm[self.chan[0]]+par[0]*self.rawbpm[self.chan[1]]))/self.numdata)
	f[0]=rms
    
    def diffsum(self,a0,a1,g=1):
	return (a0-g*a1)/(a0+g*a1)
    
    def getgxgy(self,i,chan):
	coordischan,coordis=self.getmaxposdis(i)
	#for the ratio of diff gxgy(pos=g1*diff/sum+g2)
	coorpeak=map(lambda c0:map(lambda c1:self.peaks[c0][c1],chan),coordischan)
	def getratio(g):
	    r=matrix([[self.diffsum(coorpeak[0][0],coorpeak[0][1],g),1],\
		[self.diffsum(coorpeak[1][0],coorpeak[1][1],g),1]]).getI()\
		    *matrix(coordis).getT()
	    return r.getT().tolist()[0]
	#get nearest calbration run -- get most accurate calibration const
	ratio=getratio(1)
	p=mean(self.diffsum(self.rawbpm[chan[0]],self.rawbpm[chan[1]])*ratio[0]+ratio[1])
	coordischan,coordis=self.getbestposdis(i,p)
	coorpeak=map(lambda c0:map(lambda c1:self.peaks[c0][c1],chan),coordischan)
	
	fig=figure()
	subplot(211)
	plot(self.rawbpm[chan[0]])
	subplot(212)
	plot(self.rawbpm[chan[1]])
	savefig("pic/raw%i_%i_%i.png"%(self.run,chan[0],chan[1]))
	print max(self.rawbpm[chan[0]])-min(self.rawbpm[chan[0]]),max(self.rawbpm[chan[1]])-min(self.rawbpm[chan[1]])
	
	#get rms for diff gxgy
	rms,gg=array('d',[]),array('d',[])
	glow,ghigh,step=0.5,3,0.01
	for g in arange(glow,ghigh,step):
	    ratio=getratio(g)
	    p=self.diffsum(self.rawbpm[chan[0]],self.rawbpm[chan[1]],g)*ratio[0]+ratio[1]
	    fig=figure()
	    plot(p)
	    ylim([0,2])
	    savefig("pic/gxypos%i_%i_%f.png"%(self.run,i,g))
	    m=mean(p)
	    #p=p/m
	    r=sqrt(mean(((p-m))**2)/len(p))
	    rms.append(r)
	    gg.append(g)
	    #sys.exit()
	g1=TGraph(len(gg),gg,rms)
	gxgy=gg[asarray(rms).argmin()]
	c1=TCanvas("c1","gxgy",1024,768)
	#c1.SetLogx()
	g1.Draw("AL")
	#g1.GetXaxis().SetMoreLogLabels()
	#g1.GetXaxis().SetNoExponent()
	c1.Print("pic/gxgy%i_%i.png"%(self.run,i))
	#c1.WaitPrimitive()
	return gxgy
    
    def getmaxposdis(self,c):
	#find the max pos distance from existed bpm calibration run, c is 0~3,which is bpma x,y;bpm b x,y
	maxdischan,maxdis2=[0,1],[0,0]
	maxdis=0
	for i in range(len(self.bpmpos)):
	    for j in range(i+1,len(self.bpmpos)):
		dis=abs(self.bpmpos[i][c]-self.bpmpos[j][c])
		if dis>maxdis:
		    maxdischan=[i,j]
		    maxdis=dis
		    maxdis2=[self.bpmpos[i][c],self.bpmpos[j][c]]
	#return:which 2 runs have maximum position distance
	return maxdischan,maxdis2
    
    def getbestposdis(self,c,p):
	#find the best pos group to get calibration const
	dischan,dis2=[0,1],[0,0]
	#find the most close point
	dis=1000
	for i in range(len(self.bpmpos)):
	    d=abs(self.bpmpos[i][c]-p)
	    if d<dis:
		dischan[0]=i
		dis=d
		dis2[0]=self.bpmpos[i][c]
	#find the max distance point with the most close one
	dis=0
	for i in range(len(self.bpmpos)):
	    d=abs(self.bpmpos[i][c]-dis2[0])
	    if d>dis:
		dischan[1]=i
		dis=d
		dis2[1]=self.bpmpos[i][c]
	return dischan,dis2
    
    def fit(self):
	chans=[[2,3],[0,1],[6,7],[4,5]]
	gxy,dgxy=[],[]
	for i in range(0,len(chans)):
	    gxy.append(self.getgxgy(i,chans[i]))
	    sys.exit()
	return gxy
	
#if __name__ == '__main__':self.minuitfun
#    dominuit(minuitfun,5,[3,1,0.1,0.01,0.001],[0.1,0.1,0.01,0.001,0,0001],[0,0,0,0,0],[10,10,10,10,10],14)






