#!/usr/bin/env python
import os,sys,numpy
try:import cPickle as pickle
except:import pickle
from runinfo import calibsetting

class bpmfit:
    def __init__(self,keywords,xy,real,data,centerid):
	self.realpos,self.datax,self.datay=[],[],[]
	#self.ddatax,self.ddatay=[],[]
	for i in range(len(data[0])):
	    self.realpos+=[real[i]]*len(data[0][i])
	    if i==0:
		self.datax=data[0][i]
		self.datay=data[1][i]
		#self.ddatax=ddata[0][i]
		#self.ddatay=ddata[1][i]
	    else:
		self.datax=numpy.concatenate((self.datax,data[0][i]))
		self.datay=numpy.concatenate((self.datay,data[1][i]))
		#self.ddatax=numpy.concatenate((self.ddatax,ddata[0][i]))
		#self.ddatay=numpy.concatenate((self.ddatay,ddata[1][i]))
	self.realpos=numpy.asarray(self.realpos,dtype=numpy.float32)
	self.centerid=centerid
	#get calibsetting
	self.calibconf=calibsetting(keywords)
	self.fitorder=self.calibconf.fitorder
	if xy>0:self.fitorder=[self.fitorder[1],self.fitorder[0]]
	self.maxorder=max(self.calibconf.fitorder)
	
    def fitpy(self,fun,initvalue=False,parswitch=False):
	#use pyminuit2 package,also can use pyminuit package
	from minuit2 import Minuit2 as mnt
	if not parswitch:parswitch=self.calibconf.switch
	m=mnt(fun)
	#get total var number without parswitch
	paranum=0
	masker=0
	for i in range(self.maxorder+1):
	    for j in range(i+1):
		if i-j<=self.fitorder[0] and j<=self.fitorder[1]:
		    masker+=pow(2,paranum)
		    paranum+=1
	parswitch=parswitch&masker
	nvars=0 #vars need to fit
	for i in range(paranum):
	    if (parswitch>>i)&0x1: nvars+=1
	totns=self.calibconf.scantimes #scan times
	if totns:ns=int(pow(totns,1./nvars))
	sr={} #scanrange
	defaultinitvalue=initvalue
	if not defaultinitvalue:initvalue={}
	for i in range(21):
	    fitid="a%.2i"%i
	    m.errors[fitid]=self.calibconf.steps[i] if i<len(self.calibconf.steps) else 0.0001
	    if totns:
		sr[fitid]=(fitid,ns,self.calibconf.parlow[i],self.calibconf.parhigh[i]) if i<len(self.calibconf.parlow) else (fitid,ns,-0.01,0.01)
	    if not defaultinitvalue or fitid not in initvalue.keys():
		initvalue[fitid]=self.calibconf.initvalue[i] if i<len(self.calibconf.initvalue) else 0
	m.values=initvalue
	scanparams=""
	for i in range(21):
	    if (parswitch>>i)&0x1: 
		#print "set %ird parameter to fit"%i
		m.fixed["a%.2i"%i]=False
		scanparams+="('%s',%i,%f,%f),"%sr["a%.2i"%i]
	    else:
		#print "fixed %ird parameter"%i
		m.fixed["a%.2i"%i]=True
	#scan,improve minimization
	#print "m.scan(%soutput=False)"%scanparams
	print "scan all fit range first...",
	if totns:eval("m.scan(%soutput=False)"%scanparams)
	print "done"
	try:m.migrad()
	except Exception as err:
	    print err
	return m.values,m.errors,numpy.sqrt(m.fval)
    
    def fitresult2list(self,result):
	#sort,convert dict to list
	f,e=[],[]
	for w in sorted(result[0].keys()):
	    f.append(result[0][w])
	    e.append(result[1][w])
	totpara=0
	for i in range(len(f)):
	    if sum(f[i:])==0:break
	    totpara+=1
	return f[:totpara],e[:totpara],result[2]
    
    def fun(self,a00,a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, a11,a12,a13,a14,a15,a16,a17,a18,a19,a20):
	calpos=0
	paraid=0
	for i in range(self.maxorder+1):
	    for j in range(i+1):
		if i-j<=self.fitorder[0] and j<=self.fitorder[1]:
		    calpos=calpos+vars()["a%.2i"%paraid]*pow(self.datax,i-j)*pow(self.datay,j)
		    paraid+=1
	delta=(self.realpos-calpos)
	chisq=numpy.mean(delta*delta)
	#print chisq
	return chisq
    
    def fita(self):
	#fit
	firstresult=self.fitpy(self.fun,parswitch=7)[0]
	secondresult=self.fitpy(self.fun,firstresult,parswitch=self.calibconf.switch)
	result=secondresult
	return self.fitresult2list(result)
	
    def funm(self,a00,a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, a11,a12,a13,a14,a15,a16,a17,a18,a19,a20):
	#print a,b,c,d,
	calpos=self.centerx
	paraid=0
	for i in range(self.maxorder+1):
	    for j in range(i+1):
		if i-j<=self.fitorder[0] and j<=self.fitorder[1]:
		    if i-j!=0 or j!=0:
			calpos=calpos+vars()["a%.2i"%paraid]*(pow(self.datax,i-j)*pow(self.datay,j)\
			-pow(self.centerdata[0],i-j)*pow(self.centerdata[1],j))
		    paraid+=1
	delta=(self.realpos-calpos)
	chisq=numpy.mean(delta*delta)
	#print chisq
	return chisq
	
    #force set center point match center harp pos 
    def fit(self):
	#get center pos
	self.centerdata=[self.datax[self.centerid],self.datay[self.centerid]]
	self.centerx=self.realpos[self.centerid]
	self.realpos=numpy.delete(self.realpos,self.centerid)
	self.datax=numpy.delete(self.datax,self.centerid)
	self.datay=numpy.delete(self.datay,self.centerid)
	#firstresult=self.fitpy(initvalue={'a00':self.centerx},parswitch=6)[0]
	#secondresult=self.fitpy(initvalue=firstresult,parswitch=self.calibconf.switch&1048574)
	#result=secondresult
	result=self.fitpy(self.funm,parswitch=self.calibconf.switch&1048574)
	#get c
	result[0]["a00"]=self.centerx
	result[1]["a00"]=0
	paraid=0
	for i in range(self.maxorder+1):
	    for j in range(i+1):
		if i-j<=self.fitorder[0] and j<=self.fitorder[1]:
		    if i-j!=0 or j!=0:
			result[0]["a00"]-=result[0]["a%.2i"%paraid]*pow(self.centerdata[0],i-j)*pow(self.centerdata[1],j)
			result[1]["a00"]+=result[1]["a%.2i"%paraid]*result[1]["a%.2i"%paraid]
		    paraid+=1
	result[1]["a00"]=numpy.sqrt(result[1]["a00"])
	return self.fitresult2list(result)
 
	


