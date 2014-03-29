#!/usr/bin/env python
import re,urllib,os,pickle
from odsread import odsread

#find runinfo for a run,read from period.ods and runcurr.pkl
class runinfo(odsread):
    def __init__(self,filename="period.ods"):
	link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/period.ods"
	if not os.path.exists(filename):
	    dbdir=os.getenv("BEAMDBPATH")
	    if dbdir==None:
		print "please define BEAMDBPATH in your env"
		return False
	    self.pydb=os.path.join(dbdir,"pyDB")
	    if not os.path.exists(self.pydb):os.makedirs(self.pydb)
	    filename=os.path.join(self.pydb,filename)
	    if not os.path.exists(filename):
		urllib.urlretrieve(link,filename)
	odsread.__init__(self,filename)
	self.parse("period")
	self.init=False
	
    def __whicharm(self,run):
        if run<20000:return "left"
        elif run<40000:return "right"
        else:return "third"
        
    def __s2i(self,s):
	a=re.search('\d+',s)
	if a==None:return 0
	else:return int(a.group(0))
	
    def __samerun(self,run):
	if not self.init:return False
	if self.run != run:
	    self.init=False
	    return False
	return True
    
    def __findrun(self,run):
	arm=self.__whicharm(run)
	for col in range(len(self.values[0])):
	    #print self.values[0][col],"angle" in self.values[0][col]
	    if arm in self.values[0][col]:self.runcol=col
	    elif "orbit" in self.values[0][col]:self.orbitcol=col
	    elif "field" in self.values[0][col]:
		if "angle" in self.values[0][col]:self.anglecol=col
		else:self.fieldcol=col
	    elif "energy" in self.values[0][col]:self.energycol=col
	    elif "survey" in self.values[0][col]:
		if "bpma" in self.values[0][col]:self.bpmasurveycol=col
		elif "bpmb" in self.values[0][col]:self.bpmbsurveycol=col
		elif "harp04" in self.values[0][col]:self.harp04surveycol=col
		elif "harp05" in self.values[0][col]:self.harp05surveycol=col
	for row in range(1,len(self.values)):
	    self.runrange=map(lambda a:self.__s2i(a),self.values[row][self.runcol].split("~"))
	    if self.runrange[1]==0:self.runrange[1]=100000
	    if run in range(self.runrange[0],self.runrange[1]+1):
		self.runrow=row
		break
	self.init=True
	self.run=run
	
    def __surveydatasplit(self,surveystring):
	pos,angle=tuple(surveystring.split(";"))
	pos,angle=pos.split(","),angle.split(",")
	pos,angle=map(lambda a:float(a),pos),map(lambda a:float(a),angle)
	return [pos,angle]
    #get run range for this orbit
    def runrange(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return self.runrange
    #find orbit for run
    def orbit(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return int(self.values[self.runrow][self.orbitcol])
    #find beam energy for run 
    def energy(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return float(self.values[self.runrow][self.energycol])
    #find target field for run
    def field(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return float(self.values[self.runrow][self.fieldcol])
    #find target field angle for run
    def fieldangle(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return int(self.values[self.runrow][self.anglecol])
    #find bpm A survey info for run
    def bpmasurvey(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return self.__surveydatasplit(self.values[self.runrow][self.bpmasurveycol])
    #find bpm B survey info for run
    def bpmbsurvey(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return self.__surveydatasplit(self.values[self.runrow][self.bpmbsurveycol])
    #find harp 04 survey info for run
    def harp04survey(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return self.__surveydatasplit(self.values[self.runrow][self.harp04surveycol])
    #find harp 05 survey info for run
    def harp05survey(self,run):
	if not self.__samerun(run):self.__findrun(run)
	return self.__surveydatasplit(self.values[self.runrow][self.harp05surveycol])
    #find current info for run
    def current(self,run):
	currdb=os.path.join(self.pydb,"runcurr.pkl")
	if not os.path.exists(currdb):
	    link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/runcurr.pkl"
	    try:urlretrieve(link,currdb)
	    except:"sorry can not download database,please check your network connection"
	currinfo=pickle.load(open(currdb,"rb"))
	if currinfo.has_key(run):currepics=currinfo[run]
	else:return False
	return currepics
 
#get harp peak from harppeaks.ods
class harpinfo(odsread):
    def __init__(self,filename="harppeaks.ods"):
	link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/harppeaks.ods"
	if not os.path.exists(filename):
	    dbdir=os.getenv("BEAMDBPATH")
	    if dbdir==None:
		print "please define BEAMDBPATH in your env"
		return False
	    self.pydb=os.path.join(dbdir,"pyDB")
	    if not os.path.exists(self.pydb):os.makedirs(self.pydb)
	    filename=os.path.join(self.pydb,filename)
	    if not os.path.exists(filename):
		urllib.urlretrieve(link,filename)
	odsread.__init__(self,filename)
	self.parse("harp")
	self.__arrange()
    def __arrange(self):
	#split for different periods
	self.rowperiods={}
	lastperiod=""
	for row in range(len(self.values)):
	    if "#" in self.values[row][0] and self.values[row][0][1].isdigit():
		thisperiod=self.values[row][0][1:]
		self.rowperiods[thisperiod]={}
		self.rowperiods[thisperiod]["firstrow"]=row+1
		if len(lastperiod)>0:self.rowperiods[lastperiod]["lastrow"]=row-1
		lastperiod=thisperiod
	self.rowperiods[lastperiod]["lastrow"]=len(self.values)-1
    #find bpm calibration data from ods file, including harp peak,runs for calib
    #keywords is a string list,
    def getdata(self,keywords):
	#find period
	matches={}
	for key in self.rowperiods.keys():
	    matches[key]=0
	    for keyword in keywords:
		if keyword in key:
		    matches[key]+=1
	    if matches[key]<1:del matches[key]
	if len(matches)>1:
	    matches=sorted(matches.items(),key=lambda d:d[1],reverse=True)
	    if matches[0][1]==matches[1][1]:
		print "sorry not enough keyword!!!"
		print "your keywords:",keywords
		print "avail period:",self.rowperiods.keys()
		return False
	    period=matches[0][0]
	elif len(matches)<1:
	    print "sorry can not find right data from your keywords."
	    print "your keywords:",keywords
	    print "avail period:",self.rowperiods.keys()
	    return False
	else:period=matches.keys()[0]
	rowperiod=self.rowperiods[period]
	#get right setting col
	matches={}
	choice=[]
	for col in range(8,len(self.values[rowperiod["firstrow"]-1])):
	    if len(self.values[rowperiod["firstrow"]-1][col])<1:continue
	    choice.append(self.values[rowperiod["firstrow"]-1][col])
	    matches[col]=0
	    for keyword in keywords:
		if keyword in self.values[rowperiod["firstrow"]-1][col]:
		    matches[col]+=1
	    if matches[col]<1:del matches[col]
	if len(matches)>1:
	    matches=sorted(matches.items(),key=lambda d:d[1],reverse=True)
	    if matches[0][1]==matches[1][1]:
		print "sorry not enough keyword!!!"
		print "your keywords:",keywords
		print "avail choice:",choice
		return False
	    runcol=matches[0][0]
	elif len(matches)<1:
	    print "sorry can not find right data from your keywords."
	    print "your keywords:",keywords
	    print "avail choice:",choice
	    return False
	else:runcol=matches.keys()[0]
	#for left arm or right arm
	LEFT=True
	for key in ["Right","right","RIGHT","R","r"]:
	    if key in keywords:
		LEFT=False
		break
	#get runs
	rowperiod["data"]={}
	ndata=0
	for row in range(rowperiod["firstrow"],rowperiod["lastrow"]+1):
	    try:
		if len(self.values[row][runcol])<1:continue
	    except:continue
	    allruns=[int(a) for a in re.split('[\s\\\/,]',self.values[row][runcol]) if a.isdigit()]
	    runs=[]
	    for run in allruns:
		if run<20000 and run>1000 and LEFT:runs.append(run)
		elif run>20000 and not LEFT:runs.append(run)
	    if len(runs)<1:continue
	    if len(self.values[row][1])>0:
		harp04peaks=[float(a) for a in self.values[row][1:4]]
		harp05peaks=[float(a) for a in self.values[row][5:8]]
		rowperiod["data"][ndata]={"run":runs,"harp04":harp04peaks,"harp05":harp05peaks}
		ndata+=1
	    elif "ped" in self.values[row][0]:rowperiod["data"]["pedrun"]=runs
	    elif "curr" in self.values[row][0]:rowperiod["data"]["currrun"]=runs
	rowperiod["data"]["ndata"]=ndata
	return rowperiod["data"]
 
class pkgsetting():
    def __init__(self):
	dbdir=os.path.join(os.getenv("BEAMDBPATH"),"pyDB")
	filepath=os.path.join(dbdir,"beampackage.conf")
	if not os.path.exists(filepath):
	    if not os.path.exists(dbdir):os.makedirs(dbdir)
	    link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/beampackage.conf"
	    urllib.urlretrieve(link,filepath)
	if not os.path.exists(filepath):
	    sys.exit("Error!!beampackage.conf lost,please add it first.")
	datfile=open(filepath,"r")
	for line in datfile:
	    line=line.strip().split("#")[0]
	    if len(line)<1:continue
	    line=line.split("=")
	    if line[0]=="hapunavail":
		period=re.split(",",line[1].strip())
		self.hapunavail=[]
		for p in period:
		    if p.isdigit():self.hapunavail.append(int(p))
		    else:
			ptmp=[int(x) for x in re.split("[ ~-]",p)]
			self.hapunavail+=range(ptmp[0],ptmp[1]+1)
	    elif line[0]=="filter1":self.filter1=float(line[1])
	    elif line[0]=="filter2":self.filter2=float(line[1])
	    elif line[0]=="filter3":self.filter3=float(line[1])
	    elif line[0]=="filter1type":self.filter1type=line[1].strip()
	    elif line[0]=="filter2type":self.filter2type=line[1].strip()
	    elif line[0]=="filter3type":self.filter3type=line[1].strip()
   

if __name__ == '__main__':
    period=runinfo()
    print period.current(22254)
    