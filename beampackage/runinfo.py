#!/usr/bin/env python
import re,urllib,os,glob,sys,numpy,zlib,hashlib,sqlite3,shutil
from odsread import odsread
try:import cPickle as pickle
except:import pickle
try:import sqlite3
except:
    try:from pysqlite2 import dbapi2 as sqlite3
    except:
      raise Exception("Sorry no sqlite3 module found")
      
#find runinfo for run
class runinfo(odsread):
    def __init__(self,periodfile="period.ods",sqlfile="g2p.db",pkgfile="beampackage.conf",\
      currfile="runcurr.pdt",noisefile="pedgain.pdt",happedfile="pedestal_hap.pdt"):
      self.initperiod,self.periodfile=False,periodfile
      self.initsql,self.sqlfile=False,sqlfile
      self.initpkg,self.pkgfile=False,pkgfile
      self.initconst=False
      self.pklavail={} #for pkl database available
      self.pklinfo={} #save pickled var
      self.initpkl={"curr":False,"noise":False,"happed":False}
      self.pklfile={"curr":currfile,"noise":noisefile,"happed":happedfile}
      self.bpmset={} #save bpm setting
      self.fbpmconst={} #save bpm const by file
      self.bpmconst={} #save bpmconst by run
      
    #period file
    def __initperiod(self):
      if self.initperiod:return
      self.periodfile=checkdb(self.periodfile)
      if not self.periodfile:self.periodavail=False
      else:
          odsread.__init__(self,self.periodfile)
          self.parse("period")
          self.init=False
          self.col={}
          self.periodavail=True
      self.initperiod=True
    #sqlite
    def __initsql(self):
      if self.initsql:return
      tmpdb=os.path.join("/tmp",os.path.split(self.sqlfile)[-1])
      if not os.path.exists(tmpdb):
          self.sqlfile=checkdb(self.sqlfile)
          if not self.sqlfile:
              self.sqlavail=False
              self.initsql=True
              return
          shutil.copyfile(self.sqlfile,tmpdb)
      con=sqlite3.connect(tmpdb)
      con.row_factory = sqlite3.Row
      self.cur=con.cursor()
      self.result={}
      self.sqlavail=True
      self.initsql=True
    #init pkl file
    def __initpkl(self,key):
      if self.initpkl[key]:return
      self.pklfile[key]=checkdb(self.pklfile[key])
      if not self.pklfile[key]:self.pklavail[key]=False
      else:
          self.pklinfo[key]=zload(self.pklfile[key])
          self.pklavail[key]=True
      self.initpkl[key]=True
    #init bpm constant read
    def __initconst(self):
      if self.initconst:return
      dbdir=os.getenv("BEAMDBPATH")
      if dbdir==None:
          print "please define BEAMDBPATH in your env"
          return False
      self.pydb=os.path.join(dbdir,"pyDB")
      if not os.path.exists(self.pydb):os.makedirs(self.pydb)
      self.initconst=True
    #pkgsetting
    def __initpkg(self):
      if self.initpkg:return
      self.pkgfile=checkdb(self.pkgfile)
      self.tgtz=[]
      if not self.pkgfile:self.pkgavail=False
      else:
          datfile=open(self.pkgfile,"r")
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
            elif line[0]=="autogain":
                period=re.split(",",line[1].strip())
                self.autogain=[]
                for p in period:
                  if p.isdigit():self.hapunavail.append(int(p))
                  else:
                      ptmp=[int(x) for x in re.split("[ ~-]",p)]
                      self.autogain+=range(ptmp[0],ptmp[1]+1)
            elif line[0]=="filter1":self.filter1=float(line[1])
            elif line[0]=="filter2":self.filter2=float(line[1])
            elif line[0]=="filter3":self.filter3=float(line[1])
            elif line[0]=="filter1type":self.filter1type=line[1].strip()
            elif line[0]=="filter2type":self.filter2type=line[1].strip()
            elif line[0]=="filter3type":self.filter3type=line[1].strip()
            elif "target" in line[0]:
                #target z
                splitline0=re.split("[\s;]",line[0])
                tgttype=splitline0[1].lower()
                if len(splitline0)>2:
                  tgtrange=[]
                  for rg in re.split(",",splitline0[2]):
                      rg=[int(x) for x in re.split("[-~]",rg)]
                      tgtrange+=range(rg[0],rg[1]+1)
                  tgtrange=numpy.asarray(tgtrange,numpy.int32)
                else:
                  tgtrange=xrange(60000)
                tgtz=[float(x) for x in re.split(",",line[1])]
                self.tgtz.append((tgttype,tgtrange,tgtz))
          self.pkgavail=True
      self.initpkg=True
      
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
    
    #find run period from period.ods,run could be run number or orbit
    def __findrun(self,run):
      self.__initperiod()
      if self.__samerun(run):return
      arm=self.__whicharm(run)
      for col in range(len(self.values[0])):
          if arm in self.values[0][col]:self.col["RunNumber"]=col
          elif "orbit" in self.values[0][col]:self.col["orbit"]=col
          elif "field" in self.values[0][col]:
            if "angle" in self.values[0][col]:self.col["TargetAngle"]=col
            else:self.col["TargetField"]=col
          elif "energy" in self.values[0][col]:self.col["Energy"]=col
          elif "survey" in self.values[0][col]:
            if "bpma" in self.values[0][col]:self.col["bpmasurvey"]=col
            elif "bpmb" in self.values[0][col]:self.col["bpmbsurvey"]=col
            elif "harp04" in self.values[0][col]:self.col["harp04survey"]=col
            elif "harp05" in self.values[0][col]:self.col["harp05survey"]=col
          elif "wire pos" in self.values[0][col]:
            if "bpm" in self.values[0][col]:self.col["bpmwire"]=col
            elif "harp04" in self.values[0][col]:self.col["harp04wire"]=col
            elif "harp05" in self.values[0][col]:self.col["harp05wire"]=col
          elif "type" in self.values[0][col]:
            if "harp04" in self.values[0][col]:self.col["harp04type"]=col
            elif "harp05" in self.values[0][col]:self.col["harp05type"]=col
      if run>100:
          #choose from run
          for row in range(1,len(self.values)):
            self.runrange=map(lambda a:self.__s2i(a),self.values[row][self.col["RunNumber"]].split("~"))
            if self.runrange[1]==0:self.runrange[1]=100000
            if run in range(self.runrange[0],self.runrange[1]+1):
                self.runrow=row
                break
      else:
          #choose from orbit
          for row in range(1,len(self.values)):
            orbit=int(self.values[row][self.col["orbit"]])
            if run==orbit:
                self.runrow=row
                self.runrange=map(lambda a:self.__s2i(a),self.values[row][self.col["RunNumber"]].split("~"))
                if self.runrange[1]==0:self.runrange[1]=100000
                break
      self.init=True
      self.run=run
      
    def __surveydatasplit(self,surveystring):
      try:
          pos,angle=tuple(surveystring.split(";"))
          pos,angle=pos.split(","),angle.split(",")
          pos,angle=[float(a) for a in pos],[float(a) for a in angle]
      except Exception as err:
          raise Exception("Error!please check the database file %s in pyDB!Something wrong with it!\n Error info:%s"%(self.periodfile,err))
      return [pos,angle]
    #get info from sqlite
    def getsqlinfo(self,run,key=False):
      self.__initsql()
      if not self.sqlavail:return False
      arm=self.__whicharm(run)
      if not self.result.has_key(run):
          self.cur.execute("select * from %s where RunNumber=%i"%(arm,run))
          result=self.cur.fetchone()
          if result==None:return None
          self.result[run]=result
      if key:
          if not key in self.result[run].keys():return None
          try:return float(self.result[run][key])
          except:return self.result[run][key]
      return dict(self.result[run])
    #get info automatically from sqlite or period.ods
    def autogetinfo(self,run,key):
      o=self.getsqlinfo(run,key)
      if o==None:
          self.__findrun(run)
          o=float(self.values[self.runrow][self.col[key]])
      return o
          
    #get run range for this orbit
    def getrunrange(self,run):
      self.__initperiod()
      if not self.__samerun(run):self.__findrun(run)
      return self.runrange
    #find orbit for run
    def orbit(self,run):
      return int(self.autogetinfo(run,"orbit"))
    #find beam energy for run 
    def energy(self,run):
      return self.autogetinfo(run,"Energy")
    #find target field for run
    def field(self,run):
      return self.autogetinfo(run,"TargetField")
    #find target field angle for run
    def fieldangle(self,run):
      return int(self.autogetinfo(run,"TargetAngle"))
    #find bpm A survey info for run
    def bpmasurvey(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["bpmasurvey"]])
    #find bpm B survey info for run
    def bpmbsurvey(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["bpmbsurvey"]])
    #find harp 04 survey info for run
    def harp04survey(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["harp04survey"]])
    #find harp 05 survey info for run
    def harp05survey(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["harp05survey"]])
    #find harp 04 wire pos
    def harp04wire(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["harp04wire"]])
    #find harp 05 wire pos
    def harp05wire(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["harp05wire"]])
    #find bpm wire pos
    def bpmwire(self,run):
      self.__findrun(run)
      return self.__surveydatasplit(self.values[self.runrow][self.col["bpmwire"]])
    #find harp04 type
    def harp04type(self,run):
      self.__findrun(run)
      return [int(a) for a in self.values[self.runrow][self.col["harp04type"]].split(",")]
    #find harp05 type
    def harp05type(self,run):
      self.__findrun(run)
      return [int(a) for a in self.values[self.runrow][self.col["harp05type"]].split(",")]
    #find current info for run
    def current(self,run):
      c=self.getsqlinfo(run,"Current")
      if c==None:
          self.__initpkl("curr")
          if self.pklavail["curr"] and self.pklinfo["curr"].has_key(run):
            c=self.pklinfo["curr"][run]
      return c
    #get target type
    def targettype(self,run):
      try:return self.getsqlinfo(run,"TargetType").lower()
      except:return None
    #get target z for target type
    def gettgtz(self,run):
      self.__initpkg()
      tgttype=self.targettype(run)
      tgtz=[0]
      for z in self.tgtz:
          if tgttype==z[0] and run in z[1]:
            tgtz=z[2]
            break
      return tgtz
    #check if happex avail for run, from pkgfile
    def ifhapavail(self,run):
      self.__initpkg()
      if not self.pkgavail:sys.exit("please make sure %s in your database dir"%self.pkgfile)
      if run in self.hapunavail:return False
      else:return True
    #check if run in autogain mode, from pkgfile
    def ifautogain(self,run):
      self.__initpkg()
      if not self.pkgavail:sys.exit("please make sure %s in your database dir"%self.pkgfile)
      if run in self.autogain:return True
      else:return False
      
    #get pedestal for special gain, data from pedestal runs after experiment,only happex available
    def getpedgain(self,arm,div,g1,g2):
      self.__initpkl("noise")
      if not self.pklavail["noise"]:return None
      if g1>49 or g2>49 or div>4:
          print "ped not found"
          return False
      kg1,kg2=[],[]
      aimkey="%s%i%0.2i%0.2i"%(arm,div,g1,g2)
      for key in sorted(self.pklinfo["noise"].keys()):
          if arm!=key[0]:continue
          if div!=int(key[1]):continue
          kg1.append(int(key[2:4]))
          kg2.append(int(key[4:]))
      def getclosedped(fixed):
          for i in range(30):
            if fixed==1:tmpkey="%s%i%0.2i%0.2i"%(arm,div,g1,g2+i)
            elif fixed==2:tmpkey="%s%i%0.2i%0.2i"%(arm,div,g1+i,g2)
            if self.pklinfo["noise"].has_key(tmpkey):
                if fixed==1:gp=g2+i
                elif fixed==2:gp=g1+i
                pp=self.pklinfo["noise"][tmpkey]
                break
          for i in range(30):
            if fixed==1:tmpkey="%s%i%0.2i%0.2i"%(arm,div,g1,g2-i)
            elif fixed==2:tmpkey="%s%i%0.2i%0.2i"%(arm,div,g1-i,g2)
            if self.pklinfo["noise"].has_key(tmpkey):
                if fixed==1:gm=g2-i
                elif fixed==2:gm=g1-i
                pm=self.pklinfo["noise"][tmpkey]
                break
          aaim=10**((g1+g2)/40.)
          if fixed==1:
            a0=10**((g1+gp)/40.)
            a1=10**((g1+gm)/40.)
          elif fixed==2:
            a0=10**((gp+g2)/40.)
            a1=10**((gm+g2)/40.)
          return pp+(pm-pp)*(aaim-a0)/(a1-a0)
      if g1 in kg1:
          if self.pklinfo["noise"].has_key(aimkey):
            return self.pklinfo["noise"][aimkey]
          return getclosedped(1)
      elif g2 in kg2:
          return getclosedped(2)
      return False
    
    #get pedestal from sqlite
    def pedestal(self,run,fastbus=False,nsuf="bpmped"):
      if self.ifautogain(run):return {"a":False,"b":False}
      suffix="fb" if fastbus else "hap"
      ped=[self.getsqlinfo(run,("%s_%s%i"%(nsuf,suffix,i+1))) for i in range(8)]
      if not any(a==None for a in ped[:4]):peda=ped[:4]
      else:peda=False
      if not any(a==None for a in ped[4:]):pedb=ped[4:]
      else:pedb=False
      return {"a":peda,"b":pedb}
    
    #get pedestal rms from sqlite
    def pedestalrms(self,run,fastbus=False):
      return self.pedestal(run,fastbus,"bpmpedrms")
    
    #get bpm setting like diff,div,gain,...
    def getbpmsetting(self,run):
      if self.bpmset.has_key(run):return self.bpmset[run]
      self.bpmset[run]={}
      for key in ['axdiff','aydiff','bxdiff','bydiff','axpdiv','axmdiv','aypdiv','aymdiv','bxpdiv','bxmdiv','bypdiv',\
          'bymdiv','afilter','bfilter','again','bgain','adiv','bdiv','axp1','axp2','axm1','axm2','ayp1','ayp2','aym1',\
            'aym2','bxp1','bxp2','bxm1','bxm2','byp1','byp2','bym1','bym2']:
          self.bpmset[run][key]=self.getsqlinfo(run,key)
      return self.bpmset[run]
    
    #get the most close pedestal value for run from pedestal.pdt
    def pedestalfrompkl(self,run,fastbus=False,silent=False):
      if self.ifautogain(run):
          if fastbus:
            return {"a":False,"b":False}
      suffix="fb" if fastbus else "hap"
      pedkey=suffix+"ped"
      self.__initpkl(pedkey)
      if not self.pklavail[pedkey]:return {"a":False,"b":False}
      elif not silent:print "using pedestal file %s"%self.pklfile[pedkey]
      #getdiv
      ped,pederr=[0]*8,[None]*8
      arm=self.__whicharm(run)
      arml="L" if arm=="left" else "R"
      gain=self.getbpmsetting(run)
      for i in range(8):
          ab="a" if i<4 else "b"
          xy="y" if i%4<2 else "x"
          pm="p" if i%2<1 else "m"
          abxypm="%s%s%s"%(ab,xy,pm)
          div,g1,g2=gain["%sdiv"%abxypm],gain["%s1"%abxypm],gain["%s2"%abxypm]
          if any(a==None for a in (div,g1,g2)):
            ped[i]=None
            continue
          if self.ifautogain(run) and not fastbus:
            try:
                ped[i]=self.getpedgain(arml,div,g1,g2)[i]
            except:ped[i]=None
            if not silent and ped[i]!=None:print "using %s for run %i chan %i,value %i"%(self.pklfile["noise"],run,i,ped[i])
            continue
          runs=[]
          for key in sorted(self.pklinfo[pedkey].keys()):
            keyarm=self.__whicharm(key)
            if keyarm!=arm:continue
            g=self.getbpmsetting(key)
            gdiv,gg1,gg2=g["%sdiv"%abxypm],g["%s1"%abxypm],g["%s2"%abxypm]
            if div!=gdiv:continue
            if g1!=gg1:continue
            if g2!=gg2:continue
            if not self.pklinfo[pedkey].has_key(key):continue
            if not isinstance(self.pklinfo[pedkey][key]["peaks"][i],float):continue
            runs.append(key)
          if len(runs)<1:
            ped[i]=None
            continue
          runs=numpy.asarray(runs,dtype=numpy.int32)
          argmin=numpy.abs(runs-run).argmin()
          ped[i]=self.pklinfo[pedkey][runs[argmin]]["peaks"][i]
          pederr[i]=self.pklinfo[pedkey][runs[argmin]]["rms"][i]
          if not silent:print "using pedestal run %i for run %i chan %i,value %i"%(runs[argmin],run,i,ped[i])
      if not any(a==None for a in ped[:4]):peda=ped[:4]
      else:peda=False
      if not any(a==None for a in ped[4:]):pedb=ped[4:]
      else:pedb=False
      if not any(a==None for a in pederr[:4]):pedaerr=pederr[:4]
      else:pedaerr=False
      if not any(a==None for a in pederr[4:]):pedberr=pederr[4:]
      else:pedberr=False
      return {"a":peda,"b":pedb,"aerr":pedaerr,"berr":pedberr}
    
    #read const from file,filename is file path,ab="a" or "b"
    def bpmconstreadfromfile(self,filename,ab,dbread=False):
      if self.fbpmconst.has_key(filename):
          return self.fbpmconst[filename]
      datfile=open(filename,"r")
      if not isinstance(dbread,dict):
          dbread={"period":False,"curr":False,"tgtz":False,"ped":False,"offset":False,"const":False,"fitorder":False,"keywords":False,"fval":False}
      importfile=False
      for line in datfile:
          if line[0]=="#":continue
          reline=re.split("[:\n]",line)
          data=re.split("[\s,]",reline[1].strip())
          if "keywords" in reline[0] and not dbread["keywords"]:
            dbread["keywords"]=data
          elif "run period" in reline[0] and not dbread["period"]:
            #runperiod=[int(x) for x in data]
            runperiods=re.split("[\s,]",reline[1].strip())
            runperiod=[]
            for p in runperiods:
                if p.isdigit():runperiod.append(int(p))
                else:
                  ptmp=[int(x) for x in re.split("[ ~-]",p)]
                  runperiod+=range(ptmp[0],ptmp[1]+1)
            dbread["period"]=numpy.asarray(runperiod,dtype=numpy.int32)
          elif "curr" in reline[0] and not dbread["curr"]:
            dbread["curr"]=[float(x) for x in data]
          elif "target z" in reline[0] and not dbread["tgtz"]:
            dbread["tgtz"]=[float(x) for x in data]
          elif "pedestal" in reline[0] and not dbread["ped"]:
            dbread["ped"]=[float(x) for x in data]
          elif "offset" in reline[0] and not dbread["offset"]:
            dbread["offset"]=[float(x) for x in data]
          elif "fitorder" in reline[0] and not dbread["fitorder"]:
            dbread["fitorder"]=[int(x) for x in data]
          elif "fval" in reline[0] and not dbread["fval"]:
            dbread["fval"]=[float(x) for x in data]
          elif "bpm%s"%ab in reline[0] and not dbread["const"]:
            if all(x in reline[0] for x in ["ar","gx","gy"]):rgxy=[float(x) for x in data]
            elif all(x in reline[0] for x in ["x","a","b","c"]):cx=[float(x) for x in data]
            elif all(x in reline[0] for x in ["y","a","b","c"]):cy=[float(x) for x in data]
            elif all(x in reline[0] for x in ["x","err"]):cxerr=[float(x) for x in data]
            elif all(x in reline[0] for x in ["y","err"]):cyerr=[float(x) for x in data]
          if "import" in reline[0]:
            importfile=os.path.join(os.path.split(filename)[0],re.split(" ",reline[0].strip())[1])
            break
      if locals().has_key("cy"):dbread["const"]=[tuple(rgxy),cx,cy,len(cx)]
      if locals().has_key("cyerr"):
          dbread["const"].append(cxerr)
          dbread["const"].append(cyerr)
      if importfile:return self.bpmconstreadfromfile(importfile,ab,dbread)
      else:
          self.fbpmconst[filename]=dbread
          return dbread      
    
    #read bpm constant
    def bpmconstread(self,run,forcefastbus=False,silent=False):
      if self.bpmconst.has_key(run):return self.bpmconst[run]
      self.__initconst()
      if not self.ifhapavail(run) or forcefastbus:fastbus=True
      else:fastbus=False
      runcurr=self.current(run)
      def readab(ab):
          datprefix="bpm%sfb"%ab if fastbus else "bpm%s"%ab
          datfiles=glob.glob(os.path.join(self.pydb,"%s_*.dat"%datprefix))
          calruns=[]
          for datfile in datfiles:
            try:calruns.append(int(re.split("[_.]",os.path.split(datfile)[1])[1]))
            except:continue
          calruns=sortrun(run,calruns)
          coinruns=[]
          for calrun in calruns:
            const=self.bpmconstreadfromfile(os.path.join(self.pydb,"%s_%i.dat"%(datprefix,calrun)),ab)
            if not const.has_key("period"):continue
            if run in const["period"] and any(abs(x-(100 if runcurr>100 else runcurr))<25 for x in const["curr"]):
                coinruns.append([calrun,min(abs(x-(100 if runcurr>100 else runcurr)) for x in const["curr"])])
          if len(coinruns)>0:
            calrun=sorted(coinruns,key=lambda r:r[1])[0][0]
            constfile=os.path.join(self.pydb,"%s_%i.dat"%(datprefix,calrun))
            const=self.bpmconstreadfromfile(constfile,ab)
            if not const["offset"]:const["offset"]=[0]*4
            if not silent:print "using calibration dat file %s for run %i"%(constfile,run)
            const["constfile"]=os.path.split(constfile)[1]
            #del const["period"]
            #del const["curr"]bpmconst
            return const
          else:
            if not silent: "sorry no bpm %s calibration constant available for run %i,please contact pengjia"%(ab,run)
            return False
      aconst=readab("a")
      bconst=readab("b")
      const={"a":aconst,"b":bconst}
      self.bpmconst[run]=const
      return const
    
    #get run number list for special condition
    def getruns(self,conditions=[],arm=False):
      self.__initsql()
      if not arm:arms=["left","right"]
      elif arm in ["left","L","Left","l"]:arms=["left"]
      elif arm in ["right","R","Right","r"]:arms=["right"]
      else:arms=["left","right"]
      runs=[]
      for arm in arms:
          if isinstance(conditions,str):
            conditions=[conditions]
          if len(conditions)<1:
            self.cur.execute("select RunNumber from %s"%arm)
          else:
            cmd=conditions[0]
            if len(conditions)>1:
                for c in conditions[1:]:
                  cmd+=" and %s"%c
            self.cur.execute("select RunNumber from %s where %s"%(arm,cmd))
          runs+=[x[0] for x in self.cur.fetchall()]
      return runs
    
#get harp peak from harppeaks.ods
class harpinfo(odsread):
    def __init__(self,filename="harppeaks.ods"):
      filename=checkdb(filename)
      if not filename:return False
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
            print "avail matches:",matches
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
      lr="left" if LEFT else "right"
      #get runs
      rowperiod["data"]={}
      ndata=0
      for row in range(rowperiod["firstrow"],rowperiod["lastrow"]+1):
          try:
            if len(self.values[row][runcol])<1:continue
          except:continue
          allruns=[int(a) for a in re.split('[\s\\\/\-\,]',self.values[row][runcol]) if a.isdigit()]
          runs=[]
          for run in allruns:
            if run<20000 and run>1000 and LEFT:runs.append(run)
            elif run>20000 and not LEFT:runs.append(run)
          if len(self.values[row][runcol])<4:continue
          if len(self.values[row][4])>0:
            try:harp04peaks=[float(a) for a in self.values[row][1:4]]
            except:harp04peaks=False
            harp05peaks=[float(a) for a in self.values[row][5:8]]
            rowperiod["data"][ndata]={"run":runs,"harp04":harp04peaks,"harp05":harp05peaks}
            ndata+=1
          elif "ped" in self.values[row][0]:rowperiod["data"]["pedrun"]=runs
          elif "curr" in self.values[row][0]:rowperiod["data"]["currrun"]=runs
          elif "avail_%s"%lr in self.values[row][0]:
            tmp=self.values[row][runcol].strip().split(";")
            rowperiod["data"]["availa"]=tmp[0]
            rowperiod["data"]["availb"]=tmp[-1]
      rowperiod["data"]["ndata"]=ndata
      return rowperiod["data"]
 
class calibsetting():
    #bpm calibration setting
    def __init__(self,keyword=[]):
      filepath=checkdb("bpmcalib.conf")
      if not filepath:sys.exit("Error!!bpmcalib.conf lost,please add it first.")
      datfile=open(filepath,"r")
      getfirstkeyword,getkeyword,availdata=False,False,True
      keyword=sorted(keyword)
      self.pedpeaks,self.offset,self.gxy,self.switch,self.datanum=False,False,False,False,False
      self.initvalue,self.steps,self.parlow,self.parhigh=False,False,False,False
      self.scantimes,self.fitorder=False,False
      for line in datfile:
          line=line.strip().split("#")[0]
          if len(line)<1:continue
          line=line.split("=")
          if "pedpeaks" in line[0] and availdata:
            if "peak" in line[1]:self.pedpeak=False
            else:
                try:self.pedpeaks=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
                except:pass
          elif "offset" in line[0] and availdata:
            try:self.offset=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "gxy" in line[0] and availdata:
            if "gxy" in line[1]:self.gxy =False
            else:
                try:self.gxy=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
                except:pass
          elif "switch" in line[0] and availdata:
            try:self.switch=int(line[1],2)
            except:pass
          elif "datanum" in line[0] and availdata:
            try:self.datanum=int(line[1])
            except:pass
          elif "initvalue" in line[0] and availdata:
            try:self.initvalue=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "steps" in line[0] and availdata:
            try:self.steps=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "parlimitslow" in line[0] and availdata:
            try:self.parlow=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "parlimitshigh" in line[0] and availdata:
            try:self.parhigh=[float(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "scantimes" in line[0] and availdata:
            try:
                self.scantimes=int(line[1])
                if self.scantimes<100:self.scantimes=False
            except:pass
          elif "fitorder" in line[0] and availdata:
            try:self.fitorder=[int(x) for x in re.split("[, ]",line[1]) if len(x)>0]
            except:pass
          elif "keywords" in line[0]:
            getfirstkeyword=True
            if getkeyword:break
            keywords=sorted([x for x in re.split("[ ]",line[1])])
            if keyword==keywords:
                getkeyword=True
            availdata=getfirstkeyword and getkeyword
                
class getbpmeventcut(odsread):
    def __init__(self,filename="harppeaks.ods"):
      filename=checkdb(filename)
      if not filename:return False
      odsread.__init__(self,filename)
      self.parse("eventcut")
    def getcut(self,run,fastbus=False):
      for i in range(len(self.values)):
          try:r=int(self.values[i][0])
          except:continue
          if run==r:
            if fastbus:cut=self.values[i][2]
            else:cut=self.values[i][1]
            try:
                cut=[int(c) for c in re.split("[\s,]",cut)]
                if len(cut)<2:return False
                return cut[:2]
            except:
                return False
      return False

def sortrun(frun,runs):
    diffrun=map(lambda r:(abs(frun-r),r),runs)
    runs=sorted(diffrun,key=lambda r:r[0])
    return [run[1] for run in runs if run[0]<15000]
      
def rasterconstreadfromfile(filename,dbread=False):
    datfile=open(filename,"r")
    if not isinstance(dbread,dict):
      dbread={"period":False,"tgtz":False,"clkrate":False,"rebuild":False,"slxslope":False,\
          "slyslope":False,"fstxslope":False,"fstyslope":False}
    importfile=False
    for line in datfile:
      reline=re.split("[:\n]",line)
      data=re.split(" ",reline[1].strip())
      #print data
      if "run period" in reline[0] and not dbread["period"]:
          runperiods=re.split(",",reline[1].strip())
          runperiod=[]
          for p in runperiods:
            if p.isdigit():runperiod.append(int(p))
            else:
                ptmp=[int(x) for x in re.split("[ ~-]",p)]
                runperiod+=range(ptmp[0],ptmp[1]+1)
          dbread["period"]=runperiod
      elif "fastclock rate" in reline[0] and not dbread["clkrate"]:
          dbread["clkrate"]=float(data[0])
      elif "function rebuild" in reline[0] and not dbread["rebuild"]:
          dbread["rebuild"]=[int(x)>0 for x in data]      
      elif "target z" in reline[0] and not dbread["tgtz"]:
          dbread["tgtz"]=[float(x) for x in data]
      elif "slow raster x slope" in reline[0] and not dbread["slxslope"]:
          dbread["slxslope"]=[float(x) for x in data]
      elif "slow raster y slope" in reline[0] and not dbread["slyslope"]:
          dbread["slyslope"]=[float(x) for x in data]
      elif "fast raster x slope" in reline[0] and not dbread["fstxslope"]:
          dbread["fstxslope"]=[float(x) for x in data]
      elif "fast raster y slope" in reline[0] and not dbread["fstyslope"]:
          dbread["fstyslope"]=[float(x) for x in data]
      if "import" in reline[0]:
          importfile=os.path.join(os.path.split(filename)[0],re.split(" ",reline[0].strip())[1])
          break
    if importfile:return constread(importfile,dbread)
    else:return dbread      
      
def rasterconstread1(run,happex=False):
    dbdir=os.getenv("BEAMDBPATH")
    #basename="hapraster" if happex else "raster"
    basename="raster" 
    if dbdir==None:
      print "please define BEAMDBPATH in your env"
      return False
    pydb=os.path.join(dbdir,"pyDB")
    if not os.path.exists(pydb):os.makedirs(pydb)
    datfiles=glob.glob(os.path.join(pydb,"%s_*.dat"%basename))
    calruns=[]
    for datfile in datfiles:
      calruns.append(int(re.split("[_.]",os.path.split(datfile)[1])[-2]))
    calruns=sortrun(run,calruns)
    for calrun in calruns:
      const=rasterconstreadfromfile(os.path.join(pydb,"%s_%i.dat"%(basename,calrun)))
      if run in const["period"]:
          del const["period"]
          const["constfile"]="%s_%i.dat"%(basename,calrun)
          return const
    print "sorry no raster calibration constant available for run %i,please contact pengjia"%run
    return False

def rasterconstread(run,happex=False):
    if not happex:return rasterconstread1(run)
    else:
        const_fb=rasterconstread1(run)
        const_hap=rasterconstread1(run,happex)
        info=runinfo()
        r0=info.getsqlinfo(run,"raster_hapfstratio_x") 
        r1=info.getsqlinfo(run,"raster_hapfstratio_y")
        if r0==None or r1==None:
            orbit=info.orbit(run)
            arm="left" if run<20000 else "right"
            availruns=numpy.asarray(info.getruns(["orbit=%i"%orbit,"raster_hapfstratio_x is not NULL","raster_hapfstratio_y is not NULL"],arm))
            if len(availruns)<1:
                return const_hap
            crun=availruns[numpy.argmin(abs(availruns-run))]
            r0=info.getsqlinfo(crun,"raster_hapfstratio_x") 
            r1=info.getsqlinfo(crun,"raster_hapfstratio_y")
        if not const_fb:return False
        const_hap["slxslope"][0]=const_fb["slxslope"][0]/r0
        const_hap["slyslope"][0]=const_fb["slyslope"][0]/r1
        return const_hap

#get pkl file path
class getpklpath:
    def __init__(self,rootfilepath=os.getenv("REPLAY_OUT_PATH")):
      self.rootfilepath=rootfilepath
      self.pkldir=os.getenv("BEAMPKLPATH")
      if self.pkldir==None:
         self. pkldir=os.path.join(self.rootfilepath,"pkl")
      self.pklbak=os.getenv("BEAMPKLBAK")
      
    def getdir(self):
      return self.pkldir
    
    #writemode:if writemode=1,will ignore pklbak dir,that is only save pkl data to pkldir
    def getpath(self,prefix="",varname="",run=0,writemode=0):
      filename="bpm%s_%s_%i.pkl"%(prefix,varname,run) if len(prefix)>0 else "bpm%s_%i.pkl"%(varname,run)
      path=os.path.join(self.pkldir,filename)
      if writemode:
          if not os.path.exists(self.pkldir):
            try:os.makedirs(self.pkldir)
            except:
                sys.exit("Error!!do you have permission to create directory for %s?"%self.pkldir)
          return path
      elif os.path.exists(path):return path
      else:
          if self.pklbak==None:return path
          else:
            path2=os.path.join(self.pklbak,filename)
            if os.path.exists(path2):return path2
            else:return path

def checkrunavail(path,run,checkorder=["bpm","ring","g2p"]):
    rootfiles=filter(lambda p:"%i"%run in p,os.listdir(path))
    filelist=[]
    for order in checkorder:
      pattern="%s_?[LR]?_%i.root"%(order,run)
      t=[]
      for f in rootfiles:
          if re.match(pattern,f):t.append(f)
      t.sort()
      filelist+=t
    if len(filelist)<1:
      raise Exception("sorry no file found for run %s in %s"%(run,path))
    return os.path.join(path,filelist[0])

#automatically download database file if not exists
def checkdb(filename,link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB"):
    if os.path.exists(filename):return filename
    dbdir=os.getenv("BEAMDBPATH")
    if dbdir==None:
      print "please define BEAMDBPATH in your env"
      return False
    pydb=os.path.join(dbdir,"pyDB")
    if not os.path.exists(pydb):os.makedirs(pydb)
    filename=os.path.join(pydb,filename)
    if not os.path.exists(filename):
      print "sorry no file %s exists, please update your database"%filename
      #try:urllib.urlretrieve(os.path.join(link,filename),filename)
      #except:
         # print "sorry can not download %s,please check your network connection"%filename
         # return False
    return filename
    
#md5
def md5_file(name):
    with open(name, 'rb') as f:
        try:f.seek(- 4096 * 1024, 2)# for large file,only read last 4mb 
        except:pass
        md5=hashlib.md5(f.read()).hexdigest()+hashlib.md5(name).hexdigest()
    #print name,md5
    return md5
    
def md5_files(filearray):
    md5=""
    for f in filearray:md5+=md5_file(f)
    return hashlib.md5(md5).hexdigest()[:6]

#compress pickle file by using zlib and cpickle
def zdump(value,filename):
    with open(filename,"wb",-1) as fpz:
      fpz.write(zlib.compress(pickle.dumps(value,-1),9))

#load compressed pkl file from zdump
def zload(filename):
    with open(filename,"rb") as fpz:
      value=fpz.read()
      try:return pickle.loads(zlib.decompress(value))
      except:return pickle.loads(value)
