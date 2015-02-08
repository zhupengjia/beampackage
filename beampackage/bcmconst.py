#!/usr/bin/env python
import re,urllib,os,numpy
from odsread import odsread

#keyword: sis3800,sis3801,happex,fast,slow,up,down
class bcmconst(odsread):
    def __init__(self,filename="bcm calibration.ods"):
        link="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/bcm%20calibration.ods"
        if not os.path.exists(filename):
          dbdir=os.getenv("BEAMDBPATH")
          if dbdir==None:
            print "please define BEAMDBPATH in your env"
            return False
          pydb=os.path.join(dbdir,"pyDB")
          if not os.path.exists(pydb):os.makedirs(pydb)
          filename=os.path.join(pydb,filename)
          if not os.path.exists(filename):
            urllib.urlretrieve(link,filename)
        odsread.__init__(self,filename)
        self.parse("bcm")

    def __whicharm(self,run):
        if run<20000:return ("left",[22,45])
        elif run<40000:return ("right",[47,70])
        else:return ("third",[72,83])

    #device=sis3800,sis3801,happex  updown=up,down
    def findcol(self,run,device,updown):
        arm=self.__whicharm(run)[0]
        mat=re.compile("%s\s*%s\s*%s"%(arm,device,updown),re.I)
        repstart=re.compile("start",re.I)
        repend=re.compile("end",re.I)
        cols=[]
        for row in range(4,20):
            if mat.search(self.values[row][0]):
                col=0
                for cell in self.values[row]:
                    runrange=map(lambda a:repstart.sub("0",a),cell.split("~"))
                    runrange=map(lambda a:repend.sub("100000",a),runrange)
                    try:
                        runrange=[int(runrange[0]),int(runrange[1])]
                        if run in range(runrange[0],runrange[1]+1):
                            cols.append(col)
                    except ValueError:
                        pass
                    col+=1
                break
        if len(cols)>0:return cols
        else: return False

    def getconst(self,run,device,updown,clock="fast",quiet=False):
        cols=self.findcol(run,device,updown)
        if not cols:
          if not quiet:
            print "can not find constant for run %i,device %s %s, check run number, maybe it is broken during that time"%(run,device,updown)
          return False
        if device=="happex":clk=""
        else:
            if clock=="fast":clk="fstclk"
            elif clock=="slow":clk="slowclk"
            else:
                print "please confirm to use slow or fast clock"
                return False
        mat=re.compile("%s.*%sslope.*%s"%(device,updown,clk),re.I)
        armrows=self.__whicharm(run)[1]
        consts=[]
        for row in range(armrows[0],armrows[1]):
            if mat.search(self.values[row][0]):
                for col in cols:
                    try:
                        slope=float(self.values[row][col])
                        ped=float(self.values[row+1][col])
                    except ValueError:
                        continue
                    consts.append([slope,ped])
                break
        if len(consts)==0:
          if not quiet:
            print "bcm constant not exists for run %i,%s %sclock %sstream"%(run,device,clock,updown)
          return False
        slope=0
        ped=0
        for i in range(len(consts)):
            slope+=float(consts[i][0])
            ped+=float(consts[i][1])
        slope/=len(consts)
        ped/=len(consts)
        return (slope,ped)

#const is (slope,ped),unit is uA
def getcurr(rate,const,device,clock="slow"):
    if not const:
      return numpy.nan
    if device=="happex":
      return const[0]*(rate-const[1])*875/1041.65
    elif device=="sis3801":
      return const[0]*(rate-const[1]*103700*971.65e-6)/1041.65e-6
    else:
      if clock=="fast":clockrate=103700
      else:clockrate=1024
      return const[0]*(rate-const[1]*clockrate)

#count is total count, clockcnt is total clock count(happex is total entry),unit is uC
def getcharge(count,clockcnt,const,device,clock="fast"):
    if not const:
      return numpy.nan
    if device=="happex":
      return const[0]*875e-6*(count-const[1]*clockcnt)
    else:
      return const[0]*(count-const[1]*clockcnt)

#get raw value from current
def curr2raw(curr,const,device,clock="fast"):
    if not const:
      return numpy.nan
    if device=="happex":
      return curr*1041.65/875/const[0]+const[1]
    elif device=="sis3801":
      return curr*1041.65e-6/const[0]+const[1]*103700*971.65e-6
    else:
      if clock=="fast":clockrate=103700
      else:clockrate=1024
      return curr/const[0]+const[1]*clockrate

if __name__ == '__main__':
    bcmconstfile="/home/pzhu/work/run record/bcm calibration.ods"
    const=bcmconst(bcmconstfile)

    print const.getconst(42017,"sis3800","down","fast")
    del const
    #use getconst to get calibration constant
    #first parameter is run number
    #second is which device do you use, it should be "sis3800","sis3801" or "happex"
    #third is upstream or downstream, it should be "up" or "down"
    #fourth is fast clock or slow clock, for happex it is useless. it should be "fast" or "slow"
    #you do not need to verify the bcm calibration file path, this code will automatic download it to your computer

