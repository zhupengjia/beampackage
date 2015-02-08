#!/usr/bin/env python
import time,linecache,sys,re,os
def times2r(t):
    os.environ['TZ'] = 'US/Eastern'
    time.tzset()
    if t.isdigit():
      return int(t)
    else:
      return int(time.mktime(time.strptime(t,"%Y-%m-%d %H:%M:%S")))
def timer2s(t):
    os.environ['TZ'] = 'US/Eastern'
    time.tzset()
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t))

class epicst2v:
    def __init__(self, filename):
      self.epicsfilen=filename
      try:
#         tt1=time.clock()
          self.epicsfile=open(self.epicsfilen,"rb")
          self.totline = 0
          while True:
              fbuffer = self.epicsfile.read(8192*1024)
              if not fbuffer:
                  break
              self.totline += fbuffer.count('\n')
#          self.totline=len(self.epicsfile.readlines())
#          print self.totline
#          tt2=time.clock()
      except Exception as e:
          print e
          raise Exception("file not existed or not available!!")
      tmptime=linecache.getline(self.epicsfilen,1)[0:10]
      if tmptime.isdigit():
          self.valsplit=10
          self.timestart=int(tmptime)
      else:
          self.valsplit=19
          self.timestart=times2r(linecache.getline(self.epicsfilen,1)[:self.valsplit])
 #       tt3=time.clock()
          self.timeend=times2r(linecache.getline(self.epicsfilen,self.totline)[:self.valsplit])
 #       tt4=time.clock()
 #       print "total time:",tt2-tt1,tt3-tt2,tt4-tt3

    def __del__(self):
        self.epicsfile.close()
        linecache.clearcache()

    def rawtime2value(self,t):
        t=int(t)
#      print t,self.timestart,self.timeend
        if t>self.timeend or t<self.timestart:
            print "time out of range, please use range %s to %s"%(self.timestart,self.timeend)
            #sys.exit()
            return None
        found=0
        line1=1
        line2=int(self.totline/2)
        line3=self.totline
        time1=times2r(linecache.getline(self.epicsfilen,line2)[:self.valsplit])
        while found==0:
            if t>time1:
                line1=line2
            elif t<time1:
                line3=line2
            else:
                found=1
                continue
            line2=int((line3+line1)/2)
 #           print line1,line2,line3
            if line2-line1<1:
                found=1
                continue
            elif line3-line2<2:
                found=1
                continue
            time1=times2r(linecache.getline(self.epicsfilen,line2)[:self.valsplit])
        value=linecache.getline(self.epicsfilen,line2)[self.valsplit+1:].replace("\n", "")
        if value.isdigit():
            return float(value)
        else:
            return value
    def time2value(self,t):
 #       print self.times2r(t)
        return self.rawtime2value(times2r(t))

if __name__ == '__main__':
    if len(sys.argv)<3:
        print "usage: ./epicst2v epics_path time"
        print "format for time can be raw or \"%Y-%m-%d %H:%M:%S\""
        sys.exit()
    ereadname=sys.argv[1]
    etime=sys.argv[2]
    eread=epicst2v(ereadname)
    if etime.isdigit():
        print eread.rawtime2value(int(etime))
    else:
        print eread.time2value(etime)

