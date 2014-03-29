#!/usr/bin/env python
import numpy,math,sys
from runinfo import runinfo
def rotate(pos,angle,anti=0): #angle is yaw(y),pitch(x),roll(z) in optim coordinate
    #anti=0:from original coordinate(hall) to rotated coordinate(local)
    #anti=1:from rotated coordinate(local) to original coordinate(hall)
    #rotation should be Rz(roll)Rx(pitch)Ry(yaw),each angle should anticlockwise for position when axis direct to your eye
    a=numpy.matrix(pos).T
    (yaw,pitch,roll)=(angle[0],angle[1],angle[2])
    (sy,cy)=(math.sin(yaw),math.cos(yaw))
    (sp,cp)=(math.sin(pitch),math.cos(pitch))
    (sr,cr)=(math.sin(roll),math.cos(roll))
    Myaw=numpy.matrix([[cy,0,-sy],[0,1,0],[sy,0,cy]])
    Mpitch=numpy.matrix([[1,0,0],[0,cp,sp],[0,-sp,cp]])
    Mroll=numpy.matrix([[cr,sr,0],[-sr,cr,0],[0,0,1]])
    M=Mroll*Mpitch*Myaw
    if anti==0:return numpy.array((M*a).T)[0].tolist()
    else:return numpy.array((M.I*a).T)[0].tolist()	

class harppos(runinfo):
  def __init__(self,run):
      #run is used for getting survey data
      runinfo.__init__(self)
      harp04survey=self.harp04survey(run)
      harp05survey=self.harp05survey(run)
      harp04wire=self.harp04wire(run)
      harp05wire=self.harp05wire(run)
      self.wirepos_04=harp04wire[0]
      self.wirepos_05=harp05wire[0]
      self.wiretype_04=self.harp04type(run) #order of wire for | \ /, inside the beam is first
      self.wiretype_05=self.harp05type(run)
      self.wirerot_04=map(lambda a,b:math.radians(a+b),harp04wire[1],[harp04survey[1][0],-harp04survey[1][1],harp04survey[1][2]])
      self.wirerot_05=map(lambda a,b:math.radians(a+b),harp05wire[1],[harp05survey[1][0],-harp05survey[1][1],harp05survey[1][2]])
      #    self.harppos_04=[-harp04survey[0],harp04survey[1],-876.933-harp04survey[2]]
      #    self.harppos_05=[-harp05survey[0],harp05survey[1],-876.933-harp05survey[2]]
      self.harppos_04=[-harp04survey[0][0],harp04survey[0][1],-harp04survey[0][2]]
      self.harppos_05=[-harp05survey[0][0],harp05survey[0][1],-harp05survey[0][2]]

  def __getharppos(self,peak,wirepos,wiretype,wirerot):
      pos_raw=[wirepos[wiretype[0]]-peak[wiretype[0]],\
	       ((wirepos[wiretype[1]]-wirepos[wiretype[2]])-(peak[wiretype[1]]-peak[wiretype[2]]))/2.,\
	       0]
      return rotate(pos_raw,wirerot,1)

  def getpos_04(self,peak):
      return map(lambda a,b:a+b,self.getpos_04_local(peak),self.harppos_04)

  def getpos_05(self,peak):
      return map(lambda a,b:a+b,self.getpos_05_local(peak),self.harppos_05)
    
  def getpos_04_local(self,peak):
      return self.__getharppos(peak,self.wirepos_04,self.wiretype_04,self.wirerot_04)

  def getpos_05_local(self,peak):
      return self.__getharppos(peak,self.wirepos_05,self.wiretype_05,self.wirerot_05)

class bpmpos(runinfo):
  def __init__(self,run,ab):
      #ab is "a" for bpm A or "b" for bpm B,run is used for getting survey data
      runinfo.__init__(self)
      bpmsurvey=self.bpmasurvey(run) if ab=="a" else self.bpmbsurvey(run)
      bpmwirepos=self.bpmwire(run)
      #    self.bpmpos_survey=[-bpmsurvey[0][0],bpmsurvey[0][1],-876.933-bpmsurvey[0][2]]
      self.bpmpos_survey=[-bpmsurvey[0][0],bpmsurvey[0][1],-bpmsurvey[0][2]]
      self.bpmrot=map(lambda a,b:math.radians(a+b),bpmwirepos[1],[bpmsurvey[1][0],-bpmsurvey[1][1],bpmsurvey[1][2]])
      #bpm wire is not in the center of bpm
      #wire pos=(5/2-(1.86+(8-7.03)/2)-1.525/2)*2.54
      self.bpmwire_pos=bpmwirepos[0]
      self.bpm_hardpos=map(lambda a,b:a+b,self.bpmpos_survey,rotate(self.bpmwire_pos,self.bpmrot,1))

  def getpos_hall(self,pos1,pos2):
      x=(self.bpm_hardpos[2]-pos1[2])/(pos2[2]-pos1[2])*(pos2[0]-pos1[0])+pos1[0]
      y=(self.bpm_hardpos[2]-pos1[2])/(pos2[2]-pos1[2])*(pos2[1]-pos1[1])+pos1[1]
      return [x,y,self.bpm_hardpos[2]]

  def getpos_bpm(self,pos1,pos2): #bpm local coordinate
      pos_hall=self.getpos_hall(pos1,pos2)
      pos=map(lambda a,b:a-b,pos_hall,self.bpm_hardpos)
      return rotate(pos,self.bpmrot)
    
  def getpos_rot(self,pos1,pos2): #bpm rot coordinate
      pos_hall=self.getpos_hall(pos1,pos2)
      return map(lambda a,b:a-b,pos_hall,self.bpm_hardpos)

  def posbpm2hall(self,pos_bpm):
      pos_bpm=[pos_bpm[0],pos_bpm[1],0]
      pos=rotate(pos_bpm,self.bpmrot,1)
      return map(lambda a,b:a+b,pos,self.bpm_hardpos)

  def posbpmrotate(self,pos_bpm):
      #convert bpm coordinate to rot coordinate
      pos_bpm=[pos_bpm[0],pos_bpm[1],0]
      pos=rotate(pos_bpm,self.bpmrot,1)
      return pos
  
  def posbpmhall2hall(self,pos_hall):
      #add z to hall coordinate position
      return [pos_hall[0],pos_hall[1],self.bpm_hardpos[2]]
  
  def posbpmhall2rot(self,pos_hall):
      #convert hall coordinate to rot coordinate
      if len(pos_hall)<3:pos_hall.append(self.bpm_hardpos[2])
      return map(lambda a,b:a-b,pos_hall,self.bpm_hardpos)

  def posbpmrot2hall(self,pos_bpmrot):
      #convert rot coordinate to hall coordinate
      return map(lambda a,b:a+b,pos_bpmrot,self.bpm_hardpos)
  
  def posbpmrot2bpm(self,pos_bpmrot):
      return rotate(pos_bpmrot,self.bpmrot)

  def posbpmhall2bpm(self,pos_hall):
      pos_rot=self.posbpmhall2rot(pos_hall)
      return rotate(pos_rot,self.bpmrot)
