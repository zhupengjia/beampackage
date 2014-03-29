#!/usr/bin/env python
import numpy,math,sys
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

class harppos:
  def __init__(self,harpsurvey):
    self.wirepos_04=[30.449,59.2646,69.5387]
    self.wirepos_05=[30.26544,42.88606,58.96446]
    self.wiretype_04=[2,1,0] #order of wire for | \ /, inside the beam is first
    self.wiretype_05=[1,2,0]
    self.wirerot_04=map(lambda a,b:math.radians(a+b),[0,0,0],[harpsurvey[0][1][0],-harpsurvey[0][1][1],harpsurvey[0][1][2]])
    self.wirerot_05=map(lambda a,b:math.radians(a+b),[0,0,45],[harpsurvey[1][1][0],-harpsurvey[1][1][1],harpsurvey[1][1][2]])
#    self.harppos_04=[-harpsurvey[0][0],harpsurvey[0][1],-876.933-harpsurvey[0][2]]
#    self.harppos_05=[-harpsurvey[1][0],harpsurvey[1][1],-876.933-harpsurvey[1][2]]
    self.harppos_04=[-harpsurvey[0][0][0],harpsurvey[0][0][1],-harpsurvey[0][0][2]]
    self.harppos_05=[-harpsurvey[1][0][0],harpsurvey[1][0][1],-harpsurvey[1][0][2]]

  def __getharppos(self,peak,wirepos,wiretype,wirerot):
    pos_raw=[wirepos[wiretype[0]]-peak[wiretype[0]],
	  ((wirepos[wiretype[1]]-wirepos[wiretype[2]])-(peak[wiretype[1]]-peak[wiretype[2]]))/2.,
	  0]
#    print pos_raw
#    pos_raw=[wirepos[wiretype[0]]-peak[wiretype[0]],
#	    wirepos[wiretype[1]]-peak[wiretype[1]]-(wirepos[wiretype[0]]-peak[wiretype[0]]),
#	    0]
#    print pos_raw
#    pos_raw=[wirepos[wiretype[0]]-peak[wiretype[0]],
#	    -wirepos[wiretype[2]]+peak[wiretype[2]]+(wirepos[wiretype[0]]-peak[wiretype[0]]),
#	    0]
#    print pos_raw,"\n"
#    print "should be 0:",\
#	(wirepos[wiretype[1]]+wirepos[wiretype[2]]-2*wirepos[wiretype[0]])\
#	-(peak[wiretype[1]]+peak[wiretype[2]]-2*peak[wiretype[0]])
    return rotate(pos_raw,wirerot,1)

  def getpos_04(self,peak):
    pos=self.__getharppos(peak,self.wirepos_04,self.wiretype_04,self.wirerot_04)
    return map(lambda a,b:a+b,pos,self.harppos_04)

  def getpos_05(self,peak):
    pos=self.__getharppos(peak,self.wirepos_05,self.wiretype_05,self.wirerot_05)
    #print pos
    return map(lambda a,b:a+b,pos,self.harppos_05)

class bpmpos:
  def __init__(self,bpmsurvey):
#    self.bpmpos_survey=[-bpmsurvey[0][0],bpmsurvey[0][1],-876.933-bpmsurvey[0][2]]
    self.bpmpos_survey=[-bpmsurvey[0][0],bpmsurvey[0][1],-bpmsurvey[0][2]]
    self.bpmrot=map(lambda a,b:math.radians(a+b),[0,0,45],[bpmsurvey[1][0],-bpmsurvey[1][1],bpmsurvey[1][2]])
    #bpm wire is not in the center of bpm
    #wire pos=(5/2-(1.86+(8-7.03)/2)-1.525/2)*2.54
    self.bpmwire_pos=[0,0,15.4305]
    self.bpm_hardpos=map(lambda a,b:a+b,self.bpmpos_survey,rotate(self.bpmwire_pos,self.bpmrot,1))

  def getpos_hall(self,pos1,pos2):
    x=(self.bpm_hardpos[2]-pos1[2])/(pos2[2]-pos1[2])*(pos2[0]-pos1[0])+pos1[0]
    y=(self.bpm_hardpos[2]-pos1[2])/(pos2[2]-pos1[2])*(pos2[1]-pos1[1])+pos1[1]
    return [x,y,self.bpm_hardpos[2]]

  def getpos_bpm(self,pos1,pos2): #bpm local coordinate
    pos_hall=self.getpos_hall(pos1,pos2)
    pos=map(lambda a,b:a-b,pos_hall,self.bpm_hardpos)
    return rotate(pos,self.bpmrot)

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

if __name__ == '__main__':
    #survey data: 1H04 pos,1H05 pos,1H04 rot,1H05 rot(yaw,pitch,roll in optim coordinate)
    harpsurvey=[[-0.1,0.1,8025.2],[0.4,-0.3,815.6],[0.2074,0.1556,0],[0.0316,0.035,0.1157]]
    #survey data, pos and rot
    bpmasurvey=[[0.4,-0.5,958.1],[0.1036,0.1518,0.0456]]
    bpmbsurvey=[[0.4,-0.3,692.5],[0.1101,-0.1375,-12.3058]]
    harp=harppos(harpsurvey)
    bpma=bpmpos(bpmasurvey)
    bpmb=bpmpos(bpmbsurvey)
    peak_04=[[26.9583,60.3917,68.325],
	[30.6417,56.7083,68.325],
	[21.8583,57.8417,64.3583],
	[24.125,65.775,69.4583],
	[26.1083,66.9083,71.1583],
	[27.8083,65.2083,70.875],
	[20.725,58.4083,64.075]]
    peak_05=[[28.9417,41.6917,57.5583],
	[28.9417,37.1583,48.775],
	[17.8917,38.2917,60.9583],
	[31.775,51.0417,72.575],
	[44.525,51.325,61.2417],
	[40.275,52.175,67.1917],
	[22.9917,44.525,68.6083]]
#    peak_04=[[26.454,62.1248,68.9421],
#             [29.0105,64.3954,71.401],
#             [29.0317,59.1047,68.8396],
#             [26.3632,61.3382,68.4757]]
#    peak_05=[[29.9253,42.9757,59.0408],
#             [37.4853,46.4791,58.7262],
#             [30.5925,40.3133,53.1247],
#             [24.4169,40.4344,59.127]]

#    print rotate([1,1,0],[0,0,45],1)
    for i in range(len(peak_04)):
#    for i in range(1):
	pos_harp04=harp.getpos_04(peak_04[i])
	pos_harp05=harp.getpos_05(peak_05[i])
#	print pos_harp05
	print bpma.getpos_hall(pos_harp04,pos_harp05),bpmb.getpos_hall(pos_harp04,pos_harp05)
	print bpma.getpos_bpm(pos_harp04,pos_harp05),bpmb.getpos_bpm(pos_harp04,pos_harp05)
   # print
