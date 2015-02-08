#!/usr/bin/env python
modules=["odsread","harppos","harppeak","bpmfit","bcmconst","bpmcalib","bpminsert",\
    "rasterrecon","signalfilter","epicst2v","pacupdate","runinfo"]
for m in modules:
    try:
      exec("from %s import *"%m)
    except Exception as err:
      print "Error!!",err