#!/usr/bin/env python
import time,os,sys
from urllib2 import urlopen
import socket
socket.setdefaulttimeout(30)
#def timer2s(t):
#        return time.strftime("%Y-%m-%d %H:%M:%S %Z", time.localtime(t))

def grubwikifile(link):
    import re
    localtz=os.getenv("TZ")
    filepattern="<td><a\s+href=\"(?P<link>\S+)\">(?P<filename>[\w\s]+.\w+)</a></td><td\s+align=\"right\">(?P<time>\d+-\w+-\d+ \d+:\d+)\s+</td><td\s+align=\"right\">(?P<size>[.\s\w]+)</td>"
    filere=re.compile(filepattern)
    try:web=urlopen(link,timeout=20)
    except:
	print "can not check update,please check your network"
	return False
    content=web.read()
    fileiter=filere.finditer(content)
    remotefiles={}
    os.environ['TZ']='EST+05EDT,M4.1.0,M10.5.0'
    time.tzset()
  #  timeoff=time.altzone if time.localtime().tm_isdst==1 else time.timezone
    for matcher in fileiter:
	filelink=matcher.group("link")
	filename=matcher.group("filename")
	filetime=matcher.group("time")
	filesize=matcher.group("size")
	filetime=int(time.mktime(time.strptime(filetime,"%d-%b-%Y %H:%M")))
	remotefiles[filename]=(filetime,os.path.join(link,filelink))
#	print "remote",filename,filetime,timer2s(filetime),timeoff
    if localtz:os.environ['TZ']=localtz
    return remotefiles

def syncfileswithwiki(localdir,link,comment="",force=False):
    from urllib import urlretrieve
    from select import select
    from glob import glob
    #list localfiles
    listfiles=map(lambda d:os.path.split(d)[1],glob(os.path.join(localdir,'*')))
    localfiles={}
 #   timeoff=time.altzone if time.localtime().tm_isdst==1 else time.timezone
    for listfile in listfiles:
	try:
	    filetime=int(os.path.getmtime(os.path.join(localdir,listfile)))
	    localfiles[listfile]=filetime
	except:continue
#	print "local",listfile,filetime,timer2s(filetime),timeoff
    #list remotefiles
    remotefiles=grubwikifile(link)
    if not isinstance(remotefiles,dict):return False
   # return
    #compare remote and local
    auto=os.getenv("BEAMUPDATE")
    askforupdate=False
    for k in remotefiles.keys():
	if not localfiles.has_key(k) or remotefiles[k][0]-localfiles[k]>7200:
	    changelog=getchangelog()
	    if not askforupdate:print "%s found newer version,changelog:\n%s\n"%(comment,changelog)
	    if auto=="ask" and not askforupdate and not force:
		print "do you want to update it?(Y/y,will default update after 30 secs)"
		i,o,e =select.select([sys.stdin],[],[],30)
		if i:ifupdate=sys.stdin.readline().strip()
		else:ifupdate="y"
		if not ifupdate in ["Y","y"]:
		    print "skip update"
		    return False
		else:"will update the %s"%comment
	    askforupdate=True
	    if auto!="0" or force:
		try:
		    urlretrieve(remotefiles[k][1],os.path.join(localdir,k))
		    print "%s %s updated,time old:%s new:%s,it is better to rerun script again"%(comment,k,localfiles[k],remotefiles[k][0])
		except Exception as e:
		    print "failed to update %s,please check your network connection"%comment
		    print e
		    return False
    return True

def updateconst(force=False):
    dbdir=os.getenv("BEAMDBPATH")
    if dbdir==None:
	print "please define BEAMDBPATH in your env"
	return False
    pydb=os.path.join(dbdir,"pyDB")
    if not os.path.exists(pydb):os.makedirs(pydb)
    dbwebpath="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/"
    return syncfileswithwiki(pydb,dbwebpath,"bpm const database",force)

def updatecode(force=False):
    import beampackage
    codedir=os.path.split(beampackage.__file__)[0]
    codewebpath="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/beampackage/"
    return syncfileswithwiki(codedir,codewebpath,"beam package code",force)

def getchangelog():
    logpath="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/beampackage/changelog"
    try:web=urlopen(logpath,timeout=20)
    except:return "failed to get changelog,network connection failed"
    content=web.read()
    return content.split("----------")[0]

def updatepackage():
    from threading import Thread
    p1=Thread(target=updateconst)
    p2=Thread(target=updatecode)
    p1.start()
    p2.start()


