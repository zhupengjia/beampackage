#!/usr/bin/env python
modules=["sys","os","threading","time","re"]
for m in modules:
    try:exec("import %s"%m)
    except:continue
try:
    from beampackage import *
    nopac=False
except Exception as e:
    nopac=True
    print e

def usage():
    print "usage: ./bpm.py runnum"
    print "or ./bpm.py rootfile"
    print "rootfile can be absolute path or just a file name"
    print "type \"./bpm.py help\" or \"./bpm.py -h\" for help"
    print "type \"./bpm.py update\" to update everything"
    print "type \"./bpm.py env\" will generate a environment file for your beampackage"
    print "type \"./bpm.py changelog\" will print changelog"

def getfile(path,run):
    allfiles=filter(lambda p:"%i"%run in p,os.listdir(path))
    checkorder=["g2p","bpm","ring"]
    filelist=[]
    for order in checkorder:
	pattern="%s_?[LR]?_%i.root"%(order,run)
	t=[]
	for f in allfiles:
	    if re.match(pattern,f) and f[-5:]==".root":
		t.append(f)
	t.sort()
	filelist+=t
    return map(lambda a:os.path.join(path,a),filelist)

def downloadfromwiki(localdir,link):
    import re
    from urllib2 import urlopen
    from urllib import urlretrieve
    import socket
    socket.setdefaulttimeout(30)

    filepattern="<td><a\s+href=\"(?P<link>\S+)\">(?P<filename>[\w\s]+.\w+)</a></td><td\s+align=\"right\">(?P<time>\d+-\w+-\d+ \d+:\d+)\s+</td><td\s+align=\"right\">(?P<size>[.\s\w]+)</td>"
    filere=re.compile(filepattern)
    try:web=urlopen(link,timeout=20)
    except:
	print "can not connect to g2p wiki,please check your network"
	return False
    content=web.read()
    fileiter=filere.finditer(content)
    for matcher in fileiter:
	filelink=os.path.join(link,matcher.group("link"))
	filename=matcher.group("filename")
	filepath=os.path.join(localdir,filename)
	try:
	    print "downloading %s..."%filename
	    urlretrieve(filelink,filepath)
	except:
	    print "download failed,please check your network"
	    return False
    return True

def downloadpackage():
    print "now will install beam package to your computer."
    path=os.getenv("BEAMPACPATH")
    if path==None:
	print "sorry please define BEAMPACPATH to your beampackage's parent dir"
	return False
    pacpath=path
    path=os.path.join(path,"beampackage")
    if not os.path.exists(path):os.makedirs(path)
    print "will download the package to: %s"%pacpath
    codewebpath="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/beampackage/"
    if not downloadfromwiki(path,codewebpath):return False
    print "successful downloaded the package to %s"%pacpath
    return True

def downloaddb():
    path=os.getenv("BEAMDBPATH")
    if path==None:
	print "sorry please define BEAMDBPATH to your beam database's parent dir"
	return False
    path=os.path.join(path,"pyDB")
    print "now will download beam database to %s"%path
    if not os.path.exists(path):os.makedirs(path)
    dbwebpath="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki/pyDB/"
    if not downloadfromwiki(path,dbwebpath):return False
    print "successful downloaded beam database to %s"%path
    return True

def updatemyself(force=False):
    from urllib import urlretrieve
    mypath=os.path.abspath(__file__)
    myname=os.path.split(mypath)[-1]
    mylink="http://hallaweb.jlab.org/experiment/g2p/collaborators/pzhu/wiki"
    try:remotetime=grubwikifile(mylink)[myname][0]
    except:return False
  #  timeoff=time.altzone if time.localtime().tm_isdst==1 else time.timezone
    localtime=int(os.path.getmtime(mypath))
    auto=os.getenv("BEAMUPDATE")
    if remotetime-localtime>7200:
	print "find new version of %s"%myname
	if auto=="ask" and not force:
	    print "do you want to update it?(Y/y,will default update after 30 secs)"
	    i,o,e =select.select([sys.stdin],[],[],30)
	    if i:ifupdate=sys.stdin.readline().strip()
	    else:ifupdate="y"
	    if not ifupdate in ["Y","y"]:
		print "skip update"
		return False
	    else:"will update the %s"%myname
	if os.getenv("BEAMUPDATE")!="0" or force:
	    try:
		urlretrieve(mylink+"/"+myname,mypath)
		print "file %s updated"%myname
	    except:
		print "failed to update %s"%myname
		return False
    return True

def checkenv():
    if nopac:
	makeenv()
	if sys.version_info<(2,6) or sys.version_info>(3,0):
	    raise "your python version is too low or too high,please use 2.6 or higher,but lower than 3.0. 2.7 is the best."
	    return False
        pacneeds=["sys","os","urllib","urllib2","socket","re","threading","numpy","scipy",\
	    "math","pickle","shutil","glob","array","time","linecache","zipfile","ctypes",\
		"xml.etree.ElementTree","select","ROOT","bisect","gc","multiprocessing"]
        haveallpacs=True
        misspacs=[]
        for pacneed in pacneeds:
            try:exec("import %s"%pacneed)
            except:
                misspacs.append(pacneed)
                haveallpacs=False
        if not haveallpacs:
            print "Sorry you need to install these packages first:",
            for misspac in misspacs:print misspac,
            print ""
            if "ROOT" in misspacs:print "please recompile your root with --enable-python flag"
            return False
        try:
            import sqlite3
        except:
            try:
                from pysqlite2 import dbapi2 as sqlite3
            except:
                print "Sorry no sqlite3 module found"
                return False
	d1=downloadpackage()
	d2=downloaddb()
	if d1 and d2:
	    print "successful downloaded beampackage and database to your computer,please rerun this script"
	return False
    else:
	return True

def getabspath(path):
    if path[0]=="$":
	pathenv=path.split("/")[0]
	realpath=os.getenv("%s"%pathenv[1:])
	if realpath==None:
	    print "sorry can not find %s"%pathenv
	    return False
	path=path.replace(pathenv,realpath)
    else:
	if path[0]=="~":path=os.path.expanduser(path)
	path=os.path.abspath(path)
    return path

def makeenv(forceupdate=False):
    if os.getenv("BEAMPACPATH")!=None and not forceupdate:return
    shell=os.path.split(os.getenv("SHELL"))[-1]
    if shell=="tcsh" or shell=="csh":
	shelltype,shellfile="csh",".cshrc"
	env="#!/bin/csh\n"
	setenv=("setenv"," ")
    else:
	shelltype,shellfile="sh",".bashrc"
	env="#!/bin/bash\n"
	setenv=("export","=")
    def ifenv(shelltype,envname,envin):
	if shelltype=="csh":
	    envout="if ($?%s) then\n"%envname
	    envout+="    setenv %s %s:$%s\nelse\n"%(envname,envin,envname)
	    envout+="    setenv %s %s\nendif\n"%(envname,envin)
	else:
	    envout="if [ -z \"${%s}\" ]; then\n"%envname
	    envout+="    export %s=%s\nelse\n"%(envname,envin)
	    envout+="    export %s=%s:$%s\nfi\n"%(envname,envin,envname)
	return envout
    rootpath,beampacpath,beamdbpath=False,False,False
    if os.getenv("REPLAY_OUT_PATH")==None:
	while not rootpath:rootpath=getabspath(raw_input("Please input rootfile path:\n"))
	env+="#default rootfile path\n"
	env+="%s REPLAY_OUT_PATH%s%s\n"%(setenv[0],setenv[1],rootpath)
    while not beampacpath:
	beampacpath=getabspath(raw_input("Please input the beampackage path where you want to install(parent dir):\n"))
    while not beamdbpath:
	beamdbpath=getabspath(raw_input("Please input beam database path(parent dir,for example $DB_DIR):\n"))
    env+="#beam package path(parent dir)\n"
    env+="%s BEAMPACPATH%s\"%s\"\n"%(setenv[0],setenv[1],beampacpath)
    env+="#database path(parent dir)\n"
    env+="%s BEAMDBPATH%s\"%s\"\n"%(setenv[0],setenv[1],beamdbpath)
    env+="#1 will force update,0 will not update,ask will ask before update\n"
    env+="%s BEAMUPDATE%s1\n"%setenv
    env+="#1 will check update,0 will skip check\n"
    env+="%s BEAMUPDATECHECK%s1\n"%setenv
    pythonpath="${BEAMPACPATH}"
    #for ifarm
    if "ifarm" in os.uname()[1]:
	env+="%s G2P%s\"/w/work/halla/g2p/disk1\"\n"%setenv
	env+="%s PYTHONDIR%s/u/apps/python/python-2.7.1-shared\n"%setenv
	env+=ifenv(shelltype,"PATH","${PYTHONDIR}/bin:$G2P/software/pythonlib/bin")
	env+=ifenv(shelltype,"LD_LIBRARY_PATH","${PYTHONDIR}/lib:$G2P/software/pythonlib/lib")
	pythonpath+=":$G2P/software/pythonlib:$ROOTSYS/lib"
    env+=ifenv(shelltype,"PYTHONPATH",pythonpath)
    envfilename=os.path.join(beampacpath,"beampackage.%s"%shelltype)
    if not os.path.exists(beampacpath):os.makedirs(beampacpath)
    envfile=open(envfilename,"w")
    envfile.write(env)
    print "a environment file for your beampackage saved to %s"%envfilename
    def printattention(envfilename):
	for i in range(5):
	    print "\n-------------attention %i !!!-----------------"%(5-i)
	    print "\nplease add the line below to your %s file"%shellfile
	    print "\nsource %s"%envfilename
	    print "\n----------------------------------------------"
	    time.sleep(2)
    p1=threading.Thread(target=printattention,args=(envfilename,))
    p1.start()
    os.environ["BEAMPACPATH"]=beampacpath
    os.environ["BEAMDBPATH"]=beamdbpath

if __name__ == '__main__':
    if not checkenv():sys.exit()
    if os.getenv("BEAMUPDATECHECK")=="1":
	updateme=threading.Thread(target=updatemyself)
	updateme.start()
    if len(sys.argv)<2:
        usage()
        sys.exit()
    if len(sys.argv)<3:
        if sys.argv[1]=="update":
            print "try to update everything"
            updatecode("True")
            updatemyself("True")
            updateconst("True")
            print "everything are updated,please rerun code."
            sys.exit()
        elif sys.argv[1]=="env":
	    print "will regenerate a environment file for your beampackage"
	    makeenv(True)
	    sys.exit()
	elif sys.argv[1]=="changelog":
	    logfile=open(os.path.join(os.path.join(os.getenv("BEAMPACPATH"),"beampackage"),"changelog"),"r").read().split("-----")
	    for i in range(len(logfile)-1,-1,-1):
		print logfile[i]
	    sys.exit()
	elif "help" in sys.argv[1] or "-h" in sys.argv[1]:
	    usage()
	    sys.exit()
        runs=sys.argv[1:]
    else:
        runs=sys.argv[1:]
    rootpath=os.getenv("REPLAY_OUT_PATH")
    if rootpath==None:
	print "sorry please define REPLAY_OUT_PATH to your rootfile path"
	sys.exit()
    for run in runs:
        if run.isdigit():
            files=getfile(rootpath,int(run))
        else:
            if "/" in run:files=[run]
            else:files=[os.path.join(rootpath,run)]
        if len(files)<1:print "sorry no file found for run %s in %s"%(run,rootpath)
        for f in files:
	    if "ring" in f:treename="ring"
            else:treename="T"
            print "inserting bpm info for file %s now..."%f
	    bpminsert(f,treename)
