## Beam Package

Beam Package for Jefferson Lab experiment E08-027 (g2p).


### Feature

g2p beampackage is a python package developed for position calculation, it contains the package for:

* Harp data analysis from raw harp data
* BPM calibration
* Raster calibration and reconstruction
* Decode bpm information from rootfile
* Calculate position and dump to pickle file
* Insert position to rootfile
* Read run information, including:
 * BPM calibration constant
 * BCM calibration constant
 * Any run information saved in sqlite database which synced with Ryan's Mysql database plus some additional BPM information 

Also it includes several small tools to let life easier. For the principle on how to calibrate bpm please check the technote. I will update it if I think something is changed. 


### Download and install 

There have two parts for the beampackage, source code and database. You can copy them from the g2p work disk:/w/halla-sfs62/g2p/centos62/g2ppylib Once you copied it, you need to set the environments below:

* REPLAY_OUT_PATH used to find your rootfile
* BEAMPACPATH location for beampackage source code directory
* BEAMDBPATH location for beampackage database directory
* BEAMPKLPATH The position to save your pickle file, if there is corresponding pkl file exists beampackage will read it directly and skip decode or calculate
* BEAMPKLBAK The position to save your pickle files which you don't want to redecode or recalculate and you want stay it as stable version, beampackage will check BEAMPKLPATH first then check BEAMPKLBAK
* BEAMDEBUG debug bit, used to set 0
* BEAMUPDATE automatically update beampackage, you can set to 0 to close it
* BEAMUPDATECHECK check beampackage updated, you can set to 0 to close it
* PYTHONPATH you need to add your BEAMPACPATH to PYTHONPATH 

To use the beampackage, you need to make sure the package following is installed:

* python 2.7
* root with python module compiled from cern
* scipy,numpy and matplotlib 

The following modules should be came with python and you don't need to install any more: "sys","os","urllib","urllib2","socket","re","threading","math","pickle","shutil","glob","array","time","linecache","zipfile","ctypes","xml.etree.ElementTree","select","bisect","gc","multiprocessing"


### Database 

There are several type of databases you can find in BEAMDBPATH directory, ascii txt database, ods file, c code, pkl file and sqlite database

* sqlite database you can use command "sqlite3 g2p.db" to open and use mysql like command, here is the link for command line shell of sqlite. Also you can use GUI software like sqliteman, sqlite-manager to check the data so that you don't need to remember the command.
 * g2p.db converted from Ryan's mysql database, and added some additional BPM information 
* c code
 * bpm_$1_$2.c pure c code, they are fitted function used to transport the positions from two BPMs to target(x_a,y_a,x_b,y_b -> x_tgt,y_tgt,theta_tgt,phi_tgt) when there is target magnet field exists. The Makefile is also exists in the database directory used to compile them. Those functions are fitted by using BdL method and mudifi fit program, you can find them in g2p work disk: /w/halla-sfs62/g2p/pzhu/work/beamDrift, you need to install root and clhep and gcc,fortran before using it. $1 is the beam orbit which you can find from sqlite database and $2 is the integered target z position. 
* ascii file the normal ascii txt file, you can edit it by any txt editor
 * bpma/b_$run$.dat are the bpm calibration constants files. a/b shows which bpm it used, run number in filename is just a tag came from one of the bpm calibration run, it used to check which arm is available. The constant file contains the available run period info, available target z position info, bpm calibration constant information
 * raster_$run$.dat are the raster calibration constants files. The avail run period info, avail target z info, raster calibration constant are included in the file. run number in filename is just a tag same as bpm calibration constant file. They used to transfer the raster magnet current information to position offset caused by raster, which used to get the event by event position at target.
 * beampackage.conf defined several beampackage settings.
 * bpmerr.conf contains the bpm and harp survey information
 * bpmcalib.conf contains several special calibration settings for each calibration period 
* ods file the OpenDocument Spreadsheet file, you can open it by using openoffice or libreoffice
 * bcm calibration.ods saved bcm calibration constants for each periods, each devices, each arms
 * harppeaks.ods contains harp peak informations from raw harp data and bpm calibration informations for each period, and the event cuts informations for the bpm calibration runs
 * period.ods the beamline period during experiment. You can find same information from g2p.db 
* pkl file dumped from pickle module in python, you can load it in python code to reconstruct to a original python object structure
 * pedgain.pdt pedestals for two gains both from 0 to 40, div from 0 to 4, fitted from noise study at the end of run
 * bpmraws.pdt contains the average bpm raw ADC value for each run


### Source code 

The beampackage is built by using python 2.7, also includes several c functions in order to improve the efficiency by using ctypes.

* signalfilter.py This file contains several functions and classes, the main purpose of it is decode the bpm information from rootfile, add the software filter to the raw ADC signal, calculate the current from raw, judge the beam move. Also include the function to calculate position from raw.
* bpminsert.py The main purpose of this file is used to insert the bpm information to existed rootfile's "T" tree, calculate the positions at bpms and target. Also include the function to calculate the position at target from two bpms.
* bpmcalib.py This file provides the class to calibrate bpm, combine with "bpmfit.py"
* bpmfit.py This file is used to do the fit for bpm calibration
* harppos.py This file gives the classes and functions to calculate the beam position from harp peaks, and transfer it to several coordinates and transport the position to bpm.
* harppeak.py This file is used to automatically find the harp peak from raw harp datas and save it to rootfile
* rasterrecon.py This file is used to calibrate rasters, and reconstruct the raster shape.
* pacupdate.py This file is used to update the beampackage from wiki by using urllib2
* runinfo.py This file provides a class to read a lot of run informations and bpm calibration constant from database, it also provides an api to access sqlite database.
* epicst2v.py This file provides an api to read raw epics data in python
* odsread.py This file provides an api to read data from ods file
* bcmconst.py This file is used to read bcm calibration constant from bcm calibration.ods
* functions.c This file provides additional c function that used in python code to improve efficiency
* Cavestack.c This file is used in signalfilter.py, for averiging the signal. 

There are also additional python scripts for easier usage for beampackage existed in /w/halla-sfs62/g2p/centos62/g2ppylib:

* bpm.py script to decode,calculate,insert bpm position to rootfile, also will check if you have enough python module installed. 
 

### Dumped pkl Files 

The dumped pkl files from beampackage are located in $BEAMPKLPATH, also you can copy it to $BEAMPKLBAK if you want to share with others. The bpmraw_$key_$run.pkl files are dumped from signalfilter.py, which contain raw bpm data from rootfile. The bpmpos_$key_$run.pkl; foles are dumped from bpminsert.py, which are calculated position. For g2p experiment there are two DAQ system, one is normal DAQ system which triggered by Detector, another is Happex DAQ system which triggered by helicity signal. If the key have "rbpm","sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail","hapevent" all of these pkl files are aligned with helicity trigger and have aligned array index; and if the key have "raster","clock","event","fbbpm", all of them are aligned with normal detector trigger have aligned array index. More details are below:

* bpmraw_event_$run.pkl recorded the event number for each index, same index with "raster","clock","event","fbbpm" key. event number from leaf fEvtHdr.fEvtNum
* bpmraw_raster_$run.pkl recorded fast raster and slow raster ADC information from fastbus ADC
* bpmraw_hapraster_$run.pkl recorded fast raster and slow raster ADC information from Happex ADC
* bpmraw_clock_$run.pkl fast clock
* bpmraw_fbbpm_$run.pkl bpm raw ADC from fastbus, only available when the run is in fastbus only range defined in beampackage.conf
* bpmraw_hapevent_$run.pkl recorded the event number(fEvtHdr.fEvtNum) for each index, aligned with happex DAQ, same index with "rbpm","sbpm","ssbpm","fbpm","curr","hapevent","bpmavail","sbpmavail","fbpmavail","hapevent" key.
* bpmraw_curr_$run.pkl beam current
* bpmraw_rbpm_$run.pkl bpm raw ADC from happex, or from fastbus if happex not available but aligned to happex entry first before dump.
* bpmraw_fbpm_$run.pkl bpm raw ADC added frequency=filter2 low pass filter from bpmraw_rbpm_$run.pkl, filter1/2/3 is defined in beampackage.conf.It is used for getting slow raster shape.
* bpmraw_sbpm_$run.pkl bpm raw ADC added frequency=filter1 low pass filter. It is used for getting the average beam position
* bpmraw_ssbpm_$run.pkl bpm raw ADC added frequency=filter3 low pass filter. It is used for judging the beam move.
* bpmraw_bpmavail_$run.pkl saved the bpm information available for each entry, it is used for beam cut. In code all of the entries that current<20nA and 2000 events before and after that region are set to False, others are set to True.It is judged from beamraw_curr_$run.pkl file and filter1 value.
* bpmraw_fbpmavail_$run.pkl similar with bpmraw_bpmavail_$run.pkl but with filter2 value
* bpmraw_sbpmavail_$run.pkl similar with bpmraw_bpmavail_$run.pkl but with filter3 value
* bpmpos_sbpmabpm_$run.pkl calculated position for bpm A in bpm coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmabpm_$run.pkl calculated position for bpm A in bpm coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_ssbpmabpm_$run.pkl calculated position for bpm A in bpm coordinate with frequency=filter3 low pass filter from bpmraw_ssbpm_$run.pkl
* bpmpos_sbpmarot_$run.pkl calculated position for bpm A in bpm rotated coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmarot_$run.pkl calculated position for bpm A in bpm rotated coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_sbpmahall_$run.pkl calculated position for bpm A in hall coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmahall_$run.pkl calculated position for bpm A in hall coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_sbpmbbpm_$run.pkl calculated position for bpm B in bpm coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmbbpm_$run.pkl calculated position for bpm B in bpm coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_ssbpmbbpm_$run.pkl calculated position for bpm B in bpm coordinate with frequency=filter3 low pass filter from bpmraw_ssbpm_$run.pkl
* bpmpos_sbpmbrot_$run.pkl calculated position for bpm B in bpm rotated coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmbrot_$run.pkl calculated position for bpm B in bpm rotated coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_sbpmbhall_$run.pkl calculated position for bpm B in hall coordinate with frequency=filter1 low pass filter from bpmraw_sbpm_$run.pkl
* bpmpos_fbpmbhall_$run.pkl calculated position for bpm B in hall coordinate with frequency=filter2 low pass filter from bpmraw_fbpm_$run.pkl
* bpmpos_tgt$z_$run.pkl calculated target position and angle in hall coordinate, theta is dy/dz and phi is dx/dz. $z is integered target z(mm).
* bpmpos_rms_$run.pkl rms change information for position, used to judge for beam sharp move.
* bpmpos_split_$run.pkl beampackage will split events to several relative stable region, any beam sharp move and beam slow move will be splitted. This file saved the event region for relative stable region. 


### Insert Leaves to Rootfile 

For the considering of rootfile size and insert speed, only several pkl files will insert to "T" tree in rootfile, and only z=0 and additional z which depends on the target type recorded in sql database will be inserted:

* Lrb.bpmavail from bpmraw_bpmavail_$run.pkl, used for beam cut, >0 means beam available, =0 means not
* Lrb.curr from bpmraw_curr_$run.pkl, shows the beam current.
* Lrb.tgt_$z_x beam x position at target z(mm). If raster available, the final result=average beam position+fast raster caused offset+slow raster caused offset, used for position event by event. of raster not available, only average beam position is used.
* Lrb.tgt_$z_y beam y position at target z(mm) event by event, same definition as x
* Lrb.tgt_$z_theta beam theta angle(dy/dz) at target z(mm) event by event, same definition as x
* Lrb.tgt_$z_phi beam phi angle(dx/dz) at target z(mm) event by event, same definition as x
* Lrb.tgtave_$z_x average beam x position at target z(mm), only use freq=filter1 filter bpm information.


### Contact 

If you have any question please email me: zhupengjia@gmail.com 
