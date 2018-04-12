#!/usr/bin/env python
""" SLURM usage summary 

A simple way to do this is to use Python 2.6+ for development and begin each of your Python .py files with the following:

https://www.dwheeler.com/essays/python3-in-python2.html
"""

from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division

import sys,os
import json
from  pprint import pprint
import copy
import datetime

import subprocess

__author__ = "Jan Balewski"
__email__ = "balewski@lbl.gov"


def timeStr2sec(st): #------------------
    
    x0=st.split('-')
    #print (st, len(x0 ))    
    tday=0
    if len(x0)==2:
        tday=int(x0[0])*24*3600
        x1 = x0[1].split(':')    
    else:
        x1 = x0[0].split(':')    
    
    t=0.
    try:
        if len(x1)==3:
            t=float(x1[0])*3600. +float(x1[1])*60 +float(x1[2])
        elif len(x1)==2:
            t=float(x1[0])*60 +float(x1[1])
        else:
            t=float(x1[0])
    except:
        #print('failed to convert time, assume 0, inp=',x1,t)
        a=1

    return t+tday
                    
#------------------------------
def  niceAccountTable(dataD,text1):
    keyL=['Rjob', 'Rcpu', 'Rcpu*h', 'PDjob','PDcpu']
    for x in keyL:
        print('%7s  '%x,end='')

    print('    '+text1)
    
    men=sorted(dataD.keys())
    men.append(men.pop(0))
    
    for us in men:
        if us=='  TOTAL': print('')
        for x in keyL:
            #print(x)
            if x=='Rcpu*h': print('%7.1f  '%dataD[us][x],end='')
            else:            print('%7d  '%dataD[us][x],end='')
        print('    %s'%us)

    print('')



#------------------------------
def  agregatePartitionData(byAcct):

    byPart={}
    for x in byAcct:
        #print(x, byAcct[x])
        if 'TOTAL' in x: continue
        partN=x.split()[1]
        if partN not in byPart : byPart[partN]={'Rjob':0, 'Rcpu':0, 'Rcpu*h':0 ,'PDjob':0,'PDcpu':0}
        for y in  byAcct[x]:
            byPart[partN][y]+=byAcct[x][y]
    #print('ppp',byPart)
    return byPart


#------------------------------
def  agregateAccountData(byUser):

    total0={'Rjob':0, 'Rcpu':0, 'Rcpu*h':0 ,'PDjob':0,'PDcpu':0}
    total=copy.deepcopy(total0)
    byAcct={}
    for u in byUser:
        user=byUser[u]
        acct=user['acct-part']
        #print(user,acct)
        if acct not in byAcct : byAcct[acct]=copy.deepcopy(total0)
        for x in  byAcct[acct]:
            byAcct[acct][x]+=user[x]
            total[x]+=user[x]

    byAcct['  TOTAL']=total
    #pprint(byAcct)
    return byAcct


#------------------------------
def  niceUserTable(dataD):
    keyL=['Rjob', 'Rcpu', 'Rcpu*h', 'PDjob','PDcpu']
    for x in keyL:
        print('%7s  '%x,end='')

    print('    user:account:partition')
    for us in sorted(dataD.keys()):
         for x in keyL:
            #print(x)
            if x=='Rcpu*h': print('%7.1f  '%dataD[us][x],end='')
            else:            print('%7d  '%dataD[us][x],end='')
         print('    %s'%us)

    print('')

#------------------------------
def agregateUserData( allD):
    byUser={}

    for job in allD:
        #print('ana job:',job)
        pref=job['ST']
        if pref=='CG' : continue

        user=job['USER']
        acct=job['ACCOUNT']
        part=job['PARTITION']
        #print('uu',pref,job['TIME'])
        timeSec=timeStr2sec(job['TIME'])
        name=user+' '+acct+' '+part
        if name not in byUser: byUser[name]={'Rjob':0, 'Rcpu':0, 'Rcpu*h':0 ,'PDjob':0,'PDcpu':0,'acct-part':acct+' '+part}

        nCpu=int(job['CPUS'])
        #print('aa',pref)
        byUser[name][pref+'job']+=1
        byUser[name][pref+'cpu']+=nCpu
        if pref=='R':  byUser[name]['Rcpu*h']+= nCpu*timeSec/3600.
       
        #break
    #pprint(byUser)
    return byUser    

#--------------------------------------------------
def checkEnv():
    #print('checkEnv():')
    #task = subprocess.Popen("module list  2>&1 >/dev/null ", shell=True, stdout=subprocess.PIPE)
    task = subprocess.Popen("/usr/share/modules/3.2.10/Modules/$MODULE_VERSION/bin/modulecmd tcsh list  2>&1 >/dev/null ", shell=True, stdout=subprocess.PIPE)
    data = task.stdout.read()
    #print(task.wait(),'dataByte=',data)
    assert task.wait() == 0
    #
    dataStr=data.decode('utf-8')
    #dataStr=data.encode('ascii', 'ignore')
    #print('dataStr=',dataStr)
    if 'slurm' not in dataStr:
        print('aborted, execute first:\n module load slurm')
        exit(3)
    #print('slurm environment is set up')

#--------------------------------------------------
def scanSlurmJobs():
    #print('scanSlurmJobs():')
    #task = subprocess.Popen("squeue --array", shell=True, stdout=subprocess.PIPE)
    task = subprocess.Popen('squeue --array --format="%all"', shell=True, stdout=subprocess.PIPE)
    data = task.stdout.read()
    assert task.wait() == 0
    #print('dataByte=',data)
    dataStr=data.decode('utf-8')
    #dataStr=data.encode('ascii', 'ignore')
    #print('dataStr=',dataStr)
    n=0
    #1 -> ['JOBID', 'USER', 'ACCOUNT', 'NAME', 'ST', 'REASON', 'START_TIME', 'TIME', 'TIME_LEFT', 'NODES', 'CPUS', 'PARTITION', 'PRIORITY']
    #nameL= ['USER', 'ACCOUNT', 'ST', 'TIME', 'CPUS', 'PARTITION']

    if 0: # debug on raw dump of : squeue --array
        fd=open('../../0x/squeue.log2','r')
        dataL = fd.readlines()
        print ("Read input len" ,len(dataL),dataL[:4])
        dataStr=''.join(dataL)
        #exit(3)

    outD=[]
    idxL=[8,20,0,19,38,28,41]
    #idxL=[0,1,2,4,7,10,11]
    for line in dataStr.split('\n'):
        
        if len(line)<2: continue
        if 'launch failed ' in line: continue
        n+=1
        xL=str(line).split('|')
        #print('lineL=',xL)
        if n<0: print(n,'->',xL)
        if n==1:
            nameL=[]
            for x in idxL:
                nameL.append(xL[x])
            #print('nameL=',nameL)
            continue
        rec={}
        #print('kk',xL)
        for i,x in enumerate(idxL):
            #print(i,x)
            rec[nameL[i]]=xL[x]
        #print('rr',rec)
        outD.append(rec.copy())
        #if n>10:        break
    return outD

################################
#     MAIN
################################

if __name__ == '__main__':
    #print ('Number of arguments:', len(sys.argv)-1)
    #print ('Argument List:', str(sys.argv[1:]))
    dateStop=datetime.datetime.now()
    dateNowStr=dateStop.strftime("%Y-%m-%d_%H.%M")
    #print('dump env:',os.environ)

    print ("\n %s SLURM usage, all PDSF users, ver 1.7\n"%dateNowStr)
    if  len(sys.argv)>1:
        print (" Columns:  Rjob     Rcpu   Rcpu*h    PDjob    PDcpu ")
        print (" R - running,  PD - pending ")
        print (" job - slurm job count, arrays are unrolled")
        print (" cpu - task on node count, one job may lock many tasks")
        print (" cpu*h - sum of CPU*hours used by jobs in execution")
        exit(1)


    #checkEnv()  # not working in crontab
    allD=scanSlurmJobs()
    #print('M: num accnts=',len(allD))
    #pprint(allD)
    byUser=agregateUserData( allD)    
    niceUserTable(byUser)
    byAcct=agregateAccountData( byUser)
    byPart=agregatePartitionData( byAcct)
    niceAccountTable(byAcct,'account + partition')
    niceAccountTable(byPart,'account')

'''
copy new version to
ssh -I /Library/OpenSC/lib/opensc-pkcs11.so sg-crt.nersc.gov

ssh root@mc0154-ib.nersc.gov

cd /chos/common/nsg/slurm/17.02.7/bin 
mv slusers.py old-slusers.py

scp -rp balewski@pdsf:janNersc/pdsfVaria/slusers.py .

chmod a+x slusers.py
chmod a+r slusers.py

'''
