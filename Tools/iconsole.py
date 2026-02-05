# An interactive console to
# - compile SeisJIMU [not for now]
# - write setup.in
# - run the program (only 1 process?)
# - return results

from os import path, system
import subprocess
from sys import exit
import numpy as np


#==================
# prepare setup.in
#==================
setup={}

def setup_update(new):
    setup.update(new)

def setup_print(): #formatted print
    print('setup.in will include:')
    for key in setup:
        print(key,'\t',setup[key])

#==================
# run
#==================
exedir='/home/zhouw/Codes/GitHub/SeisJIMU/exe/'

def run(#executable
        exe=None,
        verbose=2
        ):
    """
    Run 
    :exe: directory to the executable
    #:return: directory to the result folder
    """
    
    print('Writing setup.in ...')
    fp=open("setup.in","w")
    for key in setup:
        fp.write(f"{key}\t\t{setup[key]}\n")  
    fp.close()
    
    #executable
    exefull=exedir+exe
    if not path.isfile(exefull) : exit(f"Executable {exefull} does not exist! Stop now")
    
    print('Running ...')
    sp=subprocess.run([exefull+' setup.in'],shell=True,
                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if verbose==2 : print(sp.stdout.decode('utf-8'))
    if sp.returncode<0 : exit('Run with error! Stop now')
    print("Finished running")
    
    print('Finish with success.')
    
    #return setup['DIR_OUT'] if 'DIR_OUT' in setup else None


def load(fname,n=None,format='bin',dir=None):
    
    if dir==None:
        dir=setup['DIR_OUT'] if 'DIR_OUT' in setup else './results/'
    
    fname=dir+fname
    
    if fname[-2:]=='su': format='su'
    
    if format=='bin':
        data=np.transpose(np.fromfile(fname,dtype='float32').reshape((n[1],n[0])))
    
    if format=='su':
        
        print(fname)
        n=np.fromfile(fname,dtype='int16',count=58)[-1]
        data=np.fromfile(fname,dtype='float32').reshape((-1,60+n))
        data=data[:,60:60+n]
        data=np.transpose(data)
    
    return data


def clean(reserve=[]):
    list=[
        'setup.in',
        'shotlist',
        ]
    
    list.append(setup['DIR_OUT'] if 'DIR_OUT' in setup else './results/')
        
    for item in list:
        if(not item in reserve):
            print('rm -r '+item)
            system('rm -r '+item)
    
