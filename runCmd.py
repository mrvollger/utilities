#!/usr/bin/env python
import subprocess
import sys
import os
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


def exe(cmd, printme=False):
    if(printme):
        print("starting", cmd)
    
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    err = err.decode("utf-8")
    out = out.decode("utf-8")
    
    return(out, err)


