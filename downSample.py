#!/usr/bin/env python 
import os
import argparse
import h5py
import sys
import subprocess
import numpy as np
ap = argparse.ArgumentParser(description="down sample to an expected coverage for wgs")
ap.add_argument("--genomeSize", type=int ,default=200000)
ap.add_argument("fofn", help="fofn file with list of bax.h5 files to downsample")
ap.add_argument("--outdir", default="downSample")
ap.add_argument("--wgsCoverage", type=int ,default=50)
args = ap.parse_args()
genomeSize = args.genomeSize
fofn = args.fofn
outdir = args.outdir
wgsCov = args.wgsCoverage
downSampleScript="~mchaisso/projects/blasr/cpp/pbihdfutils/bin/writeHDFSubset"
debug = False 
print("parse and import over")


def runCmd(cmd):
    #print(cmd)
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    err = err.decode("utf-8")
    if(len(err) > 1):
        print(cmd)
        print(err)
    return(str(out.decode("utf-8")))



def getReadLength(fofn):
    cmd = "readLength.py " + fofn 
    rtn = runCmd(cmd)
    #print(rtn)
    check="Total number of bases"
    rtn = rtn.split("\n")
    total = []
    for line in rtn:
        line = line.split(":")
        if(line[0] == check):
            #print(line)
            num = float(line[1].strip())
            total.append(num)
    return(np.array(total))

def getHoles(bax):
    baxH5File = h5py.File(bax)
    regions = baxH5File["PulseData"]["Regions"]
    r = regions[:]
    holeNums = {}
    for i in range(len(r)):
        h = r[i][0]
        if(h not in holeNums):
            holeNums[h] = h
    rtn = np.array( holeNums.keys() )
    return(rtn)



def getIndexes(fofn, fractionToKeep):
    f = open(fofn)
    holes = []
    for bax, frac in zip(f, fractionToKeep):
        bax = bax.strip()
        hole = getHoles(bax)
        # sub selections
        numToSelect = int( frac * len(hole) ) 
        subset = np.sort(np.random.choice(hole, numToSelect , replace=False))
        print(subset)
        print(len(hole), len(subset), frac)
        holes.append(subset)
    return(holes)

def subset(holes, fofn):
    runCmd("mkdir -p " + outdir)
    files = open(fofn).read().split("\n")
    counter = 0
    for hole, f in zip(holes, files):
        counter += 1 
        print(f)
        cmd = downSampleScript +" "+ f + " " + outdir + "/" + str(counter) + ".bax.h5"
        for num in hole:
            cmd += " " + str(num)
        runCmd(cmd)

def main():
    np.random.seed(1)
    if(debug):
        fractionToKeep = np.array([0.06622517,0.06127451,0.05107252])
    else:
        totals = getReadLength(fofn)
        #print(totalLength)
        expCov = totals/genomeSize
        fractionToKeep = wgsCov * 1.0 / expCov 
    holes = getIndexes(fofn, fractionToKeep)
    subset(holes,fofn)
    
    copyfofn="cp " + fofn + " " + outdir + "/preDownSample.fofn"
    runCmd(copyfofn)
    newfofn = "cd "+outdir+" && ls *.h5 | xargs -n 1 readlink -f > downSample.fofn"
    runCmd(newfofn)




main()
print("done")



