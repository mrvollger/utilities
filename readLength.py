#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputFile', nargs='?', help="autodetects plus does regex, but the regex requires quotes around it", default=None)
parser.add_argument("--fasta", help="reads file fasta")
parser.add_argument("--h5", help="reads h5 file")
parser.add_argument("--fofn", help="file of bas/bax.h5 files")
parser.add_argument("--printn", help="print the lenght of the first n reads", default="0")
#parser.add_argument('--selfMap', help="map the bac assmbly to itsefl", action='store_true')
#parser.set_defaults(selfMap=False)
args = parser.parse_args()

#selfMap = args.selfMap
fofn = args.fofn
h5 = args.h5
fasta = args.fasta
inputFile=args.inputFile
printn = int(args.printn)

# other imports, a little slow 
import glob
import subprocess
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import h5py
import numpy as np
import multiprocessing as mp
# global varibles

# convert h5 to fasta
def runPls2Fasta(reads, outreads):
    cmd = "/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/pbihdfutils/bin/pls2fasta "
    cmd_options = reads + " " + outreads + " -trimByRegion"
    proc = subprocess.Popen(cmd + cmd_options, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    if(len(out.strip()) > 0 or len(err.strip()) > 0): 
        print(out)
        print(err)

def fastaLengths(fastareads):
    lengths = []
    names = []
    records = list(SeqIO.parse(fastareads, "fasta")) 
    #print(len(records))
    for read in records:
        #print(name)
        name = read.name
        seq = read.seq
        seqlen = len(seq)
        if(seqlen > 0):
            lengths.append(seqlen) 
            names.append(name)
    return names,  lengths


def printLengths(names, lengths):
    if(len(names) == 0 or len(lengths) == 0):
        return
    printx = min(len(names), printn)
    printLater = ""
    for i in range(printx):
        s = "{}\t{}".format(lengths[i], names[i]) + "\n"
        printLater += s
   
    # get average length   
    totalseqs = len(lengths)
    averageLen = np.mean(lengths)
    # make histogram
    
    window  = 1000 
    windows = int(max(lengths) / window) + 1 
    #print(windows, max(lengths), window)
    hist = {}
    s_hist = "Length histogram:   * => 1-2%    . => <1%  \n"
    for i in range(windows):
        start = i * window
        end = (i+1) * window
        hist[i] = 0
        for seq in lengths:
            if(seq >= start and seq < end ):
                hist[i] += 1
        # determine how many * to add
        if(hist[i] > 0): 
            num = int( hist[i]/(totalseqs/100+1) )
            s_hist += "{}-{}:{}\t\t".format(start, end, hist[i])
            if(num == 0):
                s_hist += "."
            else:
                for i in range(num): s_hist += "*"
            s_hist += "\n"

    # calculate n50
    # this sorts the by length 
    lenName = sorted(zip(lengths, names), reverse=True)
    lengths, names = zip(*lenName)
    totalseq = np.sum(lengths)  
    sofar = 0
    N50=0
    for seq in lengths:
        sofar += seq
        if(sofar >= totalseq/2 ):
            N50 = seq 
            break
    
     
    rtn =  (printLater + s_hist[:-1]) + "\n"
    rtn += ("Number of sequences:\t{}".format(totalseqs)) + "\n"
    rtn += ("Total number of bases:\t{}".format(totalseq))+ "\n"
    rtn += ("Average Length:\t{0:.2f}".format(averageLen))+ "\n"
    rtn += ("N50:\t{}".format(N50)) + "\n"
    return(rtn)

def holeGroup(data): 
    holeIdx = data[:,0]
    reads = []
    for holeNum in sorted(set(holeIdx)):
        hole = data[holeIdx == holeNum,:]
        hqRegion = hole[ hole[:,1] == 2, : ]
        # skip fi there is no hq region, idk if this ever happens
        if(len(hqRegion) == 0):
            continue
        # make sure there is only one hqregion, pretty sure it should be unique
        assert len(hqRegion) == 1
        hqRegion = hqRegion[0]
        hqStart = hqRegion[2]
        hqEnd =   hqRegion[3]
        inserts = hole[ hole[:,1]==1,: ]
        for insert in inserts:
            start = 2
            end =   3
            if( insert[start] < hqStart):
                insert[start] = hqStart
            if( insert[end]   > hqEnd):
                insert[end]   = hqEnd
            reads.append(insert)

    return(reads)

def h5Lengths(h5): 
    h5file = h5py.File(h5, 'r')
    # top level keys: PulseData ScanData
    # Scan data does not have anything useful for us
    # next level BaseCalls Regions 
    # pretty sure we want the dir PulseData/Regions
    data = h5file["PulseData"]["Regions"] 
    # creates a copy
    data = data[:] 
    h5file.close()
    trimmedInserts = holeGroup(data)
    names = []
    lengths = []
    start = 2; end = 3
    for insert in trimmedInserts:
        length = insert[end] - insert[start]
        if( length > 0 ):
            lengths.append(length)
            names.append(insert[0])

    return(printLengths(names, lengths))
    
    '''  # the slow lame way
    tempFasta = "ForReadLengthTemp.fasta"
    runPls2Fasta(h5, tempFasta)
    names, lengths = fastaLengths(tempFasta)
    printLengths(names, lengths)
    os.system("rm ForReadLengthTemp.fasta")
    '''

def detectAndRun(readFile):
    readFileL = readFile.split(".")
    length = len(readFileL)
    ftype = readFileL[length -1]
    rtn = (readFile) + "\n"
    if(ftype == "fa" or ftype == "fasta"):
        names, lengths = fastaLengths(readFile)
        rtn += printLengths(names,lengths)
    
    elif(ftype == "h5" ):
        rtn += h5Lengths(readFile)
    
    elif(ftype == "fofn"):
        print("does not supported nested FOFNs")
        #rtn += multiFiles(readFile)
    return(rtn)

def expandFofnsInList( mylist ):
    expandedList = []
    anotherFofn = False
    for item in mylist:
        ftype = item.split(".")[-1]
        if(ftype == "fofn"):
            f = open(item)
            anotherFofn = True
            for line in f:
                expandedList.append(line.strip())
        else:
            expandedList.append(item)
    if(anotherFofn):
        expandedList = expandFofnsInList(expandedList)
    return(expandedList)


def multiFiles(ListOrFofn):
    mylist = []
    # this means it is a regex  
    if( type(ListOrFofn) == type([]) ):
        mylist = ListOrFofn
    else: # this means it is an fofn
        mylist = [ListOrFofn]

    mylist = expandFofnsInList(mylist)
    #print(len(mylist))
    pool = mp.Pool()
    rtn = pool.map(detectAndRun, mylist)
    for output in rtn:
        print(output)


def main():
    if(fasta):
        print(fasta)
        names, lengths = fastaLengths(fasta)
        printLengths(names,lengths)
    if(h5):
        print(h5)
        h5Lengths(h5)
    if(fofn):
        multiFiles(fofn)

    # this means there is an input file or it might be a regex for multiple files
    if(inputFile != None):
        myglob = glob.glob(inputFile)
        multiFiles(myglob)


main()
