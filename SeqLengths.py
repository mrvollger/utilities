#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(
    description="Get the lengths of reads from fasta/fastq/(.gz) files."
)
parser.add_argument("infiles", nargs="*", help="input fasta/fastq file")
parser.add_argument(
    "-f", "--fofn", help="input FOFN of fasta/fastq files, replaces infiles argument"
)
parser.add_argument(
    "-o",
    "--out",
    help="output lengths file. Has the form: readname, length",
    type=argparse.FileType("w"),
    default=sys.stdout,
)
parser.add_argument(
    "-p",
    "--plot",
    default=None,
    help="Specify a place to plot a hsitogram of read lengths",
)
parser.add_argument(
    "--plotwith",
    default=None,
    help="If you have already calcualted readlength before you can specify plotwith to avoid rereading the fasta files",
)
parser.add_argument(
    "-t",
    "--threads",
    type=int,
    default=1,
    help="number of threads to use, only used in reading in fasta/fastq files",
)
args = parser.parse_args()


# from Bio import SeqIO
import pysam
import gzip
from multiprocessing.dummy import Pool as ThreadPool

# from tqdm import tqdm
poses = {}
pos = 0


def getDesc(f, ftype):
    desc = "\033[92mFile:{} File Type:{} Progres\033[0m".format(f.name, ftype)
    return desc


def getPos(f):
    global poses, pos
    if f.name not in poses:
        pos += 1  # position for the tqdm display.
        poses[f.name] = pos
    return poses[f.name]


def readSeqs(f, ftype):
    lengths = []
    sys.stderr.write(f"file:{f}")
    if ftype in ["fasta", "fastq"]:
        # recs = list(SeqIO.parse(f,ftype))
        # for rec in tqdm(SeqIO.parse(F,ftype), position=pos, desc=getDesc(f,ftype)):
        fasta = pysam.FastxFile(f.name)
        for rec in fasta:
            # lengths.append((rec.name, len(rec.seq)))
            lengths.append((rec.name, len(rec.sequence)))
    elif ftype in ["bam"]:
        bamfile = pysam.AlignmentFile(f, check_sq=False)
        # for read in tqdm(bamfile.fetch(until_eof=True), desc=getDesc(f,ftype), position=getPos(f)):
        for read in bamfile.fetch(until_eof=True):
            lengths.append((read.query_name, read.infer_query_length()))
    elif ftype in ["bed"]:
        # for line in tqdm(f.readlines(), desc=getDesc(f,ftype), position=getPos(f)):
        for line in f.readlines():
            line = line.split()
            name, start, end = line[0], int(line[1]), int(line[2])
            lengths.append((name, end - start))
    return lengths


# for seqfile in args.infiles:
def readFile(seqfile):
    # get file type
    ftype = None
    tags = seqfile.split(".")
    if ("fasta" in tags) or ("fa" in tags):
        ftype = "fasta"
    elif ("fastq" in tags) or ("fq" in tags):
        ftype = "fastq"
    elif (seqfile[-4:] == ".bam") or (seqfile[-4:] == ".sam"):
        ftype = "bam"
    elif seqfile[-4:] == ".bed":
        ftype = "bed"
    assert ftype is not None, "File type unknown"

    # uncompress if needed
    if seqfile[-3:] == ".gz":
        f = gzip.open(seqfile, "rt")
    else:
        f = open(seqfile)

    # get read lengths
    lengths = readSeqs(f, ftype)
    # clsoe read file
    f.close()

    print(
        "\033[92mRead {} reads in {}\033[0m".format(len(lengths), seqfile),
        file=sys.stderr,
    )
    return lengths


def NX(nums, X):
    nums = sorted(nums, key=int, reverse=True)
    datathresh = sum(nums) * (X / 100.0)
    total = 0
    for num in nums:
        total += num
        if total >= datathresh:
            return num
    return 0


flat_list = None
if args.plotwith is None:
    myfiles = args.infiles
    if args.fofn is not None:
        myfiles = []
        for line in open(args.fofn).readlines():
            myfiles.append(line.strip())

    # create pool of to thread process.
    threads = min(args.threads, len(myfiles))
    pool = ThreadPool(threads)
    print("\033[92mUsed Threads:{}\033[0m".format(threads), file=sys.stderr)
    lengths = pool.map(readFile, myfiles)
    # output buffer aorund tqdm lines
    print(pos * "\n", file=sys.stderr)

    flat_list = []
    for l in lengths:
        flat_list += l

    out = ""
    for name, length in flat_list:
        out += "{}\t{}\n".format(name, length)
    args.out.write(out)
    plotwith = args.out


if args.plot is not None:
    import pandas as pd
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set(font_scale=2)
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(16, 9))

    if args.plotwith is not None:
        df = pd.read_csv(args.plotwith, sep="\t", header=None, names=["name", "length"])
    elif flat_list is not None:
        df = pd.DataFrame(flat_list, columns=["name", "length"])
    else:
        print("No data to plot with, try --plotwith", file=sys.stderr)

    df["LogLength"] = np.log10(np.clip(df["length"], 10, None))

    # make histogram
    g = sns.distplot(df.LogLength, bins=100, kde=False, rug=False, ax=ax)
    # print(g)
    vv, ll = np.histogram(df.LogLength, bins=1000)
    for v, l in zip(vv, ll):
        if 10**l > 14000 and 10**l < 18000:
            print(v, 10**l)

    # get maxes
    myy = plt.gca().get_ylim()[1]
    mymax = max(df["length"])

    # add vertical lines
    vals = [
        df["length"].median(),
        df["length"].mean(),
        NX(df["length"], 50.0),
        NX(df["length"], 1.0),
        mymax,
    ]
    names = ["Median", "Mean", "N50", "N1", "Max"]
    divs = [0.9, 0.8, 0.7, 0.6, 0.5]
    for name, val, div in zip(names, vals, divs):
        plt.axvline(x=np.log10(val), color="darkred", linestyle="dashed")
        label = "{}={}".format(name, round(val / 1000, 1))
        plt.text(np.log10(val) + 0.01, myy * div, label, fontsize=18)

    # make x ticks
    xts = [2, 3, 4, 5, 6, 7]
    xls = ["0.1", "1", "10", "100", "1,000", "10,000"]
    if np.log10(mymax) <= 6:
        xts = xts[:-1]
        xls = xls[:-1]

    # plot
    plt.xticks(xts, xls)
    Gb = sum(df["length"]) / 1000000000
    plt.xlabel("Length (kbp), Total Gb {:.2f}".format(Gb))
    plt.ylabel("Number of reads")
    plt.savefig(args.plot)
