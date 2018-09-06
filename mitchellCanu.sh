#!/bin/bash
unset PYTHONPATH
# needed by blasr
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/net/eichler/vol5/home/mchaisso/software/hdf/lib

module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
export PATH=/net/eichler/vol2/home/mvollger/anaconda3/bin:$PATH

source activate abp-python-3
export PATH=/net/eichler/vol2/home/mvollger/anaconda3/envs/abp-python-3/bin:$PATH



INPUT=$1
THREADS=$(nproc)

canu -pacbio-raw $INPUT \
	genomeSize=1200000 \
	corOutCoverage=9999 \
	corMhapSensitivity=high \
	corMinCoverage=1 \
	gnuplotTested=true  \
	useGrid=false  \
	maxThreads=$THREADS cnsThreads=$THREADS ovlThreads=$THREADS mhapThreads=$THREADS \
	contigFilter="25 5000 1.0 .75 25" \
	-d mitchell.assembly \
	-p asm 
	
# contigFilter, min total reads | min total length | max fraction covered by a single read | fraction of contig that must have min coverage | min coverage


# params to try 
"""
minReadLength <integer=1000>
Reads shorter than this are not loaded into the assembler. Reads output by correction and trimming that are shorter than this are discarded.

Must be no smaller than minOverlapLength.

If set high enough, the gatekeeper module will halt as too many of the input reads have been discarded. Set stopOnReadQuality to false to avoid this.

minOverlapLength <integer=500>
Overlaps shorter than this will not be discovered. Smaller values can be used to overcome lack of read coverage, but will also lead to false overlaps and potential misassemblies. Larger values will result in more correct assemblies, but more fragmented, assemblies.

Must be no bigger than minReadLength.
"""


