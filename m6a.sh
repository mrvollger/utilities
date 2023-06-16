#!/bin/bash
set -euo pipefail
export bam=$1
export T=8 # threads
export TMP="/tmp/$(whoami)/m6a"
export ODIR="m6a-bb"
rm -rf $TMP/*

mkdir -p $TMP $ODIR

samtools idxstats $bam | grep -v '^*' | cut -f 1,2 | sort -u > $TMP/ref.fai

make_m6a_bed () {
  chrom=$1
  hap=$2
  if [ $hap == "pat" ]; then
    hp="-d HP:1"
  elif [ $hap == "mat" ]; then
    hp="-d HP:2"
  elif [ $hap == "unk" ]; then
    hp=$'-e ![HP]'
  fi
  #>&2 echo "${hap} ${hp}"
  samtools view -@ $T ${hp} -u ${bam} ${chrom} | ft extract -r -t $T - --m6a ${TMP}/m6a.${chrom}.hp.${hap}.tbed
}
export -f make_m6a_bed

for hap in "unk" "pat" "mat"; do
    cut -f 1 $TMP/ref.fai | parallel -n 1 "make_m6a_bed {} ${hap}"
    cat $(ls $TMP/m6a.*.hp.${hap}.tbed | sort -u) > $TMP/${hap}.m6a.${hap}.bed
    bedToBigBed -allow1bpOverlap $TMP/${hap}.m6a.${hap}.bed $TMP/ref.fai $ODIR/m6a.${hap}.bb
    rm $TMP/m6a.*.hp.${hap}.tbed $TMP/${hap}.m6a.${hap}.bed
done

