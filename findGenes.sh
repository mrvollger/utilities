#!/bin/bash

grep "$1" */ref.fasta.genes.bed -l | awk -F "/" '{print $1"/ref.fasta.bed"}' | xargs wc -l


for x in $( grep "$1" */ref.fasta.genes.bed -l | awk -F "/" '{print $1}'  ); do
	echo "ln -sf $PWD/$x testCases/$1"
done




