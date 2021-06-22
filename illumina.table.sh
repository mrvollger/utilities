basedir=/net/eichler/vol24/projects/sequencing/nextseq500/nobackups/vyqtdang/Melanesian_megapool/fastq_split
index_suffix=_I1_001.fastq.gz


for fn in `ls $basedir | grep $index_suffix`; do
	sn=$( basename $fn $index_suffix )
	sn2=(${sn//_S/ });
   	echo -e "${sn2}\tbarcode\t${sn2}\t$basedir/${sn}_R1_001.fastq.gz,$basedir/${sn}_R2_001.fastq.gz"; 
done | \
   	sort -V | \
	sed '1iwell\tbarcode\tsample_name\treads'


