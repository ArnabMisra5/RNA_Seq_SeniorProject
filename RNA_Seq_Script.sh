#!/bin/bash
input1=$(ls /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Sample_*/*.fastq.gz)
input_files=(${input1[@]})
#echo $input_files[@]}
/home/sheetalshetty/Downloads/FastQC/fastqc --extract --outdir=/media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData ${input1}
#md5sum
#Generate checksum (unique value for a files contents)

Fastqc -o /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/PreAlignmentQC -t 16 /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData/*/*.fastq.gz
#fastqc run

STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Alignment_Files \
--genomeFastaFiles /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
--sjdbGTFfile /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/gencode.v45.chr_patch_hapl_scaff.basic.annotation.gtf --sjdbOverhang 99
#initial star test run

rawsFolder=/media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData
myFolderPaths=$(ls $rawsFolder)
myFolders=($myFolderPaths[@])
for name in ${myFolders[29]}; do
	echo $rawsFolder/$name
	fastq1=$rawsFolder/$name/*_R1*.fast.gz
	fastq2=$rawsFolder/$name/*_R2*.fast.gz
	echo $fastq1
	echo $fastq2
	mkdir -p /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Outputs/ReMapped/$name
	STAR --runThreadN 32 --genomeDir /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Alignment_Files \
	--readFilesIn $fastq1 $fastq2 /
	--readFilesCOmmand gunzip -c --outFileNamePrefix /media/sheetalshetty/NewDrive!/Arnab_PNPLA7/STAR_Outputs/ReMapped/$name/ --outSAMtype BAM unsorted
done
#STAR Alignmentcode

htseq-count --format=bam /media/sheetlashetty/NewDrive1/Arnab_PNPLA7/STAR_Outputs/*.bam /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/Homo_sapiens.GRCh38.111.gtf > AllHTSeq.tsv 
#htseq-count and gff/gtf files