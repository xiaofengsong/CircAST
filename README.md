# CircAST

# 1. Introduction
CircAST v1.0.0 (CircRNAs with Alternative Spliced Transcripts) is a computational tool for circular RNA full-length assembly and quantification using RNA-Seq data with the back-spliced events detected by upstream circRNA identification tools, such as UROBORUS, CIRI2, or CIRCexplorer2. Original codes and detailed manuals are included in the downloading URL: https://github.com/xiaofengsong/CircAST/archive/v1.0.0.tar.gz.

Authors: Jing Wu (wujing@njmu.edu.cn), Xiaofeng Song (xfsong@nuaa.edu.cn) 
Maintainer: Jing Wu (wujing@njmu.edu.cn)


# 2. Prerequisites
1). Software Prerequisites:
    CircAST is implemented in Python2 under Linux system. Since several packages are required for using CircAST, we suggest that users install corresponding version of Anaconda(Anaconda 2.5.0 +, Python 2.7 version).

2). Input files:
    CircAST works with three input files. A GTF annotation file, a sorted Sequence Alignment/Map (SAM) of pair-end reads generated by tophat2 and samtools, and a circRNA junction list containing the information of the back-spliced junctions which is generated by upstream identification tools such as CIRCexplorer2/CIRI2/UROBORUS or other similar tools.

    SAM file sotsorted by name can be obtained by the following commands:
	* tophat -o results -p 12 -G /home/***/genes_homo.gtf /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome reads_1.fastq reads_2.fastq
	* samtools sort -n -o accepted_hits.sorted.bam accepted_hits.bam
	* samtools view accepted_hits.sorted.bam > accepted_hits.sorted.sam

    CircRNA lists containing back-spliced junctions can be obtained by CIRCexplorer2/CIRI2/UROBORUS,  circ_out/denovo/annotated_circ.txt, result.txt, or circRNA_list.txt. CircRNA lists transformed into such formats using other circRNA back-spliced site detection software are also also feasible. 

# 3. Usage 
	Command:
		python CircAST.py -L length -T thresh -G GTF -F accepted_hits.sorted.sam -J junction.txt
	Options:
	-L,	--length
		  reads length
	-T,	--threshold
		  set threshold for back-splicing junction supporting reads(optional, default: 10)  
	-G,	--gtf
		  gene annotation file (*.gtf file)
	-F,	--file
		  SAM alignment of pair-end reads generated by tophat2 and should be sorted by name
	-J,	--junction
		  circRNA junction list generated by CIRCexplorer2, CIRI, UROBORUS or other tools in similar format
	-H,	--help
          	  show this help information
		  
# 4. Note
	1). The accepted_hits.bam file generated by tophat2 should be sorted by name first, then be converted to accepted_hits.sorted.sam file. 
	2). CircAST is compatible with popular upstream circRNA detection tools (CIRCexplorer2, CIRI2, UROBORUS). Other circRNA junction detection tools can also be used in CircAST if their output file format is transformed to CircAST input format files.
	3). CircAST use the GTF annotation file with UCSC compatable GTF format.

# 5. Example
	python CircAST.py -L 100 -T 5 -G /home/***/circRNA/Gene/genes.gtf -F accepted_hits.sorted.sam -J circRNA_list.txt

# 6. Output file format
	The first 3 columns are the same with bed file format.
	1) Chromosome
	2) Start of circular RNA
	3) End of circular RNA
	4) Name of Circular transcript
	5) Parental gene name
	6) Strand
	7) Number of exons
	8) Exon sizes
	9) Exon offsets
	10) Number of junction reads
	11) Expression level (FPKM)
	12) Expression level (TPM)
	13) Expression level(Read counts)
	Columns of output file are split by tabs ("\t" in python).

# 7. Reference

# 8. Contact
	Please contact Jing Wu (wuijng@njmu.edu.cn) for questions and comments.

Copyright (C) 2018 Xiaofeng Song.
