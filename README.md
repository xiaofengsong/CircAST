# CircAST

# 1. Introduction
CircAST v1.0.0 (CircRNA Alternative Splicing Transcripts) is a computational tool to assemble and quantify circular transcripts using all the back-spliced events detected by circRNA identification tools. 
Original codes and detailed manuals are included in the downloading URL: ???

Authors: Jing Wu (jingwu26@163.com), Xiaofeng Song (xfsong@nuaa.edu.cn) 
Maintainer: Jing Wu (jingwu26@163.com)

# 2. Usage
	CircAST works with three input files. A GTF annotation file, a Sequence Alignment/Map (SAM) of pair-end reads generated by tophat2, and a circRNA junction list containing the information of the back-spliced junctions which is generated by upstream identification tools. 
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
		  SAM alignment of pair-end reads which is generated by tophat2 and should be sorted by name
	-J,	--junction
		  circRNA junction list generated by CIRCexplorer2, CIRI, UROBORUS or other tools in same format
	-H,	--help
          	  show this help information
		  
# 3. Note
*	1). CircAST is implemented in Python2 under Linux system. Since several packages are required for using CircAST, we suggest that users install corresponding version of Anaconda. 
*	2). The accepted_hits.bam file generated by tophat2 should be sorted by name first, then be converted to accepted_hits.sorted.sam file. 
*	3). CircAST is directly compatible with popular upstream circRNA detection tools (CIRCexplorer2, CIRI2, UROBORUS). Other circRNA junction detection tools can also be used in CircAST if their output file format is transformed to CircAST input format files.
*	4). CircAST use the GTF annotation file with UCSC compatable GTF format.

# 4. Example
	python CircAST.py -L 100 -T 5 -G /home/***/circRNA/Gene/genes.gtf -F accepted_hits.sorted.sam -J circRNA_list.txt

# 5. Prerequisites
Software Prerequisites:
The following three software should be installed in your cluster or computer before running the CircAST.py
* TopHat2
	tophat -o results -p 12 -G /home/***/genes_homo.gtf /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome reads_1.fastq reads_2.fastq
* samtools sort accepted_hits.bam file by name and then convert accepted_hits.sorted.bam to accepted_hits.sorted.sam (using samtools sort and samtools view)
* samtools sort -n -o accepted_hits.sorted.bam accepted_hits.bam
* samtools view accepted_hits.sorted.bam > accepted_hits.sorted.sam
* CIRCexplorer2/CIRI2/UROBORUS
* corresponding circRNA list: circ_out/denovo/annotated_circ.txt, result.txt, circRNA_list.txt

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
11) Expression level(FPKM)
12) Expression level(TPM)
Columns of output file are split by tabs ("\t" in python).

# 7. Reference

# 8. Contact
Please contact Jing Wu (jingwu26@163.com) for questions and comments.
