# CRISPRAnalysis
A tool to analyze the NGS data of CRISPR/Cas9 gene editing experiment.

<h2>Run the CRISPRAnalysis pipeline</h2>

#### 1. Clone the repository

        git clone https://github.com/yuht4/CRISPRAnalysis.git

The programme file is in the folder CRISPRAnalysis/ 


#### 2. Install prerequisite tools    

    
- BWA: 

    Example install by non-root user:
    
    	cd $HOME
	    git clone https://github.com/lh3/bwa.git
        cd bwa; make
	  
    Be sure to put executable 'bwa' in your PATH, for example, by adding this line to $HOME/.bashrc:
    
	    export PATH=$HOME/bwa:$PATH
 	
- samtools:

    Example install by non-root user:

	    cd $HOME
	    wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
	    tar -xvfj samtools-1.7.tar.bz2
	    ./configure
	    make
	    
    Be sure to put executable 'samtools' in your PATH, for example, by adding this line to $HOME/.bashrc:
  
	    export PATH=$HOME/samtools-1.7:$PATH
	    

- FASTX-Toolkit:

    Example install by non-root user:

	    cd $HOME
	    git clone https://github.com/agordon/fastx_toolkit
	    cd fastx_toolkit
        ./reconf
        ./configure
        make
	    
    Be sure to put executable 'fastx_toolkit' in your PATH, for example, by adding this line to $HOME/.bashrc:
  
	    export PATH=$HOME/fastx_toolkit:$PATH


- Java:

    Your running environment must have Java 1.8 or higher version installed. Run the following command to check the java version:
    	    
	    java -version

- matplotlib
    
    Your running environment must have matplotlib python package.
    Example install by non-root user:

        pip install matplotlib
	    
	    
#### 3. Preparing input files for the pipeline

- Fastq files:

    Prepare the raw fastq files by whole genome sequencing, with file extension .gz.
    For example:

        crispr_sorted_remdup_rg_1.fastq.gz
        crispr_sorted_remdup_rg_2.fastq.gz


- Reference genome:

 Prepare your reference genome fasta files:
  For example, to parepare human genome hg38.fa, run the following commands:
    
 1. Create bwa index:

        
            bwa index hg38.fa
    
    
 2. Create fasta index:
 
            samtools faidx hg38.fa

- Experiment info file:
    
    A csv file with 7 or 8 columns for: chr, start, end, cutsite, genome, sgRNA, PAM, HDR,   separated by ';'
If there is no HDR sequence involved, just not fill it in the csv file.

  This file can contain multiple rows.
    
          chr:        the crispr on-target chromosome
          start:      the -30bp position of cutsite position
          end:        the +30bp position of cutsite position
          cutsite:    the cutsite name
          sgRNA:      the sgRNA sequence
          PAM:        the PAM sequence
          HDR:        the HDR repair sequence.
    

   for example:
    
        chr3;46414647;46414707;sysu;GGCTGTCACGTGGTCAACAT;NGG;GTGCCACGTGGTGAATATTGGAGCA;

#### 4. Run the pipeline

 
- Run the pipeline:
 
   Use the script, CRISPRAnalysis.py, in CRISPRAnalysis/ to run the pipeline.
  
        Usage: CRISPRAnalysis.py [OPTIONS]
        OPTIONS:
             --csv       experiment csv file
             --genome    reference genome fasta
             --output    output folder
             --in1       first fastq data
             --in2       second fastq data
- Pipeline result
 
   The pipeline output result is in the folder, user specified. For example, in outputDir/
 

        CrisprCas9Report.txt
        IndelInfo.txt
        INDEL Figures.

        
    A. CrisprCas9Report.txt
         This file contains the crispr on-target site information and the crispr efficiency.
         For example:
     

        The crispr efficiency assessment of on-target site, csgo
        
        on-target site csgo info:
            Reference Genome : /mnt/c/Users/Jon/Documents/GitHub/crispr-wgs/demo/genome/demo.fa
            Cutsite Postion : chr14: 13724028-13724087
            gRNA : GGCTGTCACGTGGTCAACAT
            PAN : NGG
            HDR sequence : GTGCCACGTGGTGAATATTGGAGCA
        
        Analysis result :
        
        Mutation efficiency : 11.2601%  reads have indels
        Repair efficiency : 1.0724%  reads have be HDR repaired
            4 reads with HDR: GTGCCACGTGGTGAATATTGGAGCA
    
    
    B. IndelInfo.txt
        This file contains the crispr on-target site information and the indels information in the on-target site region.
        For example:
    
    
        The indel distribution analysis of on-target site, csgo
        on-target site csgo info:
            Reference Genome : /mnt/c/Users/Jon/Documents/GitHub/crispr-wgs/demo/genome/demo.fa
            Cutsite Postion : chr14: 13724028-13724087
            gRNA : GGCTGTCACGTGGTCAACAT
            PAN : NGG
            HDR sequence : GTGCCACGTGGTGAATATTGGAGCA
        
        Analysis result :
        Chromosome	Position	Type	Length	FlankingSequence(+)	ReadsNumber	Frequency	TotalReads
        chr14	13724036	Deletion	1	TCAGCCTCTG[C]CATTGGCTGT	2	0.0027	746	
        chr14	13724037	Deletion	21	CAGCCTCTGC[CATTGGCTGTCACGTGGTCAA]CATTGGAGCT	4	0.0054	746	
        chr14	13724046	Deletion	1	CCATTGGCTG[T]CACGTGGTCA	2	0.0027	746	
        chr14	13724048	Deletion	9	ATTGGCTGTC[ACGTGGTCA]ACATTGGAGC	2	0.0027	746	
        chr14	13724051	Deletion	10	GGCTGTCACG[TGGTCAACAT]TGGAGCTTTA	2	0.0027	746	
        chr14	13724052	Deletion	1	GCTGTCACGT[G]GTCAACATTG	2	0.0027	746	
        chr14	13724054	Deletion	6	TGTCACGTGG[TCAACA]TTGGAGCTTT	4	0.0054	746	
        chr14	13724055	Deletion	3	GTCACGTGGT[CAA]CATTGGAGCT	30	0.0402	746	
        chr14	13724056	Insertion	3	TCACGTGGTC[]AACATTGGAG	2	0.0027	746	
        chr14	13724057	Deletion	2	CACGTGGTCA[AC]ATTGGAGCTT	8	0.0107	746	
        chr14	13724057	Insertion	5	CACGTGGTCA[]ACATTGGAGC	2	0.0027	746	
        chr14	13724057	Insertion	7	CACGTGGTCA[]ACATTGGAGC	2	0.0027	746	
        chr14	13724057	Deletion	7	CACGTGGTCA[ACATTGG]AGCTTTAGAC	2	0.0027	746	

   C. INDEL Figures
   ![INDEL Figures](https://github.com/yuht4/CRISPRAnalysis/blob/master/C.png)
   


