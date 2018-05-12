#! /usr/bin/env python
# coding=utf-8
# --author: yuht4--

import os, sys, copy, getopt, re, argparse

PREFIX = os.path.split(os.path.realpath(__file__))[0]
PREFIX.replace('\n', '')

PICARD = PREFIX + "/lib/picard.jar"
STRELKA = PREFIX + "/lib/strelka/bin/"
BBMAP = PREFIX + "/lib/bbmap/"

def main():

    OutputDir = "/"
    RefGenomeFasta = "/"
    CSV = "/"
    CutSiteTupleList = []

    parser = argparse.ArgumentParser(description="crispr NGS data analysis")
    parser.add_argument("--output", type=str, help="output folder", required=True)
    parser.add_argument("--genome", type=str, help="reference fasta", required=True)
    parser.add_argument("--in1", type=str, help="first fastq file", required=True)
    parser.add_argument("--in2", type=str, help="second fastq file", required=True)
    parser.add_argument("--csv", type=str, help="csv file", required=True)
    args = parser.parse_args()


    FirstFastq = os.path.abspath(args.in1)
    SecondFastq = os.path.abspath(args.in2)
    OutputDir = os.path.abspath(args.output)
    RefGenomeFasta = os.path.abspath(args.genome)
    CSV = os.path.abspath(args.csv)

    FirstFastq = re.sub('\s', '', FirstFastq)
    SecondFastq = re.sub('\s', '', SecondFastq)
    OutputDir = re.sub('\s', '', OutputDir)
    CSV = re.sub('\s', '', CSV)
    RefGenomeFasta = re.sub('\s', '', RefGenomeFasta)


    #print(OutputDir);
    #print(os.path.exists(OutputDir));

    if not os.path.exists(OutputDir):
        print("The OutputDir not exist! Error\n")
        sys.exit()
    if not os.path.exists(FirstFastq) or not os.path.exists(SecondFastq):
        print("The sequencing data not exist! Error\n")
        sys.exit()
    if not os.path.exists(RefGenomeFasta):
        print("The Reference fasta not exist! Error\n")
        sys.exit()
    if not os.path.exists(CSV):
        print("The CSV file not exist! Error\n")
        sys.exit()

    FILE = open(CSV, 'r')

    while True:

        Line = FILE.readline()
        Line = re.sub('\n', '', Line)
        Line = re.sub('\r', '', Line)

        if not Line:
            break

        LineValues = Line.split(';')

        while '' in LineValues:
            LineValues.remove('')

        Chr = LineValues[0]
        Start = int(LineValues[1])
        Stop = int(LineValues[2])
        CutSite = LineValues[3]
        GRnaSequence = LineValues[4]
        PAM = LineValues[5]

        RepairSequence = ""

        if len(LineValues) == 7:
            RepairSequence = LineValues[6]

        tempTuple = (Chr, Start, Stop, CutSite, GRnaSequence, PAM, RepairSequence)

        if tempTuple not in CutSiteTupleList:
            CutSiteTupleList.append(tempTuple)

    FILE.close()

    if not os.path.exists(OutputDir + "/Logs/"):
        os.system("mkdir -p " + OutputDir + "/Logs/")
    else:
        os.system("rm -rf " + PREFIX + "/Logs/*")

    if not os.path.exists(OutputDir + "/Temp"):
        os.system("mkdir -p " + OutputDir + "/Temp")
    else:
        os.system("rm -rf " + OutputDir + "/Temp/*")

    if os.path.exists(OutputDir + "/Result"):
        os.system("rm -rf " + OutputDir + "/Result")

    if os.path.exists(OutputDir + "/cutsiteinfo.csv"):
        os.system("rm -rf " + OutputDir + "/cutsiteinfo.csv")

    if os.path.exists(OutputDir + "/CrisprCas9Report.txt"):
        os.system("rm -rf " + OutputDir + "/CrisprCas9Report.txt")

    if os.path.exists(OutputDir + "/IndelInfo.txt"):
        os.system("rm -rf " + OutputDir + "/IndelInfo.txt")

    if os.path.exists(OutputDir + "/OffTargetMotif.txt"):
        os.system("rm -rf " + OutputDir + "/OffTargetMotif.txt")

    if os.path.exists(OutputDir + "/OffTargetSites.txt"):
        os.system("rm -rf " + OutputDir + "/OffTargetSites.txt")

    if os.path.exists(OutputDir + "/variants.vcf"):
        os.system("rm -rf " + OutputDir + "/variants.vcf")


    FliterReads(FirstFastq, SecondFastq, OutputDir)
    MappingAndSortingReads(OutputDir, RefGenomeFasta)

    FILE = open(OutputDir + "/cutsiteinfo.csv", 'w+')
    for val in CutSiteTupleList:
        FILE.write(val[5] + ";" + val[6] + ";\n")
    FILE.close()


    os.system(
        "python " + PREFIX + "/UtmOn.py -o " + OutputDir + " -g " + RefGenomeFasta
        + " -c " + CSV + " -b " + OutputDir + "/Temp/WGS.rg.sort.dd.bam")

    removeTempFile(OutputDir)


def removeTempFile(OutputDir):

    if os.path.exists(OutputDir + "/Temp"):
        os.system("rm -rf " + OutputDir + "/Temp")
    if os.path.exists(OutputDir + "/Result"):
        os.system("rm -rf " + OutputDir + "/Result")
    if os.path.exists(OutputDir + "/cutsiteinfo.csv"):
        os.system("rm -rf " + OutputDir + "/cutsiteinfo.csv")
    if os.path.exists(OutputDir + "/Logs"):
        os.system("rm -rf " + OutputDir + "/Logs")

def FliterReads(FirstFastq, SecondFastq, OutputDir):

    FastXLogFile = OutputDir + "/Logs/fastx.log"
    BBMapLogFile = OutputDir + "/Logs/bbmap.log"

    os.system(
        "gunzip -c " + FirstFastq + " > " + OutputDir + "/Temp/first.fastq")
    os.system(
        "gunzip -c " + SecondFastq + " > " + OutputDir + "/Temp/second.fastq")

    os.system("fastq_quality_trimmer -Q33 -t 30 -z -i " + OutputDir +
              "/Temp/first.fastq  -o " + OutputDir +
              "/Temp/first.fliter.fastq.gz > " + FastXLogFile)
    os.system("fastq_quality_trimmer -Q33 -t 30 -z -i " + OutputDir +
              "/Temp/second.fastq  -o " + OutputDir +
              "/Temp/second.fliter.fastq.gz >> " + FastXLogFile)

    os.system("rm -rf " + OutputDir + "/Temp/first.fastq")
    os.system("rm -rf " + OutputDir + "/Temp/second.fastq")

    os.system(
        BBMAP + "/repair.sh --overwrite=t in1=" + OutputDir +
        "/Temp/first.fliter.fastq.gz in2=" + OutputDir +
        "/Temp/second.fliter.fastq.gz out1=" + OutputDir +
        "/Temp/first.fliter.repair.fastq.gz out2=" + OutputDir +
        "/Temp/second.fliter.repair.fastq.gz >> " + BBMapLogFile + " 2>&1")


def MappingAndSortingReads(OutputDir, RefGenomeFasta):

    BwaLogFile = OutputDir + "/Logs/bwa.log"
    SamLogFile = OutputDir + "/Logs/sam.log"
    PicardLogFile = OutputDir + "/Logs/picard.log"
    StrelkaLogFile = OutputDir + "/Logs/strelka.log"

    os.system(
        "bwa mem -t 4 -M " + RefGenomeFasta +
        " " + OutputDir + "/Temp/first.fliter.repair.fastq.gz " +
        OutputDir + "/Temp/second.fliter.repair.fastq.gz  > " + OutputDir +
        "/Temp/WGS.sam 2>> " + BwaLogFile)
    os.system("samtools view -Sb  " + OutputDir + "/Temp/WGS.sam > " +
              OutputDir + "/Temp/WGS.bam 2>> " + SamLogFile)

    os.system(
        "java -jar " + PICARD + " SortSam I=" + OutputDir + "/Temp/WGS.bam O="
        + OutputDir + "/Temp/WGS.sort.bam SO=coordinate 2>> " + PicardLogFile)
    os.system(
        "java -jar " + PICARD + " MarkDuplicates I=" + OutputDir +
        "/Temp/WGS.sort.bam O=" + OutputDir + "/Temp/WGS.rg.sort.bam M=" +
        OutputDir + "/Temp/WGS.rg.sort.dd.metrics REMOVE_DUPLICATES=true 2>> "
        + PicardLogFile)
    os.system(
        "java -jar " + PICARD + " AddOrReplaceReadGroups I=" + OutputDir +
        "/Temp/WGS.rg.sort.bam O=" + OutputDir +
        "/Temp/WGS.rg.sort.dd.bam SO=coordinate LB=hg19ID PL=Illumina PU=hg19PU SM=1 2>> "
        + PicardLogFile)

    os.system("samtools index " + OutputDir + "/Temp/WGS.rg.sort.dd.bam")

if __name__ == "__main__":
    main()
