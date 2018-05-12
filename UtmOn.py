#! /usr/bin/env python
# coding=utf-8
# --author: yuht4--

import getopt
import os, sys, re, copy, commands
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import pyplot as plt


PREFIX = os.path.split(os.path.realpath(__file__))[0]
PREFIX.replace('\n', '')

ABRA = PREFIX + "/lib/abra2.jar"


def usage():
    print(
        "CRISPR pipeline for assessing the efficiency of modification at the on-target site.\n"
    )
    print("USAGE: python AnalysisOnTarget.py OPTIONS\n")
    print(" Sample: python AnalysisOnTarget.py [OPTIONS]")
    print("     OPTIONS:\n")
    print("     -b    bam file\n")
    print("     -c    csv file\n")
    print("     -g    reference genome fasta\n")
    print("     -o    output folder\n")
    print("\n\n")


def main():

    OutputDir = "/"
    BAMFile = "/"
    RefGenomeFasta = "/"
    CSV = "/"

    CutSiteTupleList = []

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:b:g:c:")
    except getopt.GetoptError:
        print("Command wrong, please view the help information!\n")
        sys.exit()

    if not opts:
        usage()
        sys.exit()

    for o, a in opts:
        if o == "-o":
            OutputDir = a
        elif o == "-b":
            BAMFile = a
        elif o == "-g":
            RefGenomeFasta = a
        elif o == "-c":
            CSV = a
        elif o == "-h":
            usage()
            sys.exit()
        else:
            print("Command wrong, please view the help information!\n")
            sys.exit()

    BAMFile = os.path.abspath(BAMFile)
    OutputDir = os.path.abspath(OutputDir)
    CSV = os.path.abspath(CSV)
    RefGenomeFasta = os.path.abspath(RefGenomeFasta)

    OutputDir = re.sub('\s', '', OutputDir)
    RefGenomeFasta = re.sub('\s', '', RefGenomeFasta)
    BAMFile = re.sub('\s', '', BAMFile)
    CSV = re.sub('\s', '', CSV)

    if not os.path.exists(OutputDir):
        print("The OutputDir not exist! Error\n")
        sys.exit()
    if not os.path.exists(BAMFile):
        print("The bam file not exist! Error\n")
        sys.exit()
    if not os.path.exists(CSV):
        print("The CSV file no exist! Error\n")
        sys.exit()
    if not os.path.exists(RefGenomeFasta):
        print("The Reference fasta not exist! Error\n")
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

        tempTuple = (Chr, Start, Stop, CutSite, GRnaSequence, PAM,
                     RepairSequence)

        if tempTuple not in CutSiteTupleList:
            CutSiteTupleList.append(tempTuple)

    FILE.close()

    if not os.path.exists(OutputDir + "/Temp"):
        os.system("mkdir -p " + OutputDir + "/Temp")

    CutSiteFile = OutputDir + "/Temp/CutSite.bed"
    CutSiteTupleList.sort()

    FILE = open(CutSiteFile, "w+")

    for val in CutSiteTupleList:
        FILE.write(val[0] + "\t" + str(val[1]) + "\t" + str(val[2]) + "\t" +
                   val[3] + "\n")

    FILE.close()

    if not os.path.exists(OutputDir + "/Logs/"):
        os.system("mkdir -p " + OutputDir + "/Logs/")

    AbraLogFile = OutputDir + "/Logs/abra.log"

    if os.path.exists(OutputDir + "/CrisprCas9Report.txt"):
        os.system("rm -rf " + OutputDir + "/CrisprCas9Report.txt")
    if os.path.exists(OutputDir + "/IndelInfo.txt"):
        os.system("rm -rf " + OutputDir + "/IndelInfo.txt")
    if os.path.exists(OutputDir + "/ABRAtemp"):
        os.system("rm -rf " + OutputDir + "/ABRAtemp")
    if os.path.exists(AbraLogFile):
        os.system("rm -rf " + AbraLogFile)

    os.system("samtools view -b -L  " + CutSiteFile + " " + BAMFile + " > " +
              OutputDir + "/Temp/CutSite2.bam")
    os.system("samtools index " + OutputDir + "/Temp/CutSite2.bam")
    os.system(
        "java -jar " + ABRA + " --in " + OutputDir + "/Temp/CutSite2.bam --out "
        + OutputDir + "/Temp/CutSite.bam --ref "  + RefGenomeFasta + " --targets " + CutSiteFile +
        " --threads 4 " + AbraLogFile)
    os.system("samtools index " + OutputDir + "/Temp/CutSite.bam")
    os.system("samtools view -h  " + OutputDir + "/Temp/CutSite.bam" + " -o " +
              OutputDir + "/Temp/CutSite.sam")

    SAMFile = OutputDir + "/Temp/CutSite.sam"

    for val in CutSiteTupleList:
        AnalysisSamFile(OutputDir, SAMFile, val[0], val[1], val[2], val[3],
                        val[4], val[5], val[6], RefGenomeFasta)

    os.system("rm -rf " + OutputDir + "/ABRAtemp")


def AnalysisSamFile(OutputDir, SAMFile, Chr, Start, Stop, CutSite,
                    GRnaSequence, PAM, RepairSequence, RefGenomeFasta):

    CutSiteReadNoIndelGroupsList = []
    CutSiteReadIndelGroupsList = []
    #### 对应read pairname 和 indel

    CutSiteReadRepairGroupsList = []
    ReadNameRepairSequenceDict = {}
    #### 对应 read pairname 和 repair seq的信息

    InDelsInfoTupleList = []
    ### (start_end del length readpairname)

    ConfictReadPairNameList = []
    ### 对应 那些冲突的read pairname

    PossibleRepairSequencesList = []
    ### 对应clean sequ 和 repa seq 的可能

    CleanRepairSequence = RepairSequence

    for ch in ['(', ')', '[', ']']:
        CleanRepairSequence = CleanRepairSequence.replace(ch, '')

    ###print(RepairSequence + '\n');
    if RepairSequence != "":
        PreparationIfHaveRepairSequence(RepairSequence,
                                        PossibleRepairSequencesList)

    FILE = open(SAMFile)

    while True:

        Line = FILE.readline()

        Line = Line.replace('\n', '')

        if not Line:
            break
        if Line[0] == '@':
            continue

        LineValues = Line.split('\t')

        while '' in LineValues:
            LineValues.remove('')

        ReadPairName = LineValues[0]
        ReadChr = LineValues[2]
        ReadStart = int(LineValues[3])
        ReadCigar = LineValues[5]
        ReadSequence = LineValues[9]

        if ReadChr == Chr:
            ## print "fuck ofdfhsdfsdfsdfds";
            FinalCigar = []
            GetFinalCigar(ReadCigar, FinalCigar)
            ReadStop = GetReadStopFromReadStart(ReadStart, ReadCigar)

            ##### sam 文件中的比对 涵盖了 CutSite位点

            ##print Chr + " " + Start + " " + Stop;
            ##print ReadChr + " " +  ReadStart + " " + ReadStop;
            if ReadStart <= Start and ReadStop >= Stop:

                if 'N' in ReadCigar or 'P' in ReadCigar or 'X' in ReadCigar:
                    die()

                if 'D' in ReadCigar or 'I' in ReadCigar:

                    IndelInCutSite = 'no'
                    PostionInRead = ReadStart

                    for i in range(0, len(FinalCigar), 2):

                        if 'D' in FinalCigar[i + 1] or 'I' in FinalCigar[i
                                                                         + 1]:

                            if (PostionInRead < Start and
                                (PostionInRead + FinalCigar[i]) > Start) or (
                                    PostionInRead >= Start
                                    and PostionInRead <= Stop):

                                if (ReadPairName not in
                                        CutSiteReadNoIndelGroupsList) and (
                                            ReadPairName not in
                                            CutSiteReadRepairGroupsList):

                                    if ReadPairName not in CutSiteReadIndelGroupsList:
                                        CutSiteReadIndelGroupsList.append(
                                            ReadPairName)

                                    if FinalCigar[i + 1] == 'D':
                                        DeletionStart = PostionInRead
                                        DeletionEnd = PostionInRead + FinalCigar[i] - 1
                                        Length = FinalCigar[i]

                                        Turple = (DeletionStart, DeletionEnd,
                                                  "Deletion", Length,
                                                  ReadPairName)

                                        if Turple not in InDelsInfoTupleList:
                                            InDelsInfoTupleList.append(Turple)

                                    if FinalCigar[i + 1] == 'I':
                                        InsertionStart = PostionInRead
                                        InsertionEnd = PostionInRead - 1
                                        Length = FinalCigar[i]

                                        Turple = (InsertionStart, InsertionEnd,
                                                  "Insertion", Length,
                                                  ReadPairName)

                                        if Turple not in InDelsInfoTupleList:
                                            InDelsInfoTupleList.append(Turple)

                                    IndelInCutSite = "yes"

                                else:
                                    if ReadPairName not in ConfictReadPairNameList:
                                        ConfictReadPairNameList.append(
                                            ReadPairName)

                        PostionInRead += FinalCigar[i]

                    if IndelInCutSite == 'no':
                        ProcessReadsAsNoIndelOrRepair(
                            RepairSequence, CleanRepairSequence,
                            PossibleRepairSequencesList, ReadPairName,
                            ReadSequence, CutSiteReadIndelGroupsList,
                            CutSiteReadRepairGroupsList,
                            CutSiteReadNoIndelGroupsList,
                            ConfictReadPairNameList,
                            ReadNameRepairSequenceDict)

                else:
                    ProcessReadsAsNoIndelOrRepair(
                        RepairSequence, CleanRepairSequence,
                        PossibleRepairSequencesList, ReadPairName,
                        ReadSequence, CutSiteReadIndelGroupsList,
                        CutSiteReadRepairGroupsList,
                        CutSiteReadNoIndelGroupsList, ConfictReadPairNameList,
                        ReadNameRepairSequenceDict)

    FILE.close()

    #### REMOVE CONFICTS

    for read in ConfictReadPairNameList:
        if read in CutSiteReadIndelGroupsList:
            while read in CutSiteReadIndelGroupsList:
                CutSiteReadIndelGroupsList.remove(read)

        if read in CutSiteReadNoIndelGroupsList:
            while read in CutSiteReadNoIndelGroupsList:
                CutSiteReadNoIndelGroupsList.remove(read)

        if read in CutSiteReadRepairGroupsList:
            while read in CutSiteReadRepairGroupsList:
                CutSiteReadRepairGroupsList.remove(read)
            while ReadNameRepairSequenceDict.has_key(read):
                del ReadNameRepairSequenceDict[read]

        if InDelsInfoTupleList:
            for tup in InDelsInfoTupleList:
                if tup[4] == read:
                    while tup in InDelsInfoTupleList:
                        InDelsInfoTupleList.remove(tup)

    #### Output the File result

    ReportFile = OutputDir + "/CrisprCas9Report.txt"
    FILE = open(ReportFile, "a+")

    FILE.write("The crispr efficiency assessment of on-target site, " +
               CutSite + "\n")
    FILE.write("on-target site " + CutSite + " info:\n")
    FILE.write("    Reference Genome : " + RefGenomeFasta + "\n")
    FILE.write("    Cutsite Postion : " + Chr + ": " + str(Start) + "-" +
               str(Stop) + "\n")
    FILE.write("    gRNA : " + GRnaSequence + "\n")
    FILE.write("    PAN : " + PAM + "\n")

    if RepairSequence:
        FILE.write("    HDR sequence : " + RepairSequence + "\n")

    FILE.write("\nAnalysis result :\n")

    NumberOfIndels = len(CutSiteReadIndelGroupsList)
    NumberOfNondels = len(CutSiteReadNoIndelGroupsList)

    if CutSiteReadRepairGroupsList:
        NumberOfNondels += len(CutSiteReadRepairGroupsList)

    TotalNumberOfReadPairs = NumberOfIndels + NumberOfNondels

    if TotalNumberOfReadPairs != 0:

        EfficiencyRatio = float(NumberOfIndels) / (TotalNumberOfReadPairs)
        FILE.write("Mutation efficiency : " + ('%.4f' % (
            EfficiencyRatio * 100)) + "%  reads have indels\n")
    else:
        FILE.write("No analysis results.\n")

    if len(CutSiteReadRepairGroupsList) != 0:

        RepairsNumber = len(CutSiteReadRepairGroupsList)

        RepairFraction = float(RepairsNumber) / TotalNumberOfReadPairs

        FILE.write("Repair efficiency : " + ('%.4f' % (
            RepairFraction * 100)) + "%  reads have be HDR repaired\n")

        RepairSequenceNoIndelCountDict = {}

        for read in CutSiteReadRepairGroupsList:

            if RepairSequenceNoIndelCountDict.has_key(
                    ReadNameRepairSequenceDict[read]):
                RepairSequenceNoIndelCountDict[ReadNameRepairSequenceDict[
                    read]] += 1
            else:
                RepairSequenceNoIndelCountDict[ReadNameRepairSequenceDict[
                    read]] = 1

        if RepairSequenceNoIndelCountDict.has_key(CleanRepairSequence):

            TempFullHDRNumber = RepairSequenceNoIndelCountDict[
                CleanRepairSequence]

            FILE.write("    " + str(TempFullHDRNumber) + " reads with HDR: " +
                       CleanRepairSequence + "\n")

        for other in RepairSequenceNoIndelCountDict:
            if other != CleanRepairSequence:

                FILE.write("    " + str(RepairSequenceNoIndelCountDict[other])
                           + " reads with HDR: " + other + "\n")
    else:
        FILE.write("No analysis results.\n")

    FILE.write("\n\n")
    FILE.close()

    #############################################
    ###           Write out variants          ###
    #############################################

    IndelInfoFile = OutputDir + "/IndelInfo.txt"

    FILE = open(IndelInfoFile, "a+")

    FILE.write(
        "The indel distribution analysis of on-target site, " + CutSite + "\n")
    FILE.write("on-target site " + CutSite + " info:\n")
    FILE.write("    Reference Genome : " + RefGenomeFasta + "\n")
    FILE.write("    Cutsite Postion : " + Chr + ": " + str(Start) + "-" +
               str(Stop) + "\n")
    FILE.write("    gRNA : " + GRnaSequence + "\n")
    FILE.write("    PAN : " + PAM + "\n")

    if RepairSequence:
        FILE.write("    HDR sequence : " + RepairSequence + "\n")

    FILE.write("\nAnalysis result :\n")

    if len(InDelsInfoTupleList) == 0:
        FILE.write("No analysis results of the cutsite " + CutSite + '.\n')
    else:
        FILE.write(
            "Chromosome\tPosition\tType\tLength\tFlankingSequence(+)\tReadsNumber\tFrequency\tTotalReads\n"
        )

    ### (start, end, type, length, readname)
    InDelsInfoCountDict = {}
    # (start, end, type, length) ==> int

    for each in InDelsInfoTupleList:

        IndelStart = each[0]
        IndelEnd = each[1]
        IndelType = each[2]
        IndelLength = each[3]
        IndelReadName = each[4]
        IndelReadNumber = 1

        if not InDelsInfoCountDict.has_key(
            (IndelStart, IndelEnd, IndelType, IndelLength)):

            for second in InDelsInfoTupleList:
                if second[0] == IndelStart and second[1] == IndelEnd and second[2] == IndelType and second[3] == IndelLength and second[4] != IndelReadName:
                    IndelReadNumber += 1

            InDelsInfoCountDict[(IndelStart, IndelEnd, IndelType,
                                 IndelLength)] = IndelReadNumber

    IndelInfoList = sorted(
        InDelsInfoCountDict.items(), key=lambda item: item[0][0])
    ## [((a, b), c)]

    for each in IndelInfoList:

        IndelStart = each[0][0]
        IndelEnd = each[0][1]
        IndelType = each[0][2]
        IndelLength = each[0][3]
        IndelReadNumber = each[1]

        NumberOfReads = 2 * IndelReadNumber

        TotalReads = 2 * (len(CutSiteReadNoIndelGroupsList) + len(
            CutSiteReadIndelGroupsList) + len(CutSiteReadRepairGroupsList))

        Frequency = float(NumberOfReads) / TotalReads

        tempStart = IndelStart - 10
        tempEnd = IndelEnd + 10

        (tempInt, tempString) = commands.getstatusoutput(
            "samtools faidx " +
            RefGenomeFasta + " " + Chr + ":" + str(tempStart) + "-" + str(tempEnd))

        tempSequence = tempString.split('\n')[1]
        tempSequence = tempSequence.replace('\n', '')
        tempSequence = tempSequence.upper()

        if len(tempSequence) > 20:
            tempSequence = tempSequence[0:10] + '[' + tempSequence[10:(
                len(tempSequence) - 10)] + ']' + tempSequence[(
                    len(tempSequence) - 10):]
        else:
            tempSequence = tempSequence[0:10] + '[]' + tempSequence[10:20]

        FILE.write(Chr + "\t" + str(IndelStart) + "\t" + IndelType + "\t" +
                   str(IndelLength) + "\t" + tempSequence + "\t")
        FILE.write(
            str(NumberOfReads) + "\t" +
            ('%.4f' % Frequency) + "\t" + str(TotalReads) + '\t\n')

    FILE.write('\n\n')
    FILE.close()


    TempList = copy.deepcopy(IndelInfoList);

    InsertionPositionAndCountDict = {};
    DeletionPositionAndCountDict = {};

    InsertionSizeAndCountDict = {};
    DeletionSizeAndCountDict = {};

    small = TempList[0][0][0];


    for each in TempList:

        IndelStart = each[0][0]

        if IndelStart < small:
            small = IndelStart;

    small = small - 1;


    for each in TempList:

        IndelStart = each[0][0]
        IndelEnd = each[0][1]
        IndelType = each[0][2]
        IndelLength = each[0][3]
        IndelReadNumber = each[1]


        if IndelType == "Insertion":

            if not InsertionPositionAndCountDict.has_key(int(IndelStart)- small):

                InsertionPositionAndCountDict[int(IndelStart) - small] = int(IndelReadNumber);
            else:
                InsertionPositionAndCountDict[int(IndelStart) -small] += int(IndelReadNumber);

            if not InsertionSizeAndCountDict.has_key(int(IndelLength)):

                InsertionSizeAndCountDict[int(IndelLength)] = int(IndelReadNumber);
            else:
                InsertionSizeAndCountDict[int(IndelLength)] += int(IndelReadNumber);


        if IndelType == "Deletion":

            if not DeletionPositionAndCountDict.has_key(int(IndelStart) - small):

                DeletionPositionAndCountDict[int(IndelStart) -small ] = int(IndelReadNumber);
            else:
                DeletionPositionAndCountDict[int(IndelStart) - small] += int(IndelReadNumber);

            if not DeletionSizeAndCountDict.has_key(int(IndelLength)):

                DeletionSizeAndCountDict[int(IndelLength)] = int(IndelReadNumber);
            else:
                DeletionSizeAndCountDict[int(IndelLength)] += int(IndelReadNumber);



    sorted(DeletionSizeAndCountDict.keys());
    sorted(InsertionSizeAndCountDict.keys());
    sorted(InsertionPositionAndCountDict.keys());
    sorted(DeletionPositionAndCountDict.keys());



    plt.bar(InsertionSizeAndCountDict.keys(), InsertionSizeAndCountDict.values());

    plt.xlabel('Insertion Size (bp)')
    plt.ylabel('Reads Count')
    plt.title('Insertion Size / Reads Number distribution')

    plt.savefig(OutputDir + "/" + CutSite +"_insertionSize.png");
    plt.clf() # 清图。
    plt.cla() # 清坐标轴。
    plt.close() # 关窗口


    plt.bar(DeletionSizeAndCountDict.keys(), DeletionSizeAndCountDict.values());

    plt.xlabel('Deletion Size (bp)')
    plt.ylabel('Reads Count')
    plt.title('Deletion Size / Reads Number distribution')

    plt.savefig(OutputDir + "/" + CutSite + "_deltionSize.png");

    plt.clf() # 清图。
    plt.cla() # 清坐标轴。
    plt.close() # 关窗口

    plt.bar(InsertionPositionAndCountDict.keys(), InsertionPositionAndCountDict.values());

    plt.xlabel('Insertion Postion (the coord 0 starts from ' + str(small)  + " )");
    plt.ylabel('Reads Count')
    plt.title('Insertion Postion / Reads Number distribution')

    plt.savefig(OutputDir + "/" + CutSite + "_insertionPos.png"); 

    plt.clf() # 清图。
    plt.cla() # 清坐标轴。
    plt.close() # 关窗口

    plt.bar(DeletionPositionAndCountDict.keys(), DeletionPositionAndCountDict.values());

    plt.xlabel('Deletion Postion (the coord 0 starts from ' + str(small)  + " )");
    plt.ylabel('Reads Count')
    plt.title('Deletion Postion / Reads Number distribution')

    plt.savefig(OutputDir + "/" + CutSite +"_deletionPos.png");

    plt.clf() # 清图。
    plt.cla() # 清坐标轴。
    plt.close() # 关窗口

    #print(DeletionPositionAndCountDict);


def PreparationIfHaveRepairSequence(RepairSequence,
                                    PossibleRepairSequencesList):

    LastOpenBracket = 0
    LastClosedBracket = 0
    NucleotidesList = ["A", "C", "G", "T"]
    FinalSequencesList = []

    for i in range(len(RepairSequence)):

        if RepairSequence[i] == '(':

            ChangeString = ""
            ChangeNucleArray = []
            LastOpenBracket = i

            i += 1

            while RepairSequence[i] != ')':
                ChangeString += RepairSequence[i]
                i += 1
            LastClosedBracket = i

            ChangeNucleArray = copy.deepcopy(NucleotidesList)

            if (len(ChangeString) != 1):
                for key in range(1, len(ChangeString)):
                    NewChangeNucleArray = []

                    for change in ChangeNucleArray:
                        for nucl in NucleotidesList:
                            NewChangeNucleArray.append(change + nucl)
                    ChangeNucleArray = copy.deepcopy(NewChangeNucleArray)

            if FinalSequencesList:
                NewFinalSequencesList = []

                for FinalChange in FinalSequencesList:
                    InterChangeSequence = RepairSequence[len(FinalChange):
                                                         LastOpenBracket]

                    for mmchange in ChangeNucleArray:
                        NewFinalSeq = FinalChange + InterChangeSequence + '(' + mmchange + ')'
                        NewFinalSequencesList.append(NewFinalSeq)
                FinalSequencesList = copy.deepcopy(NewFinalSequencesList)
            else:
                for myChange in ChangeNucleArray:
                    FianlSeq = RepairSequence[0:
                                              LastOpenBracket] + '(' + myChange + ')'
                    FinalSequencesList.append(FianlSeq)

    for FinalSeqS in FinalSequencesList:

        for ch in ['(', ')', '[', ']']:
            FinalSeqS = FinalSeqS.replace(ch, '')

        FinalSeqS += RepairSequence[(LastClosedBracket + 1):]

        for ch in ['(', ')', '[', ']']:
            FinalSeqS = FinalSeqS.replace(ch, '')

        PossibleRepairSequencesList.append(FinalSeqS)

    if '(' not in RepairSequence and ')' not in RepairSequence:

        CleanRepairSequence = RepairSequence

        for ch in ['(', ')', '[', ']']:
            CleanRepairSequence = CleanRepairSequence.replace(ch, '')

        PossibleRepairSequencesList.append(CleanRepairSequence)

    PossibleRepairSequencesList = list(set(PossibleRepairSequencesList))


def GetReadStopFromReadStart(ReadStart, ReadCigar):

    ReadStop = ReadStart - 1
    FinalCigar = []
    GetFinalCigar(ReadCigar, FinalCigar)

    for i in range(0, len(FinalCigar), 2):
        if FinalCigar[i + 1] == 'M' or FinalCigar[i + 1] == 'D':
            ReadStop += int(FinalCigar[i])

    return ReadStop


def GetFinalCigar(ReadCigar, FinalCigar):

    NumberOfBases = ""
    for i in range(len(ReadCigar)):
        if ReadCigar[i].isalpha():
            FinalCigar.append(int(NumberOfBases))
            FinalCigar.append(ReadCigar[i])
            NumberOfBases = ""
        else:
            NumberOfBases += ReadCigar[i]


def ProcessReadsAsNoIndelOrRepair(
        RepairSequence, CleanRepairSequence, PossibleRepairSequencesList,
        ReadPairName, ReadSequence, CutSiteReadIndelGroupsList,
        CutSiteReadRepairGroupsList, CutSiteReadNoIndelGroupsList,
        ConfictReadPairNameList, ReadNameRepairSequenceDict):

    if RepairSequence != "":
        TempRepairSequence = CleanRepairSequence
        TempReapairSequencesList = []

        for seq in PossibleRepairSequencesList:

            seqRC = seq[::-1]

            tempDict = {
                'A': 'T',
                'a': 't',
                'T': 'A',
                't': 'a',
                'C': 'G',
                'c': 'g',
                'G': 'C',
                'g': 'c'
            }
            tempList = list(seqRC)

            for i in range(len(tempList)):
                tempList[i] = tempDict[tempList[i]]

            seqRC = "".join(tempList)

            if seq in ReadSequence or seqRC in ReadSequence:
                if seq not in TempReapairSequencesList:
                    TempReapairSequencesList.append(seq)

        if (ReadPairName not in CutSiteReadIndelGroupsList) and (
                ReadPairName not in CutSiteReadRepairGroupsList
        ) and len(TempReapairSequencesList) == 0:
            if ReadPairName not in CutSiteReadNoIndelGroupsList:
                CutSiteReadNoIndelGroupsList.append(ReadPairName)

        elif (ReadPairName not in CutSiteReadIndelGroupsList) and (
                ReadPairName not in CutSiteReadNoIndelGroupsList
        ) and len(TempReapairSequencesList) == 1:
            if ReadPairName not in CutSiteReadRepairGroupsList:
                CutSiteReadRepairGroupsList.append(ReadPairName)

            ReadNameRepairSequenceDict[
                ReadPairName] = TempReapairSequencesList[0]
        else:
            if ReadPairName not in ConfictReadPairNameList:
                ConfictReadPairNameList.append(ReadPairName)

    else:
        if ReadPairName not in CutSiteReadIndelGroupsList:
            if ReadPairName not in CutSiteReadNoIndelGroupsList:
                CutSiteReadNoIndelGroupsList.append(ReadPairName)

        else:
            if ReadPairName not in ConfictReadPairNameList:
                ConfictReadPairNameList.append(ReadPairName)


if __name__ == "__main__":
    main()
