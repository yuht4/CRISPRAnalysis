#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
stand alone filter script based on jgi-rqc-pipeline/filter/fastq_filter2.py

Command
    $ filter.py -f FASTQ.GZ -o OUT_DIR --skip-blast -html

Outputs
  - normal filter outputs + index.html

Created: March 15, 2018

Shijie Yao (syao@lbl.gov)

Revision:

"""

import os
import sys
import argparse
import datetime
import shutil
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg") ## This needs to skip the DISPLAY env var checking
import matplotlib.pyplot as plt
import mpld3

## append the pipeline lib and tools relative path:
SRC_ROOT = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SRC_ROOT + "/lib")   # common

import rqc_fastq as fastqUtil
from common import get_logger, get_status, checkpoint_step, set_colors, get_subsample_rate, append_rqc_file, append_rqc_stats
from os_utility import make_dir, change_mod, run_sh_command, get_tool_path, make_dir_p
from rqc_utility import safe_basename
from db_access import jgi_connect_db
from rqc_constants import RQCReferenceDatabases
from fastq_filter2_constants import *

from rqc_utility import get_dict_obj, pipeline_val
from readqc import do_html_body


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## global vars
VERSION = "3.2.8"
SUCCESS = 0
FAILURE = -1
SCRIPT_NAME = __file__
logLevel = "DEBUG"
# logLevel = "INFO"

## the following constants needs work - some central location ?
# TOOLS_BIN = SRC_ROOT + '/tools'

DEBUG = False
color = {}
color = set_colors(color, True)
NERSC_HOST = os.environ['NERSC_HOST']

READQC_DIR = 'filtered-readqc'

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pipeline function definitions

def log_and_print(msg):
    print(msg)
    log.info(msg)


"""
Run rqcfilter.sh

@param fastq: the raw fastq input file [input of filtering]
@param outDir: the output directory
@param prodType: product type
@param status: current status
@return outFastqFile, outRrnaFastqFile

"""
def run_rqcfilter(infastq, outDir, prodType, status, enableRmoveMicrobes, enableAggressive, disableRmoveMicrobes, disableClumpify, taxList, log):
    log_and_print("\n\n%s - RUN RQCFILTER <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))

    bbtoolsRqcfilterSh = get_tool_path("rqcfilter.sh", "bbtools")

    make_dir_p(outDir)
    opt = None

    extraOptions = ""
    optionFile = ""

    if prodType.endswith("-OLD"):
        extraOptions = " barcodefilter=f chastityfilter=f"
        prodType = prodType.replace("-OLD", "")

    if prodType in FILTER_METHODS_TXT:
        optionFile = os.path.join(SRC_ROOT, "filter_param/", FILTER_METHODS_TXT[prodType].replace(".txt", ".config"))
        log_and_print("Filter option file: %s" % optionFile)
        opt = open(optionFile, 'r').readline().rstrip()
    else:
        log_and_print("The product type, %s, is not supported yet." % prodType)
        sys.exit(2)

    assert (opt), "Null filter options."

    if prodType in ("METATRANSCRIPTOME", "MTF"):
        if infastq.endswith(".gz"):
            opt += " outribo=%s " % (os.path.basename(infastq).replace(".fastq", ".rRNA.fastq"))
        else:
            opt += " outribo=%s " % (os.path.basename(infastq).replace(".fastq", ".rRNA.fastq.gz"))

    if enableRmoveMicrobes:
        if opt.find("removemicrobes=f") != -1:
            opt = opt.replace("removemicrobes=f", "removemicrobes=t")
        opt += " removemicrobes=t "

    if disableRmoveMicrobes:
        if opt.find("removemicrobes=t") != -1:
            opt = opt.replace("removemicrobes=t", "removemicrobes=f")
        else: opt += " removemicrobes=f "

    if enableAggressive:
        opt += " aggressive=t microbebuild=3 "

    if taxList:
        opt += " taxlist=%s " % (taxList)

    ## Temp set clumpify=t for all prod types (RQC-890)
    if not disableClumpify:
        opt += " clumpify=t "
    else:
        opt += " clumpify=f "

    ## temporary
    if NERSC_HOST == "cori":
        opt += " tmpdir=null "

    opt += extraOptions


    filterLogFile = os.path.join(outDir, "filter.log")
    cmdStr = "%s in=%s path=%s %s usejni=f > %s 2>&1" % (bbtoolsRqcfilterSh, infastq, outDir, opt, filterLogFile)

    rtn = [None, None, None, None, None, None, None, None]
    outFastqFile = None
    outRrnaFastqFile = None
    outUnknownFile = None
    outFragFile = None
    outSingletonFile = None
    filterCode = None
    fileNamePrefix = None


    shFileName = "%s/filter.sh" % outDir

    def find_filtered_fastq(adir):
        outFastqFile = None
        outRrnaFastqFile = None
        outUnknownFile = None
        outFragFile = None
        outSingletonFile = None
        filterCode = None
        fileNamePrefix = None

        searchPatt = os.path.join(adir, "*.fastq.gz")
        outFastqFileList = glob.glob(searchPatt)

        assert len(outFastqFileList) >= 1, "ERROR: cannot find *.fastq.gz output file."
        for f in outFastqFileList:
            f = safe_basename(f, log)[0]
            t = f.split(".")

            if t[-3] not in ("frag", "singleton", "unknown", "rRNA", "lmp"):
                filterCode = t[-3]

            elif t[-3] == "lmp": ## nextera
                filterCode = t[-4]

            if len(t) == 7: ## ex) 12345.1.1234.ACCCC.anqdpht.fastq.gz
                fileNamePrefix = '.'.join(t[:4])

            elif len(t) == 6: ## ex) 6176.5.39297.anqrpht.fastq.gz
                fileNamePrefix = '.'.join(t[:3])

            else:
                log.warning("Unexpected filtered file name, %s", outFastqFileList)
                fileNamePrefix = '.'.join(t[:-3])
                log_and_print("Use %s as file prefix." %fileNamePrefix)
                #return return_values(rtn, False)

        assert filterCode and fileNamePrefix, "ERROR: unexpected filter file name: %s" % (outFastqFileList)

        of = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "fastq.gz"]))
        rof = os.path.join(adir, '.'.join([fileNamePrefix, "rRNA.fastq.gz"]))
        lof = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "lmp.fastq.gz"]))
        uof = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "unknown.fastq.gz"]))
        fof = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "frag.fastq.gz"]))
        sof = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "singleton.fastq.gz"]))
        rpf = os.path.join(adir, "reproduce.sh")
        cof = os.path.join(adir, os.path.basename(infastq).replace(".fastq", "").replace(".gz", "") + ".filter_cmd-%s.sh" % (prodType)) # "(seq_unit_name_prefix).filter_cmd-(prodType).sh"

        if os.path.isfile(of):
            outFastqFile = of
            if os.path.isfile(rof):
                outRrnaFastqFile = rof

        elif os.path.isfile(lof):
            outFastqFile = lof
            if os.path.isfile(uof):
                outUnknownFile = uof
            if os.path.isfile(fof):
                outFragFile = fof
            if os.path.isfile(sof):
                outSingletonFile = sof
        else:
            log.error("Cannot find fastq.gz file.")
            return return_values(rtn, False)

        return outFastqFile, rpf, outRrnaFastqFile, outFragFile, outSingletonFile, outSingletonFile, cof

    def find_filter_number(outFastqFile, outRrnaFastqFile):
        filteredReadNumRrna = 0
        filteredReadNum = fastqUtil.check_fastq_format(outFastqFile) / 4

        if outRrnaFastqFile:
            filteredReadNumRrna = fastqUtil.check_fastq_format(outRrnaFastqFile) / 4

        if filteredReadNum < 0:
            log_and_print("RUN RQCFILTER - filtered fastq format error: %s." % outFastqFile)
            return return_values(rtn, False)

        elif filteredReadNumRrna < 0:
            log_and_print("RUN RQCFILTER - rRna filtered fastq format error: %s" % outRrnaFastqFile)
            return return_values(rtn, False)

        return filteredReadNum, filteredReadNumRrna


    if STEP_ORDER[status] < STEP_ORDER[RQCFILTER_END]:

        create_shell(shFileName, (cmdStr,))

        log_and_print("rqcfilter cmd=[%s]" % cmdStr)
        log_and_print("sh file name=[%s]" % shFileName)

        stdOut, stdErr, exitCode = run_sh_command(shFileName, True, log, True) ## stdOut of 0 is success

        if exitCode != 0:
            log.error("Failed to run : %s, stdout : %s, stderr: %s", shFileName, stdOut, stdErr)
            return return_values(rtn, False)

        outFastqFile, rpf, outRrnaFastqFile, outFragFile, outSingletonFile, outSingletonFile, cof = find_filtered_fastq(outDir)


        filteredReadNum, filteredReadNumRrna = find_filter_number(outFastqFile, outRrnaFastqFile)
        log_and_print("Read counts after RQCFILTER step = %d" % filteredReadNum)

        if outRrnaFastqFile:
            log_and_print("rRNA read counts after RQCFILTER step = %d" % filteredReadNumRrna)

        if os.path.isfile(rpf):
            with open(cof, 'w') as outFH:
                with open(rpf, 'r') as rpfFH:
                    for l in rpfFH:
                        w = l.strip().split()
                        newLine = ""
                        for i in w:
                            if i.find("=") != -1:
                                t = i.split("=")
                                if t[1].find("/") != -1:
                                    t[1], exitCode = safe_basename(t[1], log)
                                newLine += '='.join(t) + ' '
                            elif i.startswith("/"):
                                i, exitCode = safe_basename(i, log)
                                newLine += i + ' '
                            else:
                                newLine += i + ' '
                        newLine += '\n'
                        outFH.write(newLine + '\n')


        else:
            log.error("Cannot find reproduce.sh file.")
            return return_values(rtn, False)

        log_and_print("RUN RQCFILTER - completed")
        checkpoint(RQCFILTER_END, status)
        status = RQCFILTER_END

        if filteredReadNum == 0:
            log.warning("No reads left after filtering")
            checkpoint(PIPE_COMPLETE, status)
            with open(STATS_LIST_FILE_NAME, 'a') as fh:
                write_stats(fh, FILTER_READ_COUNT, 0, log)
                write_stats(fh, FILTER_READ_BASE_COUNT, 0, log)
    else:
        log_and_print("No need to rerun RQCFILTER step, get filtered files and stats ... ")

        ## DEBUG : TODO uncomments !
        outFastqFile, rpf, outRrnaFastqFile, outFragFile, outSingletonFile, outSingletonFile, cof = find_filtered_fastq(outDir)
        # outFastqFile='/global/projectb/scratch/syao/standalone/filter_large/12248.8.247376.CAGAGTG-ACACTCT.anqrpht.fastq.gz'
        # rpf='/global/projectb/scratch/syao/standalone/filter_large/reproduce.sh'
        # outRrnaFastqFile=None
        # outFragFile=None
        # outSingletonFile=None
        # outSingletonFile=None
        # cof='/global/projectb/scratch/syao/standalone/filter_large/12248.8.247376.CAGAGTG-ACACTCT.filter_cmd-RNA.sh'

        filteredReadNum, filteredReadNumRrna = find_filter_number(outFastqFile, outRrnaFastqFile)
        # filteredReadNum=8048702
        # filteredReadNumRrna=0
        # print('filteredReadNum=%s' % filteredReadNum)
        # print('filteredReadNumRrna=%s' % filteredReadNumRrna)

    rtn = [outFastqFile, outRrnaFastqFile, outFragFile, outSingletonFile, outUnknownFile, filteredReadNum, filteredReadNumRrna, cof, status]
    return return_values(rtn)


"""
gzip the final filtered fastq file, and generate the FILE_LIST_FILE_NAME and STATS_LIST_FILE_NAME log files.

@param fastq: the raw fastq file
@param outDir: output dir
@param filteredFastq: the filtered fastq file
@param compressFlag: True = Gzip filtered fastq, False = leave uncompressed
@param status: where the pipeline was at by last run
@param log

@return filteredFastq or None if error
"""
def post_process(fastq, outDir, filteredFastq, rrnaFilteredFastq, fragFile, singletonFile, unknownFile, prodType, filteredReadNum, filteredReadNumRrna, cof, status, log):
    ## obtain read counts from input and filtered fastq files and save the values to STATS_LIST_FILE_NAME file;
    ## compress the filtered fastq file
    log_and_print("\n\n%s - RUN POST PROCESS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    if STEP_ORDER[status] < STEP_ORDER[POST_END]:
        checkpoint(POST_START, status)
        rawCnt = 0
        newCnt = 0
        newCntRrna = 0

        ## Get input read counts
        ## This is redundant but just make sure the input read num
        ## TODO: use "inputReads" from filterStats.txt
        log_and_print("Get read counts from input fastq file [%s]" % fastq)
        rawCnt = fastqUtil.read_count(fastq)     ## read count from the original fastq
        readCounts = {'rawCnt': rawCnt}
        log_and_print("Read counts in input fastq file = %d" % rawCnt)

        if filteredFastq and os.path.isfile(filteredFastq):
            log_and_print("Get read counts from filtered fastq file [%s]" % filteredFastq)
            if filteredReadNum >= 0:
                newCnt = filteredReadNum
            else:
                newCnt = fastqUtil.read_count(filteredFastq)
            log_and_print("Read counts in final filtered fastq file = %d" % newCnt)
        else:
            log.error("Filtered Fastq file does not exist [%s]" % filteredFastq)
            return None

        if prodType in ("METATRANSCRIPTOME", "MTF") and rrnaFilteredFastq and os.path.isfile(rrnaFilteredFastq):
            log_and_print("Get read counts from filtered fastq file [%s]" % rrnaFilteredFastq)
            if filteredReadNumRrna >= 0:
                newCntRrna = filteredReadNumRrna
            else:
                newCntRrna = fastqUtil.read_count(rrnaFilteredFastq)

            log_and_print("Read counts in final rRNA filtered fastq file = %d" % newCntRrna)


        readCounts['newCnt'] = newCnt
        readCounts['newCntRrna'] = newCntRrna


        ###################################################################
        log_and_print(" - write to file-list.txt file, %s" % FILE_LIST_FILE_NAME)
        ###################################################################
        if os.path.isfile(FILE_LIST_FILE_NAME):
            os.remove(FILE_LIST_FILE_NAME)

        refStats = {}
        cardinality = None
        bbdukVersion = None
        bbmapVersion = None

        with open(FILE_LIST_FILE_NAME, 'a') as fh:
            write_files(fh, "filtered_fastq_old", filteredFastq, log)
            if rrnaFilteredFastq:
                write_files(fh, "rrna_filtered_fastq", rrnaFilteredFastq, log)
            if fragFile:
                write_files(fh, "frag_filtered_fastq", fragFile, log)
            if singletonFile:
                write_files(fh, "singleton_filtered_fastq", singletonFile, log)
            if unknownFile:
                write_files(fh, "unknown_filtered_fastq", unknownFile, log)

            ## Write stat files generated from rqcfilter.sh to file-list.txt file
            if os.path.isfile("kmerStats.txt"):
                write_files(fh, "kmerStats", os.path.join(outDir, "kmerStats.txt"), log)
            if os.path.isfile("kmerStats1.txt"):
                write_files(fh, "kmerStats1", os.path.join(outDir, "kmerStats1.txt"), log)
            if os.path.isfile("kmerStats2.txt"):
                write_files(fh, "kmerStats2", os.path.join(outDir, "kmerStats2.txt"), log)

            if os.path.isfile("scaffoldStats.txt"):
                write_files(fh, "scaffoldStats", os.path.join(outDir, "scaffoldStats.txt"), log)
            if os.path.isfile("scaffoldStats1.txt"):
                write_files(fh, "scaffoldStats1", os.path.join(outDir, "scaffoldStats1.txt"), log)
            if os.path.isfile("scaffoldStats2.txt"):
                write_files(fh, "scaffoldStats2", os.path.join(outDir, "scaffoldStats2.txt"), log)

            if os.path.isfile("filterStats.txt"):
                write_files(fh, "filterStats", os.path.join(outDir, "filterStats.txt"), log)
            if os.path.isfile("nexteraStats.txt"):
                write_files(fh, "nexteraStats", os.path.join(outDir, "nexteraStats.txt"), log)
            if os.path.isfile(cof):
                write_files(fh, "filterCmdShFile", cof, log)

            ## With mtst enabled for MTA prod
            # if os.path.isfile("spikein.fq.gz"):
                # write_files(fh, "spikein.fq.gz", os.path.join(outDir, "spikein.fq.gz"), log)
            if os.path.isfile("scaffoldStatsSpikein.txt"):
                write_files(fh, "scaffoldStatsSpikein.txt", os.path.join(outDir, "scaffoldStatsSpikein.txt"), log)

            if os.path.isfile("khist.txt"):
                kHistDataFile = os.path.join(outDir, "khist.txt")
                write_files(fh, FILTER_READ_KHIST_TEXT, kHistDataFile, log)
                kHistPngFile, kHistHtmlPlotFile = gen_read_khist_plot(kHistDataFile, log)
                log.debug("kmer hist plotting: %s, %s, %s", kHistDataFile, kHistPngFile, kHistHtmlPlotFile)
                write_files(fh, FILTER_READ_KHIST_PLOT, kHistPngFile, log)
                write_files(fh, FILTER_READ_KHIST_D3_HTML_PLOT, kHistHtmlPlotFile, log)

            if os.path.isfile("filter.log"):
                write_files(fh, "filterLog", os.path.join(outDir, "filter.log"), log)

                filterLogStat = {}

                with open(os.path.join(outDir, "filter.log"), "r") as FLFH:
                    ## filter.log
                    ##
                    ## ...
                    ## Input:                   218162374 reads         22034399774 bases.
                    ## FTrimmed:                218162374 reads (100.00%)   218162374 bases (0.99%)
                    ## KTrimmed:                9548466 reads (4.38%)   923023698 bases (4.19%)
                    ## Trimmed by overlap:      132062 reads (0.06%)    2684566 bases (0.01%)
                    ## Result:                  209021330 reads (95.81%)    20890529136 bases (94.81%)
                    ## Unique 31-mers:          966117157
                    ## ...
                    ## Input:                   209021330 reads         20890529136 bases.
                    ## Contaminants:            1661712 reads (0.79%)   165909728 bases (0.79%)
                    ## QTrimmed:                7871517 reads (3.77%)   161911788 bases (0.78%)
                    ## Low quality discards:    171448 reads (0.08%)    16752146 bases (0.08%)
                    ## Result:                  206220610 reads (98.66%)    20545955474 bases (98.35%)
                    ## Unique 31-mers:          789287410
                    ## ...

                    isContamNumChecked = False ## Contamination will be done twice for removeribo or for MTF
                    isKtrimmedTotalRemovedNumChecked = False ## for parsing "Total Removed" after ktrimming

                    for l in FLFH:
                        if l.startswith("Input:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 2
                            if 'adaptertriminput' not in filterLogStat:
                                filterLogStat["adaptertriminput"] = {"numreads": toks[0], "numbases": toks[1]}
                            elif 'contamtriminput' not in filterLogStat:
                                filterLogStat["contamtriminput"] = {"numreads": toks[0], "numbases": toks[1]}

                        elif l.startswith("FTrimmed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["ftrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        elif l.startswith("KTrimmed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["ktrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                            isKtrimmedTotalRemovedNumChecked = True

                        ## RQCSUPPORT-1987
                        elif l.startswith("Total Removed:") and isKtrimmedTotalRemovedNumChecked:
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["ktrimmed_total_removed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                            isKtrimmedTotalRemovedNumChecked = False

                        elif l.startswith("Trimmed by overlap:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["trimmedbyoverlap"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        elif l.startswith("Result:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            if 'adaptertrimresult' not in filterLogStat:
                                filterLogStat["adaptertrimresult"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                            elif 'contamtrimresult' not in filterLogStat:
                                filterLogStat["contamtrimresult"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        elif l.startswith("Unique 31-mers:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 2 or len(toks) == 1
                            if 'adaptertrimunique31mers' not in filterLogStat:
                                if len(toks) == 2:
                                    filterLogStat["adaptertrimunique31mers"] = {"num": toks[1]}
                                else:
                                    filterLogStat["adaptertrimunique31mers"] = {"num":"0"}
                            else:
                                if len(toks) == 2:
                                    filterLogStat["contamtrimunique31mers"] = {"num": toks[1]}
                                else:
                                    filterLogStat["contamtrimunique31mers"] = {"num":"0"}

                        elif not isContamNumChecked and l.startswith("Contaminants:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["contaminants"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                            isContamNumChecked = True

                        elif l.startswith("QTrimmed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["qtrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        elif l.startswith("Short Read Discards:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["shortreaddiscards"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        elif l.startswith("Low quality discards:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["lowqualitydiscards"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        elif l.startswith("BBDuk version"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 1
                            bbdukVersion = toks[0]
                        elif l.startswith("BBMap version"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 1
                            bbmapVersion = toks[0]

                        ## BBDuk 36.12 06272016
                        elif l.startswith("Adapter Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["adaptersequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        elif l.startswith("Synthetic Contam Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["syntheticcontamsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        ## 08112016
                        elif l.startswith("Short Synthetic Contam Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["shortsyntheticcontamsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        elif l.startswith("Ribosomal Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["ribosomalsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        ## BBMap 36.12 06272016
                        elif l.startswith("Human Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["humansequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                        ## RQC-862, RQC-880
                        elif l.startswith("Microbial Sequence Removed:"):
                            toks = re.findall("(\d+.\d*)", l.rstrip())
                            assert len(toks) == 4
                            filterLogStat["microbialremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}


            ##
            ## refStats.txt format
            ##
            ## name %unambiguousReads   unambiguousMB   %ambiguousReads ambiguousMB unambiguousReads    ambiguousReads
            ## human_masked 85.24693    498.92052   0.09378 0.55290 3350692 3686
            ## mouse_masked 0.03765 0.21670 0.10802 0.63690 1480    4246
            ## cat_masked   0.01862 0.09568 0.02514 0.14820 732 988
            ## dog_masked   0.00697 0.03815 0.01384 0.08160 274 544
            ##
            if os.path.isfile("refStats.txt"):
                refStatsFile = os.path.join(outDir, "refStats.txt")
                write_files(fh, "refStats", refStatsFile, log)

                with open(refStatsFile) as RFH:
                    ## Need to report 0 if nothing matched
                    refStats['human'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                    refStats['cat'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                    refStats['dog'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                    refStats['mouse'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}

                    for l in RFH:
                        if l:
                            if l.startswith("#"):
                                continue

                            toks = l.rstrip().split()
                            assert len(toks) == 7

                            ## the number and percent of reads that map unambiguously or ambiguously to human, cat, dog.
                            ## take the sum of the two numbers (ambiguous plus unambiguous) to use as the final percentage.
                            if l.startswith("human"):
                                refStats['human'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                            if l.startswith("cat"):
                                refStats['cat'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                            if l.startswith("dog"):
                                refStats['dog'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                            if l.startswith("mouse"):
                                refStats['mouse'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}

                log.debug("refStats.txt: %s", str(refStats))

            if os.path.isfile("cardinality.txt"):
                cardinalityFile = os.path.join(outDir, "cardinality.txt")
                write_files(fh, "cardinality", cardinalityFile, log)

                with open(cardinalityFile) as CFH:
                    for l in CFH:
                        if l:
                            cardinality = l.rstrip()

                log_and_print("Cardinality: %s" % cardinality)


            ## itags
            if prodType == "ITAG":
                prefix = safe_basename(fastq, log)[0].replace(".fastq", "").replace(".gz", "")
                itaggerStatsTxtFile = os.path.join(outDir, "itagger", prefix + ".stats.txt")
                itaggerLogFile = os.path.join(outDir, "itagger", prefix + ".log.gz")
                itaggerSeqobsFile = os.path.join(outDir, "itagger", prefix + ".seqobs.gz")
                if os.path.isfile(itaggerStatsTxtFile):
                    write_files(fh, "itagger_stats_txt", itaggerStatsTxtFile, log)
                    write_files(fh, "itagger_logs_gz", itaggerLogFile, log)
                    write_files(fh, "itagger_seqobs_gz", itaggerSeqobsFile, log)


        ###########################################################
        log_and_print("Write to stats file %s" % STATS_LIST_FILE_NAME)
        ###########################################################
        if os.path.isfile(STATS_LIST_FILE_NAME):
            os.remove(STATS_LIST_FILE_NAME)

        bbtoolsVersion = None

        with open(STATS_LIST_FILE_NAME, 'a') as fh:
            for key in readCounts:
                write_stats(fh, key, readCounts[key], log)

            ## 01.06.2016 Found mismatches in filtered fastq read count
            ## Overwrite "outputReads" in filterStats.txt
            write_stats(fh, "outputReads", readCounts["newCnt"], log)
            write_stats(fh, "outputReadsRrna", readCounts["newCntRrna"], log)

            ## Write refStats to filterStats.txt file
            for key in refStats:
                for k in refStats[key]:
                    write_stats(fh, key+'_'+k, refStats[key][k], log)

            write_stats(fh, "cardinality", cardinality, log)

            ## Write refStats to filterStats.txt file
            for key in filterLogStat:
                for k in filterLogStat[key]:
                    write_stats(fh, key + '_' + k, filterLogStat[key][k], log)

            bbversionCmd = get_tool_path("bbversion.sh", "bbtools")
            cmd = "%s" % (bbversionCmd)
            stdOut, _, exitCode = run_sh_command(cmd, True, log)
            assert stdOut is not None
            bbtoolsVersion = stdOut.strip()

            ## 05112017 Now bbtools version = bbmap version
            # bbtoolsVersion = bbmapVersion if bbmapVersion else "37.xx"
            assert bbtoolsVersion is not None
            write_stats(fh, "filter_tool", "bbtools " + bbtoolsVersion, log)
            write_stats(fh, "filter", VERSION, log)


            ## Version recording
            if bbdukVersion is None: bbdukVersion = bbtoolsVersion
            if bbmapVersion is None: bbmapVersion = bbtoolsVersion
            write_stats(fh, "bbduk_version", bbdukVersion, log)
            write_stats(fh, "bbmap_version", bbmapVersion, log)


            ## itags
            if prodType == "ITAG":
                prefix = safe_basename(fastq, log)[0].replace(".fastq", "").replace(".gz", "")
                itaggerStatsTxtFile = os.path.join(outDir, "itagger", prefix + ".stats.txt")

                log.debug("itagger stat file: %s", itaggerStatsTxtFile)
                log_and_print("Adding iTags stats values to the filtering stat file.")

                with open(itaggerStatsTxtFile, 'r') as ISFH:
                    for l in ISFH:
                        if l:
                            key, val = l.rstrip().split('=')
                            if key and val:
                                if key.endswith("_count"):
                                    val = int(val) * 2 ## pair-based count ---> read-based count

                                write_stats(fh, key, val, log)


        log_and_print("Sequnit name: %s" % fastq)
        seqUnitName = os.path.basename(fastq)
        ###########################################################
        log_and_print("Creating chaff tar ball file")
        ###########################################################

        ## collect chaff file list from file-list.txt file
        chaffFileList = []
        with open(FILE_LIST_FILE_NAME, 'r') as FLFH:
            for l in FLFH:
                if l.startswith("#"):
                    continue
                if l.startswith("chaff"):
                    chaffFile = l.strip().split("=")[1]
                    if chaffFile.endswith(".fq.gz"):
                        log.debug("chaff file = %s, exist? = %s", chaffFile, os.path.isfile(chaffFile))
                        if os.path.isfile(chaffFile):
                            chaffFileList.append(chaffFile)


        CHAFF_README_HEADER = """This tarball file is a collection of removed reads from:
    - Synthetic Contamination sequence removal (synth1.fq.gz)
    - Short Synthetic Contamination sequence removal (synth2.fq.gz)
    - Human contamination removal (human.fq.gz)
    - Cat contamination removal (human.fq.gz)
    - Dog contamination removal (human.fq.gz)
    - Mouse contamination removal (human.fq.gz)

    If the tarball does not include some of the files then no reads matched the references for that type.

    Chaff files:
    """

        chaffReadMeFile = os.path.join(outputPath, "chaff_readme.txt")
        with open(chaffReadMeFile, 'w') as CRMFH:
            CRMFH.write(CHAFF_README_HEADER)
            for f in chaffFileList:
                CRMFH.write(f + '\n')

        chaffTarFileName = safe_basename(fastq, log)[0].replace(".fastq", "").replace(".gz", "") + ".chaff.tar"
        chaffTarFile = os.path.join(outputPath, chaffTarFileName)

        tarCmd = "cd %s && tar cvf %s chaff_readme.txt -C %s " % (outputPath, chaffTarFile, outputPath)

        for f in chaffFileList:
            tarCmd += f + " "

        # tarCmd += "--remove-files" ## delete chaff files

        stdOut, _, exitCode = run_sh_command(tarCmd, True, log, True)

        if exitCode != 0:
            log.error("Failed to run tar command: %s", tarCmd)
            return None
        else:
            ## remove chaff files after tar'd
            for f in chaffFileList:
                rmCmd = "rm %s" % (os.path.join(outputPath, f))
                run_sh_command(rmCmd, True, log, True)
                log.debug("chaff file removed = %s", f)

            rmCmd = "rm %s" % (chaffReadMeFile)
            run_sh_command(rmCmd, True, log, True)
            log.debug("chaff_readme file removed = %s", chaffReadMeFile)

        assert os.path.isfile(chaffTarFile), "ERROR: chaff tar file not found."

        with open(FILE_LIST_FILE_NAME, 'a') as fh:
            write_files(fh, "chaffTarFileName", chaffTarFile, log)

        log_and_print("Chaff Tar ball successfully created: %s" % chaffTarFile)



        ###########################################################
        log_and_print("Creating methods txt file")
        ###########################################################
        ## read the method description from filter/filter_desc/*.txt
        InputReads = ""
        lowQualPerc = ""
        lowQualReads = ""
        artifactPerc = ""
        artifactReads = ""
        riboRnaPerc = "0.00"
        riboRnaReads = "0"
        microbialPerc = "0.00"
        microbialReads = "0"
        humanPerc = dogPerc = catPerc = mousePerc = "0.00"
        humanReads = dogReads = catReads = mouseReads = 0
        RemainReads = ""
        inputBases = ""
        outputBases = ""
        ktrimmed_total_removed_percreads = "0.00"
        ktrimmed_total_removed_percbases = "0.00"
        ktrimmed_total_removed_numbases = 0
        ktrimmed_total_removed_numreads = 0

        with open(STATS_LIST_FILE_NAME, 'r') as statsListFH:
            for l in statsListFH:
                if len(l) > 1:
                    toks = l.strip().split("=")
                    if toks[0] == 'rawCnt':
                        InputReads = toks[1]
                        # lowQualTotalReads = InputReads
                        # artifactTotalReads = InputReads
                        # riboRnaTotalReads = InputReads
                        # humanTotalReads = dogTotalReads = catTotalReads = mouseTotalReads = InputReads

                    elif toks[0] == 'newCnt':
                        RemainReads = toks[1]
                    elif toks[0].lower() == 'lowqualitydiscards_percreads':
                        lowQualPerc = toks[1]
                    elif toks[0].lower() == 'lowqualitydiscards_numreads':
                        lowQualReads = toks[1]
                    elif toks[0].lower() == 'contaminants_percreads':
                        artifactPerc = toks[1]
                    elif toks[0].lower() == 'contaminants_numreads':
                        artifactReads = toks[1]
                    elif toks[0].lower() == 'ribosomalsequenceremoved_percreads':
                        riboRnaPerc = toks[1]
                    elif toks[0].lower() == 'ribosomalsequenceremoved_numreads':
                        riboRnaReads = toks[1]
                    elif toks[0].lower() == 'microbialremoved_percreads':
                        microbialPerc = toks[1]
                    elif toks[0].lower() == 'microbialremoved_numreads':
                        microbialReads = toks[1]

                    elif toks[0].lower() == 'human_unambiguousreads':
                        humanReads = toks[1]
                    elif toks[0].lower() == 'cat_unambiguousreads':
                        dogReads = toks[1]
                    elif toks[0].lower() == 'dog_unambiguousreads':
                        catReads = toks[1]
                    elif toks[0].lower() == 'mouse_unambiguousreads':
                        mouseReads = toks[1]

                    elif toks[0].lower() == 'human_unambiguousreadsperc':
                        humanPerc = "%.2f" % (float(toks[1]))
                    elif toks[0].lower() == 'cat_unambiguousreadsperc':
                        dogPerc = "%.2f" % (float(toks[1]))
                    elif toks[0].lower() == 'dog_unambiguousreadsperc':
                        catPerc = "%.2f" % (float(toks[1]))
                    elif toks[0].lower() == 'mouse_unambiguousreadsperc':
                        mousePerc = "%.2f" % (float(toks[1]))

                    ## RQCSUPPORT-1987
                    ## need to report the number of adapter trimmed reads and bases in the *.filtered-report.txt
                    elif toks[0].lower() == 'ktrimmed_total_removed_percreads':
                        ktrimmed_total_removed_percreads = "%.2f" % (float(toks[1]))
                    elif toks[0].lower() == 'ktrimmed_total_removed_percbases':
                        ktrimmed_total_removed_percbases = "%.2f" % (float(toks[1]))
                    elif toks[0] == 'ktrimmed_total_removed_numbases':
                        ktrimmed_total_removed_numbases = toks[1]
                    elif toks[0] == 'ktrimmed_total_removed_numreads':
                        ktrimmed_total_removed_numreads = toks[1]

                    elif toks[0] == 'inputBases':
                        inputBases = toks[1]
                    elif toks[0] == 'outputBases':
                        outputBases = toks[1]


                    percReadsRemoved = "n/a"
                    try:
                        percReadsRemoved = "%.2f" %  (100.0 - (100.0 * float(RemainReads) / float(InputReads)))
                    except:
                        pass

                    percBasesRemoved = "n/a"
                    try:
                        percBasesRemoved = "%.2f" % (100.0 - (100.0 * float(outputBases) / float(inputBases)))
                    except:
                        pass

        filteredSeqUnitName = os.path.basename(filteredFastq)
        statsReportFileName = os.path.basename(fastq).replace(".fastq", "").replace(".gz", "") + ".filtered-report.txt"
        statsMethodFileName = os.path.basename(fastq).replace(".fastq", "").replace(".gz", "") + ".filtered-methods.txt"

        methodsText = ""
        seqUnitNameToWrite = seqUnitName if seqUnitName.endswith(".gz") else seqUnitName + ".gz"
        fullSeqUnitName = fastq if fastq.endswith(".gz") else fastq + ".gz"

        import io
        with io.open(os.path.join(SRC_ROOT, "filter_desc/", FILTER_METHODS_TXT[prodType]), "r", encoding="utf-8") as MFH:
            for l in MFH:
                if l.startswith("#"):
                    continue

                methodsText += l
        libraryName = libraryType = 'N.A.'
        sequencingRunMode = sequencerName = 'N.A.'
        # print("DEBUG statsReportFileName=%s" % statsReportFileName)
        # print("DEBUG statsMethodFileName=%s" % statsMethodFileName)
        methodsText = methodsText.replace("<library_name>", libraryName)\
                                 .replace("<library_type>", libraryType)\
                                 .replace("<library type>", libraryType)\
                                 .replace("<sequencing_run_mode>", sequencingRunMode)\
                                 .replace("<sequencer_name>", sequencerName)\
                                 .replace("<raw_read_count>", InputReads)\
                                 .replace("<raw_base_count>", inputBases)\
                                 .replace("<bbtools_version>", bbtoolsVersion)\
                                 .replace("<bbtools version>", bbtoolsVersion)\
                                 .replace("<contams_list>", "contamination")\
                                 .replace("<filtered_read_count>", RemainReads)\
                                 .replace("<filtered_base_count>", outputBases)\
                                 .replace("<rqc_filter_command>", safe_basename(cof, log)[0])\
                                 .replace("<filter_command_file>", safe_basename(cof, log)[0])\
                                 .replace("<filter_chaff_fastq_file>", chaffTarFileName)

        with io.open(statsMethodFileName, "w", encoding="utf-8") as statsMethodFileNameFH:
            statsMethodFileNameFH.write(methodsText)

        if os.path.isfile(statsMethodFileName):
            with open(FILE_LIST_FILE_NAME, 'a') as fh:
                write_files(fh, "filterMethodsFileName", os.path.join(outputPath, statsMethodFileName), log)

            log_and_print("Filter methods file successfully created: %s" % os.path.join(outputPath, statsMethodFileName))

        else:
            log.error("Filter method file does not exist [%s]", statsMethodFileName)
            return None



        ###########################################################
        log_and_print("Creating report txt file")
        ###########################################################

        reportText = ""

        ## Read report template file
        with open(os.path.join(SRC_ROOT, "filter_desc/filter_report.txt"), 'r') as RFH:
            for l in RFH:
                if l.startswith("#"):
                    continue
                reportText += l

        reportText = reportText.replace("<timestamp>", datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))\
                               .replace("<seq_unit_base_name>", seqUnitNameToWrite)\
                               .replace("<library_name>", libraryName)\
                               .replace("<filter_fastq_file>", filteredSeqUnitName)\
                               .replace("<input_read_count>", InputReads)\
                               .replace("<input_base_count>", inputBases)\
                               .replace("<low_quality_read_count_pct>", lowQualPerc)\
                               .replace("<low_quality_read_count>", lowQualReads)\
                               .replace("<artifact_read_count_pct>", artifactPerc)\
                               .replace("<artifact_read_count>", artifactReads)\
                               .replace("<rrna_read_count_pct>", riboRnaPerc)\
                               .replace("<rrna_read_count>", riboRnaReads)\
                               .replace("<microbial_read_count_pct>", microbialPerc)\
                               .replace("<microbial_read_count>", microbialReads)\
                               .replace("<human_read_count_pct>", humanPerc)\
                               .replace("<human_read_count>", str(humanReads))\
                               .replace("<dog_read_count_pct>", dogPerc)\
                               .replace("<dog_read_count>", str(dogReads))\
                               .replace("<cat_read_count_pct>", catPerc)\
                               .replace("<cat_read_count>", str(catReads))\
                               .replace("<mouse_read_count_pct>", mousePerc)\
                               .replace("<mouse_read_count>", str(mouseReads))\
                               .replace("<ktrimmed_total_removed_percreads>", ktrimmed_total_removed_percreads)\
                               .replace("<ktrimmed_total_removed_percbases>", ktrimmed_total_removed_percbases)\
                               .replace("<ktrimmed_total_removed_numreads>", str(ktrimmed_total_removed_numreads))\
                               .replace("<ktrimmed_total_removed_numbases>", str(ktrimmed_total_removed_numbases))\
                               .replace("<filter_read_count>", RemainReads)\
                               .replace("<filter_base_count>", outputBases)\
                               .replace("<filter_read_count_pct>", percReadsRemoved)\
                               .replace("<filter_base_count_pct>", percBasesRemoved)\
                               .replace("<rqc_filter_command>", safe_basename(cof, log)[0])\
                               .replace("<filter_command_file>", safe_basename(cof, log)[0])

        with open(statsReportFileName, 'w') as statsReportFileNameFH:
            statsReportFileNameFH.write(reportText)

        if os.path.isfile(statsReportFileName):
            with open(FILE_LIST_FILE_NAME, 'a') as fh:
                write_files(fh, "filterReportFileName", os.path.join(outputPath, statsReportFileName), log)

            log_and_print("Filter report file successfully created: %s" % os.path.join(outputPath, statsReportFileName))

            with open(STATS_LIST_FILE_NAME, 'a') as fh:
                write_stats(fh, "filter_report_input_read_count", InputReads, log)
                write_stats(fh, "filter_report_input_base_count", inputBases, log)
                write_stats(fh, "filter_report_low_quality_read_count_pct", lowQualPerc, log)
                write_stats(fh, "filter_report_low_quality_read_count", lowQualReads, log)
                write_stats(fh, "filter_report_artifact_read_count_pct", artifactPerc, log)
                write_stats(fh, "filter_report_artifact_read_count", artifactReads, log)
                write_stats(fh, "filter_report_rrna_read_count_pct", riboRnaPerc, log)
                write_stats(fh, "filter_report_rrna_read_count", riboRnaReads, log)
                write_stats(fh, "filter_report_microbial_read_count_pct", microbialPerc, log)
                write_stats(fh, "filter_report_microbial_read_count", microbialReads, log)
                write_stats(fh, "filter_report_human_read_count_pct", humanPerc, log)
                write_stats(fh, "filter_report_human_read_count", str(humanReads), log)
                write_stats(fh, "filter_report_dog_read_count_pct", dogPerc, log)
                write_stats(fh, "filter_report_dog_read_count", str(dogReads), log)
                write_stats(fh, "filter_report_cat_read_count_pct", catPerc, log)
                write_stats(fh, "filter_report_cat_read_count", str(catReads), log)
                write_stats(fh, "filter_report_mouse_read_count_pct", mousePerc, log)
                write_stats(fh, "filter_report_mouse_read_count", str(mouseReads), log)
                write_stats(fh, "filter_report_filter_read_count", RemainReads, log)
                write_stats(fh, "filter_report_filter_base_count", outputBases, log)
                write_stats(fh, "filter_report_filter_read_count_pct", percReadsRemoved, log)
                write_stats(fh, "filter_report_filter_base_count_pct", percBasesRemoved, log)

        else:
            log.error("Filter report file does not exist [%s]", statsReportFileName)
            return None

        checkpoint(POST_END, status)
        status = POST_END
    else:
        log_and_print('No need to do post processing.')

    return filteredFastq, status

"""
Gen read kmer hist plots

"""
def gen_read_khist_plot(kHistFile, log):
    log_and_print("\n\n - RUN gen_read_khist_plot\n")

    log.debug("khist file: %s", kHistFile)

    ## Gen kmer hist Plot
    if not os.path.isfile(kHistFile):
        log.error("Cannot find khist.txt file.")
        return None, None

    ## Load data from txt
    rawDataMatrix = np.loadtxt(kHistFile, delimiter="\t", comments="#", skiprows=1)

    try:
        assert len(rawDataMatrix[0][:]) == 2, "ERROR: read kmer histogram text file format error in %s: %s" % (kHistFile, str(rawDataMatrix))
    except:
        if len(rawDataMatrix) != 2:
            return None, None

    ## Data
    ## #Depth   Count
    ## 1    1915198
    ## 2    498682
    ## 3    183976
    ## 4    144168
    ## 5    113947

    fig, ax = plt.subplots()

    markerSize = 3.5
    lineWidth = 1.0
    try:
        p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 1], "r", marker="o", markersize=markerSize, linewidth=lineWidth, label="depth")
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=list(rawDataMatrix[:, 1])))
    except:
        p1 = ax.plot(rawDataMatrix[0], rawDataMatrix[1], "r", marker="o", markersize=markerSize, linewidth=lineWidth, label="depth")

    ax.set_xlabel("Depth (log)", fontsize=12, alpha=0.5)
    ax.set_ylabel("Count (log)", fontsize=12, alpha=0.5)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(color="gray", linestyle=":")

    pngPlotFile = os.path.join(kHistFile + ".read_len_hist.png")
    htmlPlotFile = os.path.join(kHistFile + ".read_len_hist.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlPlotFile)

    ## Save Matplotlib plot in png format
    plt.savefig(pngPlotFile, dpi=fig.dpi)


    return pngPlotFile, htmlPlotFile


##==============================================================================
## Helper functions

def clean_up(fList, log):
    log_and_print("\n\n%s - CLEAN UP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % ( color['pink'], color['']))

    for f in fList:
        if os.path.isfile(f):
            log_and_print("Removing %s ... ", f)
            os.remove(f)

    log_and_print("CLEAN UP - completed")


def create_shell(shName, cmdArray):
    with open(shName, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("set -e\n")
        fh.write("set -o pipefail\n")
        #fh.write("module unload blast+; module load blast+\n")
        for item in cmdArray:
            fh.write("%s\n" % item)
        os.chmod(shName, 0755)  #-rwxr-xr-x


def return_values(vals, flag=True):
    if flag:
        return vals
    else:
        return [None for _ in range(len(vals))]


## checkpoint logging
def checkpoint(status, fromStatus=PIPE_START):
    if status == PIPE_START or STEP_ORDER[status] > STEP_ORDER[fromStatus]:
        checkpoint_step(STATUS_LOG_FNAME, status)


def write_stats(fh, k, v, log):
    line = "%s=%s" % (k, v)
    fh.write("%s\n" % line)
    log_and_print("Added '%s' in the stats file" % line)


def write_files(fh, k, v, log):
    line = "%s=%s" % (k, v)
    fh.write("%s\n" % line)
    log_and_print("Added '%s' in the files file" %line)


def file_name_trim(fname):
    return fname.replace(".gz", "").replace(".fastq", "").replace(".fasta", "")

def read_qc(odir, fastq, status):
    log_and_print("\n\n%s - RUN QC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    if not os.path.isfile(fastq):
        return False, status

    qcdir = os.path.join(odir, READQC_DIR)

    if STEP_ORDER[status] < STEP_ORDER[QC_END]:
        checkpoint(QC_START, status)

        qcexe = os.path.join(os.path.dirname(__file__), 'readqc.py')
        cmd = '%s -o %s -f %s --skip-blast' % (qcexe, qcdir, fastq)
        # print('DEBUG : %s' % cmd)
        stdOut, stdErr, exitCode = run_sh_command(cmd, True)
        if exitCode != 0:
            print('ERROR : %s' % stdErr)
            return False, qcdir, status

        checkpoint(QC_END, status)
        status = QC_END
    else:
        log_and_print("No need to do qc step.")

    return True, qcdir, status

def do_html(odir, qcdir, rawFastq, filteredFastq, status):
    exedir = os.path.dirname(os.path.abspath(__file__))
    log_and_print("\n\n%s - Create HTML file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    print('Output dir - %s' % odir)
    fname = os.path.basename(rawFastq)
    print('Create HTML page in %s for %s ..' % (odir, fname))
    temp = os.path.join(exedir, 'template/filter_template.html')
    print('open the template file %s' % temp)

    stats = get_dict_obj(os.path.join(odir, 'filterStats.txt'))

    tok_map = {
            'filter_report_input_read_count' : {'token' : '[_RAW-READ-CNT_]', 'type': 'bigint'},
            'filter_report_input_base_count' : {'token' : '[_RAW-BASE-CNT_]', 'type': 'bigint'},
            'filter_report_filter_read_count' : {'token' : '[_FILTERED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_filter_base_count' : {'token' : '[_FILTERED-BASE-CNT_]', 'type': 'bigint'},
            'filter_report_filter_read_count_pct' : {'token' : '[_REMOVED-READ-PCT_]', 'type': 'raw'},
            'filter_report_filter_base_count_pct' : {'token' : '[_REMOVED-BASE-PCT_]', 'type': 'raw'},

            'filter_report_low_quality_read_count' : {'token' : '[_LOW-QUAL-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_low_quality_read_count_pct' : {'token' : '[_LOW-QUAL-REMOVED-PCT_]', 'type': 'raw'},

            'filter_report_artifact_read_count' : {'token' : '[_ARTI-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_artifact_read_count_pct' : {'token' : '[_ARTI-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_rrna_read_count' : {'token' : '[_RRNA-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_rrna_read_count_pct' : {'token' : '[_RRNA-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_microbial_read_count' : {'token' : '[_MICROBE-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_microbial_read_count_pct' : {'token' : '[_MICROBE-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_human_read_count' : {'token' : '[_HUMAN-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_human_read_count_pct' : {'token' : '[_HUMAN-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_dog_read_count' : {'token' : '[_DOG-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_dog_read_count_pct' : {'token' : '[_DOG-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_cat_read_count' : {'token' : '[_CAT-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_cat_read_count_pct' : {'token' : '[_CAT-REMOVED-PCT_]', 'type': 'raw'},
            'filter_report_mouse_read_count' : {'token' : '[_MOUSE-REMOVED-READ-CNT_]', 'type': 'bigint'},
            'filter_report_mouse_read_count_pct' : {'token' : '[_MOUSE-REMOVED-PCT_]', 'type': 'raw'},
    }

    with open(temp, 'r') as fh:
        html = fh.read()

        ## do the place-holder replacement !!
        html = html.replace('[_REPORT-DATE_]', '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
        html = html.replace('[_INPUT-FILE-NAME_]', os.path.basename(rawFastq))
        html = html.replace('[_RAW-FILE-LOCATION_]', rawFastq)
        html = html.replace('[_FILTERED-FILE-LOCATION_]', filteredFastq)
        fsize = format(os.stat(rawFastq).st_size / (1024*1024), ',')
        html = html.replace('[_RAW-FILE-SIZE_]', fsize)
        fsize = format(os.stat(filteredFastq).st_size / (1024*1024), ',')
        html = html.replace('[_FILTERED-FILE-SIZE_]', fsize)

        for key in tok_map:
            dat = tok_map[key]
            html = html.replace(dat['token'], pipeline_val(key, dat, stats))

        # readqc on the filter file
        fdir = os.path.join(odir, READQC_DIR)
        hbody = do_html_body(fdir, odir)
        html = html.replace('[_FILTERED-READ-QC_]', hbody)

        ## write the html to file
        idxfile = os.path.join(odir, 'index.html')
        with open(idxfile, 'w') as fh2:
            fh2.write(html)
        print('HTML index file written to %s' % idxfile)

        # copy the css file
        cssdir = os.path.join(exedir, 'css')
        todir = os.path.join(odir, 'css')
        if os.path.isdir(todir):
            shutil.rmtree(todir)
        shutil.copytree(cssdir, todir, False, None)

        # copy the image file
        imgdir = os.path.join(exedir, 'images')
        todir = os.path.join(odir, 'images')
        if os.path.isdir(todir):
            shutil.rmtree(todir)
        shutil.copytree(imgdir, todir, False, None)

    return status

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## main program
if __name__ == "__main__":
    ## Parse options
    usage = "* Filter Pipeline, version %s\n" % (VERSION)
    origCmd = ' '.join(sys.argv)

    ## command line options
    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f", "--fastq", dest="fastq", help="Set input fastq file (full path to fastq)", required=True)
    parser.add_argument("-o", "--output-path", dest="outputPath", help="Set output path to write to", required=True)
    parser.add_argument("-p", "--prod-type", dest="prodType", help="Set product type: [DNA | FUNGAL | METAGENOME | \
VIRAL-METAGENOME |SAG | ISO | SAG | CELL-ENRICHMENT | PLANT-2x150 | PLANT-2x250 | RNA | SMRNA | METATRANSCRIPTOME \
(or MTF) | LFPE | CLRS | CLIP-PE | NEXTERA-LMP (or NEXTERA) | ITAG | MICROTRANS | BISULPHITE | 3PRIMERNA | CHIPSEQ \
| RNAWOHUMAN]", required=True)
    parser.add_argument("-t", "--taxlist", dest="taxList", help="A list of taxid(s) to exclude in CSV format", required=False)
    parser.add_argument("-ap", "--ap", dest="apNum", help="Set AP (Analysis Project) ID. Ex) -ap 123 or -ap 123,456,789", required=False)
    parser.add_argument("-at", "--at", dest="atNum", help="Set AT (Analysis Task) ID. Ex) -at 123 or -at 123,456,789", required=False)
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    ## switches
    parser.add_argument("-a", "--aggressive", dest="enableAggressive", action="store_true", help="Enable aggressive=t and microbebuild=3", default=False)
    parser.add_argument("-b", "--skip-blast", dest="skipBlastFlag", action="store_true", help="Skip Blast search", default=False)
    parser.add_argument("-c", "--skip-cleanup", dest="skipCleanUp", action="store_true", help="Skip file clean up after pipeline complete", default=False)
    parser.add_argument("-d", "--debug", dest="doDebug", action="store_true", help="Enable debug mode")
    parser.add_argument("-l", "--disable-clumpify", dest="disableClumpify", action="store_true", help="Disable clumpify", default=False)
    parser.add_argument("-m", "--skip-microbes-removal", dest="disableRmoveMicrobes", action="store_true", help="Skip microbes removal", default=False)
    parser.add_argument("-r", "--contam", dest="enableRmoveMicrobes", action="store_true", help="Enable removemicrobes=t", default=False)
    parser.add_argument("-s", "--skip-subsampleqc", dest="skipSubsampleQc", action="store_true", help="Skip the subsample and qc step", default=False)
    parser.add_argument("-pl", "--print-log", dest="print_log", default = False, action = "store_true", help = "print log to screen")

    ## produce html when processing is done
    parser.add_argument("-html", "--html", action="store_true", help="Create html file", dest="html", default=False, required=False)

    skipCleanUp = False
    skipSubsampleQc = False
    outputPath = None ## output path, defaults to current working directory
    fastq = None ## full path to input fastq
    # logLevel = "DEBUG"
    enableRmoveMicrobes = False
    disableRmoveMicrobes = False
    disableClumpify = False
    enableAggressive = False
    skipBlastFlag = False
    apNum = None
    atNum = None

    taxList = ""

    options = parser.parse_args()
    print_log = options.print_log

    if options.outputPath:
        outputPath = options.outputPath

    if not outputPath:
        outputPath = os.getcwd()

    ## create output_directory if it doesn't exist
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)

    outputPath = os.path.realpath(outputPath)
    outputPath = outputPath.replace("/chos", "") if outputPath.startswith("/chos") else outputPath

    ## initialize my logger
    logFile = os.path.join(outputPath, "rqc_filter_pipeline.log")

    print "Started filtering pipeline with %s, writing log to: %s" % (SCRIPT_NAME, logFile)

    log = get_logger("filter", logFile, logLevel, print_log, True)

    if options.doDebug:
        DEBUG = True

    if options.skipCleanUp:
        skipCleanUp = True

    if options.skipBlastFlag:
        skipBlastFlag = True

    if options.skipSubsampleQc:
        skipSubsampleQc = True

    if options.enableRmoveMicrobes:
        if options.disableRmoveMicrobes:
            log.error("Conflict in option parameters: cannot set skip-contam with skip-microbes-removal.")
            sys.exit(0)

        enableRmoveMicrobes = True

    if options.disableRmoveMicrobes:
        if options.enableRmoveMicrobes:
            log.error("Conflict in option parameters: cannot set skip-contam with skip-microbes-removal.")
            sys.exit(0)

        disableRmoveMicrobes = True

    if options.enableAggressive:
        enableAggressive = True

    if options.disableClumpify:
        disableClumpify = True

    if options.fastq:
        fastq = options.fastq

    if options.taxList:
        taxList = options.taxList

    if options.apNum:
        apNum = options.apNum
    if options.atNum:
        atNum = options.atNum


    log_and_print("%s" % '#' * 80)
    log_and_print("  Filtering pipeline (version %s)" % VERSION)
    log_and_print("%s" % '#' * 80)


    prodType = options.prodType.upper()
    skipSubsampleQcFlag = options.skipSubsampleQc ## run subsample and qc or not

    if not os.path.isdir(outputPath):
        log.error("Cannot work with directory: %s", outputPath)

    ## check for fastq file
    if fastq:
        if not os.path.isfile(fastq):
            log.error("Input fastq file, %s not found, abort!", fastq)
    else:
        log.error("No fastq defined, abort!")

    fastq = os.path.realpath(fastq)
    fastq = fastq.replace("/chos", "") if fastq.startswith("/chos") else fastq


    ##--------------------------------
    ## init log
    log_and_print("Started pipeline, writing log to: %s" % logFile)
    log_and_print("CMD: %s" % origCmd)
    log_and_print("Fastq file: %s" % fastq)
    log_and_print("Output path: %s" % outputPath)

    os.chdir(outputPath)

    cycle = 0
    cycleMax = 1
    bIsPaired = False

    status = get_status(STATUS_LOG_FNAME, log)
    log_and_print("Starting pipeline at [%s]" % status)

    if status == PIPE_START:
        checkpoint(PIPE_START)

    ## main loop: retry upto cycleMax times
    while cycle < cycleMax:

        cycle += 1

        log_and_print("ATTEMPT [%d]" % cycle)

        filesToRemove = []  # list of intermediate files for clean up
        lastFastq = fastq   # lastFastq : fastq produced by each step, init to input
        rRnaFilterFile = None # MTF generates this 2nd filtered output file
        subsampledFastq = None
        bIsPaired = None
        filteredReadNum = -1
        filteredReadNumRrna = -1
        FragFile = SingletonFile = UnknownFile = None

        if cycle > 1:
            status = get_status(STATUS_LOG_FNAME, log)

        if status != PIPE_COMPLETE:
            ##
            ## Run rqcfilter.sh
            ##
            lastFastq, rRnaFilterFile, fragFile, singletonFile, unknownFile, filteredReadNum, filteredReadNumRrna, cof, status = run_rqcfilter(lastFastq, outputPath, prodType, status, enableRmoveMicrobes, enableAggressive, disableRmoveMicrobes, disableClumpify, taxList, log) ## Only MTF type generates rRnaFilterFile
            if filteredReadNum == 0:
                break

            ##--------------------------------
            ## Run post processing
            if lastFastq is not None and lastFastq != -1:
                lastFastq, status = post_process(fastq, outputPath, lastFastq, rRnaFilterFile, fragFile, singletonFile, unknownFile, prodType, filteredReadNum, filteredReadNumRrna, cof, status, log)
            else:
                print "Failed @ sketch"

            ##--------------------------------
            ## Clean up
            if lastFastq is not None and lastFastq != -1:

                ## run readQC on the filtered fastq
                rtn, qcdir, status = read_qc(outputPath, lastFastq, status)

                ## create html file on the readQC results
                if rtn:
                    status = do_html(outputPath, qcdir, fastq, lastFastq, status)

                exit(0)

                checkpoint(PIPE_COMPLETE)

                if not skipCleanUp:
                    clean_up(filesToRemove, log)
                else:
                    log_and_print("SKIP CLEANUP")

                log_and_print("Pipeline Completed")
                cycle = cycleMax + 1

            else:
                print "Failed @ postprocess"

        else:
            cycle = cycleMax + 1
            log_and_print("Pipeline already completed")


    print "Done."
    sys.exit(0)


## EOF
