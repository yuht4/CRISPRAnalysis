#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
 Read qc utilities

 Created: Sep 12 2016
 sulsj

"""

## the following constatns are specific to THIS pipeline
STATUS_LOG_FNAME = 'status.log'

PIPE_START        = 'start'
RQCFILTER_START   = 'start rqcfilter'

## used by rqcfilter.sh
CLUMPIFY_START     = 'clumpify start'
CLUMPIFY_END       = 'clumpify finish'
KTRIM_START        = 'ktrim start'
KTRIM_END          = 'ktrim finish'
FILTER_START       = 'filter start'
FILTER_END         = 'filter finish'
DELETE_TMP_START   = 'delete temp files start'
DELETE_TMP_END     = 'delete temp files finish'
SHORT_FILTER_START = 'short filter start'
SHORT_FILTER_END   = 'short filter finish'
REMOVECOMMONMICROBES_START = 'removecommonmicrobes start'
REMOVECOMMONMICROBES_END   = 'removecommonmicrobes finish'
DEHUMANIZE_START   = 'dehumanize start'
DEHUMANIZE_END     = 'dehumanize finish'
MERGE_START        = 'merge start'
MERGE_END          = 'merge finish'
KHIST_START        = 'khist start'
KHIST_END          = 'khist finish'
TAXFILTER_START    = 'taxfilter start'
TAXFILTER_END      = 'taxfilter finish'
##

RQCFILTER_END     = 'end rqcfilter'
RQCFILTER_END2    = 'rqcfilter complete'

SUBSAMPLEQC_START = 'start subsampleqc'
SUBSAMPLEQC_END   = 'end subsampleqc'
GC_DIV_START      = 'start gcdiv'
GC_DIV_END        = 'end gcdiv'
ITAGGER_START     = 'start itagger'
ITAGGER_END       = 'end itagger'
MERANAL_START     = 'start uniq mer analysis'
MERANAL_END       = 'end uniq mer analysis'
BLAST_NT_START    = 'start blast search vs nt'
BLAST_NT_END      = 'end blast search vs nt'
KMER_DEPTH_START  = 'start kmer depth analysis'
KMER_DEPTH_END    = 'end kmer depth analysis'
SKETCH_START      = 'start sketch'
SKETCH_END        = 'end sketch'
POST_START        = 'start post process'
POST_END          = 'end post process'
QC_START        = 'start qc'
QC_END          = 'end qc'
HTML_START        = 'start html generation'
HTML_END          = 'end html generation'
PIPE_COMPLETE     = 'complete'

STEP_ORDER = {
    PIPE_START        : 0,

    RQCFILTER_START   : 10,
        CLUMPIFY_START    : 20,
        CLUMPIFY_END      : 30,
        KTRIM_START       : 40,
        KTRIM_END         : 50,
        FILTER_START      : 60,
        FILTER_END        : 70,
        DELETE_TMP_START  : 80,
        DELETE_TMP_END    : 90,
        REMOVECOMMONMICROBES_START: 100,
        REMOVECOMMONMICROBES_END  : 110,
        DEHUMANIZE_START  : 120,
        DEHUMANIZE_END    : 130,
        MERGE_START       : 140,
        MERGE_END         : 150,
        KHIST_START       : 160,
        KHIST_END         : 170,
        TAXFILTER_START   : 180,
        TAXFILTER_END     : 190,

    RQCFILTER_END     : 300,
    RQCFILTER_END2    : 301,

    SUBSAMPLEQC_START : 310,
    SUBSAMPLEQC_END   : 311,

    GC_DIV_START      : 321,
    GC_DIV_END        : 322,

    ITAGGER_START     : 331,
    ITAGGER_END       : 332,

    MERANAL_START     : 341,
    MERANAL_END       : 342,

    BLAST_NT_START    : 351,
    BLAST_NT_END      : 352,

    KMER_DEPTH_START  : 361,
    KMER_DEPTH_END    : 362,

    SKETCH_START      : 371,
    SKETCH_END        : 371,

    POST_START        : 381,
    POST_END          : 382,

    QC_START          : 400,
    QC_END            : 410,

    HTML_START        : 420,
    HTML_END          : 430,

    PIPE_COMPLETE     : 999
}


## stats for subsample and qc
FILTER_READ_GC_MEAN = "read GC mean"
FILTER_READ_GC_STDEV = "read GC std"
FILTER_READ_GC_MED = "read GC median"
FILTER_READ_GC_MODE = "read GC mode"
FILTER_READ_GC_TEXT = "read GC text hist"
FILTER_READ_GC_PLOT = "read GC plot"
FILTER_READ_GC_D3_HTML_PLOT = "FILTER_READ_GC_D3_HTML_PLOT"
FILTER_READ_BASE_COUNT = "read base count"
FILTER_READ_COUNT = "read count"
FILTER_READ_SAMPLED_COUNT = "read sampled count"


## for khist option enabled
FILTER_READ_KHIST_TEXT = "FILTER_READ_KHIST_TEXT"
FILTER_READ_KHIST_PLOT = "FILTER_READ_KHIST_PLOT"
FILTER_READ_KHIST_D3_HTML_PLOT = "FILTER_READ_KHIST_D3_HTML_PLOT"

## GAA-1738
FILTER_READ_LENGTH_HIST_TEXT = "FILTER_READ_LENGTH_HIST_TEXT"
FILTER_READ_LENGTH_HIST_PLOT = "FILTER_READ_LENGTH_HIST_PLOT"
FILTER_READ_LENGTH_HIST_D3_HTML_PLOT = "FILTER_READ_LENGTH_HIST_D3_HTML_PLOT"

FILTER_READ_Q20_READ1 = "read q20 read1"
FILTER_READ_Q20_READ2 = "read q20 read2"

FILTER_READ_QUAL_POS_QRPT_MERGED = "FILTER_READ_QUAL_POS_QRPT_MERGED"
FILTER_READ_QUAL_POS_PLOT_MERGED = "read qual pos plot merged"
FILTER_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT = "FILTER_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT"

FILTER_READ_BASE_COUNT_TEXT_1 = "read base count text 1"
FILTER_READ_BASE_COUNT_TEXT_2 = "read base count text 2"
FILTER_READ_BASE_COUNT_PLOT_1 = "read base count plot 1"
FILTER_READ_BASE_COUNT_D3_HTML_PLOT_1 = "FILTER_READ_BASE_COUNT_D3_HTML_PLOT_1" # cyclone nucleotide comp plot
FILTER_READ_BASE_COUNT_PLOT_2 = "read base count plot 2"
FILTER_READ_BASE_COUNT_D3_HTML_PLOT_2 = "FILTER_READ_BASE_COUNT_D3_HTML_PLOT_2" # cyclone nucleotide comp plot

FILTER_READ_INSERT_SIZE_HISTO_PLOT = "FILTER_READ_INSERT_SIZE_HISTO_PLOT"
FILTER_READ_INSERT_SIZE_HISTO_DATA = "FILTER_READ_INSERT_SIZE_HISTO_DATA"
FILTER_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT = "FILTER_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT"

FILTER_READ_INSERT_SIZE_TOTAL_TIME = "FILTER_READ_INSERT_SIZE_TOTAL_TIME"
FILTER_READ_INSERT_SIZE_NUM_READS = "FILTER_READ_INSERT_SIZE_NUM_READS"
FILTER_READ_INSERT_SIZE_JOINED_NUM = "FILTER_READ_INSERT_SIZE_JOINED_NUM"
FILTER_READ_INSERT_SIZE_JOINED_PERC = "filtered_insert_joined_pct"

FILTER_READ_INSERT_SIZE_AMBIGUOUS_NUM = "FILTER_READ_INSERT_SIZE_AMBIGUOUS_NUM"
FILTER_READ_INSERT_SIZE_AMBIGUOUS_PERC = "FILTER_READ_INSERT_SIZE_AMBIGUOUS_PERC"
FILTER_READ_INSERT_SIZE_NO_SOLUTION_NUM = "FILTER_READ_INSERT_SIZE_NO_SOLUTION_NUM"
FILTER_READ_INSERT_SIZE_NO_SOLUTION_PERC = "FILTER_READ_INSERT_SIZE_NO_SOLUTION_PERC"
FILTER_READ_INSERT_SIZE_TOO_SHORT_NUM = "FILTER_READ_INSERT_SIZE_TOO_SHORT_NUM"
FILTER_READ_INSERT_SIZE_TOO_SHORT_PERC = "FILTER_READ_INSERT_SIZE_TOO_SHORT_PERC"

FILTER_READ_INSERT_SIZE_AVG_INSERT2 = "filtered_insert_mean"
FILTER_READ_INSERT_SIZE_STD_INSERT2 = "filtered_insert_stdev"
FILTER_READ_INSERT_SIZE_MODE_INSERT2 = "filtered_insert_mode"

FILTER_READ_INSERT_SIZE_INSERT_RANGE_START = "FILTER_READ_INSERT_SIZE_INSERT_RANGE_START"
FILTER_READ_INSERT_SIZE_INSERT_RANGE_END = "FILTER_READ_INSERT_SIZE_INSERT_RANGE_END"

FILTER_READ_INSERT_SIZE_90TH_PERC = "FILTER_READ_INSERT_SIZE_90TH_PERC"
FILTER_READ_INSERT_SIZE_50TH_PERC = "FILTER_READ_INSERT_SIZE_50TH_PERC"
FILTER_READ_INSERT_SIZE_10TH_PERC = "FILTER_READ_INSERT_SIZE_10TH_PERC"
FILTER_READ_INSERT_SIZE_75TH_PERC = "FILTER_READ_INSERT_SIZE_75TH_PERC"
FILTER_READ_INSERT_SIZE_25TH_PERC = "FILTER_READ_INSERT_SIZE_25TH_PERC"

GC_DIVERGENCE_CSV_FILE = "GC_DIVERGENCE_CSV_FILE"
GC_DIVERGENCE_PLOT_FILE = "GC_DIVERGENCE_PLOT_FILE"
GC_DIVERGENCE_COEFFICIENTS_CSV_FILE = "GC_DIVERGENCE_COEFFICIENTS_CSV_FILE"

GC_DIVERGENCE_COEFF_R1_AT = "GC_DIVERGENCE_COEFF_R1_AT"
GC_DIVERGENCE_COEFF_R1_ATCG = "GC_DIVERGENCE_COEFF_R1_ATCG"
GC_DIVERGENCE_COEFF_R1_CG = "GC_DIVERGENCE_COEFF_R1_CG"
GC_DIVERGENCE_COEFF_R2_AT = "GC_DIVERGENCE_COEFF_R2_AT"
GC_DIVERGENCE_COEFF_R2_ATCG = "GC_DIVERGENCE_COEFF_R2_ATCG"
GC_DIVERGENCE_COEFF_R2_CG = "GC_DIVERGENCE_COEFF_R2_CG"

FILTER_20MER_UNIQUENESS_TEXT = "read 20mer uniqueness text"
FILTER_20MER_UNIQUENESS_PLOT = "read 20mer uniqueness plot"
FILTER_20MER_UNIQUENESS_D3_HTML_PLOT = "ILLUMINA_READ_20MER_UNIQUENESS_D3_HTML_PLOT"
FILTER_20MER_SAMPLE_SIZE = "read 20mer sample size"
FILTER_20MER_PERCENTAGE_STARTING_MERS = "read 20mer percentage starting mers"
FILTER_20MER_PERCENTAGE_RANDOM_MERS = "read 20mer percentage random mers"

STATS_LIST_FILE_NAME = 'filterStats.txt'   ## the key=value file, need to be matched in pipeline config
FILE_LIST_FILE_NAME = 'file-list.txt'      ## the key=filename file, need to be matched in pipeline config

FILTER_METHODS_TXT = {"DNA": "dna.txt",
                      "FUNGAL": "fungal.txt",
                      "METAGENOME": "metagenome.txt",
                      "VIRAL-METAGENOME": "viral-metagenome.txt",
                      "ISO": "iso.txt",
                      "SAG": "sag.txt",
                      "CELL-ENRICHMENT": "cell-enrichment.txt",
                      "PLANT-2X150": "plant-2x150.txt",
                      "PLANT-2X250": "plant-2x250.txt",
                      "RNA": "rna.txt",
                      "RNAWOHUMAN": "rnawohuman.txt",
                      "3PRIMERNA": "3primerna.txt",
                      "SMRNA": "smrna.txt",
                      "METATRANSCRIPTOME": "mtaa.txt",
                      "MTF": "mtaa.txt",
                      "LFPE": "lfpe.txt",
                      "CLRS": "clrs.txt",
                      "CLIP-PE": "clip-pe.txt",
                      "NEXTERA-LMP": "nextera-lmp.txt",
                      "NEXTERA": "nextera-lmp.txt",
                      "NEXTSEQ": "nextseq.txt",
                      "ITAG": "itag.txt",
                      "MICROTRANS": "microtrans.txt",
                      "BISULPHITE": "bisulphite.txt",
                      "CHIPSEQ": "chip-seq.txt"}

FILTER_METADATA_JSON = {"METATRANSCRIPTOME": "filter_mt.json",
                        "MTF": "filter_mt.json",
                        "NEXTERA-LMP": "filter_nextera.json",
                        "NEXTERA": "filter_nextera.json"}

# unique mer analysis
MER_SAMPLE_MER_SIZE = 25
MER_SAMPLE_REPORT_FRQ = 25000


## Constants for reporting
class FilterQcStats:
    FILTER_READ_MATCHING_HITS = "read matching hits of"
    FILTER_READ_PARSED_FILE = "read parsed file of"
    FILTER_READ_TOP_HITS = "read top hits of"
    FILTER_READ_TOP_100_HITS = "read top 100 hits of"
    FILTER_READ_TAX_SPECIES = "read tax species of"
    FILTER_READ_TAXLIST_FILE = "read taxlist file of"
    FILTER_READ_TOPHIT_FILE = "read tophit file of"
    FILTER_READ_TOP100HIT_FILE = "read top100hit file of"
    FILTER_READ_KHIST_LOG = "khist log file"
    FILTER_READ_KHIST_TXT = "khist txt file"

## EOF
