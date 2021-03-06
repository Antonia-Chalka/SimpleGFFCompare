"""
Created on 14 May 2019
@author: Antonia Chalka
"""

import os
import argparse
import re

# Setting up Arguments
parser = argparse.ArgumentParser(description='Compare two GFF files. Currently checks only CDS between files.')
parser.add_argument('SampleFileName', metavar='S',
                    help='Name of sample GFF file')
parser.add_argument('TemplateFileName', metavar='T',
                    help='Name of template GFF file')
parser.add_argument('-i', '--IterationNumber', metavar='I', type=int, default=10,
                    help='Max times to run fuzzy filter. Each iteration increases range by 3. Default: 10')
parser.add_argument('-m', '--mismatch_permissibility', metavar='M', type=int, default=3,
                    help='The max length of mismatches that will NOT be logged in mismatch file. Default: 3')
parser.add_argument('-p' '--GeneNameListFile', metavar='P', default=str(os.getcwd() + "/" + "proteinList.txt"),
                    help='Directory of file with gene names (1 per line) to be used in generating the counts table. ' +
                         'Default is proteinList.txt.')
argument = parser.parse_args()
IterationNumber = argument.IterationNumber 
mismatch_permissibility = argument.mismatch_permissibility

# Return real, absolute, file path (eg ~/directory -> /home/user/directory)
def make_path_sane(path):
    path = os.path.expanduser(path)
    path = os.path.normpath(path)
    path = os.path.realpath(path)
    path = os.path.abspath(path)
    return path

SampleFile = open(make_path_sane(argument.SampleFileName), "r")
TemplateFile = open(make_path_sane(argument.TemplateFileName), "r")

# Set up custom names for output files: strainName_mismatch.txt
<<<<<<< HEAD
StrainName = os.path.basename(argument.TemplateFileName)
LogFile = open(StrainName + "_" + "Log.txt", "w")
MismatchFile = open(StrainName + "_" + "Mismatch.txt", "w")
StrainStatsFile = open(StrainName + "_" + "StrainStats.csv", "w")
=======
os.makedirs(os.path.dirname("./GFFCompare_Output/"), exist_ok=True) # Make output folder if i does not exist
StrainName = os.path.basename(argument.SampleFileName)
LogFile = open("./GFFCompare_Output/" + StrainName + "_" + "Log.txt", "w")
MismatchFile = open("./GFFCompare_Output/" + StrainName + "_" + "Mismatch.txt", "w")
StrainStatsFile = open("./GFFCompare_Output/" + StrainName + "_" + "StrainStats.txt", "w")
ProteinOutputFile = open("./GFFCompare_Output/" + StrainName + "_" + "GeneCounts.txt", "w")
MismatchStatsOutputFile = open("./GFFCompare_Output/" + StrainName + "_" + "MismatchStats.csv", "w")

# Stats (False Positive, False Negative, True Positive)
def calc_stats_lenght(sample_start, sample_stop, template_start, template_stop):
    FP=0
    FN=0
    # mismatch = sample length - template lenght
    mismatch = ((sample_stop - sample_start) + 1) - ((template_stop - template_start) + 1)
    if mismatch < 0: # template longer
        FN = abs(mismatch)
    elif mismatch > 0: # sample longer
        FP = mismatch

    TP = (template_stop - template_start) - FN
    return FP, FN, TP

def calc_stats_position(sample_start, sample_stop, template_start, template_stop):
    FP = 0
    FN = 0
    
    start_mismatch = template_start - sample_start
    if start_mismatch < 0: # template longer
        FN  += abs(start_mismatch)
    elif start_mismatch > 0: # sample longer
        FP += start_mismatch

    stop_mismatch = sample_stop - template_stop
    if stop_mismatch < 0: # template longer
        FN += abs(stop_mismatch)
    elif stop_mismatch > 0: # sample longer
        FP += stop_mismatch
    #TP = template lenght - FN total
    TP = ((template_stop - template_start) +1) - FN
    return FP, FN, TP

FP_length = 0
FN_length = 0
TP_length = 0
FP_position = 0
FN_position = 0
TP_position = 0
>>>>>>> e3350d9ecddf1109638f6779f37f9970da81c0b5

# Reading Sample file
LogFile.write("Reading " + argument.SampleFileName + " as sample file... \n")
SampleCDSs = []
SamplemRNAs = []
SampleGenes = []
SampleMiscFeatures = []
SampleOther = []

with SampleFile:
    for line in SampleFile:
        if (line.startswith("##") is False) and (not line.strip() == ''):  # Do not process comment or empty lines
            fields = tuple(line.split("\t"))
            # Assign each feature info to appropriate list (can be modified)
            if fields[2] == "CDS":
                SampleCDSs.append(fields)
            elif fields[2] == "mRNA":
                SamplemRNAs.append(fields)
            elif fields[2] == "gene":
                SampleGenes.append(fields)
            elif fields[2] == "misc_feature":
                SampleMiscFeatures.append(fields)
            else:
                SampleOther.append(fields)
SampleFile.close()

# Log stats of sample file
LogFile.write("GGF Sample file read. \n")
sample_total_lines = len(SampleCDSs) + len(SamplemRNAs) + len(SampleGenes) + len(SampleMiscFeatures) + len(SampleOther)
LogFile.write("Total lines: " + str(sample_total_lines) + "\n" +
              "CDC lines: " + str(len(SampleCDSs)) + "\n" +
              "mRNA lines: " + str(len(SamplemRNAs)) + "\n" + 
              "Gene lines: " + str(len(SampleGenes)) + "\n" +
              "Misc. Features lines: " + str(len(SampleMiscFeatures)) + "\n" +
              "Other lines: " + str(len(SampleOther)) + "\n \n")

# Reading template file
LogFile.write("Reading " + argument.TemplateFileName + " as template file... \n")

TemplateCDSs = []
TemplatemRNAs = []
TemplateGenes = []
TemplateMiscFeatures = []
TemplateOther = []

with TemplateFile:
    for line in TemplateFile:
        if (line.startswith("##") is False) and (not line.strip() == ''):  # Do not process comment or empty lines
            fields = tuple(line.split("\t"))
            # Assign each feature info to appropriate list (can be modified)
            if fields[2] == "CDS":
                TemplateCDSs.append(fields)
            elif fields[2] == "mRNA":
                TemplatemRNAs.append(fields)
            elif fields[2] == "gene":
                TemplateGenes.append(fields)
            elif fields[2] == "misc_feature":
                TemplateMiscFeatures.append(fields)
            else:
                TemplateOther.append(fields)
TemplateFile.close()

# Log stats of template file
LogFile.write("GFF Template File read.\n")
template_total_lines = len(TemplateCDSs) + len(TemplatemRNAs) + len(TemplateGenes) + len(TemplateMiscFeatures) + len(TemplateOther)
LogFile.write("Total lines: " + str(template_total_lines) + "\n" +
              "CDC lines: " + str(len(TemplateCDSs)) + "\n" +
              "mRNA lines: " + str(len(TemplatemRNAs)) + "\n" +
              "Gene lines: " + str(len(TemplateGenes)) + "\n" +
              "Misc. Features lines: " + str(len(TemplateMiscFeatures)) + "\n" +
              "Other lines: " + str(len(TemplateOther)) + "\n \n")

# Comparison & Detection of mismatches (start position is 3 in list, stop is 4)
LogFile.write("Beginning comparison of CDS features..." + "\n")
# For stats later:
perfect_matches = 0
single_matches = 0
fuzzy_matches = 0
no_matches = 0
<<<<<<< HEAD

true_positives = 0
false_negatives = 0
false_positives = 0

=======
>>>>>>> e3350d9ecddf1109638f6779f37f9970da81c0b5
mismatches = []  # List to put lists of possible matches in and filter later on

for CDS1 in SampleCDSs:  # Iterate for each cds in sample file
    LogFile.write("\nSearching for matching Sample CDS: " + str(CDS1) + "\n")

    found_pair = False  # Assume no match at beginning
    mismatch = [CDS1]  # List to store original mismatches

    for CDS2 in TemplateCDSs:
        if (CDS1[3] == CDS2[3]) and (CDS1[4] == CDS2[4]):  # Perfect match
            LogFile.write("Perfect match with Template CDS: " + str(CDS2) + "\n" +
                          "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" + 
                          "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
            found_pair = True  # To skip all other mismatch checks
            # TODO Check if calculated correctly
            true_positives += (CDS1[4] - CDS1[3])

            mismatch = []  # Empty mismatch list
            perfect_matches += 1
            
            FP, FN, TP = calc_stats_lenght(int(CDS1[3]), int(CDS1[4]), int(CDS2[3]), int(CDS2[4]))
            FP_length += FP
            FN_length += FN
            TP_length += TP

            FP, FN, TP = calc_stats_position(int(CDS1[3]), int(CDS1[4]), int(CDS2[3]), int(CDS2[4]))
            FP_position += FP
            FN_position += FN
            TP_position += TP

            break

    # If no perfect match, search for single start/stop perfect matches.
    # Allow for multiple possible matches to be detected (no breaks)
    if found_pair is False:
        LogFile.write("No exact match for both start and stop positions found." +
                      "Iterating for start or stop only matches...\n")
        for CDS2 in TemplateCDSs:        
            if CDS1[3] == CDS2[3] and CDS1[4] != CDS2[4]:  # Mismatch at stop
                LogFile.write("Found match with starting position: " + str(CDS2) + "\n" +
                              "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" +
                              "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
                mismatch.append(CDS2)
                found_pair = True
            elif CDS1[3] != CDS2[3] and CDS1[4] == CDS2[4]:  # Mismatch at start
                LogFile.write("Found match with ending position: " + str(CDS2) + "\n" +
                              "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" +
                              "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
                mismatch.append(CDS2)
                found_pair = True
        if found_pair is True:
            single_matches += 1

            
    if found_pair is False:  # No perfect matches for either start and/or stop positions
        LogFile.write("No exact match for either start and stop positions found. Implementing fuzzy filter...\n")
        
        for i in range(1, IterationNumber+1):  # Iterate filter to user-specified number (default=10)
            nts = 3 * i  # Mismatch length based on filter iteration
            LogFile.write("Filter allowance :" + str(i) + " codons (" + str(nts) + " nucleotides)" + "\n")
                
            for CDS2 in TemplateCDSs:
                # Check if CDS2 starting position is within range
                StartMatch = int(CDS1[3])-nts <= int(CDS2[3]) <= int(CDS1[3]) + nts
                # Check if CDS2 stopping position is within range
                EndMatch = int(CDS1[4])-nts <= int(CDS2[4]) <= int(CDS1[4]) + nts
                
                if StartMatch is True and EndMatch is True:
                    LogFile.write("Found start and stop position within range, with CDS: " + str(CDS2) + "\n" +
                                  "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" + 
                                  "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
                    mismatch.append(CDS2)
                    found_pair = True
                    
                elif StartMatch is True and EndMatch is False:
                    LogFile.write("Found matching start position within range, with CDS: " + str(CDS2) + "\n" +
                                  "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" + 
                                  "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
                    mismatch.append(CDS2)
                    found_pair = True
                    
                elif StartMatch is False and EndMatch is True:
                    LogFile.write("Found matching stop position within range, with CDS: " + str(CDS2) + "\n" +
                                  "Start positions: " + str(CDS1[3]) + " " + str(CDS2[3]) + "\n" + 
                                  "Stop positions: " + str(CDS1[4]) + " " + str(CDS2[4]) + "\n")
                    mismatch.append(CDS2)
                    found_pair = True

            if found_pair is True:
                fuzzy_matches += 1
                break
    # if others cds have been added to initial list: [CDS1, CDS2] /[CDS1, CDS2, CDS2,...] -> len >1).
    # No matches (len=1), have already been added to mismatch file in above if statement
    if len(mismatch) > 1:
        mismatches.append(mismatch)

    if found_pair is False:  # No matches found
        LogFile.write("Could not find any match, even with fuzzy filter. Moving on to next CDS in file...\n")
        MismatchFile.write("WARNING: No matches found for:" + str(CDS1) + ".\n \n") 
        no_matches += 1

        
# Matching stats
LogFile.write("\nComparison has been completed.\n" +
              "Perfect matches (both positions):" + str(perfect_matches) + "\n" +
              "Perfect match (one position only):" + str(single_matches) + "\n" +
              "Matches (fuzzy filter):" + str(fuzzy_matches) + "\n" +
              "No matches:" + str(no_matches) + "\n \n")


# Filter mismatches
LogFile.write("Filtering mismatches...\n")
MismatchFile.write("Logging mismatches bigger than " + str(mismatch_permissibility) + "\n")
MismatchesLogged = 0
MismatchesDiscarded = 0
TotalMismatchLength = 0
LoggedMismatches = []  # To be used in gene count table output

for mismatch in mismatches:
    # CDS1 = mismatch[0]
    # CDS2 = mismatch[1+]
    log_mismatch = True  # Whether to log mismatch in file, assume yes as default
    CDS1 = mismatch[0]  # for ease of reference
    LogFile.write("\nFiltering mismatches for " + str(CDS1) + "\n")
    
    for i in range(1, len(mismatch)):
        CDS2 = mismatch[i]  # for ease of reference
        mismatch_length_start = abs(int(CDS1[3]) - int(CDS2[3]))
        mismatch_length_stop = abs(int(CDS1[4]) - int(CDS2[4]))

        LogFile.write("Checking mismatch with " + str(CDS2) + "\n" +
                      "Start position mismatch length: " + str(mismatch_length_start) + "\n"
                      "Stop position mismatch length: " + str(mismatch_length_stop) + "\n")

        # if mismatch is below threshold, do not log to mismatch file
        if (mismatch_length_start <= mismatch_permissibility) and (mismatch_length_stop <= mismatch_permissibility):
            # TODO ADD TO TP/FP/FN totals

            log_mismatch = False
            LogFile.write("Mismatches lower than threshold (" + str(mismatch_permissibility) + "). Discarding...\n")
            
    if log_mismatch is True:  # Start logging mismatches above threshold
        MismatchesLogged += 1
        LoggedMismatches.append(CDS1)
        LogFile.write("Mismatches bigger than threshold (" + str(mismatch_permissibility) +
                      "). Logging to mismatch file...\n")
        MismatchFile.write("\nMismatch Origin\t" + CDS1[8].split(";")[0] + "\t" + CDS1[3] + "\t" + CDS1[4] + "\t\n")
        for i in range(1, len(mismatch)):
            CDS2 = mismatch[i]
            mismatch_length_start = int(CDS1[3]) - int(CDS2[3])
            mismatch_length_stop = int(CDS1[4]) - int(CDS2[4])
            # TODO TEST IF FP/FN WORKS
            if mismatch_length_start > 0:  # negative
                false_negatives += mismatch_length_start
            elif mismatch_length_start < 0:  # positive
                false_positives += mismatch_length_start

            if mismatch_length_stop > 0:  # negative
                false_negatives += mismatch_length_stop
            elif mismatch_length_stop < 0:  # positive
                false_positives += mismatch_length_stop
            # TODO IMPLEMENT TP
            # true_positives =

            MismatchFile.write("Possible Match\t" + CDS2[8].split(";")[0] + "\t" + CDS2[3] + "\t" + CDS2[4] + "\t\n" +
                               "Mismatch Length \tStart: " + str(mismatch_length_start) + "\tStop: " +
                               str(mismatch_length_stop) + "\n")
            TotalMismatchLength += mismatch_length_start + mismatch_length_stop

            # TODO FP FN AND TP
            FP, FN, TP = calc_stats_lenght(int(CDS1[3]), int(CDS1[4]), int(CDS2[3]), int(CDS2[4]))
            FP_length += FP
            FN_length += FN
            TP_length += TP

            FP, FN, TP = calc_stats_position(int(CDS1[3]), int(CDS1[4]), int(CDS2[3]), int(CDS2[4]))
            FP_position += FP
            FN_position += FN
            TP_position += TP


    else:
        MismatchesDiscarded += 1

# Stats for mismatch filtering
LogFile.write("\nTotal mismatches: " + str(MismatchesDiscarded + MismatchesLogged) + "\n" +
              "Mismatches Discarded (below " + str(mismatch_permissibility) + " threshold ): " +
              str(MismatchesDiscarded) + "\n" +
              "Mismatches logged: " + str(MismatchesLogged) + "\n")

# Strain Stats
StrainStatsFile.write("Strain\t" +
                      "Num_CDS_Automated\t" +
                      "Num_CDS_Deposited\t" +
                      "Num_Genes_Automated\t" +
                      "Num_Genes_Deposited\t" +
                      "Num_mRNA_Automated\t" +
                      "Num_mRNA_Deposited\t" +
                      "Other/MiscFeatures_Automated\t" +
                      "Other/MiscFeatures_Deposited\t" +
                      "Total_Features_Automated\t" + 
                      "Total_Features_Deposited\t" +
                      "Perfect_Matches\t" + 
                      "No_Matches\t" +
                      "Total_Mismatches\t" +
                      "Mismatches_Discarded\t" +
                      "Mismatches_Logged\t" +
                      "Total_Mismatch_length(Logged_Only)\t" +
                      "Avg_Mismatch_Length\n"
                      )
def avg_mismatch_length_calc(n, d):
    return n / d if d else 0
StrainStatsFile.write(StrainName + "\t" +
                      str(len(SampleCDSs)) + "\t" +
                      str(len(TemplateCDSs)) + "\t" +
                      str(len(SampleGenes)) + "\t" +
                      str(len(TemplateGenes)) + "\t" +
                      str(len(SamplemRNAs)) + "\t" +
                      str(len(TemplatemRNAs)) + "\t" +
                      str(len(SampleMiscFeatures) + len(SampleOther)) + "\t" +
                      str(len(TemplateMiscFeatures) + len(TemplateOther)) + "\t" +
                      str(sample_total_lines) + "\t" +
                      str(template_total_lines) + "\t" +
                      str(perfect_matches) + "\t" +
                      str(no_matches) + "\t" +
                      str(single_matches + fuzzy_matches) + "\t" +
                      str(MismatchesDiscarded) + "\t" +
                      str(MismatchesLogged) + "\t" +
                      str(TotalMismatchLength) + "\t" +
                      str(avg_mismatch_length_calc(TotalMismatchLength, MismatchesLogged)) + "\n"
                      )

# Protein Table Output. Populate protein list with protein names. Does NOT check for duplicates
LogFile.write("\nBeginning gene counter...\n")
LogFile.write("Getting gene names from file...\n")# 
GeneListFile = open(os.path.dirname(os.path.abspath(__file__)) + "/proteinList.txt", "r")

GeneNames = []
with GeneListFile:
    for line in GeneListFile:
        GeneNames.append(line.strip())

# Counter method
def feature_counter(feature_list, count_gene):
    counter = 0
    for feature in feature_list:
        '''Regex explanation:
        gene : the gene to search for, eg RL1
        (?<![A-Z]|[a-z]|[0-9]| ): make sure that the gene is not *preceded* by letters or characters, eg dont match 
        with 1RL1 or ARL1. Space is prevent 'comment' fields to be captured (see below)
        (?<!Note=): make sure match is not preceded by 'Note=' in gff file, as some ncbi files include the 'parent' 
        gene which can create fake duplicate pattern matching and inflate count of parent gene.
        (?![A-Z]|[a-z]|[0-9]): as above, make sure that the gene is not *followed* by letters or characters, 
        eg dont match with RL11 or RL1A
        
        Overall when searching for RL1, the following should NOT match:
        1RL1; ARL1; ;RL11 ;RL1A Note=RL1_ BRL1A etc
        The following should match:
        ;RL1_ =RL1_ =RL1;
        '''
        regex = "(?<![A-Z]|[a-z]|[0-9]| )(?<!Note=)" + count_gene + "(?![A-Z]|[a-z]|[0-9]| )"
        if re.search(regex, feature[8]) is not None:  # if no hit, will return None
            counter += 1
    return counter

# File Writing
LogFile.write("Writing headers to gene counter output file...\n")
ProteinOutputFile.write("Gene\t" +  # Headers
                        "Num_CDS_Automated\t" +
                        "Num_CDS_Deposited\t" +
                        "Num_Genes_Automated\t" +
                        "Num_Genes_Deposited\t" +
                        "Num_mRNA_Automated\t" +
                        "Num_mRNA_Deposited\t" +
                        "Other/MiscFeatures_Automated\t" +
                        "Other/MiscFeatures_Deposited\t" +
                        "Logged_Mismatches(Sample)\n")
for gene in GeneNames:
    LogFile.write("Writing " + gene + " counter...\n")
    ProteinOutputFile.write(gene + "\t" +
                            str(feature_counter(SampleCDSs, gene)) + "\t" +
                            str(feature_counter(TemplateCDSs, gene)) + "\t" +
                            str(feature_counter(SampleGenes, gene)) + "\t" +
                            str(feature_counter(TemplateGenes, gene)) + "\t" +
                            str(feature_counter(SamplemRNAs, gene)) + "\t" +
                            str(feature_counter(TemplatemRNAs, gene)) + "\t" +
                            str(feature_counter(SampleMiscFeatures, gene) +
                                feature_counter(SampleOther, gene)) + "\t" +
                            str(feature_counter(TemplateMiscFeatures, gene) +
                                feature_counter(TemplateOther, gene)) + "\t" +
                            str(feature_counter(LoggedMismatches, gene)) + "\n")

# Mismatch Output File
MismatchStatsOutputFile.write("Strain\t" +  # Headers
                        "FP_Lenght\t" +
                        "FN_Lenght\t" +
                        "TP_Lenght\t" +
                        "Sensitivity_Length\t" +
                        "FDR_Length\t" +
                        "FP_Position\t" +
                        "FN_Position\t" +
                        "TP_Position\t" +
                        "Sensitivity_Position\t" +
                        "FDR_Position\t" + "\n")

MismatchStatsOutputFile.write(StrainName + "\t" +
                            str(FP_length) + "\t" +
                            str(FN_length) + "\t" +
                            str(TP_length) + "\t" +
                            str(TP_length / (TP_length + FN_length)) + "\t" +
                            str(FP_length / (FP_length + TP_length)) + "\t" +
                            str(FP_position) + "\t" +
                            str(FN_position) + "\t" +
                            str(TP_position) + "\t" +
                            str(TP_position / (TP_position + FN_position)) + "\t" +
                            str(FP_position / (FP_position + TP_position)) + "\n")


# File handling
ProteinOutputFile.close()
LogFile.close()
MismatchFile.close()
MismatchStatsOutputFile.close()

# $filePath = Get-ChildItem "C:\Users\s1438773\Onedrive - University of Edinburgh\Msc\GFFCompare\data\resultGracy" -Filter *.fasta
# cd  "C:\Users\s1438773\Onedrive - University of Edinburgh\Msc\GFFCompare\"
#  ForEach($file in $filePath){
#    Write-Host $file
#    python GFFCompare.py .\data\resultGRACy\$file .\data\NCBI\NCBI_GFF\$file -m 0
#}
# copy *.csv merged.csv 




