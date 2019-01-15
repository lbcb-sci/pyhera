#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

import paramsparser

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))

from fastqparser import read_fastq
import utility_sam

# Parameter definitions for paramparser
paramdefs = {'--version' : 0,
             '-v' : 0}


# Setting run names for various tools used for analysisw
MINIMAP2 = os.path.join(SCRIPT_PATH, 'minimap2', 'minimap2')
SAMSCRIPTS = os.path.join(SCRIPT_PATH, 'samscripts')
FASTQFILTER = os.path.join(SAMSCRIPTS, 'src', 'fastqfilter.py')
GEPARD_JAR = '/home/kkrizanovic/Downloads/Gepard-1.40.jar'
GEPARD_MATRIX = '/home/kkrizanovic/src/gepard/resources/matrices/edna.mat'
PYTHON = 'python'

default_MM2options = '-ax map-pb'

def check_tools():
    global SCRIPT_PATH, BRIDGER, MINIMAP2, EZRA, PYTHON
    if not os.path.exists(SCRIPT_PATH):
        sys.stderr.write('\nChecking tools: folder %s does not exist!\n' % SCRIPT_PATH)
        return False
    elif not os.path.exists(MINIMAP2):
        sys.stderr.write('\nChecking tools: Minmap2 executable (%s) does not exist!\n' % MINIMAP2)
        return False
    elif not os.path.exists(FASTQFILTER):
        sys.stderr.write('\nChecking tools: FastqFilter script (%s) does not exist!\n' % FASTQFILTER)
        return False
    elif not os.path.exists(GEPARD_JAR):
        sys.stderr.write('\nChecking tools: Gepard dot plotter jar file (%s) does not exist!\n' % GEPARD_JAR)
        return False

    (status, output) = commands.getstatusoutput(PYTHON + ' --verion')
    if not output.startswith('Python 2.7'):
        PYTHON = 'python2'
        (status, output) = commands.getstatusoutput(PYTHON + ' --version')
        if not output.startswith('Python 2.7'):
            sys.stderr.write('\nThis script requires python 2.7 to run! Cannot find appropriate python!\n')
            return False

    return True

def print_version():
    sys.stdout.write('\nScaffold analysis script, version 1.0\n');

# Analyze scaffolds by comparing them to the reference, in several steps
# 1. Map all scaffolds against reference using Minimap2. Produce SAM (TODO: PAF) to easily load using samscripts
# 2. Print scaffolds that do not map to any part of reference and print scaffolds that map to different parts of the reference
# 3. Run Gepard on all pairs of scaffold-reference fro which a valid mapping exists
def scara_analyze(scaffolds_file, reference_file, output_folder):
    sys.stderr.write('\nSTARTING SCAFFOLDING ANALYSIS SCRIPT')

    output_folder_path = os.path.join(os.getcwd(), output_folder)

    ### STEP 0. Checking paths and folders
    if not os.path.exists(scaffolds_file):
        sys.stderr.write('\nScaffolds file does not exist (%s)! Exiting ...' % scaffolds_file)
        return
    elif not os.path.exists(reference_file):
        sys.stderr.write('\nReference file does not exist (%s)! Exiting ...' % reference_file)
        return
    elif not os.path.exists(output_folder):
        sys.stderr.write('\nOutput folder does not exist (%s)! Creating it ...' % output_folder)
        os.mkdir(output_folder_path)

    ### STEP 1. Running Minimap2
    sys.stderr.write('\nCALCULATING MAPPINGS BETWEEN SCAFFOLDS AND REFERENCE!')
    minimap2_output_file = os.path.join(output_folder_path, 'scaffolds2reference.sam')
    if os.path.exists(minimap2_output_file):
        sys.stderr.write('\nMapping file already present! Skipping ...!')
    else:
        cmd = '%s %s %s %s > %s' % (MINIMAP2, default_MM2options, reference_file, scaffolds_file, minimap2_output_file)
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)
        logfile = os.path.join(output_folder_path, 'Minimap2_r2r.log')
        with open(logfile, 'w') as lfile:
            lfile.write(output)

    ### STEP 2. Load and analyze Minimap2 file
    # Loading SAM file into a dictionary
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    sys.stderr.write('\nANALYZING MAPPINGS!')
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(minimap2_output_file, qnames_with_multiple_alignments)
    
    # Load scaffolds
    [theaders, sseqs, squals] = read_fastq(scaffolds_file)
    # Cutting headers at first space
    sheaders = []
    for theader in theaders:
        sheader = theader[:theader.find(' ')]
        sheaders.append(sheader)
    
    # Load reference
    [theaders, rseqs, rquals] = read_fastq(reference_file)
    # Cutting headers at first space
    rheaders = []
    for theader in theaders:
        rheader = theader[:theader.find(' ')]
        rheaders.append(rheader)
    
    # Analyze SAM
    scaffold_mappings = {}      # A dictionary that for each scaffold that is mapped to a reference contains
                                # a list of reference parts (chromosome) to which the scaffold is mapped
    reference_mappings = {}     # A dictionary that for each reference that is mapped to a scaffold contains
                                # a list of scaffolds mapped to it
    for sheader in sheaders:
        scaffold_mappings[sheader] = []
    for rheader in rheaders:
        reference_mappings[rheader] = []
    for (qname, sam_lines) in sam_hash.iteritems():
        for samline in sam_lines:
            # Skip samlines with invalid CIGAR
            if samline.cigar == '*':
                continue
            if qname in sheaders:
                sname = qname
                rname = samline.rname
            elif qname in rheaders:
                sname = samline.rname
                rname = qname
            else:
                sys.stderr.write('\nERROR: Invalid query name in mappings file (%s)!' % qname)
                return

            smappings = scaffold_mappings[sname]
            if rname not in smappings:
                scaffold_mappings[sname].append(rname)
            
            rmappings = reference_mappings[rname]
            if sname not in rmappings:
                reference_mappings[rname].append(sname)

    # Print scaffold-reference mappings (both dictionaries)
    found_double_mappings = False
    found_zero_mappings = False
    mapping_analysis_file = os.path.join(output_folder_path, 'mapping_analysis.txt')
    with open(mapping_analysis_file, 'w') as mafile:
        mafile.write('SCAFFOLD: REFERENCE LIST\n')
        for sname, rname_list in scaffold_mappings.iteritems():
            if len(rname_list) > 1:
                found_double_mappings = True
            if len(rname_list) == 0:
                found_zero_mappings = True
            mafile.write('%s: %s\n' % (sname, ', '.join(rname_list)))

        mafile.write('REFERENCE: SCAFFOLD LIST\n')
        for rname, sname_list in reference_mappings.iteritems():
            mafile.write('%s: %s\n' % (rname, ', '.join(sname_list)))

    if found_double_mappings:
        sys.stderr.write('\nWARNING: Found scaffolds mapped to multiple references!')
    if found_zero_mappings:
        sys.stderr.write('\nWARNING: Found unmapped scaffolds!')

    ### STEP 3. Generate gepard dot plots for all mappings between scaffolds and references
    sys.stderr.write('\nGENERATING DOT PLOTS FOR SCAFFOLDS!')
    # Create separate fasta files for each scaffold
    scaffolds_folder = os.path.join(output_folder_path, 'scaffolds')
    os.mkdir(scaffolds_folder)
    for i in xrange(len(sheaders)):
        sheader =sheaders[i]
        sseq = sseqs[i]
        sfilename = os.path.join(scaffolds_folder, sheader + '.fasta')
        with open(sfilename, 'w') as sfile:
            sfile.write('>%s\n%s\n' % (sheader, sseq))

    # Create separate fasta file for each reference
    referencess_folder = os.path.join(output_folder_path, 'references')
    os.mkdir(referencess_folder)
    for i in xrange(len(rheaders)):
        rheader = rheaders[i]
        rseq = rseqs[i]
        rfilename = os.path.join(referencess_folder, rheader + '.fasta')
        with open(rfilename, 'w') as rfile:
            rfile.write('>%s\n%s\n' % (rheader, rseq))

    # Generate dot plots
    gepard_folder = os.path.join(output_folder_path, 'scaff_gepard')
    os.mkdir(gepard_folder)
    for sname, rname_list in scaffold_mappings.iteritems():
        for rname in rname_list:
            sfilename = os.path.join(scaffolds_folder, sname + '.fasta')
            if not os.path.exists(sfilename):
                sys.stderr.write('\nERROR: Scaffold fasta file not found: %s' % sfilename)
            rfilename = os.path.join(referencess_folder, rname + '.fasta')
            if not os.path.exists(rfilename):
                sys.stderr.write('\nERROR: Reference fasta file not found: %s' % rfilename)
            gepard_file = os.path.join(gepard_folder, '%s_%s.png' % (sname, rname))
            cmd = 'java -cp %s org.gepard.client.cmdline.CommandLine -seq1 %s -seq2 %s -matrix %s -outfile %s' \
                % (GEPARD_JAR, sfilename, rfilename, GEPARD_MATRIX, gepard_file)
            sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
            (status, output) = commands.getstatusoutput(cmd)


    ### STEP 3.1 Generate additional dotplots, one for each reference against all scaffolds mapped to it
    sys.stderr.write('\nGENERATING DOT PLOTS FOR REFERENCES!')
    gepard_folder2 = os.path.join(output_folder_path, 'ref_gepard')
    os.mkdir(gepard_folder2)
    for rname, sname_list in reference_mappings.iteritems():
        if len(sname_list) > 0:
            sfilename = os.path.join(gepard_folder2, '%s_scaffolds.fasta' % rname)
            rfilename = os.path.join(referencess_folder, rname + '.fasta')
            with open(sfilename, 'w') as sfile:
                for i in xrange(len(sheaders)):
                    sheader =sheaders[i]
                    sseq = sseqs[i]
                    if sheader in sname_list:
                        sfile.write('>%s\n%s\n' % (sheader, sseq))
            gepard_file2 = os.path.join(gepard_folder2, '%s_scaffolds.png' % rname)
            cmd = 'java -cp %s org.gepard.client.cmdline.CommandLine -seq1 %s -seq2 %s -matrix %s -outfile %s' \
                % (GEPARD_JAR, sfilename, rfilename, GEPARD_MATRIX, gepard_file2)
            sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
            (status, output) = commands.getstatusoutput(cmd)
  
    return



def verbose_usage_and_exit():
    sys.stderr.write('\nscara_analyze - analyze the results of scaffolding by comparing them to the reference.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [scaffolds file] [reference file] [output folder]\n' % sys.argv[0])
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[1:])
        if '-v' in paramdict or '--version' in paramdict:
            print_version()
        verbose_usage_and_exit()

    scaffolds_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_folder = sys.argv[3]

    # Currently not needed, but setting up the papramsparser anyways
    pparser = paramsparser.Parser(paramdefs)
    paramdict = pparser.parseCmdArgs(sys.argv[4:])
    paramdict['command'] = ' '.join(sys.argv)

    scara_analyze(scaffolds_file, reference_file, output_folder)
