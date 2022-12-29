# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 09:53:22 2020

@author: zhangguoqing

E-mail: zhangguoqing84@westlake.edu.cn

"""
__author__ = ("Ju Feng",
              "Zhang Guoqing")
__version__ = '1.0.0'
__date__ = '2020.01.13'

import sys, os
import argparse as ap
import subprocess as subp

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

# The script path
pwdpath = os.path.dirname(os.path.abspath(__file__))

# Default database path
arg_db = "/".join([pwdpath, "ARGfams_v0.1.hmm"])
arg_large_db = "/".join([pwdpath, "Pfam-Tigrfam.hmm"])

# Parameter Info
def read_params(args):
    p = ap.ArgumentParser(description=
        "DESCRIPTION\n"
        "FunGeneTyper version: " + __version__ + "\n"
        "Detailed introducion \n"
        "Usage \n",
        formatter_class=ap.RawTextHelpFormatter,
        add_help=False)
    arg = p.add_argument

    arg('-i', '--input', dest='input_file', metavar='INPUT_FILE', type=str,
        required=True, default=None, help="the input file ORFs.faa")

    arg('-o', '--output', dest='output_file_name', metavar='OUTPUT_FILE_NAME',
        type=str, default='FunGeneTyperout', nargs='?',
        help='the outputfile prefix name: eg. PRIFEX.tlout')

    g = p.add_argument_group('Required arguments')
    arg = g.add_argument

#    arg('-t', '--type', dest='type', metavar='Analysis_Type ',
#        choices=['MGE', 'ARG', 'Both'], type=str, default='Both',
#        help='Selective analysis type, This tools provides three types of analysis\n'
#        )

    arg('-db', dest = 'sub_db', metavar='Sub_Database', type=str, 
        default=arg_db, help='Sub database; Default Antibiotic Resistance Genes')



    arg('-DB', dest='large_db', metavar='synthesis_database', type=str,
        default=arg_large_db, help='synthesis database, Default Antibiotic Resistance Genes database')

    arg(dest='hmmtype', nargs='?', choices=['--cut_ga', '--cut_nc', '--cut_tc'],
        default='--cut_ga', help='hmm type; chose from {--cut_ga, --cut_nc, --cut_tc} \n default:cut_ga')

    arg('-n', '--nproc', metavar="N", dest='cpucore', type=int, default=2,
        help="The number of CPUs to use for parallelizing the mapping [default 1]")

    g = p.add_argument_group('Other arguments')
    arg = g.add_argument

    arg('--check', action='store_true',
        help="Only checks if the Default ARG DB is installed and installs it if not.")
    arg('-v', '--version', action='version',
        version="FunGeneTyper version {} ({})".format(__version__, __date__),
        help="Prints the current MetaPhlAn version and exit")
    arg("-h", "--help", action="help", help="show this help message and exit")

    pars = p.parse_args()
    return pars


##################################################################################
### HMM mapping function
### inputfile: ORF file(faa)
### ouputfile: outputfile + '.hmm.tlout'
##################################################################################

def hmm_subdatabase(inputfile, outputfile, cpucore, database, hmmtype):
    try:
        subp.check_call('hmmscan ' + hmmtype + ' --cpu ' + str(cpucore) + ' --domtblout '
                        + outputfile + ' ' + database + ' ' + inputfile + ' >> protein_seq_c_ORFs.log', shell=True)
    except subp.CalledProcessError:
        print("hmmscan has some unexpected problems, please please check whether the parameters are accurate;\n \
        If you are sure that the parameters are correct, please check whether hmmscan is working properly in the analysis environment.")

def hmm_largedatabase(inputfile, outputfile, cpucore, database, hmmtype):
    try:
        subp.check_call('hmmscan ' + hmmtype + ' --cpu ' + str(cpucore) + ' --domtblout '
                        + outputfile + ' ' + database + ' ' + inputfile + ' >> protein_seq_c_ORFs.log', shell=True)
    except subp.CalledProcessError:
        print("hmmscan has some unexpected problems, please please check whether the parameters are accurate;\n \
        If you are sure that the parameters are correct, please check whether hmmscan is working properly in the analysis environment.")


# HMM result split
def hmm_parser(line):
    lis1 = line.rstrip().split(' ')
    lis2 = [x for x in lis1 if x != '']
    return lis2


##################################################################################
### extract sequence from large database
### inputfile: subdatabase hmm mapping result; ORF file(faa)
### ouputfile: hmm_output.replace('.tlout','')+'.extracted.fa'
### 
##################################################################################

def extract_seqs(hmm_result, faafile, extract_output):
    ###Enter the full name of .tlout file
    ###Enter the full name of the fasta contains the query sequences
    global n
    fileinput = open(hmm_result, 'r')  # full name of .tlout file
    fileoutput = open(extract_output, 'w')  # full name of the fasta contains the query sequences
    print("Extracting sequence")

    a = {}
    m, k = 0, 0
    for line in open(hmm_result, 'r'):
        if line.startswith('#'):
            continue
        lis = hmm_parser(line)
        m += 1
        try:
            percent = float(lis[2]) / float(lis[5])
            if percent < 0:  ## Filter condition ('if percent < 0' means No filtering)
                k += 1
                continue
            else:
                a[lis[3]] = m  ## store query ids as the key of a dic
        except IndexError:
            print(line)
    print(len(a), 'unique ids in ' + hmm_result)
    print(str(k) + 'No data filtering')

    Num = 0
    for line in open(faafile, 'r'):

        Num += 1
        if Num % 100000 == 0:
            print(Num, 'sequences have been searched!')

        if line.startswith('>'):
            ID = str(line.rstrip().split(' ')[0][1:])
            n = 0
            try:
                j = a[ID]
                fileoutput.write(line)
                n += 1
            except KeyError:
                continue
        else:
            if n == 1:
                fileoutput.write(line)
            else:
                continue

    fileinput.close()
    fileoutput.close()
    print('OK,', 'Extract Sequence Finished!')


##################################################################################
### hmm_compare
### inputfile: subdatabase hmm mapping result; large database hmm mapping result
### ouputfile: hmm_output.replace('.tlout','')+'.extracted.fa'
##################################################################################

def compare_result(arghmmoutfile, arglargedatahmmoutfile, outputfile):
    a, b = {}, {}
    ### hmm output from the whole database (PFAMTIGR) and sub-databases (Resfams)
    for line in open(arglargedatahmmoutfile, 'r'):
        if line[0] == '#':
            continue
        else:
            lis = hmm_parser(line)
            ID = lis[3]
            PFAMTIGR = [ID, lis[0], lis[1], lis[13]]
            try:
                if a[lis[3]] >= float(lis[13]):
                    continue
                else:
                    a[lis[3]] = float(lis[13])
                    b[lis[3]] = PFAMTIGR
            except KeyError:
                a[lis[3]] = float(lis[13])
                b[lis[3]] = PFAMTIGR

    print(str(len(b)) + ' unique hits')

    # outputfile
    f = open(outputfile + '.compare.csv', 'w')
    f1 = open(outputfile + '.check.csv', 'w')

    f.write(','.join(['Database', 'ORF_ID', 'Resfams_name', 'Resfams_ID',
                      'Resfams_bitscore', 'ORF_ID', 'PFAMTIGR_name', 'PFAMTIGR_ID', 'PFAMTIGR_bitscore']) + '\n')
    f1.write(','.join(['Resfams_name', 'Resfams_ID', 'PFAMTIGR_name', 'PFAMTIGR_ID', 'Frequency']) + '\n')
    ### hmm output from subdatabase (e.g., Resfams)
    c, d = {}, []
    for line in open(arghmmoutfile, 'r'):  # prodigal_result_ORFs_cut_ga.hmm.domain.tlout.txt
        if line[0] == '#':
            continue
        else:
            lis = hmm_parser(line)
            ID = lis[3]
            Resfams = [ID, lis[0], lis[1], lis[13]]
            try:
                if c[lis[3]] >= float(lis[13]):
                    continue
                else:
                    c[lis[3]] = float(lis[13])
                    try:
                        PFAMTIGR = b[lis[3]]
                        if float(Resfams[3]) >= float(PFAMTIGR[3]):
                            f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')
                        else:
                            if Resfams[1] == PFAMTIGR[1] or Resfams[2] == PFAMTIGR[2]:
                                f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')
                            else:
                                f.write(','.join(['PFAMTIGR'] + Resfams + PFAMTIGR) + '\n')
                                d.append('__'.join(
                                    Resfams[1:3] + PFAMTIGR[1:3]))  ### check whether it is the same HMM profiles
                    except KeyError:
                        print(lis[3], 'not found in the PFAMTIGR')
                        PFAMTIGR = [ID, '', '', '']
                        f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')

            except KeyError:
                c[lis[3]] = float(lis[13])
                try:
                    PFAMTIGR = b[lis[3]]
                    if float(Resfams[3]) >= float(PFAMTIGR[3]):
                        f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')
                    else:
                        if Resfams[1] == PFAMTIGR[1] or Resfams[2] == PFAMTIGR[2]:
                            f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')
                        else:
                            f.write(','.join(['PFAMTIGR'] + Resfams + PFAMTIGR) + '\n')
                            d.append('__'.join(
                                Resfams[1:3] + PFAMTIGR[1:3]))  ### check whether it is the same HMM profiles
                except KeyError:
                    print(lis[3] + ' not found in the PFAMTIGR')
                    PFAMTIGR = [ID, '', '', '']
                    f.write(','.join(['Resfams'] + Resfams + PFAMTIGR) + '\n')
    f.close()
    d1 = list(set(d))
    d1.sort()

    for item in d1:
        f1.write(','.join(item.split('__') + [str(d.count(item))]) + '\n')

    f1.close()
    print('Compare_result_DONE!')

def classifier(metadata, compareout, finalized):
    Resfams = {}
    for line in open(metadata, 'r'):  # Resfams_metadata_v0.1.txt
        lis = line.strip().split('\t')
        Resfams[lis[0]] = [lis[1], lis[2], lis[4], lis[8]]
    print(len(Resfams), 'records in the Resfams database')

    m, n = 0, 0
    f2 = open(finalized, 'w')  # MarianaTrench_ARGs.finalized.csv
    f2.write(','.join(
        ['Query_ID', 'Subject_ID', 'Database', 'ARG_type', 'ARG_subtype', 'Full annotation', 'Resistance_mechanisms',
         'Bitscore']) + '\n')

    for line in open(compareout, 'r'):  # compare output file
        lis = line.strip().split(',')
        if lis[0] == 'Resfams':
            try:
                annotation1 = Resfams[lis[3]]
                outlist = [lis[1], lis[3], lis[0]] + annotation1 + [lis[4]]
                f2.write(','.join(outlist) + '\n')
            except KeyError:
                continue
            m += 1
        else:
            try:
                annotation2 = Resfams[lis[7]]
                outlist = [lis[1], lis[7], lis[0]] + annotation2 + [lis[8]]
                f2.write(','.join(outlist) + '\n')
                m += 1
            except KeyError:
                n += 1
                continue
    print(m, 'genes were finalized as ARGs')
    print(n, 'genes were discarded')
    f2.close()

def develop():
    # get pars
    pars = read_params(sys.argv)

    metadata_arg = "/".join([pwdpath, 'ARGfams_metadata_v0.1.txt'])
    hmm_output = pars.output_file_name + '.hmm.tlout'
    hmm_large_output = pars.output_file_name + '_bdb.hmm.tlout'
    extract_arg = hmm_output.replace('.tlout', '') + '.extracted.fa'
    compare_output = pars.output_file_name  # outputfile + '_ARGfams.PFAMTIGRs.compare.csv'

    classifierinput = pars.output_file_name + '.compare.csv'
    finalized_arg = pars.output_file_name + '.finalized.csv'

    ### check database
    #   if expression:
    #       pass
    # DELETE

    hmm_subdatabase(pars.input_file, hmm_output, pars.cpucore, pars.sub_db, pars.hmmtype)
    extract_seqs(hmm_output, pars.input_file, extract_arg)
    hmm_largedatabase(extract_arg, hmm_large_output, pars.cpucore, pars.large_db, pars.hmmtype)
    compare_result(hmm_output, hmm_large_output, compare_output)
    classifier(metadata_arg, classifierinput, finalized_arg)
if __name__ == '__main__':
    develop()

# updata Database

# 不切片
# 保留compare
# 20种耐药基因，(TETA, card)
# 
