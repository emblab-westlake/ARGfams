#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 20:23:13 2022
@author: xinyu
"""

__author__ = ("Ju Feng",
              "Huang xinyu"
              "Zhang Guoqing")
__version__ = '0.5.0'
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
        help="Prints the current MGEfams tool version and exit")
    arg("-h", "--help", action="help", help="show this help message and exit")

    pars = p.parse_args()
    return pars

def hmm_subdatabase(inputfile, outputfile, cpucore, database, hmmtype):
    try:
        subp.check_call('hmmscan ' + hmmtype + ' --cpu ' + str(cpucore) + ' --domtblout '
                        + outputfile + ' ' + database + ' ' + inputfile + ' >> protein_seq_c_ORFs.log', shell=True)
    except subp.CalledProcessError:
        print("hmmscan has some unexpected problems, please please check whether the parameters are accurate;\n \
        If you are sure that the parameters are correct, please check whether hmmscan is working properly in the analysis environment.")

















def hmm_parser(line):
    lis1 = line.rstrip().split(' ')
    lis2 = [x for x in lis1 if x != '']
    return lis2

def Read_hmm_result(hmm_output):
    ###Enter the full name of .tlout file
    ###Enter the full name of the fasta contains the query sequences
    global n
    fileinput = open(hmm_output, 'r') 

    a,b = {},{}

    for line in open(hmm_output,'r'):
        if line[0]=='#':
            continue
        else:
            lis = hmm_parser(line)
            ID = lis[3] 
            PFAMTIGR = [ID, lis[0], lis[1],lis[13]]
            try:
                if a[lis[3]] >= float(lis[13]):
                    continue
                else:
                    a[lis[3]] = float(lis[13])
                    b[lis[3]] = PFAMTIGR
            except KeyError:
                 a[lis[3]] = float(lis[13])
                 b[lis[3]] = PFAMTIGR
                 
    print(str(len(b))+' unique hits')

    fileinput.close()

    c,d = {},[]
    f = open('process_file.csv','w')

    f.write(','.join(['Database','ORF_ID','Resfams_name','Resfams_ID','Resfams_bitscore','ORF_ID','PFAMTIGR_name','PFAMTIGR_ID','PFAMTIGR_bitscore'])+'\n')



### hmm output from subdatabase (e.g., Resfams)
#for line in open('CFUs_bins_renamed_merged.Resfams.hmm.domain.cov0.95.besthit.txt','r'):

    fileinput = open(hmm_output, 'r')

    for line in open(hmm_output,'r'):
        if line[0]=='#':
            continue
        else:
            lis = hmm_parser(line)
            ID = lis[3] 
            Resfams = [ID, lis[0], lis[1],lis[13]]
            try:
                if c[lis[3]] >= float(lis[13]):
                    continue
                else:
                    c[lis[3]] = float(lis[13])
                    try:
                        PFAMTIGR = b[lis[3]]
                        if float(Resfams[3]) >= float(PFAMTIGR[3]):
                            f.write(','.join(['Resfams'] + Resfams + PFAMTIGR)+'\n')
                        else:
                            if Resfams[1] == PFAMTIGR[1] or Resfams[2] == PFAMTIGR[2]:
                                f.write(','.join(['Resfams'] + Resfams + PFAMTIGR)+'\n')
                            else:
                                f.write(','.join(['PFAMTIGR']  +Resfams + PFAMTIGR)+'\n')
                                d.append('__'.join(Resfams[1:3]+PFAMTIGR[1:3]))  ### check whether it is the same HMM profiles  
                    except KeyError:
                        print(lis[3]), 'not found in the PFAMTIGR'
                        PFAMTIGR = [ID, '', '','']
                        f.write(','.join(['Resfams'] +Resfams + PFAMTIGR)+'\n')
                    
            except KeyError:
                c[lis[3]] = float(lis[13])
                try:
                    PFAMTIGR = b[lis[3]]
                    if float(Resfams[3]) >= float(PFAMTIGR[3]):
                        f.write(','.join(['Resfams'] + Resfams + PFAMTIGR)+'\n')
                    else:
                        if Resfams[1] == PFAMTIGR[1] or Resfams[2] == PFAMTIGR[2]:
                            f.write(','.join(['Resfams'] + Resfams + PFAMTIGR)+'\n')
                        else:
                            f.write(','.join(['PFAMTIGR']  +Resfams + PFAMTIGR)+'\n')
                            d.append('__'.join(Resfams[1:3]+PFAMTIGR[1:3]))  ### check whether it is the same HMM profiles      
                except KeyError:
                    print(lis[3]), 'not found in the PFAMTIGR'
                    PFAMTIGR = [ID, '', '','']
                    f.write(','.join(['Resfams'] +Resfams + PFAMTIGR)+'\n')
    f.close()
    d1 = list(set(d))
    d1.sort()

    fileinput.close()
       
#    print('DONE!')

    Resfams = {}
    for line in open('ARGfams_metadata_v0.1.txt','r'):
        lis = line.strip().split('\t')
        Resfams[lis[0]] = [lis[1],lis[2],lis[4],lis[8]]
#    print(len(Resfams),'records in the Resfams database')

    m, n = 0, 0
    f2=open('Finalized_Rawdata.txt','w')
    f2.write('\t'.join(['Query_ID','Subject_ID','Database','ARG_type','ARG_subtype','Full annotation', 'Resistance_mechanisms', 'Bitscore'])+'\n')

    for line in open('process_file.csv','r'):
        lis = line.strip().split(',')
        if lis[0] == 'Resfams':
            try:
                annotation1 = Resfams[lis[3]]
                outlist = [lis[1], lis[3], lis[0]] + annotation1 + [lis[4]]
                f2.write('\t'.join(outlist)+'\n')
            except KeyError:
                continue
            m+=1
        else:
            try:
                annotation2 = Resfams[lis[7]]
                outlist = [lis[1], lis[7], lis[0]] + annotation2 + [lis[8]]
                f2.write('\t'.join(outlist)+'\n')
                m+=1
            except KeyError:
                n+=1
                continue 

    f2.close()



    r7 = open('Finalized_bit50.txt','w')
    r7.write('\t'.join(['Query_ID','Subject_ID','Database','ARG_type','ARG_subtype','Full annotation', 'Resistance_mechanisms', 'Bitscore'])+'\n')

    f7 = open('Finalized_Rawdata.txt','r')
    next(f7)

    for line7 in f7:
        bitscore = float(line7.split('\t')[7])
        if bitscore >= 50:
            r7.write(line7)
    r7.close()
    f7.close()


    Types_and_sybtypes = []
    f3 = open('ARGfams_metadata_v0.1.txt')
    next(f3)

    for line3 in f3:
        Type = line3.split('\t')[1]
        SubType = line3.split('\t')[2]
        Type_and_sybtype = Type + '___' + SubType
        Types_and_sybtypes.append(Type_and_sybtype)

    f3.close()

    Types_and_sybtypes = list(set(Types_and_sybtypes))
#    print(len(Types_and_sybtypes))


    MAGs = []
    f4 = open('Finalized_bit50.txt', 'r')
    next(f4)
    for line4 in f4:
        Raw_MAG_info = line4.split('\t')[0]
        MAG_info = '_'.join(Raw_MAG_info.split('_')[:-1])  
        MAGs.append(MAG_info)

    f4.close()
    MAGs = list(set(MAGs))



    matrix = {}
    for MAG in MAGs:
        matrix[MAG] = {}
        for Types_and_sybtype in Types_and_sybtypes:
            matrix[MAG][Types_and_sybtype] = 0


    #Delect_Type = ['multidrug']
    f3 = open('Finalized_bit50.txt', 'r')
    next(f3)
    for line3 in f3:
        TypeInfo = line3.split('\t')[3]
        SubTypeInfo = line3.split('\t')[4]
        Type_and_subtype_Info = TypeInfo + '___' + SubTypeInfo
        Raw_MAG_info = line3.split('\t')[0]
        MAG_info = '_'.join(Raw_MAG_info.split('_')[:-1])  
        matrix[MAG_info][Type_and_subtype_Info] += 1
    f3.close()


    r2 = open('ARG_Summary.txt','w')
    r2.write('\t')

    for MAG in MAGs:
        r2.write('\t' + MAG)
    r2.write('\n')
    for Types_and_sybtype in Types_and_sybtypes:
        Type_Info = Types_and_sybtype.split('___')[0]
        SubType_Info = Types_and_sybtype.split('___')[1]
        r2.write(Type_Info + '\t' + SubType_Info + '\t')
        for MAG in MAGs:
            r2.write(str(matrix[MAG][Types_and_sybtype])+'\t')
        r2.write('\n')
    r2.close()




    r8 = open('ARG_Summary_without_Multdrug.txt','w')
    f8 = open('ARG_Summary.txt', 'r')
    for line8 in f8:
        if line8.startswith('multidrug'):
            continue
        else:
            r8.write(line8)

    f8.close()
    r8.close()



def develop():
    # get pars
    pars = read_params(sys.argv)
    hmm_output = pars.output_file_name + '.hmm.tlout'


    hmm_subdatabase(pars.input_file, hmm_output, pars.cpucore, pars.sub_db, pars.hmmtype)
    Read_hmm_result(hmm_output)

    

if __name__ == '__main__':
    develop()



