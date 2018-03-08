#!/usr/bin/env python3

"""
This software is intended to analyze sequencing data from NGS RNA-seq. It 
performs quality check, trimming, and alignment to reference geneome, followed 
by metrics collection and count of mapped reads.

It uses the following resources:

    environment requirements:

        python: release 3.0 or higher
        fastqc
        trim galore
        STAR
        samtools
        metrics
        IGVtools
        htseq-count (HTSeq-0.5.4p1)

    reference files:

        reference genome      --   used for mapping
        gtf file              --   used for gene annotations

"""
################################################################################
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   #
# POSSIBILITY OF SUCH DAMAGE.                                                  #
################################################################################

__author__     = "Awni Mousa <amousa@rockefeller.edu>"
__version__    = "$ Rev: 1.1 $"
__date__       = "$ Date: 2018-02-01 (Thu, 01 Feb 2018) $"
__maintainer__ = "Awni Mousa"

################################################################################
# IMPORT
# common modules

import os
import subprocess
import argparse
import sqlite3 as lite
import stat
import re
from difflib import SequenceMatcher
import logging
import logging.handlers
from tools import *
import multiprocessing
import sys
import time
import pwd

################################################################################
# Dictionary
job_list = {    # input           process
                'qCheck' :  1,  # perform quality check on reads
                  'trim' :  2,  # trim fastq files
                 'align' :  3,  # align to genome and create sam/bam files
               'metrics' :  4,  # sequencing metrics to analyze mapped reads location
            'igvConvert' :  5,  # create bai/tdf files for igv tools
            'countReads' :  6,  # count mapped reads to the genome
                 'clean' :  7   # perform cleanup
           }

################################################################################

def main():
    """
    A method to analyze fastaq files.
    Performs Quality check for reads, trimming, and alignment, followed by
    counting reads.
    For each step a log line is created and stored in a log file.
    """
    # time stamp
    t1 = time.time()
    
    # determine host
    import socket
    HOST = socket.gethostname()
    # set paths on server
    DATA = "/data5/" if HOST.split('.')[0] == "nhgalaxy2" else "/data1/"
    dataBase = "/data5/public/all_log.sqlite" if HOST.split('.')[0] == \
               "nhgalaxy2" else "/data2/raw/public/all_log.sqlite"

    # instantiate global debugger
    Glogger = logging.getLogger('logdebug')
    Ghdlr = logging.FileHandler(os.path.join(os.getcwd() , 'Global.log'))
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    Ghdlr.setFormatter(formatter)
    Glogger.addHandler(Ghdlr)
    Glogger.setLevel(logging.INFO)

##### Collect global parameters
    
    PATHtoGENOMES = os.path.join("/","mnt","heintz-bambi1",
                                 "Genomes_and_Annotations")

    STAR_genome_config = {
      'mm10_sj'        :'%sgenomes/STARgenomeDir_mm10_sj' %DATA,
      'mm10_GENCODE_sj':'%sgenomes/STARgenomeDir_mm10_GENCODE_sj' %DATA,
      'mm10'           :'%sgenomes/STARgenomeDir_mm10' %DATA,
      'hg38'           :'%sgenomes/STAR.UCSC_hg38_Dec2013_GRCh38_assembly' %DATA,
      'monkey'         :'%sgenomes/STARgenomeDir_ChlSab1.1_sj' %DATA,
      'rat'            :'%sgenomes/STARgenomeDir_rn6_sj' %DATA
        }

    igv_genome_config = {
        'mm10' :'%sgenomes/IGV_genomes/mm10.chrom.sizes' %DATA,
        'hg38' :'%sgenomes/IGV_genomes/hg38.chrom.sizes' %DATA,
       'monkey':'%sgenomes/IGV_genomes/monkey.chrom.sizes' %DATA,
        'rat'  :'%sgenomes/IGV_genomes/rn6.chrom.sizes' %DATA
       }

    gtf_genome_config = {
    'mm10'  :   os.path.join( PATHtoGENOMES,
                             'GTFs','mm10_knownGene_exon_updated.gtf'),
    'hg38'  :   os.path.join( PATHtoGENOMES,
                             'GTFs','hg38_ucsc_exon_updated.gtf'),
    'monkey':   os.path.join( PATHtoGENOMES,
                             'GTFs','Chlorocebus_sabaeus.ChlSab1.1.90.chr.gtf'),
    'rat'   :   os.path.join( PATHtoGENOMES,
                             'GTFs','Rattus_norvegicus.Rnor_6.0.90.chr.gtf')
    }
    
    picard_genome_config = {
    'mm10':  ' '.join([
                 os.path.join( PATHtoGENOMES,
                              'Picard_mm10_refFlat.xxu.txt'),
                 os.path.join( PATHtoGENOMES,
                              'Picard_mm10_rrna_intervalList.xxu.txt')
                     ]),
    'hg38':  ' '.join([
                 os.path.join( PATHtoGENOMES,'Annotations','human',
                              'hg38_v22_refflat.txt'),
                 os.path.join( PATHtoGENOMES,'Annotations',
                              'human','hg38_rrna_intervallist.txt')
                     ]),
  'monkey':  ' '.join([
                 os.path.join( PATHtoGENOMES,'Annotations',
                              'ChlSab1.1','ChlSab1.1_refFlat.txt'),
                 os.path.join( PATHtoGENOMES,'Annotations','ChlSab1.1',
                              'ribosomalRNAs_ChlSab1.1.90.interval_list')
                     ]),
    'rat' :  ' '.join([
                 os.path.join( PATHtoGENOMES,'Annotations','rat','rn6_refFlat.txt'),
                 os.path.join( PATHtoGENOMES,'Annotations',
                              'rat','ribosomalRNAs_rn6.interval_list')
                     ])
    }
        
    # set Global Parameters
    GlobalParameters = {
        'pDir'     : args.parentDirectory,
        'sDir'     : args.sourceDirectory,
        'align'    : STAR_genome_config[args.genome_config],
        'metrics'  : picard_genome_config[args.genome_config.split('_')[0]],
        'igv'      : igv_genome_config[args.genome_config.split('_')[0]],
        'gtf'      : gtf_genome_config[args.genome_config.split('_')[0]],
        'gtf_f'    : args.gtf_feature,
        'gtf_attr' : args.gtf_attribute
        }
    
    if args.star_genomePath is not None:
        GlobalParameters['align'] = args.star_genomePath

    if args.igv_genomePath is not None:
        GlobalParameters['igv'] = args.igv_genomePath

    if args.gtf is not None:
        GlobalParameters['gtf'] = args.gtf
        
    gp = {'pDir'    :'Parent Directory',
          'sDir'    :'Source Directory',
          'align'   :'STAR ref genome',
          'metrics' :'sequencing metrics',
          'igv'     :'IGV tools ref genome',
          'gtf'     :'Path to GTF',
          'gtf_f'   :'Feature type',
          'gtf_attr':'Attribute type'}

    gp_names = ['pDir','sDir','align','metrics','igv','gtf','gtf_f','gtf_attr']
    gp_str=[''.join([gp[i].rjust(20),' | ',GlobalParameters[i].ljust(20),'\n'])
            for i in gp_names]
    gp_str = ''.join(gp_str)

    Glogger.info('*' * 64 )
    Glogger.info('Program:\t' + str(__file__.split('/')[-1]) )
    Glogger.info('Author:\t' + str(__author__) )
    Glogger.info('Version:\t' + str(__version__) )
    Glogger.info('*' * 64 + '\n\n\t global parameters:\n\t %s\n'
                 %('-'*len('global parameters:')) + str(gp_str) )

##### collect fastq file names
    fQnames = []

    # for single end fQnames input
    if args.inputFile:
        fQnames += args.inputFile

    # for paired end fQnames input
    if args.Mates:
        fQnames += args.Mates

    # for batch fQnames stored in a file
    if args.lbatch:
        # read file's content and add fQnames to fQnames placeholder
        try:
            f = open( args.lbatch , 'r' )
            lines = [','.join(l.strip().split()) for l in f.readlines()]
            fQnames += lines
            f.close()
        except Exception as e:
            Glogger.error('batch file: ' + str(e) )

    # remove duplicated items and trailing suffixes
    fQnames = sorted(set(fQnames))
    # derive the directory names from the fQnames
    dNames = {}
    fqToRemove = []
    for fq in fQnames:
        ids = fq.split(",")
        if len(ids) == 1:
            keyDN = re.sub(".f[ast]*q[.tar]*.gz$","",ids[0])
            #keyDN = re.sub('.fastq$','',
             #           re.sub('.tar$','',
              #              re.sub('.gz$','',ids[0])))
            dNames[keyDN] = ids
        elif len(ids) == 2:
            s = SequenceMatcher(None, ids[0], ids[1])
            for i,j,n in s.get_matching_blocks():
                keyDN = re.sub('_R$','_paired',ids[0][i:n])
                break
            # make sure Mate1 comes first in order
            try:
                Mate1,Mate2,S = PE_check(ids,GlobalParameters['sDir'])
                dNames[keyDN] = [Mate1,Mate2]
            except Exception as e:
                Glogger.error('MATES: ' + str(e) )
        else:
            fqToRemove += [fq]
    if fqToRemove:
        Glogger.error('error fastq file names!\n'+'\n'.join(fqToRemove))
        [fQnames.remove(i) for i in fqToRemove]

##### collect processes
    processes = []
    # capture user input for running jobs if flagged
    if args.p:
        # empty processes
        processes = Inp('count reads') # capture and sort for user selected jobs
        try:
            assert(processes)
        except:
            Glogger.error('no jobs were selected!')
            sys.exit()

    # capture command line arguments if flagged
    elif args.comm:
        try:                        # need to manage typos quietly
            # convert character job_list to numeric based on dictionary
            processes = [job_list[i] for i in args.comm]
            assert( processes )
        except Exception as e:
            Glogger.error('unrecognized job name: ' + str(e) )
            print('ERROR: unrecognized job name: ' + str(e) )
            sys.exit()

    # remove duplicated items
    processes = sorted(set(processes))

    # quit if selected 'quit'
    if 8 in processes:
        Glogger.info('quit by user!')
        sys.exit('quit!')
    
    # when choosing the cleanup option only, cleanup will be carried out on
    # all directories. There is no need to provide fQnames
    if processes and not fQnames:
        if len(processes) == 1 and 7 in processes:
            run_job( (clean , (GlobalParameters['pDir'], Glogger)))
        else:
            Glogger.error('missing fastq file names!')
            sys.exit('missing fastq file names!')
            
    if fQnames and not processes:
        # default selection of jobs 1...7
        processes = list(range(1,8))

    # Intended for logging purposes:
    # will print to logger which jobs were selected to process
    j = [[key for key,value in job_list.items() if value == process][0] 
              for process in processes]
    Glogger.info('selected jobs | ' + ' '.join(j) )
    
    globalPar_filename='globalPar.%s.out' %multiprocessing.current_process().pid
    f = open(globalPar_filename , 'w')
    f.write('*' * 80)
    f.write('\nProgram:'.ljust(20) + str(__file__.split('/')[-1]) )
    f.write('\nAuthor:'.ljust(20) + str(__author__) )
    f.write('\nVersion:'.ljust(20) + str(__version__) )
    f.write('\nDate:'.ljust(20) + time.strftime("%c") )
    f.write('\n' + '*' * 80 + '\n\n\t global parameters:\n\t %s\n'
                 %('-'*len('global parameters:')) + str(gp_str) )
    f.write('\n' + 'selected jobs'.rjust(20) + ' | ' + ' '.join(j) )
    f.write('\n' + '*' * 80 + '\n')
    f.close()

    print('selected jobs: %s' %processes)

##### START PROCESSING EACH FASTQ FILE
    if dNames:
        for key,value in dNames.items():
            Glogger.info('inputFiles: %s' %'\t'.join(value))
    
            # add entry to summary table
            if not os.path.exists(dataBase):
                db( dataBase )  # create database table
            # connect to database table
            con = lite.connect( dataBase )
            with con:
                cur = con.cursor()
                t = time.asctime().strip().split()
                t = ' '.join([t[1],t[2],t[-1]])
                to_insert = "INSERT OR IGNORE INTO log_data ( ID, UsrID," +\
                            " notes ) VALUES (\'" + key +\
                            "\',\'" + pwd.getpwuid(os.getuid())[0] +\
                            "\',\'" + t + "\')"
                # run sqlite command
                cur.execute( to_insert )
            # close connection
            con.commit()
    
        # set the cap of maximum multiprocesses
        n = len(dNames) if len(dNames) < 11 else 10
        s = multiprocessing.Semaphore(n)
        # map processes
        
############        

        d = [ multiprocessing.Process(name = key , target = worker,
                                      args = (s, key, value, processes,
                                              GlobalParameters, Glogger,
                                              globalPar_filename, dataBase)
                                      )
             for key,value in dNames.items() ]
        # start all processes
        for task in d:
            task.start()
        # add processes to the pool
        for task in d:
            task.join()
    
    # time stamp
    t2 = time.time()

    # calculate total processing time
    Glogger.info('Total processing time: %s hrs\n' %round((t2-t1)/3600,2) +
                 '-' * 80 + '\n')
    os.remove(globalPar_filename)

################################################################################

# this will be executed when the script is run
if __name__ == '__main__':

    # parser object for managing input options
    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '''
        WORK FLOW TO PROCESS FASTq DATA FILES:
        -------------------------------------------------
            1) Check the reads for quality and trim contaminants
            2) Perform alignment to target genome and generate metrics
            3) Count mapped reads
        ''' ,
        epilog = '''
        All data created by this script will be stored in target directories. 
        By the end of run, each fastq file will be analyzed based on a variety 
        of input values and a number of cutoff thresholds. The end point of run 
        will be creating a counts of reads of RNAseq data mapped to the target 
        genome.
        ''' )

    parser.add_argument( '-i' , dest = 'inputFile' , nargs = '*' ,
        help = 'Space separated input fastq files. For paired-end fastq files,\
               use "-pe" parameter.')
    
    parser.add_argument( '-pe' , dest = 'Mates' , nargs = '*' ,
        help = 'Comma separated paired-end fastq files, in the form of: \
               mate1,mate2. Multiple PE fastq files can be provided separated \
               by white space: [mate1,mate2 [mate1,mate2 .... [mate1,mate2]]].')

    parser.add_argument( '-l', dest = 'lbatch' , nargs = '?' ,
        help = 'Filename with a list of fastq files to be processed. Paired-end \
               fastq files should be listed as mate1 mate2 in the same row, \
               separated by white space.' )

    parser.add_argument( '-pd', dest = 'parentDirectory' , nargs = '?' ,
            default = os.getcwd() ,
        help = 'Parent Directory path where the new processed fastq files will \
                be stored in child directories. Default: current working \
                directory' )

    parser.add_argument( '-sd', dest = 'sourceDirectory' , nargs = '?' ,
            default = '/mnt/heintz-bambi1/FastQ/' , const = os.getcwd() ,
        help = 'Source directory of the input fastq. Default: bambi1/FastQ.\
               If this option is encountered with no command-line argument \
               following it, it will default to the current working directory' )

    parser.add_argument( '--p', action = 'store_true' ,
        help = 'Select job from a list to processes. To be used with inputFile \
                argument.' )

    parser.add_argument( '-c', dest = 'comm' , nargs = '*' ,
        help = 'Command line argument to run jobs: qCheck trim align \
                metrics igvConvert countReads clean' )

    parser.add_argument( '--genome-config' , required = False ,
                         choices = ['mm10_sj','mm10_GENCODE_sj',
                                    'mm10','hg38','monkey','rat'] ,
                         default = 'mm10_GENCODE_sj',
        help = 'Mus musculus mouse reference genome configuration. Default: \
                mm10_GENCODE_sj; Sets the evironmental variables of reference \
                genomes and gtf to mm10 configuration. Can be switched\
                to any of the choice configurations.' )

    parser.add_argument( '--star-genomePath' , required = False ,
        help = 'Set path to STAR reference genome. Overrides reference genome \
                         configuration defaults.')

    parser.add_argument( '--igv-genomePath' , required = False ,
        help = 'Set path to igv reference genome. Overrides reference genome \
                         configuration defaults.')

    parser.add_argument( '--gtf' , required = False ,
        help = 'Set path to gtf. Overrides reference genome configuration \
                         defaults.')

    parser.add_argument( '--gtf-feature' , required = False , default = 'exon',
        help = 'Set the feature type (3rd column in GTF file) to be used. \
                         default: exon')

    parser.add_argument( '--gtf-attribute' , required = False ,
                         default = 'gene_id',
        help = 'Set the attribute type (9th column in GTF file) to be used. \
                         default: gene_id')
    
    # load the inputs
    args = parser.parse_args()

    # defaults to reading the system args
    main()

