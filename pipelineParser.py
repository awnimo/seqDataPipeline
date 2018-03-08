#!/usr/bin/env python3

################################################################################

__author__  = "Awni Mousa <amousa@rockefeller.edu>"
__version__ = "$ Rev: 0.2 $"
__date__    = "$ Date: 2018-02-01 (Thu, 01 Feb 2018) $"
__maintainer__ = "Awni Mousa"

################################################################################
import os, sys
import subprocess
import multiprocessing
import argparse
import logging
import datetime
from parserTools import *
import glob
import stat
################################################################################

def main():
    
##### check environment
    os.chdir("/mnt/heintz-bambi1/processing/RNASEQ_jobs")
    wDir = os.getcwd()
    today = datetime.date.today().strftime('%d-%b-%Y')
    # instantiate a logger
    logger = logging.getLogger('logdebug')
    hdlr = logging.FileHandler('temp/%s.log' %today)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    try:
        assert args.inputArgs
    except:
        sys.exit( "missing input!" )
        
    filenames = args.inputArgs
    # check if 'form' exists; if not omit from the list and report missing form
    i = ''
    toRemove = []
    # ENV
    for i in filenames:
        try:
            assert( os.path.exists(i) )
        except Exception as e:
            logger.error("ENV:\n"+str(e))
            logger.info(i)
            # remove flagged forms from list of processing
            f = open('temp/%s.notProcessed.txt' %today,'a+')
            f.write(i + "\n")
            f.close()
            filenames.remove( i )
    
##### parse the forms
    IDS = {}
    CONTACT = {}
    DICT = {}
    # Keys to DICT one of the following:
    #   RNA_mouse
    #   RNA_human
    #   RNA_mouse_ex
    #   RNA_human_ex
    #   RNA_mouse_nuclear
    #   RNA_human_nuclear
    #   RNA_mouse_nuclear_ex
    #   RNA_human_nuclear_ex
    #   CHIP_CLIP_mouse
    #   CHIP_CLIP_human
    #   CHIP_CLIP_mouse_ex
    #   CHIP_CLIP_human_ex
    #   GDNA_mouse
    #   GDNA_human
    for filename in filenames:
        filenameA = '@'+filename
        formInputs = Read_form(filename=filename)
        arguments = formInputs.parseForm()
        print(filename)
        print(arguments)
        changeName = ['mv',filename,filenameA]
        subprocess.call( changeName,stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,universal_newlines = True )
        # aggregate data arguments for batch processing
        if not arguments['userName'] in CONTACT:
            CONTACT[arguments['userName']]=' '.join(
                                                 reversed(arguments['contact']))
        # if application not specified assign RNA
        if arguments['application'].lower() not in \
                             ["rna-seq","gdna-seq","chip/clip-seq","atac-seq"]:
            logger.error("ERROR: application not specified!!")
        else:
            KEY = arguments['application'].split("-")[0].upper()
            KEY = '_'.join([re.sub('/','_',KEY),arguments['species'].lower()])
            if arguments['sampleType'].lower() != 'Nuclear'.lower():
                arguments['sampleType'] = ''
            KEY = '_'.join([KEY,arguments['sampleType'].lower()]).strip('_')
            if arguments['external'].lower() == 'y':
                KEY = KEY+'_ex'
            # assign key to dictionary
            if KEY not in DICT:
                DICT[KEY] = {'-i':[] , '-pe':[] , 'PARAM': []}
            if arguments['external'].lower() == 'y':
                if arguments['external_type'].upper() == "GEO":
                    DICT[KEY]['PARAM'] += ['-sd /mnt/heintz-bambi2/other_fastq/GEO']
                    DICT[KEY]['PARAM'] = list(set(DICT[KEY]['PARAM']))
                    # call for function to download, extract, and rename the fastq
                    fq = download_GEO(arguments['fq'],arguments['seqSampleID'],wDir,
                                      "/mnt/heintz-bambi2/other_fastq/GEO",logger)
                    arguments['fq'] = fq[0]
                    # safely set arguments['mates'] parameter
                    arguments['mates'] = fq[1]
                else:
                    DICT[KEY]['PARAM'] += ['-sd /mnt/heintz-bambi2/external_fastq']
                    DICT[KEY]['PARAM'] = list(set(DICT[KEY]['PARAM']))
                    fq = arguments['fq'].split(',')
                    # rename the fastq
                    for i in range(len(fq)):
                        try:
                            nfq = arguments['seqSampleID']+"_"+fq[i]
                            os.rename('/mnt/heintz-bambi2/external_fastq/' + fq[i],
                                      '/mnt/heintz-bambi2/external_fastq/' + nfq)
                            fq[i] = nfq
                        except Exception as e:
                            logger.error('External Fastq File Rename: ' + str(e))
                    arguments['fq'] = ','.join(fq)
                    # safely set arguments['mates'] parameter
                    if len(fq) == 1:
                        arguments['mates'] = 'SE'
                    else:
                        arguments['mates'] = 'PE'
            if arguments['mates'] == 'PE':
                DICT[KEY]['-pe'] += [arguments['fq']]
            elif arguments['mates'] == 'SE':
                DICT[KEY]['-i'] += [arguments['fq']]
            if not arguments['userName'] in IDS:
                IDS[arguments['userName']] = [arguments['seqSampleID']]
            else:
                IDS[arguments['userName']] += [arguments['seqSampleID']]
            # assign other arguments
            for k in DICT:
                # set genome configuration to human or rat
                if 'human' in k:
                    DICT[k]['PARAM'] += ['--genome-config hg38']
                if 'rat' in k:
                    DICT[k]['PARAM'] += ['--genome-config rat']
                if 'monkey' in k:
                    DICT[k]['PARAM'] += ['--genome-config monkey']
                # set gtf configuration for nuclear data
                PATHtoGTFs = os.path.join("/","mnt","heintz-bambi1",
                                          "Genomes_and_Annotations","GTFs")
                human_nuclear = os.path.join( PATHtoGTFs,'hg38_ucsc_wg_updated.gtf')
                mouse_nuclear = os.path.join( PATHtoGTFs,
                                             'mm10_knownGene_wg_updated.gtf')
                if 'RNA_human_nuclear' in k:
                    DICT[k]['PARAM'] += [' '.join(['--gtf',human_nuclear,
                                                   '--gtf-feature','gene'])]
                if 'RNA_mouse_nuclear' in k:
                    DICT[k]['PARAM'] += [' '.join(['--gtf',mouse_nuclear,
                                                   '--gtf-feature','gene'])]
                if 'RNA_rat_nuclear' in k:
                    DICT[k]['PARAM'] += ['--gtf-feature gene']
                if 'RNA_monkey_nuclear' in k:
                    DICT[k]['PARAM'] += ['--gtf-feature gene']
                DICT[k]['PARAM'] = list(set(DICT[k]['PARAM']))
##### prepare the command line arguments
    CMDs = []
    print("DICT:",DICT)
    for i in DICT:
        if DICT[i]['-i'] or DICT[i]['-pe']:
            if i.split("_")[0] == 'GDNA':
                PROGRAM = ''
                print("Application: ", i)
                print("Module for gDNA data not supported yet!")
                # report IDs
                f = open('temp/%s.notProcessed.txt' %today,'a+')
                r = [i+"\t"+n.split("_")[0] for n in DICT[i]['-i']+DICT[i]['-pe']]
                f.write('\n'.join(r) + "\n")
                f.close()
            elif i.split("_")[0] == 'CHIP':
                PROGRAM = ['mainChIP']
            elif i.split("_")[0] == 'ATAC':
                PROGRAM = ['mainATAC']
            else:
                 PROGRAM = ['mainScript']
        else:
            PROGRAM = ''
        if PROGRAM:
            PAR = [' '.join(DICT[i]['PARAM'])] if DICT[i]['PARAM'] else []
            FQ = ['-i'] + DICT[i]['-i'] if DICT[i]['-i'] else []
            PE = ['-pe'] + DICT[i]['-pe'] if DICT[i]['-pe'] else []
            if FQ:
                CMDs += [((PROGRAM , FQ , PAR),(i,len(FQ)-1))]
            if PE:
                CMDs += [((PROGRAM , PE , PAR),(i,len(PE)-1))]
##### multiprocess
    print("CMDs:",CMDs)
    if CMDs:
        n = len(CMDs) if len(CMDs) < 7 else 6
        s = multiprocessing.Semaphore(n)
        d = [ multiprocessing.Process(
                name = CMD[1][0] , target = runJobs , 
                args = (s,CMD[0],CMD[1][1],wDir,IDS,CONTACT,logger))
              for CMD in CMDs ]
        # start all processes
        for task in d:
            task.start()
        # add processes to the pool
        for task in d:
            task.join()
            

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument( 'inputArgs', nargs = '*' )
    args = parser.parse_args()
    
    main()

