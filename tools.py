import os
import stat
import subprocess
import glob
import re
import numpy as np

def db( dataBase='all_log.sqlite' ):
    '''
    a simple function to create sql database that holds progress information
    for easy tracking
    '''
    import sqlite3 as lite
    con = lite.connect( dataBase )
    with con:

        cur = con.cursor()
        # create log table
        cur.execute("CREATE TABLE IF NOT EXISTS log_data(ID TEXT " +\
                    "PRIMARY KEY, q_check TEXT, trim TEXT, align TEXT, " +\
                    "picard TEXT, igv TEXT, count TEXT, callPeaks TEXT, " +\
                    "clean TEXT, UsrID TEXT, notes TEXT)" )
    # close connection
    con.commit()
    subprocess.call('chmod 0777 %s' %dataBase, shell = True)


def all_log( i , log_proc , dataBase='all_log.sqlite', fail = '' ):
    '''
    Add log information to log database
    '''
    import sqlite3 as lite
    index = {  1 : ['q_check',  'q'],
               2 : ['trim',     't'],
               3 : ['align',    'a'],
               4 : ['picard',   'p'],
               5 : ['igv',      'i'],
               6 : ['count',    'c'],
               8 : ['callPeaks','m'],
               9 : ['callPeaks','m'],
               7 : ['clean',   'Ok']    }
    if fail:
        col = fail
    else:
        col = index[ log_proc ][1]
    # open connection to log file
    con = lite.connect( dataBase )
    with con:
        cur = con.cursor()
        # update table
        to_insert = "UPDATE log_data SET " + index[ log_proc ][0] + " = \'" + \
                    col + "\' where ID = \'" + i + "\'"
        # run sqlite command
        cur.execute( to_insert )
    # close connection
    con.commit()


def Inp( i ):
    '''
    A small function that captures user's input choice of selected jobs
    '''
    n = 5       # 5 attempts increment
    text = '\n> select job number/s from the list:\n\n' +\
                      ' 1. quality check\n' +\
                      ' 2. trim\n' +\
                      ' 3. align\n' +\
                      ' 4. Metrics Tools\n' +\
                      ' 5. igv_converter\n' +\
                      ' 6. %s\n' %i +\
                      ' 7. clean\n' +\
                      ' 8. quit\n\n-->   ' 
    while n > 0:
        try:
            Type = input(text)

            # manage errors quietly
            if not Type:
                raise ValueError
            t = list(map(int, re.findall(r'\d+',Type)))   # keep only integers

            if not t:
                print (64*'-','\nInvalid selection!')
                raise ValueError

            # for selection out of range
            notlisted = [ n for n in t if n not in range(1,9) ]
            if notlisted:
                print (64*'-','\n Selection out of range!')
                raise ValueError

            return t
            break

        except:
            n = n-1
            if n > 0 :
                print ( '\nTry Again:' )


def qcheck( i, path=os.getcwd(), sd=os.getcwd(), logger='' ):
    '''
    A simple method to perform quality check on fastq files using fastqc
    '''
    logger.info('start: quality check...')
    
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    if len(i) >= 3 or not i:
        logger.error('Missing/Too many input fastq files')
        s = False
    else:
        try:
            for inFile in i:
                assert( os.stat( os.path.join(sd,inFile) )[stat.ST_SIZE] )
        except:
            logger.error( 'impaired fastq files in %s!' %sd)
            s = False
        else:
            inputfiles = [os.path.join(sd,inFile) for inFile in i]
            t = '2' if len(inputfiles) == 2 else '1'
            cmd = ['fastqc','-q','-t',t,'-o',path] + inputfiles
            # run quality check
            try:
                print(' '.join(cmd))
                p = subprocess.Popen( cmd,stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          universal_newlines = True )
                outs, errs = p.communicate()
                if errs:
                    logger.error(str(errs))
                if outs:
                    logger.info(str(outs))
                s = True
            except Exception as e:
                logger.error( 'quality check: ' + str(e) )
                s = False
    # rollback working directory
    os.chdir( wDir )
    logger.info('end: quality check...')
    return [s,'']


def trim( i, QC=False, path=os.getcwd(), sd=os.getcwd(), logger='' ):
    '''
    A simple function to perform trimming using trim galore.
    As input it asks for a fastq file name
    '''
    logger.info('start: trimming...')

    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    
    if len(i) >= 3 or not i:
        logger.error('Missing/Too many input fastq files')
        s = False
    else:
        try:
            for inFile in i:
                assert( os.stat( os.path.join(sd,inFile) )[stat.ST_SIZE] )
        except:
            logger.error( 'missing or impaired fastq files in %s!' %sd)
            s = False
        else:
            inputfiles = [os.path.join(sd,inFile) for inFile in i]
            _paired,t = [['--paired'],'2'] if len(inputfiles) == 2 else [[],'1']
            if QC:
                cmdQC=['--fastqc','--fastqc_args',' '.join(['-q','-t',t,
                       '-o',os.getcwd()])]
            else:
                cmdQC = []
            cmd = ['trim_galore','--stringency','3','--gzip'] + _paired +\
                   cmdQC + inputfiles
            # run trim galore
            try:
                print(' '.join(cmd))
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines = True)
                outs, errs = p.communicate()
                f = open("trimming_report.txt","w")
                f.write(outs)
                f.write(errs)
                f.close()
                s = True
            except Exception as e:
                logger.error( 'trimming: ' + str(e) )
                s = False
        finally:
            # rollback working directory
            os.chdir( wDir )
            logger.info('end: trimming...')
    return [s,'']


def bowtie2( ref, path=os.getcwd(), logger='' ):
    '''
    A simple function to perform alignment using bowtie2.
    As input it asks a fastq file name/s and a reference genome.
    Refer to manual for further information.
    '''
    print('start: mapping...')
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    threads = 32
    newFiles = glob.glob('*.fq.gz')
    try:
        assert( len(newFiles) < 3 and newFiles )
    except:
        logger.error('Empty/Too many fastq files to align')
        s = False
    else:
        try:
            for inFile in newFiles:
                assert( os.stat(inFile)[stat.ST_SIZE] )
        except:
            logger.error( 'missing or impaired fastq files in %s!' %path)
            s = False
        else:
            # run alignment
            try:
                if len(newFiles) == 2:
                    Mate1,Mate2,S = PE_check(newFiles)
                    assert(S)
                    logger.info('PAIRED: paired end fastq files detected...')
                    MATES = ['-1',Mate1,'-2',Mate2]
                else:
                    MATES = ['-U ' , newFiles[0]]
            except Exception as e:
                logger.error("PAIRED: " + str(e))
                logger.error("PAIRED: paired end fastq files don't match...")
                s = False
                MATES = None
            else:
                logger.info( '\nPerforming alignment...please wait\n' )
                # call for command line
                cmd = [ 'bowtie2', '-x', ref, '-p', str(threads), '-X', '2000', 
                        '--no-mixed', '--no-discordant', '--mm', '-S',
                        'alignedReads.sam'] + MATES
                cmd_grep = ['grep','-v','chrM','alignedReads.sam']
                cmd_awk=['awk','BEGIN {FS=OFS="\t"}; {if ($1 ~ /^@/) print '+
                         '$0; else if (($9>-100) && ($9<100)) print $0}']
                cmd_sam1 = ['samtools','view','-@','32','-bh','-F4']
                cmd_sam2 = ['samtools','sort','-@','32','-o',
                            'Aligned.sortedByCoord.out.bam']
                cmd_fragSize = ['awk','BEGIN {FS=OFS="\t"}; '+
                                '{print $9}', 'alignedReads.sam']
                
                logger.info("bowtie2 command:\n" + ' '.join(cmd))
                try:
                    print(' '.join(cmd))
                    p1 = subprocess.Popen( cmd,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE,
                                           universal_newlines=True)
                    p1.wait()
                    p5 = subprocess.Popen( cmd_fragSize,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE,
                                           universal_newlines=True)
                    p2 = subprocess.Popen( cmd_grep,
                                           #stdin=p1.stdout,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                    p3 = subprocess.Popen( cmd_awk,
                                           stdin=p2.stdout,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                    p4 = subprocess.Popen( cmd_sam1, 
                                           stdin=p3.stdout, 
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                    subprocess.call( cmd_sam2, stdin=p4.stdout)
                    p1.stdout.close()
                    p2.stdout.close()
                    p3.stdout.close()
                    p4.stdout.close()
                    outs = p1.stderr.read()
                    fragSize = p5.stdout.read()
                    p1.wait()
                    p2.wait()
                    p3.wait()
                    p4.wait()
                    if outs:
                        if 'error' in str(outs).lower():
                            logger.error("bowtie2:" + str(outs))
                            raise Exception("bowtie2 error occured!")
                        else:
                            logger.info("bowtie2 summary: " + str(outs))
                    
                    # fragment length metrics
                    fragmentLen = [abs(int(i.strip()))
                                   for i in fragSize.strip().split("\n")]
                    fragmentLenArr = np.array(fragmentLen,dtype=int)
                    fragmentLenArrBins = np.bincount(fragmentLenArr)
                    ii = np.nonzero(fragmentLenArrBins)[0]
                    fragmentLenArr = np.vstack((ii,fragmentLenArrBins[ii])).T
                    freqArr = fragmentLenArr[:,1]/sum(fragmentLenArr[:,1])
                    freqArr = np.resize(freqArr,(freqArr.size,1))
                    fragmentLenArr = np.hstack((fragmentLenArr,freqArr))
                    f = open("fragmentLenArr.txt","w")
                    f.write('\t'.join(["fLength","freq","density"])+'\n')
                    for i in fragmentLenArr:
                        row = [str(int(float(i[0]))),
                               str(int(float(i[1]))),
                               str("{0:.6f}".format(float(i[2])))]
                        f.write('\t'.join(row)+'\n')
                    f.close()

                    # rename bam file and cleanup sam
                    BAM = os.path.basename(path)+'.sortedByCoord.under100bp.bam'
                    cmd_rmdup = ['samtools','rmdup',
                                 'Aligned.sortedByCoord.out.bam',BAM]
                    subprocess.call( cmd_rmdup )
                    os.remove('Aligned.sortedByCoord.out.bam')
                    os.remove('alignedReads.sam')
                    cmd = ['samtools','index',BAM]
                    subprocess.call( cmd )
                    for IN in newFiles:
                        os.remove(IN)
                    logger.info('Alignment finished successfully')
                    s = True
                except Exception as e:
                    logger.error( 'mapping: ' + str(e) )
                    s = False
    finally:
        # rollback working directory
        os.chdir( wDir )
        logger.info('end: mapping...')
    return [s,'']


def align( ref, path=os.getcwd(), logger='', CALLPEAK=False, ATAC=False ):
    '''
    A simple function to perform alignment using STAR 2.4.2a.
    As input it asks a fastq file name/s and a reference genome.
    Refer to manual for further information.
    '''
    if ATAC:
        s,e = bowtie2( ref, path, logger )
        return [s,'']
    else:
        logger.info('start: mapping...')

        # set working directory
        wDir = os.getcwd()
        os.chdir( path )
    
        threads = 16
        newFiles = glob.glob('*.fq.gz')
    
        try:
            assert( len(newFiles) < 3 and newFiles )
        except:
            logger.error('Empty/Too many fastq files to align')
            s = False
        else:
            try:
                for inFile in newFiles:
                    assert( os.stat(inFile)[stat.ST_SIZE] )
            except:
                logger.error( 'missing or impaired fastq files in %s!' %path)
                s = False
            else:
            # run alignment
            # set the command line for alignment
    
            #                   command line for STAR
            #                   ---------------------
            #                  --genomeDir|path to Genome Directory
            #                --readFilesIn|fastq file name/s
            #                 --runThreadN|number of threasds to run
            # --outFilterMismatchNoverLmax|int: alignment will be output only if
            #                             | its ratio of mismatches to mapped
            #                             | length is less than this value
            #          --outSAMstrandField|intronMotif : strand derived from the
            #                             | intron motif. Reads with inconsistent
            #                             | and/or non-canonical introns are
            #                             | filtered out
                try:
                    if len(newFiles) == 2:
                        Mate1,Mate2,S = PE_check(newFiles)
                        assert(S)
                        logger.info('PAIRED: paired end fastq files detected...')
                        MATES = [Mate1,Mate2]
                    else:
                        MATES = [newFiles[0]]
                except Exception as e:
                    logger.error("PAIRED: " + str(e))
                    logger.error("PAIRED: paired end fastq files don't match...")
                    s = False
                    MATES = None
                else:
                    print ( '\nPerforming alignment...please wait\n' )
                    logger.info('Alignment process started.')
                    # call for command line
                    if CALLPEAK:
                        cmd_args = [ '--outFilterScoreMinOverLread','0', 
                                     '--outFilterMatchNminOverLread','0',
                                     '--outFilterMatchNmin','40',
                                     '--alignIntronMax','1',
                                     '--alignEndsType', 'EndToEnd']
                    else:
                        cmd_args = [ '--outFilterMismatchNoverLmax','0.06',
                                     '--outSAMstrandField', 'intronMotif', 
                                     '--outFilterType','BySJout' ]
                    cmd = [ 'STAR','--genomeDir',ref,'--runThreadN',
                             str(threads),'--limitBAMsortRAM','20000000000',
                             '--genomeLoad','LoadAndRemove','--outSAMtype','BAM',
                             'SortedByCoordinate','--readFilesCommand','zcat',
                             '--readFilesIn' ] + MATES + cmd_args
                    logger.info("STAR command:\n" + ' '.join(cmd))
                    try:
                        print(' '.join(cmd))
                        p = subprocess.Popen( cmd, stdout=subprocess.PIPE,
                                                   stderr=subprocess.PIPE,
                                                   universal_newlines = True)
                        outs, errs = p.communicate()
                        if outs:
                            logger.info("STAR STDOUT: " + str(outs))
                        if errs:
                            logger.error("STAR STDERR: " + str(errs))
                            raise Exception("STAR error occured!")
                        # rename bam file
                        BAM = path.split("/")[-1] + '.sortedByCoord.bam'
                        os.renames('Aligned.sortedByCoord.out.bam',BAM)
                        cmd = ['samtools','index', BAM]
                        subprocess.call( cmd )
                        for IN in newFiles:
                            os.remove(IN)
                        logger.info('Alignment finished successfully')
                        subprocess.call(['rm','-rf','_STARtmp'])
                        s = True
                    except Exception as e:
                        logger.error( 'mapping: ' + str(e) )
                        s = False
        finally:
            # rollback working directory
            os.chdir( wDir )
            logger.info('end: mapping...')
    return [s,'']


def Metrics( path=os.getcwd(),logger='',REF_FLAT='',
            RIBOSOMAL_INTERVALS='',ATAC=False):
    '''
    A simple function to location of mapped reads
    '''
    logger.info('start: Metrics tools...')
    
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    
    sorted_bams = glob.glob( '*.bam' )
    try:
        assert(sorted_bams)
        for sorted_bam in sorted_bams:
            assert( os.stat( sorted_bam )[stat.ST_SIZE] )
    except:
        logger.error( 'missing or impaired bam files!' )
        s = 'failed'
    else:
        for sorted_bam in sorted_bams:
            if ATAC:
                cmd = ['Rscript','/home/amousa/fragSizeDist.AM.R',sorted_bam]
            else:
                CHART_OUTPUT = sorted_bam + '.pdf'
                METRICS_OUTPUT = sorted_bam + '_metrics.txt'
                RIB_INT = ['RIBOSOMAL_INTERVALS=%s' %RIBOSOMAL_INTERVALS] \
                        if RIBOSOMAL_INTERVALS else []
                cmd = ['java','-jar','/usr/picard/CollectRnaSeqMetrics.jar',
                       'REF_FLAT=%s' %REF_FLAT,'STRAND_SPECIFICITY=NONE',
                       'CHART_OUTPUT=%s' %CHART_OUTPUT,'INPUT=%s' %sorted_bam,
                       'OUTPUT=%s' %METRICS_OUTPUT] + RIB_INT
            try:
                print(' '.join(cmd))
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines = True)
                outs, errs = p.communicate()
                if outs:
                    print('STDOUT: CollectRnaSeqMetrics\n' + str(outs))
                if errs:
                    print('STDERR: CollectRnaSeqMetrics\n' + str(errs))
                s = ''
            except Exception as e:
                logger.error( 'location metrics: ' + str(e) )
                s = 'failed'
    finally:
        logger.info('end: Metrics tools...')
        os.chdir( wDir )
    return [True,s]


def igv_tools( path=os.getcwd(), genome='', logger='', w='25'):
    '''
    A simple function to create indexed bam and tdf files
    '''
    logger.info('start: igv_converter...')
    
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    
    sorted_bams = glob.glob( '*.bam' )
    try:
        assert(sorted_bams)
        for sorted_bam in sorted_bams:
            assert( os.stat( sorted_bam )[stat.ST_SIZE] )
    except:
        logger.error( 'missing or impaired bam files!' )
        s = 'failed'
    else:
        for sorted_bam in sorted_bams:
            TDF = sorted_bam + '.tdf'
            # cmd = ['samtools','view','-@','32','-c','-f', '1',sorted_bam]
            # p = subprocess.Popen( cmd,
                                  # stdout=subprocess.PIPE,
                                  # stderr=subprocess.PIPE,
                                  # universal_newlines=True)
            # outs,errs = p.communicate()
            # if int(outs.strip()) > 0:
                # cmd=['igvtools','count','--pairs','-w',w,sorted_bam,TDF,genome]
            # else:
                # cmd=['igvtools','count','-w',w,sorted_bam,TDF,genome]
            cmd=['igvtools','count','-w',w,sorted_bam,TDF,genome]    
            try:
                print(' '.join(cmd))
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines = True)
                outs, errs = p.communicate()
                if errs:
                    logger.error('STDERR:\n' + str(errs))
                    logger.error('check for errors!')
                    s = 'failed'
                else:
                    if outs:
                        logger.info('STDOUT:\n' + str(outs))
                    s = ''
            except Exception as e:
                s = 'failed'
                break
        os.remove("igv.log")
    finally:
        logger.info('end: igv_converter...')
        os.chdir( wDir )
    return [True,s]


def callpeak( path=os.getcwd(), genome='', logger='', nomodel='', q='0.01'):
    '''
    A simple function to create indexed bam and tdf files
    '''
    if logger:
        logger.info('start: macs callpeak...')
    
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )
    
    sorted_bams = glob.glob( '*.bam' )
    try:
        assert(sorted_bams)
        for sorted_bam in sorted_bams:
            assert( os.stat( sorted_bam )[stat.ST_SIZE] )
    except:
        if logger:
            logger.error( 'missing or impaired bam files!' )
        else:
            print('missing or impaired bam files!'),
        s = False
    else:
        for sorted_bam in sorted_bams:
            ChIPname = sorted_bam
            cmd = ['samtools','view','-@','32','-c','-f', '1',sorted_bam]
            p = subprocess.Popen( cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  universal_newlines=True)
            outs,errs = p.communicate()
            BAM="BAMPE" if int(outs.strip()) > 0 else "BAM"
            cmd = ['macs2','callpeak','-t',sorted_bam,'-f',BAM,'-B',
                   '-g',genome,'-n',ChIPname,'-q',q]
            if nomodel:
                cmd += [nomodel]
                
            try:
                print(' '.join(cmd))
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines = True)
                outs, errs = p.communicate()
                if errs:
                    if logger:
                        logger.error('STDERR:\n' + str(errs))
                    else:
                        print('STDERR:\n' + str(errs)),
                if outs:
                    if logger:
                        logger.info('STDOUT:\n' + str(outs))
                    else:
                        print('STDOUT:\n' + str(outs)),
                s = True
            except Exception as e:
                if logger:
                    logger.error( 'macs2 callpeak: ' + str(e) )
                else:
                    print('macs2 callpeak: ' + str(e)),
                s = False
                break
    finally:
        if logger:
            logger.info('end: macs2 callpeak...')
        else:
            print('end: macs2 callpeak...'),
        os.chdir( wDir )
    return [s,'']


def count(path=os.getcwd(),i='',gtf='',logger='',gtf_f='exon',gtf_attr='gene_id'):
    '''
    call script to run count of reads; wait until job ends
    '''
    logger.info('start: count reads...')
    
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )    
    
    # locate BAMs in current directory
    BAMs = glob.glob('*.bam')
    if not BAMs:
        # if same data were processed before and the user interested in counting
        # reads using different options
        BAMs = glob.glob(os.path.join('/mnt/heintz-bambi2/BAM',i,'*.bam'))
    try:
        assert(BAMs)
    except:
        logger.error('missing alignment BAM files')
        s = False
    else:
        try:
            for BAM in BAMs:
                suff = '.' + gtf_f if gtf_f != 'exon' else ''
                outfile = BAM + suff + '.counts.txt'
                cmd = ['htseq-count','-q','-s','no','-r','pos','-i',gtf_attr,
                       '-t',gtf_f,'-f','bam',BAM,gtf]
                logger.info("HTSeq count command:\n"+' '.join(cmd))
                print(' '.join(cmd))
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines = True)
                outs, errs = p.communicate()
                if errs:
                    logger.error(str(errs))
                    if errs.startswith("Error"):
                        raise Exception
                with open(outfile,'w') as f:
                    f.write(outs)
                f.close()
                logger.info("htseq-count successfully completed!")
                s = True
        except Exception as e:
            logger.error('htseq-count:\n'+str(errs))
            s = False
    finally:
        logger.info('end: count reads...')
        os.chdir( wDir )
    return [s,'']


def clean( path , Glogger , i = '' ):
    '''
    A simple function to remove unnecessary files from the working environment
    For safest operation of this function, it requires that the user provides
    an identifying name of directory (fastq file name) to clean unnecessary
    files
    If the user does not provide a name, this method will retrieve the directory
    contents downwards the parent directory
    '''
    Glogger.info('start: cleanup...')
    print('moving files...')
        
    # set working directory
    wDir = os.getcwd()
    os.chdir( path )    
    
    if i:
        dirs = [os.path.join(path,i)]

    if not i:
        # optional for batch cleanup starting from top parent directory
        for root, d, files in os.walk( path ):
            if root == path:
                dirs = d
        # timed out prompt
        import signal
        signal.alarm(30)
        prompt = input('\n\n\t>>> Are you sure you want to cleanup all '+
            'directories (you have 30 seconds to respond before termination)?\n' +
            'a (abort), c (continue) >> ')
        if prompt == 'a':
            s = False
            sys.exit('Cleanup aborted!!')
        if prompt != 'c':
            s = False
            sys.exit('Unrecognized choice!!')

        # disable the alarm after success
        signal.alarm(0)
        
    for directory in dirs:
        cmd = 'rsync --remove-source-files --chmod=ug+rwx,o-rwx -r ' +\
               directory + ' /mnt/heintz-bambi2/BAM/'
        Glogger.info("cleanup:\n"+cmd)
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             universal_newlines = True , shell = True )
        outs, errs = p.communicate()
        try:
            assert(not errs)
            os.rmdir( directory )
            Glogger.info( 'rsync outs:\n' + str(outs) )
            s = True
        except:
            Glogger.error('rsync STDERR:\n' + str(errs))
            s = False
            
    os.chdir( wDir )
    return [s,'']
    

def PE_check(Mates, sd=''):
    try:
        import gzip
        f = gzip.open(os.path.join(sd,Mates[0]),'rt')
        line1 = f.readline().strip().split()
        f.close()
        f = gzip.open(os.path.join(sd,Mates[1]),'rt')
        line2 = f.readline().strip().split()
        f.close()
    except Exception as e:
        s = False
    else:
        try:
            # if data came from SRR public source, paired-end fastq files
            # are encoded differently
            if '@SRR' in line1[0]:
                # determine if paired
                pe1 = line1[0].split('.')[-1]
                pe2 = line2[0].split('.')[-1]
                assert( pe1 != pe2)
                if pe1 == '1' and pe2 == '2':
                    Mate1 = Mates[0]
                    Mate2 = Mates[1]
                elif pe1 == '2' and pe2 == '1':
                    Mate2 = Mates[0]
                    Mate1 = Mates[1]
                else:
                    raise ValueError
                s = True
            else:
                # make sure they have identical identifiers
                assert( line1[0] == line2[0] )
                # determine if paired
                pe1 = line1[1].split(':')[0]
                pe2 = line2[1].split(':')[0]
                assert( pe1 != pe2)
                if pe1 == '1' and pe2 == '2':
                    Mate1 = Mates[0]
                    Mate2 = Mates[1]
                elif pe1 == '2' and pe2 == '1':
                    Mate2 = Mates[0]
                    Mate1 = Mates[1]
                else:
                    raise ValueError
                s = True
        except Exception:
            Mate1 = ''
            Mate2 = ''
            s = False
    return [Mate1,Mate2,s]


def run_job( args ):
    '''
    args[0] - function (job) to run
    args[1] - input arguments to function
    result  - boolean: True or False
    '''
    result = args[0]( *args[1] )
    return result


def worker( s='',key='',values='',processes='',GlobalParameters='',  
            Glogger='', globalPar='', dataBase='', CALLPEAK=False, ATAC=False):
    '''
    A simple function to run jobs on fastq files.
        Inputs:
    s          - number of workers in pool waiting to be consumed by worker
    threads    - number of threads allowed in a single process
    key        - name of directory
    values     - fastq file names to be processed
    processes  - jobs to be run on a fastq file
    GlobalParameters - a dictionary of global parameters:
                       parent dir, source dir, path to STAR ref genome,
                       path to igv ref genome, path to gtf file
    Glogger    - global logger to Global.log
    globalPar  - global environment parametrs placeholder

    This function will be called by a multiprocessing object;
    each 'current_process' is identified by 'name' and 'ID' that can be tracked
    in terminal.
    '''
    with s:
        import multiprocessing
        p = multiprocessing.current_process()
        if key:
            print('--> processing:' , p.name, p.pid)
            # p.name and p.pid will be logged to allow easy tracking of processes in
            # terminal ($ ps aux | grep {p.name} OR $ ps p.pid)
            Glogger.info('process name: ' + p.name + ' process ID: ' + str(p.pid))
            # create a path to make a special directory to store output files
            # jobs will be stored recursively, in directories identified by
            # the fQnames.
            pDirectory = os.path.join(GlobalParameters['pDir'],key)
    
            # create special directory
            try:
                os.makedirs( pDirectory , exist_ok = True )
                os.chmod( pDirectory , mode=0o775 )
            except Exception as e:
                print(str(e))
                m = re.search('(?<=expected mode )\d*', str(e))
                MODE = int(m.group(0), 8) # convert octal to decimal
                os.chmod( pDirectory , mode=MODE )
            finally:
                # write global parameters to file
                import shutil
                shutil.copy2(globalPar, pDirectory)
                
                # instantiate individual debugger
                import logging
                import logging.handlers
                logger = logging.getLogger('Alogdebug.' + key)
                fdebug = os.path.join(pDirectory, 'debug.log')
                hdlr = logging.FileHandler( fdebug )
                formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
                hdlr.setFormatter(formatter)
                logger.addHandler(hdlr)
                logger.setLevel(logging.INFO)
        
                logger.info( 'Start processing ' + key )

                try:
                    cmd=['/usr/bin/sqlite3','-header','-column',
                         '/data4/amousa/processing/software.sqlite',
                         'select software,version,release_date,link from '+
                         'software_version where status==\'current\'']
                    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             universal_newlines = True)
                    outs, errs = p.communicate()
                    with open(os.path.join(pDirectory,"software_summary.txt"),
                              'w') as f:
                        f.write(outs)
                        f.close()
                except Exception as e:
                    logger.error("software_summary error!\n" + str(e))
                    logger.error("errs output:\n" + str(errs))

                # run job/s stored in processes
                tasks = {
#    function    parameters to function    #
#    --------    ----------------------    #
  2:(qcheck     , (values,pDirectory,GlobalParameters['sDir'],logger)),
  1:(trim       , (values,False,pDirectory,GlobalParameters['sDir'],logger)),
  3:(align      , (GlobalParameters['align'],pDirectory,logger,CALLPEAK,ATAC)),
  5:(igv_tools  , (pDirectory,GlobalParameters['igv'],logger)),
  7:(clean      , (GlobalParameters['pDir'],Glogger,key)),
                }
                if 6 in processes:
                    tasks[6] = (count, (pDirectory,key,
                                        GlobalParameters['gtf'],logger,
                                        GlobalParameters['gtf_f'],
                                        GlobalParameters['gtf_attr']))
                if 8 in processes:
                    tasks[8] = (callpeak, (pDirectory,
                                           GlobalParameters['callpeak'],logger))
                if 9 in processes:
                    tasks[9] = (callpeak, (pDirectory,
                                           GlobalParameters['callpeak'],logger,
                                           '--nomodel','0.05'))
                if ATAC:
                    tasks[4] = (Metrics, (pDirectory,logger,"","",True))
                    tasks[5] = (igv_tools, (pDirectory,
                                           GlobalParameters['igv'],logger,'5'))
                else:
                    tasks[4] = (Metrics,(pDirectory,logger,
                                         GlobalParameters['metrics'].split()[0],
                                         GlobalParameters['metrics'].split()[1]))
                if all(i in processes for i in [1,2]):
                    processes.remove(1)
                    tasks[2] = (trim , (values,True,pDirectory,
                                        GlobalParameters['sDir'],logger))
                processes = iter(processes)
                job = next(processes)
                while job:
                    success = run_job(tasks[job])
                    # break loop if failed in job
                    try:
                        assert( success[0] )
                        # update database
                        all_log( key , job , dataBase, success[1] )
                        try:
                            job = next(processes)
                        except:
                            break
                    except:
                        Glogger.error( 
                            "AssertionError: %s\nrun_job(tasks[%s] returned %s)" 
                            %(p.name,job,success) )
                        # update database
                        all_log( key , job , dataBase , success[1] )
                        job = None
                        break
        else:
            Glogger.error("ERROR WORKER missing key!")
            Glogger.error('ERROR WORKER: process name '+p.name+' process ID:'+str(p.pid))           

