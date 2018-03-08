#!/usr/bin/env python3
import re
import os
import glob
import subprocess
import stat
import itertools

class Read_form():
    '''
    This method will parse submitted form
    '''
    def __init__(self, filename=None, string=None):
        if bool(filename):
            self.filename = filename
            if os.path.isfile(self.filename):
                self.parseForm(filename=self.filename)
        else:
            self.string = string
            assert(self.string)
            self.parseForm(string=self.string)
            
    def __str__(self):
        return print(self.string)
        
    def parseForm(self,filename=None,string=None):
        if bool(filename):
            f = open(filename,'r')
            self.string = f.read()
            f.close()
        elif bool(string):
            self.string = string
        dataString = self.extractLines(self.string)
        return dataString
        
    def extractLines(self,string=None):
        '''
        This method will extract lines and parse the parameters from string 
        input of submitted form for DE
        '''
        if not string:
            string = self.string
        # parse ControlGenes
        FastQFiles = re.findall('<FastQFile>(.*?)</FastQFile>',string, re.DOTALL)
        # parse the remaining parameters        
        rows = string.strip().split('\n')
        # unpack rows
        contact,userName,application,species,sampleType,mates,external,\
        external_type,seqSampleID = self.unpack(*rows)
        FastQFile = FastQFiles[0].strip().split()
        FastQFile = [i.strip().split("/")[-1] for i in FastQFile]
        fq = ','.join(i for i in FastQFile)
        data = {
            'contact':contact,
            'userName':userName,
            'application':application,
            'species':species,
            'sampleType':sampleType,
            'mates':mates,
            'external':external,
            'external_type':external_type,
            'seqSampleID':seqSampleID,
            'fq':fq
        }
        return data
        
    def unpack(self, *rows):
        contact,application,species,sampleType,mates,external,external_type,\
        seqSampleID, *fq = rows 
        contact = contact.strip().split()[1:]
        userName = contact[-1].split('@')[0]
        contact = contact[:2]
        application = '_'.join(application.strip().split()[1:])
        species = species.strip().split()[-1]
        sampleType = '_'.join(sampleType.strip().split()[1:])
        mates = mates.strip().split()[-1]
        external = external.strip().split()[-1]
        external_type = external_type.strip().split()[-1]
        seqSampleID = seqSampleID.strip("_").split()[-1]
        rows = [
            contact,
            userName,
            application,
            species,
            sampleType,
            mates,
            external,
            external_type,
            seqSampleID
        ]
        return rows


def constructEmail(completed='',not_completed='',userName='', contact=''):
    '''
    Automatic email notification constructor
    '''
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    
    COMPLETED = ' '.join(completed)
    NON_COMPLETED = ' '.join(not_completed)
    # send email notification with job name and link when done
    bioinformaticsRU = 'bioinformaticsru@rockefeller.edu'
    # Create message container - the correct MIME type is multipart/alternative.
    msg = MIMEMultipart('alternative')
    msg['Subject'] = "Sequencing Data Notification Center"
    msg['From'] = bioinformaticsRU
    userEmailArddress = [userName+'@rockefeller.edu','amousa@rockefeller.edu']
    msg['To'] = ', '.join(userEmailArddress)
    # Create the body of the message (a plain-text and an HTML version).
    text1 = "Dear %s,\n\n" %contact
    text2 = "The following sequencing data were successfully processed and "\
            "completed:\n%s\n" %COMPLETED
    text3 = "The following sequencing data were NOT processed:\n%s\nPlease "\
            "contact Awni Mousa mailto:amousa@rockefeller.edu for more "\
            "information.\n" %NON_COMPLETED
    text4 = "\n\nPlease do not reply to this email.\nAll replies "\
            "should be addressed to Awni Mousa mailto:amousa@rockefeller.edu"
    
    html1 = """"""
    html2 = """"""
    text = text1
    if COMPLETED:
        text = text + text2
        html1 = """<br>
        The following sequencing data were successfully processed and completed:
                   <br>
                   %s<br>""" %COMPLETED
    if NON_COMPLETED:
        text = text + text3
        html2 = """<br>
        The following sequencing data were NOT processed:
                   <br>
                   %s<br>
        Please contact <a href="mailto:amousa@rockefeller.edu">Awni Mousa</a>
        for more information.<br>""" %NON_COMPLETED
        
    text = text + text4
    
    html = """\
    <html>
      <head></head>
      <body>
        <p>Dear <b>%s</b>,<br>
        <br>
        %s
        %s
        <br>
        </p>
      </body>
      <footer><font size="3">
        Please do not reply to this email.<br>
        All replies should be addressed to 
        <a href="mailto:amousa@rockefeller.edu">Awni Mousa</a></font>
      </footer>
    </html>
    """ %(contact,html1,html2)
    # Record the MIME types of both parts - text/plain and text/html.
    part1 = MIMEText(text, 'plain')
    part2 = MIMEText(html, 'html')
    # Attach parts into message container.
    # According to RFC 2046, the last part of a multipart message, in this case
    # the HTML message, is best and preferred.
    msg.attach(part1)
    msg.attach(part2)
    # Send the message via local SMTP server.
    s = smtplib.SMTP_SSL('smtp.rockefeller.edu')
    user = bioinformaticsRU.split('@')[0]
    password = 'f75-hTT-4Jw-y6nq'
    s.login(user,password)
    # sendmail function takes 3 arguments: sender's address, recipient's address
    # and message to send - here it is sent as one string.
    s.sendmail(bioinformaticsRU, userEmailArddress, msg.as_string())
    s.quit()


def runJobs(s,COMMAND,m,path,IDS,CONTACT,logger):
    with s:
        # don't attempt to run more that 100 jobs at a time
        CMDs = []
        while m > 100:
            FQ = COMMAND[1][:101]
            # get the fastq ids in command
            fq_ids = [i.split("_")[0] for i in FQ]
            fq_ids = sorted(set(fq_ids))
            CMDs+=[[' '.join(COMMAND[0]+ COMMAND[1][:101]+ COMMAND[2]),fq_ids]]
            [COMMAND[1].remove(Pr) for Pr in FQ[1:]]
            m = len(COMMAND[1][1:])
        # get the fastq ids in command
        fq_ids = [i.split("_")[0] for i in COMMAND[1][1:]]
        fq_ids = sorted(set(fq_ids))
        CMDs += [[' '.join(COMMAND[0] + COMMAND[1] + COMMAND[2]),fq_ids]]
        for CMD in CMDs:
            # run the command
            logger.info("\n\n"+CMD[0])
            os.chdir('/data4/amousa/processing')
            p = subprocess.Popen( CMD[0],stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  universal_newlines = True, 
                                  shell = True )
            outs, errs = p.communicate()
            if errs:
                logger.error( 'STDERR in Command:\n%s\n%s' %(CMD[0],errs))
            if outs:
                logger.info('STDOUT in Command:\n%s\n%s' %(CMD[0],outs))
            COMPLETED = []
            COMP = []
##### locate the completed files and send notification
            for ID,seqSampleID in IDS.items():
                seqSampleID = sorted(set(seqSampleID))
                for nfq in CMD[1]:
                    if nfq in seqSampleID:
                        COMPLETED += [nfq]
                if COMPLETED:
                    logger.info("\nID: "+ID+
                                "\nContact: "+CONTACT[ID]+
                                "\nfq_Ids:\n"+' '.join(COMPLETED)+"\n\n")
                    # check the integrity of the transferred files
                    for fqID in COMPLETED:
                        fqOUTS=glob.glob(os.path.join(
                                      '/mnt/heintz-bambi2/BAM',
                                            fqID+'_*/*'))
                        for inFile in fqOUTS:
                            if os.path.exists(inFile) and\
                               os.stat( inFile )[stat.ST_SIZE] != 0:
                                COMP += [fqID]
                    # send email notification with job name and link when done
                    if COMP:
                        COMP = sorted(set(COMP))
                        constructEmail(COMP,'',ID,CONTACT[ID])
                        logger.info("send notification\n"+
                                    "\nID: "+ID+
                                    "\nContact: "+CONTACT[ID]+
                                    "\nfq_Ids:\n"+' '.join(COMP)+"\n\n")
                        toMove = ' '.join(["@" + i + "_*" for i in COMP])
                        # change back to wDir
                        os.chdir( path )
                        subprocess.call("mv "+ toMove +" finished_jobs", 
                                        shell=True)
                    # empty placeholder to avoid redundant email notifications
                    COMPLETED = []
                    COMP = []                    


def download_GEO(SRA,ID='',wDir=os.getcwd(),cwDir=os.getcwd(),logger=''):
    import os
    import subprocess
    import glob
    import re
    # download SRA data file
    os.chdir(cwDir)
    # clean directory of any SRR with the same SRR ID
    old = glob.glob("*%s*" %re.sub(".sra","",SRA))
    for i in old:
        os.remove(i)
    CMD = 'wget -q ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/'+\
          'ByRun/sra/SRR/%s/%s/%s.sra' %(SRA[:6],SRA,SRA)
    p = subprocess.Popen( CMD,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                          universal_newlines = True, shell = True )
    outs, errs = p.communicate()
    fq = []
    paired = 'SE'
    if errs:
        if logger:
            logger.error( 'STDERR in Command:\n%s\n%s' %(CMD,errs))
        else:
            print("ERROR: " , CMD)
            print(errs)
    else:
        if outs:
            if logger:
                logger.info('STDOUT in Command:\n%s\n%s' %(CMD,outs))
            else:
                print("INFO: " , CMD)
                print(outs)
        # determine if the file is readable and whether the data is paired- or
        # single-end 
        CMD = "fastq-dump -I -X 1 -Z --split-spot %s.sra" %SRA
        p = subprocess.Popen( CMD,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                              universal_newlines = True, shell = True )
        outs, errs = p.communicate()
        if not outs:
            if logger:
                logger.error( 'STDERR in Command:\n%s\n%s' %(CMD,errs))
            else:
                print("ERROR: " , CMD)
                print(errs)
        else:
            first_read = outs.strip().split()
            first_read = [i for i in first_read if '@'+SRA in i]
            pairs = [int(i.split('.')[-1]) for i in first_read]
            if len(pairs) == 1 and pairs[0] == 1:
                if logger:
                    logger.info("%s: single-end" %SRA)
                else:
                    print("%s: single-end" %SRA)
                CMD = "fastq-dump -I --gzip --outdir %s %s.sra" %(cwDir,SRA)
                p = subprocess.Popen( CMD,stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True,shell=True )
                outs, errs = p.communicate()
                if errs:
                    if logger:
                        logger.error('STDERR in Command:\n%s\n%s' %(CMD,errs))
                    else:
                        print("ERROR: " , CMD)
                        print(errs)
                else:
                    if outs:
                        if logger:
                            logger.info('STDOUT in Command:\n%s\n%s' %(CMD,outs))
                        else:
                            print('STDOUT in Command:\n%s\n%s' %(CMD,outs))
                    fq += [[SRA+'.fastq.gz',
                            "_".join([ID,SRA+'.fastq.gz']).strip("_")]]
            elif len(pairs) == 2 and pairs[0] != pairs[1]:
                if logger:
                    logger.info("%s: paired-end" %SRA)
                else:
                    print("%s: paired-end" %SRA)
                # split SRA into 2 separate mates
                CMD = "fastq-dump -I --gzip --outdir %s --split-files %s.sra" \
                      %(cwDir,SRA)
                p = subprocess.Popen( CMD,stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True,shell=True )
                outs, errs = p.communicate()
                if errs:
                    if logger:
                        logger.error('STDERR in Command:\n%s\n%s' %(CMD,errs))
                    else:
                        print('STDERR in Command:\n%s\n%s' %(CMD,errs))
                else:
                    if outs:
                        if logger:
                            logger.info('STDOUT in Command:\n%s\n%s' %(CMD,outs))
                        else:
                            print('STDOUT in Command:\n%s\n%s' %(CMD,outs))
                    fq += [[SRA+'_1.fastq.gz',
                            "_".join([ID,SRA+'_R1.fastq.gz']).strip("_")]]
                    fq += [[SRA+'_2.fastq.gz',
                            "_".join([ID,SRA+'_R2.fastq.gz']).strip("_")]]
            else:
                if logger:
                    logger.error( 'STDERR in Command:\n%s\n%s' %(CMD,outs))
                else:
                    print('STDERR in Command:\n%s\n%s' %(CMD,outs))
            # cleanup binaries
            os.remove("%s.sra" %SRA)
    # rename the fastq files
    if fq:
        FQ = []
        for IN,OUT in fq:
            try:
                os.rename(cwDir+"/"+IN , cwDir+"/"+OUT)
                if logger:
                    logger.info('External Fastq File %s Renamed Successfully' \
                                %SRA)
                else:
                   print('External Fastq File %s Renamed Successfully' %SRA)
                FQ += [OUT]
            except Exception as e:
                if logger:
                    logger.error('External Fastq File Rename: %s\n%s' \
                                 %(SRA,str(e)))
                else:
                    print('External Fastq File Rename: %s\n%s' %(SRA,str(e)))
        if len(FQ) == 2:
            paired = 'PE'
        fq = ','.join(FQ)
    os.chdir(wDir)
    return [fq,paired]
    

def printSi():
    L = [32,45,45,45,45,45,45,45,45,45,45,45,45,13,10,124,32,65,119,110,105,32,
         77,111,117,115,97,32,45,32,66,105,111,105,110,102,111,114,109,97,116,
         105,99,115,32,83,112,101,99,105,97,108,105,115,116,13,10,124,32,84,104,
         101,32,82,111,99,107,101,102,101,108,108,101,114,32,85,110,105,118,101,
         114,115,105,116,121,13,10,32,45,45,45,45,45,45,45,45,45,45,45,45]
    I = ''.join([chr(i) for i in L])
    return(I)
