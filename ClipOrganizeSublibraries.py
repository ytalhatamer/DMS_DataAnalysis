import re, sys
import pandas as pd
import gzip
global readcounter, wt
from pathlib import Path

import globalvars
import readFastQ

# cmd='python '+scriptpath+'/ClipOrganizeSublibraries.py '+inputfolder+' '+filename+' '+primers2Use
def ClipOrganizeReads(inp):
    scriptPath=globalvars.scriptpath
    inputfolder, filename , primers2Use = inp
    seqFileName=inputfolder / 'dnafiles' / filename
    name=filename.split('_')[0]
    dnafilename=name+'.dna'
    dnaFileName=inputfolder / 'dnafiles' / dnafilename

    print('Script Directory: ',scriptPath,\
          '\nInput Foldername: ',inputfolder,
          '\nOutput Filename:',name,\
          '\nSequence Filename:',filename,\
          '\nSequence Type: ', primers2Use)

    setnum=int(name.split('-Set-')[1])
    print("\nFile : "+seqFileName.as_posix()+" is running ... \nSetNum= "+str(setnum))
    

    def checkzipped(filename):
        if filename.as_posix()[-2:]=='gz':
            print('gzip.open command opens the sequence file')
            return gzip.open(filename,'rt')
        else:
            print('open command opens the sequence file')
            return open(filename,'rt')
    def curatesequence(seq,qsc):
        curseq=''
        if len(wt)==len(seq):
            for i in range(len(wt)):
                if i==(len(wt)-1):
                    if seq[i]==wt[i]:
                        curseq+=seq[i]
                    else:
                        if qsc[i]<10:
                            if seq[i-1]==wt[i-1]:
                                curseq+=wt[i]
                            else:
                                curseq+=seq[i]
                        else:
                            curseq+=seq[i]
                elif i==0:
                    if seq[i]==wt[i]:
                        curseq+=seq[i]
                    else:
                        if qsc[i]<10:
                            if seq[i+1]==wt[i+1]:
                                curseq+=wt[i]
                            else:
                                curseq+=seq[i]
                        else:
                            curseq+=seq[i]

                else:
                    if seq[i]==wt[i]:
                        curseq+=seq[i]
                    else:
                        if qsc[i]<10:
                            if seq[i-1]==wt[i-1] and seq[i+1]==wt[i+1]:
                                curseq+=wt[i]
                            else:
                                curseq+=seq[i]
                        else:
                            curseq+=seq[i]
        else:
            return None
        return curseq

    def translate(dna):
        protein=''
        codon2aa=globalvars.codon2aa
        if len(dna)%3==0:
            for i in range(0,len(dna),3):
                if dna[i:i+3] in codon2aa:
                    protein+=codon2aa[dna[i:i+3]]
                else:
                    return None
        else:
            return None
        return protein

    primersfile=scriptPath / 'WTSequencesPrimers.xlsx'
    WT=pd.read_excel(primersfile, primers2Use).set_index('SetNum')
    wt=WT['Sequence'][setnum].strip().replace('\ufeff','')
    primerfwd=WT['Fwd-Pattern'][setnum]
    primerrev=WT['Rev-Pattern'][setnum]
    print ('FWD primer : '+primerfwd)
    print ('REV primer : '+primerrev)
    readcounter=0
    validreadcounter=0
    numlenproblem=0
    numtranslationproblem=0
    numreadwritten=0
    shortread=0

    with open(dnaFileName,'w') as txtfile :
        with checkzipped(seqFileName) as currFile:
            AllReads=readFastQ.read_fastq('Author',currFile)
            for i in AllReads:
                readcounter+=1
                if len(i[1])<len(wt):
                    shortread+=1
                    continue
                qsc_now=[ord(j)-33 for j in i[2]]
                match=re.search(r'(?P<before>\S*'+primerfwd+')(?P<on>\S+)(?P<after>'+\
                                                  primerrev+'\S*)',str(i[1]))

                if match is not None:
                    if len(match.groups())==3:
                        validreadcounter+=1
                        strTxt=match.groups()
                        lenstrTxt=[len(strTxt[0]),len(strTxt[1])+len(strTxt[0]),\
                                   len(strTxt[2])]

                        seq_clipped=strTxt[1]
                        qsc_clipped=qsc_now[lenstrTxt[0]:-len(strTxt[-1])]
                        curseq=curatesequence(strTxt[1],qsc_clipped)
                        if curseq is not None:
                            proteinseq=translate(curseq)
                            if proteinseq is not None:
                                txtfile.write(curseq+'\n')
                                numreadwritten+=1
                            else:
                                numtranslationproblem+=1
                                continue
                        else:
                            numlenproblem+=1
                            continue
                else:
                    continue

    ReadStatFile = inputfolder / 'ReadStats.csv'
    readstats=pd.read_csv(ReadStatFile).set_index('SampleName')
    readstats['SampleType'][name] = name.split('-Set-')[0]
    readstats['SetNum'][name] = int(name.split('-Set-')[1])
    readstats['#InFile'][name]=readcounter
    readstats['#lenRead<lenRef'][name] = shortread
    readstats['#PrimerMatch'][name] = validreadcounter
    readstats['#lenRead!=lenRef'][name] = numlenproblem
    readstats['#TranslationProblem'][name] = numtranslationproblem
    readstats['#ValidRead'][name] = numreadwritten

    print('*'*40)
    print('Sample Name: ',name)
    print ('Total Read Count: ', readcounter)
    print ('Number of Reads with length shorter than WT: ',shortread)
    print ('Total Read Number with Primer Match: ', validreadcounter)
    print ('Number of Reads with length different than WT: ', numlenproblem)
    print ('Number of Reads with translation problem: ', numtranslationproblem)
    print ('Number of Reads written in DNA file: ', numreadwritten)
    print('*'*40)
    return readstats
