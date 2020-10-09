from multiprocessing import Pool, Lock
import os, sys, glob
import pandas as pd
import numpy as np
from pathlib import Path
import globalvars
from ClipOrganizeSublibraries import ClipOrganizeReads
from DNACallMutations import CallMutations

def FlashMe(inp):
    treatmentname, fwd, rev, inputfolder = inp
    flashpath = workpath / 'readjoiners/flash'
    os.chdir(inputfolder)
    print('Entering to input folder: ',os.getcwd())

    treatmentfolder=inputfolder / treatmentname
    if not treatmentfolder.exists():
        os.makedirs(treatmentfolder)

    cmd=flashpath.as_posix().replace(' ','\ ') + ' -z -M 100 -m 40 -o '+\
        treatmentname+'/'+treatmentname+' '+fwd+' '+rev+' | tee '+treatmentname+'.log'
    print(cmd)
    os.system(cmd)
    os.chdir(workpath)
    print('Back to work path : ',os.getcwd())
    
def SymbLinkOutput(inp):
    treatmentname, fwd, rev, inputfolder = inp
    seqDate=inputfolder.stem
    outputfolder=analysesfolder / seqDate
    tfilename=treatmentname+'.extendedFrags.fastq.gz'
    targetfile=inputfolder / treatmentname / tfilename
    dfilename=treatmentname+'_flashed.fastq.gz'
    destinationfile=outputfolder / 'dnafiles' / dfilename
    destinationfile.symlink_to(targetfile)


if __name__=='__main__':


    print('In this folder :',os.getcwd())

    num_cores=int(sys.argv[1])
    func2run = sys.argv[2]
    seqDate=sys.argv[3]
    scriptpath=globalvars.scriptpath
    workpath=globalvars.workpath
    datapath=globalvars.datapath
    analysesfolder=globalvars.analysesfolder
    inputfolder=analysesfolder / seqDate
    dnafilesfolder= inputfolder / 'dnafiles'
    seqdatafolder = datapath / 'SequencingData' / seqDate

    if not inputfolder.exists():
        os.makedirs(inputfolder)
    if not dnafilesfolder.exists():
        os.makedirs(dnafilesfolder)

    samplesTablefile = scriptpath / 'SamplesTable.xlsx'
    dfSequencing = pd.read_excel(samplesTablefile, 'SequencingRunParams') \
        .set_index('SequencingDate')
    dfSamples = pd.read_excel(samplesTablefile, seqDate).set_index('Index')

    primers2Use = dfSequencing['PrimerSheet2Use'][int(seqDate)]


    ReadStatFile= inputfolder / 'ReadStats.csv'
    ReadStatA = inputfolder / 'ReadStatsA.csv'
    ReadStatB = inputfolder / 'ReadStatsB.csv'
    if not ReadStatFile.exists():
        samplenames=dfSamples['SampleName'].values  #[i.split('_')[0] for i in dfSamples['FileForward'].values]
        dfReadStats=pd.DataFrame(np.nan,columns=['SampleType','SetNum','#InFile','#lenRead<lenRef',\
                                            '#PrimerMatch','#lenRead!=lenRef',\
                                            '#TranslationProblem','#ValidRead','#TrueWT','#SynonymousWT',\
                                            '#TotalWT','#TotalMutants'],index=samplenames)
        dfReadStats.index.name='SampleName'
        dfReadStats['SampleType']=[i.split('-Set-')[0] for i in dfReadStats.index.values]
        dfReadStats['SetNum']=[i.split('-Set-')[1] for i in dfReadStats.index.values]
        dfReadStats.to_csv(ReadStatFile,',')
    ## Number of cores to parallelize..
    print( 'Total Number of Cores used: ' + str(num_cores) )

    p=Pool(processes=num_cores)#,initializer=init,initargs=(globalvars.ll,))
    joblist=[]
    if func2run=='Flash':
        for i in dfSamples.index:
            treatment=dfSamples['SampleName'][i]#['FileForward'][i].split('_')[0]
            fwd=dfSamples['FileForward'][i]
            rev=dfSamples['FileReverse'][i]
            joblist.append((treatment,fwd,rev,seqdatafolder))
        p.map(FlashMe,joblist)
        p.close()
        p.join()
        print('Flash job completed!')
        [SymbLinkOutput(i) for i in joblist]
    elif func2run=='ClipOrganize':
        for i in dfSamples.index:
            name=dfSamples['SampleName'][i]#dfSamples['FileForward'][i].split('_')[0]
            print(i,len(Path(dfSamples['FileForward'][i]).suffix))
            if Path(dfSamples['FileForward'][i]).suffix.strip() == '.gz' :
                filename=name+'_flashed.fastq.gz'
            else:
                filename=name+'_flashed.fastq'
            print(i, filename)
            joblist.append((inputfolder,filename,primers2Use))

        results=p.map(ClipOrganizeReads,joblist)
        p.close()
        p.join()
        finaltable=pd.concat(results,axis=0).dropna(thresh=4)
        print(finaltable)
        finaltable.to_csv(ReadStatA)

    elif func2run=='CallMutations':
        for i in dfSamples.index:
            sample=dfSamples['SampleName'][i]#dfSamples['FileForward'][i].split('_')[0]
            dnafilename=sample + '.dna'
            filename=dnafilesfolder / dnafilename
            joblist.append((inputfolder,filename,primers2Use))
        results=p.map(CallMutations,joblist)
        p.close()
        p.join()
        finaltable = pd.concat(results, axis=0).dropna(thresh=4)
        print(finaltable)
        ReadStatB= inputfolder / 'ReadStatsB.csv'
        finaltable.to_csv(ReadStatB)
