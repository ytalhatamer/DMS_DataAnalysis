#!/usr/bin/env python3
import numpy as np
import pandas as pd
import globalvars

#cmd='python '+scriptpath+'/DNACallMutations.py '+' '+inp+' '+primers2Use
def CallMutations(inp):
    scriptpath=globalvars.scriptpath
    inputfolder, filename , primers2Use =inp
    codon2aa=globalvars.codon2aa

    name=filename.stem
    setnum=int(name.split('-Set-')[-1])
    primersfile = scriptpath / 'WTSequencesPrimers.xlsx'
    WT = pd.read_excel(primersfile, primers2Use).set_index('SetNum')
    WT_sequence=WT['Sequence'][setnum].strip().replace('\ufeff','')
    WTcodons=[WT_sequence[i:i+3] for i in range(0,len(WT_sequence),3)]

    print("\nFile : "+filename.as_posix()+" is running ... \nSetNum= "+str(setnum))


    print('Script Directory: ',scriptpath,\
          ' Input Filename:',filename,\
          'Sequence Type: ', primers2Use)

    startresidue=int(WT['StartResidues'][setnum])
    endresidue=int(WT['EndResidues'][setnum])
    orderedaalist=['X','R','K','D','E','Q','N','H',\
                   'S','T','Y','C','W','A','I','L',\
                   'M','F','V','P','G','SynWT']

    counts=pd.DataFrame(0.0, index=range(startresidue,endresidue),columns=orderedaalist)

    for i,j in zip(WTcodons,range(startresidue,endresidue)):
        counts[codon2aa[i]][j]=np.nan

    trueWT=0
    synWT=0

    m=0
    dnafile2read= inputfolder / filename
    with open(dnafile2read,'rt') as reader:
        counter=0
        for line in reader:
            seq=line.strip()
            counter+=1
            if seq==WT_sequence:
                trueWT+=1
                continue
            dS=[]
            dN=[]
            codons=[seq[i:i+3] for i in range(0,len(seq),3)]
            for query,reference,residnum in zip(codons,WTcodons,range(startresidue,endresidue)):
                if query==reference:
                    continue
                else:
                    if codon2aa[query]==codon2aa[reference]:
                        dS.append([codon2aa[reference],residnum,reference,query])
                    else:
                        dN.append([reference,residnum,query,codon2aa[query]])
            if not dN:
                if not dS:
                    trueWT+=1
                elif len(dS)==1:
                    counts['SynWT'][dS[0][1]]+=1
                    synWT+=1
            else:
                if len(dN)==1:
                    if not dS:
                        counts[dN[0][3]][dN[0][1]]+=1

    counts.index.name='Residues'

    outputfile=name+'-Counts.csv'
    outputloc= inputfolder / outputfile
    counts.to_csv(outputloc,',',na_rep=np.nan)
    readstatsfile=inputfolder / 'ReadStats.csv'
    readstats = pd.read_csv(readstatsfile).set_index('SampleName')
    readstats['#TrueWT'][name] = trueWT
    readstats['#SynonymousWT'][name] = synWT
    readstats['#TotalWT'][name] = trueWT + synWT
    readstats['#TotalMutants'][name] = counts.sum().sum()

    print('Mutation counting for '+name+' is finished!')
    return readstats