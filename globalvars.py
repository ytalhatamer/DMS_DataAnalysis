from pathlib import Path
global scriptpath, workpath, datapath, analysesfolder, codon2aa

workpath = Path('/work/greencenter/ytamer/')
scriptpath = workpath / 'TolC_SM_Analysis_Scripts'
datapath = workpath / 'TolC_SM_Data'
analysesfolder= datapath / 'StandardizedAnalyses'
codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
            "TAA": "X", "TAC": "Y", "TAG": "X", "TAT": "Y",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
            "TGA": "X", "TGC": "C", "TGG": "W", "TGT": "C",
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

sepFac=0.1 # Min distinguishable read count to add when calculating enrichment or fitness
