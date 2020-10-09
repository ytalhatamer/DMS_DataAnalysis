# SaturationMutagenesisAnalysis
Data Analysis scripts that are used to analyze TolC saturation mutagenesis pipeline outlined in Tamer et al 2020. doi: 

For detailed methods please refer to  Tamer et al 2020. doi: 

In short, analysis pipeline has four steps.

1. Merging the paired end reads using Flash tool(run-Flash.sh). https://ccb.jhu.edu/software/FLASH/
2. Clips the primer sequences provided in WTSequencesPrimers.xlsx by user (run-ClipOrganize.sh).
- Writes the sequence in between two primers in to \*.dna file
3. Aligns clipped sequence to WT sequence and finds mismatches and makes up a lookup table for each mutation count (run-CallMutations.sh)
4. Calculates frequency and enrichment of each mutation along with relative fitness value. (PopulateDataCalculateEnrichment.py)

