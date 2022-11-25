# Wastewater-SARS-CoV-2-Analysis-Pipeline
SARS-CoV-2 from Wastewater samples Analysis Pipeline
Wastewater SAR-CoV-2 analysis is a bash-based bioinformatics pipeline for the analysis of either long-read whole genome sequencing data of DNA extracted from wastewater. It was developed for SARS-CoV2 and its variants.

## The process includes the following:
1. Alignment of reads to the reference via BWA
2. Processing of alignment results via samtools
3. Detection of variant positions with ivar
4. Generation of consensus and sequencing statistics
5. Determine composition of variants via freyja

## Dependencies

User provided:
- [BWA](https://github.com/lh3/bwa)
- [SamTools](https://github.com/samtools/samtools)
- [ivar](https://github.com/andersen-lab/ivar)
- [pangolin](https://github.com/cov-lineages/pangolin)
- [freyja](https://github.com/andersen-lab/Freyja)


## Usage
Install all the required tool for the pipeline and run "bash Wastewater_SARS-CoV-2_Script.sh"
