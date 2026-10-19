# CARP snakemake with CENSOR WU-BLAST RepeatMasker
This is a Snakemake pipeline for CARP that uses CENSOR, WU-BLAST, and RepeatMasker. Due to software licensing restrictions, this project does not provide a container. Users are required to install the necessary software dependencies themselves.

```bash
conda create -n snkm9 -c bioconda snakemake=9.27.0
conda activate snkm9
snakemake --version
```


```bash
git clone https://github.com/Luan-Zhong/carp_snkm_with_censor_wublast_repeatmasker
```

### Download databases
Make sure to move the downloaded databases to the lib directory and update the file names in the "# Reference libraries section" of `config.yaml`
- RepBase FASTA format</br>
  Download RepBase differnt releases from [Dfam Repbase Archive](https://dfam.org/repbase/releases).
- Swissprot library </br>
  Download uniprot protein data set from [UniProt download](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
  ```bash
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  gzip -d uniprot_sprot.fasta.gz
  ```
- Download GB_TE dataset </br>
  First install EDirect from EDirect download page [Entrez Direct: E-utilities on the Unix Command Line](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Installation)).

  Then download the data, this taking a bit time due to the file will be large
  ```bash
  esearch -db protein -query "reverse transcriptase or transposon or repetitive element or RNA-directed DNA polymerase or pol protein or non-LTR retrotransposon or mobile gag or polymerase (pol)" | efetch -format fasta > GB_TE.2026-10-18.fa
  ```
  if this not working try
  ```bash
  NQUIRE_HELPER=wget esearch -db protein -query "reverse transcriptase OR transposon OR repetitive element OR RNA-directed DNA polymerase OR pol protein OR non-LTR retrotransposon OR mobile gag OR polymerase (pol)" | NQUIRE_HELPER=wget efetch -format fasta > GB_TE.2026-10-18.fa
  ```
- Download retrovirus datasets </br>
  Download retrovirus datasets from your browser: [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Viruses,taxid:11632)


## Quick start

```bash
# Dry run
conda activate snkm9
snakemake --snakefile Snakefile_censor_wublast_repeatmasker --use-singularity --printshellcmds --configfile config_censor_wublast_repeatmasker.yaml --cores 4 -n
```
NB: You need a container to run this, or if you set up your local environment, then you can discard the `--use-singularity`
