# ARGfams: 

Fast and robust identification of antibiotic resistance genes (ARGs) from genomic and metagenomic assemblies of high-throughput DNA sequencing. The ARGfams is a high-quality, manually curated and structured subdatabase of profile Hidden Markov Models (HMM) for ARGs annotation.

This subdatabase consists of HMM models for ARGs  were built from Pfam(v34.0) and TIGRFAMs(v15.0), based on string match in their functional annotations to one of the keywords (Details are in the citations below). [Download](https://doi.org/10.6084/m9.figshare.21610416.v1 )

For ARG-like ORFs with domain bit-score of best hits greater than 50 were finalized as an ARG. The identified ARGs were then classified into 12 types, including aminoglycoside, beta-lactam, bleomycin, chloramphenicol, daunorubicin, macrolide-lincosamide-streptogramin (MLS), multidrug, quinolone, tetracycline, trimethoprim, vancomycin, unclassified, followed by further classification  into 161 subtypes (Details are in the citations below).

We provide two similar ARG annotation strategies but slightly different in annotation speed.

The first strategy is a two-step scan. The advantage of this annotation strategy is that when you update the version of the large database of the protein HMM model, you can continuously update and optimize the sub-ARG database by manually verifying the difference annotation results in the two-step scan.

However, large protein hmm model databases are not updated frequently. So if the large database has not been updated within a certain period, using the sub-ARG databases for a one-step scan can obtain reliable annotation results in a shorter time, which is the second strategy.

---

### Development Record

ARGfams created by Feng Ju in Oct 2018
ARGfams v0.1 - created by Xinyu Huang and Guoqing Zhang based and Pfam (v34.0) and TIGRFAMs (v15.0)


### Dependence

python: >=3.6
HMMER3: >=3.3.2

Update Notes
ARGfams v0.1 was constructed based on the lastest version of Pfam (v34.0) and TIGRFAMs (v15.0), and will be regularly updated.

Dependence
python: >=3.6
HMMER3: >=3.3.2

### Usage
DESCRIPTION ARGfams version: 0.5.0 Detailed introducion

optional arguments:
**-i INPUT_FILE, --input INPUT_FILE**
the input file ORFs.faa  

**-o [OUTPUT_FILE_NAME], --output [OUTPUT_FILE_NAME]**
the outputfile prefix name: eg. PRIFEX.tlout  

Required arguments:  
**-db ARGfams_Database**
ARGfams_Database; Default Antibiotic Resistance Genes  

**-DB synthesis_database**
Synthesis database, Default Pfam and Tigrfam

**{--cut_ga,--cut_nc,--cut_tc}**
hmm type; chose from {--cut_ga, --cut_nc, --cut_tc default:cut_ga  

**-n N, --nproc N**
The number of CPUs to use for parallelizing the mapping [default 1]  

Other arguments:  

**--check **
Only checks if the Default ARG DB is installed and installs it if not.  

**-v, --version**
Prints the current MetaPhlAn version and exit  

**-h, --help**
show this help message and exit



```bash
# One-step scan
python ARGfams_OneStep_v0.1.py -i <protein_c_ORFs>.faa -o <OUTPUT_NAME> -db ARGfams_V0.1/ARGfams_v0.1.hmm -n 2

# Two-step scan
python ARGfams_TwoSteps_v0.1.py -i <protein_c_ORFs>.faa -o <OUTPUT_NAME> -db ARGfams_V0.1/ARGfams_v0.1.hmm -DB Pfam-Tigrfam.hmm -n 2
```



**Citation:** He L, Huang X, Zhang G, Yuan L, Shen E, Zhang L, et al. Distinctive signatures of pathogenic and antibiotic resistant potentials in the hadal microbiome. Environmental Microbiome. 2022;17(1):19. https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-022-00413-5
