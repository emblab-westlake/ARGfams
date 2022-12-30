# ARGfams: 

ARGfams leverages the power of hmmscan search against a high-quality, manually curated, and structured sub-database of profile Hidden Markov Models (HMM) for fast annotation of Antibiotic Resistance Genes (ARGs) from genomic and metagenomic assemblies of high-throughput DNA sequencing data. This annotation tool identifies protein domains of known ARGs from genome or metagenome-assembled open reading frames (ORFs) or protein-coding genes (PCGs) and excels classic sequence alignment-based approaches for predicting relative remote homologues of known ARGs from environmental microbiome (e.g., soil, water, and sediment) in which new or less-homologous ARGs commonly occur.

**Database resource:** The structured sub-database of ARGs consists of 197 HMM models extracted from the full database of Pfam(v34.0), TIGRFAMs(v15.0) and Resfams (Full - v1.2), based on string match in their functional annotations to one of the indicative keywords of ARGs (Dataset S2 of Environmental Microbiome. 2022;17(1):19), followed by manual validation. The ARGs in the sub-database (v0.1) were classified into 11 types based on the class of antibiotics to which they confer resistance, including aminoglycoside, beta-lactam, bleomycin, chloramphenicol, daunorubicin, macrolide-lincosamide-streptogramin (MLS), multidrug, quinolone, tetracycline, trimethoprim, vancomycin, followed by further classification into 161 subtypes based on protein families. More details on the structured sub-database are available in the Methods of the Citation below).

**Annotation strategy:** two alternative strategies of ARG annotation that in principle generated the same results.

Strategy 1 (default): ARGfams (ARGfams_v0.1.py) performs one-step scan of query ORF sequences against the ARGfams sub-database (-db), and significant best hits with a domain bit-score greater than 50 are predicted as ARG.

Strategy 2 (optional): ARGfams (ARGfamsPlus_v0.1.py) performs two-step scan of query ORF sequences against the structured sub-database ARGfams (-db) and user-defined HMM database (-DB). The 2nd scan outputs are compared against those of 1st scan outputs to check whether certain models in the user-defined database (-DB) might generate higher-confidence alignments of ARGs. If yes, the users can expand and update the current version of structured sub-database of ARGfams by incorporating these models from user-defined database.

---

### Development Record

ARGfams conceived and initially created by Feng Ju in Oct 2018. Structured sub-database ARGfams v0.1.hmm created by Xinyu Huang and Guoqing Zhang based on Pfam (v34.0), TIGRFAMs (v15.0) and Resfams(Full - v1.2)


### Dependence

python: >=3.6  
HMMER3: >=3.3.2



### Usage
DESCRIPTION ARGfams version: 0.1 Detailed introduction

#### optional arguments:
**-i INPUT_FILE, --input INPUT_FILE**  
the input file ORFs.faa  

**-o [OUTPUT_FILE_NAME], --output [OUTPUT_FILE_NAME]**  
the outputfile prefix name: eg. PRIFEX.tlout   

#### Required arguments:  
**-db ARGfams_Database**  
ARGfams_Database; Default Antibiotic Resistance Genes  

**-DB synthesis_database**  
Synthesis database, user-defined database

**{--cut_ga,--cut_nc,--cut_tc}**  
hmm type; chose from  --cut_ga, --cut_nc, --cut_tc [default: cut_ga] 

**-n N, --nproc N**  
The number of CPUs to use for parallelizing the mapping [default 1]  

#### Other arguments:  

**--check**  
Only checks if the Default ARG DB is installed and installs it if not.  

**-v, --version**  
Prints the current version

**-h, --help**  
show this help message



```bash
# Strategy 1 (default) with One-step scan (faster)
python ARGfams_v0.1.py -i <INPUT_FILE> -o <OUTPUT_Prefix> -db ARGfams_V0.1/ARGfams_v0.1.hmm -n 2

# Strategy 2 (optional) with two-step can (slower)
python ARGfamsPlus_v0.1.py -i <INPUT_FILE> -o <OUTPUT_Prefix> -db ARGfams_V0.1/ARGfams_v0.1.hmm -DB user_defined.hmm -n 2
```



**Citation:** He L, Huang X, Zhang G, Yuan L, Shen E, Zhang L, et al. Distinctive signatures of pathogenic and antibiotic resistant potentials in the hadal microbiome. Environmental Microbiome. 2022;17(1):19. https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-022-00413-5
