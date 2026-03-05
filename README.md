# fast-GWAS
*Note that the project is currently a work in progress. The description below is not (yet) accurate. Files of interest at this time are: - GWAS_draft2.ipynb, which reads a genotype, phenotype, and optional covariate file and runs GWAS by linear regression
- Augment_data.ipynb, which was used to subsample 1000 Genomes data and simulate our phenotypic controls


## TODO
- Implement MP (accelerate over multiple cores)
- Convert to python script and support command line arguments 
- Generate larger subsample
- Benchmark runtime against PLINK for small vs large datasets
- Analyze/interpret results of positive vs negative control for our program versus PLINK


## Installation
This tool can be installed and run from the command line directly. A singularity definition file and requirements.txt are included, but dependency versioning is not strict. See next section for examples on how to run the script.
'''
John_Doe@ubuntu:~$ git clone https://github.com/Siddharth-Gaywala/GWAS.git .
'''

## Arguments
There are two mandatory arguments:  
-g: path to Genotype file, .vcf format  
-p: path to Phenotype file, .phen format  

There are two optional arguments:  
-c: Covariates file, csv/tsv format  
-mp: Integer value (>1) to denote number of cores for mp. If no value is set, only one core will be used.  

## Example run
'''
John_Doe@ubuntu:~$ /Users/John_Doe/FGWAS/FGWAS.py -g /mnt/data/genotypes.vcf -p /mnt/data/phenotypes.phen -c /mnt/data/covariates.cov -mp 4
'''
