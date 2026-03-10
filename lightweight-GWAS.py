import numpy as np
import pandas as pd
from cyvcf2 import VCF
import os
from scipy import stats
from qqman import qqman
import matplotlib.pyplot as plt
import argparse
import os
from sklearn.decomposition import PCA
import tracemalloc
import time

def gwas_linear(G, y):
    """
    G: (n_samples, n_snps)
    y: (n_samples,)
    Returns: beta, pvals
    """
    n = G.shape[0]

    # Center variables (important!)
    G_centered = G - G.mean(axis=0)
    y_centered = y - y.mean()

    # Compute beta for each SNP
    numerator = np.sum(G_centered * y_centered[:, None], axis=0)
    denominator = np.sum(G_centered**2, axis=0)
    beta = numerator / denominator

    # Predicted values
    y_hat = G_centered * beta

    # Residual variance
    rss = np.sum((y_centered[:, None] - y_hat)**2, axis=0)
    sigma2 = rss / (n - 2)

    # Standard error
    se = np.sqrt(sigma2 / denominator)

    # t-statistic
    t_stat = beta / se

    # Two-sided p-value
    pvals = 2 * stats.t.sf(np.abs(t_stat), df=n-2)

    return beta, pvals

def gwas_linear_cov(G, y, C):
    """
    G: (n_samples, n_snps)
    y: (n_samples,)
    C: (n_samples, n_covariates)

    Returns:
        beta (n_snps,)
        pvals (n_snps,)
    """

    n = G.shape[0]

    # Add intercept to covariates
    C = np.column_stack((np.ones(n), C))

    # Precompute projection matrix components
    CtC_inv = np.linalg.inv(C.T @ C)
    P = C @ CtC_inv @ C.T   # projection onto covariates

    # Residualize phenotype
    y_res = y - P @ y

    # Residualize genotypes
    G_res = G - P @ G

    # Now run simple regression on residuals
    numerator = np.sum(G_res * y_res[:, None], axis=0)
    denominator = np.sum(G_res**2, axis=0)
    beta = numerator / denominator

    # Residual variance
    y_hat = G_res * beta
    rss = np.sum((y_res[:, None] - y_hat)**2, axis=0)

    df = n - C.shape[1] - 1
    sigma2 = rss / df

    se = np.sqrt(sigma2 / denominator)

    t_stat = beta / se
    pvals = 2 * stats.t.sf(np.abs(t_stat), df=df)

    return beta, pvals

def gwas_fxn_cov(vcf_file: str, phenotype_file: str, out_name: str, covariates_matrix = None):
    vcf = VCF(vcf_file)
    
    # read phenotypes
    phen_df = pd.read_csv(phenotype_file, delim_whitespace=True, header=None)
    phen_df.columns = ["FID", "IID", "PHENO"]
    
    # read genotypes
    sample_ids = np.array(vcf.samples)
    n_samples = len(sample_ids)
    print("Number of individuals in VCF:", n_samples)
    
    genotypes = []
    chroms = []
    positions = []

    for variant in vcf:
        # variant.genotypes gives [ [0,1,True], ... ]
        # first two entries are alleles
        gt = np.array([g[0] + g[1] if g[0] != -1 else np.nan 
                       for g in variant.genotypes])
        genotypes.append(gt)
        chroms.append(variant.CHROM)
        positions.append(variant.POS)

    G = np.array(genotypes).T  # transpose to get individuals x SNPs
    
    print("# individuals x # SNPs", flush=True)
    print(G.shape, flush=True) # individuals x SNPs
    
    # Align Phenotypes
    phenotypes = phen_df["PHENO"].values
    phenotype_ids = phen_df["IID"].values
    id_to_pheno = dict(zip(phenotype_ids, phenotypes))

    aligned_phenotypes = np.array([id_to_pheno[s] for s in sample_ids])
    aligned_phenotypes = np.array(aligned_phenotypes, dtype=float)
    
    if covariates_matrix is None:
        beta, pvals = gwas_linear(G, aligned_phenotypes)
    else:
        beta, pvals = gwas_linear_cov(G, aligned_phenotypes, covariates_matrix)
    
    assoc_df = pd.DataFrame({
        "CHR": chroms,
        "BP": positions,
        "BETA": beta,
        "P": pvals
    })
    assoc_df["CHR"] = pd.to_numeric(assoc_df["CHR"], errors="coerce")
    assoc_df = assoc_df.sort_values(["CHR", "BP"]).reset_index(drop=True)

    print(assoc_df.head(), flush=True)
    
    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    fig.set_size_inches((15, 5))
    qqman.manhattan(assoc_df, ax=ax0)
    qqman.qqplot(assoc_df, ax=ax1)
    
    print("Saving visualization to: " + out_name+"_manhattan.png", flush=True)
    fig.savefig(out_name+"_manhattan.png", dpi=150, bbox_inches="tight")
   
    print("Saving analysis results to: " + out_name+"_analysis_results.csv", flush=True)
    assoc_df.to_csv(out_name+"_analysis_results.csv", index=False)
    
if __name__ == "__main__":
    # runtime and memory profiling
    tracemalloc.start()
    start = time.time()

    parser = argparse.ArgumentParser(description="Run simple-GWAS from CLI")
    parser.add_argument("--vcf", type=str, required=True, help="Path to valid VCF file")
    parser.add_argument("--phenotype", type=str, required=True, help="Path to valid phenotype file")                    
    parser.add_argument("--covariates", type=str, required=False, default=None, help="Path to covariates file")
    parser.add_argument("--pca_covariates", type=str, required=False, default="False", help="Boolean(as a string) to indicate whether PCA should be run on genotypes as an alternative covariates method")
    parser.add_argument("--out", required=True, help="Prefix for analysis results")
    args = parser.parse_args()
    
    if args.pca_covariates == "False" or args.pca_covariates == "false":
        to_pca = False
    elif args.pca_covariates == "True" or args.pca_covariates == "true":
        to_pca = True
    else:
        raise ValueError("String passed to pca_covariates is not false/False or true/True")
    
    # Read in covariates file if it is provided
    covariates_matrix = args.covariates
    if covariates_matrix is not None:
        if to_pca == True:
            raise ValueError("Both covariates matrix and pca_covariates are passed as arguments, only one method can be run at a time")
        print("Running with consideration to passed covariates.", flush=True)    
        covar_df = pd.read_csv(args.covariates, delim_whitespace=True, header=None)
        covariates_matrix = covar_df.iloc[:, 2:].values
    
    if to_pca == True:
        print("Running PCA on genotypes to determine top 4 components; alternative covariate method", flush=True)
        
        # read genotypes
        vcf = VCF(args.vcf)
        sample_ids = np.array(vcf.samples)
        n_samples = len(sample_ids)

        genotypes = []
        chroms = []
        positions = []

        for variant in vcf:
            # variant.genotypes gives [ [0,1,True], ... ]
            # first two entries are alleles
            gt = np.array([g[0] + g[1] if g[0] != -1 else np.nan 
                           for g in variant.genotypes])
            genotypes.append(gt)
            chroms.append(variant.CHROM)
            positions.append(variant.POS)

        # transpose to get individuals x snps
        G = np.array(genotypes).T  

        # Run PCA to get top 4 components
        pca = PCA(n_components=4)
        pcs = pca.fit_transform(G)  
        covariates_matrix = pcs
                       
    gwas_fxn_cov(vcf_file=args.vcf, phenotype_file=args.phenotype, out_name=args.out, covariates_matrix=covariates_matrix)
    
    elapsed = time.time() - start
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print(f"Wall time: {elapsed:.2f}s")
    print(f"Peak memory: {peak / 1024 / 1024:.2f} MB")