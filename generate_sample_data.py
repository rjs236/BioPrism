import os
import numpy as np
import pandas as pd

np.random.seed(42)

sample_data_dir = os.path.join(os.path.dirname(__file__), "static", "sample_data")
os.makedirs(sample_data_dir, exist_ok=True)

# 20 samples: 10 Tumor, 10 Normal
tumor_samples = [f"Tumor_{i:02d}" for i in range(1, 11)]
normal_samples = [f"Normal_{i:02d}" for i in range(1, 11)]
all_samples = tumor_samples + normal_samples

# Create Microarray Metadata
meta_df = pd.DataFrame({
    'Sample': all_samples,
    'Group': ['Tumor'] * 10 + ['Normal'] * 10
})
meta_df.to_csv(os.path.join(sample_data_dir, "demo_microarray_metadata.csv"), index=False)

# Gene lists
upregulated_genes = [
    "TP53", "EGFR", "ERBB2", "MYC", "VEGFA", "CCND1", "CDK4", "MKI67", "PCNA",
    "E2F1", "BCL2", "HIF1A", "STAT3", "AKT1", "PIK3CA", "KRAS", "BRAF", "MET",
    "CD274", "CTLA4", "MMP9", "PLAU", "AURKA", "FOXM1", "CDK2"
]

downregulated_genes = [
    "PTEN", "BRCA1", "BRCA2", "RB1", "CDKN1A", "CDKN2A", "BAX", "CASP3",
    "CASP8", "ATM", "CHEK2", "CDH1", "APC", "SMAD4", "VHL", "FOXO3", "DUSP1",
    "CAV1", "GATA3", "RAD51"
]

housekeeping_genes = [
    "ACTB", "GAPDH", "B2M", "HPRT1", "RPLP0", "TBP", "PGK1", "PPIA", "GUSB",
    "SDHA", "YWHAZ", "ALDOA", "TUBB", "LDHA"
]

background_genes = [f"GENE_{i:03d}" for i in range(1, 31)]

all_genes = upregulated_genes + downregulated_genes + housekeeping_genes + background_genes

# Generate log2 normalized microarray expression (range 4 to 14)
expr_data = {}
for g in all_genes:
    if g in upregulated_genes:
        # Higher in tumor: mean 10.5 vs 7.2 (log2FC ~ +3.3)
        tumor_vals = np.random.normal(loc=10.5, scale=0.45, size=10)
        norm_vals = np.random.normal(loc=7.2, scale=0.40, size=10)
    elif g in downregulated_genes:
        # Lower in tumor: mean 6.8 vs 10.2 (log2FC ~ -3.4)
        tumor_vals = np.random.normal(loc=6.8, scale=0.42, size=10)
        norm_vals = np.random.normal(loc=10.2, scale=0.45, size=10)
    elif g in housekeeping_genes:
        # Stable in both: mean ~ 12.0
        tumor_vals = np.random.normal(loc=12.2, scale=0.15, size=10)
        norm_vals = np.random.normal(loc=12.2, scale=0.15, size=10)
    else:
        # Random non-significant background
        base = np.random.uniform(5.0, 9.0)
        tumor_vals = np.random.normal(loc=base, scale=0.5, size=10)
        norm_vals = np.random.normal(loc=base, scale=0.5, size=10)

    expr_data[g] = list(np.round(tumor_vals, 2)) + list(np.round(norm_vals, 2))

expr_df = pd.DataFrame(expr_data, index=all_samples).T
expr_df.index.name = 'Gene'
expr_df.to_csv(os.path.join(sample_data_dir, "demo_microarray_expression.csv"))

# Also generate RNA-seq Count matrix demo
# Non-negative integers
counts_data = {}
for g in all_genes:
    if g in upregulated_genes:
        tumor_vals = np.random.negative_binomial(n=20, p=0.015, size=10) + 1200
        norm_vals = np.random.negative_binomial(n=10, p=0.04, size=10) + 200
    elif g in downregulated_genes:
        tumor_vals = np.random.negative_binomial(n=10, p=0.04, size=10) + 150
        norm_vals = np.random.negative_binomial(n=20, p=0.015, size=10) + 1400
    elif g in housekeeping_genes:
        tumor_vals = np.random.negative_binomial(n=50, p=0.01, size=10) + 4000
        norm_vals = np.random.negative_binomial(n=50, p=0.01, size=10) + 4000
    else:
        tumor_vals = np.random.negative_binomial(n=15, p=0.03, size=10) + 300
        norm_vals = np.random.negative_binomial(n=15, p=0.03, size=10) + 300
        
    counts_data[g] = list(tumor_vals) + list(norm_vals)

counts_df = pd.DataFrame(counts_data, index=all_samples).T
counts_df.index.name = 'Gene'
counts_df.to_csv(os.path.join(sample_data_dir, "demo_rnaseq_counts.csv"))

rnaseq_meta_df = pd.DataFrame({
    'Sample': all_samples,
    'Group': ['Tumor'] * 10 + ['Normal'] * 10
})
rnaseq_meta_df.to_csv(os.path.join(sample_data_dir, "demo_rnaseq_metadata.csv"), index=False)

# Generate independent validation dataset (8 new samples: 4 Tumor, 4 Normal)
indep_samples = [f"Validation_Tumor_{i}" for i in range(1, 5)] + [f"Validation_Normal_{i}" for i in range(1, 5)]
indep_meta_df = pd.DataFrame({
    'Sample': indep_samples,
    'Group': ['Tumor'] * 4 + ['Normal'] * 4
})

indep_expr = {}
for g in all_genes:
    if g in upregulated_genes:
        t_v = np.random.normal(loc=10.3, scale=0.5, size=4)
        n_v = np.random.normal(loc=7.4, scale=0.45, size=4)
    elif g in downregulated_genes:
        t_v = np.random.normal(loc=7.0, scale=0.45, size=4)
        n_v = np.random.normal(loc=10.0, scale=0.5, size=4)
    else:
        base = np.random.uniform(5.0, 9.0)
        t_v = np.random.normal(loc=base, scale=0.5, size=4)
        n_v = np.random.normal(loc=base, scale=0.5, size=4)
    indep_expr[g] = list(np.round(t_v, 2)) + list(np.round(n_v, 2))

indep_df = pd.DataFrame(indep_expr, index=indep_samples)
indep_combined = pd.concat([indep_meta_df.set_index('Sample'), indep_df], axis=1).reset_index()
indep_combined.to_csv(os.path.join(sample_data_dir, "independent_validation_data.csv"), index=False)

print("Generated all sample datasets successfully!")
