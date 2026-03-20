import allel
import numpy as np
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt

# Load VCF
callset = allel.read_vcf('new_filtered_variants.vcf')
samples = callset['samples']
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# Convert genotypes to alt allele count
genotype_alt = genotypes.to_n_alt()
variant_filter = genotype_alt.std(axis=1) > 0
genotype_alt_filt = genotype_alt[variant_filter]
geno_matrix = genotype_alt_filt.T

# Impute missing values
imputer = SimpleImputer(strategy='mean')
geno_matrix = imputer.fit_transform(geno_matrix)

# PCA
pca = PCA(n_components=2)
coords = pca.fit_transform(geno_matrix)

# Load descendant groups
with open('sd_samples.txt') as f:
    sd_samples = set(line.strip() for line in f)

with open('nat_samples.txt') as f:
    nat_samples = set(line.strip() for line in f)

# Assign groups
group_labels = []
for s in samples:
    if s == '892B':
        group_labels.append('892B ancestor')
    elif s == '942':
        group_labels.append('942 ancestor')
    elif s == 'ERR1308867':
        group_labels.append('ERR1308867 ancestor')
    elif s in sd_samples:
        group_labels.append('892B descendants')
    elif s in nat_samples:
        group_labels.append('ERR1308867 descendants')
    else:
        group_labels.append('Other')

# Define styles: darker for ancestors, lighter for descendants
group_styles = {
    '892B ancestor': {'color': 'darkred', 'marker': 'o', 'size': 100},
    '942 ancestor': {'color': 'navy', 'marker': 's', 'size': 100},
    'ERR1308867 ancestor': {'color': 'darkgreen', 'marker': 'D', 'size': 100},
    '892B descendants': {'color': 'lightcoral', 'marker': 'o', 'size': 50},
    'ERR1308867 descendants': {'color': 'lightgreen', 'marker': 'D', 'size': 50},
    'Other': {'color': 'gray', 'marker': 'x', 'size': 40}
}

# Plot
plt.figure(figsize=(12, 10))
unique_groups = set(group_labels)

for group in unique_groups:
    idx = [i for i, g in enumerate(group_labels) if g == group]
    style = group_styles[group]

    # Scatter plot
    plt.scatter(coords[idx, 0], coords[idx, 1],
                label=group,
                c=style['color'],
                marker=style['marker'],
                s=style['size'],
                alpha=0.85)

    # Add text labels for each point
    for i in idx:
        plt.text(coords[i, 0] + 0.005, coords[i, 1] + 0.005, samples[i],
                 fontsize=8, color='black', alpha=0.9)


plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
plt.title("PCA of Genotypes by Ancestral Group")
plt.legend(title="Group")
plt.grid(True)
plt.tight_layout()
plt.show()
