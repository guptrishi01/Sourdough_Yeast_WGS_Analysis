if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("ggtree")
install.packages("ape")
install.packages("ggplot2")
install.packages("ggtree")
install.packages("crayon")

library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# STEP 1: Convert VCF to GDS
vcf.fn <- "all_chromosomes_combined.vcf"
sample_ancestors <- "sample_ancestors.csv"
gds.fn <- "all_chromosomes_combined.gds"
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")

# STEP 2: Open GDS file
genofile <- snpgdsOpen(gds.fn)

# STEP 3: Compute IBS distance matrix
ibs <- snpgdsIBS(genofile, num.thread=4)
ibs.dist <- 1 - ibs$ibs

# STEP 4: Build Neighbor-Joining tree
nj.tree <- nj(as.dist(ibs.dist))
# Set sample IDs to tip
nj.tree$tip.label <- ibs$sample.id
# Root tree based on sample ID
rooted_tree <- root(nj.tree, outgroup = "ERR1308867")

# Metadata for label
samplecategories = read.csv(sample_ancestors,
                            sep = ",",
                            col.names = c("Sample Name","Ancestor"),
                            header = FALSE,
                            stringsAsFactors = FALSE
                            )
dd = as.data.frame(samplecategories)

# STEP 6: Plot with ggtree in radial layout
p <- ggtree(rooted_tree, layout="circular", branch.length="none") +
  geom_tiplab(size=1.5, align=FALSE) +
  theme_tree2() +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p %<+% dd + 
  geom_tiplab(aes(fill = factor(Ancestor)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) +
  theme(legend.position = c(-0.8, 0), 
        legend.title = element_blank(), # no title
        legend.key = element_blank())

ggsave()
ggsave("neighbor_tree.png", p)