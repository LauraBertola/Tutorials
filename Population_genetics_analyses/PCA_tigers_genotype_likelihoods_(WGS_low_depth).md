**PCA on the 🦁 data from the WGS low depth tutorial**  

## Population structure with PCAngsd

Apart from knowing how many SNPs we have, we'd actually also like to know what the data look like, for example when plotting a PCA. There is a special tool, [PCAngsd](https://www.popgen.dk/software/index.php/PCAngsd), which allows you to estimate the covariance matrix and individual allele frequencies for low-depth data, directly from the beagle file you created in the previous step.

First we need to activate the correct environment:
```
conda activate pcangsd
```

Then, run this:
```
pcangsd -b all_minind5.beagle.gz -o all_minind5 -t 8 --iter 1000
```

You tell it the input file (-b), the output files prefix (-o), how many threads to use for computation (-t) and a number of iterations to run (--iter). When the run is finished, and hopefully converged (it will tell you on the screen), you should have the following files: .cov and .log. Some versions of PCAngsd also produce .eigenvec and .eigenval files, but we can easily compute them ourselves from the covariance matrix (.cov), so we don't really need them.

Now, download the .cov file to your computer, as well as the all_bams.list file, so we can look at the results in R.

In R studio, run the following:
```
# Read in eigenvectors
cov <- as.matrix(read.table("all_minind5.cov"))
eig <- eigen(cov)
eigval <- eig$values
eigvec <- eig$vectors
write.table(eig$values, "all_minind5.eigenval.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(eig$vectors, "all_minind5.eigenvec.txt",
            quote = FALSE, row.names = FALSE, col.names = paste0("PC", 1:ncol(eig$vectors)))
variance_explained <- eig$values / sum(eig$values)
print(variance_explained[1:6])  # first 6 PCs

# Load sample names (from your bam list)
samples <- basename(readLines("all_bams.list"))

# Species PCA
# Assign species manually
species <- c(
  rep("Somalia", 1),
  rep("India", 1),
  rep("Zambia", 1),
  rep("Kenya", 1),
  rep("Zambia", 1),
  rep("DRC", 1),
  rep("Benin", 1),
  rep("South_Africa", 1),
  rep("Namibia", 1),
  rep("Cameroon", 1)
)

# Colors for species
species_colors <- c("Benin"="#9DB3F1",
                    "Cameroon"="#244932",   
                    "DRC"="#495C24", 
                    "Somalia"="#652F24", 
                    "Kenya"="#8f3a13", 
                    "Zambia"="#FF429A", 
                    "South_Africa"="#9848A8", 
                    "Namibia"="#C97FA2",
                    "India"="#FFC000"
                    )

# Make a data.frame
pca_df <- data.frame(
  Sample = samples,
  Species = species,
  PC1 = eigvec[,1],
  PC2 = eigvec[,2]
)

png("PCA.png", width = 1500, height = 1500, res = 300)

xpad <- diff(range(pca_df$PC1)) * 0.05
ypad <- diff(range(pca_df$PC2)) * 0.05

# Plot
plot(pca_df$PC1, pca_df$PC2, 
     col = species_colors[pca_df$Species],
     pch = 19, 
     cex = 1.5,
     xlab = paste0("PC1 (", round(eigval[1]/sum(eigval)*100, 1), "%)"),
     ylab = paste0("PC2 (", round(eigval[2]/sum(eigval)*100, 1), "%)"),
     xlim = range(pca_df$PC1) + c(-xpad, xpad),
     ylim = range(pca_df$PC2) + c(-ypad, ypad))

text(pca_df$PC1, pca_df$PC2, 
     labels = pca_df$Species, 
     pos = 4,  # position: above points
     cex = 0.7, 
     col = species_colors[as.character(pca_df$Species)])

dev.off()
```

This code will save two tables, one with eigenvectors and one with eigenvalues, calculated from the covariance matrix we got from PCAngsd. It also plots a PCA, with the population labels and colors defined in the code. It should look something like this:
<img src="./Images/pca.png" alt="pca" width="70%">
