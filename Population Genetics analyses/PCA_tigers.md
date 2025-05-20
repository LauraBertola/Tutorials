**PCA on the :tiger: data from the WGS (high depth) tutorial**  
To get an idea of how our filtering impacts actual downstream results, let's run a couple of PCAs. This is to illustrate that choices about filtering are extremely important, and not thinking about this carefully can lead to wrong conclusions about your study system!

Because the dataset we've been using so far only has 9 individuals, we'll use another, pre-filtered vcf with more populations. Download this datafile:
```
wget https://zenodo.org/records/15173226/files/machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz
```

Check how many individuals are in this file:
```
/softwares/bcftools1.12/bcftools query -l machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz | wc -l
```

And how many SNPs:
```
/softwares/bcftools1.12/bcftools view -H machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz | wc -l
```

Maybe we should plot the amount of missing data here as well:
```
/softwares/bcftools1.12/bcftools query -l machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz > sample_names_machali.txt
```

```
paste sample_names_machali.txt <( \
  /softwares/bcftools1.12/bcftools query -f '[%GT\t]\n' machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz | \
  awk '{
    for (i=1; i<=NF; i++) if ($i == "./.") missing[i]++;
  }
  END {
    for (i=1; i<=length(missing); i++) print missing[i];
  }'
) | sort -k2 -nr > missing_data_machali.txt
```

We'll plot this again in R. Do the following:
```
/softwares/R-4.2.3/bin/R
```

We're now inside of R. Run the following code here:
```
# Load libraries
library(ggplot2)

# Read the data
missing_data <- read.table("missing_data_machali.txt", header=FALSE, col.names=c("Sample", "Missing"))

# Plot
ggplot(missing_data, aes(x = reorder(Sample, -Missing), y = Missing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Missing Genotypes per Sample", x = "Sample", y = "Count of Missing Genotypes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

Quit R with (saving your workspace is not necessary):
```
quit()
```

These data are already filtered for missingness, so we can start exploring the results with further filtering. First we need to convert the vcf files to plink format, and then we can use plink to calculate the eigenvectors.
```
/softwares/plink/plink --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz \
  --make-bed --double-id --allow-extra-chr --out plink_file
```

This creates the following files, which are needed for a lot of downstream analyses:
| File                  | Description                                   | Contents Summary                              |
|-----------------------|-----------------------------------------------|-----------------------------------------------|
| `plink_file.bed`       | Binary genotype data file                      | Genotype matrix (SNPs × individuals) stored in compact binary format |
| `plink_file.bim`       | Variant information file (text)                | One SNP per line: chromosome, SNP ID, genetic distance, physical position, allele 1, allele 2 |
| `plink_file.fam`       | Sample information file (text)                 | One individual per line: family ID, individual ID, paternal/maternal ID, sex, phenotype |
| `plink_file.nosex`     | Optional sample IDs file without sex info     | List of sample IDs when sex info is missing or irrelevant           |

Now we run the PCA:
```
/softwares/plink/plink --bfile plink_file --pca 10 --allow-extra-chr --out plink_pca
```

We are going to plot this result in R:
```
/softwares/R-4.2.3/bin/R
```
```
library(ggplot2)      

eigenvec_data <- read.table("plink_pca.eigenvec", header=FALSE)
colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep=""))
head(eigenvec_data)
ggplot(eigenvec_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot Missing01: PC1 vs PC2") +
  theme_minimal()
```
```
quit()
```

