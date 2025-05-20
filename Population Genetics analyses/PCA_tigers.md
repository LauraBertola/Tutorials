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
/softwares/bcftools1.12/bcftools query -l machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz | \
sed 's|output_files/||; s/_aligned_reads_deduplicated.bam$//' > sample_names_cleaned.txt
```

```
paste sample_names_cleaned.txt <( \
  /softwares/bcftools1.12/bcftools query -f '[%GT\t]\n' machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz | \
  awk '{
    for (i=1; i<=NF; i++) if ($i == "./.") missing[i]++;
  }
  END {
    for (i=1; i<=length(missing); i++) print missing[i];
  }'
) | sort -k2 -nr > missing_data_machali.txt
```

We don't have a lot of samples here, but it's good to plot these results, so we have a visual if there are samples which are poorer in quality than others. We do this in R. This requires a slightly different syntax than bash, as we'll be working in another environment. Do the following:
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

Let's also do some filtering for different levels of missingness:
```
/softwares/bcftools1.12/bcftools +fill-tags -Oz machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz.vcf.gz -- -t F_MISSING  | \
/softwares/bcftools1.12/bcftools view -i 'INFO/F_MISSING<0.1' -Oz -o machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing01.vcf.gz
/softwares/bcftools1.12/bcftools view -H machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing01.vcf.gz | wc -l
```
```
/softwares/bcftools1.12/bcftools +fill-tags -Oz machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz.vcf.gz -- -t F_MISSING  | \
/softwares/bcftools1.12/bcftools view -i 'INFO/F_MISSING<0.25' -Oz -o machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing025.vcf.gz
/softwares/bcftools1.12/bcftools view -H machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing025.vcf.gz | wc -l
```
```
/softwares/bcftools1.12/bcftools +fill-tags -Oz machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz.vcf.gz -- -t F_MISSING  | \
/softwares/bcftools1.12/bcftools view -i 'INFO/F_MISSING<0.5' -Oz -o machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing05.vcf.gz
/softwares/bcftools1.12/bcftools view -H machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing05.vcf.gz | wc -l
```
```
/softwares/bcftools1.12/bcftools +fill-tags -Oz machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz.vcf.gz -- -t F_MISSING  | \
/softwares/bcftools1.12/bcftools view -i 'INFO/F_MISSING<0.75' -Oz -o machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing075.vcf.gz
/softwares/bcftools1.12/bcftools view -H machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3.recode.vcf.gz_missing075.vcf.gz | wc -l
```

Now, let's do some popgen, but seeing how a PCA for this dataset looks like. First we need to convert the vcf file to plink format, and then we can use plink to calculate the eigenvectors.
```
plink --vcf machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile.recode.vcf \
  --make-bed --double-id --allow-extra-chr --out plink_file
```
```
plink --bfile plink_file --pca 10 --allow-extra-chr --out plink_pca
```

We are going to plot this result in R:
```
/softwares/R-4.2.3/bin/R
```
```
library(ggplot2)      
eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE)
colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep=""))
head(eigenvec_data)
ggplot(eigenvec_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") +
  theme_minimal()
ggsave("pca_plot.png")
```
```
quit()
```

