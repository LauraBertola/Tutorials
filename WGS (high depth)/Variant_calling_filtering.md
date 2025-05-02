## Variant calling
Now we all the genomes mapped, we'd like to filter out the SNPs that we'll be using for downstream analyses. There are different variant callers out there, notably [bcftools](https://samtools.github.io/bcftools/bcftools.html) and [strelka](https://github.com/Illumina/strelka). Here, we'll use bcftools. Run the following:
```
bcftools mpileup -Ou -f reference.fasta .bam > output.bcf
```

