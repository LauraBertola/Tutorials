## Variant calling
Now we have all the genomes mapped, with reads sorted and deduplicated, and all .bam files indexed. We can now proceed to call the variant positions that we'll be using for downstream analyses. There are different variant callers out there, but here we'll use [bcftools] (https://samtools.github.io/bcftools/bcftools.html) and [strelka](https://github.com/Illumina/strelka). Run the following:
```
/softwares/bcftools1.12/bcftools mpileup -Ou -f input_files/reference/GCA_021130815.1_PanTigT.MC.v3_genomic.fna output_files/*deduplicated.bam | \
/softwares/bcftools1.12/bcftools call -f GQ,GP -mv -Ov -o unfiltered_variants.vcf
```

The command has quite a few layers, so let's unpack those here. The `mpileup` command basically collects and summarized the data from all the .bam files, using the reference (indicated by the -f flag). The flag -Ou specifies that the output file should be uncompressed. This file then gets piped (|) into the next command, which is `call`. You ask it to use a multiallelic caller and output variant sites only (-mv). You can also ask for more information fields by using the -f flag with the `call` command (note that in this context, the -f flag means something else than the -f flag with the `mpileup` command!). Further, you specify that the output file should be an uncompressed vcf file (-Ov) and you specify the where the output file should be saved (-o).

You can also output a bcf file. The difference between vcf and bcf is similar to sam and bam, as we saw previously. vcf is human-readable, but slowerw to process. bcf is binary and not human-readable, but therefore smaller in file size and faster to process in pipelines. Either way, this step will take some time... :hourglass:

A vcf file contains a lot of information, and there are many ways of adjusting what the output should look like. More information about the format of vcf files, as well as additional flags to use, can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Our vcf file will contain the following information fields:
| FORMAT Tag | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `GT`       | Genotype — encoded like `0/0`, `0/1`, `1/1`, etc.                           |
| `GQ`       | Genotype Quality — Phred-scaled confidence in the genotype call            |
| `DP`       | Read Depth — total number of reads covering the site for the sample        |
| `PL`       | Phred-scaled genotype Likelihoods — likelihoods for `0/0`, `0/1`, and `1/1` |

In the fields above, 0 stands for the reference allele, whereas 1 stands for the derived allele. Let's take a closer look at our vcf file. First, we'll start with the header, by using -h:
```
/softwares/bcftools1.12/bcftools view -h unfiltered_variants.vcf
```

Lots of info here, but the last line is quite important, because it tells you what the format is of the vcf file:  
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT 

Now, let's take a look at the actual data, instead of the header, and use `less` -S so it doesn't print all the variants on the screen:
```
/softwares/bcftools1.12/bcftools view -H unfiltered_variants.vcf | less -S
```

You should see something like this:
![vcf](Images/vcf.png)

You can scroll through it, see what the reference allele and the alternative alleles are for each of the positions, and what information is available per sample. There is a lot of information in the INFO column:  
| INFO Tag   | Description                                                                                   |
|------------|--------------------------------------------------------------------------------------------|
| `DP`       | Total read depth at this position (from all samples)                                      |
| `SGB`      | Segregation-based score used internally by `bcftools` to penalize low-quality calls       |
| `MQ0F`     | Fraction of reads with mapping quality zero — `0` means all reads were well-mapped        |
| `AC`       | Allele count for the ALT allele                                                           |
| `AN`       | Allele number — total number of called alleles (e.g., diploid × number of samples)        |
| `DP4`      | Strand-specific read depth: REF-forward, REF-reverse, ALT-forward, ALT-reverse            |
| `MQ`       | Average mapping quality of reads covering this site                                       |

Finally, you may be curious how many variants are in your vcf file. Do:
```
/softwares/bcftools1.12/bcftools view -v snps unfiltered_variants.vcf | grep -vc "^#"
```

## Variant filtering
Our vcf file has a lot of variants, but we don't want to use everything. For example, there will be a lot of noise in the file, e.g. positions which do not have reliable information, or having missing data for most samples. We therefore need to filter our vcf file.


region: /softwares/bcftools1.12/bcftools view -H -r chr2:100000-200000 variants.vcf
genotypes: /softwares/bcftools1.12/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' variants.vcf | head
basecalls: /softwares/bcftools1.12/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' variants.vcf | head
