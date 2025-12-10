## Exercise 1 (FastQC_Trimming)

You've downloaded two samples already. But to make the SNP calling a little bit more interesting we should have a few more. Softlinking to my folder should now work, so try this again:
```
ln -s /home/uramakri/laurabertola/Tutorial_WGS_HD/input_files input_files
```

If you look inside this folder, you'll see that I've downloaded files for seven more individuals. Follow the tutorial, running FastQC and Trimming.

By the end of this exercise you should have a trimmed version for each of the downloaded samples in your output_files folder, nine in total. You can quikcly check this, by looking how many trimming reports you have in your folder. Each file should have two trimming reports, one of the forward reads, one for the reverse reads. If you don't feel like counting lines on your screen, let the computer do it:
```
ls *trimming_report.txt | wc -l
```

`ls` will list all the files ending with trimming_report.txt, but instead of displaying this, we'll pipe the result to the next command. `wc` is short for 'word count', and the -l flag tells it to count the lines.

## Exercise 2 (FastQC_Trimming)

Next, we will map the samples to a reference genome. We will talk about it more during our next session, but to prepare for it, you can do the following two steps.

Make a folder called 'reference' within your input_files folder first, and navigate into is with `cd`. Then download the reference genome, using `wget`:
```
wget https://zenodo.org/records/14258052/files/GCA_021130815.1_PanTigT.MC.v3_genomic.fna
```

For a reference genome to be used efficiently, it needs to be indexed. But indexing takes some time (~2 hours in my case). Your internet connection may drop, or maybe you feel like closing your laptop and going home. Both of these things will interrupt the run, which is really annoying when you're running something which takes long (hours, days, weeks...). For that, you can use `screen`.
'screen' is a terminal multiplexing tool that allows you to create, manage, and switch between multiple shell sessions from a single terminal window. Let's try to run the indexing from within screen.

First, start a new screen session, which you call 'index':
```
screen -S index
```

You're now inside your screen session. Now let's start the indexing, by doing the following:
```
bwa index GCA_021130815.1_PanTigT.MC.v3_genomic.fna
```

You can exit your screen session by typing CTRL+A and D. Nothing will show on your screen, it will just go back to your previous terminal window. You can safely close your computer and go offline, the process will continue to run in the screen session on the server.

To see if you have screen session running, do:
```
screen -ls
```

It should show you that you have a screen session called 'index'. To attach to the screen session again, for example to check if the run is finished yet, do:
```
screen -r index
```

If it is finished, you don't need this screen session anymore, and you can kill it. Do CTRL+A and K. It will ask you "Really kill this window [y/n]", and you can type y. If you now do:
```
screen -ls
```

It should say: "No sockets found", which means you don't have active screen sessions anymore.

The `bwa index` command has created the following files:
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.bwt
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.pac
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.ann
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.amb
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.sa

Check if these files indeed are in your folder now, using `ls`.

## Exercise 3 (Mapping)

Next time, we'll used the mapped data in the .bam files to call variants. However, calling variants from two files and then doing some popgen is not very interesting. Therefore, if you have downloaded and trimmed more files (see Exercise 1), process them in the same way. You should end up with an indexed bam file for all nine samples. Note that mapping takes some time (~15 minutes per sample), so it is advisable to use `screen` in this exercise as well.

## Exercise 4 (Variant_Calling_Filtering)

In the last session, we created a couple of PCA plots from a dataset with various levels of missingness. Compare this to results if you were to make a PCA directly from raw variants, without filtering. Use the following files to compare the results:
```
variants_snps.vcf
variants_snps_qual30.vcf
variants_snps_qual30_maf05_gq30_alldp18-60.vcf 
```

This is to illustrate that choices about filtering are extremely important, and not thinking about this carefully can lead to wrong conclusions about your study system!

Maybe you're also curious to know what the **real** data look like. You can download a more complete, pre-filtered test datafile here:
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

Explore this dataset for missingness per sample, and run a PCA to get an idea of the real 🐯 population structure.
