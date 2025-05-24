# Empirical data & Quality Control (QC)

We will be reanalysing RAD-Seq data from cheetahs (*Acinonyx jubatus*) sampled from across their distribution in Africa and Iran and published in [Prost *et al.* 2022](https://onlinelibrary.wiley.com/doi/10.1111/mec.16577). The study used various datatypes, including whole genome sequencing (WGS), mitochondrial data, MHC data, minisatellites and RADseq data. For this workshop, we will focus only on (part of ) the RADseq data, which consists of 23 individuals from 6 populations, and one outgroup (puma; *Puma concolor*). The data were generated using a double-digest restriction-site associated DNA (ddRAD) sequencing approach [Peterson *et al.*, 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135). Note that raw reads have been randomly downsampled to 125,000 reads per sample, in order to create a dataset that will be computationally tractable with the expectation of finishing in a reasonable time. 

To keep things organized, first make a folder in your directory for this Tutorial, using `mkdir`, check if it is there with `ls`, and then go to you're newly created folder, using `cd`.
```
mkdir Tutorial_RADseq
```
```
ls
```
```
cd Tutorial_RADseq
```

The downsampled cheetah data are in my input_files folder. So, now create a symbolic link to my folder by doing the following:
```
ln -s /home/uramakri/laurabertola/Tutorial_RADseq/input_files input_files
```

You should now have something which *looks* like a folder, called input_files, in your directory, but actually it teleports you to *my* folder when you enter it.

If you do `ls` now, you'll see the files of the individual samples:
```
(ipyrad) osboxes@osboxes:~$ ls
subset-R1-raws subset-R1-raws.tgz
```
```
(ipyrad) osboxes@osboxes:~$ ls subset-R1-raws/
```
![png](images/ls_raws.png)

> **Special Note:** You'll see that every file contains `_R1_`. Most of the time, the data you will be working on are paired-end, meaning that each sample has a `_R1_` and `_R2_` file. For this workshop, and to ensure that the steps run quickly, we will only use `_R1_`. 

## The `FASTQ` data file format

Before we get started with the assembly, let's take a look at what the raw data
looks like. We can use `zcat` and `head` to do this.

Here we have our first look at a fastq formatted file. Each sequenced read is spread over four lines, one of which contains sequence and another the quality scores stored as ASCII characters. The other two lines are used as headers to store information about the read.

```bash
## zcat: unZip and conCATenate the file to the screen
## head -n 20: Just take the first 20 lines of input

(ipyrad) osboxes@osboxes:~/ipyrad-workshop$ zcat subset-R1-raws/SRR19760910_R1_.fastq.gz | head -n 20
```
```
@SRR19760910.1 1 length=110
CATGCACGTGCAGCATATAAGAAGGATGTTTGTCATGCATTATCTTATTTGATGTTTACGGAAGCCCCATGGTTATCCCCATTTTAGGGATGAAGAAACGCCACAGAGAT
+SRR19760910.1 1 length=110
BFFFFFF<FBFFFFFF<FB//////<<FFFBFF//<FFFFFFBF/FBFFFFFFFFFFFFBB<F/BFFFFFFFFBFF/<<</BFBBFF/<FF<FF<7FFFF/7B/FF/B<7
@SRR19760910.2 2 length=110
CATGCAACTCTTGGTCTCGGGGTCTTGAGTTCGAGCCCCACGTTGGATTAGAGATTACTTAAATAAATAAAGTTCAAAAGTTTTAGAATGTTATCATTTTCTTTAACAGT
+SRR19760910.2 2 length=110
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR19760910.3 3 length=110
CATGCCATTTCCCATGGGCAAGGATCTCAGGCTGTGCTCATTCCCAAGGACAAGACCAAGCCAATTCCCAATCCCCATATTTAAGGAGCTGCTTCCTGGGACCAATTCTG
+SRR19760910.3 3 length=110
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFF
@SRR19760910.4 4 length=110
CATGCAACTCTTGATCTCAGGGTCATGAGTTCAAGCCCCACATTGGGTATGGATCCTACTGAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGACAAG
+SRR19760910.4 4 length=110
FFBFFFFFFFFFFFFFBFFF<FFBFFBBFBFBF/B/BBFFBF<FFFB</BFBFB<BBFFFFFFBFFFBFFFFFF<FFF/FFFBFFF</FFBFFFFFFBFFFFFFFFFFFF
@SRR19760910.5 5 length=110
CATGCATTTGTGTTTGCTTCTATTTGTATGAAGAGTCGAGAAACCAGAAGCTAATACAAAGGGTTGCCCTTGGTAGGGGATGCTGACTGGATGGCTTTGGGGCAGGAGGA
+SRR19760910.5 5 length=110
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFBFFFFFFFFFFBFFFFFFFFFFF<FFFFFFFFFBFFFFFFFFFFFB/7
```

The first is the name of the read (its location on the plate). The second line
contains the sequence data. The third line is unused. And the fourth line is
the quality scores for the base calls. The [FASTQ wikipedia](https://en.wikipedia.org/wiki/FASTQ_format)
page has a good figure depicting the logic behind how quality scores are encoded.

## FastQC for quality control
The first step of any RAD-Seq assembly is to inspect your raw data to estimate overall quality. At this stage you can then attempt to improve your dataset by identifying and removing samples with failed sequencing. Another key QC procedure involves inspecting average quality scores per base position and trimming read edges, which is where low quality base-calls tend to accumulate. In this figure, the X-axis shows the position on the read in base-pairs and the Y-axis depicts information about [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score) per base for all reads, including median (center red line), IQR (yellow box), and 10%-90% (whiskers). As an example, here is a very clean base sequence quality report for a 75bp RAD-Seq library. These reads have generally high quality across their entire length, with only a slight (barely worth mentioning) dip toward the end of the reads:

![png](images/fastqc-high-quality-example.png)

In contrast, here is a somewhat typical base sequence quality report for R1 of a 300bp paired-end Illumina run of another RADseq dataset:

![png](images/fastqc-quality-example.png)

This figure depicts a common artifact of current Illumina chemistry, whereby quality scores per base drop off precipitously toward the ends of reads, with the effect being magnified for read lengths >150bp. The purpose of using FastQC to examine reads is to determine whether and how much to trim our reads to reduce sequencing error interfering with basecalling. In the above figure, as in most real dataset, we can see there is a tradeoff between throwing out data to increase overall quality by trimming for shorter length, and retaining data to increase value obtained from sequencing with the result of increasing noise toward the ends of reads.

### Running FastQC on the cheetah data
In preparation for running FastQC on our raw data we need to make an output directory to keep the FastQC results organized:

```
(ipyrad) osboxes@osboxes:~/ipyrad-workshop$ mkdir fastqc-results
```
Now run fastqc on one of the samples:
```
(ipyrad) osboxes@osboxes:~/ipyrad-workshop$ fastqc -o fastqc-results subset-R1-raws/SRR19760910_R1_.fastq.gz
```
> **Note:** The `-o` flag tells fastqc where to write output files. **Especially Notice** the *relative path* to the raw file. The difference between *relative* and *absolute* paths is an important one to learn. Relative paths are specified with respect to the current working directory. Since I am in `/home/ipyrad-workshop`, and this is the directory the `subset-R1-raws` directory is in, I can simply reference it directly. If I was in any other directory I could specify the *absolute path* to the target fastq.gz file which would be `/home/ipyrad-workshop/subset-R1-raws/SRR19760910_R1_.fastq.gz`. Absolute paths are always more precise, but also always (often _much_) longer.

FastQC will indicate its progress in the terminal. This toy data will run quite quickly, but real data can take somewhat longer to analyse (10s of minutes).

![png](images/fastqc-run.png)

If you feel so inclined you can QC all the raw data using a wildcard substitution:
```
(ipyrad) osboxes@osboxes:~/ipyrad-workshop$ fastqc -o fastqc-results subset-R1-raws/*
```
> **Note:** The `*` here is a special command line character that means "Everything that matches this pattern". So here `subset-R1-raws/*` matches _everything_ in the raws directory. Equivalent (though more verbose) statements are: `ls subset-R1-raws/*.gz`, `ls subset-R1-raws/*.fastq.gz`, `ls subset-R1-raws/*_R1_.fastq.gz`. All of these will list all the files in the `raws` directory. **Special Challenge:** Can you construct an `ls` command using wildcards that only lists samples in the `subset-R1-raws` directory that include the digit 5 in their sample name?

Examining the output directory you'll see something like this (assuming that you ran FastQC on all files):

![png](images/fastqc-folder.png)

Now we have output files that include html and images depicting lots of information about the quality of our reads, but we can't inspect these from inside our little black window. But we can go back to the Notebook server, and navigate to the results by first clicking on `ipyrad-workshop` and then `fastqc-results`. Open one of the html files by double clicking on them.

![png](images/notebook-fastqc.png)

### Instpecting and Interpreting FastQC Output

Just taking a random one, lets spend a moment looking at the results from `SRR19760910_R1__fastqc.html`. Opening up this html file, on the left you'll see a summary of all the results, which highlights areas FastQC indicates may be worth further examination. We will only look at a few of these.

![png](images/fastqc-summary.png)

Lets start with Per base sequence quality, because it's very easy to interpret, and often times with RAD-Seq data results here will be of special importance.

![png](images/fastqc-perbasequal.png)

For the cheetah data the sequence quality per base is uniformly quite high, with minor dips only in the first and last few bases (again, this is typical for Illumina reads). Based on information from this plot we can see that the cheetah data doesn't need trimming, which is good.

Now lets look at the `Per base sequece content`, which FastQC highlights with a scary red **X**.
![png](images/fastqc-perbasecontent.png)

The squiggles indicate base composition per base position averaged across the reads. It looks like the signal FastQC is concerned about here is related to the *extreme* base composition bias of the first 5 positions. We happen to know this is a result of the restriction enzyme overhang present in all reads (`CATGC` in this case for the SphI enzyme used), and so it is in fact of no concern. 

All in all, the data look good and we can proceed with the ipyrad analysis.

# References
Prost, Stefan, Ana Paula Machado, Julia Zumbroich, Lisa Preier, Sarita Mahtani-Williams, Rene Meissner, Katerina Guschanski, et al. 2022. Genomic Analyses Show Extremely Perilous Conservation Status of African and Asiatic Cheetahs (Acinonyx Jubatus). Molecular Ecology 31 (16): 4208–23.

Peterson, Brant K., Jesse N. Weber, Emily H. Kay, Heidi S. Fisher, and Hopi E. Hoekstra. 2012. Double Digest RADseq: An Inexpensive Method for de Novo SNP Discovery and Genotyping in Model and Non-Model Species. PloS One 7 (5): e37135.

