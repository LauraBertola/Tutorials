# Raw data

We'll provide raw :tiger: data for three samples to illustrate the steps of quality control, trimming and mapping. The raw data come in .fastq.gz (or .fq.gz) format. Let's take a look at the data first.

Go to the folder where the data are stored, using `cd`, and use `ls` to display the contents of the folder. Most of the time, you'll be working with paired-end data, meaning that each sample has two files. These are usually identified by _R1 and _R2, or _1 and _2. Those two files contain the forward and reverse reads, respectively. For more information about paired-end Illumina sequencing, watch [this video](https://www.youtube.com/watch?v=fCd6B5HRaZ8).

Now let's take a look at what's inside these files. Type the following command:
```
user@cluster:~$ zcat yourfile.fq.gz | less -S
```

The `zcat` command is used to unzip the file (note it ends with .gz), and `less` allows you to view it. Here we add `-S` to chop off long lines. Otherwise it wraps around and becomes messy. Try scrolling right and down, using the arrow keys. To quit `less` type q. 
If you'd like to see the messy format, with long lines wrapping, try using `head -n 20` instead of the `less` command. `head` shows the first part of the file, and `-n 20` tells it to show the first 20 lines. 

Now to the actual data. Each sequenced read is spread over four lines, one of which contains sequence and another the quality scores stored as ASCII characters. The other two lines are used as headers to store information about the read.
It'll look something like this:
![fastq.gz](Images/fastq.gz.png)

The first is the name of the read, with information about its location on the plate, or in this case the identified from NCBI, where the data were downloaded from. The second line contains the sequence data. The third line is unused (identified with +). And the fourth line is the quality scores for the base calls. The (FASTQ wikipedia page)[https://en.wikipedia.org/wiki/FASTQ_format] has a good figure depicting the logic behind how quality scores are encoded.



