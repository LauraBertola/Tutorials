# Raw data

We'll provide raw :tiger: data for three samples to illustrate the steps of quality control, trimming and mapping. The raw data come in .fastq.gz (or .fq.gz) format. Let's take a look at one of these files first.

```
user@cluster:~$ zcat yourfile.fq.gz | less -S
```

The `zcat` command is used to unzip the file (note it ends with .gz), and `less` allows you to view it. Here we add `-S` to chop off long lines. Otherwise it wraps around and becomes messy. If you'd like to see that, try using `head -n 20` instead of the `less` command. `head` shows the first part of the file, and `-n 20` tells it to show the first 20 lines. 
