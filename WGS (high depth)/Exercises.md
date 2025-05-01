## Exercise 1

You've downloaded two samples already. But to make the SNP calling a little bit more interesting download the following samples as well, and process them in the same way as we did with the first two.
```
wget https://zenodo.org/records/14258052/files/BEN_CI18_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_CI18_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW10_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW10_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW12_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW12_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW13_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW13_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI18_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI18_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI9_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI9_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI19_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI19_sub_2.fq.gz

```

By the end of this exercise you should have a trimmed version for each of the downloaded samples in your output_files folder.

## Exercise 2

Next, we will map the samples to a reference genome. We will talk about it more during our next session, but to prepare for it, please do the following two steps.

First download the reference genome, the same way you downloaded the other data. Make a folder called 'reference' within your input_files folder first, and navigate into is with `cd`.
Then do:
```
wget https://zenodo.org/records/14258052/files/GCA_021130815.1_PanTigT.MC.v3_genomic.fna
```

For a reference genome to be used efficiently, it needs to be indexed. But indexing takes some time (~2 hours in my case). You're internet connection may drop, or maybe you feel like closing your laptop and going home. Both of these things will interrupt the run, which is really annoying when you're running something which takes long (hours, days, weeks...). For that, you can use `screen`.
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
screen ls
```

It should show you that you have a screen session called 'index'. To attach to the screen session again, for example to check if it is finished yet, do:
```
screen -r index
```

If it is finished, you don't need this screen session anymore, and you can kill it. Do CTRL+A and K. It will ask you "Really kill this window [y/n]", and you can type y. If you now do:
```
screen ls
```

It should say: "No sockets found", which means you don't have active screen sessions anymore.

The `bwa index` command has created the following files:
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.bwt
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.pac
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.ann
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.amb
- GCA_021130815.1_PanTigT.MC.v3_genomic.fna.sa

Check if these files indeed are in your folder now, using `ls`.
