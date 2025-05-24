## Branching

It is likely that you want to explore your data with a number of different filters, to get a better understanding. In the tutorial we used `min_samples_locus` = 4, which means that only SNPs which are succesfully called in 4 or more samples are retained. Let's try some other values, using branching.

```
ipyrad -p params-cheetah.txt -b minsamples10
```
```
  loading Assembly: cheetah
  from saved path: ~/ipyrad-workshop/cheetah.json
  creating a new branch called 'minsamples10' with 24 Samples
  writing new params file to params-minsamples10.txt
```

This creates a new params file (as it says) which you should edit and modify to update the following parameter:

```
10              ## [21] [min_samples_locus]: Min # samples per locus for output
```

Now you can run step 7 again to generate the new output files with this new
`min_samples_locus` setting:

```
ipyrad -p params-minsamples10.txt -s 7 -c 4
```

This will create a new set of output files in `minsamples10_outfiles` which have only retained loci present in 10 or more samples. Look at the stats file
to see how many loci are retained in this dataset. You can try a number of different thresholds.

Another case where branching comes in handy is when you decide to remove some individuals, for example because they have a lot of missing data. Imagine, we'd like to remove the following samples:
- SRR19760947
- SRR19760949
- SRR19760954

You can do this with the following command:
```
ipyrad -p params-cheetah.txt -b 3inds_removed - SRR19760947 SRR19760949 SRR19760954
```
```
  loading Assembly: cheetah
  from saved path: ~/ipyrad-workshop/cheetah.json
  creating a new branch called '3inds_removed' with 21 Samples
  writing new params file to params-3inds_removed.txt
```

Note that here, we branched off the original assembly, with `min_samples_locus` at 4. Also note that with the latter example the samples are immediately removed and you don't have to edit the params file afterwards.

Proceed to rerunning step 7, and see how this affects the number of loci retained.
```
ipyrad -p params-3inds_removed.txt -s 7 -c 4
```
