# Genomic data, the fastq format and basic quality control
Here, we will familiarize ourselves with genomic data, typically arriving in the 'fastq' format. We will run some analysis which allow us to get an idea of the quality of the data and if we need to do any pre-processing before moving on to the next step, like mapping the data to a reference genome.
This section will use data from lions from [Bertola et al. (2022)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08510-y). And below is a young lion I met on my last trip to Kruger, as few years ago.

:camera: Lion in Kruger NP, South Africa
![lion2](./Images/.jpg)
©Laura Bertola

## Raw data
Imagine you're doing a project on lions and would like to include some previously published data in your analyses. For example, because it allows you to compare your population to other populations elsewhere in Africa. You may want to check [NCBI](https://www.ncbi.nlm.nih.gov/) or [ENA](https://www.ebi.ac.uk/ena/browser/home) if they have anything interesting online.
