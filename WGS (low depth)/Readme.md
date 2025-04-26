# WGS (low depth)
This part will deal with processing of whole genome sequencing data (WGS) of low sequencing depth. Here, I'm referring to 1-5X genome-wide average sequencing depth, which does not allow for calling heterozygous sites reliably using hard calls. Instead, we'll used the ANGSD framework to obtain genotype likelihoods, which leverages information from all samples by looking at overall allele frequencies.

This section will use data from lions from [Bertola et al. (2022)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08510-y). Note that the data are downsampled, so these steps can be run quickly, for educational purposes only.
