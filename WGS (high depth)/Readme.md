# WGS (high depth)

This part will deal with processing of whole genome sequencing data (WGS) of medium-to-high sequencing depth. Here, I'm referring to >6X genome-wide average sequencing depth, which allows for reliably calling heterozygous sites following this pipeline.

This section will use data from Bengal tigers from [Khan et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8813985/). Note that the data are downsampled, so these steps can be run quickly, for educational purposes only. Here, we kind of _pretend_ that we have high depth data, while the downsampled data set actually has more resemblance to a low depth data set.

You can find the steps for quality control, trimming and mapping [here](FastQC_Trimming_Mapping.md).
Steps for variant calling and filtering are [here](Variant_Calling_Filtering.md).
Note that there are exercises with additional information [here](Exercises.md).

:camera: Tiger in Chitwan National Park, Nepal
![tiger](./Images/DSC_8741b.jpg)
©Laura Bertola
