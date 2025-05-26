# RADseq

This part will deal with processing of a Restriction site Associated DNA sequencing (RADseq), which is a form of reduced representation sequencing. This allows for cost-effective analyses of a larger number of samples, making this a useful strategy when funding is limited (and it always is :slightly_frowning_face:).

This section relies heavily on resources created for [RADcamp](https://radcamp.github.io/), so a big thank you to Dr. Isaac Overcast and all RADcamp contributors over the years! We will use data from cheetahs from [Prost et al. (2022)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16577). Note that the data are downsampled, so these steps can be run quickly, for educational purposes only.

You can find more information about the raw data and steps for quality control [here](Data_FastQC.md).  
The steps for using the actual Ipyrad pipeline are [here](Ipyrad.md).  
Note that there are exercises with additional information [here](Exercises.md).  

:camera: Cheetah coalition in the Masai Mara, Kenya
![cheetahs](./Images/DSC_3251.jpg)
©Laura Bertola

## Important note on why not to use a WGS pipeline on RADseq data

You technically can, but it's usually a bad idea unless you're very careful. Here's why:

1. WGS tools assume even, genome-wide coverage
Tools in standard WGS pipelines (e.g., GATK best practices):
* Expect reads to be distributed across the whole genome.
* Use population-level models for genotype likelihoods and filtering.
* Expect low missing data and higher coverage consistency.
  
→ RADseq violates all of these assumptions. Note that there is a difference between _randomly_ missing data (e.g. low depth WGS) and _systematically_ missing data (e.g. RADseq).

2. High missingness confuses variant callers
Variant callers like GATK HaplotypeCaller, FreeBayes, etc., assume:
* Every sample will have coverage at most variant sites.  
* Depth can be used as a proxy for genotype confidence.  

→ In RADseq, many loci are missing completely from many individuals, and coverage can be low and variable, which leads to:

  Poor genotype calling:
* Inflated missingness in VCF.
* Overfiltered or miscalled variants.

3. Duplicate handling & local realignment are overkill or misleading
WGS pipelines include steps like:
* Marking duplicates.
* Local realignment around indels.  

→ These aren't helpful (and may be harmful) for RADseq because:
- RADseq reads from the same fragment often look like duplicates but aren't PCR duplicates.
- Local realignment makes little sense with RAD loci that don't span structural variants.


**However,** You can use a reference genome to map RADseq reads (e.g., with bwa mem) if you want:
* Better locus ID consistency across individuals.
* Reference-anchored SNP positions.

But after mapping, you should still use a RADseq-aware pipeline (e.g., ipyrad, Stacks) to cluster loci and call variants.

## In summary
| **Feature**                      | **Low-Depth WGS Pipelines**                                                                 | **RADseq Pipelines**                                                           |
|----------------------------------|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| **Typical Tools**                | GATK, ANGSD, bcftools, FreeBayes                                                             | ipyrad, Stacks, pyRAD                                                           |
| **Data Structure**               | Reads spread across genome (sparse but random)                                               | Reads concentrated at restriction sites (locus-level)                          |
| **Missing Data Pattern**        | Randomly distributed                                                                          | Systematic (entire loci missing in individuals)                                |
| **Causes of Missingness**       | Low coverage, sequencing errors                                                              | Restriction site dropout, uneven digestion, low coverage                       |
| **Genotype Calling Approach**    | Genotype likelihoods (GLs), population priors (e.g., in ANGSD, GATK's HaplotypeCaller)       | Clustering-based consensus calling per locus                                   |
| **Handling of Missing Data**     | Statistical models infer genotypes/frequencies even with low depth                           | Filtering by minimum samples per locus (e.g., `min_samples_locus`)            |
| **Assumption About Coverage**    | Generally assumes uniform or randomly missing                                                | Explicitly expects high missingness and dropout                                |
| **Effect of Hard Filtering**     | Can be tuned by depth/quality thresholds                                                     | Over-stringent filtering removes most loci (bad)                               |
| **SNP Position Resolution**      | Precise, genome-wide                                                                         | Clustered within RAD loci (may lack genomic coordinates if de novo)            |


Of course, there are different arguments for using a specific pipeline, and consistency with previously generated datasets is one of them. But RADseq pipelines exist for a reason.
