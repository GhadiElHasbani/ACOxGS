# ACOxGS

## An Optimization of Gene-Level Subnetwork Identification Inspired by Ant Colony and GeneSurrounder

### Abstract
High-throughput experimental technologies can provide deeper insights into pathway perturbations in biomedical studies. Accordingly, their usage is central to the identification of molecular targets and the development of suitable treatments for various diseases. Classical interpretations of generated data, such as differential gene expression and pathway analyses, disregard interconnections between studied genes when looking for gene-disease associations. Given that these interconnections are central to cellular processes, there has been a recent interest in incorporating them in such studies. The latter allows the detection of gene modules that underlie complex phenotypes in gene interaction networks. Existing methods either impose radius-based restrictions or grow modules freely at the expense of a statistical bias towards large modules. We propose a heuristic method, inspired by Ant Colony Optimization, to apply gene-level scoring and module identification with distance-based search constraints and penalties rather than radius-based restrictions. We test and compare our results to other approaches using three different neurodegenerative diseases, namely Alzheimer’s, Parkinson’s, and Huntington’s, over three independent experiments. We report the outcomes of enrichment analyses and concordance of gene-level scores for each disease. Results indicate that the proposed approach generally shows superior stability in comparison to existing methods. It produces stable and meaningful enrichment results in all three datasets which have different case to control proportions, sample sizes, and preprocessing steps. 

### Guide
1. Download any missing dependecies using `packages_download.R`
2. Perform module identification with ACOxGS, GS, and/or LEAN using your dataset of choice through sample code in `testrun.R` (includes KEGG network construction steps) and `networkX.ipynb` (for largest connected component extraction)
3. Perform enrichment and/or cross-study concordance analyses using `enrichment.R` and `concordance.R`, respectively
