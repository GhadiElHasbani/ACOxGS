# ACOxGS

## Ant colony optimization for the identification of dysregulated gene-level subnetworks from expression data

### Abstract
**Background:** High-throughput experimental technologies can provide deeper insights into pathway perturbations in biomedical studies. Accordingly, their usage is central to the identification of molecular targets and the subsequent development of suitable treatments for various diseases. Classical interpretations of generated data, such as differential gene expression and pathway analyses, disregard interconnections between studied genes when looking for gene-disease associations. Given that these interconnections are central to cellular processes, there has been a recent interest in incorporating them in such studies. The latter allows the detection of gene modules that underlie complex phenotypes in gene interaction networks. Existing methods either impose radius-based restrictions or freely grow modules at the expense of a statistical bias towards large modules. We propose a heuristic method, inspired by Ant Colony Optimization, to apply gene-level scoring and module identification with distance-based search constraints and penalties, rather than radius-based constraints.

**Results:** We test and compare our results to other approaches using three datasets of different neurodegenerative diseases, namely Alzheimer’s, Parkinson’s, and Huntington’s, over three independent experiments. We report the outcomes of enrichment analyses and concordance of gene-level scores for each disease. Results indicate that the proposed approach generally shows superior stability in comparison to existing methods. It produces stable and meaningful enrichment results in all three datasets which have different case to control proportions and sample sizes.

**Conclusion:** The presented network-based gene expression analysis approach
successfully identifies dysregulated gene modules associated with a certain disease. Using a heuristic based on Ant Colony Optimization, we perform a distance-based search with no radius constraints. Experimental results support the effectiveness and stability of our method in prioritizing modules of high relevance. Our tool is publicly available at github.com/GhadiElHasbani/ACOxGS.git.

**Keywords:** gene expression, enrichment analysis, gene interaction network, ant colony
optimization

### Guide
1. Download any missing dependecies using `packages_download.R`
2. Perform module identification with ACOxGS, GS, and/or LEAN using your dataset of choice through sample code in `testrun.R` (includes KEGG network construction steps) and `networkX.ipynb` (for largest connected component extraction)
3. Perform enrichment and/or cross-study concordance analyses using `enrichment.R` and `concordance.R`, respectively
