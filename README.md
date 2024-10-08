# Research Compendium for: "A macroevolutionary analysis of European Late Upper Palaeolithic stone tool shape using a Bayesian phylodynamic framework"

David N. Matzig¹,✉ ORCID: [0000-0001-7349-5401](http://orcid.org/0000-0001-7349-5401)  
Ben Marwick² ORCID: [0000-0001-7879-4531](http://orcid.org/0000-0001-7879-4531)  
Felix Riede¹ ORCID: [0000-0002-4879-7157](http://orcid.org/0000-0002-4879-7157)    
Rachel C.M. Warnock³ ORCID: [0000-0002-9151-4642](http://orcid.org/0000-0002-9151-4642)    

¹ Department of Archaeology and Heritage Studies, Aarhus University, Denmark  
² Department of Anthropology, University of Washington, USA   
³ GeoZentrum Nordbayern, Friedrich-Alexander-University Erlangen, Germany

✉ Correspondence: David N. Matzig <david.matzig@cas.au.dk>  


### Compendium DOI:

[![DOI](https://zenodo.org/badge/511112726.svg)](https://zenodo.org/doi/10.5281/zenodo.10693325)

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/yesdavid/stone_tool_evolution_article_2024> are the development versions and may have changed since the paper was published.

### Maintainer of this repository:

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7349--5401-green.svg)](http://orcid.org/0000-0001-7349-5401) David N. Matzig (<david.matzig@cas.au.dk>) 

### Published in:
Matzig, D.N., Marwick, B., Riede, F., & Warnock, R.C.M. (2024). A macroevolutionary analysis of European Late Upper Palaeolithic stone tool shape using a Bayesian phylodynamic framework. R. Soc. Open Sci. 11: 240321. [https://doi.org/10.1098/rsos.240321](https://doi.org/10.1098/rsos.240321).  

### Abstract:
Phylogenetic models are commonly used in paleobiology to study the patterns and processes of organismal evolution. In the human sciences, phylogenetic methods have been deployed for reconstructing ancestor-descendant relationships using linguistic and material culture data. Within evolutionary archaeology specifically, phylogenetic analyses based on maximum parsimony and discrete traits dominate, which sets limitations for the downstream role cultural phylogenies, once derived, can play in more elaborate analytical pipelines. Moreover, the use of discrete character traits in these efforts prevails, which in turn sets a number of non-trivial challenges. Recent methodological advances in computational paleobiology, however, now allow us to infer Bayesian phylogenies using continuous characters. Capitalizing on these developments, we here present an exploratory analysis of cultural macroevolution of projectile point shape evolution in the European Final Palaeolithic and earliest Mesolithic (~15,000-11,000 BP) using a time-scaled Bayesian phylogeny and a fossilised birth-death sampling process model. This model-based approach leaps far beyond the application of parsimony, in that it not only produces a tree, but also divergence times, and diversification rates, which we compare to the pronounced climatic changes that occurred during this timeframe. While common in cultural evolutionary analyses of language, the extension of Bayesian phylodynamic models to archaeology represents a major methodological breakthrough.

### Keywords:
Cultural macroevolution; Bayesian phylogenies; phylogenetic comparative methods; geometric morphometrics; archaeology; Late Upper Palaeolithic; stone tools

### Required software packages and their versions:
The Bayesian phylogenies were inferred using `BEAST v2.6.6`, and the following packages: `ORC v1.0.3`, `BDMM-Prime v0.0.32`, `BDSKY v1.4.8`, `MODEL_SELECTION v1.5.3`, `BEASTLabs v1.9.7`, `SA v2.0.2`, `FastRelaxedClockLogNormal v1.1.1`, and `contraband v0.0.1` (which is provided in the `./2_scripts/contraband.jar` file).

All further analyses were conducted in `R v4.3.2` using the following packages: `beastio (>= 0.3.3)`, `coda (>= 0.19-4)`, `cowplot (>= 1.1.1)`, `data.table (>= 1.14.8)`, `dplyr (>= 1.1.2)`, `forcats (>= 1.0.0)`, `ggplot2 (>= 3.4.3)`, `ggpubr (>= 0.6.0)`, `ggrepel (>= 0.9.3)`, `ggthemes (>= 4.2.4)`, `ggtree (>= 3.6.2)`, `magrittr (>= 2.0.3)`, `Momocs (>= 1.4.0)`, `outlineR (>= 0.1.0)`, `raster (>= 3.6-20)`, `rcarbon (>= 1.5.0)`, `readr (>= 2.1.4)`, `RevGadgets (>= 1.1.0)`, `rgeos (>= 0.6-2)`, `rworldmap (>= 1.3-6)`, `sp (>= 1.6-0)`, `splitstackshape (>= 1.4.8)`, `tibble (>= 3.2.1)`, `tidyr (>= 1.3.0)`, `treeio (>= 1.22.0)`.

These particular package versions can be downloaded and installed from the Posit Package Manager via the `install_packages.R` script provided in this repository. Alternatively, users can download and employ the `StoneToolEvol2024.sif` Singularity/Apptainer container provided on [10.5281/zenodo.11186033](https://doi.org/10.5281/zenodo.11186033). To use the Singularity file interactively, use the following command: `singularity exec StoneToolEvol2024.sif R`.

### How to reproduce:

It is important to execute the scripts (found in `1_scripts`) in a certain order. First, scripts `10_C14_dates_v3.R`, `11_subset_outlines.R`, and `12_fill_BEAST2_scripts.R` have to be run in R. Script `12_fill_BEAST2_scripts.R` produces .xml files which have to be run with BEAST2 and the contraband package. BEAST2 with the above specified packages has to be installed on the computer/HPC, but the .xml files have to be run using the `contraband.jar` file (e.g. `java -jar contraband.jar 32_09_ULNC_FBD_skyline_age_uncertainty_underPrior.xml`). Afterwards, the results can be merged, analysed, and visualised in R using the R scripts `13_LOCAL_move_converged_runs.R`, `14_logfiles.R`, `15_plot_Skyline.R`, `15_plot_Skyline_ULNC.R`, `15_plot_Skyline_underPrior.R`, and `15_plot_Tree.R`. A short description of what each script does is given in `2_scripts/README.md`.

### Licenses:

Code: MIT <http://opensource.org/licenses/MIT> year: 2024, copyright holder: David Nicolas Matzig

