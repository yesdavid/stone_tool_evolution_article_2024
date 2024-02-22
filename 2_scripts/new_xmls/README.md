# File guide and notes

1. ./FBD_constant_rate_BLANK.xml - constant FBD with fixed ages
2. ./FBD_constant_rate_age_uncertainty_BLANK.xml - constant FBD with estimated ages, except the youngest (= 0.638)
3. ./FBD_constant_rate_age_uncertainty_woffset_BLANK.xml  - constant FBD with estimated ages, inc. the youngest
4. ./FBD_skyline_BLANK.xml - variation in rates, fixed fossil ages
5. ./FBD_skyline_age_uncertainty_BLANK.xml - skyline FBD with estimated ages, except the youngest (= 0.638)

Unforunately it's not possible to combine trees with the offset parameter with the skyline model, so we have to fix the youngest sample age. 

Useful links:

* [Constant rate FBD SA xml template](https://github.com/CompEvol/sampled-ancestors/blob/master/examples/bears.xml)
* [Taming the BEAST FBD tutorial](https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf) (alternative parameterisation, fossil age uncertainty and offset)
* [Tree with offset xml template](https://bitbucket.org/bjoelle/fbd_study_code/src/master/brachiopods/empirical_input/beast_analyses/brachiopods_FBD_interval_ages.xml)
* [SA FBD skyline template](https://github.com/gavryushkina/sampled-ancestor-bdsky/blob/master/templates/SAFBDSKY.xml)
* [Taming the BEAST BDSky tutorial](https://taming-the-beast.org/tutorials/Skyline-analyses-for-macroevolution/)
* [Contraband xml template](https://github.com/fkmendes/contraband/blob/master/examples/testing/BMPruneLikelihood.xml) (multiple traits, BM model)
* [Contraband xml template](https://github.com/fkmendes/contraband/blob/master/examples/testing/BMMVNShiftLikelihoodOneTrait_FBDTree_RateCatClock.xml) (multiple clock rate categories, nCat = 2)

# Folder contents

The `*_BLANK.xml` files are used by the R script `../12_fill_BEAST2_scripts.R` to create the BEAST2 xml-scripts with the chosen data ready to be run. `run_BLANK.sh` and `resume_BLANK.sh` are the scripts to either start or resume these scripts as jobs on a HPC using Slurm. Depending on the user's HPC environment, the run and resume scripts have to be manually edited. 

The final xml-scripts used in the article can be found in the sub-folder `./BEAST2_contraband` for each taxa-trait combination. The sub-folder `./done` contains all MCC-trees resulting from the scripts in `./BEAST2_contraband`. The sub-folder `./underPrior` contains the scripts for running the models under the prior. 
