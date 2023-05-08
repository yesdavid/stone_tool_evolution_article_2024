# File guide and notes

1. FBD_constant_rate.xml - constant FBD with fixed ages
2. FBD_constant_rate_age_uncertainty.xml - constant FBD with estimated ages, except the youngest (= 0.638)
3. FBD_constant_rate_age_uncertainty_woffset.xml  - constant FBD with estimated ages, inc. the youngest
4. FBD_skyline.xml - variation in rates, fixed fossil ages5. FBD_skyline_age_uncertainty.xml - skyline FBD with estimated ages, except the youngest (= 0.638=Unforunately it's not possible to combine trees with the offset parameter with the skyline model, so we have to fix the youngest sample age. 

Useful links:

* [Constant rate FBD SA xml template](https://github.com/CompEvol/sampled-ancestors/blob/master/examples/bears.xml)
* [Taming the BEAST FBD tutorial](https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf) (alternative parameterisation, fossil age uncertainty and offset)
* [Tree with offset xml template](https://bitbucket.org/bjoelle/fbd_study_code/src/master/brachiopods/empirical_input/beast_analyses/brachiopods_FBD_interval_ages.xml)
* [SA FBD skyline template](https://github.com/gavryushkina/sampled-ancestor-bdsky/blob/master/templates/SAFBDSKY.xml)
* [Contraband xml template](https://github.com/fkmendes/contraband/blob/master/examples/testing/BMPruneLikelihood.xml) (multiple traits, BM model)
* [Contraband xml template](https://github.com/fkmendes/contraband/blob/master/examples/testing/BMMVNShiftLikelihoodOneTrait_FBDTree_RateCatClock.xml) (multiple clock rate categories, nCat = 2)


