# Code for the paper "REFINE2: A tool to evaluate real-world performance of machine-learning based effect estimators for molecular and clinical studies"

Code to reproduce results in the paper [REFINE2: A tool to evaluate real-world performance of machine-learning based effect estimators for molecular and clinical studies](https://arxiv.org/abs/2105.13148) by Xiang Meng, Jonathan Huang

Xiang Meng 2021.06.11


[![DOI](https://zenodo.org/badge/373040701.svg)](https://zenodo.org/badge/latestdoi/373040701)

---

## How to use the app?
1. Clone [the folder](https://github.com/mengeks/drml-plasmode/tree/master/REFINE2) on your local machine  
2. Convert your data set as a CSV file. Make sure to name the outcome as `Y` and the binary treatment as `A`  
3. Make sure you installed the package **Shiny** in R  
4. In R session, run `shiny::runApp(PATH)`, where `PATH` is the path of the folder  
   - Example: `shiny::runApp("~/Desktop/HuangGroup/cvtmle_plasmode/Code/REFINE2")`  
5. Set the “Path of the data” to your data path  
6. Set 4 models and other parameters  
   - Models should be in the format of R formula: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/formula.html  
7. Click the **Run** button  

### Demo Screenshot
![REFINE2 Screenshot](REFINE2/REFINE2-Screenshot.png)


⚠️ **Note (Updated Aug 2026):** The app has been primarily tested on **macOS (OSX)**, however it should also be fully functional on PC/Windows thanks to fixes proposed by @sbresnahan. If Windows users are still encountering issues, please let us know!


---


## 1. Code Structure
* [2020309-Plasmode_Sim_CV-TMLE.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/2020309-Plasmode_Sim_CV-TMLE.R): Master code to first generate dataset and then produce estimation. It has dependencies on all other codes. 
* [20200803-Sims-Function.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200803-Sims-Function.R): Code to generate simulation dataset and plasmode dataset
* [20200226-PlasmodeCont_Revised.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200226-PlasmodeCont_Revised.R): Function to generate Plasmode dataset
* [20200904-run-sim-code.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200904-run-sim-code.R): Code to parallelize estimation. 
* [20200720-Algos-code.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200720-Algos-code.R): Function that contains all estimators
* [20200705-DCDR-Functions.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200705-DCDR-Functions.R): Function for Double Cross-Fit which is modified from Breskin and Zivich (2020).
* [20200816-Result-Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200816-Result-Summary.R): A script containing all summary functions.
* [3.1 Single-Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/3.1%20Single-Summary.R): Code that produces summary for a single result file (need to manually change path)
* [3.Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/3.Summary.R): Code that produces a list of summaries and plots for the paper.
* [OutputCoeff.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/OutputCoeff.R): Code to produce the actual simulation coefficients as shown in the appendix.

## 2. Results in Paper and Corresponding Code
As briefly mentioned above, parameters in [2020309-Plasmode_Sim_CV-TMLE.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/2020309-Plasmode_Sim_CV-TMLE.R) controls all the simulations. To run for different simulation scenario, one should just change parameters here. eg. "idx.handpick" is the indices of covariates picked; "interact" means whether we add first order interactions between covariates; "interact.w.exp" means we add first order interactions between covariates and the exposure variable. 


## 3. The REFINE2 Shiny app
The REFINE2 tool allows the analyst to check whether the chosen effect estimation strategy, with or without the use of machine learning algorithms, will reliably estimate an unbiased average treatment effect (ATE) and appropriate confidence intervals given their specific dataset and user-provided parameters. Namely, the user provides: 
* the dataset  
* number of simulation draws
* **Simulation** models representing the "true" relationships between treatment, outcome, and covariates (in R syntax)
* **Estimation** models representing how the effect would be estimated
* chosen estimator (*e.g.* IPW, DC-TMLE)
* library used to fit the estimator (smooth, non-smooth)

Allowing the user to specify both Simulation and Estimator models ensures maximum flexibility on degree of misspecification. For example, variables could be used to simulate the truth but be omitted or mismeasured in the estimation model to test sensitivity. Using the exact same Simulation and Estimation models will test the finite sample / bias convergence properties of the algorithms on the dataset at hand.

## 4. R Package Dependencies

REFINE2 automatically installs required dependencies when first launched. Key packages include:

**Core Dependencies:**
- `shiny` - Web application framework
- `SuperLearner` - Ensemble learning (includes `glmnet`, `randomForest`, `xgboost`, `polspline`)
- `tidyverse` - Data manipulation and visualization (includes `dplyr`, `ggplot2`, `readr`, etc.)

**Additional Packages:**
- `mgcv`, `splines` - Statistical modeling
- `doParallel`, `foreach` - Parallel computing
- `readxl`, `haven` - Data import
- `table1` - Descriptive statistics

All dependencies are checked and installed automatically on first run.
