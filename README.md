# BNPclassification

This repostitory contains the files need to run the methods outlined in the manuscript "A Bayesian Nonparametric Model for Classification of Longitudinal Profiles" by Jeremy Gaskins, Claudio Fuentes, and Rolando de la Cruz (https://arxiv.org/abs/1711.01512).

This contains the following files:
<ul>
  <li> <b>BNCLD_sim1_1.RData</b> : This is a representative dataset for the methodology.  It comes from the simulation study described in Section 4 of the manuscript.   </li>
  <li> <b>BNPclassificiation - two component model.R</b> : This code runs the MCMC analysis of the two-component version of the model.  The MCMC output from this analysis is also used to initialize the adaptive Metropolis parameters in the "full BNP model.R" code. </li>
  <li> <b>BNPclassificiation - two component MCMC output analysis.R</b> : This code is used for the analysis of the MCMC output produced by the "two component model.R" code.  It will produce parameter summaries, graphical diagnositics, and within sample error rates. </li>
  <li> <b>BNPclassificiation - full BNP model.R</b> : This code runs the MCMC analysis of the full BNP model that mixes over the patient clusterings.  This model requires the MCMC output from the "two component model.R" analysis to initialize the adaptive Metropolis parameters. </li>
  <li> <b>BNPclassificiation - BNP stage 2 model.R</b> : This code performs the stage 2 BNP model that estimates the optimal clustering and runs MCMC with the fixed choice of the patient clusterings.  This model requires the MCMC output from the "full BNP model.R" analysis to estimate the optimal clustering. </li>  <li> <b>BNPclassificiation - MCMC output analysis.R</b> : This code is used for the analysis of the MCMC output produced by "full BNP model.R" and "BNP stage model.R" codes.  It will produce parameter summaries, graphical diagnositics, and within sample error rates. </li>
</ul>
