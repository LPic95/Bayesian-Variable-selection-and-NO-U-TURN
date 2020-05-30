[![HitCount](http://hits.dwyl.com/LPic95/Bayesian-Variable-selection-and-NO-U-TURN.svg)](http://hits.dwyl.com/LPic95/Bayesian-Variable-selection-and-NO-U-TURN)


Bayesian methods for sparse modeling: no “turning back” on the path to variable selection. 
 ----------------
<p align="justify">
In the statistical model definition and in the related building process, the selection of variables from a large group of regressors has recently become widespread. The crucial role is further expressed in “high dimensional contexts”, whose large datasets are so difficult to cope with. Therefore, sparsity is focal to recover the underlying signals in complex relationships between variables.
 </p>
 <p align="justify">
The main aims of this project are to compare salient Bayesian techniques with the most consolidated and well-known Lasso methodology and consequently, to jointly evaluate a juxtaposition between MCMC, Gibbs Sampling algorithm, and Hamiltonian Monte Carlo, specifically the No-U-Turn Sampler. 
This thesis is made up of two sections, with the first pinpointing and summarising the theoretical underpinnings of the main Bayesian methodologies, and the second focusing on empirical evidences. More in detail, the second chapter shifts its focus from modelling to simulation, with the purpose of studying performances in various simulation designs and consequently detecting the most appropriate methodology according to each context, e.g. predictors drawn from asymmetric and heavy-tailed distributions, multicollinearity and high dimensional data.
</p>
 <p align="justify">
In light of the  twofold purpose of this research, the “open” interpretation of the title is now clear: statistical learning can no longer disregard the selection of variables and, at the same time, the No-U-Turn Sampler is needful to reduce the Gibbs Sampling’s “random walk”.
</p>


 Research questions 
 ----------------

|`1`| Among the main Bayesian feature selection techniques is there a preferable methodology in forecasting and inferential terms with high degree of sparsity, multicollinearity and high dimensional data?|
|-|-|


|`2`| Does changing the procedure for approximating the posterior distribution improve forecasting accuracy? |
|-|-|


Methodologies
----------------

 **Sparse modelling** 
 
 * Least Absolute Shrinkage and Selection Operator
 * Spike and Slab 
   * continuous priors and SVS
   * rescaled Spike and Slab with a stability study 
 * Spike and Slab Lasso
 
**Sampler** 

* Gibbs sampling (Monte Carlo)
* No-U-Turn (Hamiltonian Monte Carlo)
<p align="justify">
The entire project is developed considering all the possible combinations of the sparse methodologies and posterior sampler algorithms shown. 
<p align="justify">
 
Simulation study
----------------
<p align="justify">
The data generator model is: 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?Y_{i}\&space;=\boldsymbol{x_{i}}^{t}&space;\boldsymbol{\beta}&space;&plus;&space;\varepsilon&space;_i,&space;\,\,\,\,\,&space;i=1,...,n," title="Y_{i}\ =\boldsymbol{x_{i}}^{t} \boldsymbol{\beta} + \varepsilon _i, \,\,\,\,\, i=1,...,n," />
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7Bx%7D_i%20%3D%20%28x_%7Bi%2C1%7D%2C...%2Cx_%7Bi%2C50%7D%29%5Et)
is the 50 dimensional vector of covariates, ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D%20%3D%20%28%5Cbeta_1%2C...%2C%5Cbeta_%7B50%7D%29%5Et) is the vector of regression coefficients and the random error ![](https://latex.codecogs.com/gif.latex?%5Cvarepsilon_i%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%280%2C4%20%5Cright%20%29%2C%20%5Cforall_i%3D1%2C..%2Cn).

The original simulation design is defined as follows:
* ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7BX%7D%20%5Csim%20%5Cmathbb%7BN%7D%5Cleft%20%28%20%5Ctextbf%7B0%7D%2C%5Ctextbf%7BI%7D%20%5Cright%20%29)
* the degree of sparsity ranges from 0.96 to 0.1
* sample size adopted (![](https://latex.codecogs.com/gif.latex?n)) is equal to 1000.
<p align="justify">
 
Variants on the simulation design
----------------
 <p align="justify">
 
1. ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D%20%5Csim%20%5Cmathbb%7BN%7D%5Cleft%20%28%20%5Ctextbf%7B0%7D%2C%5Csigma%5E2%5Ctextbf%7BI%7D%20%5Cright%20%29) with *![](https://latex.codecogs.com/gif.latex?%5Csigma) firstly equal to 3 and then 15.* 

2. Draw coefficients from *heavy-tailed distribution, Student’s T, and asymmetric probability densities, Skew Normal*.

3. Move beyond the *orthogonality of predictors*: ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7BX%7D%5Csim%20%5Cmathbb%7BN%7D%5Cleft%20%28%20%5Ctextbf%7B0%7D%2C%5Cboldsymbol%7B%5CSigma%7D%20%5Cright%20%29) with a covariance block-diagonal matrix (![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5CSigma%7D)) where each 10 x 10 size sub-matrix has unit variance and fixed covariance. *The relative correlation ranges from 0.90 to 0.10*.

4. Variation in the *number of observations (![](https://latex.codecogs.com/gif.latex?n)): the focus lies on high-dimensional context, e.g. 20 observational units and 50 predictors*.

5. In the last simulation design, *the cross product between ![](https://latex.codecogs.com/gif.latex?n%20%3D%20%5C%7B20%2C%2030%2C%2040%2C%2050%5C%7D) and ![](https://latex.codecogs.com/gif.latex?corr%20%3D%20%5C%7B0.90%2C%200.65%2C%200.4%2C%200.2%2C%200.1%5C%7D) vectors* is considered. The applied degree of sparsity is 0.8.

</p>

Results
----------------

| :memo:        | Among the main Bayesian feature selection techniques is there a preferable methodology in forecasting and inferential terms with high degree of sparsity, multicollinearity and high dimensional data?|
|---------------|:------------------------|


* SSL seems to be robust and to perform adequately both from the inferential and forecasting standpoint, regardless the degree of sparsity, especially with a strong correlation among predictors and high dimensionality contexts.

 <p align="center"> <sub>Test mse values of different estimation methods and different sampling procedures over 100 simulated test datasets varying the sample size.</sub></p>

<table style="border-collapse:collapse;border-spacing:0;table-layout: fixed; width: 568px" class="tg"><colgroup><col style="width: 77px"><col style="width: 77px"><col style="width: 65px"><col style="width: 73px"><col style="width: 69px"><col style="width: 69px"><col style="width: 69px"><col style="width: 69px"></colgroup><thead><tr><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal" colspan="3">  Methodology</th><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">1000 obs.</th><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">100 obs.</th><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">50 obs.</th><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">30 obs.</th><th style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;font-weight:normal;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">20 obs.</th></tr></thead><tbody><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal" colspan="2" rowspan="2">Bayesian Lasso</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">NUTS</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.188943</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">6.081969</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">14.72189</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">26.52780</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">15.43429</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">G.S</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.251924</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">9.259484</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal">22.30086</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">27.46183</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">37.31440</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal" colspan="2" rowspan="2">George and McCulloch</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">NUTS</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.198545</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">7.434634</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">17.08161</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">57.90704</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">42.52329</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">G.S</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.093047</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">5.172765</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">20.64758</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">45.62200</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">31.60901</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal" colspan="2" rowspan="2">Ishwaran and Rao</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">NUTS</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.220462</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">7.524676</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">18.21790</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">33.33155</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">53.73723</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">G.S</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.103748</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.562725</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">23.67782</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">29.20373</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">31.46803</span></td></tr><tr><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal" colspan="2">SSL</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">G.S</td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.083747</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">4.410538</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">11.23710</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">15.42313</span></td><td style="border-color:#ffffff;border-style:solid;border-width:0px;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal"><span style="font-weight:normal;font-style:normal;text-decoration:none">28.48944</span></td></tr></tbody></table>

This statement is further supported by the inferential metrics (sensitivity, specificity, accuracy, Hamming distance) and by the following graphics.


<p align="center">
<img src="Grafici%20github/20_obs.png" width="350" height="200"/> <img src="Grafici%20github/30_obs.png" width="380" height="200"/> 
 </p>
 <p align="center">
<img src="Grafici%20github/50.png" width="350" height="200"/> <img src="Grafici%20github/1000_obs.png" width="380" height="200"/>
 </p>


 
  <p align="center"> <sub> Graphical comparison among estimates and true coefficients by varying the number of observations over 50 predictors.</sub></p>


| :memo:        | Does changing the procedure for approximating the posterior distribution improve forecasting accuracy?       |
|---------------|:------------------------|
<p align="justify">
 
* NO-U-TURN algorithm results in lower computational and running time, thanks to only 8 chains compared to the 200 involved in Gibbs sampling procedure and thanks to the specific programming language, Stan(C++).

* NUTS Bayesian Lasso is the best method in high dimensionality setting with moderate and fairly high correlation among predictors.

* In all the simulation designs, the implementation of the Bayesian Lasso methodology with the NO-U-TURN algorithm results in more accurate Lasso predictions than the respective ones obtained through Gibbs Sampling.
<p>
