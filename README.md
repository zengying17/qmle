# qmle
Matlab program to perform Quasi-maximum likelihood estimation for the linear-in-means peer effects model with group effects. The estimation strategy and the properties of the estimator are discussed in "Kuersteiner, G. M., Prucha, I. R., & Zeng, Y. (2021). Efficient Peer Effects Estimators with Random Group Effects. ArXiv:2105.04330 [Econ, Stat]. http://arxiv.org/abs/2105.04330". The program includes one main function and four dependent functions:
- limpe_re_unc: main function for the QMLE. 
- lnf_unc, lnf_c, mom34fun, ghfun, VCfun: dependent functions called by the main function to calculate the objective function and variance-covariance matrix.

Below is the syntax of the main function.

      [a,b]=limpe_re_unc(datatable, yvar, groupvar, xivar, xpvar, hetvar, options, x0);

The outputs are a (results) and b (sample), which can be omited.
- a includes all results. For example, a.restabl is the table of estimates and standard errors for the parameters. a.lambda is the estimate for endogenous peer effect lambda. 
- b is the final sample used for estimation. 

The inputs include the following:
- datatable is the data table;
- yvar is the dependent variable;
- groupvar is the variable for group id; 
- xivar is the list of variables for indiviudal charateristics, i.e., x1 in the paper;   
- xpvar is the list of variables whose leave-out-mean are added as indepednet variables, i.e., x2 in the paper;  
- hetvar is the categorical variable that define group types, i.e., Dr in the paper;  
- options are optimization options in matlab, which can be left empty;  
- x0 is the starting value of the optimization, which can be left empty.  
