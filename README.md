# viLDA
Varitional Inference Latent Dirichlet Allocation. viLDA is an R package that estimates posterior distributions for LDA models.

## Project Description
viLDA is an R package that can estimate LDA models using Markov Chain Monte Carlo (MCMC) methods and Variational Inference methods. In particular, we implement the Gibbs Sampler designed by Griffiths and Steyvers (2004) and the Stochastic Variational Inference (SVI) algorithm developed by Hoffman et al (2013). The Gibbs sampler was mainly used as a baseline comparison and proof of concept, and the main focus of the package was to write optimized R and C++ code using the Rcpp interface to run SVI to estimate LDA models. We compare the performance of our package to the topicmodels package, which runs the Variational Expectation Maximization algorithm designed by Blei et al. (2003).

Speed is a high priority for viLDA, which is why we developed one of the fastest known algorithms for LDA posterior distributions (SVI) and wrote the performance critical code in C++. A challenge that we are currently still working on is scaling up SVI to work on data sets with very large vocabularies.

## Using Our Project
To install our package, you can visit the colab page at the URL http://bit.ly/viLDA and follow the steps to install the package using dev tools and run the main function.

## Acknowledgements
We would like to thank Dr. Hyun Min Kang for teaching requisite knowledge for us to undertake this project and for advising our project as part of the Statistical Computing course at the University of Michigan Ann Arbor.

## References
David M. Blei, Andrew Y. Ng, and Michael I. Jordan. Latent dirichlet allocation. J. Mach. Learn. Res., 3: 993–1022, March 2003, https://www.jmlr.org/papers/volume14/hoffman13a/hoffman13a.pdf
Hoffman, Matthew D., et al. Stochastic Variational Inference.. Journal of Machine Learning Research 14 (2013) 1303-1347, https://www.jmlr.org/papers/volume14/hoffman13a/hoffman13a.pdf
Thomas L. Griffiths and Mark Steyvers. Finding scientific topics. PNAS April 6, 2004 101 (suppl 1) 5228-5235; https://doi.org/10.1073/pnas.0307752101.
