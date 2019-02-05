# CBPS and Entropy Balancing: A Simulation Study

See .pdf for paper. Alternatively, see [blog post](https://chansoosong.com/2018/12/11/cbps-and-entropy-balancing-simulation-study/).

A challenge in the application of propensity scores for matching is that the propensity score is unknown and must be estimated. Whether the propensity score achieves good balance depends on correct specification of the propensity score model. However, the correct specification is unknown and we check the balance properties to evaluate the propensity score model specification. This is what [Imai (2008)](https://gking.harvard.edu/files/abs/matchse-abs.shtml) refers to as the 'propensity score tautology': the propensity score works when it works. 

In this simulation study, I analyze two approaches that seek to bypass this 'propensity score tautology': [Covariate Balancing Propensity Score](https://imai.fas.harvard.edu/research/files/CBPS.pdf) and [Entropy Balancing](https://web.stanford.edu/~jhain/Paper/PA2012.pdf).

