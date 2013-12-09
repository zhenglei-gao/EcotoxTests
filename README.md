EcotoxTests
===========


### 2.1  Statistical Power Analysis

The power of a statistical test is the probability of rejecting a false null hypothesis or accepting a true alternative hypothesis. As the power increases, we have more chance to detect an effect if there exists an effect and the chances of a Type II error occurring decrease. The probability of a Type II error occurring is referred to as the false negative rate ($\beta$). Therefore power is equal to 1 − $\beta$, which is also known as the sensitivity.

Power analysis can be used to calculate the minimum sample size required so that one can be reasonably likely to detect an effect of a given size. Power analysis can also be used to calculate the minimum effect size that is likely to be detected in a study using a given sample size. In addition, the concept of power is used to make comparisons between different statistical testing procedures: for example, between a parametric and a nonparametric test of the same hypothesis.

### 2.2  William's test and its power

The biological activity (i.e. toxicity) of a substance may be investigated by experiments in which the treatments are a series of monotonically increasing (or decreasing) doses of the substance and an overall dose related trend is to be expected. One of the aims of such a study is to determine the minimum effective dose (the lowest dose at which there is activity). Williams (1971, 1972) describes such a test. William's test is a multiple comparison test and is very sensitive in estimating a lower LOEC than other available test because it takes into account the order of concentrations according to their increase or decrease. The test used maximum likelihood procedure.

The null hypothesis is:

$$
H_0: \mu_0=\mu_1=...=\mu_k
$$

The alternative is:

$$
H_1: \mu_0\leq\mu_1\leq...\leq\mu_k
$$

with at least one strict inequality. 

The MLE of the treatment means $\mu_i$ are obtained by the sample means $\bar{x}_i$ and the number of observations in the samples $n_i$, by the formular

$$
\mu_i^*=\max_{l\leq u\leq i}\min_{i\leq v \leq k}\sum_{j=u}^{v}n_j \bar{x}_j/\sum_{j=u}^{v}n_j
$$

The test statistic is:

$$
T_i=\frac{\mu_i^*-\bar{x}_0}{\sqrt{s^2/n_i+s^2/n_0}}$$

where, $s$ is an unbiased estimated of $\sigma$, the within group standard deviation that is independent from $\bar X_i$. The null hypothesis is rejected and the fact that the $i$-th dose level is the minimum effectie dose is concluded if 

$$T_j> t_{j,\alpha}, \mbox{for all } j \geq i$$

where,  $t_{j,\alpha}$ is the upper $\alpha$th percentile of the distribution $T_j$. These critical values were partly tabulated by Williams(1972) and are obtained in my `EcotoxTests` package either using inverse interpolation routines or by simulations. 


The power of the Williams test can be calculated by approximation method or by simulation based method.

### Method 1: Approximation based:

The following approximate function for power is given by Chow et al. (2008) page 288.

$$
1-\beta=1-\Phi(t_{K,\alpha}-\frac{\Delta}{\sigma\sqrt{2/n}})
$$

where $\Delta$ is the clinically meaningful minimal difference. In the program, calculated by the difference between the NOEC and the cotrol level.

The approximate method assumes that there are $G+1$ level of dose concentrations and for each concentration $n$ samples are collected. The alternative hypothesis is at the $K$-th($K=1,...G$) concentration level, the difference in the means is differently from 0 by $\Delta$(almost the whole alternative space). 


### Method 2: Simulation based:

We also implement the simulation based power analysis, which makes more sense from my point of view. We pick an alternative hypothesis that at each dose level, the difference is a true difference. Then conditional on this alternative hypothesis being true, we generate data from this alternative model, perform the Williams' test on these simulated data. The statistical power is the proportion of p-values that are lower than a specified $\alpha$-level, which is the probability of accepting the alternative. 
