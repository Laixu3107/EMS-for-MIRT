# Latent Variable Selection in MIRT Models using the EMS Algorithm

## Introduction

The latent variable selection for multidimensional item response theory (MIRT) models
is equivalent to the statistical model selection in the presence of missing data.
In our paper (Xu et al., 2021), we consider the multidimensional two-parameter logistic (M2PL) model,
and apply the expectation model selection (EMS; Jiang et al., 2015) algorithm to find the optimal model (i.e., the incidence matrix)
and the parameter estimates (including the item discrimination and difficulty parameters) under optimal model which results in the smallest BIC value.

There are three files in the **R-code** directory:

- **MIRT_EML1.R** implements the EM-based L1 regularization for MIRT (Sun et al., 2016).
- **MIRT_EMS.R** implements the EMS algorithm for MIRT (Xu et al., 2021).
- **MIRT_EMS_sp.R** implements the accelerating version of EMS for MIRT (Xu et al., 2021).

Note that, all the implementations above are based on the R package *glmnet* (Friedman et al., 2010). In addition, each of them contains a simple example at the bottom of the file, which can be run directly for the users to see the details of the entire iteration procedure and the results including the parameter estimates under the optimal model, the computing time and etc.



## Citation
To cite these codes in publications, please use the following reference:

Xu, P. F., Shang, L., Zheng, Q. Z., Shan, N., & Tang, M. L. (2021). Latent variable selection in multidimensional item response theory models using the EMS algorithm. Submitted for publication. https://github.com/Laixu3107/EMS-for-MIRT.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2010). “Regularization paths for generalized linear models via coordinate descent.” *Journal of Statistical Software*, 33(1), 1–22. https://www.jstatsoft.org/v33/i01/.

Jiang, J., Nguyen, T., & Rao, J. S. (2015). The E-MS algorithm: model selection with incomplete data. *Journal of the American Statistical Association*, 110(511), 1136-1147.

Sun, J., Chen, Y., Liu, J., Ying, Z., \& Xin, T. (2016). Latent variable selection for multidimensional item response theory models via $L_{1}$ regularization. *Psychometrika*, 81(4), 921-939.

Xu, P. F., Shang, L., Zheng, Q. Z., Shan, N., & Tang, M. L. (2021). Latent variable selection in multidimensional item response theory models using the EMS algorithm. Submitted for publication. 


