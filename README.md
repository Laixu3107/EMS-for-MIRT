# Latent Variable Selection in MIRT Models using the EMS Algorithm

## Introduction
The latent variable selection for multidimensional item response theory (MIRT) models
is equivalent to the statistical model selection in the presence of missing data.
In our paper (Xu et al., 2021), we consider the multidimensional two-parameter logistic (M2PL) model,
and apply the expectation model selection (EMS; Jiang et al., 2015) algorithm to find the optimal model (i.e., the incidence matrix)
and the parameter estimates (including the item discrimination and difficulty parameters) under optimal model which results in the smallest BIC value.

There are six files in the **M2PL_Rcode** directory:

- **M2PL_EML1.R**  implements the EM-based L1 regularization for MIRT when $\Sigma$ is known (Sun et al., 2016).
- **M2PL_EMS.R**   implements the EMS algorithm for MIRT when $\Sigma$ is known (Xu et al., 2021).
- **M2PL_EMSSP.R** implements the accelerating version of EMS for MIRT when $\Sigma$ is known (Xu et al., 2021).
- **M2PL_EMS_sigma_unknown.R**  implements the EMS algorithm for MIRT when $\Sigma$ is unknown (Xu et al., 2021).
- **M2PL_EMSSP_sigma_unknown.R** implements the accelerating version of EMS algorithm for MIRT when $\Sigma$ is unknown (Xu et al., 2021). 
- **M2PL_Esti_sigma.cpp** is used for updating $\Sigma$ in each MS-step when $\Sigma$ is unknown. This code is written based on that in R-package lvmcomp (Zhang el al., 2020).

Note that, all the first five implementations above are based on the R package *glmnet* (Friedman et al., 2010).
In addition, each of them contains a simple example at the bottom of the file,
which can be run directly for the users to see the details of the entire iteration procedure and the results including the parameter estimates under the optimal model, the computing time and etc.

## Citation
To cite these codes in publications, please use the following reference:

Xu, P. F., Shang, L., Zheng, Q. Z., Shan, N., & Tang, M. L. (2021). Latent variable selection in multidimensional item response theory models using the EMS algorithm. Submitted for publication. https://github.com/Laixu3107/EMS_for_MIRT.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2010). “Regularization paths for generalized linear models via coordinate descent.” *Journal of Statistical Software*, 33(1), 1–22. https://www.jstatsoft.org/v33/i01/.

Jiang, J., Nguyen, T., & Rao, J. S. (2015). The E-MS algorithm: model selection with incomplete data. *Journal of the American Statistical Association*, 110(511), 1136-1147.

Sun, J., Chen, Y., Liu, J., Ying, Z., \& Xin, T. (2016). Latent variable selection for multidimensional item response theory models via $L_{1}$ regularization. *Psychometrika*, 81(4), 921-939.

Xu, P. F., Shang, L., Zheng, Q. Z., Shan, N., & Tang, M. L. (2021). Latent variable selection in multidimensional item response theory models using the EMS algorithm. Submitted for publication. 

Zhang, S., Chen, Y., \& Liu, Y. (2020). An improved stochastic EM algorithm for large-scale full-information item factor analysis. *British Journal of Mathematical and Statistical Psychology*, 73(1), 44-71.
