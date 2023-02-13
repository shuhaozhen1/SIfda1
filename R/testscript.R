coef_list <- list(
  coef1 = function(x){0*x},
  coef2 = function(x){sin(2*pi*x)},
  coef3 = function(x){x},
  coef4 = function(x){-x},
  coef5 = function(x){x^2}
)

VCMdata <- VCM_process_values_error(n=100, m=10, p = 5, domain=c(0,1), mean_list = rep(list(function(x) {0*x}), 5),
                                    coef_list=coef_list, covariancef=function(x,y){exp(-abs(x-y))},
                                    distribution = 'normal', snr = 5, sig = NULL, depend = FALSE,
                                    trans = FALSE, transfmatrix = NULL, num_basis = 1000)

sin(2*pi*0.2)
