
source('common.R')
source('jhu_data.R')


incompatible_curves = c(42,  43,  80, 104, 186)


#accumulated new cases per 100,000 inhabitants
cov_data = normalized[-incompatible_curves,]
y_len = dim(cov_data)[2]
country_count = dim(cov_data)[1]
arg_range = c(-2,2)
x = seq(arg_range[1], arg_range[2], length.out = y_len)

nbasis=5
wbasis = create.hermite.basis(arg_range, nbasis=nbasis)
cvec0 = matrix(0, nbasis, country_count)
Wfd0 = fd(cvec0, wbasis)


Lfdobj    <- 10        #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
covfdPar <- fdPar(Wfd0, Lfdobj, lambda)
fdCov_data <-smooth.monotone(x, t(cov_data), covfdPar)


Dvalues = eval.monfd(x, fdCov_data$Wfdobj, 1)
colnames(Dvalues) <- rownames(cov_data)
fd_eval = eval.monfd(x, fdCov_data$Wfdobj)

fd_values = fdCov_data$beta[1,] + fdCov_data$beta[2,] * t(fd_eval)
rownames(fd_values) <- rownames(cov_data)

plot(fd_values['Lithuania',], type='l')
plot(Dvalues[,'Lithuania'], col='blue', type='l')

matplot(t(fd_values), type="l")
matplot(Dvalues, type="l")


#Recovered

cov_data = normalized_recov
ii = c()
for (i in 1:254) {
  res <- tryCatch({
    y = t(cov_data)[,i]
    y = cumsum(y)
    x = seq(arg_range[1], arg_range[2], length.out = length(y))
    cur = smooth.monotone(x, y, covfdPar)
  }, warning = function(w) {
  }, error = function(e) {
    print(e)
    ii = c(ii, i)
    print('Failed')
    return(NA)
  }, finally = {
  })
  if (length(res)<2) {
    ii = c(ii, i)
    print(ii)
  }
}




incompatible_curves = c(90, 172, 175, 221)

rec_data = normalized_recov[-incompatible_curves,]
y_len = dim(rec_data)[2]
country_count = dim(rec_data)[1]
arg_range = c(-2,2)
x = seq(arg_range[1], arg_range[2], length.out = y_len)

nbasis=5
wbasis = create.hermite.basis(arg_range, nbasis=nbasis)
cvec0 = matrix(0, nbasis, country_count)
Wfd0 = fd(cvec0, wbasis)


Lfdobj    <- 10        #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
recfdPar <- fdPar(Wfd0, Lfdobj, lambda)
fdRec_data <-smooth.monotone(x, t(rec_data), recfdPar)


Dvalues_rec = eval.monfd(x, fdRec_data$Wfdobj, 1)
colnames(Dvalues_rec) <- rownames(rec_data)
fd_eval_rec = eval.monfd(x, fdRec_data$Wfdobj)

fd_values_rec = fdRec_data$beta[1,] + fdRec_data$beta[2,] * t(fd_eval_rec)
rownames(fd_values_rec) <- rownames(rec_data)

plot(fd_values_rec['Lithuania',], type='l')
plot(Dvalues_rec[,'Lithuania'], col='blue', type='l')

matplot(t(fd_values_rec), type="l")
matplot(Dvalues_rec, type="l")




#Deaths

cov_data = normalized_deaths
ii = c()
for (i in 1:length(rownames(cov_data))) {
  res <- tryCatch({
    cur = smooth.monotone(x, t(cov_data)[,i], covfdPar)
  }, warning = function(w) {

  }, error = function(e) {
    print(e)
    ii = c(ii, i)
    print('Failed')
    return(NA)
  }, finally = {

  })

  if (length(res)<2) {
    ii = c(ii, i)
    print(ii)
  }
}



incompatible_curves = c(11,  27,  38,  42, 43,
                        47,  49,  74,  78,  79,
                        80,  84,  87, 102, 103,
                        104, 106, 112, 124, 127,
                        135, 141, 161, 176, 182,
                        186, 213, 215, 216, 222,
                        227, 242, 251, 256, 264)

dth_data = normalized_deaths[-incompatible_curves,]
y_len = dim(dth_data)[2]
country_count = dim(dth_data)[1]
arg_range = c(-2,2)
x = seq(arg_range[1], arg_range[2], length.out = y_len)

nbasis=5
wbasis = create.hermite.basis(arg_range, nbasis=nbasis)
cvec0 = matrix(0, nbasis, country_count)
Wfd0 = fd(cvec0, wbasis)


Lfdobj    <- 10       #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
dthfdPar <- fdPar(Wfd0, Lfdobj, lambda)
fdDth_data <-smooth.monotone(x, t(dth_data), dthfdPar)


Dvalues_dth = eval.monfd(x, fdDth_data$Wfdobj, 1)
colnames(Dvalues_dth) <- rownames(dth_data)
fd_eval_dth = eval.monfd(x, fdDth_data$Wfdobj)

fd_values_dth = fdDth_data$beta[1,] + fdDth_data$beta[2,] * t(fd_eval_dth)
rownames(fd_values_dth) <- rownames(dth_data)


matplot(t(fd_values_dth), type="l")
matplot(Dvalues_dth, type="l")
