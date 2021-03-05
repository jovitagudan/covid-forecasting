# Restricted OLS described here http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xegbohtmlnode18.html
source('jhu_data.R')
source('common.R')


nbasis = 7
country = 'Lithuania'

lt_acc_data = t(normalized[country,])
lt_recov_data = t(normalized_recov[country,])
lt_deaths_data = t(normalized_deaths[country,])


y = lt_acc_data - lt_recov_data - lt_deaths_data

x = seq(-2, 2, length.out = length(y))
# evaluated basis
X = hermite(x, nbasis = nbasis)

point_count = length(y)
restriction_idx = c(
  1, point_count
)

beta = solve(t(X) %*% X) %*% t(X) %*% y

R = hermite(x[restriction_idx], nbasis=nbasis)
r = y[restriction_idx]

XX_inv = solve(t(X) %*% X)
p1 = t(R) %*% solve( R %*% XX_inv %*% t(R)   )
p2 = (R %*% beta - r)

beta_rls = beta - XX_inv %*% p1 %*% p2

plot(x, y, type='l')
lines(x, X %*% beta, col='red')
lines(x, X %*% beta_rls, col='blue')
