# Caompare complexity effect models: 
# LOO: (out-of-sample) predictive performance for new observations that are produced by the same data-generating process. 
# Predictive performance can be estimated by calculating the expected predictive performance (Gelman et al., 2014b; Vehtari and Ojanen, 2012).
# leave-one-out cross-validation (LOO-CV; Geisser and Eddy (1979))
library("loo")

m_normal <- readRDS(file="stanout/normalbycond.rda")
m_lognormal <- readRDS(file="stanout/lognormalbycond.rda")
m_mog <- readRDS(file="stanout/mogbycond.rda")
m_moglog <- readRDS(file="stanout/mogbycondlognormal.rda")
m_mogeachcond <- readRDS(file="stanout/mogeachcondnormal.rda")


# Extract pointwise log-likelihood and compute LOO
# the efficient PSIS-LOO approximation to exact LOO-CV
# The smaller the elpd, the better the fit

(models <- ls(pattern = "m_" ))

m <- m_normal
log_likout <- extract_log_lik(m)
looout <- loo(log_likout)
loo_mnormal <- looout

m <- m_lognormal
log_likout <- extract_log_lik(m)
looout <- loo(log_likout)
loo_mlognormal <- looout

m <- m_mog
log_likout <- extract_log_lik(m)
looout <- loo(log_likout)
loo_mmog <- looout



plot(loo_mnormal, label_points = FALSE)
print(loo_mnormal)

plot(loo_mlognormal, label_points = FALSE)
print(loo_mlognormal)

plot(loo_mmog, label_points = FALSE)
print(loo_mmog)

# Compare: diff elpd (positive difference prefers second model)
diff_mn.mln <- compare(loo_mnormal, loo_mlognormal)
print(diff_mn.mln)
diff_mln.mog <- compare(loo_mnormal, loo_mmog)
print(diff_mln.mog)

# WAIC
#(waic1 <- waic(log_lik_1))
#(waic2 <- waic(log_lik_2))
#print(compare(waic1, waic2), digits = 2)
