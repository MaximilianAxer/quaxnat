
# Load quaxnat package:
library(quaxNat)

# Load and check data
data("quax_data")
#After loading the quaxNat package, we are now ready to apply the quax function,
#which we do for the .999th quantile and five dispersal kernels implemented in
#quaxnat.

# Estimate regeneration potential based on various dispersal kernels:
tau <- 0.995
f1 <- quax(oak_regen ~ distance_oak, data, tau=tau, fun=k_t)
f2 <- quax(oak_regen ~ distance_oak, data, tau=tau, fun=k_weibull)
f3 <- quax(oak_regen ~ distance_oak, data, tau=tau, fun=k_lognormal)
f4 <- quax(oak_regen ~ distance_oak, data, tau=tau, fun=k_power)
f5 <- quax(oak_regen ~ distance_oak, data, tau=tau,
           fun=k_exponential_power)
f <- list(`Spatial t`=f1, Weibull=f2, Lognormal=f3, Power=f4,
          `Exp. Power`=f5)

# Plot regeneration density as a function of the distance to
# the nearest seed tree:
plot(oak_regen ~ distance_oak, data, xlim=c(0,1500), cex=.8)

# Pick some colors:
col <- hcl.colors(length(f), palette="Dynamic")

# Add estimated functions to diagram:
for (i in seq_along(f)){
  curve(f[[i]](x), add=TRUE, n=10000, col=col[i], lwd=2)}

# Add legend:
legend("topright", title=paste0(tau," quantile"),
       legend=names(f), lty=1, lwd=2, col=col)

# Compare quality of fits:
sapply(f, summary)


