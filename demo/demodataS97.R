#read example data
data(dataS97)

#example with a few competing covariates
dataS97.bvs <- bvs.lm(formula = "gr56092~gdpsh60l+lifee060+p60",
                      data = dataS97, n.keep = 2^3, imp.seed = 123)

#print the result
dataS97.bvs

#summary of the result
summary(dataS97.bvs)

#plots for the result
plot(dataS97.bvs)
