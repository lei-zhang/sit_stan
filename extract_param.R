a <- summary(out1)
a <- a$summary

# get the column-names and row-names
colnames(a)
rownames(a)

parMean <- as.matrix(a[,1])
parMedn <- as.matrix(a[,6])


## calculating is a bit tricky
library(modeest)





#### extract mcmc ####
lr_mu <- extract (out1, "lr_mu")$lr_mu  # returns an array









