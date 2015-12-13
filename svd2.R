rm(list=ls())
#install.packages("softImpute")
library(softImpute)

#x = simulate_matrix(5,1000,1000)$M
#system.time(total(x, 0.7,5))

#randomly select the base matrix to cover
base_matrix = function(x, base)
{ 
  n=nrow(x)
  p=ncol(x)
  np=n*p
  ix=seq(np)
  notmiss=sample(ix,np*base,replace=FALSE)
  xna=matrix(NA, nrow = n, ncol = p)
  xna[notmiss]=x[notmiss]
  return (list(xna, notmiss))
}

#add the additional part for query ofr svd
query = function(base_xna, x, add, notmiss,r, image = FALSE)
{ n=nrow(x)
  p=ncol(x)
  np=n*p
  ix=seq(np) 
  error = c()
  j = 1
  i = add/100
  xna = base_xna 
  notmiss2=sample(ix[-notmiss],np*i,replace=FALSE)
  xna[notmiss2]=x[notmiss2]
  fit1=softImpute(xna,rank=r)
  comp = complete(xna,fit1) 
  error = (sum((comp - x)^2))/(sum(x^2))
  j = j+1
  if (image == TRUE)
  {
    ynn = t(comp[nrow(comp):1,])
    image(ynn, col = grey(seq(0, 1, length = 256)))
  }
  return(list(error = error, xna = xna, comp = comp))
}

library("pixmap")
setwd("~/Desktop/sta289/final/")
library("png")
library("colorspace")
x <- readPNG("image.png")[,,1]
#x = simulate_matrix(2,600,400)$M
base_result = base_matrix(x, base = 0.30045067987837814627)
result = query(base_result[[1]], x, 0, base_result[[2]],5, image = TRUE)
cover_svd = result$ynn
miss_svd = result$xna
error_svd = result$error
basematrix = base_result[[1]]
save(miss_svd, x,file = "svd.RData")
#matrix_havemiss = result$xna

#calculate error and plot
error = c()
mismatrix = list()
for (i in 0:15)
{ result = query(base_result[[1]], x, i, base_result[[2]],2)
  error[i+1] = result[[1]]
  mismatrix[[i+1]] = result$xna
}

plot(35:50,error, type = "l")
save(error,file = "error_svd_xingxing.RData")
basematrix = base_result[[1]]
save(mismatrix, x, file = "matrix.RData")

#install.packages("pheatmap")
yn = t(result[[2]][nrow(x):1,])
image(yn, col = grey(seq(0, 1, length = 256)))

#import x to server
#inserver
#pmf method
load("matrix.RData")
pmfmatrix = lapply(mismatrix, function(xx) {xx[is.na(xx)] = 0; xx})
source("/home/yinz/libpmf-1.41/R/libpmf.R")
pmf_error = function(matrix, true_matrix){
model = pmf.train.matrix(matrix, '-k 2')
comp1 = model$W %*% t(model$H)
comp = comp1 + matrix
error = (sum((comp - true_matrix)^2))/(sum(true_matrix^2))
return(list(error = error, comp = comp))
}

pmf_error = sapply(pmfmatrix, function(xx) pmf_error(xx, x))
save(pmf_error, file = "pmf_error.RData")

#my computer
load("pmf_error.RData")
b = 35:50
plot(x = b, y = error, type = "o", ylim = range(0, 0.45), col = "red", ylab = "Error", xlab = "Querry percent(%)")
lines(x = b, y = pmf_error, col = "blue", type = "o")
lines(x = b, y = qerror, col = "green", type = "o")
legend("topright", legend = c("SVP", "PMF", "Querry"), col = c("red", "blue", "green"), lty = 1, cex = 0.7)










