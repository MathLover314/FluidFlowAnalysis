library(rmatio)
library(imager)

fluids = read.mat("DATA\\FLUIDS\\CYLINDER_ALL.mat")

fluids_velocity = as.matrix(fluids$UALL)
svd_decomp = svd(fluids_velocity)
cutoff = (4/sqrt(3)) * sqrt(151)
r = 50

egien_flow = diag(svd_decomp$d[1:r]) %*% t(svd_decomp$v[,1:r])
plot(as.cimg(egien_flow[,1:20]))
W = egien_flow[,1:150]
Wp = egien_flow[,2:151]
svd_w = svd(W)
A = Wp %*% svd_w$v %*% solve(diag(svd_w$d)) %*% t(svd_w$u)

W_res = matrix(0, nrow = r, ncol=151)
W_res[,1] = W[, 1]
for (i in 2:151) {
  W_res[,i] = predict_flow(i)
}
norm(egien_flow - W_res)
plot(as.cimg(W_res[,1:20]))

predict_flow <- function(k) {
  return(A^(k-1) %*% W[,1])
}

decomp_rank <- function(r) {
  return(svd_decomp$u[,1:r] %*% diag(svd_decomp$d[1:r]) %*% t(svd_decomp$v[,1:r]))
}

B = numeric(100)
for (i in 2:100) {
  B[i] = norm(fluids_velocity - decomp_rank(i)) 
}

fluids_cut = decomp_rank(10)
get_state <- function(t) {
  return(fluids_cut[(1+199*t):(199+199*t),])
}
plot(as.cimg(get_state(1)))

for (i in 1:100) {
  B[i] = cumsum(svd_decomp$d[1:i])/ sum(svd_decomp$d[1:i])
}

par(mfrow = c(2, 1), mar = rep(1, 4))