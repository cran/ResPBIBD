
#' Affine resolvable PBIB designs
#'
#' @param v Total number of treatments
#' @param s Positive integer (>=2)
#' @param p Positive integer (>=2)
#' @param q Positive integer (>=2)
#'
#' @return This function generates a PBIB design and its parameters, variance factors and efficiency factor.

#'@description This function generates a new series of affine resolvable Partially Balanced Incomplete Block Designs (PBIBDs) and its parameters (v, b1, b2, r, k1, k2), canonical efficiency factors, variance factor between associates and average variance factors.
#' @examples
#' library(ResPBIBD)
#' PBIBD2(12, 2, 3, 2)
#' @export


PBIBD2=function(v, s, p, q){
if(v==s*p*q && p>=2 && q>=2 && s>=2 && p>=q){
b1=p
b2=q
r=2
k1=s*q
k2=s*p
a1=matrix(1:v, nrow = b1, ncol = k1, byrow = T)
a2=matrix(1:v, nrow = b2, ncol = k2, byrow = F)
a3=matrix(0, nrow = b1, ncol = k2-k1)
a=rbind(cbind(a1,a3),a2)
##################################################
row.names(a) = c(1:nrow(a))
ss = 1
while (ss <= nrow(a)) {
  rownames(a)[ss] <- paste("Block", as.character(ss), sep = "")
  ss = ss + 1
}
################################################
design=rbind(cbind(a1,a3),a2)
k=ncol(design)
b = b1+b2
N_mat = matrix(0, v, b)
for (i in 1:b) {
  for (j in 1:k) {
    N_mat[design[i, j], i] = N_mat[design[i, j], i] + 1
  }
}
K = diag(colSums(N_mat), b, b)
R = diag(rowSums(N_mat), v, v)
kvec = colSums(N_mat)
Kinv = diag(1/kvec, nrow = b, ncol = b)
C_mat = R - N_mat %*% Kinv %*% t(N_mat)
E = eigen(C_mat, only.values = T)
E1 = unlist(E)
E_positive = E1[E1 >= 1e-09]
n = length(E_positive)
C_Efficiency = n/(r * sum(c(1/E_positive)))
p_matrix = matrix(, nrow = 0, ncol = v)
i = 1
j = 1
while (i <= (choose(v, 2))) {
  j = i + 1
  while (j <= v) {
    p1 = matrix(0, nrow = 1, ncol = v)
    p1[i] = 1
    p1[j] = -1
    p_matrix = rbind(p_matrix, p1)
    j = j + 1
  }
  i = i + 1
}
p_invC_Pprme = (p_matrix) %*% MASS::ginv(C_mat) %*% t(p_matrix)
var = diag(p_invC_Pprme)
var1 = round(var, digits = 4)
var2 = unique(var1)
Average_var = mean(var)
A1 = c("Number of treatments", "First set of blocks", "Second set of blocks",
       "Number of replications", "Size b1 blocks","Size of b2 blocks" )
A2 = c("v", "b1","b2", "r", "k1", "k2")
A3 = c(v, b1, b2, r, k1, k2)
A = cbind(A1, A2, A3)
###########################################################
message("\n", "Design parameters")
prmatrix(A, rowlab = , collab = rep("", ncol(A)),
         quote = FALSE, na.print = "")

message("\n"," Affine Resolvable PBIB Design")
a[a == 0] <- NA
prmatrix(a, rowlab = , collab = rep("", ncol(a)),
         quote = FALSE, na.print = "")

message("\n", "Canonical efficiency factor")
print(round(C_Efficiency,4), quote = F)

message("\n", "Variance factor between associates" )
B1 <- c("variance factor between first associates",
        "variance factor between second associates",
        "variance factor between third associates",
        "variance factor between fourth associates")
B2 <- c(var2[1], var2[2], var2[3], var2[4])
if (p == q) {
  B1 <- B1[1:3]
  B2 <- B2[1:3]
}
B <- cbind(B1, B2)
prmatrix(B, rowlab = , collab = rep("", ncol(B)), quote = FALSE, na.print = "")

message("\n", "Average variance factor" )
print(round(Average_var,4), quote = F)
}
else {
  message("Please enter v = s*p*q, where s,p,q >=2 and p>=q")
}
}

###########################################################





