
# Set-up ##### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(pacman)
p_load(glasso, huge 
       , dplyr, ggplot2, tidyr
       , cglasso )

head.matrix = function(mat, n=5){
  print(mat[1:n, 1:n])
}

# MAR #### 
N = 1000 
d = 64
rho = 0.3 

Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probna = 0.05)
out <- cglasso(. ~ ., data = Z)



b0 <- runif(p)
B <- matrix(runif(q * p), nrow = q, ncol = p)
X <- matrix(rnorm(n * q), nrow = n, ncol = q)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))
Z <- rcggm(n = n, b0 = b0, X = X, B = B, Sigma = Sigma, probna = 0.05)
out <- cglasso(. ~ ., data = Z)
out

.Z = Z
.Z$X = NULL 

.test = cglasso(. ~ ., data=.Z)
.rho = mean(.test$rho)
.testY = impute(.test ,"mar", rho.new = .rho)

.test = cglasso(. ~ ., data=Z)

















# function from huge package 
# 
# generate_impute = function (n = 200, d = 50, graph = "random", v = NULL, u = NULL, 
#                             g = NULL, prob = NULL, vis = FALSE, verbose = TRUE) {
#   gcinfo(FALSE)
#   if (verbose) 
#     cat("Generating data from the multivariate normal distribution with the", 
#         graph, "graph structure....")
#   if (is.null(g)) {
#     g = 1
#     if (graph == "hub" || graph == "cluster") {
#       if (d > 40) 
#         g = ceiling(d/20)
#       if (d <= 40) 
#         g = 2
#     }
#   }
#   if (graph == "random") {
#     if (is.null(prob)) 
#       prob = min(1, 3/d)
#     prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 - 
#                                                      0.5 * prob)) * (prob >= 0.5)
#   }
#   if (graph == "cluster") {
#     if (is.null(prob)) {
#       if (d/g > 30) 
#         prob = 0.3
#       if (d/g <= 30) 
#         prob = min(1, 6 * g/d)
#     }
#     prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 - 
#                                                      0.5 * prob)) * (prob >= 0.5)
#   }
#   g.large = d%%g
#   g.small = g - g.large
#   n.small = floor(d/g)
#   n.large = n.small + 1
#   g.list = c(rep(n.small, g.small), rep(n.large, g.large))
#   g.ind = rep(c(1:g), g.list)
#   rm(g.large, g.small, n.small, n.large, g.list)
#   gc()
#   theta = matrix(0, d, d)
#   if (graph == "band") {
#     if (is.null(u)) 
#       u = 0.1
#     if (is.null(v)) 
#       v = 0.3
#     for (i in 1:g) {
#       diag(theta[1:(d - i), (1 + i):d]) = 1
#       diag(theta[(1 + i):d, 1:(d - 1)]) = 1
#     }
#   }
#   if (graph == "cluster") {
#     if (is.null(u)) 
#       u = 0.1
#     if (is.null(v)) 
#       v = 0.3
#     for (i in 1:g) {
#       tmp = which(g.ind == i)
#       tmp2 = matrix(runif(length(tmp)^2, 0, 0.5), length(tmp), 
#                     length(tmp))
#       tmp2 = tmp2 + t(tmp2)
#       theta[tmp, tmp][tmp2 < prob] = 1
#       rm(tmp, tmp2)
#       gc()
#     }
#   }
#   if (graph == "hub") {
#     if (is.null(u)) 
#       u = 0.1
#     if (is.null(v)) 
#       v = 0.3
#     for (i in 1:g) {
#       tmp = which(g.ind == i)
#       theta[tmp[1], tmp] = 1
#       theta[tmp, tmp[1]] = 1
#       rm(tmp)
#       gc()
#     }
#   }
#   if (graph == "random") {
#     if (is.null(u)) 
#       u = 0.1
#     if (is.null(v)) 
#       v = 0.3
#     tmp = matrix(runif(d^2, 0, 0.5), d, d)
#     tmp = tmp + t(tmp)
#     theta[tmp < prob] = 1
#     rm(tmp)
#     gc()
#   }
#   if (graph == "scale-free") {
#     if (is.null(u)) 
#       u = 0.1
#     if (is.null(v)) 
#       v = 0.3
#     out = .Call("_huge_SFGen", 2, d)
#     theta = matrix(as.numeric(out$G), d, d)
#   }
#   diag(theta) = 0
#   omega = theta * v
#   diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
#   sigma = cov2cor(solve(omega))
#   omega = solve(sigma)
#   x = mvrnorm(n, rep(0, d), sigma)
#   sigmahat = cor(x)
#   if (vis == TRUE) {
#     fullfig = par(mfrow = c(2, 2), pty = "s", omi = c(0.3, 
#                                                       0.3, 0.3, 0.3), mai = c(0.3, 0.3, 0.3, 0.3))
#     fullfig[1] = image(theta, col = gray.colors(256), main = "Adjacency Matrix")
#     fullfig[2] = image(sigma, col = gray.colors(256), main = "Covariance Matrix")
#     g = graph.adjacency(theta, mode = "undirected", diag = FALSE)
#     layout.grid = layout.fruchterman.reingold(g)
#     fullfig[3] = plot(g, layout = layout.grid, edge.color = "gray50", 
#                       vertex.color = "red", vertex.size = 3, vertex.label = NA, 
#                       main = "Graph Pattern")
#     fullfig[4] = image(sigmahat, col = gray.colors(256), 
#                        main = "Empirical Matrix")
#     rm(fullfig, g, layout.grid)
#     gc()
#   }
#   if (verbose) 
#     cat("done.\n")
#   rm(vis, verbose)
#   gc()
#   sim = list(data = x, sigma = sigma, sigmahat = sigmahat, 
#              omega = omega, theta = Matrix(theta, sparse = TRUE), 
#              sparsity = sum(theta)/(d * (d - 1)), graph.type = graph)
#   class(sim) = "sim"
#   return(sim)
# }