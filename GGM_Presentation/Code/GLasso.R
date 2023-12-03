# Set-Up ####
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  library(pacman)
  p_load(glasso, huge 
         , dplyr, ggplot2, tidyr)

# General Function #### 
  glasso.perf = function(N, d, L
                         , .cor.str = "band"
                         , .pkg = "huge"
                         , .band = 1
                         , rho = 0.2 # default values to recreate SLS Fig 9.5
                         , return.obj = F
                         , return.mat = F
                         , .edge.p
                         , .edge.thresh = 0
                         , exact = T){
    if(missing(L)) L = 2*sqrt(log(d)/N)
      
    # if(missing(.cor.str)) .cor.str = "random"; #warning("No .cor.str provided, setting cor structure to banded/AR(1)")
  
    if(tolower(.cor.str) %in% c("band", "banded", "AR", "auto", "autoregressive")){
      
      # Theta = matrix(rep(NA, d*d), nrow=d)
      # for(i in 1:nrow(Theta)){
      #   for(j in 1:ncol(Theta)){
      #     # Theta[i,j] = ifelse(abs(i-j)<=.band, rho, 0) # banded/AR(.band)
      #     Theta[i,j] = ifelse(abs(i-j)<=.band, rho^abs(i-j), 0) # banded/AR(.band)
      #   }
      # }
      
      Theta = huge.generator(n = N, d = d, g = .band, v=rho, graph = "band"
      , verbose = F)$omega
      
    }else if(tolower(.cor.str) %in% c("unstructured", "random", "constructed")){ # not ER, my own 
      
      # Theta = matrix(rep(1, d*d), nrow=d)
      # for(i in 1:nrow(Theta)){
      #     for(j in i:ncol(Theta)){
      #       Theta[i,j] = rbinom(1, 1, prob=.edge.p)*rho
      #       Theta[j,i] = Theta[i,j] # enforce symmetry 
      #     }
      #   }
      
      Theta = huge.generator(n=N, d=d
                             , graph="random"
                             , v = rho 
                             , prob = .edge.p
                             , verbose=F
                             )$omega 
      
    }else if(tolower(.cor.str) %in% c("hub", "star")){
      Theta = huge.generator(n = N, d = d, v=rho, graph = "hub"
                             , verbose = F)$omega
    }
    
    
    # diag(Theta) = rep(1, d) # cleaning/check 
    
    S = solve(Theta)
    if(.pkg == "glasso"){
      g = glasso(S, rho=L, approx = T)
      wi.thr = g$wi
    }else{
      if(exact == T){
        gh = huge::huge.glasso(S, lambda = L, nlambda = 1, verbose = F)
        wi.thr = gh[["icov"]][[1]]
      }else{
        gh = huge(x = S, lambda = L, nlambda = 1, verbose=F, method="mb")
        # gh = huge(x = S, nlambda = 1, verbose=F, method="mb")
        wi.thr = as.matrix(gh$path[[1]])
        }
    }
    
    wi.thr[abs(wi.thr)<.edge.thresh] = 0

    .op = norm(wi.thr - Theta, "2")
    .frob = norm(wi.thr - Theta, "F")
    .acc = sum(wi.thr>1e-10 & abs(Theta)>1e-10) / sum(abs(Theta)>1e-10)
    .neg.acc = sum(abs(wi.thr)<1e-10 & abs(Theta)<1e-10) / sum(abs(Theta)<1e-10)
    # accuracy, diagonal corrected (no self edgess)
    wi.thr.off = wi.thr[row(wi.thr)!=col(wi.thr)]
    Theta.off = Theta[row(Theta)!=col(Theta)]
    
    .acc = sum(wi.thr.off>1e-10 & abs(Theta.off)>1e-10) / sum(abs(Theta.off)>1e-10)
    .neg.acc = sum(abs(wi.thr.off)<1e-10 & abs(Theta.off)<1e-10) / sum(abs(Theta.off)<1e-10)
    
    .output = tibble::lst(.op
                          , .frob
                          , .acc 
                          , .neg.acc
                          , .df.op = data.frame(N, d, .op, .frob, .acc
                                                , .neg.acc
                                                )
                          )
    .output
    
    if(return.obj) .output = tibble::lst(.output, wi.thr)
    if(return.mat) .output = tibble::lst(.output, Theta)
    
    return(.output)
  }
  
# Banded Results #### 
band.glasso.FixN = function(.band, rho, exact = T){
  n = c(10, 50, 100, 250, 300, 500, 1000, 1500, 2000)#, 5e3, 1e5)
  
  d64b2 = lapply(n, glasso.perf, d = 64, .cor.str="band", .band = .band, rho = rho
                 , exact = exact)
  d100b2 = lapply(n, glasso.perf, d = 128, .cor.str="band", .band = .band, rho = rho
                  , exact = exact)
  d225b2 = lapply(n, glasso.perf, d = 256, .cor.str="band", .band = .band, rho = rho
                  , exact = exact)
  
  .df = bind_rows(lapply(d64b2, "[[", ".df.op")
            , lapply(d100b2, "[[", ".df.op")
            , lapply(d225b2, "[[", ".df.op")
  ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc, .neg.acc)
                 , names_to = "metric") %>% 
    filter(metric %in% c(".acc", ".neg.acc", ".op"))
    
  .plot = .df %>% 
    # mutate(tmp = N / (4*log(d))) %>% 
    # ggplot(aes(x=tmp, y=value, color=as.factor(d))) + 
    ggplot(aes(x=N, y=value, color=as.factor(d))) +
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_wrap(~ metric
               , scale = "free_y"
               , labeller = as_labeller(c(.acc = "TPR"
                                          , .neg.acc = "TNR"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
    ylab("Value") + xlab("N")
  
  return(list(.df, .plot))
}

.r = 0.4
set.seed(245)
b1 = band.glasso.FixN(1, rho = .r, exact = T) 
b1mb = band.glasso.FixN(1, rho = .r, exact = F) 
b2 = band.glasso.FixN(2, rho = .r, exact = T) 
b2mb = band.glasso.FixN(2, rho = .r, exact = F) 
b4 = band.glasso.FixN(4, rho = .r, exact = T) 
b4mb = band.glasso.FixN(4, rho = .r, exact = F) 
b8 = band.glasso.FixN(8, rho = .r, exact = T) 
b8mb = band.glasso.FixN(8, rho = .r, exact = F) 


ggsave("glasso_complete_fixN_b1.png", plot = b1[[2]], scale = 0.6)
ggsave("glasso_complete_fixN_b2.png", plot = b2[[2]], scale = 0.6)
ggsave("glasso_complete_fixN_b4.png", plot = b4[[2]], scale = 0.6)
ggsave("glasso_complete_fixN_b8.png", plot = b8[[2]], scale = 0.6)

ggsave("glasso_complete_fixN_b1mb.png", plot = b1mb[[2]], scale=0.6)
ggsave("glasso_complete_fixN_b2mb.png", plot = b2mb[[2]], scale = 0.6)
ggsave("glasso_complete_fixN_b4mb.png", plot = b4mb[[2]], scale = 0.6)
ggsave("glasso_complete_fixN_b8mb.png", plot = b8mb[[2]], scale = 0.6)


# ER Results #### 

band.erdos.FixN = function(.edge.p, .band=1, rho=0.3, exact = T){
  n = c(10, 50, 100, 250, 300, 500, 1000, 1500, 2000)#, 5e3, 1e5)
  
  d64b2 = lapply(n, glasso.perf, d = 64
                 , .cor.str="random", .edge.p = .edge.p
                 , rho = rho
                 , exact = exact)
  d100b2 = lapply(n, glasso.perf, d = 128
                  , .cor.str="random", .edge.p = .edge.p
                  , rho = rho
                  , exact = exact)
  d225b2 = lapply(n, glasso.perf, d = 256
                  , .cor.str="random", .edge.p = .edge.p
                  , rho = rho
                  , exact = exact)
  
  .df = bind_rows(lapply(d64b2, "[[", ".df.op")
                  , lapply(d100b2, "[[", ".df.op")
                  , lapply(d225b2, "[[", ".df.op")
  ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc, .neg.acc)
                 , names_to = "metric") %>% 
    filter(metric %in% c(".acc", ".neg.acc", ".op"))
  
  .plot = .df %>% 
    # mutate(tmp = N / (4*log(d))) %>% 
    # ggplot(aes(x=tmp, y=value, color=as.factor(d))) + 
    ggplot(aes(x=N, y=value, color=as.factor(d))) +
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_wrap(~ metric
               , scale = "free_y"
               , labeller = as_labeller(c(.acc = "TPR"
                                          , .neg.acc = "TNR"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
    ylab("Value") + xlab("N")
  
  return(list(.df, .plot))
}

set.seed(245)

er01 = band.erdos.FixN(.edge.p = 0.01, rho=0.8, exact = T) 
er01mb = band.erdos.FixN(.edge.p = 0.01, rho=0.8, exact = F) 

er05 = band.erdos.FixN(.edge.p=0.05, exact = T) 
er05mb = band.erdos.FixN(.edge.p=0.05, exact = F) 

er1 = band.erdos.FixN(.edge.p=0.1, exact = T) 
er1mb = band.erdos.FixN(.edge.p=0.1, exact = F) 

er2 = band.erdos.FixN(.edge.p=0.2, exact = T) 
er2mb = band.erdos.FixN(.edge.p=0.2, exact = F) 



ggsave("glasso_complete_ER_FixN_01.png", er01[[2]], scale=0.6)
ggsave("glasso_complete_ER_FixN_05.png", er05[[2]], scale=0.6)
ggsave("glasso_complete_ER_FixN_1.png", er1[[2]], scale=0.6)
ggsave("glasso_complete_ER_FixN_2.png", er2[[2]], scale=0.6)

ggsave("glasso_complete_ERmb_FixN_01.png", er01mb[[2]], scale=0.6)
ggsave("glasso_complete_ERmb_FixN_05.png", er05mb[[2]], scale=0.6)
ggsave("glasso_complete_ERmb_FixN_1.png", er1mb[[2]], scale=0.6)
ggsave("glasso_complete_ERmb_FixN_2.png", er2mb[[2]], scale=0.6)









#############################
  ## Fix d #### 
  d = c(100, 250, 500, 1000)
  n100 = lapply(d, glasso.perf, N=100, .cor.str="band", .band = 2, rho = 0.6
               , .pkg = "huge")
  n500 = lapply(d, glasso.perf, N=500, .cor.str="band", .band = 2, rho = 0.6
               , .pkg = "huge")

  bind_rows(lapply(n100, "[[", ".df.op")
            , lapply(n500, "[[", ".df.op")) %>% 
    mutate(OpBound = 2*sqrt(log(d)/N)) %>% 
  pivot_longer(cols = c(.op, .frob, .acc, .neg.acc)
               , names_to = "metric") %>% 
    filter(metric!=".frob") %>%
    ggplot(aes(x=d, y=value, color=as.factor(N))) + 
    geom_point() + geom_line() + 
    theme_minimal() + 
    # coord_cartesian(ylim=c(0, 1.3)) + 
    facet_wrap(~ metric
               , scales = "free_y"
               , labeller = as_labeller(c(.acc = "TPR"
                                          , .neg.acc = "TNR"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="N") + 
    ylab("Value") + xlab("d")
  
# non-sparse graph 
# 
# glasso.perf(N=2000, d=20, .cor.str = "random", .edge.p = 0.5
#             , L = 2*sqrt(log(10)/100))
# 
# 
# 
# 
# 
# 
# .edge.p = 0.7
#   Theta = matrix(rep(1, d*d), nrow=d)
#   for(i in 1:nrow(Theta)){
#     for(j in i:ncol(Theta)){
#       Theta[i,j] = rbinom(1, 1, prob=.edge.p)*runif(1, 0.3, 5)
#       Theta[j,i] = Theta[i,j] # enforce symmetry 
#     }
#   }
# 
# 
# diag(Theta) = rep(1, d) # cleaning/check 
# 
# S = solve(Theta)
# gh = huge(x = S, method="glasso")
# gh = huge::huge.glasso(S, lambda = 2*sqrt(log(d)/N), nlambda = 1, verbose = F)
# g = glasso(S, 2*sqrt(log(d)/N))
# 
# wi.thr[abs(wi.thr)<.edge.thresh] = 0
# 
# .op = norm(wi.thr - Theta, "2")
# .frob = norm(wi.thr - Theta, "F")
# .acc = sum(as.logical(wi.thr) == as.logical(Theta)) / length(as.logical(Theta))
# 
# .output = tibble::lst(.op
#                       , .frob
#                       , .acc 
#                       , .df.op = data.frame(N, d, .op, .frob, .acc)
# )
# .output