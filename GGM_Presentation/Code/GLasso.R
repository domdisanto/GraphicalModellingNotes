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
                         , exact = F){
    if(missing(L)) L = 2*sqrt(log(d)/N)
      
    # if(missing(.cor.str)) .cor.str = "random"; #warning("No .cor.str provided, setting cor structure to banded/AR(1)")
  
    if(tolower(.cor.str) %in% c("band", "banded", "AR", "auto", "autoregressive")){
      
      Theta = matrix(rep(NA, d*d), nrow=d)
      for(i in 1:nrow(Theta)){
        for(j in 1:ncol(Theta)){
          # Theta[i,j] = ifelse(abs(i-j)<=.band, rho, 0) # banded/AR(.band)
          Theta[i,j] = ifelse(abs(i-j)<=.band, rho^abs(i-j), 0) # banded/AR(.band)
        }
      }
      
    }else if(tolower(.cor.str) %in% c("unstructured", "random", "constructed")){ # not ER, my own 
      
      Theta = matrix(rep(1, d*d), nrow=d)
      for(i in 1:nrow(Theta)){
          for(j in i:ncol(Theta)){
            Theta[i,j] = rbinom(1, 1, prob=.edge.p)*runif(1, 0.3, 0.6)
            Theta[j,i] = Theta[i,j] # enforce symmetry 
          }
        }
      }
    
    
    diag(Theta) = rep(1, d) # cleaning/check 
    
    S = solve(Theta)
    if(.pkg == "glasso"){
      g = glasso(S, rho=L, approx = T)
      wi.thr = g$wi
    }else{
      if(exact == T){
        gh = huge::huge.glasso(S, lambda = L, nlambda = 1, verbose = F)
        # gh = huge::huge.glasso(S, nlambda = 1, verbose = F)
        wi.thr = gh[["icov"]][[1]]
      }else{
        gh = huge(x = S, lambda = L, nlambda = 1, verbose=F, method="mb")
        wi.thr = as.matrix(gh$path[[1]])
        }
    }
    
    wi.thr[abs(wi.thr)<.edge.thresh] = 0
    
    .op = norm(wi.thr - Theta, "2")
    .frob = norm(wi.thr - Theta, "F")
    .acc = sum(wi.thr>0 & Theta>0) / sum(Theta>0)
    # .acc = sum(as.logical(wi.thr) == as.logical(Theta)) / length(as.logical(Theta))
    # .edge.recovery = sum(as.logical(wi.thr) ) / sum(as.logical(Theta))
    
    .output = tibble::lst(.op
                          , .frob
                          , .acc 
                          # , .edge.recovery
                          , .df.op = data.frame(N, d, .op, .frob, .acc
                                                # , .edge.recovery
                                                )
                          )
    .output
    
    if(return.obj) .output = tibble::lst(.output, wi.thr)
    if(return.mat) .output = tibble::lst(.output, Theta)
    
    return(.output)
  }

  
  .test = glasso.perf(N=50, d=64, .cor.str="band", .pkg="huge"
              , .band = 1
              , rho = 0.2
              , return.mat = F
              , return.obj = F
              , exact = T)

  .test
  
# Banded Structure (Acc, Op) #### 
  ## Fix n #### 
  ### Banded(2) #### 
  # n = seq(1e2, 1e3, length.out=10)
  n = c(10, 50, 100, 250, 300, 500, 1e3)
  
  d64b2 = lapply(n, glasso.perf, d = 64, .cor.str="band", .band = 1, rho = 0.2)
  d100b2 = lapply(n, glasso.perf, d = 128, .cor.str="band", .band = 1, rho = 0.2)
  d225b2 = lapply(n, glasso.perf, d = 256, .cor.str="band", .band = 1, rho = 0.2)
  
  
  b3.plot = bind_rows(lapply(d64b2, "[[", ".df.op")
            , lapply(d100b2, "[[", ".df.op")
            , lapply(d225b2, "[[", ".df.op")
            ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc)
                 , names_to = "metric") %>% 
    filter(metric %in% c(".acc", ".op")) %>%
    # mutate(tmp = N / (4*log(d))) %>% 
    # ggplot(aes(x=tmp, y=value, color=as.factor(d))) + 
    ggplot(aes(x=N, y=value, color=as.factor(d))) +
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_wrap(~ metric
               , scale = "free_y"
               , labeller = as_labeller(c(.acc = "Proportion of Edges Selected"
                                 , .op = "Operator 2-Norm"))
               ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
    ylab("Value") + xlab("N")
  
  
  ggsave("glasso_complete_fixN_b3.png"
         , plot = b3.plot
         , scale = 0.6)


  ### Banded(2) #### 
  # n = seq(1e3, 1e4, length.out=10)
  n = c(10, 50, 100, 250, 500, 1e3)
  
  d64b8 = lapply(n, glasso.perf, d = 64, .cor.str="band", .band = 8, rho = 0.3
               , .pkg = "huge")
  d100b8 = lapply(n, glasso.perf, d = 128, .cor.str="band", .band = 8, rho = 0.3)
  d225b8 = lapply(n, glasso.perf, d = 256, .cor.str="band", .band = 8, rho = 0.3)
  
  
  b8.plot = bind_rows(lapply(d64b8, "[[", ".df.op")
            , lapply(d100b8, "[[", ".df.op")
            , lapply(d225b8, "[[", ".df.op")
  ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc)
                 , names_to = "metric") %>% 
    filter(metric!=".frob") %>%
    ggplot(aes(x=N, y=value, color=as.factor(d))) + 
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_wrap(~ metric
               , scales = "free_y"
               , labeller = as_labeller(c(.acc = "Accuracy"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
    ylab("Value") + xlab("N")

  ggsave("glasso_complete_fixN_b8.png"
         , plot = b8.plot
         , scale = 0.6)

  ## Fix d #### 
  d = c(100, 250, 500, 1000, 5000)
  n50 = lapply(d, glasso.perf, N=50, .cor.str="band", .band = 3, rho = 0.3
               , .pkg = "huge")
  n100 = lapply(d, glasso.perf, N=100, .cor.str="band", .band = 3, rho = 0.3
               , .pkg = "huge")

  bind_rows(lapply(n50, "[[", ".df.op")
            , lapply(n100, "[[", ".df.op")) %>% 
    mutate(OpBound = 2*sqrt(log(d)/N)) %>% 
  pivot_longer(cols = c(.op, .frob, .acc)
               , names_to = "metric") %>% 
    filter(metric!=".frob") %>%
    ggplot(aes(x=d, y=value, color=as.factor(N))) + 
    geom_point() + geom_line() + 
    theme_minimal() + 
    coord_cartesian(ylim=c(0, 1.3)) + 
    facet_grid(cols = vars(metric)
               , scales = "free"
               , labeller = as_labeller(c(.acc = "Accuracy"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
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