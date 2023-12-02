# Set-Up ####
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  library(pacman)
  p_load(glasso, huge 
         , dplyr, ggplot2, tidyr)

# General Function #### 
  glasso.perf = function(N, d, L
                         , .cor.str = "band"
                         , .pkg = "glasso"
                         , .band = 1
                         , rho = 0.2 # default values to recreate SLS Fig 9.5
                         , return.obj = F
                         , return.mat = F
                         , .edge.p
                         , .edge.thresh = 0){
    if(missing(L)) L = 2*sqrt(log(d)/N)
      
    # if(missing(.cor.str)) .cor.str = "random"; #warning("No .cor.str provided, setting cor structure to banded/AR(1)")
  
    if(tolower(.cor.str) %in% c("band", "banded", "AR", "auto", "autoregressive")){
      
      Theta = matrix(rep(NA, d*d), nrow=d)
      for(i in 1:nrow(Theta)){
        for(j in 1:ncol(Theta)){
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
      g = glasso(S, rho=L)
      wi.thr = g$wi
    }else{
      gh = huge::huge.glasso(S, lambda = L, nlambda = 1, verbose = F)
      # gh = huge::huge.glasso(S, verbose=F)
      wi.thr = gh[["icov"]][[1]]
    }
    
    wi.thr[abs(wi.thr)<.edge.thresh] = 0
    
    .op = norm(wi.thr - Theta, "2")
    .frob = norm(wi.thr - Theta, "F")
    .acc = sum(as.logical(wi.thr) == as.logical(Theta)) / length(as.logical(Theta))

    .output = tibble::lst(.op
                          , .frob
                          , .acc 
                          , .df.op = data.frame(N, d, .op, .frob, .acc)
                          )
    .output
    
    if(return.obj) .output = tibble::lst(.output, g)
    if(return.mat) .output = tibble::lst(.output, Theta)
    
    return(.output)
  }

# Banded Structure (Acc, Op) #### 
  n = seq(1e3, 1e4, length.out=10)
  
  d64 = lapply(n, glasso.perf, d = 64, .cor.str="band", .band = 2, rho = 0.3)
  d100 = lapply(n, glasso.perf, d = 128, .cor.str="band", .band = 2, rho = 0.3)
  d225 = lapply(n, glasso.perf, d = 256, .cor.str="band", .band = 2, rho = 0.3)
  
  
  bind_rows(lapply(d64, "[[", ".df.op")
            , lapply(d100, "[[", ".df.op")
            , lapply(d225, "[[", ".df.op")
            ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc)
                 , names_to = "metric") %>% 
    filter(metric!=".frob") %>%
    ggplot(aes(x=N, y=value, color=as.factor(d))) + 
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_grid(cols = vars(metric)
               , labeller = as_labeller(c(.acc = "Accuracy"
                                 , .op = "Operator 2-Norm"))
               ) + 
    scale_color_brewer(type = "qual", palette = 2
                       , name="d") + 
    ylab("Value") + xlab("N")
    


  
  d64 = lapply(n, glasso.perf, d = 64, .cor.str="band", .band = 8, rho = 0.3
               , .pkg = "huge")
  d100 = lapply(n, glasso.perf, d = 128, .cor.str="band", .band = 8, rho = 0.3)
  d225 = lapply(n, glasso.perf, d = 256, .cor.str="band", .band = 8, rho = 0.3)
  
  
  bind_rows(lapply(d64, "[[", ".df.op")
            , lapply(d100, "[[", ".df.op")
            , lapply(d225, "[[", ".df.op")
  ) %>% 
    pivot_longer(cols = c(.op, .frob, .acc)
                 , names_to = "metric") %>% 
    filter(metric!=".frob") %>%
    ggplot(aes(x=N, y=value, color=as.factor(d))) + 
    geom_point() + geom_line() + 
    theme_minimal() + 
    facet_grid(cols = vars(metric)
               , labeller = as_labeller(c(.acc = "Accuracy"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "seq", palette = "BuGn"
                       , name="d") + 
    ylab("Value") + xlab("N")


  d = c(100, 250, 500, 1000)
  n50 = lapply(d, glasso.perf, N=50, .cor.str="band", .band = 8, rho = 0.3
               , .pkg = "huge")
  n100 = lapply(d, glasso.perf, N=100, .cor.str="band", .band = 8, rho = 0.3
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
    facet_grid(cols = vars(metric)
               , labeller = as_labeller(c(.acc = "Accuracy"
                                          , .op = "Operator 2-Norm"))
    ) + 
    scale_color_brewer(type = "seq", palette = "BuGn"
                       , name="d") + 
    ylab("Value") + xlab("d")
  
    np.test = glasso.perf(N=20, d=1000, .cor.str="band", .band=2, rho=0.4, .pkg="huge")

# non-sparse graph 

glasso.perf(N=2000, d=20, .cor.str = "random", .edge.p = 0.5
            , L = 2*sqrt(log(10)/100))






.edge.p = 0.7
  Theta = matrix(rep(1, d*d), nrow=d)
  for(i in 1:nrow(Theta)){
    for(j in i:ncol(Theta)){
      Theta[i,j] = rbinom(1, 1, prob=.edge.p)*runif(1, 0.3, 5)
      Theta[j,i] = Theta[i,j] # enforce symmetry 
    }
  }


diag(Theta) = rep(1, d) # cleaning/check 

S = solve(Theta)
gh = huge(x = S, method="glasso")
gh = huge::huge.glasso(S, lambda = 2*sqrt(log(d)/N), nlambda = 1, verbose = F)
g = glasso(S, 2*sqrt(log(d)/N))

wi.thr[abs(wi.thr)<.edge.thresh] = 0

.op = norm(wi.thr - Theta, "2")
.frob = norm(wi.thr - Theta, "F")
.acc = sum(as.logical(wi.thr) == as.logical(Theta)) / length(as.logical(Theta))

.output = tibble::lst(.op
                      , .frob
                      , .acc 
                      , .df.op = data.frame(N, d, .op, .frob, .acc)
)
.output