# Generating three types of graphs 
  # ER 
  # 


# adapting code from GI-JOE repo 

p = 4

precision_generation <- function(p
                                 , edge_val_range,eigmin,g_type,target_nodes = NULL, signal = NULL
                                 , n.blocks, ...){
  ## generate random p by p precision matrix for a given graph type, with specified seed
  ## edge_val_range: for a given edge, generate precision entry from a uniform distribution in edge_val_range
  ## eigmin: minimum eigenvalue of the generated precision martix
  ## g_type: graph type
  ## target_nodes and signal: if provided, set the precision entry corresponding to the target node pair as the value of signal
  ## additional input are graph parameters like number of stars in a multi-star graph, connected probability in a random graph
  ## return: precision matrix and adjacency matrix
  input <- list(...)
  Adj_mat <- matrix(rep(0,p^2),p,p); Theta <- matrix(rep(0,p^2),p,p)
  # GTYPE BANDED
  if(g_type == "banded"){
    for(i in 1 : p){
      nb <- i + c(-3 : -1, 1 : 3); nb <- nb[nb >= 1 & nb <= p]
      Adj_mat[i, nb] <- 1
      Theta[i, nb] <- 0.8^(abs(nb - i))
    }
  # All other G types 
  else if(g_type == "block"){
      K = n.blocks
      for(i in 1:(p/K)){
        Adj_mat[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)] <- 1
      }
      diag(Adj_mat) <- 0
    }
    
    # G TYPE Circle cleanup 
    # for(i in 1:(p-1)){
    #   for(j in (i+1):p){
    #     if(Adj_mat[i,j]==1){
    #       Theta[i,j] <- runif(1,edge_val_range[1],edge_val_range[2])
    #     }
    #   }
    # }
    # Theta <- Theta + t(Theta)
  
    
  if(!is.null(target_nodes) & !is.null(signal)){
    Theta[target_nodes[1], target_nodes[2]] <- signal
    Theta[target_nodes[2], target_nodes[1]] <- signal
  }
    
  Theta <- Theta + (eigmin-min(eigen(Theta)[[1]]))*diag(p)
  return(list(Theta,Adj_mat))
  }
  
#----------------------------------------------------#
data_generation_pairwise <- function(Sigma, N, seed){
  ## generate pairwise measured data from from N(0,Sigma), with specified seed 
  ## N: p by p pairwise sample size matrix specified by the user; only upper triangular elements are used
  ## return: data matrix "X", observation pattern "Obs", unbiased entry-wise estimate of covariance "Sigma_hat", and a list of list "Ind" that specifies the row indices of X corresponding to each pair (i,j)
  p <- dim(Sigma)[1]
  ind <- upper.tri(N, diag = FALSE)
  total_n <- sum(N[ind]);
  set.seed(seed)
  X <- matrix(rep(0, total_n * p), total_n, p);
  Obs <- matrix(rep(0, total_n * p), total_n, p);
  start_n <- 1;
  diag_n <- rep(0, p); diag_Sigma <- rep(0, p);
  Sigma_hat <- matrix(rep(0, p^2), p, p)
  Ind <- list(); 
  for(i in 1:(p-1)){
    Ind[[i]] <- list(); 
    for(j in (i+1):p){
      end_n <- start_n + N[i,j] - 1;
      Obs[start_n:end_n, c(i,j)] <- 1;
      sqrt_Sigma_tmp <- msqrt(Sigma[c(i,j), c(i,j)])[[1]];
      X[start_n : end_n, c(i, j)] <- matrix(rnorm(N[i, j] * 2, 0, 1), N[i, j], 2) %*% sqrt_Sigma_tmp; 
      Ind[[i]][[j]] <- start_n : end_n #Ind only includes the sample sizes for i < j
      Sigma_hat[i,j] <- sum(X[start_n : end_n, i] * X[start_n : end_n, j])/N[i, j];
      diag_Sigma[i] <- diag_Sigma[i] + sum(X[start_n : end_n, i]^2);
      diag_n[i] <- diag_n[i] + N[i, j];
      diag_Sigma[j] <- diag_Sigma[j] + sum(X[start_n : end_n,j]^2);
      diag_n[j] <- diag_n[j] + N[i, j];
      start_n <- end_n + 1;
    }
  }
  Sigma_hat <- Sigma_hat + t(Sigma_hat);
  diag(Sigma_hat) <- diag_Sigma / diag_n; 
  X <- X * Obs;
  return(list(X = X, Obs = Obs, Sigma_hat = Sigma_hat, Ind = Ind))
}

generate_Obs <- function(p, n_total, Adj_mat, measurement_scenario, Obs_seed){
  ## generate observational patterns used in the comparative studies of the paper
  ## p: number of nodes
  ## n_total: total number of incomplete samples
  ## Adj_mat: adjacency matrix of the considered graph
  ## measurement_scenario: a number from 1 to 5, indication five different erose measurement types
  ## Obs_seed: seed for generating random observational patterns
  set.seed(Obs_seed)
  if(measurement_scenario == 1){
    #missing independently with small, medium, large sample sizes
    prob <- sqrt(c(rep(0.1, floor(p / 3)), rep(0.5, floor(p / 3)), rep(0.9, p - floor(p / 3) * 2)))
    Obs <- matrix(rep(0, p * n_total), n_total, p)
    for(j in 1 : p){
      Obs[, j] <- as.numeric(runif(n_total, 0, 1) <= prob[j])
    }
  }else if(measurement_scenario == 2){
    #missing independently, missing probability exponentially depending on degrees
    degree <- apply(Adj_mat, 1, sum)
    prob <- 1 - 0.95 ^ (4 * degree)
    Obs <- matrix(rep(0, p * n_total), n_total, p)
    for(j in 1 : p){
      Obs[, j] <- as.numeric(runif(n_total, 0, 1) <= prob[j])
    }
  }else if(measurement_scenario == 3){
    #size-constrained measurements; each time randomly sample a subset, 
    #each node being picked with probability depending on degree
    degree <- apply(Adj_mat, 1, sum)
    prob <- 1 - 0.95 ^ (4 * degree)
    Obs <- matrix(rep(0, p * n_total), n_total, p)
    for(i in 1 : n_total){
      Obs[i, sample(p, size = 20, replace = FALSE, prob = prob)] <- 1
    }
  }else if(measurement_scenario == 4){
    #size-constrained measurements; each time first randomly sample one center, 
    #then sample the rest of the nodes with probability weights depending on their distances to the selected center
    Obs <- matrix(rep(0, p * n_total), n_total, p)
    for(i in 1 : n_total){
      center <- sample(p, size = 1);
      prob <- 0.95 ^ (abs((1 : p)[-center] - center))
      Obs[i, c(center, sample((1 : p)[-center], size = 19, replace = FALSE, prob = prob))] <- 1
    }
  }else if(measurement_scenario == 5){
    #graph-quilting setting
    Obs <- matrix(rep(0, p * n_total), n_total, p)
    overlap <- 10; nblocks <- 10
    blocksize <- floor((p + (nblocks - 1) * overlap) / nblocks)
    nsample <- floor(n_total / nblocks)
    for(i in 1 : nblocks){
      if(i < nblocks){
        Obs[(nsample * (i - 1) + 1) : (nsample * i), ((blocksize - overlap) * (i - 1) + 1) : ((blocksize - overlap) * (i - 1) + blocksize)] <- 1
      }else{
        Obs[(nsample * (i - 1) + 1) : n_total, ((blocksize - overlap) * (i - 1) + 1) : p] <- 1
      }
    }
  }
  return(Obs)
}