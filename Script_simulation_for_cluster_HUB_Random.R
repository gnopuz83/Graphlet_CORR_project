library(MASS)
library(huge)
library(Matrix)
library(parallel)
library(orca)
library(batchtools)
library(matrixcalc)

Hub_graph = function(d, s = 20){
  
  sub_hub = matrix(data = rep(0, s*s), nrow = s, ncol = s)
  sub_hub[1,] = 1
  sub_hub[,1] = 1
  
  number_of_hubs = round(d/s)
  my_matrix = bdiag(replicate(number_of_hubs, sub_hub, simplify = FALSE))
  return(as.matrix(my_matrix))
}


simulate_HUB_graph = function(n, d ,rho_hub = 0.22, s = 20, seed = 123, K = 100, scale_cov = TRUE){
  ####### create the random graph
  set.seed(seed)
  #ro = 2/(s-1)
  ro = rho_hub
  Random_hub = Hub_graph(d = d, s = s)
  
  
  ##### set the values in the inverse as a uniform with values -1, 1
  I_C = Random_hub*ro
  diag(I_C) = 1
  #B = seq(1000,0.1, -0.0001)
  #BS = binary_search(B, length(B), K,function_condition_number, I_C)
  #diag(I_C) = B[BS]
  if(scale_cov==TRUE){
    S = cov2cor(solve(I_C))
  }else{
    S = solve(I_C)
  }
  S2 = t(S)
  S[upper.tri(S)] = S2[upper.tri(S2)] 
  #print(sum(S==t(S)))
  #S = solve(I_C)
  S1 = solve(S)
  x = mvrnorm(n, rep(0, d), S)
  print(dim(I_C))
  condition_number = kappa(I_C)
  res = list("simulation_data" = x, "real_graph" = Random_hub, "real_inverse" = I_C, "inverse_of_cov" = S1, "condition_number" = condition_number)
  return(res)
}


GCD_oracle_random = function(True_graph, pulsar_obj, GC_fun, GCF = F){
  
  path_1 = c()
  
  #interval = pulsar_obj$stars$lb.index - pulsar_obj$stars$ub.index
  
  G1 = True_graph
  diag(G1) = 0
  if(GCF ==T){
    GCM1 = GC_fun(G1)
  }else{
    GCM1 = gcvec_extended_random3(G1)
  }
  
  opt_lambda = pulsar_obj$stars$opt.index
  
  for(i in 1:length(pulsar_obj$est$lambda)){
    pulsar_obj$stars$opt.index = i
    res1 = refit(pulsar_obj)
    GCM2 = GC_fun(res1$refit$stars)
    dist_gcd = GCD(GCM2, GCM1)
    #print(dist_gcd)
    path_1 = c(dist_gcd, path_1)
  }
  pulsar_obj$stars$opt.index = opt_lambda
  return(path_1)
}

GCD_oracle_random_nothing = function(True_graph, pulsar_obj, GC_fun, GCF = F){
  
  path_1 = c()
  
  #interval = pulsar_obj$stars$lb.index - pulsar_obj$stars$ub.index
  
  G1 = True_graph
  diag(G1) = 0
  if(GCF ==T){
    GCM1 = GC_fun(G1)
  }else{
    GCM1 = gcvec_extended_random1(G1)
  }
  
  opt_lambda = pulsar_obj$stars$opt.index
  
  for(i in 1:length(pulsar_obj$est$lambda)){
    pulsar_obj$stars$opt.index = i
    res1 = refit(pulsar_obj)
    GCM2 = GC_fun(res1$refit$stars)
    dist_gcd = GCD(GCM2, GCM1)
    #print(dist_gcd)
    path_1 = c(dist_gcd, path_1)
  }
  pulsar_obj$stars$opt.index = opt_lambda
  return(path_1)
}



Hamming_distance = function(X,Y){
  
  X1 = X[upper.tri(X)]
  Y1 = Y[upper.tri(Y)]
  
  X2 = as.vector(X1)
  Y2 = as.vector(Y1)
  
  Res = sum(X2!=Y2)
  
  return(Res)
}

Hamming_distance_path = function(pulsar_obj, true_graph){
  
  path_1 = c()
  
  #interval = pulsar_obj$stars$lb.index - pulsar_obj$stars$ub.index
  
  #print(pulsar_obj$stars$ub.index)
  #print(pulsar_obj$stars$lb.index)
  
  
  opt_lambda = pulsar_obj$stars$opt.index
  for(i in 1:length(pulsar_obj$est$lambda)){
    pulsar_obj$stars$opt.index = i
    res1 = refit(pulsar_obj)
    H_D = Hamming_distance(res1$refit$stars, true_graph)
    #print(H_D)
    path_1 = c(H_D, path_1)
  }
  pulsar_obj$stars$opt.index = opt_lambda
  return(path_1)
}


Function_GCD_data = function(pulsar_obj, real_G){
  mean_gcd1 = c()
  mean_gcd2 = c()
  mean_gcd3 = c()
  mean_gcd4 = c()
  HD_mean = c()
  
  tot_df1 = list()
  K = 1
  for(lambda_sub in rev(pulsar_obj$gcd_new$graph_subsambling)){
    sub_df1 = c()
    sub_df2 = c()
    sub_df3 = c()
    sub_df4 = c()
    
    tot_GCM = list()
    
    HD_df1 = c()
    
    GCM1 = gcvec_extended_random1(real_G)
    GCM12 = gcvec_extended_random1(real_G)
    #GCM13 = gcvec_extended_random2(real_G)
    GCM13 = gcvec_extended_random3(real_G)
    
    C = 1
    
    for(graph1 in lambda_sub){
      
      GCM2 = gcvec_extended_random1(graph1)
      GCM3 = gcvec_extended_random2(graph1)
      GCM4 = gcvec_extended_random3(graph1)
      
      dist_gcd1 = GCD(GCM2, GCM1)
      dist_gcd2 = GCD(GCM2, GCM13)
      dist_gcd3 = GCD(GCM3, GCM1)
      dist_gcd4 = GCD(GCM3, GCM13)
      
      HD1 = Hamming_distance(graph1, real_G)
      #print(dist_gcd)
      sub_df1 = c(sub_df1, dist_gcd1)
      sub_df2 = c(sub_df2, dist_gcd2)
      sub_df3 = c(sub_df3, dist_gcd3)
      sub_df4 = c(sub_df4, dist_gcd4)
      HD_df1 = c(HD_df1, HD1)
      tot_GCM$GCM1[[C]] = GCM2
      tot_GCM$GCM2[[C]] = GCM3
      tot_GCM$GCM3[[C]] = GCM4
      C = C+1
    }
    #tot_df1$GCD1[[K]] = sub_df1
    #tot_df1$GCD2[[K]] = sub_df2
    #tot_df1$GCD3[[K]] = sub_df3
    #tot_df1$GCD4[[K]] = sub_df4
    #tot_df1$GCMs[[K]] = tot_GCM
    
    
    K = K+1
    mean_gcd1 = c(mean_gcd1, mean(sub_df1))
    mean_gcd2 = c(mean_gcd2, mean(sub_df2))
    mean_gcd3 = c(mean_gcd3, mean(sub_df3))
    mean_gcd4 = c(mean_gcd4, mean(sub_df4))
    HD_mean = c(HD_mean, mean(HD_df1))
  }
  
  tot_df1$mean_gcd1 = mean_gcd1
  tot_df1$mean_gcd2 = mean_gcd2
  tot_df1$mean_gcd3 = mean_gcd3
  tot_df1$mean_gcd4 = mean_gcd4
  tot_df1$HD_mean_subsampling = HD_mean
  
  return(tot_df1)
}




binary_search = function(A, n, T1, F1, M, Thresh = 0.1){
  # Binary search algorithm for the condition number
  L = 0
  R = n-1
  while(L<=R){
    m = floor((L+R)/2)
    if(F1(M,A[m])<(T1-Thresh)){
      L = m+1
    }else if(F1(M,A[m])>(T1+Thresh)){
      R = m-1
    }else if(F1(M,A[m])>=(T1-Thresh)&F1(M,A[m])<=(T1+Thresh)){
      return(m)
    }else{
      print("not_found")
      return(FALSE)}
  }
  return(FALSE)
  print("not_found")
}
function_condition_number = function(Matrix, scale_diag){
  diag(Matrix) = scale_diag
  C_N = kappa(Matrix)
  return(C_N)
}

# Function for set the element of the inverse covariance matrix [-1,1]
generate_minus1_one = function(adj_m, values_interval, values){
  V = as.integer(values/2)
  tot_seq1 = runif(V , values_interval[1], values_interval[2])
  if((V+V)!=values){
    V = V+1
  }
  tot_seq2 = runif(V , values_interval[3], values_interval[4])
  tot_seq = c(tot_seq1, tot_seq2)
  sub_tri = adj_m[upper.tri(adj_m)]
  sub_tri[sub_tri!=0] = tot_seq
  adj_m[upper.tri(adj_m)] = sub_tri
  t_ad_2 = t(adj_m)
  t_ad_2[upper.tri(t_ad_2)] = sub_tri
  return(t_ad_2)
}

# Function for simulate the random graph and the data
simulate_random_graph = function(n, d, interval, K, seed = 123, scale_cov = FALSE){
  ####### create the random graph
  set.seed(seed)
  G = Prob(d)
  p_choose2 = d*(d-1)/2
  p = 3/d
  number_of_edges = as.integer(p_choose2*p)
  values_ones = as.integer(sum(G)/2)
  if(values_ones==number_of_edges){
    ##### set the values in the inverse as a uniform with values -1, 1
    I_C = generate_minus1_one(G,interval, number_of_edges)
    diag(I_C) = 1
    while(is.positive.definite(I_C) == FALSE){
      I_C = generate_minus1_one(G,interval, number_of_edges)
      B = seq(2000,0.1, -0.0001)
      BS = binary_search(B, length(B), K,function_condition_number, I_C)
      if(BS==FALSE){
        I_C = I_C
      }else{
        diag(I_C) = B[BS]}
      #print(is.positive.definite(I_C))
    }
    #B = seq(1000,0.1, -0.0001)
    #BS = binary_search(B, length(B), K,function_condition_number, I_C)
    #diag(I_C) = B[BS]
    if(scale_cov == TRUE){
      S = cov2cor(solve(I_C))
    }else{
      S = solve(I_C)
    }
    S1 = solve(S)
    x = mvrnorm(n, rep(0, d), S)
    return(list("simulation_data" = x, "real_graph" = G, "real_inverse" = I_C, "inverse_of_cov" = S1))
  }else{print("something wrong in the number of edges")}
}

# Function for construct the random adjacency matrix
Prob = function(d){
  p_choose2 = d*(d-1)/2
  p = 3/d
  number_of_edges = as.integer(p_choose2*p)
  random_matrix_uniform = runif(d*(d-1)/2, 0,1)
  Matrix_random = matrix(rep(0, d*d), d,d)
  Matrix_random[upper.tri(Matrix_random)] = random_matrix_uniform
  T_MR = t(Matrix_random)
  T_MR[upper.tri(T_MR)] = random_matrix_uniform
  T_MR[T_MR<p] = 1
  T_MR[T_MR<1] = 0
  diag(T_MR) = 0
  L = T_MR
  tot_edge = sum(L)/2
  V = sign(tot_edge-number_of_edges)
  if(V==1){
    remove_n = tot_edge-number_of_edges
    set_to_0 = sample(which(L[upper.tri(L)] == 1), remove_n)
    zeros = as.vector(L[upper.tri(L)])
    zeros[set_to_0] = 0
    L[upper.tri(L)] = zeros
    Matrix_random = matrix(rep(0, d*d), d,d)
    Matrix_random[upper.tri(Matrix_random)] = as.vector(L[upper.tri(L)])
    T_m = t(Matrix_random)
    T_m[upper.tri(T_m)] = as.vector(L[upper.tri(L)])
    L = T_m
    if(sum(L==t(L))==d*d){
      return(L)
    }else{
      print("no simmetric")
    }
  }else if(V==0){
    if(sum(L==t(L))==d*d){
      return(L)
    }else{
      print("no simmetric")
    }
  }else if(V==-1){
    while(V!=0) {
      rows = seq(1,d)
      r1 = sample(rows, 1)
      c1 = sample(rows,1)
      L[r1,c1] = 1
      Matrix_random = matrix(rep(0, d*d), d,d)
      Matrix_random[upper.tri(Matrix_random)] = as.vector(L[upper.tri(L)])
      T_m = t(Matrix_random)
      T_m[upper.tri(T_m)] = as.vector(L[upper.tri(L)])
      L = T_m
      diag(L) = 0
      tot_edge = sum(L)/2
      number_of_edges = as.integer(p_choose2*p)
      V = sign(tot_edge-number_of_edges)
    }
    if(sum(L==t(L))==d*d){
      return(L)
    }else{
      print("no symmetric")
    }
  }else{}
  return(L)
}






gcvec_extended_random1 = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],runif(length(orbind),0, 1)), method='spearman'))
  gcor[upper.tri(gcor)]
}



gcvec_extended_random = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],runif(length(orbind),0, 1)), method='spearman'))
  gcor[upper.tri(gcor)]
}


gcvec_extended_random1 = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],runif(length(orbind),0, 1)), method='spearman'))
  gcor[upper.tri(gcor)]
}



gcvec_extended_random = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  
  after1 = rbind(gcount[,orbind],1)
  after2 = rbind(gcount[,orbind],runif(length(orbind),0,1))
  #after2 = rbind(after2,runif(length(orbind),0,1))
  
  gcor1 <- suppressWarnings(cor(after1, method='spearman'))
  gcor2 = suppressWarnings(cor(after2, method='spearman'))
  gcor = gcor1[upper.tri(gcor1)]
  gcor_random = gcor2[upper.tri(gcor2)]
  
  return(list("corr_1" = gcor1, "res_1" = gcor, "gcount_1" = gcount, "corr_2" = gcor2, "gcor2" = gcor2, "res2" = gcor_random, "after1" = after1, "after2"=after2))
}


gcvec_extended_random3 = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  
  for(i in colnames(gcount)){
    if(sum(gcount[,i])==0){
      random_vector = runif(dim(gcount)[1],0, 0.001)
      gcount[,i] = random_vector
    }
  }
  
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],runif(length(orbind),0, 1)), method='spearman'))
  gcor[upper.tri(gcor)]
}




gcvec_extended_random2 = function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  add_1_minus1 = gcount[,orbind]
  minus_1_1 = matrix(runif(length(orbind)*dim(gcount)[1],0,0.01), nrow = dim(gcount)[1], ncol = length(orbind))
  add_1_minus1 = add_1_minus1 + minus_1_1
  gcor <- suppressWarnings(cor(add_1_minus1, method='spearman'))
  gcor[upper.tri(gcor)]
}












#' pulsar: batch mode
#'
#' Run pulsar using stability selection, or another criteria, to select an undirected graphical model over a lambda-path.
#'
#' @param wkdir set the working directory if different than \code{\link{getwd}}
#' @param regdir directory to store intermediate batch job files. Default will be a tempory directory
#' @param init text string appended to basename of the regdir path to store the batch jobs for the initial StARS variability estimate (ignored if `regdir` is NA)
#' @param conffile path to or string that identifies a \code{\link[batchtools:batchtools-package]{batchtools}} configuration file. This argument is passed directly to the \code{name} argument of the \code{\link[pulsar]{findConfFile}} function. See that help for detailed explanation.
#' @param job.res named list of resources needed for each job (e.g. for PBS submission script). The format and members depends on configuration and template. See examples section for a Torque example
#' @param cleanup Flag for removing batchtools registry files. Recommended FALSE unless you're sure intermediate data shouldn't be saved.
#' @param refit Boolean flag to refit on the full dataset after pulsar is run. (see also \code{\link{refit}})
#' @return an S3 object of class \code{\link{batch.pulsar}} with a named member for each stability criterion/metric. Within each of these are:
#' \itemize{
#'    \item summary: the summary criterion over \code{rep.num} graphs at each value of lambda
#'    \item criterion: the stability metric
#'    \item merge: the raw criterion merged over the \code{rep.num} graphs (constructed from \code{rep.num} subsamples), prior to summarization
#'    \item opt.ind: index (along the path) of optimal lambda selected by the criterion at the desired threshold. Will return \eqn{0} if no optimum is found or \code{NULL} if selection for the criterion is not implemented.
#'   }
#' If \code{stars} is included as a criterion then additional arguments include
#' \itemize{
#'    \item lb.index: the lambda index of the lower bound at \eqn{N=2} samples if \code{lb.stars} flag is set to TRUE
#'    \item ub.index: the lambda index of the upper bound at \eqn{N=2} samples if \code{ub.stars} flag is set to TRUE
#'}
#' @return reg: Registry object. See \code{batchtools::makeRegistry}
#' @return id: Identifier for mapping graph estimation function. See \code{batchtools::batchMap}
#' @return call: the original function call
#' @examples
#' \dontrun{
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 400 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(.2, .01, len=40)
#' hugeargs  <- list(lambda=lams, verbose=FALSE)
#'
#' ## Run batch.pulsar using snow on 5 cores, and show progress.
#' options(mc.cores=5)
#' options(batchtools.progress=TRUE, batchtools.verbose=FALSE)
#' out <- batch.pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                  rep.num=20, criterion='stars', conffile='snow')

#' ## Run batch.pulsar on a Torque cluster
#' ## Give each job 1gb of memory and a limit of 30 minutes
#' resources <- list(mem="1GB", nodes="1", walltime="00:30:00")
#' out.p <- batch.pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                  rep.num=100, criterion=c('stars', 'gcd'), conffile='torque'
#'                  job.res=resources, regdir=file.path(getwd(), "testtorq"))
#' plot(out.p)

#' ## take a look at the default torque config and template files we just used
#' file.show(findConfFile('torque'))
#' file.show(findTemplateFile('simpletorque'))
#' }
#' @references Müller, C. L., Bonneau, R., & Kurtz, Z. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv https://arxiv.org/abs/1605.07072
#' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. Proceedings of the Twenty-Third Annual Conference on Neural Information Processing Systems (NIPS).
#' @references Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. The Journal of Machine Learning Research, 13, 1059–1062.
#' @references Michel Lang, Bernd Bischl, Dirk Surmann (2017). batchtools: Tools for R to work on batch systems. The Journal of Open Source Software, 2(10). URL https://doi.org/10.21105/joss.00135.
#' @inheritParams pulsar
#' @importFrom Matrix mean triu
#' @seealso \code{\link{pulsar}} \code{\link{refit}}
#' @export
batch.pulsar <- function(data, fun=huge::huge, fargs=list(),
                         criterion=c("stars"), thresh = 0.1, subsample.ratio = NULL,
                         lb.stars=FALSE, ub.stars=FALSE, rep.num = 20, seed=NULL,
                         wkdir=getwd(), regdir=NA, init="init", conffile='',
                         job.res=list(), cleanup=FALSE, refit=TRUE, known_graph = FALSE, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE) {
  
  if (!requireNamespace('batchtools', quietly=TRUE)) {
    stop("'batchtools' package required to run 'batch.pulsar'")
  }
  gcinfo(FALSE)
  if (!is.na(regdir))
    if (file.exists(regdir)) stop('Registry directory already exists')
  
  n <- nrow(data)
  p <- ncol(data)
  # min requirements for function args
  knowncrits <- c("stars", "gcd", "estrada", "sufficiency", "gcd_new")
  .lamcheck(fargs$lambda)
  .critcheck0(criterion, knowncrits)
  subsample.ratio <- .ratcheck(subsample.ratio, n)
  nlams    <- length(fargs$lambda)
  conffile <- findConfFile(conffile)
  
  if (!is.null(seed)) set.seed(seed)
  ind.sample <- replicate(rep.num, sample(c(1:n),
                                          floor(n*subsample.ratio), replace=FALSE), simplify=FALSE)
  if (refit) {
    tmp <- 1L:n
    attr(tmp, 'full') <- TRUE
    ind.sample <- c(list(tmp), ind.sample)
  }
  if (!is.null(seed)) set.seed(NULL)
  
  ## build the estimator function that takes the randomized sample index
  ## and the full data
  estFun <- function(ind.sample, fargs, data, fun) {
    tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
    if (!('path' %in% names(tmp)))
      stop('Error: expected data stucture with \'path\' member')
    
    if (isTRUE(attr(ind.sample, 'full')))
      return(tmp)
    else
      return(tmp$path)
  }
  
  
  est <- list()
  reduceargs <- list() ; reduceGCDargs <- list()
  if (lb.stars) {
    if (!("stars" %in% criterion)) # || length(criterion)>1)
      stop('Lower/Upper bound method must be used with StARS')
    minN   <- 2 + refit
    if (!is.na(regdir)) regdir <- paste(regdir, init, sep="_")
  } else
    minN <- rep.num + refit
  
  isamp <- ind.sample[1:minN]
  out <- batchply(data, estFun, fun, fargs, isamp,
                  wkdir, regdir, conffile, job.res)
  reg <- out$reg ; id  <- out$id
  ## jobs w/ no errors
  doneRun <- batchtools::waitForJobs(reg=reg, id)
  jdone   <- batchtools::findDone(reg=reg, id)
  pulsar.jobs <- intersect((1+refit):minN, jdone$job.id)
  
  if (refit) {
    fullmodel <- batchtools::loadResult(id=1, reg=reg)
    minN <- minN - 1L
  } else {
    fullmodel <- NULL
  }
  ## Use this for stars crits
  starsaggfun <- function(res, aggr)
    lapply(1:length(aggr), function(i) aggr[[i]] + res[[i]])
  
  if (lb.stars) {
    est$init.reg <- reg ; est$init.id  <- id
    
    if (!doneRun)
      stop('Errors in batch jobs for computing initial stability')
    
    # collect initial results
    lb.starsmerge <- batchtools::reduceResults(reg=reg, ids=pulsar.jobs, fun=starsaggfun)
    lb.est <- stars.stability(NULL, thresh, minN, p, lb.starsmerge)
    gc(FALSE)
    # compute initial gcd if selected
    if ('gcd' %in% criterion || 'gcd_new' %in% criterion) {
      aggfun <- function(job, res) lapply(res, gcvec)
      lb.gcdpremerge <- do.call(batchtools::reduceResultsList,
                                c(list(reg=reg, ids=pulsar.jobs, fun=aggfun), reduceGCDargs))
    }
    if (cleanup) unlink(reg$file.dir, recursive=TRUE)
    if (lb.est$opt.index == 1)
      warning("Accurate lower bound could not be determined with the first 2 subsamples")
    if (ub.stars) {
      # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
      pmean <- sapply(lb.est$merge, function(x) sum(x)/(p*(p-1)))
      ub.summary <- cummax(4*pmean*(1-pmean))
      tmpub      <- .starsind(ub.summary, thresh, 1)
      if (any(ub.summary == 0)) ## adjust upper bound to exclude empty graphs
        ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
      else
        ub.index <- max(tmpub, 1)
      
    } else ub.index <- 1
    ## select middle of the lambda path
    fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
    nlams        <- length(fargs$lambda)
    reduceargs   <- list(init=lb.starsmerge[ub.index:lb.est$opt.index])
    
    ## process initial estimate if gcd is selected
    if ('gcd' %in% criterion || 'gcd_new' %in% criterion) {
      reduceGCDargs <- list(init=lapply(lb.gcdpremerge,
                                        function(gcdpm) gcdpm[ub.index:lb.est$opt.index]))
    }
    regdir <- gsub(paste("_", init, sep=""), "", regdir)
    ## Run graph estimation on the rest of the subsamples
    isamp <- ind.sample[-(1:minN)]
    out   <- batchply(data, estFun, fun, fargs, isamp,
                      wkdir, regdir, conffile, job.res)
    reg <- out$reg ; id <- out$id
    doneRun <- batchtools::waitForJobs(reg=reg, id)
    jdone   <- batchtools::findDone(reg=reg, id)
    pulsar.jobs <- intersect((1+refit):rep.num, jdone$job.id)
  }
  rep.num <- length(pulsar.jobs) # adjust denominator to successfull jobs
  if (lb.stars) rep.num <- rep.num + minN
  if (!doneRun)
    warning(paste("Only", jdone, "jobs completed... proceeding anyway"))
  
  for (i in 1:length(criterion)) {
    crit <- criterion[i]
    if (crit == "stars") {
      ## Reduce results, include init estimate from N=2 if lb/ub is used ##
      starsmerge <- do.call(batchtools::reduceResults,
                            c(list(reg=reg, ids=pulsar.jobs, fun=starsaggfun), reduceargs))
      est$stars  <- stars.stability(NULL, thresh, rep.num, p, starsmerge)
    }
    
    if (crit == "gcd") {
      gcdaggfun   <- function(res) lapply(res, gcvec)
      gcdpremerge <- c(reduceGCDargs$init,
                       batchtools::reduceResultsList(reg=reg, ids=pulsar.jobs, fun=gcdaggfun))
      gcdmerge    <- lapply(1:nlams, function(i) dist(t(sapply(1:rep.num, function(j) gcdpremerge[[j]][[i]]))))
      est$gcd <- gcd.stability(NULL, thresh, rep.num, p, nlams, gcdmerge)
    }
    else if (crit == "gcd_new") {
      gcdaggfun   <- function(res) lapply(res, gcvec)
      gcdpremerge <- c(reduceGCDargs$init,
                       batchtools::reduceResultsList(reg=reg, ids=pulsar.jobs, fun=gcdaggfun))
      gcdmerge    <- lapply(1:nlams, function(i) dist(t(sapply(1:rep.num, function(j) gcdpremerge[[j]][[i]]))))
      est$gcd_new <- gcd.stability_extended(NULL, thresh, rep.num, p, nlams, gcdmerge, known_graph = known_graph, orbind = orbind, five = five)
      
    }else if (crit == "estrada") {
      if (!("stars" %in% criterion))
        warning('Need StaRS for computing Estrada classes... not run')
      else
        est$estrada <- estrada.stability(est$stars$merge, thresh, rep.num, p, nlams)
    }
    
    else if (crit == "sufficiency") {
      if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
      else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
    }
    
  }
  
  if (lb.stars) {
    ## split indices of init and full stars estimate
    pind <- ub.index:lb.est$opt.index
    pinv <- setdiff(1:length(lb.est$summary), pind)
    ## stitch back together init and full stars summaries
    tmpsumm       <- vector('numeric', length(lb.est$summary))
    tmpsumm[pinv] <- lb.est$summary[pinv]
    tmpsumm[pind] <- est$stars$summary
    est$stars$summary <- tmpsumm
    ## stitch back together init and full stars merges
    tmpmerg <- vector('list', length(lb.est$summary))
    tmpmerg[pinv]   <- lb.est$merge[pinv]
    tmpmerg[pind]   <- est$stars$merge
    est$stars$merge <- tmpmerg
    ## record stars-related indices
    est$stars$lb.index  <- lb.est$opt.index
    est$stars$ub.index  <- ub.index
    est$stars$opt.index <- est$stars$opt.index + ub.index - 1
  }
  
  if (cleanup) unlink(reg$file.dir, recursive=TRUE)
  est$id  <- id
  est$reg <- reg
  
  if ("stars" %in% criterion) {
    if (est$stars$opt.index == 1) {
      direction <- if (any(est$stars$summary >= .1)) "larger" else "smaller"
      warning(paste("Optimal lambda may be", direction,
                    "than the supplied values"))
    }
  }
  est$call  <- match.call()
  est$est   <- fullmodel
  est$envir <- parent.frame()
  return(structure(est, class = c("batch.pulsar", "pulsar")))
}

#' @keywords internal
.get_batchtools_conffile <- function(conffile) {
  ## try to eval batchtools' default for makeRegistry
  ## otherwise use pulsar's findConfFile
  defconffile <- batchtools::findConfFile()
  if (length(defconffile)==0 || is.null(defconffile))
    defconffile <- findConfFile(conffile)
  defconffile
}

#' @keywords internal
batchply <- function(data, estFun, fun, fargs, ind.sample, wkdir, regdir,
                     conffile, job.res) {
  reg <- batchtools::makeRegistry(file.dir=regdir, work.dir=wkdir,
                                  conf.file=findConfFile(conffile))
  args <- list(fargs=fargs, data=data, fun=fun)
  id   <- batchtools::batchMap(estFun, ind.sample, more.args=args, reg=reg)
  doneSub <- batchtools::submitJobs(reg=reg, resources=job.res)
  return(list(reg=reg, id=id))
}

# batchtools_config_template.R
#' find config file
#'
#' Find a default config file. First calls \code{batchtools::findConfFile} and then find a pulsar default.
#'
#' @param name name of default config or path to config file.
#' @examples
#' ## Default config file provided by pulsar runs code in interactive mode
#' ## This is for testing purposes and executes serially.
#' findConfFile()
#' ## Use the parallel package
#' ## slower than providing the 'ncores' argument to pulsar function, due to
#' ## the overhead of creating the batchtools registry.
#' findConfFile('parallel')
#'
#' ## Use the snow package to register/execute batch jobs on socket clusters.
#' findConfFile('snow')
#' ## Use a TORQUE / PBS queing system. Requires brew template file.
#' findConfFile('torque')
#' findTemplateFile('simpletorque')
#'
#' @details
#' See the batchtools functions \code{batchtools::findConfFile} and \code{batchtools::makeRegistry}. When calling \code{batch.pulsar}, we attempt to use batchtool's default lookup for a config file before calling \code{pulsar::findConfFile}.
#'
#' For clusters with a queuing submission system, a template file, for
#' defining worker node resources and executing the batch R code, will need to
#' be defined somewhere on the system. See \code{\link{findTemplateFile}}.
#' @seealso \code{\link{findTemplateFile}}
#' @export
findConfFile <- function(name='') {
  ## if x is not a file
  ## look for config file using batchtools rules,
  ## otherwise, look in the pulsar system package
  
  conffile <- batchtools::findConfFile()
  if (!is.na(conffile)) return(conffile)
  
  if (checkmate::testFileExists(name, access = "r"))
    return(fs::path_real(name))
  
  ## append type to file extension for default config files
  if (nchar(name)==0) name <- '.R'
  else name <- paste0('.', tools::file_path_sans_ext(name), '.R')
  
  conffile <- fs::path_real(system.file('config',
                                        sprintf('batchtools.conf%s', name), package='pulsar'))
  # }
  if (checkmate::testFileExists(conffile, access = "r")) return(conffile)
  else return(character(0))
}

#' find template file
#'
#' Find a config file from batchtools or default file from pulsar
#'
#' @param name name of default template or path to template file.
#' @examples
#'  \dontrun{
#'  cluster.functions = batchtools::makeClusterFunctionsTORQUE(
#'                      template=pulsar::findTemplateFile('simpletorque'))
#'  }
#' @seealso findConfFile
#' @details
#' See the batchtools functions \code{batchtools::findTemplateFile}, \code{batchtools::makeClusterFunctionsTORQUE}, \code{batchtools::makeClusterFunctionsSGE}, etc, to employ batchtools' default lookup scheme for template files. Supply the output of this function to the \code{template} argument to override batchtools' default.
#'
#' In this case we look for "[name].tmpl" in the pulsar installation directory in the subfolder "templates".
#' @export
findTemplateFile <- function(name) {
  ## get non exported function
  #  x   <- tryCatch(.batchtools_findTemplateFile(name), error=function(e) '')
  #  if (checkmate::testFileExists(x, access = "r")) return(fs::path_real(x))
  #  else {
  x <- system.file("templates", sprintf("%s.tmpl", name), package = "pulsar")
  if (checkmate::testFileExists(x, access = "r")) return(fs::path_real(x))
  else {
    stop(sprintf('Argument \'template\' (=\\"%s\\") must point to a template file or contain the template itself as string (containing at least one newline', name))
  }
  #  }
}

# Generics.R
#' Print a \code{pulsar.refit} S3 object
#'
#' Print information about the model, path length, graph dimension, criterion and optimal indices and graph sparsity.
#'
#' @param x a \code{pulsar.refit}. output from \code{refit}
#' @param ... ignored
#' @importFrom utils capture.output
#' @export
print.pulsar.refit <- function(x, ...) {
  cat("Pulsar-selected refit of", capture.output(print(x$fun)), "\n")
  cat("Path length:", length(x$est$path), "\n")
  cat("Graph dim:  ", ncol(x$est$path[[1]]), "\n")
  crits <- names(x$refit)
  if (length(crits) > 0) {
    critext  <- ifelse(length(crits) > 1, "Criteria:", "Criterion:")
    critext2 <- lapply(crits, function(cr) {
      sp <- sum(x$refit[[cr]]) / ncol(x$refit[[cr]])^2
      optext <- paste(cr, "... sparsity ", signif(sp, 3), sep="")
      paste("  ", optext, sep="")
    })
    cat(critext, "\n", paste(critext2, collapse="\n"), "\n", sep="")
  }
}


#' Print a \code{pulsar} and \code{batch.pulsar} S3 object
#'
#' Print information about the model, path length, graph dimension, criterion and optimal indices, if defined.
#'
#' @param x a fitted \code{pulsar} or \code{batch.pulsar} object
#' @param ... ignored
#' @export
print.pulsar <- function(x, ...) {
  fin <- getArgs(getCall(x), getEnvir(x))
  mode <- ifelse(fin$ncores > 1, "parallel", "serial")
  cat("Mode:", mode)
  .print.pulsar(x, fin)
}

#' @rdname print.pulsar
#' @export
print.batch.pulsar <- function(x, ...) {
  fin <- getArgs(getCall(x), getEnvir(x))
  cat("Mode: batch")
  .print.pulsar(x, fin)
}

#' @keywords internal
.print.pulsar <- function(x, fin) {
  if (fin$lb.stars) {
    cat("... bound index: lower ", x$stars$lb.index,
        ", upper ", x$stars$ub.index, "\n", sep="")
  } else cat("\n")
  cat("Path length:", length(fin$fargs$lambda), "\n")
  cat("Subsamples: ", fin$rep.num, "\n")
  cat("Graph dim:  ", ncol(fin$data), "\n")
  critext  <- ifelse(length(fin$criterion) > 1, "Criteria:", "Criterion:")
  critext2 <- lapply(fin$criterion, function(cr) {
    opt.ind <- x[[cr]]$opt.ind
    optext  <- ifelse(is.null(opt.ind), "",
                      paste("... opt: index ", opt.ind, ", lambda ",
                            signif(fin$fargs$lambda[opt.ind], 3), sep=""))
    paste("  ", cr, optext, sep="")
  })
  cat(critext, "\n", paste(critext2, collapse="\n"), "\n", sep="")
}

#' Plot a \code{pulsar} S3 object
#'
#' @param x a \code{pulsar} or \code{batch.pulsar} object
#' @param scale Flag to scale non-StARS criterion to max StARS value (or 1)
#' @param invlam Flag to plot 1/lambda
#' @param loglam Flag to plot log[lambda]
#' @param legends Flag to plot legends
#' @param ... ignored
#'
#' @details If both invlam and loglam are given, log[1/lambda] is plotted
#' @export
plot.pulsar <- function(x, scale=TRUE, invlam=FALSE, loglam=FALSE, legends=TRUE, ...) {
  .plot.pulsar(x, scale, invlam, loglam, legends)
}

#' @importFrom graphics plot points legend
#' @keywords internal
.plot.pulsar <- function(x, scale=TRUE, invlam=FALSE, loglam=FALSE, legends=TRUE) {
  fin  <- getArgs(getCall(x), getEnvir(x))
  lams <- fin$fargs$lambda
  xlab <- "lambda"
  if (invlam) {lams <- 1/lams ; xlab <- paste("1/", xlab, sep="")}
  if (loglam) {lams <- log(lams) ; xlab <- paste("log[ ", xlab, " ]", sep="")}
  
  nlam  <- length(lams)
  crits <- fin$criterion
  n     <- length(crits)
  if (scale) {
    ylab <- "summary (scaled)"
    if ("stars" %in% crits)
      ymax <- max(x$stars$summary)
    else ymax <- 1
  } else {
    ylab <- "summary"
    ymax <- max(unlist(lapply(crits, function(c) x[[ c ]]$summary)))
  }
  
  yrange <- c(0, ymax)
  plot(lams, seq(yrange[1], yrange[2], length.out=nlam),
       xlab=xlab, ylab=ylab, type='n')
  if (!is.null(x$stars$lb.index)) {
    ilams <- 1:length(lams)
    range1 <- ilams < x$stars$ub.index
    range2 <- ilams > x$stars$lb.index
    range  <- !(range1 | range2)
    ccol   <- vector('numeric', n+1)
    ltys   <- vector('numeric', n+1)
    legs   <- vector('numeric', n+1)
  } else {
    range1 <- rep(FALSE, nlam) ; range2 <- range1
    range  <- !range1
    ccol   <- vector('numeric', n)
    ltys   <- vector('numeric', n)
    legs   <- vector('numeric', n)
  }
  
  i <- 1 ; lcol <- 1
  optcrits <- c() ; optcols <- c()
  for (cr in crits) {
    summs <- x[[ cr ]]$summary
    optind <- opt.index(x, cr)
    if (scale && cr != "stars") summs <- ymax*summs/max(summs)
    if (length(summs) == nlam) {
      points(lams[range],  summs[range],  type='b', col=lcol)
      points(lams[range1], summs[range1], type='b', col=lcol, lty=2)
      points(lams[range2], summs[range2], type='b', col=lcol, lty=2)
      optind2 <- optind
      
      if (any(range1 | range2)) {
        ccol[i:(i+1)] <- c(lcol,lcol)
        ltys[i:(i+1)] <- c(2,1)
        legs[i:(i+1)] <- c(paste("b-", cr, sep=""), cr)
        i <- i+1
      } else {
        ccol[i] <- lcol
        ltys[i] <- 1
        legs[i] <- cr
      }
    } else {
      points(lams[range], summs, type='b', col=lcol)
      optind2 <- optind-which(range)[1]+1
      ccol[i] <- lcol
      ltys[i] <- 1
      legs[i] <- cr
    }
    
    if (!is.null(optind)) {
      points(lams[optind], summs[optind2], type='p', cex=1.5, pch=16, col=lcol)
      optcrits <- c(optcrits, cr)
      optcols  <- c(optcols , lcol)
    }
    lcol <- lcol + 1 ; i <- i + 1
  }
  
  if (legends) {
    legend('bottomleft', legs, col=ccol, pch=1, lty=ltys, cex=1.4)
    if (length(optcrits) > 0)
      legend('topright', optcrits, pch=16, col=optcols, cex=1.5, title="opt lambda")
  }
}

#' Update a pulsar call
#'
#' Update a pulsar model with new or altered arguments. It does this by extracting the call stored in the object, updating the call and (by default) evaluating it in the environment of the original \code{pulsar} call.
#'
#' @param object a n existing pulsar or batch.pulsar object
#' @param ... arguments to \code{pulsar} to update
#' @param evaluate Flag to evaluate the function. If \code{FALSE}, the updated call is returned without evaluation
#' @details The \code{update} call is evaluated in the environment specified by the \code{pulsar} or \code{batch.pulsar} object, so if any variables were used for arguments to the original call, unless they are purposefully updated, should not be altered. For example, if the variable for the original data is reassigned, the output of \code{update} will not be on the original dataset.
#' @return If \code{evaluate = TRUE}, the fitted object - the same output as \code{pulsar} or \code{batch.pulsar}. Otherwise, the updated call.
#' @examples
#' \dontrun{p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(getMaxCov(dat$data), .01, len=20)
#'
#' ## Run pulsar with huge
#' hugeargs <- list(lambda=lams, verbose=FALSE)
#' out.p <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion='stars')
#'
#' ## update call, adding bounds
#' out.b <- update(out.p, lb.stars=TRUE, ub.stars=TRUE)
#' }
#' @importFrom stats update
#' @seealso \code{eval}, \code{\link{update}}, \code{\link{pulsar}}, \code{\link{batch.pulsar}}
#' @export
update.pulsar <- function(object, ..., evaluate=TRUE) {
  extras <- match.call(expand.dots=FALSE)$...
  .update.pulsar(object, extras, evaluate)
}

#' @importFrom stats getCall
#' @keywords internal
.update.pulsar <- function(object, extras, evaluate) {
  call <- getCall(object)
  if (is.null(getEnvir(object))) object$envir <- parent.frame()
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, getEnvir(object))
  else call
}

#' Get calling environment
#'
#' Generic S3 method for extracting an environment from an S3 object. A getter for an explicitly stored environment from an S3 object or list... probably the environment where the original function that created the object was called from. The default method is a wrapper for \code{x$envir}.
#'
#' @param x S3 object to extract the environment
#' @seealso \code{getCall}, \code{environment}, \code{parent.env}, \code{eval}
#' @export
getEnvir <- function(x) {
  UseMethod("getEnvir")
}

#' @rdname getEnvir
#' @export
getEnvir.default <- function(x) {
  getElement(x, "envir")
}

# graphFeatures.R
### Supporting functions to compute features of a graph ##
##   represented as an adjacency matrix (can be sparse) ##
##


#' Graph dissimilarity
#'
#' Dissimilarity matrix of a graph is here defined as the number of neighbors shared by any two nodes.
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph.
#' @param sim Flag to return Graph similarity instead (1-dissimilarity)
#' @param loops Flag to consider self loops
#'
#' @return a \eqn{p*p} dissimilarity matrix
#' @references Bochkina, N. (2015). Selection of the Regularization Parameter in Graphical Models using a Priori Knowledge of Network Structure, arXiv: 1509.05326.
#' @export
graph.diss <- function(G, sim=FALSE, loops=FALSE) {
  dmat <- GraphDiss2(G)
  dmat[is.na(dmat)] <- 1
  if (!loops) diag(dmat) <- 0
  if (sim) dmat <- 1-dmat
  dmat
}

#' @keywords internal
GraphDiss2 <- function(G) {
  Gprod <- G %*% G
  Gdiag   <- Matrix::diag(Gprod)
  degProd <- Gdiag %*% t(Gdiag)
  1 - (Gprod / sqrt(degProd))
}


#' Natural Connectivity
#'
#' Compute the natural connectivity of a graph
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph. Ignored if \code{eig} is given
#' @param eig precomputed list of eigen vals/vectors (output from \code{eigen}). If NULL, compute for \code{G}.
#' @param norm should the natural connectivity score be normalized
#'
#' @details The natural connectivity of a graph is a useful robustness measure of complex networks, corresponding to the average eigenvalue of the adjacency matrix.
#' @return numeric natural connectivity score
#' @references Jun, W., Barahona, M., Yue-Jin, T., & Hong-Zhong, D. (2010). Natural Connectivity of Complex Networks. Chinese Physics Letters, 27(7), 78902. doi:10.1088/0256-307X/27/7/078902
#' @export
natural.connectivity <- function(G, eig=NULL, norm=TRUE) {
  if (is.null(eig)) {
    eig <- eigen(G)
  }
  estrind <- exp(eig$values)
  nc <- log(mean(estrind))
  if (norm) {
    n <- length(estrind)
    nc <- nc / (n - log(n))
  }
  return(nc)
}


#' @keywords internal
.adj2elist <- function(G) {
  if (inherits(G, "sparseMatrix")) {
    G <- Matrix::triu(G, k=1)
    return(Matrix::summary(G)[,-3])
  } else {
    p <- ncol(G)
    return(arrayInd(which(as.logical(triu(G))), c(p,p)))
  }
}

#' Graphlet correlation vector
#'
#' Compute graphlet correlations over the desired orbits (default is 11 non-redundant orbits of graphlets of size <=4) for a single graph \code{G}
#'
#' @param G a \eqn{p*p} adjacency matrix (dense or sparse) of a graph.
#' @param orbind index vector for which orbits to use for computing pairwise graphlet correlations. Default is from Yaveroğlu et al, 2014 (see References), but 1 offset needed for R-style indexing.
#'
#' @references Hočevar, T., & Demšar, J. (2014). A combinatorial approach to graphlet counting. Bioinformatics (Oxford, England), 30(4), 559–65. doi:10.1093/bioinformatics/btt717
#' @references Yaveroğlu, Ö. N., Malod-Dognin, N., Davis, D., Levnajic, Z., Janjic, V., Karapandza, R., … Pržulj, N. (2014). Revealing the hidden language of complex networks. Scientific Reports, 4, 4547. doi:10.1038/srep04547
#' @importFrom stats cor
#' @importFrom methods as
#' @export
gcvec <- function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1) {
  if (length(orbind) < 2) stop("Only one orbit selected, need at least two to calculate graphlet correlations")
  if (any(orbind > 15))   stop("Only 15 orbits, from 4-node graphlets, can be selected")
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  gcount <- orca::count4(Elist)
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],1), method='spearman'))
  gcor[upper.tri(gcor)]
}

#' @keywords internal
subgraph.centrality <- function(Graph, eigs=NULL, rmdiag=FALSE) {
  ## Code from estrada, for undirected graph represented by adjacency matrix M
  ## also return odd/even contributions and the first eigen vector
  if (rmdiag) diag(Graph) <- 0
  # Calculate the subgraph centrality
  #  if (inherits(Graph, 'Matrix'))
  #    eigs <- eigs_sym(Graph, ncol(Graph)-1)
  #  else
  if (is.null(eigs))
    eigs <- eigen(Graph)
  l <- eigs$value
  v <- eigs$vector
  v2 <- v^2
  dl  <- l
  edl <- exp(dl)
  fb  <- v2 %*% edl #A vector of sugraph centralities
  
  # Partition into odd and even contributions of the subgraph centrality
  sinhl <- sinh(dl)
  fbodd <- v2 %*% sinhl
  
  coshl <- cosh(dl)
  feven <- v2 %*% coshl
  out <- list(central=fb, odd=fbodd, even=feven, evec=v, evals=l)
  class(out) <- 'subgraph.centrality'
  return(out)
}

#' @keywords internal
.SMA <- function(x) ((mean(abs(x))))


#' Estrada class
#'
#' Estrada proposes that graphs can be classified into four different classes. We call this the Estrada class.
#' These are:
#'   I. Expander-like
#'  II. Cluster
#' III. Core-Periphery
#' IV.  Mixed.
#' @param G a \eqn{p*p} adjacency matrix of a Graph
#' @param evthresh tolerance for a zero eigenvalue
#' @return Estrada class (\eqn{1-4})
#' @references Estrada, E. (2007). Topological structural classes of complex networks. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(1), 1-12. doi:10.1103/PhysRevE.75.016103
#' @export
estrada.class <- function(G, evthresh=1e-3) {
  if (!inherits(G, "subgraph.centrality"))
    G <- subgraph.centrality(G)
  
  ev1  <- G$evec[,1]
  eval <- G$evals[1]
  
  if (length(unique(sign(ev1))) == 1 && G$evals[2] > 0) { ## try the second eigen vector
    ev1  <- G$evec[,2]
    eval <- G$evals[2]
  }
  subgodd <- G$odd
  Evratio <- pmax(ev1^2 * sinh(eval) / subgodd, evthresh)
  Evratio[is.nan(Evratio)] <- 0
  if (sum(Evratio==evthresh) > (2/3)*length(Evratio)) return(0)
  delLogEv1 <- log10(sqrt(Evratio))
  delSplit <- split(delLogEv1, sign(delLogEv1))
  Devs  <- lapply(delSplit, .SMA)
  if (is.null(Devs$`-1`) && is.null(Devs$`1`)) return(0)
  if (is.null(Devs$`-1`))  Devs$`-1` <- 0
  else if (is.null(Devs$`1`))   Devs$`1` <-  0
  
  if (length(Devs) != 2) return(0) ## will we ever get here?
  if (log10(Devs$`1`+ 1e-3)  > -2.1) eclass <- c(3,4)    else eclass <- c(1,2)
  if (log10(Devs$`-1`+1e-3) > -2.1)  eclass <- eclass[2] else eclass <- eclass[1]
  return(eclass)
}

# mcPulsarSelect.R
#' @keywords internal
.lamcheck <- function(lams) {
  if (is.null(lams)) {
    stop(paste('Error: missing members in fargs:',
               paste(c('lambda')[c(is.null(lams))])))
  } else {
    if (!all(lams == cummin(lams)))
      warning("Are you sure you don't want the lambda path to be monotonically decreasing")
    if (length(lams) < 2)
      warning("Only 1 value of lambda is given. Are you sure you want to do model selection?")
  }
}

#' @keywords internal
.ratcheck <- function(subsample.ratio, n) {
  if (is.null(subsample.ratio)) {
    if (n > 144)
      return(10 * sqrt(n)/n)
    else
      return(0.8)
  } else return(subsample.ratio)
}

#' @keywords internal
.critcheck0 <- function(criterion, knowncrits) {
  if (!all(criterion %in% knowncrits)) {
    stop(paste('Unknown criterion', paste(criterion[!(criterion %in% knowncrits)],
                                          collapse=", "), sep=": "))
  }
  starsdepend <- c("estrada", "sufficiency")
  if (any(starsdepend %in% knowncrits)) {
    if (any(starsdepend %in% criterion) && !("stars" %in% criterion)) {
      stop(paste('Criterion: ', paste(starsdepend[starsdepend %in% criterion],
                                      collapse=", "), ' cannot be run unless stars is also a selected criterion', sep=""))
    }
  }
  
}

#' @importFrom parallel mclapply
#' @keywords internal
.tlist <- function(li, n, m) {
  if (missing(m)) m <- length(li)
  lapply(1L:n, function(i) {
    lapply(1L:m, function(j) li[[j]][[i]])
  })
}

#' @keywords internal
.try_mclapply <- function(X, FUN, mc.cores = getOption("mc.cores", 2L),
                          pass.errors=TRUE, ...) {
  ## capture errors/warnings from mclapply
  warn <- NULL
  env <- environment()
  withCallingHandlers({out <- mclapply(X, FUN, mc.cores=mc.cores, ...)},
                      warning=function(w) {
                        # why doesn't assignment mode work when ncores>1
                        assign('warn', c(warn, w$message), env)
                        invokeRestart("muffleWarning")
                      })
  ## handle errors
  errors <- sapply(out, inherits, what="try-error")
  if (any(errors)) {
    error <- table(trimws(sapply(out[errors], '[', 1)))
    msg   <- paste(sprintf('\n%s job%s failed with: "%s"',
                           error, ifelse(error>1, 's', ''), names(error)), collapse="" )
    if (!pass.errors || all(errors)) {
      stop(msg, call.=FALSE)
    } else {
      warning(msg, call.=FALSE)
      ## continue with successful jobs
      out <- out[!errors]
    }
  } else {
    ## will only invoke if ncore=1, otherwise mclapply suppresses warnings
    if (length(warn) > 0) {
      twarn <- table(trimws(sapply(warn, '[', 1)))
      msg   <- paste(sprintf('%s job%s had warning: "%s"',
                             twarn, ifelse(twarn>1, 's', ''), names(twarn)),
                     collapse="\n")
      warning(msg, call.=FALSE)
    }
  }## no errors detected, continue
  
  # reset previous warn option
  # options(warn=warn.opt)
  attr(out, 'errors') <- errors
  return(out)
}


#' pulsar: serial or parallel mode
#'
#' Run pulsar using StARS' edge stability (or other criteria) to select an undirected graphical model over a lambda path.
#'
#' @param data A \eqn{n*p} matrix of data matrix input to solve for the \eqn{p*p} graphical model
#' @param fun pass in a function that returns a list representing \eqn{p*p} sparse, undirected graphical models along the desired regularization path. The expected inputs to this function are: a data matrix input and a sequence of decreasing lambdas and must return a list or S3 object with a member \emph{named} \code{path}. This should be a list of adjacency matrices for each value of \code{lambda}. See \code{\link{pulsar-function}} for more information.
#' @param fargs arguments to argument \code{fun}. Must be a named list and requires at least one member \code{lambda}, a numeric vector with values for the penalty parameter.
#' @param criterion A character vector of selection statistics. Multiple criteria can be supplied. Only StARS can be used to automatically select an optimal index for the lambda path. See details for additional statistics.
#' @param thresh threshold (referred to as scalar \eqn{\beta} in StARS publication) for selection criterion. Only implemented for StARS. \code{thresh=0.1} is recommended.
#' @param subsample.ratio determine the size of the subsamples (referred to as \eqn{b(n)/n}). Default is 10*sqrt(n)/n for n > 144 or 0.8 otherwise. Should be strictly less than 1.
#' @param rep.num number of random subsamples \eqn{N} to take for graph re-estimation. Default is \eqn{N=20}, but more is recommended for non-StARS criteria or if using edge frequencies as confidence scores.
#' @param seed A numeric seed to force predictable subsampling. Default is NULL. Use for testing purposes only.
#' @param lb.stars Should the lower bound be computed after the first \eqn{N=2} subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path.
#' @param ub.stars Should the upper bound be computed after the first \eqn{N=2} subsamples (should result in considerable speedup and only implemented if stars is selected). If this option is selected, other summary metrics will only be applied to the smaller lambda path. This option is ignored if the lb.stars flag is FALSE.
#' @param ncores number of cores to use for subsampling. See \code{batch.pulsar} for more parallelization options.
#' @param refit Boolean flag to refit on the full dataset after pulsar is run. (see also \code{\link{refit}})
#'
#' @return an S3 object of class \code{pulsar} with a named member for each stability metric run. Within each of these are:
#' \itemize{
#'    \item summary: the summary statistic over \code{rep.num} graphs at each value of lambda
#'    \item criterion: the stability criterion used
#'    \item merge: the raw statistic over the \code{rep.num} graphs, prior to summarization
#'    \item opt.ind: index (along the path) of optimal lambda selected by the criterion at the desired threshold. Will return \eqn{0} if no optimum is found or \code{NULL} if selection for the criterion is not implemented.
#'   }
#' If \code{stars} is included as a criterion then additional arguments include
#' \itemize{
#'    \item lb.index: the lambda index of the lower bound at \eqn{N=2} samples if \code{lb.stars} flag is set to TRUE
#'    \item ub.index: the lambda index of the upper bound at \eqn{N=2} samples if \code{ub.stars} flag is set to TRUE
#'}
#' @return call: the original function call
#' @details
#' The options for \code{criterion} statistics are:
#' \itemize{
#'    \item stars (Stability approach to regularization selection)
#'    \item gcd   (Graphet correlation distance, requires the \pkg{orca} package) see \code{\link{gcvec}}
#'    \item diss  (Node-node dissimilarity) see \code{\link{graph.diss}}
#'    \item estrada (estrada class) see \code{\link{estrada.class}}
#'    \item nc  (natural connectivity) see \code{\link{natural.connectivity}}
#'    \item sufficiency (Tandon & Ravikumar's sufficiency statistic)
#' }
#' @examples
#'\dontrun{
#' ## Generate the data with huge:
#' library(huge)
#' p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(getMaxCov(dat$data), .01, len=20)
#'
#' ## Run pulsar with huge
#' hugeargs <- list(lambda=lams, verbose=FALSE)
#' out.p <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion='stars')
#'
#' ## Run pulsar in bounded stars mode and include gcd metric:
#' out.b <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion=c('stars', 'gcd'),
#'                 lb.stars=TRUE, ub.stars=TRUE)
#' plot(out.b)
#' }
#' @importFrom Matrix mean triu
#' @importFrom parallel mclapply
#' @references Müller, C. L., Bonneau, R., & Kurtz, Z. (2016). Generalized Stability Approach for Regularized Graphical Models. arXiv. https://arxiv.org/abs/1605.07072
#' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. Proceedings of the Twenty-Third Annual Conference on Neural Information Processing Systems (NIPS).
#' @references Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. The Journal of Machine Learning Research, 13, 1059–1062.
#' @seealso \code{\link{batch.pulsar}} \code{\link{refit}}
#' @export
pulsar <- function(data, fun=huge::huge, fargs=list(),
                   criterion=c("stars"),
                   thresh = 0.1, subsample.ratio = NULL,
                   rep.num = 20, seed=NULL,
                   lb.stars=FALSE, ub.stars=FALSE,
                   ncores = 1, refit=TRUE, known_graph = FALSE, five = FALSE, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1)  {
  gcinfo(FALSE)
  n <- nrow(data)
  p <- ncol(data)
  # min requirements for function args
  .lamcheck(fargs$lambda)
  nlams <- length(fargs$lambda)
  knowncrits <- c("stars", "diss", "estrada", "gcd", "nc", "sufficiency", "gcd_new") #"vgraphlet", "egraphlet",
  .critcheck0(criterion, knowncrits)
  subsample.ratio <- .ratcheck(subsample.ratio, n)
  
  if (!is.null(seed)) set.seed(seed)
  ind.sample <- replicate(rep.num,
                          sample(c(1L:n), floor(n * subsample.ratio),
                                 replace = FALSE), simplify=FALSE)
  if (refit) {
    tmp <- 1L:n
    attr(tmp, 'full') <- TRUE
    ind.sample <- c(list(tmp), ind.sample)
  }
  if (!is.null(seed)) set.seed(NULL)
  ## wrap estimator
  estFun <- function(ind.sample, fargs) {
    tmp <- do.call(fun, c(fargs, list(data[ind.sample,])))
    if (!('path' %in% names(tmp)))
      stop('Error: expected data stucture with \'path\' member')
    
    if (isTRUE(attr(ind.sample, 'full')))
      return(tmp)
    else
      return(tmp$path)
  }
  
  if (lb.stars) {
    if (!("stars" %in% criterion))
      stop('Lower/Upper bound method must be used with StARS')
    minN <- 2L + refit # Hard code for now
  } else minN <- rep.num + refit
  
  isamp <- ind.sample[1L:minN]
  ## don't pass on errors if lb.stars = TRUE
  premerge <- .try_mclapply(isamp, estFun, fargs = fargs, mc.cores = ncores,
                            mc.preschedule = FALSE, pass.errors = !lb.stars)
  errors <- attr(premerge, 'errors')
  # Adjust rep.num for failed jobs
  rep.num <- rep.num - ifelse(refit, sum(errors[-1]), sum(errors))
  
  if (refit) {
    fullmodel <- premerge[[1]]
    premerge  <- premerge[-1]
    minN <- minN - 1
  } else fullmodel <- NULL
  
  if (lb.stars) {
    lb.premerge       <- premerge
    lb.premerge.reord <- lapply(1L:nlams, function(i)
      lapply(1L:minN, function(j) lb.premerge[[j]][[i]]))
    
    lb.est <- stars.stability(lb.premerge.reord, thresh, minN, p)
    if (lb.est$opt.index == 1)
      warning("Accurate lower bound could not be determined with the first 2 subsamples")
    if (ub.stars) {
      # upper bound is determined by equivilent of MaxEnt of Poisson Binomial
      pmean      <- sapply(lb.est$merge, function(x) { sum(x)/(p*(p-1)) })
      ub.summary <- cummax(4*pmean*(1-pmean))
      tmpub      <- .starsind(ub.summary, thresh, 1)
      if (any(ub.summary == 0))  ## adjust upper bound to exclude empty graphs
        ub.index <- max(tmpub, max(which(ub.summary == 0))+1)
      else
        ub.index <- max(tmpub, 1)
    } else ub.index <- 1
    # reselect lambda between bounds
    fargs$lambda <- fargs$lambda[ub.index:lb.est$opt.index]
    nlams <- length(fargs$lambda)
    lb.premerge  <- lapply(lb.premerge,
                           function(ppm) ppm[ub.index:lb.est$opt.index])
    isamp <- ind.sample[-(1L:(minN+refit))]
    #    tmp   <- mclapply(isamp, estFun, fargs=fargs, mc.cores=ncores, mc.preschedule = FALSE)
    tmp <- .try_mclapply(isamp, estFun, fargs = fargs, mc.cores = ncores,
                         mc.preschedule = FALSE)
    # Adjust rep.num for failed jobs
    rep.num <- rep.num - sum(attr(tmp, 'errors'))
    premerge <- c(lb.premerge, tmp)
  }
  
  premerge.reord <- .tlist(premerge, nlams, rep.num)
  rm(premerge) ; gc()
  est <- list()
  
  for (i in 1L:length(criterion)) {
    crit <- criterion[i]
    
    if (crit == "stars")
      est$stars <- stars.stability(premerge.reord, thresh, rep.num, p)
    
    else if (crit == "diss")
      est$diss <-  diss.stability(premerge.reord, thresh, rep.num, p, nlams)
    
    else if (crit == "estrada") {
      if (!("stars" %in% criterion))
        warning('Need StaRS for computing Estrada classes... not run')
      else
        est$estrada <- estrada.stability(est$stars$merge,thresh,rep.num,p,nlams)
    }
    
    else if (crit == "sufficiency") {
      if (!("stars" %in% criterion)) warning('Need StaRS for computing sufficiency... not run')
      else  est$sufficiency <- sufficiency(est$stars$merge, rep.num, p, nlams)
    }
    
    #      else if (crit == "egraphlet")
    #        est$egraphlet <- egraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)
    
    #      else if (crit == "vgraphlet")
    #        est$vgraphlet <- vgraphlet.stability(premerge.reord, thresh, rep.num, p, nlams)
    
    else if (crit == "gcd") {
      est$gcd <- gcd.stability(premerge.reord, thresh, rep.num, p, nlams)
      est$gcd$opt.index = est$gcd$opt.index + ub.index
      
    }else if (crit == "gcd_new"){
      est$gcd_new <- gcd.stability_extended(premerge.reord, thresh, rep.num, p, nlams, known_graph = known_graph, orbind = orbind, five = five)
      est$gcd_new$opt.index = est$gcd_new$opt.index + ub.index
      
    }else if (crit == "nc")
      est$nc <- nc.stability(premerge.reord, thresh, rep.num, p, nlams)
    
  }
  
  if (lb.stars) {
    find <- 1:length(lb.est$summary)
    pind <- ub.index:lb.est$opt.index
    pinv <- setdiff(find, pind)
    tmpsumm <- vector('numeric', length(lb.est$summary))
    tmpsumm[pinv] <- lb.est$summary[pinv]
    tmpsumm[pind] <- est$stars$summary
    est$stars$summary   <- tmpsumm
    
    tmpmerg <- vector('list', length(lb.est$summary))
    tmpmerg[pinv]   <- lb.est$merge[pinv]
    tmpmerg[pind]   <- est$stars$merge
    est$stars$merge <- tmpmerg
    
    est$stars$lb.index  <- lb.est$opt.index
    est$stars$ub.index  <- ub.index
    est$stars$opt.index <- est$stars$opt.index + ub.index - 1L
  }
  
  if ("stars" %in% criterion) {
    if (est$stars$opt.index == 1) {
      direction <- if (any(est$stars$summary >= .1)) "larger" else "smaller"
      warning(paste("Optimal lambda may be", direction, "than the supplied values"))
    }
  }
  est$call  <- match.call()
  est$envir <- parent.frame()
  est$est   <- fullmodel
  return(structure(est, class="pulsar"))
}

#' @keywords internal
.starsind <- function(summary, thresh, offset=1) {
  max(which.max(summary >= thresh)[1] - offset, 1)
}

#' @keywords internal
stars.stability <- function(premerge, stars.thresh, rep.num, p, merge=NULL) {
  if (is.null(stars.thresh)) stars.thresh <- 0.1
  est <- list()
  
  
  # do.call(batchtools::reduceResults,
  #                  c(list(reg=reg, fun=starsaggfun), reduceargs))
  
  if (is.null(merge)) {
    est$merge <- lapply(premerge, function(x) Reduce("+", x))
    gc() # flush
  } else est$merge <- merge
  est$summary <- rep(0, length(est$merge))
  
  for (i in 1:length(est$merge)) {
    est$merge[[i]] <- est$merge[[i]]/rep.num
    est$summary[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]])) / (p * (p - 1))
  }
  ## monotonize variability
  est$summary   <- cummax(est$summary)
  est$opt.index <- .starsind(est$summary, stars.thresh)
  est$criterion <- "stars.stability"
  est$thresh    <- stars.thresh
  return(est)
}

#' @keywords internal
sufficiency <- function(merge, rep.num, p, nlams) {
  ## Merge solution from StARS
  est <- list()
  est$merge <- sapply(merge, function(x) apply(x*(1-x), 2, max))
  est$summary <- colMeans(est$merge)
  est$criterion <- 'sufficiency'
  return(est)
}

#' @keywords internal
.sumsq <- function(x2,y) x2 + y^2

#' @keywords internal
diss.stability <- function(premerge, diss.thresh, rep.num, p, nlams) {
  est <- list()
  disslist  <- lapply(premerge, function(pm) lapply(pm, graph.diss))
  est$merge <- lapply(disslist, function(dissmat) Reduce("+", dissmat)/rep.num)
  mergesq   <- lapply(disslist, function(dissmat)
    Reduce(.sumsq, dissmat[-1], init=dissmat[[1]]^2)/rep.num)
  
  gc() # flush
  est$summary <- rep(0, length(est$merge))
  for (i in 1:length(est$merge)) {
    tmp <- mergesq[[i]] - est$merge[[i]]^2
    est$summary[i] <- sum(triu(tmp))  / (p * (p - 1))
  }
  est$mergesq <- mergesq
  est$criterion <- "diss.stability"
  return(est)
}


#estrada.stability <- function(premerge, thresh, rep.num, p, nlams) {
#    est <- list()
#    estrlist  <- lapply(premerge, function(pm) lapply(pm, estrada.class))
#    est$merge <- lapply(estrlist, function(x) table(unlist(x)))

##    gc() # flush
#    est$summary <- rep(0, length(est$merge))
#    for (i in 1:length(est$merge)) {
#        est$summary[i] <- 1-max(est$merge[[i]])/rep.num
#    }
#    ## monotonize variability
##    est$summary <- cummax(est$summary)
#    if (!is.null(thresh))
#      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
#    else
#      est$opt.index <- 0

#    est$criterion <- "estrada.stability"
#    return(est)
#}

#' @keywords internal
estrada.stability <- function(merge, thresh, rep.num, p, nlams) {
  est <- list()
  est$summary <- unlist(lapply(merge, function(x) estrada.class(x >= .05)))
  if (!is.null(thresh))
    est$opt.index <- max(which.max(est$summary >= thresh)[1] - 1, 1)
  else
    est$opt.index <- 0
  
  est$criterion <- "estrada.stability"
  return(est)
}

#' @keywords internal
nc.stability <- function(premerge, thresh, rep.num, p, nlams) {
  est <- list()
  est$merge <- sapply(premerge, function(x) sapply(x, natural.connectivity))
  est$summary <- colMeans(est$merge)
  est$criterion <- "nc.stability"
  return(est)
}

#' @importFrom stats dist
#' @keywords internal
gcd.stability <- function(premerge, thresh, rep.num, p, nlams, merge=NULL) {
  est <- list()
  if (is.null(merge))
    est$merge <- lapply(premerge, function(pm) dist(t(sapply(pm, gcvec))))
  else
    est$merge <- merge
  
  est$summary <- vector('numeric', nlams)
  for (i in 1:nlams) est$summary[i] <- mean(est$merge[[i]])
  est$criterion <- "graphlet.stability"
  est$opt.index <- max(which.min(est$summary)[1] - 1, 1)
  return(est)
}

#egraphlet.stability <- function(premerge, thresh, rep.num, p, nlams, norbs=12) {
#    est <- list()
##    estrlist    <- lapply(premerge.reord, function(pm) lapply(pm, estrada))
##    est$merge <- lapply(estrlist, function(estrvec) Reduce("+", estrvec)/rep.num)
#    collect <- lapply(premerge, function(pm) lapply(pm, egraphletlist))
#    collect.reord <- lapply(collect, function(colam) lapply(1:norbs, function(i)
#                      lapply(1:rep.num, function(j) colam[[j]][[i]])))
#    rm(list=c('collect')) ; gc()
#    est$merge <- vector('list', nlams)
#    est$summary <- Matrix(0, length(est$merge), norbs)
#    for (i in 1:nlams) {
#      est$merge[[i]] <- vector('list', norbs)
#      for (j in 1:norbs) {
#        collam <- collect.reord[[i]][[j]]
#        EX2 <- Reduce(.sumsq, collam[-1], init=collam[[1]]^2)/rep.num
#        EX  <- Reduce('+', collam)/rep.num
#        est$merge[[i]][[j]] <- EX2 - EX^2
#        est$summary[i,j] <- 2 * sum(est$merge[[i]][[j]]) / (p * (p - 1))
#      }
#    }

##    if (!is.null(thresh))
##      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
##    else
##      est$opt.index <- 0

#    est$criterion <- "egraphlet.stability"
#    return(est)
#}


#vgraphlet.stability <- function(premerge, thresh, rep.num, p, nlams, norbs=15) {

#    est <- list()
#    collect <- lapply(premerge, function(pm) lapply(pm, vgraphletlist))
##    rm(list=c('collect')) ; gc()
#    est$merge <- vector('list', nlams)
#    est$summary <- matrix(0, nlams, norbs)
#    for (i in 1:nlams) {
##      est$merge[[i]] <- vector('list', norbs)
##      for (j in 1:norbs) {
##        collam <- collect.reord[[i]][[j]]
#        EX2 <- Reduce(.sumsq, collect[[i]][-1], init=collect[[i]][[1]]^2)/rep.num
#        EX  <- Reduce('+', collect[[i]])/rep.num
#        est$merge[[i]] <- EX2 - EX^2
#        est$summary[i,] <- 2*Matrix::colMeans(est$merge[[i]])
##      }
#    }
##    if (!is.null(thresh))
##      est$opt.index    <- max(which.max(est$summary >= thresh)[1] - 1, 1)
##    else
##      est$opt.index <- 0
#    est$criterion <- "graphlet.stability"
#    return(est)
#}

# refit.R
#' Refit pulsar model
#'
#' Run the supplied graphical model function on the whole dataset and refit with the selected lambda(s)
#'
#' @param obj a fitted \code{pulsar} or \code{batch.pulsar} object
#' @param criterion a character vector of criteria for refitting on full data. An optimal index must be defined for each criterion or a message will displayed. If missing (no argument is supplied), try to refit for all pre-specified criteria.
#' @details The \code{refit} call is evaluated in the environment specified by the \code{pulsar} or \code{batch.pulsar} object, so if any variables were used for arguments to the original call, unless they are purposefully updated, should not be altered. For example, if the variable for the original data is reassigned, the output of \code{refit} will not be on the original dataset.
#' @return a \code{pulsar.refit} S3 object with members:
#' \itemize{
#'   \item est: the raw output from the graphical model function, \code{fun}, applied to the full dataset.
#'   \item refit: a named list of adjacency matrices, for each optimal criterion in \code{obj} or specified in the \code{criterion} argument.
#'   \item fun: the original function used to estimate the graphical model along the lambda path.
#'}
#' @examples
#'
#' ## Generate the data with huge:
#' \dontrun{
#' library(huge)
#' set.seed(10010)
#' p <- 40 ; n <- 1200
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#' lams  <- getLamPath(getMaxCov(dat$data), .01, len=20)
#'
#' ## Run pulsar with huge
#' hugeargs <- list(lambda=lams, verbose=FALSE)
#' out.p <- pulsar(dat$data, fun=huge::huge, fargs=hugeargs,
#'                 rep.num=20, criterion='stars')
#'
#' fit  <- refit(out.p)
#' }
#' @seealso \code{\link{pulsar}} \code{\link{batch.pulsar}}
#' @export
refit <- function(obj, criterion) {
  UseMethod("refit")
}

#' @export
refit.pulsar <- function(obj, criterion) {
  .refit.pulsar(obj, criterion)
}

#' @keywords internal
.refit.pulsar <- function(obj, criterion) {
  est <- vector('list', 2)
  names(est) <- c('est', 'refit')
  fin <- getArgs(getCall(obj), getEnvir(obj))
  ## call est function on original dataset
  if (length(obj$est)) {
    est$est <- obj$est
  } else {
    est$est <- do.call(eval(fin$fun), c(fin$fargs, list(fin$data)))
  }
  
  if (missing(criterion)) criterion <- eval(fin$criterion)
  est$refit <- vector('list', length(criterion))
  names(est$refit) <- criterion
  for (crit in criterion) {
    optind <- obj[[crit]]$opt.index
    if (!is.null(optind)) {
      est$refit[[crit]] <- est$est$path[[optind]]
    } else {
      est$refit[[crit]] <- NULL
      if (crit %in% names(obj)) {
        message(paste('No optimal index selected for', crit, 'criterion', sep=" "))
      } else
        warning(paste('Unknown criterion', crit, sep=" "), call.=FALSE)
    }
  }
  
  ## TODO: if fun is null, get formal arg of obj
  est$fun <- obj$call$fun
  if (is.null(est$fun))
    est$fun <- formals(class(obj))$fun
  
  structure(est, class='pulsar.refit')
}



# utilities.R
## Supporting and generic functions

#' @keywords internal
getArgs <- function(call, envir=parent.frame()) {
  fin    <- lapply(call, eval, envir=envir)
  forms  <- formals(fin[[1]])
  iscall <- sapply(forms, class) == 'call'
  iscall <- iscall & !(names(forms) %in% c('regdir', 'regid'))
  forms[iscall] <- lapply(forms[iscall], eval)
  c(forms[!(names(forms) %in% names(fin))], fin)
}

#' @keywords internal
.pcheck <- function(obj) {
  if (!inherits(obj, 'pulsar'))
    stop("obj must be pulsar output")
}

#' @keywords internal
.critcheck <- function(obj, criterion=NULL) {
  if (!(criterion %in% names(obj)))
    warning('desired criterion was not used in the pulsar run')
}


#' Get or evaluate an optimal index
#'
#' If the optimal index for the lambda path is not already assigned, then use a validated method to
#' select the optimal index of the lambda path for alternate criteria  (i.e. other than StARS).
#'
#' @param obj the pulsar/batch.pulsar object to evaluate
#' @param criterion a character argument for the desired summary criterion
#' @param ... Ignored
#' @details Automated optimal index selection is [currently] only implemented for \code{gcd} (graphlet stability).
#'
#' Criterion:
#' \itemize{
#'  \item gcd: Select the minimum gcd summary score within the lower and upper StARS bounds.
#' }
#' @return index of the lambda path
#' @seealso \code{\link{opt.index}}
#' @export
get.opt.index <- function(obj, criterion="gcd", ...) {
  optind <- opt.index(obj, criterion)
  if (!is.null(optind)) optind
  if (criterion == 'gcd' || criterion == 'gcd_new') {
    if (is.null(obj$stars$lb.index) || !obj$stars$lb.index)
      stop('Lower bound needed for gcd metric (run with lb.stars=TRUE)')
    gcdind <- which.min(getElement(obj, criterion)$summary)
    gcdind <- gcdind + obj$stars$ub.index - 1
    names(gcdind) <- criterion
    return(gcdind)
  } else {
    stop("Currently, gcd is the only supported criterion")
  }
}

#' Optimal index
#'
#' Get or set the optimal index of the lambda path, as determined by a given criterion. \code{value} must be a numeric/int.
#'
#' @param obj a pulsar or batch.pulsar object
#' @param criterion a summary statistic criterion for lambda selection. If value is not named, default to gcd.
#' @seealso \code{\link{get.opt.index}}
#' @export
opt.index <- function(obj, criterion='gcd') {
  .pcheck(obj)
  .critcheck(obj, criterion)
  getElement(obj, criterion)$opt.index
}

#' @param value Integer index for optimal lambda by criterion
#' @rdname opt.index
#' @export
"opt.index<-" <- function(obj, criterion=names(value), value) {
  .pcheck(obj)
  fin <- getArgs(obj$call, obj$envir)
  .critcheck(obj, criterion)
  if (length(criterion) > 1) stop("Select one criterion")
  if (is.null(criterion)) criterion <- 'gcd'
  if (!is.null(value)) {
    if (!is.numeric(value) || value < 1 || value >= length(fin$fargs$lambda))
      stop('Index value must be positive int within range length of lambda path')
  }
  obj[[ criterion ]]$opt.index <- value
  obj
}

#' Lambda path
#'
#' Generate a lambda path sequence in descending order, equally or log-spaced.
#'
#' @param max numeric, maximum lambda value
#' @param min numeric, minimum lambda value
#' @param len numeric/int, length of lambda path
#' @param log logical, should the lambda path be log-spaced
#' @return numeric vector of lambdas
#' @examples
#' ## Generate the data with huge:
#' library(huge)
#' set.seed(10010)
#' p <- 40 ; n <- 100
#' dat   <- huge.generator(n, p, "hub", verbose=FALSE, v=.1, u=.3)
#'
#' ## Theoretical lamda max is the maximum abs value of the empirical covariance matrix
#' maxCov <- getMaxCov(dat$data)
#' lams   <- getLamPath(maxCov, 5e-2*maxCov, len=40)
#'
#' @seealso \code{\link{getMaxCov}}
#' @export
getLamPath <- function(max, min, len, log=FALSE) {
  if (max < min) stop('Did you flip min and max?')
  if (log) { min <- log(min) ; max <- log(max) }
  lams  <- seq(max, min, length.out=len)
  if (log) exp(lams)
  else lams
}

#' Max value of cov
#'
#' Get the maximum [absolute] value of a covariance matrix.
#'
#' @param x A matrix/Matrix of data or covariance
#' @param cov Flag if \code{x} is a covariance matrix, Set False is \code{x} is an nxp data matrix. By default, if \code{x} is symmetric, assume it is a covariance matrix.
#' @param abs Flag to get max absolute value
#' @param diag Flag to include diagonal entries in the max
#' @details This function is useful to determine the theoretical value for lambda_max - for Gaussian data, but may be a useful starting point in the general case as well.
#' @seealso \code{\link{getLamPath}}
#' @export
getMaxCov <- function(x, cov=isSymmetric(x), abs=TRUE, diag=FALSE) {
  if (!cov) x <- cov(x)
  tmp <- Matrix::triu(x, k=if (diag) 0 else 1)
  tmp <- if (abs) abs(tmp) else tmp
  max(tmp)
}


gcvec_extended <- function(G, orbind=c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)+1, five = FALSE){
  
  diag(G) = 0
  
  Elist <- .adj2elist(G)
  n <- length(orbind)
  if (ncol(Elist) < 1 || nrow(Elist) < 1) {
    return(rep(0, n*(n-1)/2))
  }
  
  p <- ncol(G)
  if(five == TRUE){
    gcount <- orca::count5(Elist)
  }else{
    gcount <- orca::count4(Elist)
  }
  ## expand missing nodes
  buffer <- matrix(0, nrow=p-nrow(gcount), ncol=ncol(gcount))
  gcount <- rbind(gcount, buffer)
  # deprecate direct call to count4 for CRAN submission
  #  gcount <- .C("count4", Elist, dim(Elist),
  #      orbits = matrix(0, nrow = max(Elist), ncol = 15), PACKAGE="orca")$orbits
  ##  # warnings here are due to std dev == 0. This almost always occurs for a completely connected
  ## or completely empty graph and can be safely suppressed.
  gcor <- suppressWarnings(cor(rbind(gcount[,orbind],1), method='spearman'))
  gcor[upper.tri(gcor)]
  
}

GCD <- function(gcv1, gcv2){
  res = dist(rbind(gcv1,gcv2))[1]
  return(res)
}

gcd.stability_extended <- function(premerge, thresh, rep.num, p, nlams, known_graph, orbind = FALSE, five = FALSE, merge=NULL) {
  est <- list()
  
  if(is.null(merge)){
    for(j in 1:length(premerge)){
      res_tot = c()
      estimated_graph_subsampling = list()
      for(i in 1:length(premerge[[j]])){
        #print(length(premerge[[j]]))
        r1 = gcvec_extended(premerge[[j]][[i]], orbind = orbind, five = five)
        r2 = gcvec_extended(known_graph, orbind = orbind, five = five)
        res = GCD(r1, r2)
        res_tot = c(res_tot, res)
        estimated_graph_subsampling[[i]] = premerge[[j]][[i]]
      } 
      est$graph_subsambling[[j]] = estimated_graph_subsampling
      est$merge[[j]] <- res_tot
      est$g1[[j]] = r1
      est$g2 = r2
    }
  }else{
    est$merge <- merge
  }
  
  est$summary <- vector('numeric', nlams)
  for (i in 1:nlams) est$summary[i] <- mean(est$merge[[i]])
  est$criterion <- "gcd_new"
  est$opt.index <- max(which.min(est$summary)[1] -1, 1)
  return(est)
}


#gcd.stability <- function(premerge, thresh, rep.num, p, nlams, merge=NULL) {
#    est <- list()
#    if (is.null(merge))
#        est$merge <- lapply(premerge, function(pm) dist(t(sapply(pm, gcvec))))
#    else
#        est$merge <- merge
#
#    est$summary <- vector('numeric', nlams)
#    for (i in 1:nlams) est$summary[i] <- mean(est$merge[[i]])
#    est$criterion <- "graphlet.stability"
#    est$opt.index <- max(which.min(est$summary)[1] - 1, 1)
#    return(est)
#}




Data1_HUB = simulate_HUB_graph(200, 400,rho_hub = 0.2, seed = 123)
lmax <- getMaxCov(Data1_HUB$simulation_data)
lams <- getLamPath(lmax, .001, len=100)
hugeargs <- list(lambda=lams, method = "glasso", verbose = FALSE)

ub.index = 0
lb.index = 100

GCD_Oracle_results_strategy1 = list()
GCD_Oracle_results_strategy2 = list()
GCD_Oracle_results_strategy3 = list()
GCD_Oracle_results_strategy4 = list()

Hamming_dist_res = list()

GCD_GIP_results_strategy1 = list()
GCD_GIP_results_strategy2 = list()
GCD_GIP_results_strategy3 = list()
GCD_GIP_results_strategy4 = list()

Stars_01_res = list()
Stars_005_res = list()

GCD_GIP_res_strategy1 = list()
GCD_GIP_res_strategy2 = list()
GCD_GIP_res_strategy3 = list()
GCD_GIP_res_strategy4 = list()

GCD_Oracle_res_strategy1 = list()
GCD_Oracle_res_strategy2 = list()
GCD_Oracle_res_strategy3 = list()
GCD_Oracle_res_strategy4 = list()

for(i in 2:11){
  
  Data1_HUB = simulate_HUB_graph(200, 400,rho_hub = 0.2, seed = i)
  
  pulsar_gip_400_HUB_gl = pulsar(Data1_HUB$simulation_data, fun = huge::huge, seed = 1, thresh = 0.05, fargs = hugeargs, criterion = c("stars", "gcd_new"), rep.num = 30, known_graph = Data1_HUB$real_graph)
  print(i)
  #GCD oracle calculation
  GCD_Oracle_400_HUB_scheme_nothing = GCD_oracle_random_nothing(Data1_HUB$real_graph, pulsar_gip_400_HUB_gl,gcvec_extended_random1)
  
  GCD_Oracle_400_HUB_scheme_random_noise = GCD_oracle_random_nothing(Data1_HUB$real_graph, pulsar_gip_400_HUB_gl,gcvec_extended_random2)
  
  GCD_Oracle_400_HUB_scheme_random_vector = GCD_oracle_random(Data1_HUB$real_graph, pulsar_gip_400_HUB_gl,gcvec_extended_random3)
  
  GCD_Oracle_400_HUB_scheme_both = GCD_oracle_random(Data1_HUB$real_graph, pulsar_gip_400_HUB_gl,gcvec_extended_random2)
  
  #GCD GIP calculation
  
  GCD_gip_hub_400_gl = Function_GCD_data(pulsar_gip_400_HUB_gl, Data1_HUB$real_graph)
  
  #HD oracle calculation
  HD_Hub_400 = Hamming_distance_path(pulsar_gip_400_HUB_gl, Data1_HUB$real_graph)
  
  #GCD oracle paths
  GCD_Oracle_results_strategy1 = append(GCD_Oracle_results_strategy1, list(GCD_Oracle_400_HUB_scheme_nothing))
  GCD_Oracle_results_strategy2 = append(GCD_Oracle_results_strategy2, list(GCD_Oracle_400_HUB_scheme_random_noise))
  GCD_Oracle_results_strategy3 = append(GCD_Oracle_results_strategy3, list(GCD_Oracle_400_HUB_scheme_random_vector))
  GCD_Oracle_results_strategy4 = append(GCD_Oracle_results_strategy4, list(GCD_Oracle_400_HUB_scheme_both))
  
  #GCD oracle better results
  GCD_Oracle_res_strategy1 = append(GCD_Oracle_res_strategy1, list(HD_Hub_400[which.min(GCD_Oracle_400_HUB_scheme_nothing)]))
  GCD_Oracle_res_strategy2 = append(GCD_Oracle_res_strategy2, list(HD_Hub_400[which.min(GCD_Oracle_400_HUB_scheme_random_noise)]))
  GCD_Oracle_res_strategy3 = append(GCD_Oracle_res_strategy3, list(HD_Hub_400[which.min(GCD_Oracle_400_HUB_scheme_random_vector)]))
  GCD_Oracle_res_strategy4 = append(GCD_Oracle_res_strategy4, list(HD_Hub_400[which.min(GCD_Oracle_400_HUB_scheme_both)]))
  
  #HD res
  
  Hamming_dist_res = append(Hamming_dist_res, list(HD_Hub_400))
  
  #Stars res
  
  v05 = 101-pulsar_gip_400_HUB_gl$stars$opt.index
  v1 = 101-(which.min(pulsar_gip_400_HUB_gl$stars$summary<=0.1)-1)
  Stars_01_res = append(Stars_01_res, list(HD_Hub_400[v1]))
  Stars_005_res = append(Stars_005_res, list(HD_Hub_400[v05]))
  
  #GCD GIP paths
  
  GCD_GIP_results_strategy1 = append(GCD_GIP_results_strategy1, list(GCD_gip_hub_400_gl$mean_gcd1))
  GCD_GIP_results_strategy2 = append(GCD_GIP_results_strategy2, list(GCD_gip_hub_400_gl$mean_gcd2))
  GCD_GIP_results_strategy3 = append(GCD_GIP_results_strategy3, list(GCD_gip_hub_400_gl$mean_gcd3))
  GCD_GIP_results_strategy4 = append(GCD_GIP_results_strategy4, list(GCD_gip_hub_400_gl$mean_gcd4))
  
  #GCD GIP res
  
  GCD_GIP_res_strategy1 = append(GCD_GIP_res_strategy1, list(HD_Hub_400[which.min(GCD_gip_hub_400_gl$mean_gcd1)]))
  GCD_GIP_res_strategy2 = append(GCD_GIP_res_strategy2, list(HD_Hub_400[which.min(GCD_gip_hub_400_gl$mean_gcd2)]))
  GCD_GIP_res_strategy3 = append(GCD_GIP_res_strategy3, list(HD_Hub_400[which.min(GCD_gip_hub_400_gl$mean_gcd3)]))
  GCD_GIP_res_strategy4 = append(GCD_GIP_res_strategy4, list(HD_Hub_400[which.min(GCD_gip_hub_400_gl$mean_gcd4)]))
  
}



#SAVE the HUB results
vector_GCD_GIP_res_strategy1 = c()
vector_GCD_GIP_res_strategy2 = c()
vector_GCD_GIP_res_strategy3 = c()
vector_GCD_GIP_res_strategy4 = c()
for(i in GCD_GIP_res_strategy1){
  vector_GCD_GIP_res_strategy1 = c(vector_GCD_GIP_res_strategy1,i)
}

for(i in GCD_GIP_res_strategy2){
  vector_GCD_GIP_res_strategy2 = c(vector_GCD_GIP_res_strategy2,i)
}

for(i in GCD_GIP_res_strategy3){
  vector_GCD_GIP_res_strategy3 = c(vector_GCD_GIP_res_strategy3,i)
}

for(i in GCD_GIP_res_strategy4){
  vector_GCD_GIP_res_strategy4 = c(vector_GCD_GIP_res_strategy4,i)
}

vector_GCD_Oracle_res_strategy1 = c()
vector_GCD_Oracle_res_strategy2 = c()
vector_GCD_Oracle_res_strategy3 = c()
vector_GCD_Oracle_res_strategy4 = c()

for(i in GCD_Oracle_res_strategy1){
  vector_GCD_Oracle_res_strategy1 = c(vector_GCD_Oracle_res_strategy1,i)
}

for(i in GCD_Oracle_res_strategy2){
  vector_GCD_Oracle_res_strategy2 = c(vector_GCD_Oracle_res_strategy2,i)
}

for(i in GCD_Oracle_res_strategy3){
  vector_GCD_Oracle_res_strategy3 = c(vector_GCD_Oracle_res_strategy3,i)
}

for(i in GCD_Oracle_res_strategy4){
  vector_GCD_Oracle_res_strategy4 = c(vector_GCD_Oracle_res_strategy4,i)
}



vector_Stars_01_res = c()
vector_Stars_005_res = c()

for(i in Stars_01_res){
  vector_Stars_01_res = c(vector_Stars_01_res,i)
}

for(i in Stars_005_res){
  vector_Stars_005_res = c(vector_Stars_005_res,i)
}


Best_HD = c()

for(i in Hamming_dist_res){
  R = min(i)
  Best_HD = c(Best_HD,R)
}


res_tot_hub = cbind(vector_GCD_GIP_res_strategy1, vector_GCD_GIP_res_strategy2, vector_GCD_GIP_res_strategy3, vector_GCD_GIP_res_strategy4,
                       vector_GCD_Oracle_res_strategy1, vector_GCD_Oracle_res_strategy2, vector_GCD_Oracle_res_strategy3, vector_GCD_Oracle_res_strategy4,
                       vector_Stars_005_res, vector_Stars_01_res, Best_HD)

res_tot_hub_df = as.data.frame(res_tot_hub)




res_HD_1_hub = c()
for(i in Hamming_dist_res){
  res_HD_1_hub = cbind(res_HD_1_hub, i)
}

res_HD_1_hub_df = as.data.frame(res_HD_1_hub)
colnames(res_HD_1_hub_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_Oracle_S1 = c()
for(i in GCD_Oracle_results_strategy1){
  H_res_GCD_Oracle_S1 = cbind(H_res_GCD_Oracle_S1, i)
}

H_res_GCD_Oracle_S1_df = as.data.frame(H_res_GCD_Oracle_S1)
colnames(H_res_GCD_Oracle_S1_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_Oracle_S2 = c()
for(i in GCD_Oracle_results_strategy2){
  H_res_GCD_Oracle_S2 = cbind(H_res_GCD_Oracle_S2, i)
}

H_res_GCD_Oracle_S2_df = as.data.frame(H_res_GCD_Oracle_S2)
colnames(H_res_GCD_Oracle_S2_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_Oracle_S3 = c()
for(i in GCD_Oracle_results_strategy3){
  H_res_GCD_Oracle_S3 = cbind(H_res_GCD_Oracle_S3, i)
}

H_res_GCD_Oracle_S3_df = as.data.frame(H_res_GCD_Oracle_S3)
colnames(H_res_GCD_Oracle_S3_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_Oracle_S4 = c()
for(i in GCD_Oracle_results_strategy4){
  H_res_GCD_Oracle_S4 = cbind(H_res_GCD_Oracle_S4, i)
}

H_res_GCD_Oracle_S4_df = as.data.frame(H_res_GCD_Oracle_S4)
colnames(H_res_GCD_Oracle_S4_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")



H_res_GCD_GIP_S1 = c()
for(i in GCD_GIP_results_strategy1){
  H_res_GCD_GIP_S1 = cbind(H_res_GCD_GIP_S1, i)
}

H_res_GCD_GIP_S1_df = as.data.frame(H_res_GCD_GIP_S1)
colnames(H_res_GCD_GIP_S1_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_GIP_S2 = c()
for(i in GCD_GIP_results_strategy2){
  H_res_GCD_GIP_S2 = cbind(H_res_GCD_GIP_S2, i)
}

H_res_GCD_GIP_S2_df = as.data.frame(H_res_GCD_GIP_S2)
colnames(H_res_GCD_GIP_S2_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")



H_res_GCD_GIP_S3 = c()
for(i in GCD_GIP_results_strategy3){
  H_res_GCD_GIP_S3 = cbind(H_res_GCD_GIP_S3, i)
}

H_res_GCD_GIP_S3_df = as.data.frame(H_res_GCD_GIP_S3)
colnames(H_res_GCD_GIP_S3_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


H_res_GCD_GIP_S4 = c()
for(i in GCD_GIP_results_strategy4){
  H_res_GCD_GIP_S4 = cbind(H_res_GCD_GIP_S4, i)
}

H_res_GCD_GIP_S4_df = as.data.frame(H_res_GCD_GIP_S4)
colnames(H_res_GCD_GIP_S4_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")








write.csv(res_tot_hub_df, "tot_results_HUB_simulation.csv")

write.csv(res_HD_1_hub_df, "simulations_HD_HUB.csv")

write.csv(H_res_GCD_Oracle_S1_df, "simulations_GCD_Oracle_S1_HUB.csv")
write.csv(H_res_GCD_Oracle_S2_df, "simulations_GCD_Oracle_S2_HUB.csv")
write.csv(H_res_GCD_Oracle_S3_df, "simulations_GCD_Oracle_S3_HUB.csv")
write.csv(H_res_GCD_Oracle_S4_df, "simulations_GCD_Oracle_S4_HUB.csv")

write.csv(H_res_GCD_GIP_S1_df, "simulations_GCD_GIP_S1_HUB.csv")
write.csv(H_res_GCD_GIP_S2_df, "simulations_GCD_GIP_S2_HUB.csv")
write.csv(H_res_GCD_GIP_S3_df, "simulations_GCD_GIP_S3_HUB.csv")
write.csv(H_res_GCD_GIP_S4_df, "simulations_GCD_GIP_S4_HUB.csv")




















Data1_Random = simulate_random_graph(200, 400, c(-1, -0.07, 0.07, 1), K = 10, seed = 123, scale_cov = FALSE)
lmax <- getMaxCov(Data1_Random$simulation_data)*1.2
lams <- getLamPath(lmax, .01, len=100)
hugeargs <- list(lambda=lams, method = "glasso", verbose = FALSE)

ub.index = 0
lb.index = 100

R_GCD_Oracle_results_strategy1 = list()
R_GCD_Oracle_results_strategy2 = list()
R_GCD_Oracle_results_strategy3 = list()
R_GCD_Oracle_results_strategy4 = list()

R_Hamming_dist_res = list()

R_GCD_GIP_results_strategy1 = list()
R_GCD_GIP_results_strategy2 = list()
R_GCD_GIP_results_strategy3 = list()
R_GCD_GIP_results_strategy4 = list()

R_Stars_01_res = list()
R_Stars_005_res = list()

R_GCD_GIP_res_strategy1 = list()
R_GCD_GIP_res_strategy2 = list()
R_GCD_GIP_res_strategy3 = list()
R_GCD_GIP_res_strategy4 = list()

R_GCD_Oracle_res_strategy1 = list()
R_GCD_Oracle_res_strategy2 = list()
R_GCD_Oracle_res_strategy3 = list()
R_GCD_Oracle_res_strategy4 = list()

for(i in 2:11){
  
  Data1_Random = simulate_random_graph(200, 400, c(-1, -0.07, 0.07, 1), K = 10, seed = i, scale_cov = FALSE)
  
  pulsar_gip_400_Random_gl = pulsar(Data1_Random$simulation_data, fun = huge::huge, seed = 1, thresh = 0.05, fargs = hugeargs, criterion = c("stars", "gcd_new"), rep.num = 30, known_graph = Data1_Random$real_graph)
  print(i)
  #GCD oracle calculation
  R_GCD_Oracle_400_Random_scheme_nothing = GCD_oracle_random_nothing(Data1_Random$real_graph, pulsar_gip_400_Random_gl,gcvec_extended_random1)
  
  R_GCD_Oracle_400_Random_scheme_random_noise = GCD_oracle_random_nothing(Data1_Random$real_graph, pulsar_gip_400_Random_gl,gcvec_extended_random2)
  
  R_GCD_Oracle_400_Random_scheme_random_vector = GCD_oracle_random(Data1_Random$real_graph, pulsar_gip_400_Random_gl,gcvec_extended_random3)
  
  R_GCD_Oracle_400_Random_scheme_both = GCD_oracle_random(Data1_Random$real_graph, pulsar_gip_400_Random_gl,gcvec_extended_random2)
  
  #GCD GIP calculation
  
  R_GCD_gip_hub_400_gl = Function_GCD_data(pulsar_gip_400_Random_gl, Data1_Random$real_graph)
  
  #HD oracle calculation
  R_HD_Hub_400 = Hamming_distance_path(pulsar_gip_400_Random_gl, Data1_Random$real_graph)
  
  #GCD oracle paths
  R_GCD_Oracle_results_strategy1 = append(R_GCD_Oracle_results_strategy1, list(R_GCD_Oracle_400_Random_scheme_nothing))
  R_GCD_Oracle_results_strategy2 = append(R_GCD_Oracle_results_strategy2, list(R_GCD_Oracle_400_Random_scheme_random_noise))
  R_GCD_Oracle_results_strategy3 = append(R_GCD_Oracle_results_strategy3, list(R_GCD_Oracle_400_Random_scheme_random_vector))
  R_GCD_Oracle_results_strategy4 = append(R_GCD_Oracle_results_strategy4, list(R_GCD_Oracle_400_Random_scheme_both))
  
  #GCD oracle better results
  R_GCD_Oracle_res_strategy1 = append(R_GCD_Oracle_res_strategy1, list(R_HD_Hub_400[which.min(R_GCD_Oracle_400_Random_scheme_nothing)]))
  R_GCD_Oracle_res_strategy2 = append(R_GCD_Oracle_res_strategy2, list(R_HD_Hub_400[which.min(R_GCD_Oracle_400_Random_scheme_random_noise)]))
  R_GCD_Oracle_res_strategy3 = append(R_GCD_Oracle_res_strategy3, list(R_HD_Hub_400[which.min(R_GCD_Oracle_400_Random_scheme_random_vector)]))
  R_GCD_Oracle_res_strategy4 = append(R_GCD_Oracle_res_strategy4, list(R_HD_Hub_400[which.min(R_GCD_Oracle_400_Random_scheme_both)]))
  
  #HD res
  
  R_Hamming_dist_res = append(R_Hamming_dist_res, list(R_HD_Hub_400))
  
  #Stars res
  
  R_v05 = 101-pulsar_gip_400_Random_gl$stars$opt.index
  R_v1 = 101-(which.min(pulsar_gip_400_Random_gl$stars$summary<=0.1)-1)
  R_Stars_01_res = append(R_Stars_01_res, list(R_HD_Hub_400[R_v1]))
  R_Stars_005_res = append(R_Stars_005_res, list(R_HD_Hub_400[R_v05]))
  
  #GCD GIP paths
  
  R_GCD_GIP_results_strategy1 = append(R_GCD_GIP_results_strategy1, list(R_GCD_gip_hub_400_gl$mean_gcd1))
  R_GCD_GIP_results_strategy2 = append(R_GCD_GIP_results_strategy2, list(R_GCD_gip_hub_400_gl$mean_gcd2))
  R_GCD_GIP_results_strategy3 = append(R_GCD_GIP_results_strategy3, list(R_GCD_gip_hub_400_gl$mean_gcd3))
  R_GCD_GIP_results_strategy4 = append(R_GCD_GIP_results_strategy4, list(R_GCD_gip_hub_400_gl$mean_gcd4))
  
  #GCD GIP res
  
  R_GCD_GIP_res_strategy1 = append(R_GCD_GIP_res_strategy1, list(R_HD_Hub_400[which.min(R_GCD_gip_hub_400_gl$mean_gcd1)]))
  R_GCD_GIP_res_strategy2 = append(R_GCD_GIP_res_strategy2, list(R_HD_Hub_400[which.min(R_GCD_gip_hub_400_gl$mean_gcd2)]))
  R_GCD_GIP_res_strategy3 = append(R_GCD_GIP_res_strategy3, list(R_HD_Hub_400[which.min(R_GCD_gip_hub_400_gl$mean_gcd3)]))
  R_GCD_GIP_res_strategy4 = append(R_GCD_GIP_res_strategy4, list(R_HD_Hub_400[which.min(R_GCD_gip_hub_400_gl$mean_gcd4)]))
  
}




#SAVE the Random results
vector_R_GCD_GIP_res_strategy1 = c()
vector_R_GCD_GIP_res_strategy2 = c()
vector_R_GCD_GIP_res_strategy3 = c()
vector_R_GCD_GIP_res_strategy4 = c()
for(i in R_GCD_GIP_res_strategy1){
  vector_R_GCD_GIP_res_strategy1 = c(vector_R_GCD_GIP_res_strategy1,i)
}

for(i in R_GCD_GIP_res_strategy2){
  vector_R_GCD_GIP_res_strategy2 = c(vector_R_GCD_GIP_res_strategy2,i)
}

for(i in R_GCD_GIP_res_strategy3){
  vector_R_GCD_GIP_res_strategy3 = c(vector_R_GCD_GIP_res_strategy3,i)
}

for(i in R_GCD_GIP_res_strategy4){
  vector_R_GCD_GIP_res_strategy4 = c(vector_R_GCD_GIP_res_strategy4,i)
}

vector_R_GCD_Oracle_res_strategy1 = c()
vector_R_GCD_Oracle_res_strategy2 = c()
vector_R_GCD_Oracle_res_strategy3 = c()
vector_R_GCD_Oracle_res_strategy4 = c()

for(i in R_GCD_Oracle_res_strategy1){
  vector_R_GCD_Oracle_res_strategy1 = c(vector_R_GCD_Oracle_res_strategy1,i)
}

for(i in R_GCD_Oracle_res_strategy2){
  vector_R_GCD_Oracle_res_strategy2 = c(vector_R_GCD_Oracle_res_strategy2,i)
}

for(i in R_GCD_Oracle_res_strategy3){
  vector_R_GCD_Oracle_res_strategy3 = c(vector_R_GCD_Oracle_res_strategy3,i)
}

for(i in R_GCD_Oracle_res_strategy4){
  vector_R_GCD_Oracle_res_strategy4 = c(vector_R_GCD_Oracle_res_strategy4,i)
}



vector_R_Stars_01_res = c()
vector_R_Stars_005_res = c()

for(i in R_Stars_01_res){
  vector_R_Stars_01_res = c(vector_R_Stars_01_res,i)
}

for(i in R_Stars_005_res){
  vector_R_Stars_005_res = c(vector_R_Stars_005_res,i)
}


Best_HD_R = c()

for(i in R_Hamming_dist_res){
  R = min(i)
  Best_HD_R = c(Best_HD_R,R)
}


res_tot_random = cbind(vector_R_GCD_GIP_res_strategy1, vector_R_GCD_GIP_res_strategy2, vector_R_GCD_GIP_res_strategy3, vector_R_GCD_GIP_res_strategy4,
                       vector_R_GCD_Oracle_res_strategy1, vector_R_GCD_Oracle_res_strategy2, vector_R_GCD_Oracle_res_strategy3, vector_R_GCD_Oracle_res_strategy4,
                       vector_R_Stars_005_res, vector_R_Stars_01_res, Best_HD_R)

res_tot_random_df = as.data.frame(res_tot_random)




res_HD_1 = c()
for(i in R_Hamming_dist_res){
  res_HD_1 = cbind(res_HD_1, i)
}

res_HD_1_df = as.data.frame(res_HD_1)
colnames(res_HD_1_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_Oracle_S1 = c()
for(i in R_GCD_Oracle_results_strategy1){
  res_GCD_Oracle_S1 = cbind(res_GCD_Oracle_S1, i)
}

res_GCD_Oracle_S1_df = as.data.frame(res_GCD_Oracle_S1)
colnames(res_GCD_Oracle_S1_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_Oracle_S2 = c()
for(i in R_GCD_Oracle_results_strategy2){
  res_GCD_Oracle_S2 = cbind(res_GCD_Oracle_S2, i)
}

res_GCD_Oracle_S2_df = as.data.frame(res_GCD_Oracle_S2)
colnames(res_GCD_Oracle_S2_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_Oracle_S3 = c()
for(i in R_GCD_Oracle_results_strategy3){
  res_GCD_Oracle_S3 = cbind(res_GCD_Oracle_S3, i)
}

res_GCD_Oracle_S3_df = as.data.frame(res_GCD_Oracle_S3)
colnames(res_GCD_Oracle_S3_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_Oracle_S4 = c()
for(i in R_GCD_Oracle_results_strategy4){
  res_GCD_Oracle_S4 = cbind(res_GCD_Oracle_S4, i)
}

res_GCD_Oracle_S4_df = as.data.frame(res_GCD_Oracle_S4)
colnames(res_GCD_Oracle_S4_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")



res_GCD_GIP_S1 = c()
for(i in R_GCD_GIP_results_strategy1){
  res_GCD_GIP_S1 = cbind(res_GCD_GIP_S1, i)
}

res_GCD_GIP_S1_df = as.data.frame(res_GCD_GIP_S1)
colnames(res_GCD_GIP_S1_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_GIP_S2 = c()
for(i in R_GCD_GIP_results_strategy2){
  res_GCD_GIP_S2 = cbind(res_GCD_GIP_S2, i)
}

res_GCD_GIP_S2_df = as.data.frame(res_GCD_GIP_S2)
colnames(res_GCD_GIP_S2_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")



res_GCD_GIP_S3 = c()
for(i in R_GCD_GIP_results_strategy3){
  res_GCD_GIP_S3 = cbind(res_GCD_GIP_S3, i)
}

res_GCD_GIP_S3_df = as.data.frame(res_GCD_GIP_S3)
colnames(res_GCD_GIP_S3_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")


res_GCD_GIP_S4 = c()
for(i in R_GCD_GIP_results_strategy4){
  res_GCD_GIP_S4 = cbind(res_GCD_GIP_S4, i)
}

res_GCD_GIP_S4_df = as.data.frame(res_GCD_GIP_S4)
colnames(res_GCD_GIP_S4_df) = c("sim1", "sim2", "sim3", "sim4", "sim5","sim6", "sim7", "sim8", "sim9", "sim10")








write.csv(res_tot_random_df, "tot_results_random_simulation.csv")

write.csv(res_HD_1_df, "simulations_HD_random.csv")

write.csv(res_GCD_Oracle_S1_df, "simulations_GCD_Oracle_S1_random.csv")
write.csv(res_GCD_Oracle_S2_df, "simulations_GCD_Oracle_S2_random.csv")
write.csv(res_GCD_Oracle_S3_df, "simulations_GCD_Oracle_S3_random.csv")
write.csv(res_GCD_Oracle_S4_df, "simulations_GCD_Oracle_S4_random.csv")

write.csv(res_GCD_GIP_S1_df, "simulations_GCD_GIP_S1_random.csv")
write.csv(res_GCD_GIP_S2_df, "simulations_GCD_GIP_S2_random.csv")
write.csv(res_GCD_GIP_S3_df, "simulations_GCD_GIP_S3_random.csv")
write.csv(res_GCD_GIP_S4_df, "simulations_GCD_GIP_S4_random.csv")


