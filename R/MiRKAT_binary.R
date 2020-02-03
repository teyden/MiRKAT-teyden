MiRKAT_binary = function(y, X = NULL, Ks, family= "binomial", nperm = 999, method = "davies"){
  n <- length(y)
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
                  
  qX1 <- qr(X1)
  
  ## Take care of aliased variables and pivoting in rhs
  X1 <- X1[, qX1$pivot, drop=FALSE]
  X1 <- X1[, 1:qX1$rank, drop=FALSE]
  options(warn=2)  # make sure this model is correct
  mod <- glm(y ~ X1-1, family = binomial)
  options(warn=1)
  
  px  = NCOL(X1)
  mu  = mod$fitted.values
  residual = y - mu  
  
  w   = mu*(1-mu)
  D0  =  sqrt(w)  
  DX12 = D0 * X1
  P0 = diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)

  if (method == "davies"){    
    if (n < 50){
      warning("For binary outcome and n < 50, p-value using davies method can be inaccurate at tails, permutation is recommended.")
    }
    S = sapply(Ks, getIndivP_binary, residual,  D0, px, P0)
    ps = as.numeric(unlist(S[3,]))
    if (length(Ks) ==1){
      return(indivP = ps)
    }
    eP0 = c(rep(1, n-px), rep(0, px))
    Qs = unlist(S[1,])
 
    q_sim = get_q_sim(Ks, n, nperm, residual)
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  

    return(list(indivP = ps , omnibus_p = p_final))    
  
  } else if (method == "moment"){
    S = sapply(Ks, getIndivP_hm, residual, mu, D0, P0)
    ps = as.numeric(unlist(S[1,]))
    Qs = unlist(S[2,])
    
    if (length(Ks) == 1){
      return(indivP = ps)
    }

    q_sim = get_q_sim(Ks, n, nperm, residual)
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  

    return(list(indivP = ps, omnibus_p = p_final))

  } else if (method == "permutation"){
    Qs = lapply(Ks, getQ, residual, s2 = 1)
    q_sim = get_q_sim(Ks, n, nperm, residual)

    if (length(Ks) == 1){
     p_perm = (sum(q_sim > Qs)+ 1)/(nperm + 1) 
     return(indivP = p_perm)
    }
  
    q_sim = t(q_sim)
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1)   
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm + 1)
    
    return(list(indivP = p_perm , omnibus_p = p_final))

  }
}

get_q_sim <- function(Ks, n, nperm, residual) {
  q_sim = sapply(1:nperm, function(i){
    ind <- sample(n)
    p1 = sapply(1:length(Ks), function(j) {
      Q1 = as.numeric(residual %*% Ks[[j]][ind, ind] %*% residual) # the adjusted is zero in this case 
      return(Q1)
    })  
  })
  return(q_sim)
}