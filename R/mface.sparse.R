mface.sparse <- function(data, newdata = NULL, center = TRUE, argvals.new = NULL, 
                         knots = 7, knots.option="equally-spaced", 
                         p = 3, m = 2, lambda = NULL, lambda_mean = NULL, 
                         lambda_bps = NULL, 
                         search.length = 14, lower = -3, upper = 10, 
                         calculate.scores = FALSE, pve = 0.99)
{
  # require(face)
  # require(matrixcalc)
  # require(Matrix)
  # require(splines)
  # require(mgcv)
  
  #########################
  ####step 0: read in data
  #########################
  
  if(!is.list(data)){
    stop("'data' should be a list of data frames with three variables:argvals,subj and y")
  }
  tnew <- argvals.new
  if(is.null(tnew)){
    if(any(unlist(lapply(data, function(x){is.null(x$argvals)})))){
      stop("'data' should be a list of data frames with three variables:argvals,subj and y")
    }
    temp_t = unlist(lapply(data, function(x){range(x$argvals)}))
    tnew <- seq(min(temp_t), max(temp_t), length=100)
  }
  p.m = p ## this is splines degree
  p <- length(data) ## dimension of the multivariate functional data
  
  object <- vector("list", p)
  Theta0 <- vector("list", p)
  Bnew <- vector("list", p)
  C <- vector("integer", p)
  # sigma2 <- vector("numeric", p)
  var.error.pred <- vector("list", p)
  var.error.hat <- vector("list", p)
  var.error.new <- vector("list", p)
  mu.pred <- vector("list", p)
  Bi <- vector("list", p)
  Bi_pred <- vector("list", p)
  B <- vector("list", p)
  
  W <- data
  
  #############################################
  ####step 1: univariate fpca for each response
  #############################################
  
  for(k in 1:p){
  temp <- data[[k]]
  temp_new = NULL
  if(!is.null(newdata)) {
    temp_new <- newdata[[k]]
  }
  object[[k]] <- face.sparse.inner(data = temp, newdata = temp_new, center = center, 
                 argvals.new = tnew, knots = knots, 
                 knots.option = knots.option, p = p.m, m = m, 
                 lambda = lambda, lambda_mean = lambda_mean, 
                 search.length = search.length, lower = lower, upper = upper, 
                 calculate.scores = calculate.scores, pve = pve)
  Theta0[[k]] <- object[[k]]$Theta0
  if(k==1){
    G <- crossprod(object[[k]]$Bnew) / nrow(object[[k]]$Bnew)
    eig_G <- eigen(G, symmetric = T)
    G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
    G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
  }
  Bnew[[k]] <- object[[k]]$Bnew %*% G_invhalf
  C[k] <- ncol(object[[k]]$Bnew)
  # sigma2[k] <- object[[k]]$sigma2
  var.error.hat[[k]] <- object[[k]]$var.error.hat
  var.error.new[[k]] <- object[[k]]$var.error.new
  if(!is.null(newdata) || calculate.scores == T) {
    mu.pred[[k]] <- object[[k]]$mu.pred
    var.error.pred[[k]] <- object[[k]]$var.error.pred
    Bi[[k]] <- lapply(1:length(object[[k]]$Bi),function(x) object[[k]]$Bi[[x]]%*%G_invhalf)
    Bi_pred[[k]] <- lapply(1:length(object[[k]]$Bi_pred),function(x) object[[k]]$Bi_pred[[x]]%*%G_invhalf)
    B[[k]] <- object[[k]]$B %*% G_invhalf
  }
  if(center){
    W[[k]][, "y"] <- W[[k]][, "y"] - object[[k]]$mu.hat # demean each response
    }
  }
  names(object) <- names(data)
  
  ###############################################
  ####step 2: compute empirical cross-covariances
  ###############################################

  xprod <- vector("list", p*(p-1)/2)
  N2 <- vector("list", p*(p-1)/2)
  usubj <- unique(W[[1]]$subj)
  for(i in 1:length(usubj)){
    crossCov <- vector("list", p*(p-1)/2)
    ind <- 1
    for(k in 1:(p-1))
      for(j in (k+1):p){
        temp_k <- W[[k]][which(W[[k]]$subj == usubj[i]),]
        temp_j <- W[[j]][which(W[[j]]$subj == usubj[i]),]
        N2[[ind]][i] <- nrow(temp_k) * nrow(temp_j)
        crossCov[[ind]] <- as.vector(temp_k[,"y"]%*%t(temp_j[,"y"]))
        xprod[[ind]][[i]] <- data.frame(arg1 = rep(temp_k$argvals, nrow(temp_j)),
                                        arg2 = rep(temp_j$argvals, each = nrow(temp_k)),
                                        crossCov = crossCov[[ind]])  
        ind <- ind + 1
      }
  }

  ####################################
  ####step 3: smooth cross-covariance
  ####################################

  # smooth.cov <- vector("list",p*(p-1)/2)
  Theta <- vector("list", p*(p-1)/2)
  bps.lambda <- matrix(NaN, p*(p-1)/2,2)
  for(k in 1:(p*(p-1)/2)) {
  	temp <- do.call(rbind, xprod[[k]])
  	fit <-  bps(temp$crossCov, temp$arg1, temp$arg2, knots=knots, p = p.m, m = m, 
  	            lower = lower, upper = upper, search.length = search.length,
  	            lambda = lambda_bps, knots.option=knots.option, N2 = N2[[k]])
  	Theta[[k]] <- fit$Theta
  	bps.lambda[k,] <- fit$lambda
  }
  
 Theta_all <- do.call(bdiag, Theta0)
 C_idx <- c(0, cumsum(C))
 ind <- 1
 for(k in 1:(p-1)){
   for(j in (k+1):p){
     sel1 <- C_idx[k] + 1:C[k]
     sel2 <- C_idx[j] + 1:C[j]
     Theta_all[sel1,sel2] <- Theta[[ind]]
     Theta_all[sel2,sel1] <- t(Theta[[ind]])
     ind <- ind + 1
   }
 }
 Theta_all <- as.matrix(Theta_all)
 
 
 G_halfmat <-do.call(bdiag, lapply(1:p, function(x) G_half))
 Theta_all <- G_halfmat %*% Theta_all %*% G_halfmat
 
 ####################################################
 ####step 4: calculate estimated covariance function
 ####################################################
 
 ##make covariance matrix psd
 # Eig <- eigen(cov_all)
 Eig <- eigen(Theta_all)
 sel <- which(Eig$values>0)
 if(length(sel)>1) Theta_all <- matrix.multiply(Eig$vectors[,sel], Eig$values[sel])%*%t(Eig$vector[,sel])
 if(length(sel)==1) Theta_all <- Eig$vectors[,sel]%*%t(Eig$vector[,sel])*Eig$values[sel]
   
 Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
 Chat <- Bnew%*%Theta_all%*%t(Bnew)
 
 Chat_diag = as.vector(diag(Chat))  
 Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
 
 npc <- which.max(cumsum(Eig$values[sel])/sum(Eig$values[sel])>pve)[1]
 
 eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                         ncol=min(npc, p*length(tnew))) #/ sqrt((max(tnew) - min(tnew))/(length(tnew) - 1))
 eigenvalues = Eig$values[1:min(npc, p*length(tnew))] #* ((max(tnew) - min(tnew))/(length(tnew) - 1))

 ##############################
 ####step 5: calculate variance
 ##############################
 
 # var.error.hat <- do.call(cbind, var.error.hat)
 var.error.new <- do.call(cbind, var.error.new)
 
 Chat.raw.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta_all),Bnew)) + diag(var.error.new) 
 Chat.raw.diag.new = as.vector(diag(Chat.raw.new)) 
 Cor.raw.new = diag(1/sqrt(Chat.raw.diag.new))%*%Chat.raw.new%*%diag(1/sqrt(Chat.raw.diag.new))
 
 #######################
 ####step 6: prediction
 #######################
 if(is.null(newdata) && calculate.scores==T){
   newdata = data
 }
 if(!is.null(newdata)){
 
   subj.pred = lapply(newdata, function(x) {x$subj})
   subj_unique.pred = unique(subj.pred[[1]])
   y.pred = lapply(newdata, function(x) {x[, "y"]})
   se.pred = lapply(y.pred, function(x) {0*x})
   
   Bi <- unlist(Bi)
   Bi_pred <- unlist(Bi_pred)
   B <- do.call(bdiag, B) 
   
   scores = list(subj=subj_unique.pred,
                 scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                 u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta_all)),
                 Sig = list())
   
   for(i in 1:length(subj_unique.pred)){
     sel.pred = lapply(subj.pred, function(x){
       which(x==subj_unique.pred[i])})
     pred.points <- mapply(function(newdata, sel.pred) {
       newdata$argvals[sel.pred]}, newdata = newdata, sel.pred = sel.pred, 
       SIMPLIFY = F)
     mu.predi <- mapply(function(mu.pred, sel.pred) {
       mu.pred[sel.pred]}, mu.pred = mu.pred, sel.pred = sel.pred, 
       SIMPLIFY = F)
     var.error.predi <- mapply(function(var.error.pred, sel.pred) {
       var.error.pred[sel.pred]}, var.error.pred = var.error.pred,
       sel.pred = sel.pred, SIMPLIFY = F)
     
     y.predi = mapply(function(y.pred, sel.pred, mu.predi) {
       y.pred[sel.pred] - mu.predi}, y.pred = y.pred, 
       sel.pred = sel.pred, mu.predi = mu.predi, SIMPLIFY = F)
     sel.pred.obs = lapply(y.predi, function(x){
       unname(which(!is.na(x)))})
     obs.points <- mapply(function(pred.points, sel.pred.obs) {
       pred.points[sel.pred.obs]}, pred.points = pred.points, 
       sel.pred.obs = sel.pred.obs, SIMPLIFY = F)
     
     if(!any(unlist(lapply(obs.points, function(x){is.null(x)}))==T)){
       var_sel <- mapply(function(var.error.predi, sel.pred.obs) {
        var.error.predi[sel.pred.obs]
       }, var.error.predi = var.error.predi, 
       sel.pred.obs = sel.pred.obs, SIMPLIFY = F)
       var = unlist(lapply(var_sel, function(x) {mean(x)}))
       len_sel.pred.obs <- unlist(lapply(sel.pred.obs, function(x){
         length(x)}))
       if(sum(var)==0&sum(len_sel.pred.obs) < npc)
         stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
              cannot be estimated.")
       
       B3i.pred = Bi_pred[[i]]
       B3i = Bi[[i]]
       for(j in 1:(p-1)){
         B3i.pred = bdiag(B3i.pred, Bi_pred[[i+j*length(subj_unique.pred)]])
         B3i = bdiag(B3i, Bi[[i+j*length(subj_unique.pred)]])
       }
       Chati = tcrossprod(B3i%*%Theta_all,B3i)
       
       if(sum(len_sel.pred.obs)==1) Ri = unlist(var_sel)
       if(sum(len_sel.pred.obs)>1)  Ri = diag(unlist(var_sel))
       Vi.inv = as.matrix(solve(Chati + Ri))
       Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta_all,B3i.pred))
       Hi = as.matrix(B3i.pred%*%tcrossprod(Theta_all,B3i)%*%Vi.inv)
       y.predi_sel = unname(unlist(mapply(function(y.predi, sel.pred.obs){
         y.predi[sel.pred.obs]}, y.predi = y.predi, 
         sel.pred.obs = sel.pred.obs, SIMPLIFY = F)))
       ui =tcrossprod(Theta_all,B3i)%*%Vi.inv %*% y.predi_sel
       Si = tcrossprod(Theta_all,B3i)%*%Vi.inv%*%t(tcrossprod(Theta_all,B3i))
       scores$u[i,] = as.vector(ui)
       temp0 = as.numeric(Hi%*%y.predi_sel) + unlist(mu.predi)
       len_sel.pred = unlist(lapply(sel.pred, function(x){length(x)}))
       ind_sel.pred = unname(cumsum(len_sel.pred))
       ind_temp <- vector("list", p)
       for(j in 1:p){
         if(j==1){
           ind_temp[[j]] = 1:ind_sel.pred[1]
         }else{
           ind_temp[[j]] = (1+ind_sel.pred[j-1]):ind_sel.pred[j]
         }
       }
       temp = as.matrix(B3i.pred%*%tcrossprod(Theta_all,B3i))
       for(j in 1:p){
         y.pred[[j]][sel.pred[[j]]] <- temp0[ind_temp[[j]]]
         if(sum(len_sel.pred.obs) >1){
           se.pred[[j]][sel.pred[[j]]] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))[ind_temp[[j]]]
         }
         if(sum(len_sel.pred.obs) ==1){
           se.pred[[j]][sel.pred[[j]]] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
         }
       }

       ## predict scores
       if(calculate.scores==TRUE){ 
         temp = matrix(t(Eig$vectors[,1:npc]),nrow=npc)%*%ui
         temp = as.matrix(temp)
         scores$scores[i,1:npc] = temp[,1]
         scores$Sig[[i]] = diag(eigenvalues) - 
           (matrix(t(Eig$vectors[,1:npc]),nrow=npc)) %*% as.matrix(Si) %*% matrix(Eig$vectors[,1:npc],ncol=npc)
       }
     }
   }
   
   Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all), B)))
   
 }
 if(is.null(newdata)){
   y.pred=NULL
   mu.pred = NULL
   var.error.pred = NULL
   Chat.diag.pred = NULL
   # cov.pred = NULL
   se.pred = NULL
   scores=NULL
   
 }
 
 res <- list(fit = object, Theta = Theta_all, 
             Chat.new = Chat, Cor.new = Cor, 
             npc = npc, eigenfunctions = eigenfunctions, 
             eigenvalues = eigenvalues, 
             var.error.hat = var.error.hat, 
             var.error.new = var.error.new, Chat.raw.diag.new = Chat.raw.diag.new, 
             Cor.raw.new = Cor.raw.new,
             y.pred = y.pred, mu.pred = mu.pred, var.error.pred = var.error.pred, 
             Chat.diag.pred = Chat.diag.pred, 
             se.pred = se.pred, rand_eff = scores,
             G_invhalf = G_invhalf, bps.lambda = bps.lambda, 
             U = Eig$vectors[,1:npc], 
             argvals.new = tnew,
             center=center, knots=knots, knots.option = knots.option,
             p = p.m, m = m,
             lower=lower,upper=upper,search.length=search.length,
             calculate.scores=calculate.scores,pve=pve)
 class(res) <- "mface.sparse"
 return(res)
  
}
