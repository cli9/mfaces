predict.mface.sparse <- function(object,newdata, 
                                 calculate.scores=T, ...){
  
  ## check inputs
  if(!inherits(object, "mface.sparse")) stop("'fit' has to be a mface.sparse object")
  p <- length(newdata)
  Bnew <- vector("list", p)
  var.error.pred <- vector("list", p)
  Bi <- vector("list", p)
  Bi_pred <- vector("list", p)
  B <- vector("list", p)
  mu.pred <- vector("list", p)
  uy.pred <- vector("list", p)
  uscores.pred <- vector("list", p)
  ucov.pred <- vector("list", p)
  for(k in 1:p){
    temp <- newdata[[k]]
    pred <- predict.face.sparse.inner(object$fit[[k]], temp)
    # if(type==1){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    Bi[[k]] <- lapply(1:length(pred$Bi),function(x) pred$Bi[[x]]%*%object$G_invhalf)
    Bi_pred[[k]] <- lapply(1:length(pred$Bi_pred),function(x) pred$Bi_pred[[x]]%*%object$G_invhalf)
    B[[k]] <- pred$B %*% object$G_invhalf
    # }else{
    #   Bnew[[k]] <- Matrix(object$fit[[k]]$eigenfunctions)
    #   Bi[[k]] <- pred$Pi
    #   Bi_pred[[k]] <- pred$Pi_pred
    #   B[[k]] <- pred$P
    # }
    
    var.error.pred[[k]] <- pred$var.error.pred
    mu.pred[[k]] <- pred$mu.pred
    uy.pred[[k]] <- pred$y.pred
    uscores.pred[[k]] <- pred$rand_eff$scores
    ucov.pred[[k]] <- pred$cov.pred
  }

  subj.pred = lapply(newdata, function(x) {x$subj})
  subj_unique.pred = unique(subj.pred[[1]])
  y.pred = lapply(newdata, function(x) {x[, "y"]})
  se.pred = lapply(y.pred, function(x) {0*x})
  
  Theta = object$Theta
  npc = object$npc
  
  Bnew <- do.call(bdiag, Bnew)
  Bi <- unlist(Bi)
  Bi_pred <- unlist(Bi_pred)
  B <- do.call(bdiag, B) 
  uscores.pred <- unlist(uscores.pred)
  # cov.pred <- matrix(0,length(unlist(mu.pred)),length(unlist(mu.pred)))
  
  scores = list(subj=subj_unique.pred,
                scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta)),
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
      Chati = tcrossprod(B3i%*%Theta,B3i)
      
      if(sum(len_sel.pred.obs)==1) Ri = unlist(var_sel)
      if(sum(len_sel.pred.obs)>1)  Ri = diag(unlist(var_sel))
      Vi.inv = as.matrix(solve(Chati + Ri))
      Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta,B3i.pred))
      Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
      y.predi_sel = unname(unlist(mapply(function(y.predi, sel.pred.obs){
        y.predi[sel.pred.obs]}, y.predi = y.predi, 
        sel.pred.obs = sel.pred.obs, SIMPLIFY = F)))
      ui =tcrossprod(Theta,B3i)%*%Vi.inv %*% y.predi_sel
      Si = tcrossprod(Theta,B3i)%*%Vi.inv%*%t(tcrossprod(Theta,B3i))
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
      temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
      # idx_1 = lapply(sel.pred, function(x){x[1]-1})
      # idx_n = lapply(sel.pred, function(x){x[length(x)]})
      # idx = (sum(unlist(idx_1))+1):sum(unlist(idx_n))
      for(j in 1:p){
        y.pred[[j]][sel.pred[[j]]] <- temp0[ind_temp[[j]]]
        if(sum(len_sel.pred.obs) > 1){
          # cov.pred[idx,idx] = Vi.pred - temp%*%Vi.inv%*%t(temp)
          se.pred[[j]][sel.pred[[j]]] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))[ind_temp[[j]]]
        }
        
        if(sum(len_sel.pred.obs) == 1){
          # cov.pred[idx,idx] = Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp)
          se.pred[[j]][sel.pred[[j]]] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
        }
      }
      ## predict scores
      if(calculate.scores==TRUE){ 
        temp = matrix(t(object$U),nrow=npc)%*%ui
        temp = as.matrix(temp)
        scores$scores[i,1:npc] = temp[,1]
        scores$Sig[[i]] = diag(object$eigenvalues) - 
          (matrix(t(object$U),nrow=npc)) %*% as.matrix(Si) %*% matrix(object$U,ncol=npc)
      }
    }
  }
  
  Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
  
  
  return(list(object=object,newdata=newdata,B=B,
              y.pred = y.pred,mu.pred=mu.pred,var.error.pred=var.error.pred,
              rand_eff = scores, 
              # cov.pred = cov.pred,
              se.pred = se.pred, uscores.pred = uscores.pred,
              Chat.diag.pred = diag(Chat.pred),
              ucov.pred = ucov.pred, uy.pred = uy.pred))  
}
