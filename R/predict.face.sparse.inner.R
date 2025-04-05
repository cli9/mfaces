predict.face.sparse.inner <- function(object,newdata,calculate.scores=T,...){
  
  ## check inputs
  if(!inherits(object, "face.sparse")) stop("'fit' has to be a face.sparse object")
  check.data(newdata,type="predict") 

  subj.pred = newdata$subj
  subj_unique.pred = unique(subj.pred)
  y.pred = newdata$y
  
  Theta = object$Theta
  knots = object$knots
  npc = object$npc
  p = object$p
  
  center = object$center
  #const.var.error = object$const.var.error
  fit_mean = object$fit_mean
  #fit_var_error = object$fit_var_error
  sigma2 = object$sigma2
  
  mu.pred <- rep(0,nrow(newdata))
  var.error.pred <- rep(max(sigma2,0.000001),nrow(newdata))
  # cov.pred <-  matrix(0,length(mu.pred),length(mu.pred))
  se.pred <- rep(0, nrow(newdata))
    
  if(center) {mu.pred <- predict(fit_mean,newdata$argvals) }
  #if(!const.var.error){
  #  var.error.pred <- predict.pspline(fit_var_error,newdata$argvals)
  #  var.error.pred <-  sapply(var.error.pred,function(x) max(x,0))### added Sept 8, 2015 by Luo
  #}

  
  Bnew = spline.des(knots=knots, x=object$argvals.new, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  scores = list(subj=subj_unique.pred,
                scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
                )
  
  Bi_pred <- vector("list", length(subj_unique.pred))
  Bi <- vector("list", length(subj_unique.pred))
  Pi_pred <- vector("list", length(subj_unique.pred))
  Pi <- vector("list", length(subj_unique.pred))
  
  # interpolating
  h <- function(t1){
    est <- sapply(1:ncol(object$eigenfunctions),
                  function(x) spline(object$argvals.new,
                                     object$eigenfunctions[,x],xout=t1)$y)
    return(est)
  }
  
  for(i in 1:length(subj_unique.pred)){
    sel.pred = which(subj.pred==subj_unique.pred[i])
    pred.points <- newdata$argvals[sel.pred]
    mu.predi <- mu.pred[sel.pred]
    var.error.predi <- var.error.pred[sel.pred]
    
    y.predi = y.pred[sel.pred] - mu.predi
    sel.pred.obs = which(!is.na(y.predi))
    obs.points <- pred.points[sel.pred.obs]
     
    if(!is.null(obs.points)){
      var <- mean(var.error.predi[sel.pred.obs])
      if(var==0&length(sel.pred.obs) < npc)
        stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
      B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      Bi_pred[[i]] <- B3i.pred   
      Pi_pred[[i]] <- Matrix(h(pred.points)) # <--- 
      
      Bi[[i]] <- B3i 
      Pi[[i]] <- Matrix(h(obs.points)) # <---
      Chati = tcrossprod(B3i%*%Theta,B3i)
      
      if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
      if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
      Vi.inv = as.matrix(solve(Chati + Ri))
      Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta,B3i.pred))
      Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
      ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
      scores$u[i,] = as.vector(ui)
      y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
      temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
      if(length(sel.pred.obs) > 1){
        se.pred[sel.pred] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))
      }
      if(length(sel.pred.obs) == 1){
        se.pred[sel.pred] = sqrt(Vi.pred[1, 1] - Vi.inv[1, 1]*temp%*%t(temp))
      }
      ## predict scores
      if(object$calculate.scores==TRUE || calculate.scores==TRUE){ 
        temp = matrix(t(object$eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(object$eigenfunctions[,1]^2)
        temp = as.matrix(temp)
        scores$scores[i,1:npc] = temp[,1]
      }
    }
  }
  
  B = spline.des(knots=knots, 
                 x=newdata[!is.na(newdata$y), 'argvals']
                 ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
  P = Matrix(h(newdata$argvals))  # <---
  
  return(list(object=object,newdata=newdata,y.pred = y.pred,
              mu.pred=mu.pred,var.error.pred=var.error.pred,
              rand_eff = scores, 
              # cov.pred = cov.pred, 
              se.pred = se.pred,
              Chat.pred = Chat.pred,Chat.diag.pred = diag(Chat.pred), 
              Bi = Bi, Bi_pred = Bi_pred, B = B, 
              Pi = Pi, Pi_pred = Pi_pred, P = P))  
}
