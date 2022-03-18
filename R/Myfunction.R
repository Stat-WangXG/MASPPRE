#??Predicting the survival probability of individual by using the model averaging method for survival data.
#'
#'This function is applied to the partly linear additive Cox model, and offers the weight estimators of model averaging and the prediction of survival probability of individual given the time point for right censoring data, where the weight selection criteria employ the inverse probability of censoring
#'weights and the non-parametric parts of this model are approximated by B-splines method.
#'
#'@param t  numeric, the given time point.
#'@param sim.data  a data.frame, the column names of the data.frame are U.1, U.2, U.3 and U.4 (Continuous covariate), X (The binary covariate), time (Observation time), status (Censor indicator).
#'@param p  numeric, the dimension of continuous covariate.
#'@param l  numeric, the dimension of discrete covariate.
#'@param h  numeric, the number of splines basis functions.
#'@param alpha numeric, The partition ratio of training set and test set
#'
#'@return MASPPRE returns a list containing the following components:
#'@return \item{weights}{the estimators of weights.}
#'@return \item{SP_train}{the prediction of survival probability in the train sets.}
#'@return \item{SP_test}{the prediction of survival probability in the test sets.}
#'@export
#'
#'@examples
#'data("sim.data")
#'h <- 5
#'p <- 4
#'l <- 1
#'alpha <- 0.7
#'t <- mean(sim.data$time)
#'## estimate weights and predict the survival probability
#'# MASPPRE(t, sim.data, p, l, h, alpha)

MASPPRE <- function(t, sim.data, p, l, h, alpha ){

  #### an approximation with B-spline

  e1 <- p+l

  bsdata<-list()

  for ( a in 1:p) {

    bsp <- splines::bs(sim.data[,a],df=h)
    bsdata[[a]]<-bsp

  }


  data<-NULL
  for(a1 in 1:p){

    data <- cbind(data,bsdata[[a1]])

  }

  colnames(data) <- paste('p',1:(h*p),sep='.')
  data <- cbind(sim.data[,1:e1], data, time=sim.data$time, status=sim.data$status)

  train_sub <- sample(nrow(data), alpha*nrow(data))
  train.data <- data[train_sub,]
  test.data <- data[-train_sub,]

  n1 <- dim(train.data)[1]
  n2 <- dim(test.data)[1]

  time <- train.data$time
  status <- train.data$status

  time2 <- test.data$time
  status2 <- test.data$status

  cstatus <- 1-status


  Est.eff <- matrix(0,nrow =n1,ncol = p)       ###the effects in the train set
  pEst.eff <- matrix(0,nrow =n2,ncol = p)    ###the effects in the test set

  for (m in 1:p) {

    formula_m=paste("survival::Surv(time,status)~",paste(colnames(data)[1:e1][-m],sep="",collapse = "+"),"+",paste(colnames(data)[(e1+(m-1)*h+1):(e1+(m)*h)],sep="",collapse = "+"),sep = "")
    formula_m=as.formula(formula_m)
    fit_m <- survival::coxph(formula_m,data=train.data)
    coef_m <- as.matrix(fit_m$coefficients)

    rc_m=cbind(train.data[,1:e1][,-m],train.data[,(e1+(m-1)*h+1):(e1+(m)*h)])
    rc_m=as.matrix(rc_m)
    esti_m <- rc_m%*%coef_m
    Est.eff[,m] <- esti_m

    prc_m=cbind(test.data[,1:e1][,-m],test.data[,(e1+(m-1)*h+1):(e1+(m)*h)])
    prc_m=as.matrix(prc_m)
    pesti_m <- prc_m%*%coef_m
    pEst.eff[,m] <- pesti_m

  }


  ### survival function
  SF <- function(tm,Phi,time,status){

    Phiw <- as.matrix(Phi)
    H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=tm))
    Stx <- exp( -H0t*exp(Phiw) )
    return(Stx)

  }

  St <- function(tmm)
  {
    return(SF(tmm,Est.eff[,i],time,status))
  }

  #####the estimated survival probabilities of submodels
  ESSP <- list()
  for(i in 1:p)
  {
    essp <- matrix(unlist(lapply(time,St)), nrow=n1,ncol=n1) #matrix:the first col is the vector when tm takes tm[1]
    ESSP[[i]] <- essp
  }

  formula_cens=paste("survival::Surv(time,cstatus)~",paste(colnames(data)[1:e1],sep="",collapse = "+"))
  formula_cens=as.formula(formula_cens)
  fit_cens <- survival::coxph(formula_cens,data=train.data,x=TRUE,y=TRUE,iter.max=200)
  SP_cens <- pec::predictSurvProb(fit_cens,newdata=train.data,times=time)###predict the probablity of censor data in the train set


  ####model averaging  weights

  GYZ <- diag(SP_cens)

  D_ma <- list()
  d_ma <- list()

  for( im in 1:n1){
    ### Compute I(Y>t) and W(t,G,X)

    It <- 1*(train.data$time>time[im])
    WGt <- ((1-It)*train.data$status)/GYZ+It/SP_cens[,im]
    WGt[is.nan(WGt)] <- 0

    esSP <- matrix(0,n1,p)
    for (im1 in 1:p) {

      esSP1 <- ESSP[[im1]][,im]
      esSP[,im1] <- esSP1

    }

    Dmat1_ma <- t(esSP)%*%diag(WGt)%*%esSP
    dvec1_ma <- as.vector(t(esSP)%*%diag(WGt) %*% It )
    D_ma[[im]]<-Dmat1_ma
    d_ma[[im]]<-dvec1_ma

  }

  S1_ma <- 0
  H_ma <- 0

  for(i2 in 1:n1){

    S1_ma <- S1_ma+D_ma[[i2]]
    H_ma <- H_ma+d_ma[[i2]]

  }

  Dmat_ma <- S1_ma
  dvec_ma <- H_ma
  Amat_ma <- cbind(rep(1,p), diag(p))
  bvec_ma <- c(1, rep(0,p))
  sol_ma  <- quadprog::solve.QP(Dmat_ma,dvec_ma,Amat_ma,bvec_ma,meq=1) # sol$solution/unconstrained.solution
  ### collect the estimated w for every tm's in each column
  w_ma <- sol_ma$solution


  Est.sp <- matrix(0,nrow =n1,ncol = p)       ###the estimated survival probability in the train set given the time point
  pEst.sp <- matrix(0,nrow =n2,ncol = p)      ###the estimated survival probability effects in the test set given the time point

  for (m2 in 1:p) {

    est.sp <- SF(t,Est.eff[,m2],time,status)
    pest.sp <- SF(t,pEst.eff[,m2],time2,status2)
    Est.sp[,m2] <- est.sp
    pEst.sp[,m2] <- pest.sp
  }

  ESP_ma <- Est.sp%*%w_ma
  PESP_ma <- pEst.sp%*%w_ma

  result <- list(weight=w_ma, SP_train=ESP_ma, SP_test=PESP_ma )

  return(result)
}

