# Databricks notebook source
#CRAN Libraries needed on cluster
#partykit
#tree
#RobustGaSP
#doParallel
#DiceDesign
#gridExtra
#tgp

library(ggplot2)
library(gridExtra)
library(randomForest)
library(partykit)
library(sparklyr)
library(rpart)
library(tree)
library(RobustGaSP)
library(DiceDesign)
library(foreach)
library(stringr)
library(DiceDesign)
library(stats)
library(tgp)
#library(brms)

work_directory <- "/dbfs/databricks/rstudio/srijit/"
load(paste0(work_directory,"gaussian/ATT49681.RData"))
load(paste0(work_directory,"gaussian/DPP4_both.rda"))
load(paste0(work_directory,"gaussian/modified_fns.RData"))

# COMMAND ----------

run.GPbrm_mod.fn <- function(InputData, nobs=40, randseed=18, nbasfns=10,combfeat = TRUE,
                               thin=2, chains=4, warmup=2000, iter=5000)
  {
    # This function sets up a call to brm to perform GP calculations
    # INPUT
    #   InputData = a list (usually produced by GPinut1.fn) containing
    #                train = a data frame whose last column contains response values (y) and first 
    #                        columns contain values of one or more features. Input data for GP calcs.
    #                 test = a data frame like train used to test validity of predictions based on train
    #                  nfc = number of continuous-valued features in train and test (these will
    #                        occupy the first nfc columns)  
    #                  nfb = number of categorical-valued features (possibly 0) in train and test
    #                        (these will occupy nfb columns after the first nfc columns)
    #        nobs = number of rows of the training dataset to be input to the calculations; 
    #               if 0, all rows will be included; if > 0 then nsamp rows will be sampled 
    #               from the training dataset.  
    #    combfeat = TRUE if the continuous and categorical features are to be combined as
    #               nominal continuous features, FALSE if they are to be considered
    #               separately.
    #    randseed = seed for random no. generator (for reproducibility)
    #               FALSE if the continuous features are drawn from a unit cube.
    #     nbasfns = the brm implementation of GP can approximate the response surface by a 
    #               weighted sum of basis functions, which speeds up the calculations.  
    #               Setting nbasfns = NA causes calculation of the actual response surface.
    #               combinations and the corresponding responses.
    #        thin = thinning parameter for STAN
    #      chains = number of parallel MCMC runs to be performed
    #      warmup = number of warmup iterations before obtaining actual posterior simulations
    #               (STAN parameter)
    #        iter = number of post-warmup iterations for each chain
    #    randseed = starting value for random number generator, to facilitate reproducibility
    #
    # OUTPUT
    #        call = calling sequence (arguments, etc.)
    #        date = time and date when the function was executed
    #  Elapsed.Time = time required for calculations
    #   pred_both = list consisting of predicted values for input array (pred_train) and 
    #               validation array (pred_test), and rmse values for both arrays
    #               (rmse_train, rmse_test) 
    #    pred_new = predicted values corresponding to the n_newx features selected using a
    #               Latin Hypercube design
    #       plots = plots of observed response values versus predicted values for the 
    #               training and testing data
    #      brmout = output from brm (long list)
    #                                            A. L. Gould  September 2020, January 2021
    
    step.fn <- function(xy, nbasfns=10, thin=2, chains=4, warmup=2000, iter=5000)
    {
      nx <- ncol(xy) - 1
      cn <- colnames(xy)[1:nx]
      str1 <- paste0("brm(y~gp(",paste(cn,collapse=","),",k=",nbasfns,",c=5/4)")
      str2 <- paste0(",xy,chains=",chains,",warmup=",warmup,",iter=",iter+warmup,",thin=",thin,")")
      
      print(paste0("BRMS function call: ", str1, str2))
      
      brmout <- eval(parse(text=paste0(str1,str2)))             # Execute brm
      
      print("BRMS call complete")
      return(list(brmout=brmout, sumry_diags=Sumry_Diags_brm.fn(brmout)))
    }    
    ############################  run.GPbrm.fn starts here   ############################
  
    t1 <- proc.time()[3]
    plots <- vector("list",length=2)
    
    ###################  Get key input  ##########################################
    
    xy.train <- InputData$train
    xy.test <- InputData$test 
    new_x <- InputData$new_x
    nfc <- InputData$nfc
    nfb <- InputData$nfb
    ifelse(nobs > 0, ntrpts <- min(nrow(xy.train), nobs), ntrpts <- nrow(xy.train))                          
    if (ntrpts < nrow(xy.train)) 
      xy.train <- sample_xy.fn(xy.train,nsamp=ntrpts, randseed=randseed)$samp.xy
    nx <- nfc + nfb
    ifelse(combfeat, nf <- nfc + nfb, nf <- nfc)
    
    ####################  Execute the selected algorithm  ########################
    
    z <- step.fn(xy.train, nbasfns=nbasfns, thin=thin, chains=chains, warmup=warmup, iter=iter)
    
    ####################  Training set predictions  ##############################
    
    z1 <- posterior_summary(posterior_epred(z$brmout))[,c(1:4)] 
    rmse_train <- sqrt(mean((xy.train$y-z1[,1])^2))
    pred_train <- cbind(xy.train, z1)
    names(pred_train)[nf+(2:3)] <- c("Yhat","sdYhat")
    plots[[1]] <- plot_obsVpred.fn(pred_train$Yhat, pred_train$y, titl="Train [GPbrm]")
    
    ####################  Testing set predictions  ###############################
    
    if(!is.null(xy.test)) 
    {
      z2 <- posterior_summary(posterior_epred(z$brmout, xy.test[,1:nx]))[,c(1:4)]
      rmse_test <- sqrt(mean((xy.test$y-z2[,1])^2))
      pred_test <- cbind(xy.test, z2)
      names(pred_test) <- names(pred_train)
      plots[[2]] <- plot_obsVpred.fn(pred_test$Yhat, pred_test$y, titl="Test [GPbrm]")
    }
    else pred_test <- rmse_test <- NULL
    
    ####################  Predictions for new points  #############################
    
    if (!is.null(new_x)) 
    {
      z3 <- posterior_summary(posterior_epred(z$brmout,new_x))[,c(1,3,4)]
      colnames(z3)[1] <- "Yhat"
      pred_new <- cbind(new_x, z3)
    } 
    else pred_new <- NULL
    
    return(list(call=sys.call(), date=date(), Elapsed.Time=Elapsed.time.fn(t1),
                train.pred=pred_train, test.pred=pred_test, new.pred=pred_new, train.rmse=rmse_train,
                test.rmse=rmse_test, brmout=z$brmout))
  }

# COMMAND ----------

pred_sumry_mod.fn <- function(outlist, nbest=3, bestquant=0.95, xi=0.1, eps=0.001)
{
# Summarizes the findings from separate runs of GP programs    
#   
  get.optEI.fn <- function(outlist, nbest=3, bestquant=0.95, xi=0.1, eps=0.001)
  {
    Yhat <- outlist$new.pred$Yhat
    if (bestquant < 0) ymax <- max(outlist$train.pred$y)
    else ymax <- quantile(outlist$train.pred$y,probs=bestquant)
    eval.pred.diff <- Yhat - ymax
    rmse <- sqrt(mean(eval.pred.diff^2))
    sig.pred <- outlist$new.pred$sdYhat
    z.diff <- eval.pred.diff/pmax(eps, sig.pred)
    EI <- (sig.pred > eps)*(eval.pred.diff*pnorm(z.diff) + sig.pred*dnorm(z.diff))
    pred <- as.data.frame(cbind(outlist$new.pred, EI))
    xyopt <- pred[order(pred$EI,decreasing=T),][1:min(nbest,nrow(pred)),]
    return(list(new_x=xyopt, pred=pred, rmse=rmse)) 
  }    
  
  notest <- is.null(outlist[[1]]$test.pred)
  sd2.te <- Yhat.te.comb <- rmse.tr.comb <- 0
  nsub <- length(outlist)
  rmse.tr <- rmse.te <- vector('numeric',nsub)
  EIopt <- vector("list",nsub)
  plts <- vector("list",nsub+1)
  comb.newx <- sd.tr <- Yhat.tr <- train.pred <- NULL
  #nv <- ncol(outlist[[1]]$test.pred) - 2
  nv <- ncol(outlist[[1]]$train.pred) - 2
  #print(paste0("NV value: ", nv))
  #print(outlist[[1]])
  
  for (i in 1:nsub)
  {
    rmse.tr[i] <- outlist[[i]]$train.rmse
    rmse.tr.comb <- rmse.tr.comb + rmse.tr[i]^2
    Yhat.tr <- c(Yhat.tr, outlist[[i]]$train.pred[,nv-1])
    sd.tr <- c(sd.tr, outlist[[i]]$train.pred[,nv])
    train.pred <- rbind(train.pred, outlist[[i]]$train.pred[,1:nv])
    plts[[i]] <- with(outlist[[i]]$train.pred, plot_obsVpred.fn(Yhat, y ,titl=paste0("Train, subset",i)))
    if (!notest)
    {
      rmse.te[i] <- outlist[[i]]$test.rmse
      sd2.te <- sd2.te + outlist[[i]]$test.pred$sdYhat^2
      Yhat.te.comb <- Yhat.te.comb + outlist[[i]]$test.pred$Yhat 
    }
    if (!is.null(outlist[[1]]$new_x)) EIopt[[i]] <- NULL
    else 
    {
      EIopt[[i]] <- get.optEI.fn(outlist[[i]], nbest=nbest, bestquant=bestquant, xi=xi, eps=eps)
      comb.newx <- rbind(comb.newx, EIopt[[i]]$new_x)
    }
  }
  Q2.5 <- Yhat.tr - qnorm(0.975)*sd.tr
  Q97.5 <- Yhat.tr + qnorm(0.975)*sd.tr
  train.pred <- cbind(train.pred,Q2.5,Q97.5)
  rmse.tr.comb <- sqrt(rmse.tr.comb/nsub)
  plts[[nsub+1]] <- plot_obsVpred.fn(train.pred$Yhat, train.pred$y, titl="Train, average prediction")
  if (notest){
    rmse.te <- sd2.te <- Yhat.te.comb <- rmse.te.comb <- Q2.5.te.comb <- Q97.5.te.comb <- test.pred.comb <- NULL
  }
  else      
  {
    sd2.te <- sqrt(sd2.te/nsub)
    Yhat.te.comb <- Yhat.te.comb/nsub
    rmse.te.comb <- sqrt(mean(rmse.te^2))
    Q2.5.te.comb <- Yhat.te.comb - qnorm(0.975)*sd2.te
    Q97.5.te.comb <- Yhat.te.comb + qnorm(0.975)*sd2.te
    test.pred.comb <- cbind(outlist[[1]]$test.pred[,1:nv], Yhat.te.comb, sd2.te, Q2.5.te.comb, Q97.5.te.comb)
    names(test.pred.comb) <- names(outlist[[1]]$test.pred)
    plts[[nsub+2]] <- plot_obsVpred.fn(test.pred.comb$Yhat, test.pred.comb$y, titl="Test, average prediction")
  }
  if (!is.null(outlist[[1]]$new_x)) 
    comb.newx <- comb.newx[order(comb.newx$EI,decreasing=T),][1:min(nbest,nrow(comb.newx)),]

  return(list(call=sys.call(), date=date(), rmse.tr=rmse.tr, rmse.tr.comb=rmse.tr.comb, rmse.te=rmse.te, 
              rmse.te.comb=rmse.te.comb, train.pred=train.pred, test.pred.comb=test.pred.comb, 
              newx.comb=comb.newx, EIopt=EIopt#, plts=plts
             ))
}


# COMMAND ----------

read_data.fn <- function(fname){
  #print(paste0("Reading:", fname))
  d <- readRDS(fname)
  #print(paste0("Deleting:", fname))
  file.remove(fname)
  return(d)
}

# COMMAND ----------

compute.fn <- function(i, partInput, params){
    
  calltxt2 <- paste0("(InputData=partInput[[i]], nobs=", nobs, ", randseed=", randseed)

  fnlib <- params[[1]]
  .GlobalEnv$run.GPtgp.fn  <- params[[2]]
  .GlobalEnv$run.GPhetGP.fn <- params[[3]]
  .GlobalEnv$run.GPrgasp.fn <- params[[4]]
  .GlobalEnv$run.GPfit.fn <- params[[5]]
  .GlobalEnv$run.GPLVGP.fn <- params[[6]]
  .GlobalEnv$run.GPlasvdGP.fn <- params[[7]]
  .GlobalEnv$run.GPbrm.fn <- params[[8]]
  .GlobalEnv$run.GPlaGP.fn <- params[[9]]
  .GlobalEnv$plot_obsVpred.fn <- params[[10]]
  .GlobalEnv$Elapsed.time.fn <- params[[11]]
  .GlobalEnv$vec2mat.fn <- params[[12]]
  .GlobalEnv$sample_xy.fn <- params[[13]]
  file_path <- params[[14]]
  outname <- params[[15]]

  switch(whichfn,
         st2 <- paste0("run.GPtgp.fn", calltxt2, ", fnopt=", fnopt,", verb=", verb, ", catonly4tree=", catonly4tree,")") ,# tgp
         st2 <- paste0("run.GPhetGP.fn", calltxt2, ")"),                                                                  # hetGP
         st2 <- paste0("run.GPrgasp.fn", calltxt2, ", nugget.est=", nugget.est, ", lower_bound=",lower_bound, ")"),       # RobustGaSP
         st2 <- paste0("run.GPfit.fn", calltxt2, ")"),                                                                    # GPfit
         st2 <- paste0("run.GPLVGP.fn", calltxt2, ")"),                                                                   # LVGP
         st2 <- paste0("run.GPlasvdGP.fn", calltxt2, ", nthread=", nthread,")" ),                                         # lasvdGP 
         st2 <- paste0("run.GPbrm_mod.fn", calltxt2, ", nbasfns=", nbasfns, ", thin=", thin, ", chains=", chains, ", warmup=", warmup, ", iter=",iter, ")"), # brms
         st2 <- paste0("run.GPlaGP.fn", calltxt2, ")")                                                                    # laGP
  )

  print(paste0("Method signature for ", i, " : ",st2))
  t1 <- proc.time()[3]

  result <- tryCatch({
      r <- eval(parse(text=st2))
      print(paste0("Completed method execution for ", i))
      r
    },
    error = function(e){
      print(paste0("Error while executing ", i))
      print(e)
      return(list(date=date(), call=sys.call(), Elapsed.Time=Elapsed.time.fn(t1), train.pred=NULL,
            test.pred=NULL, new.pred=NULL, train.rmse=NULL, test.rmse=NULL, plots=NULL))
    }
  )

  sl <- object.size(result)      
  print("Result size in executor:")  
  print(sl, units = "auto")

  file_name <- paste0(file_path, "res_",outname,"_",i,"_",as.numeric(Sys.time()))
  
  print(paste0("Result file:",file_name))
  
  #saveRDS(result, file=file_name)
  #return(file_name)  
  return(result)  
}

# COMMAND ----------

runGP_mod.fn <- function(outname, InputData="DPP4_both", TestData="NULL", indx.cont="confeat", indx.cat="catfeat", 
         nobs=100, n_newx=20, normXform=TRUE, UseMedResp=FALSE, combfeat=TRUE, fnlib="tgp", fnopt=3,
         catonly4tree=TRUE, cnflev=0.95, randseed=241, parallel=FALSE, progress=FALSE, eps=0.01,
         corropt="ma", pwr=2, nu=5/2, nthread=1, nbasfns=10, thin=5, nugget.est=TRUE,
         lower_bound=FALSE, chains=4, warmup=2000, verb=0, partsize=1000, iter=5000)
{
  # This function is a driver for parallelized GP computations.
  #
  # INPUT
  #      outname = text string giving the name of the list that will contain the computational results
  #    InputData = name of the data file containing the training dataset 
  #     TestData = name of the data file (if any) containing the testing dataset 
  #    indx.cont = character string identifying the name of a variable that contains the
  #                indices of the categorical features in the datasets in InputData corresponding 
  #                to continuous features to include in the calculations (or NULL).
  #     indx.cat = character string identifying the name of a variable that contains the indices
  #                categorical features in the datasets in InputData corresponding to the
  #         nobs = number of observations to be taken for GP calculations from each partition 
  #                in the training dataset contained in InputData 
  #       n_newx = number of new feature vectors to be selected using LHS sampling for predicting response
  #    normXform = TRUE if the continuous feature values are standard normal deviates
  #                FALSE if the values are selections from a unit hypercube
  #   UseMedResp = TRUE if redundant feature combinations are to be combined with a common
  #     combfeat = TRUE if categorical parameters are to be mapped to latents and
  #                treated as continuous
  #        fnlib = name of package containing the GP function that is to be run
  #                response equal to the median of the corresponding responses, FALSE if not
  #        fnopt = option for choosing which function in tgp is to be run.  Ignored if not tgp.
  # catonly4tree = TRUE if the trees constructed using the bgp and bgpllm functions are based
  #                only on categorical-valued features, FALSE if all features can be used
  #     randseed = starting value for random number generator (allows for reproducibility)
  #                categorical features to include in the calculations (or NULL)
  #     progress = parameter for LVGP
  #     parallel = parameter for LVGP
  #      nthread = parameter for lasvdGP
  #      corropt = parameter for GPfit
  #         verb = 0 means no verbose output for tgp
  #     partsize = size of each partition when partitioning the InputData training set
  #  nbasfns, thin, warmup, iter = parameters for GPbrm
  #
  # OUTPUT written directly to workspace
  #                                           A.L.Gould September, December 2021
  t1 <- proc.time()[3]
  
  ############  Choose a GP function library ####################
  
  fnlibs <- c("tgp","hetGP","RobustGaSP","GPfit","LVGP","DynamicGP","brms","laGP")
  whichfn <- grep(fnlib,fnlibs)
  
  ###########  Prepare the input dataset  #######################					 
  
  InDta <- eval(parse(text=InputData))                      
  TstDta <- eval(parse(text=TestData))
  indx_cont <- eval(parse(text=indx.cont))
  indx_cat <- eval(parse(text=indx.cat))  
  In_Dta <- GPinput.fn(InDta,TstDta, indx.cont=indx_cont, indx.cat=indx_cat, n_newx=n_newx, 
                       normXform=normXform, combfeat=combfeat, UseMedResp=UseMedResp, eps=eps)
  
  ############  Partition the input dataset  ####################  
  
  # Randomly partition the rows of the 'train' member of the InputData list
  # dataset into about equal subsets of approximately partsize  rows.  The 
  # result is a list of lists similar to InputData with smaller 'train' subsets.
  # The other members of the InputData list are copied to each of these
  # component lists.  This allows parallelization of the calculations with
  # an appreciable saving of computation time.

  mm <- simple_partition.fn(nrow(In_Dta$train), partsize, rseed=randseed)
  
  print(paste0("Partition Size:", partsize))
  print(paste0("Number of Rows:", nrow(In_Dta$train)))
  print(paste0("Number of partitions:", length(mm)))
  
  partInput <- vector("list", length=length(mm))
  for (i in 1:length(mm))
  {
    partInput[[i]] <- In_Dta
    partInput[[i]]$train <- In_Dta$train[mm[[i]],]
  }
  
  eval(paste0("library(",fnlib,")"))
  
  sc = sparklyr::spark_connect(method = "databricks")
  registerDoSpark(sc)
  
  paramsToPass = list(fnlib, run.GPtgp.fn,run.GPhetGP.fn,run.GPrgasp.fn,
                     run.GPfit.fn,run.GPLVGP.fn,run.GPlasvdGP.fn,run.GPbrm_mod.fn,
                     run.GPlaGP.fn,plot_obsVpred.fn, Elapsed.time.fn,vec2mat.fn, sample_xy.fn, work_directory, outname)
    
  
  res <- foreach(i=1:length(partInput)) %dopar% compute.fn(i, partInput, paramsToPass)
  #res_files <- foreach(i=1:length(partInput)) %dopar% compute.fn(i, partInput, paramsToPass)
  print("Results obtained")
  
  #res <- foreach(i=1:length(res_files)) %do% read_data.fn(res_files[[i]])
  #print("Calculation complete")
  #print(res)
  #return(res)
    
  output <- tryCatch(
    pred_sumry_mod.fn(res),
    error = function(e){
      print("Error calling pred_sumry_mod.fn:")
      print(e)
      return(list(date=date(), call=sys.call(), Elapsed.Time=Elapsed.time.fn(t1)))
    }
  )
  
  return(output)
}



# COMMAND ----------

catfeat <- c(1, 9, 10, 16, 17, 19, 20, 22, 23, 33, 37)


# COMMAND ----------

#Testing with a small batch size to make sure everything is working fine
start_time <- Sys.time()
print(start_time)
res_gasp <- runGP_mod.fn("DPP4a_rgasp_1000_nocat.all",InputData="DPP4_both", TestData="NULL",
         fnlib="RobustGaSP", nobs=0, indx.cont="NULL", indx.cat="zero", n_newx=0,
         combfeat=FALSE, normXform=TRUE, nugget.est=T, lower_bound=F, randseed=rs,partsize=50)

print(paste0("Run time (minutes):", (as.numeric(Sys.time() - start_time , units = "mins"))))

# COMMAND ----------

#Run with batch size 500
start_time <- Sys.time()
print(start_time)
res_gasp <- runGP_mod.fn("DPP4a_rgasp_1000_nocat.all",InputData="DPP4_both", TestData="NULL",
         fnlib="RobustGaSP", nobs=0, indx.cont="NULL", indx.cat="zero", n_newx=0,
         combfeat=FALSE, normXform=TRUE, nugget.est=T, lower_bound=F, randseed=rs,partsize=500)

print(paste0("Run time (minutes):", (as.numeric(Sys.time() - start_time , units = "mins"))))

# COMMAND ----------

# MAGIC %md
# MAGIC The above command will keep running forever, even after the Spark job is completed.. but if you look at driver logs you can see
# MAGIC ```
# MAGIC OpenJDK 64-Bit Server VM warning: ignoring option MaxPermSize=512m; support was removed in 8.0
# MAGIC Mon Oct 10 13:11:34 2022 Python shell started with PID  3976  and guid  d4341cd6d41a4ce38abee7addad9f383
# MAGIC Mon Oct 10 13:11:34 2022 Initialized gateway on port 45159
# MAGIC Mon Oct 10 13:11:37 2022 Python shell executor start
# MAGIC (22/10/10 13:11:43 ERROR sparklyr: Gateway (65529) failed with exception ,java.io.EOFException)
# MAGIC java.lang.OutOfMemoryError
# MAGIC 	at java.io.ByteArrayOutputStream.hugeCapacity(ByteArrayOutputStream.java:123)
# MAGIC 	at java.io.ByteArrayOutputStream.grow(ByteArrayOutputStream.java:117)
# MAGIC 	at java.io.ByteArrayOutputStream.ensureCapacity(ByteArrayOutputStream.java:93)
# MAGIC 	at java.io.ByteArrayOutputStream.write(ByteArrayOutputStream.java:153)
# MAGIC 	at java.io.DataOutputStream.write(DataOutputStream.java:107)
# MAGIC 	at java.io.FilterOutputStream.write(FilterOutputStream.java:97)
# MAGIC 	at sparklyr.Serializer$.writeBytes(serializer.scala:199)
# MAGIC 	at sparklyr.Serializer.writeObject(serializer.scala:449)
# MAGIC 
# MAGIC ```

# COMMAND ----------


