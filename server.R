###########################
# MSP_DSC
# version 1.0.0 Fevrier 2018. 
# version 2.0.0 October 2026. DSC alpha cases 1-4
# version 2.1.0 October 2026. DSC alpha case 5 (prior centred on [a_min, a_max])
# version 2.2.0 October 2026. Paired bootstrap comparison of the five DSC cases
###########################


source("MSP_func.R")



function(input, output) {
  # loading measurements in the selected file
  meas <- reactive({
    inFile <- input$ifile
    if (is.null(inFile)) return(NULL)
    read.csv(inFile$datapath, header = TRUE, sep = ",",comment.char = "#",blank.lines.skip=TRUE)
  })
  
  
  # N'ouvre les menus que si un fichier de donnees est ouvert et bien lu.
  observeEvent(is.null(meas()), {
    hideTab(inputId = "tabs", target = "Data explorer")
    hideTab(inputId = "tabs", target = "Analysis and Plots")
  })
  
  observeEvent(input$ifile, {
    showTab(inputId = "tabs", target = "Data explorer")
    showTab(inputId = "tabs", target = "Analysis and Plots")
  })
  
  
  
  # choose the protocol
  proto <- reactive({
    input$Protocol
  })
  
  # lower and upper bounds for NRM fraction
  #sel_fraction <- reactive({
  #  input$sliderF
  #})
  
  # lower and upper bounds for NRM fraction
  #critangle <- reactive({
  #  input$Angle
  #})
  
  # alpha parameter of case 1: fixed value (no slider)
  alpha <- reactive({
    alpha_case1
  })
  
  # Number of bootstrapped values
  R.rep <- reactive({
    input$NBoot
  })

  # Cases (1 to 5) of MSP-DSC already computed with the current data and settings.
  # The comparison of the five cases (tab 'Compare cases') is run only when all
  # of them are available.
  done <- reactiveValues(c1 = FALSE, c2 = FALSE, c3 = FALSE, c4 = FALSE, c5 = FALSE)
  observeEvent(list(input$ifile, input$sliderF, input$Angle, input$NBoot), {
    done$c1 <- FALSE; done$c2 <- FALSE; done$c3 <- FALSE; done$c4 <- FALSE; done$c5 <- FALSE
  }, ignoreInit = TRUE, priority = 100)
  observeEvent(input$kSigma, { done$c5 <- FALSE }, ignoreInit = TRUE, priority = 100)
  
  #browser()
  ###############
  # QDB
  ###############
  
  fit.DB <- reactive({
    B <- as.numeric(meas()$Field)
    #browser()    
    ndata <- nrow(meas())
    QDB <- (meas()$m1 - meas()$m0) / meas()$m0 # Eq. 15 in Fabian & Leonhardt, 2010
    
    # Ordinary least squares regresion with possible outliers
    B.QDB <- data.frame(Name = meas()$name)
    B.QDB$x <- B  
    B.QDB$y <- QDB
    model_DB <- lm(y ~ x, data = B.QDB) 

    # Testing for the presence of influential points (eg Outliers)
    B.QDB$Cook <- cooks.distance(model_DB)
    threshold_Cook <- 3 * mean(B.QDB$Cook)
    DB_influential.data <- which(B.QDB$Cook > threshold_Cook)
    quality_DB <- rep("Good",length(B))
    quality_DB[DB_influential.data] <- "Outlier"
    B.QDB$quality <- quality_DB
    nDB_influential <- length(DB_influential.data)
   
    # OLS without Influential data
    model_DB <- lm(y ~ x, data = subset(B.QDB,quality == "Good")) 
    
    # Bootstrapping for uncertainty
    data_sel <- subset(B.QDB,quality == "Good", select=c(x,y))
    reg.model <- boot(data_sel,reg.boot,R=R.rep())
    tt <- boot.ci(reg.model,conf = conf_level, type = "norm")
    
# build the data frame for statistical tests and plots
    b.val <- reg.model$t
    b.prot <- rep("DB",R.rep())
    dat <- data.frame(b.prot,b.val)
    
    res <- data.frame("pi_val" = reg.model$t0)
    res$l68 <- tt$normal[3]
    res$u68 <- tt$normal[5]
    res$l95 <- tt$normal[4]
    res$u95 <- tt$normal[6]
    res$ndataT <- ndata
    res$ndataS <- ndata - nDB_influential
    res$r2  <- summary(model_DB)[[8]]
    
    # Build dataframe for the best fit lines
    best_fit <- data.frame("Field" = c(0,max(B) + 3))
    best_fit$DB <- model_DB[[1]][2] * best_fit$Field + model_DB[[1]][1]
    
    new <- data.frame(x = seq(0,max(B) + 3,length.out = 200))
    yv <- predict(model_DB,newdata=new,interval="confidence",level=0.95)
    error_DB <- data.frame("Field" = c(new$x,rev(new$x)))
    error_DB$y_err <- c(yv[,2],rev(yv[,3]))
     
    # create a list of data for use in rendering
    info <- list(QDB=B.QDB,best_fit=best_fit,err_DB=error_DB,res=res,dat=dat)
    return(info)
    
  })
  
  ###############
  # QFC
  ###############
  # Geometry shared by FC and DSC (no bootstrap): fraction, angle, Qfc, corrected m2
  fc.base <- reactive({
    B <- as.numeric(meas()$Field)
    
    # transform Inc, Dec in deg to rad
    rD0 <- deg2rad(meas()$D0); rI0 <- deg2rad(meas()$I0)
    rD1 <- deg2rad(meas()$D1); rI1 <- deg2rad(meas()$I1)
    rD2 <- deg2rad(meas()$D2); rI2 <- deg2rad(meas()$I2)

    #browser()

    ## array of vectors in cartesian coordinates
    ndata <- nrow(meas())
    X0 <- sph2cart(matrix(c(rD0,rI0,meas()$m0),nrow = ndata))
    X1 <- sph2cart(matrix(c(rD1,rI1,meas()$m1),nrow = ndata))
    X2 <- sph2cart(matrix(c(rD2,rI2,meas()$m2),nrow = ndata))
   
    # calculate the fraction
    demag <- (X1 + X2) / 2
    norm_demag <- apply(demag,1,Norm) # func Norm from pracma
    fraction <- (meas()$m0 - norm_demag) / meas()$m0
    
    #     Angle between m0 and demag
    # and Angle between m0 and m2
    angle_demag <- vector(mode = "numeric", length = ndata)
    angle_m0m2 <- vector(mode = "numeric", length = ndata)
    
    for (i in 1:ndata) {
      angle_demag[i] <- func_angulardiff(X0[i,],demag[i,])
      angle_m0m2[i] <- func_angulardiff(X0[i,],X2[i,])
    }
    
    # Reverse m2 if the angle between m2 and m0 is greater than pi/2
    mm2 <- meas()$m2
    ind_reverse <- which(angle_m0m2 > pi / 2)
    mm2[ind_reverse] <- -mm2[ind_reverse]
    
    # Eq. 16 in Fabian and Leonhardt, 2010.
    QFC <- 2 * (meas()$m1 - meas()$m0) / (2 * meas()$m0 - meas()$m1 - mm2)
    B.QFC <- data.frame(Name = meas()$name)
    B.QFC$x <- B 
    B.QFC$y <- QFC
    B.QFC$f <- fraction
    B.QFC$angle <- rad2deg(angle_demag)

    list(B = B, B.QFC = B.QFC, mm2 = mm2)
  })

  fit.FC <- reactive({
    base  <- fc.base()
    B     <- base$B
    B.QFC <- base$B.QFC
    mm2   <- base$mm2
    QFC   <- B.QFC$y

    # datapoint selection (select "Good" data)
    ind_sel <- c(which(B.QFC$f  > input$sliderF[1]
                     & B.QFC$f  < input$sliderF[2]
                     & B.QFC$angle < input$Angle))
    
    
    
    quality_FC <- rep("Bad",length(B))
    if (!is.null(ind_sel)) {
    quality_FC[ind_sel] <- "Good"
    }
    B.QFC$quality <- quality_FC
     
    # OLS anchored to the point (0,-1)
    model_FC <- lm((y + 1) ~ x - 1, data = subset(B.QFC,quality=="Good")) 
    
    ## Testing for the presence of influential points
    B.QFC$Cook[B.QFC$quality == "Good"] <- cooks.distance(model_FC)
    threshold_Cook <- 3 * mean(subset(B.QFC, quality == "Good" )$Cook)
    
    FC_influential.data <- which(B.QFC$Cook > threshold_Cook)
    B.QFC$quality[FC_influential.data] <- "Outlier"
    nFC_influential <- length(FC_influential.data)
    
    ## OLS anchored to the point (0,-1) without Influential data
    data_sel <- subset(B.QFC,quality == "Good", select=c(x,y))
    data_sel$y <- data_sel$y + 1
    data_sel$w <- 1
    model_FC <- lm(y ~ x - 1, data = data_sel)
      
    # bootstrapping for uncertainty
    reg.model <- boot(data_sel,reg.boot_anch,R=R.rep())
    tt <- boot.ci(reg.model, conf = conf_level, type = "norm")

    # build the data frame for statistical tests and plots
    b.val <- reg.model$t
    b.prot <- rep("FCi",R.rep())
    dat <- data.frame(b.prot,b.val)
    
    res <- data.frame("pi_val" = reg.model$t0)
    res$l68 <- tt$normal[3]
    res$u68 <- tt$normal[5]
    res$l95 <- tt$normal[4]
    res$u95 <- tt$normal[6]
    res$ndataT <- length(na.omit(QFC))
    res$ndataS <- length(data_sel$x)
    res$r2  <- summary(model_FC)[[8]]
    
    # Build dataframe for the best fit lines
    best_fit <- data.frame("Field" = c(0,max(B) + 3))
    best_fit$FC <- coef(model_FC) * best_fit$Field - 1
    
    # Build dataframe for the uncertainties
    new <- data.frame(x = c(0,max(B) + 3))
    yv <- predict(model_FC,newdata=new,interval="confidence",level=0.95)
    error_FC <- data.frame("Field" = c(new$x,rev(new$x)))
    error_FC$y_err <- c(yv[,2]-1,rev(yv[,3]-1))
    
    # create a list of data for use in rendering
    info <- list(QFC=B.QFC,best_fit=best_fit,err_FC=error_FC,res=res,dat=dat,mm2=mm2)
    return(info)
    
  })
  
  # ----------------------------------------------------------
  # MSP-DSC protocol, 4 alpha strategies (input$caseDSC)
  #   1 : constant alpha (slider)
  #   2 : single optimized alpha [0, 0.8]
  #   3 : per-specimen alpha, bounds [0, 0.8]
  #   4 : per-specimen alpha, bounds [a_min, a_max] from the file
  # ----------------------------------------------------------

  # Specimen selection, shared by the four cases
  dsc.sel <- reactive({
    f     <- fc.base()$B.QFC$f
    angle <- fc.base()$B.QFC$angle
    m4m3  <- (meas()$m4 - meas()$m3) / meas()$m4
    m1m3  <- (meas()$m1 - meas()$m3) / meas()$m1
    good  <- f > input$sliderF[1] & f < input$sliderF[2] &
             angle < input$Angle & m4m3 > m4m3threshold & m1m3 > m1m3threshold
    good[is.na(good)] <- FALSE
    list(f = f, angle = angle, m4m3 = m4m3, m1m3 = m1m3, good = good)
  })

  # Cases 4 and 5 : a_min / a_max = columns "a_min"/"a_max" if present,
  # otherwise the last (a_min) and second-to-last (a_max) columns of the file
  dsc.bounds <- reactive({
    m  <- meas()
    nc <- ncol(m)
    if (all(c("a_min", "a_max") %in% names(m))) {
      lo <- m[["a_min"]]; hi <- m[["a_max"]]
    } else {
      validate(need(nc >= 2, "Cases 4 and 5 need a_min and a_max in the last two columns of the file."))
      lo <- m[[nc]]; hi <- m[[nc - 1]]   # header order: ..., a_max, a_min
    }
    list(a_min = suppressWarnings(as.numeric(lo)),
         a_max = suppressWarnings(as.numeric(hi)))
  })

  # ---- Case 1 : constant alpha (original code) ----
  fit.DSC1 <- reactive({
    B   <- as.numeric(meas()$Field)
    mm2 <- fc.base()$mm2
    sel <- dsc.sel()

    # Eq. XX in Fabian and Leonhardt, 2010.
    QDSC <- 2*((1+alpha())*meas()$m1-meas()$m0-alpha()*meas()$m3)/(2*meas()$m0-meas()$m1-mm2)

    B.QDSC <- data.frame(Name = meas()$name)
    B.QDSC$x <- B
    B.QDSC$y <- QDSC

    # alteration estimation
    B.QDSC$dQ <- calc_dQ(meas()$m0, meas()$m1, meas()$m3, meas()$m4, mm2, alpha())

    # datapoint selection
    B.QDSC$f     <- sel$f
    B.QDSC$angle <- sel$angle
    B.QDSC$m4m3  <- sel$m4m3
    B.QDSC$m1m3  <- sel$m1m3

    quality_DSC <- rep("Bad", length(B))
    quality_DSC[sel$good] <- "Good"
    B.QDSC$quality <- quality_DSC
    B.QDSC$Cook <- NA_real_

    # WLS on selected points including potential outliers
    model_DSC <- lm((y + 1) ~ x - 1, data = subset(B.QDSC,quality=="Good"), weights=1/dQ)

    # Testing for the presence of influential points
    B.QDSC$Cook[B.QDSC$quality == "Good"] <- cooks.distance(model_DSC)
    threshold_Cook <- 3 * mean(subset(B.QDSC, quality == "Good" )$Cook)
    DSC_influential.data <- which(B.QDSC$Cook > threshold_Cook)
    B.QDSC$quality[DSC_influential.data] <- "Outlier"
    nDSC_influential <- length(DSC_influential.data)

    # WLS anchored to the point (0,-1) without Influential data
    model_DSC <- lm((y + 1) ~ x - 1, data = subset(B.QDSC,quality=="Good"), weights=1/dQ)

    # bootstrapping for uncertainty
    data_sel <- subset(B.QDSC,quality == "Good", select=c(x,y,dQ))
    data_sel$y <- data_sel$y + 1
    data_sel$dQ <- 1/data_sel$dQ

    reg.model <- boot(data_sel,reg.boot_anch,R=R.rep())
    tt <- boot.ci(reg.model, conf = conf_level, type = "norm")

# build the data frame for statistical tests and plots
    b.val <- reg.model$t
    b.prot <- rep("DS",R.rep())
    dat <- data.frame(b.prot,b.val)

    res <- data.frame("pi_val" = reg.model$t0)
    res$l68 <- tt$normal[3]
    res$u68 <- tt$normal[5]
    res$l95 <- tt$normal[4]
    res$u95 <- tt$normal[6]
    res$ndataT <- length(na.omit(QDSC))
    res$ndataS <- length(data_sel$x)
    res$r2  <- summary(model_DSC)$r.squared

    # alpha used (last column of the table)
    B.QDSC$alpha <- alpha()

    # Build dataframe for the best fit lines
    best_fit <- data.frame("Field" = c(0,max(B) + 3))
    best_fit$DSC <- coef(model_DSC) * best_fit$Field - 1

    # Build dataframe for the uncertainties
    new <- data.frame(x = c(0,max(B) + 3))
    yv <- predict(model_DSC,newdata=new,interval="confidence",level=0.95)
    error_DSC <- data.frame("Field" = c(new$x,rev(new$x)))
    error_DSC$y_err <- c(yv[,2]-1,rev(yv[,3]-1))

    # create a list of data for use in rendering
    ans <- list(QDSC = B.QDSC, best_fit = best_fit, err_DSC = error_DSC, res = res, dat = dat,
         case_label = sprintf("case 1 (alpha = %.2f)", alpha()),
         alpha_text = sprintf("%.2f (constant)", alpha()),
         plot_note = sprintf(", alpha = %.2f", alpha()))
    isolate(done$c1 <- TRUE)
    ans
  })

  # ---- Cases 2, 3, 4, 5 : optimized alpha ----
  dsc.opt.fun <- function(cs) {
    m   <- meas()
    n   <- nrow(m)
    B   <- as.numeric(m$Field)
    sel <- dsc.sel()
    mm2 <- fc.base()$mm2
    Qfc <- fc.base()$B.QFC$y

    c_i <- 2 * (m$m1 - m$m3) / (2 * m$m0 - m$m1 - mm2)
    dQ  <- calc_dQ(m$m0, m$m1, m$m3, m$m4, mm2, alpha_dQ)

    # alpha bounds, ridge centre and ridge strength of each specimen
    f_min <- f_max <- NULL      # intervals read in the file (cases 4 and 5)
    a_ctr_all <- lam_all <- NULL
    k <- NA_real_
    if (cs %in% c(4, 5)) {
      bnd   <- dsc.bounds()
      f_min <- bnd$a_min
      f_max <- bnd$a_max
    }
    if (cs == 4) {
      a_min  <- f_min
      a_max  <- f_max
      centre <- "mid"
    } else if (cs == 5) {
      # alpha stays in [0, 0.8]; the interval of the file defines a Gaussian prior
      k      <- if (is.null(input$kSigma)) sigma_k_default else input$kSigma
      a_min  <- rep(alpha_lo, n)
      a_max  <- rep(alpha_hi, n)
      pri    <- case5_prior(f_min, f_max, k)
      a_ctr_all <- pri$ctr
      lam_all   <- pri$lam
      centre <- "col"
    } else {
      a_min  <- rep(alpha_lo, n)
      a_max  <- rep(alpha_hi, n)
      centre <- alpha_pop
    }

    use <- sel$good & is.finite(Qfc) & is.finite(c_i) & is.finite(dQ) & dQ > 0
    validate(need(sum(use) >= 4,
                  "Fewer than 4 specimens pass the selection criteria: alpha cannot be optimized."))
    if (cs %in% c(4, 5)) {
      validate(need(
        all(is.finite(f_min[use]), is.finite(f_max[use])) &&
          all(f_min[use] >= 0 & f_min[use] <= f_max[use] & f_max[use] <= 1),
        sprintf("Case %d: the last two columns of the file must be a_min and a_max (0 <= a_min <= a_max <= 1).", cs)))
    }

    d <- data.frame(id = seq_len(n), Field = B, m0 = m$m0, m1 = m$m1, m3 = m$m3,
                    m2_corr = mm2, Qfc = Qfc, dQ = dQ,
                    a_min = a_min, a_max = a_max)
    if (cs == 5) {
      d$a_ctr <- a_ctr_all
      d$lam_i <- lam_all
    }
    d <- d[use, , drop = FALSE]

    field_grid <- seq(0, max(B) + 3, length.out = 200)
    nb  <- R.rep()

    # cases 3 and 4: lambda selected by leave-one-out cross-validation
    # (one lambda per case: different shrinkage target and bounds)
    # case 5: lambda_i = 1 / sigma_i^2 is given by the intervals of the file (no LOOCV)
    loocv <- NULL
    lam   <- NA_real_
    if (cs %in% c(3, 4)) {
      loocv <- select_lambda_loocv(d, lambda_grid, centre)
      lam   <- loocv$lambda[which.min(loocv$loocv_error)]
    }

    res <- withProgress(message = "Optimizing alpha", value = 0, {
      prog <- function(x) setProgress(value = x)
      if (cs == 2) {
        run_global_alpha(d, alpha_lo, alpha_hi, nb, field_grid, prog)
      } else {
        run_perspecimen_alpha(d, lam, centre, nb, field_grid, prog)
      }
    })

    # alpha and Qdsc for ALL specimens: those not used in the final fit
    # (Bad, Cook outliers) keep the global alpha (case 2) or their ridge centre
    alpha_all <- if (cs == 2) rep(res$alpha_global, n) else resolve_alpha_center(centre, a_min, a_max, a_ctr_all)
    alpha_all[res$clean_ids] <- res$alpha_clean
    Q_all <- Qfc + alpha_all * c_i

    quality <- rep("Bad", n)
    quality[d$id] <- "Good"
    quality[res$outlier_ids] <- "Outlier"
    Cook <- rep(NA_real_, n)
    Cook[d$id] <- res$cook

    B.QDSC <- data.frame(Name = m$name, x = B, y = Q_all, dQ = dQ,
                         f = sel$f, angle = sel$angle, m4m3 = sel$m4m3, m1m3 = sel$m1m3,
                         quality = quality, Cook = Cook, alpha = alpha_all)

    out <- data.frame("pi_val" = 1 / res$slope)
    out$l68 <- res$boot$ci68[1]
    out$u68 <- res$boot$ci68[2]
    out$l95 <- res$boot$ci95[1]
    out$u95 <- res$boot$ci95[2]
    out$ndataT <- length(na.omit(Q_all))
    out$ndataS <- length(res$clean_ids)
    out$r2  <- res$r2

    best_fit <- data.frame("Field" = c(0, max(B) + 3))
    best_fit$DSC <- res$slope * best_fit$Field - 1

    error_DSC <- data.frame("Field" = c(field_grid, rev(field_grid)))
    error_DSC$y_err <- c(res$boot$band_lower, rev(res$boot$band_upper))

    dat <- data.frame(b.prot = rep("DS", length(res$boot$pi)), b.val = res$boot$pi)

    ac <- res$alpha_clean
    lam_rng <- if (cs == 5) range(lam_all[use]) else c(NA_real_, NA_real_)
    switch(as.character(cs),
      "2" = { lab <- sprintf("case 2 (alpha = %.3f)", res$alpha_global)
              txt <- sprintf("%.3f (single, optimized in [0, 0.8])", res$alpha_global) },
      "3" = { lab <- "case 3 (per-specimen alpha, [0, 0.8])"
              txt <- sprintf("per-specimen in [0, 0.8]: mean %.2f, range %.2f-%.2f (lambda = %s)",
                             mean(ac), min(ac), max(ac), format(lam)) },
      "4" = { lab <- "case 4 (per-specimen alpha, [a_min, a_max])"
              txt <- sprintf("per-specimen in [a_min, a_max]: mean %.2f, range %.2f-%.2f (lambda = %s)",
                             mean(ac), min(ac), max(ac), format(lam)) },
      "5" = { lab <- "case 5 (per-specimen alpha, prior centred on [a_min, a_max])"
              txt <- sprintf("per-specimen in [0, 0.8], prior centred on [a_min, a_max]: mean %.2f, range %.2f-%.2f (sigma = width/%g, lambda_i = %.0f-%.0f)",
                             mean(ac), min(ac), max(ac), k, lam_rng[1], lam_rng[2]) })

    ans <- list(QDSC = B.QDSC, best_fit = best_fit, err_DSC = error_DSC, res = out, dat = dat,
         case_label = lab, alpha_text = txt, k = k,
         plot_note = if (cs == 2) sprintf(", alpha = %.2f", res$alpha_global) else "",
         alpha_global = res$alpha_global,
         lambda = lam, lambda_range = lam_rng, loocv = loocv,
         # intervals of the file for cases 4 and 5 (case 5: informative, alpha is searched in [0, 0.8])
         bounds = if (cs == 5) data.frame(a_min = f_min, a_max = f_max)
                  else data.frame(a_min = a_min, a_max = a_max))
    isolate(done[[paste0("c", cs)]] <- TRUE)
    ans
  }

  fit.DSC.c2 <- reactive({ dsc.opt.fun(2) })
  fit.DSC.c3 <- reactive({ dsc.opt.fun(3) })
  fit.DSC.c4 <- reactive({ dsc.opt.fun(4) })
  fit.DSC.c5 <- reactive({ dsc.opt.fun(5) })

  # ---- Dispatcher used by the tables, plots and summary ----
  fit.DSC <- reactive({
    req(input$caseDSC)
    switch(as.character(input$caseDSC),
           "1" = fit.DSC1(), "2" = fit.DSC.c2(),
           "3" = fit.DSC.c3(), "4" = fit.DSC.c4(), "5" = fit.DSC.c5())
  })
  
  # Print Measurement file 
  output$contents <- DT::renderDataTable({
    DT::datatable(meas(),class="display",style="bootstrap",options = list(lengthMenu = c(10, 20, 50), pageLength = 10))
  })
  
  table.calc <- reactive({
    switch(input$Protocol,
           DB = {
             bq <- fit.DB()$QDB
             colnames(bq) <- c("Sample","B","Qdb","Cook","Quality")
             bq <- bq[
               with(bq, order(B, Sample)),
               ]
             return(bq)
           },
           FC = {
             bq <- fit.FC()$QFC
             colnames(bq) <- c("Sample","B","Qfc","f","angle","Quality","Cook")
             bq <- bq[
               with(bq, order(B, Sample)),
               ]
             bq <- bq[c(1,2,3,4,5,7,6)]
             return(bq)
           },
           DSC = {
             bq <- fit.DSC()$QDSC
             colnames(bq) <- c("Sample","B","Qdsc","dQ","f","angle","m4m3","m1m3", "Quality","Cook","alpha")
             bq <- bq[
              with(bq, order(B, Sample)),
             ]
             bq <- bq[c(1,2,3,11,4,5,6,7,8,10,9)]
             return(bq)
           }
    )
  })
  
  # Print Q values
  output$text1 <-renderTable({ 
    if (is.null(input$ifile)) { return() }
    return(table.calc())
    })
 
  output$download <- downloadHandler(filename = function(){
        paste("data-",Sys.Date(), ".csv", sep="")},
        content = function(file){ write.csv(table.calc(),file)
  })
  
  # Print the PI and statistical parameters
  output$summary <- renderPrint({
    if (is.null(input$ifile)) { return() }
    switch(input$Protocol,
     DB =
        {
        cat(" Protocol \t\tDB\n",
        "N \t \t \t",format(fit.DB()$res$ndataT,justify = "right"),"\n",
        "n \t \t \t",format(fit.DB()$res$ndataS,justify = "right"),"\n",
        "PI \t \t \t",format(fit.DB()$res$pi_val, digits = 3,nsmall=1),"\n",
        "68% Conf.Int \t \t[",format(fit.DB()$res$l68, digits =3,nsmall=1),format(fit.DB()$res$u68,digits = 3,nsmall=1),"]\n",
        "95% Conf.Int \t \t[",format(fit.DB()$res$l95, digits =3,nsmall=1),format(fit.DB()$res$u95,digits = 3,nsmall=1),"]\n",
        "r-squared \t \t",format(fit.DB()$res$r2,digits = 4,nsmall=4))
        },
     FC =
        {
        cat(" Protocol \t\tFC\n",
        "N \t \t \t",format(fit.FC()$res$ndataT,justify = "right"),"\n",
        "n \t \t \t",format(fit.FC()$res$ndataS,justify = "right"),"\n",
        "PI \t \t \t",format(fit.FC()$res$pi_val, digits = 3,nsmall=1),"\n",
        "68% Conf.Int \t \t[",format(fit.FC()$res$l68, digits =3,nsmall=1),format(fit.FC()$res$u68,digits = 3,nsmall=1),"]\n",
        "95% Conf.Int \t \t[",format(fit.FC()$res$l95, digits =3,nsmall=1),format(fit.FC()$res$u95,digits = 3,nsmall=1),"]\n",
        "r-squared \t \t",format(fit.FC()$res$r2,digits = 4,nsmall=4))
        },
     DSC =
        {
        cat(" Protocol \t\tDSC\n",
        "alpha \t \t \t",fit.DSC()$alpha_text,"\n",
        "N \t \t \t",format(fit.DSC()$res$ndataT,justify = "right"),"\n",
        "n \t \t \t",format(fit.DSC()$res$ndataS,justify = "right"),"\n",
        "PI \t \t \t",format(fit.DSC()$res$pi_val, digits = 3,nsmall=1),"\n",
        "68% Conf.Int \t \t[",format(fit.DSC()$res$l68, digits =3,nsmall=1),format(fit.DSC()$res$u68,digits = 3,nsmall=1),"]\n",
        "95% Conf.Int \t \t[",format(fit.DSC()$res$l95, digits =3,nsmall=1),format(fit.DSC()$res$u95,digits = 3,nsmall=1),"]\n",
        "r-squared \t \t",format(fit.DSC()$res$r2,digits = 4,nsmall=4))
        }
    )
  })    
  
  # Q vs B plot
  Plot1 <- reactive({
    if (is.null(input$ifile)) { return() }
    switch(input$Protocol,
     DB =
       {
       f <- fit.DB()
       msp_plot(f$QDB, data.frame(Field = f$best_fit$Field, fit = f$best_fit$DB),
                f$err_DB, f$res, ylab = "Qdb", anchor = FALSE)
       },
     FC =
       {
       f <- fit.FC()
       msp_plot(f$QFC, data.frame(Field = f$best_fit$Field, fit = f$best_fit$FC),
                f$err_FC, f$res, ylab = "Qfc")
       },
     DSC =
       {
       f <- fit.DSC()
       msp_plot(f$QDSC, data.frame(Field = f$best_fit$Field, fit = f$best_fit$DSC),
                f$err_DSC, f$res, ylab = "Qdsc", errbars = TRUE, note = f$plot_note)
       }
    )
  })

  output$QDB <- renderPlot({
    print(Plot1())
  })
  
  output$downplot <- downloadHandler(
      filename = function(){
        paste("MSP-plot",input$formatF,sep=".")},
      content = function(file){
        switch(input$formatF,
         png = {png(file, width = 700, height = 400, units ="px")},
         pdf = {pdf(file, width = 7, height = 4)},
         tiff = {tiff(file, width = 700, height = 400, units ="px")}
        )
        print(Plot1())
        dev.off()
      }
  )
  
  # ----------------------------------------------------------
  # Function for plotting cooks distances and threshold
  # ----------------------------------------------------------
  Plot2 <- reactive({
    if (is.null(input$ifile)) { return() }
    switch(input$Protocol,
    DB =
      {
       q <- fit.DB()$QDB
       cook_plot(q, 3 * mean(q$Cook))
      },
    FC =
      {
       q <- fit.FC()$QFC
       cook_plot(q, 3 * mean(subset(q, quality != "Bad")$Cook))
      },
    DSC =
      {
       q <- fit.DSC()$QDSC
       cook_plot(q, 3 * mean(subset(q, quality != "Bad")$Cook))
      }
    )
  })

  output$Pcook <- renderPlot({
    print(Plot2())
  })
  
  output$downplot2 <- downloadHandler(
    filename = function(){
      paste("Cook-plot",input$formatF2,sep=".")},
    content = function(file){
      switch(input$formatF2,
             png = {png(file, width = 700, height = 400, units ="px")},
             pdf = {pdf(file, width = 7, height = 4)},
             tiff = {tiff(file, width = 700, height = 400, units ="px")}
      )
      print(Plot2())
      dev.off()
    }
  )
  
  # ----------------------------------------------------------
  # Function for plotting the kernel density curve
  # ----------------------------------------------------------
  Plot3 <- reactive({
    if (is.null(input$ifile)) { return() }
    switch(input$Protocol,
      DB  = boot_plot(fit.DB()$dat,  fit.DB()$res, xr = input$xrangeBoot),
      FC  = boot_plot(fit.FC()$dat,  fit.FC()$res, xr = input$xrangeBoot),
      DSC = boot_plot(fit.DSC()$dat, fit.DSC()$res, fit.DSC()$plot_note, xr = input$xrangeBoot)
    )
  })

  output$Pboot <- renderPlot({
    print(Plot3())
  })
  
  output$downplot3 <- downloadHandler(
    filename = function(){
      paste("Boot-plot",input$formatF3,sep=".")},
    content = function(file){
      switch(input$formatF3,
             png = {png(file, width = 700, height = 400, units ="px")},
             pdf = {pdf(file, width = 7, height = 4)},
             tiff = {tiff(file, width = 700, height = 400, units ="px")}
      )
      print(Plot3())
      dev.off()
    }
  )
  
  ######################
  # Alpha values (cases 3 and 4) and case 4 vs case 3
  ######################

  # 'Alpha values': DSC, case 3, 4 or 5.
  # 'Case 4 vs case 3' and 'Lambda' (LOOCV): DSC, case 3 or 4.
  observe({
    isdsc  <- identical(input$Protocol, "DSC") && !is.null(input$caseDSC)
    show   <- isdsc && input$caseDSC %in% c("3", "4")
    show_a <- isdsc && input$caseDSC %in% c("3", "4", "5")
    if (show_a) showTab(inputId = "resTabs", target = "Alpha values")
    else        hideTab(inputId = "resTabs", target = "Alpha values")
    for (tb in c("Case 4 vs case 3", "Lambda")) {
      if (show) showTab(inputId = "resTabs", target = tb)
      else      hideTab(inputId = "resTabs", target = tb)
    }
    if (isdsc) showTab(inputId = "resTabs", target = "Compare cases")
    else       hideTab(inputId = "resTabs", target = "Compare cases")
    show35 <- isdsc && input$caseDSC %in% c("3", "5")
    if (show35) showTab(inputId = "resTabs", target = "Case 5 vs case 3")
    else        hideTab(inputId = "resTabs", target = "Case 5 vs case 3")
    if (isTRUE((!show_a && input$resTabs == "Alpha values") ||
               (!isdsc && input$resTabs == "Compare cases") ||
               (!show35 && input$resTabs == "Case 5 vs case 3") ||
               (!show && input$resTabs %in% c("Case 4 vs case 3", "Lambda")))) {
      updateTabsetPanel(inputId = "resTabs", selected = "Q values")
    }
  })

  # alpha value displayed in the left panel: fixed (case 1) or optimized (case 2)
  output$alphaInfo <- renderText({
    switch(as.character(input$caseDSC),
           "1" = sprintf("%.2f (fixed)", alpha_case1),
           "2" = sprintf("%.3f (optimized in [%g, %g])", fit.DSC.c2()$alpha_global, alpha_lo, alpha_hi),
           "")
  })

  # lambda selected by LOOCV for the current case (3 or 4), shown in the left panel
  output$lambdaInfo <- renderText({
    f <- switch(as.character(input$caseDSC), "3" = fit.DSC.c3(), "4" = fit.DSC.c4(), NULL)
    if (is.null(f)) return("")
    edge <- if (f$lambda %in% range(lambda_grid)) " - edge of the search range" else ""
    sprintf("%.2f (search range %.2f - %.2f%s)", f$lambda, min(lambda_grid), max(lambda_grid), edge)
  })

  # Default x-range of the bootstrap plot: range of the laboratory fields of the file
  observeEvent(meas(), {
    B  <- suppressWarnings(as.numeric(meas()$Field))
    lo <- max(0,   floor(min(B, na.rm = TRUE) / 5) * 5)
    hi <- min(100, ceiling(max(B, na.rm = TRUE) / 5) * 5)
    if (is.finite(lo) && is.finite(hi) && lo < hi) {
      updateSliderInput(inputId = "xrangeBoot", value = c(lo, hi))
    }
  })

  # lambda_i of case 5 (fixed by the intervals of the file), shown in the left panel
  output$lambdaInfo5 <- renderText({
    f <- fit.DSC.c5()
    sprintf("lambda_i = 1/sigma_i^2, from %.0f to %.0f", f$lambda_range[1], f$lambda_range[2])
  })

  # Table of the alpha values of the selected case (3, 4 or 5)
  alpha.table <- reactive({
    validate(need(identical(input$Protocol, "DSC") && input$caseDSC %in% c("3", "4", "5"),
                  "Select protocol DSC and case 3, 4 or 5."))
    f   <- fit.DSC()
    tab <- data.frame(Sample  = f$QDSC$Name,
                      B       = f$QDSC$x,
                      a_min   = f$bounds$a_min,
                      a_max   = f$bounds$a_max,
                      alpha   = f$QDSC$alpha,
                      Quality = f$QDSC$quality)
    tab[order(tab$B, tab$Sample), ]
  })

  output$alphaTable <- renderTable({
    alpha.table()
  }, digits = 3)

  output$downloadAlpha <- downloadHandler(
    filename = function() {
      paste0("alpha-case", input$caseDSC, "-", Sys.Date(), ".csv")},
    content = function(file) { write.csv(alpha.table(), file, row.names = FALSE) }
  )

  # ---------------------------------------------------------------------
  # Comparison of the five MSP-DSC cases (paired bootstrap).
  # Run on demand, and only when the five cases have been computed.
  # ---------------------------------------------------------------------
  output$compareStatus <- renderText({
    st <- c(done$c1, done$c2, done$c3, done$c4, done$c5)
    paste0("Cases computed: ",
           paste0(1:5, ifelse(st, " (yes)", " (no)"), collapse = ", "),
           if (all(st)) "" else "\nSelect each missing case in 'MSP-DSC: alpha strategy' (left panel) to compute it.")
  })

  compare.res <- eventReactive(input$runCompare, {
    validate(need(identical(input$Protocol, "DSC"), "Select protocol DSC."))
    validate(need(all(c(done$c1, done$c2, done$c3, done$c4, done$c5)),
                  "The five cases must be computed first: select each case 1 to 5 in 'MSP-DSC: alpha strategy'."))

    # the five fits are already computed (cached reactives)
    f <- list(fit.DSC1(), fit.DSC.c2(), fit.DSC.c3(), fit.DSC.c4(), fit.DSC.c5())

    m   <- meas()
    n   <- nrow(m)
    B   <- as.numeric(m$Field)
    sel <- dsc.sel()
    mm2 <- fc.base()$mm2
    Qfc <- fc.base()$B.QFC$y
    c_i <- 2 * (m$m1 - m$m3) / (2 * m$m0 - m$m1 - mm2)
    dQ  <- calc_dQ(m$m0, m$m1, m$m3, m$m4, mm2, alpha_dQ)
    use <- sel$good & is.finite(Qfc) & is.finite(c_i) & is.finite(dQ) & dQ > 0

    # common set of specimens: those passing the selection criteria
    U <- data.frame(id = seq_len(n), Field = B, m0 = m$m0, m1 = m$m1, m3 = m$m3,
                    m2_corr = mm2, Qfc = Qfc, dQ = dQ,
                    a_min = rep(alpha_lo, n), a_max = rep(alpha_hi, n))[use, , drop = FALSE]
    nU <- nrow(U)
    keep <- lapply(f, function(x) x$QDSC$quality[use] == "Good")

    # case 1 : constant alpha, weights 1/dQ (as in the fit of case 1)
    U1 <- U
    U1$dQ <- calc_dQ(m$m0, m$m1, m$m3, m$m4, mm2, alpha())[use]
    U1$Q1 <- (Qfc + alpha() * c_i)[use]
    fit1 <- function(db) { w <- 1 / db$dQ; sum(w * db$Field * (db$Q1 + 1)) / sum(w * db$Field^2) }
    # case 2 : single alpha re-optimized at each replicate
    fit2 <- function(db) {
      ob <- optimize(function(a) profile_slope(a, db)$SSR, interval = c(alpha_lo, alpha_hi))
      profile_slope(ob$minimum, db)$slope
    }
    # case 3 and 4 : lambda selected on the complete data set, kept fixed
    lam3 <- f[[3]]$lambda
    lam4 <- f[[4]]$lambda
    fit3 <- function(db) fit_block_alpha(db, lam3, alpha_pop)$slope
    U4 <- U
    U4$a_min <- f[[4]]$bounds$a_min[use]
    U4$a_max <- f[[4]]$bounds$a_max[use]
    fit4 <- function(db) fit_block_alpha(db, lam4, "mid")$slope
    # case 5 : prior (centre and lambda_i) of each specimen
    pri <- case5_prior(f[[5]]$bounds$a_min, f[[5]]$bounds$a_max, f[[5]]$k)
    U5 <- U
    U5$a_ctr <- pri$ctr[use]
    U5$lam_i <- pri$lam[use]
    fit5 <- function(db) fit_block_alpha(db, NA_real_, "col")$slope

    cases <- list(list(data = U1, keep = keep[[1]], fit = fit1),
                  list(data = U,  keep = keep[[2]], fit = fit2),
                  list(data = U,  keep = keep[[3]], fit = fit3),
                  list(data = U4, keep = keep[[4]], fit = fit4),
                  list(data = U5, keep = keep[[5]], fit = fit5))

    pi_boot <- withProgress(message = "Paired bootstrap of the five cases", value = 0, {
      paired_bootstrap_cases(nU, cases, R.rep(), progress = function(x) setProgress(value = x))
    })

    pi_hat <- sapply(f, function(x) x$res$pi_val)
    labs   <- paste("Case", 1:5)
    est <- data.frame(Case = labs,
                      PI = pi_hat,
                      CI95_low = sapply(f, function(x) x$res$l95),
                      CI95_up  = sapply(f, function(x) x$res$u95),
                      Specimens = sapply(f, function(x) x$res$ndataS),
                      stringsAsFactors = FALSE)
    diffs <- compare_cases(pi_boot, pi_hat, labs)
    diffs$Significant <- ifelse(diffs$p_Holm < 0.05, "yes", "no")
    list(est = est, diffs = diffs, n_rep = nrow(pi_boot),
         n_failed = colSums(!is.finite(pi_boot)))
  })

  output$compareEst <- renderTable({
    r <- compare.res()$est
    names(r) <- c("Case", "PI (uT)", "CI 95% low", "CI 95% up", "Specimens")
    r
  }, digits = 2)

  output$compareDiff <- renderTable({
    r <- compare.res()$diffs
    r <- data.frame(Comparison = r$Comparison,
                    "Difference (uT)" = r$Difference,
                    "CI 95% low" = r$CI_low, "CI 95% up" = r$CI_up,
                    "p" = signif(r$p_value, 3), "p (Holm)" = signif(r$p_Holm, 3),
                    "Significant (Holm, 5%)" = r$Significant,
                    check.names = FALSE, stringsAsFactors = FALSE)
    r
  }, digits = 2)

  output$compareNote <- renderText({
    r <- compare.res()
    paste0(r$n_rep, " paired bootstrap replicates; replicates with a failed fit per case: ",
           paste(r$n_failed, collapse = ", "), ".")
  })

  Plot6 <- reactive({
    r <- compare.res()$diffs
    r$Comparison <- factor(r$Comparison, levels = rev(r$Comparison))
    ggplot(r, aes(x = Difference, y = Comparison, colour = Significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_pointrange(aes(xmin = CI_low, xmax = CI_up), linewidth = 1, size = 1) +
      scale_colour_manual(values = c(yes = "darkred", no = "steelblue"),
                          breaks = c("yes", "no"),
                          name = "Significant\n(Holm, 5%)") +
      labs(x = "Difference of paleointensity (\u00b5T)", y = NULL) +
      theme_minimal() +
      theme(axis.title   = element_text(size = 20),
            axis.text    = element_text(size = 18),
            legend.title = element_text(size = 16),
            legend.text  = element_text(size = 14))
  })
  output$Pcompare <- renderPlot({ print(Plot6()) })

  output$downplot6 <- downloadHandler(
    filename = function() paste("Compare-cases", input$formatF6, sep = "."),
    content = function(file) {
      switch(input$formatF6,
             png = {png(file, width = 1000, height = 650, units = "px")},
             pdf = {pdf(file, width = 10, height = 6.5)},
             tiff = {tiff(file, width = 1000, height = 650, units = "px")})
      print(Plot6())
      dev.off()
    }
  )

  output$downloadCompare <- downloadHandler(
    filename = function() paste0("compare-cases-", Sys.Date(), ".csv"),
    content = function(file) write.csv(compare.res()$diffs, file, row.names = FALSE)
  )

  # Alpha of case 4 (constrained) vs alpha of case 3 (full range), specimens
  # whose alpha was actually optimized in both cases
  alpha.compare <- reactive({
    a3 <- fit.DSC.c3()$QDSC
    a4 <- fit.DSC.c4()$QDSC
    cmp <- data.frame(Sample = a3$Name, B = a3$x,
                      alpha_c3 = a3$alpha, alpha_c4 = a4$alpha,
                      q3 = a3$quality, q4 = a4$quality)
    cmp <- cmp[cmp$q3 == "Good" & cmp$q4 == "Good", ]
    validate(need(nrow(cmp) > 0, "No specimen is retained in both case 3 and case 4."))
    cmp
  })

  Plot4 <- reactive({
    validate(need(identical(input$Protocol, "DSC") && input$caseDSC %in% c("3", "4"),
                  "Select protocol DSC and case 3 or 4."))
    cmp <- alpha.compare()
    alpha_compare_plot(cmp, "alpha_c3", "alpha_c4",
                       expression(paste("Optimized ", alpha, " (Full range [0, 0.8])")),
                       expression(paste("Optimized ", alpha, " (MDF-constrained intervals)")))
  })

  output$Palpha <- renderPlot({
    print(Plot4())
  })

  output$downplot4 <- downloadHandler(
    filename = function(){
      paste("Alpha-case4-vs-case3",input$formatF4,sep=".")},
    content = function(file){
      switch(input$formatF4,
             png = {png(file, width = 600, height = 600, units ="px")},
             pdf = {pdf(file, width = 6, height = 6)},
             tiff = {tiff(file, width = 600, height = 600, units ="px")}
      )
      print(Plot4())
      dev.off()
    }
  )

  # Alpha of case 5 (prior on the individual intervals) vs alpha of case 3 (full
  # range), specimens whose alpha was optimized in both cases
  alpha.compare5 <- reactive({
    a5 <- fit.DSC.c5()$QDSC
    a3 <- fit.DSC.c3()$QDSC
    cmp <- data.frame(Sample = a5$Name, B = a5$x,
                      alpha_c3 = a3$alpha, alpha_c5 = a5$alpha,
                      q3 = a3$quality, q5 = a5$quality)
    cmp <- cmp[cmp$q3 == "Good" & cmp$q5 == "Good", ]
    validate(need(nrow(cmp) > 0, "No specimen is retained in both case 3 and case 5."))
    cmp
  })

  Plot7 <- reactive({
    validate(need(identical(input$Protocol, "DSC") && input$caseDSC %in% c("3", "5"),
                  "Select protocol DSC and case 3 or 5."))
    cmp <- alpha.compare5()
    alpha_compare_plot(cmp, "alpha_c3", "alpha_c5",
                       expression(paste("Optimized ", alpha, " (Full range [0, 0.8])")),
                       expression(paste("Optimized ", alpha, " (Gaussian prior on MDF-constrained intervals)")))
  })

  output$Palpha5 <- renderPlot({
    print(Plot7())
  })

  output$downplot7 <- downloadHandler(
    filename = function(){
      paste("Alpha-case5-vs-case3",input$formatF7,sep=".")},
    content = function(file){
      switch(input$formatF7,
             png = {png(file, width = 600, height = 600, units ="px")},
             pdf = {pdf(file, width = 6, height = 6)},
             tiff = {tiff(file, width = 600, height = 600, units ="px")}
      )
      print(Plot7())
      dev.off()
    }
  )

  # LOOCV error as a function of lambda (case 3 or 4), selected lambda = vertical line
  Plot5 <- reactive({
    validate(need(identical(input$Protocol, "DSC") && input$caseDSC %in% c("3", "4"),
                  "Select protocol DSC and case 3 or 4."))
    f <- fit.DSC()
    ggplot(f$loocv, aes(x = lambda, y = loocv_error)) +
      geom_line(colour = "steelblue", linewidth = 0.8) +
      geom_point(colour = "steelblue", size = 2.5) +
      geom_vline(xintercept = f$lambda, colour = "darkred", linewidth = 1) +
      # label next to the line, on the side where there is room
      annotate("text", x = f$lambda, y = Inf, label = sprintf("lambda == %.2f", f$lambda),
               parse = TRUE, vjust = 1.5, size = 6, colour = "grey20",
               hjust = if (f$lambda > mean(range(lambda_grid))) 1.15 else -0.15) +
      scale_x_continuous(breaks = seq(min(lambda_grid), max(lambda_grid), by = 0.25)) +
      labs(x = expression(lambda), y = "LOOCV error") +
      theme_minimal() +
      theme(axis.title = element_text(size = 20),
            axis.text  = element_text(size = 18))
  })

  output$Plambda <- renderPlot({
    print(Plot5())
  })

  output$downplot5 <- downloadHandler(
    filename = function(){
      paste0("Lambda-LOOCV-case", input$caseDSC, ".", input$formatF5)},
    content = function(file){
      switch(input$formatF5,
             png = {png(file, width = 700, height = 400, units ="px")},
             pdf = {pdf(file, width = 7, height = 4)},
             tiff = {tiff(file, width = 700, height = 400, units ="px")}
      )
      print(Plot5())
      dev.off()
    }
  )

  ######################
  # Readme Tab
  ######################
  
  output$downloadExample <- downloadHandler(
    filename <- function() {
      paste("example1", "csv", sep=".")
    },
      content <- function(file) {
      file.copy("data/example.csv", file)
    }
  )
}
