# this is the OpenMx script to fit the Bivariate SEM-PGS model. all the starting values for the off diagonal elements are set to 0.
# Author: Xuanyu Lyu
# Date: 05/31/2024

fitBiSEMPGS_m2 <- function(data_path){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance","1e-4")
    mxOption(NULL,"Optimality tolerance","1e-7")

    #mxOption(NULL,"Number of Threads","4")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        #VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,0.06,.84), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + 2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f) + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")
        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, name="VY_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        #VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')

    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.5,.6), label=c("a11", "a22"),    name="a", lbound = c(.2,.3))     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.02,0.02,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.03,0.03,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.10,0.15,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    # add a constraint to make k=j
        j_Algebra <- mxAlgebra(k, name = "j_Algebra")
        j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0.15,0.1,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.1) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.5 * (itlo + t(itol)), name="ic_Algebra") # ic should be full
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.05) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.05) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.05) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")

    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, 
                    #VF, 
                    VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, 
                    VF_Algebra,  Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    #VF_Constraint, 
                    #VE_Constraint,
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    j_constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

        fitModel1 <- mxTryHardWideSearch(Model1, extraTries = 30, OKstatuscodes = c(0,1), intervals=T, silent=F, showInits = F, exhaustive = F, jitterDistrib = "rnorm", loc=.5, scale = .1)
        return(summary(fitModel1))

}

# this is an Openmx fit function that assuming the genetic variances are known from other models (RDR, sibling)
fitBiSEMPGS_m2_fixH2 <- function(data_path){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance","1e-6")
    mxOption(NULL,"Number of Threads","4")
    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    # Load the simulated data for this demonstration:
        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.3,0.3,.2), label=c("VF11", "VF12", "VF12","VF22"), name="VF") # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.7,.1,0.1,.5), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(F,F), values=c(0.4949747468,0.3346640106), label=c("a11", "a22"),    name="a", lbound = -.05)     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.4,0.10,0.15,.2), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        #adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        #adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")
    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.10,0,-0.01,.25), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.1) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.15,0.1,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.15,0.1,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.15,0.1,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.1,0.1,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.5 * (itlo + t(itol)), name="ic_Algebra") # ic should be full
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.2) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.2) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.2) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")
    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, VF, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v, 
                    VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra,  gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    VF_Constraint, 
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

        #fitModel1 <- mxTryHard(Model1, extraTries = 5, intervals=T, silent=F)
        fitModel1 <- mxTryHardWideSearch(Model1, extraTries = 8, intervals=T, silent=F,showInits = T, exhaustive = T,  loc=.2, scale = .05)
        return(summary(fitModel1))

}

fitBiSEMPGS_m2_fixH2 <- function(data_path,feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 5, exhaustive = T){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))

    #mxOption(NULL,"Number of Threads","4")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        #VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,.06,.65), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        #VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')

    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.3,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(F,F), values=c(0.4949747468,0.3346640106), label=c("a11", "a22"),    name="a", lbound = -.05)     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.10,0.15,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    # add a constraint to make k=j
        j_Algebra <- mxAlgebra(k, name = "j_Algebra")
        j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0.15,0.1,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.1) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.5 * (itlo + t(itol)), name="ic_Algebra") # ic should be full
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.05) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.05) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.05) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")

    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY_Algebra,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY_Algebra,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY_Algebra,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, 
                    #VF, 
                    VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, 
                    VF_Algebra,  Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    #VF_Constraint, 
                    #VE_Constraint,
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    j_constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

        fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0), intervals=T, silent=F, showInits = F, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
        return(summary(fitModel1))

}

fitBiSEMPGS_m2_tol <- function(data_path,feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))

    #mxOption(NULL,"Number of Threads","4")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables
        #VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=FALSE, values=NA, label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        #VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,0.06,.4), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.2) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")
        #VE_Algebra <- mxAlgebra(VY - (2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra), name="VE_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        #VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
        #VE_Constraint    <- mxConstraint(VE == VE_Algebra,       name='VE_Constraint')

    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.2) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.8,.55), label=c("a11", "a22"),    name="a", lbound = c(-.2,-.2))     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.02,0.02,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.2)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,T,T,F),nrow = 2,ncol = 2), values=c(.5,0.03,0.03,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.3)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.2) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.10,0.15,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.2) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    # add a constraint to make k=j
        j_Algebra <- mxAlgebra(k, name = "j_Algebra")
        j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0,0.05,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.2) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.2)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.2)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.2)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.2)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.2) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.2)
        ic   <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic21", "ic12","ic22"), name="ic", lbound = -.2) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.5 * (itlo + t(itol)), name="ic_Algebra") # ic should be full
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.3) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.2) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.2) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")

    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY_Algebra,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY_Algebra,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY_Algebra,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, 
                    #VF, 
                    VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, 
                    VF_Algebra,  
                    #VE_Algebra,
                    Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    #VF_Constraint, 
                    #VE_Constraint,
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    j_constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)
        #print(mxCheckIdentification(Model1,details = FALSE)$status)
        fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, showInits = F, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
        #return(fitModel1)
        return(summary(fitModel1, verbose = TRUE))

}

fitBiSEMPGS_m2_tol_fixH2 <- function(data_path, avalue,feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
        # load the packages
    library(OpenMx)
    library(data.table)
    library(stringr)

    # Specify Options:
        mxOption(NULL,"Calculate Hessian","Yes")
        mxOption(NULL,"Standard Errors","Yes")
        mxOption(NULL,"Default optimizer","NPSOL")
        #mxOption(NULL,"mxByRow","TRUE")

    # some optimizer options - adapted from Yongkong's script
    
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))

    #mxOption(NULL,"Number of Threads","4")
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    #mxOption(NULL,"Analytic Gradients","No")

    options()$mxOptions$'Feasibility tolerance'
    #options()$mxOptions$'Analytic Gradients'
    options()$mxOptions$'Gradient step size'  #1e-7
    options()$mxOptions$'Optimality tolerance'  #1e-7
    #mxOption(NULL,"Analytic Gradients","No")

        Example_Data  <- fread(data_path, header = T)

        #cov(Example_Data, use="pairwise.complete.obs")

    # Create variables and define the algebra for each variables

        VY    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(2,.4,.4,1.5), label=c("VY11", "VY12", "VY12","VY22"), name="VY", lbound = -.05) # Phenotypic variance
        #VF    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.20,0.06,0.06,.04), label=c("VF11", "VF12", "VF12","VF22"), name="VF", lbound = -.1) # Variance due to VT
        VE    <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,.06,0.06,.84), label=c("VE11", "VE12", "VE12","VE22"), name="VE", lbound = -.05) # Residual variance

        VY_Algebra <- mxAlgebra(2 * delta %*% t(Omega) + 2 * a %*% t(Gamma) + w %*% t(delta) + v %*% t(a) + VF_Algebra + VE, name="VY_Algebra")
        VF_Algebra <- mxAlgebra(2 * f %*% VY %*% t(f) + f %*% VY %*% mu %*% VY %*% t(f) + f %*% VY %*% t(mu) %*% VY %*% t(f), name="VF_Algebra")

        VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
        #VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')

    # Genetic effects:
        delta <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.3), label=c("delta11", "delta22"),name="delta", lbound = -.05) # Effect of PGS on phen
        a     <- mxMatrix(type="Diag", nrow=2, ncol=2, free=c(F,F), values=avalue, label=c("a11", "a22"),    name="a", lbound = -.05)     # Effect of latent PGS on phen
        k     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("k11", "k12", "k12","k22"),    name="k", lbound = -.05)     # PGS variance (if no AM)
        j     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=matrix(c(F,F,F,F),nrow = 2,ncol = 2), values=c(.5,0.05,0.05,.5), label=c("j11", "j12", "j12","j22"),    name="j", lbound = -.05)     # Latent PGS variance (if no AM)
        Omega <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.6,0.15,0.1,.5), label=c("Omega11", "Omega21", "Omega12","Omega22"),name="Omega", lbound = -.05) # Within-person PGS-Phen covariance
        Gamma <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.5,0.10,0.15,.3), label=c("Gamma11", "Gamma21", "Gamma12","Gamma22"),name="Gamma", lbound = -.05) # Within-person latent PGS-Phen covariance

        Omega_Algebra <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + delta %*% k + 0.5 * w , name="Omega_Algebra") # E.g., cov(Yp, NTOp)
        Gamma_Algebra <- mxAlgebra(2 * a %*% hc + 2 * delta %*% t(ic) + a %*% j + 0.5 * v, name="Gamma_Algebra") # E.g., cov(Yp, NTLp)

        Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
        Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
        
        adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
        adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    # add a constraint to make k=j
        j_Algebra <- mxAlgebra(k, name = "j_Algebra")
        j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects:
        mu    <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.15,0,0.05,.3), label=c("mu11", "mu21", "mu12","mu22"), name="mu", lbound = -.1) # AM co-path coefficient
        gt     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.05,0.01,.1), label=c("gt11", "gt21", "gt12","gt22"),  name="gt", lbound = -.05)  # Increase in cross-mate PGS (co)variances from AM
        ht     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.015,0.01,.02), label=c("ht11", "ht21", "ht12","ht22"),  name="ht", lbound = -.05)  # Increase in cross-mate latent PGS (co)variances from AM
        gc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.14,0.1,0.1,.1), label=c("gc11", "gc12", "gc12","gc22"),   name="gc", lbound = -.05)  # Increase in within-mate PGS (co)variances from AM
        hc     <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(0.04,0.01,0.01,.015), label=c("hc11", "hc12", "hc12","hc22"),  name="hc", lbound = -.05)  # Increase in within-mate latent PGS (co)variances from AM
        gt_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Omega, name="gt_Algebra") # E.g., cov(TPO, TMO)
        ht_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Gamma, name="ht_Algebra") # E.g., cov(TPL, TML)
        gc_Algebra <- mxAlgebra(0.5 * (gt + t(gt)), name="gc_Algebra") # gc should be symmetric
        hc_Algebra <- mxAlgebra(0.5 * (ht + t(ht)), name="hc_Algebra") # hc should be symmetric
        gchc_constraint_Algebra <- mxAlgebra( hc * (2*delta%*%k%*%t(delta)/(2*a%*%j%*%t(a))), name = "gchc_constraint_Algebra") # g and h are equally proportional to a and delta

        itlo  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.05,0.04,.03), label=c("itlo11", "itlo21", "itlo12","itlo22"), name="itlo", lbound = -.05) 
        itol  <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.05,0.03,0.036,.03), label=c("itol11", "itol21", "itol12","itol22"), name="itol", lbound = -.05)
        ic   <- mxMatrix(type="Symm", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.07,0.05,0.05,.05), label=c("ic11", "ic12", "ic12","ic22"), name="ic", lbound = -.05) 
        itlo_Algebra <- mxAlgebra(t(Gamma) %*% mu %*% Omega, name="itlo_Algebra") # E.g., cov(TPO, TML)
        itol_Algebra <- mxAlgebra(t(Omega) %*% mu %*% Gamma, name="itol_Algebra") # E.g., cov(TPL, TMO)
        ic_Algebra <- mxAlgebra(.5 * (itlo + t(itol)), name="ic_Algebra") # ic should be full
        
        

        gt_constraint <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
        ht_constraint <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
        gc_constraint <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
        hc_constraint <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
        gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
        itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
        itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
        ic_constraint <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
        f     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.16,0.04,0.11,.09), label=c("f11", "f21","f12","f22"),  name="f", lbound = -.05) # Vertical Transmission effect
        w     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.3,0.15,0.1,.51), label=c("w11", "w21", "w12","w22"),  name="w", lbound = -.05) # Genetic nurture
        v     <- mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), values=c(.2,0.1,0.07,.08), label=c("v11", "v21", "v12","v22"),  name="v", lbound = -.05) # Latent nurture
        w_Algebra     <- mxAlgebra(2 * f %*% Omega + f %*% VY %*% mu %*% Omega + f %*% VY %*% t(mu) %*% Omega, name="w_Algebra")    
        v_Algebra     <- mxAlgebra(2 * f %*% Gamma + f %*% VY %*% mu %*% Gamma + f %*% VY %*% t(mu) %*% Gamma, name="v_Algebra")    
        wv_constraint_algebra <- mxAlgebra((w * sqrt(2*delta%*%k%*%t(delta))/sqrt(2*a%*%j%*%t(a))), name='wv_constraint_algebra')

        v_constraint <- mxConstraint(v == v_Algebra, name='v_constraint')
        w_constraint <- mxConstraint(w == w_Algebra, name='w_constraint')
        wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')
    # Between-people covariances
        thetaNT <- mxAlgebra(2 * delta %*% gc + 2 * a %*% ic + .5 * w, name="thetaNT")
        thetaT  <- mxAlgebra(1 * delta %*% k + thetaNT, name="thetaT")
        Yp_PGSm <- mxAlgebra(VY %*% mu %*% Omega, name="Yp_PGSm")
        Ym_PGSp <- mxAlgebra(VY %*% t(mu) %*% Omega, name="Ym_PGSp")
        Yp_Ym   <- mxAlgebra(VY %*% mu %*% VY,    name="Yp_Ym")
        Ym_Yp   <- mxAlgebra(VY %*% t(mu) %*% VY, name="Ym_Yp")
        Yo_Yp   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% t(mu) %*% VY + a %*% t(Gamma) %*% t(mu) %*% VY + f %*% VY + f %*% VY %*% t(mu) %*% VY, name = "Yo_Yp")
        Yo_Ym   <- mxAlgebra(delta %*% t(Omega) + a %*% t(Gamma) + delta %*% t(Omega) %*% mu %*% VY + a %*% t(Gamma) %*% mu %*% VY + f %*% VY + f %*% VY %*% mu %*% VY, name = "Yo_Ym")

    # Expected covariances matrix
        CovMatrix <- mxAlgebra(rbind(
            #     Yp1 Yp2|   Ym1 Ym2|   Yo1 Yo2|    Tp1 Tp2|   NTp1 NTp2|    Tm1 Tm2| NTm1 NTm2
            cbind(VY_Algebra,        Yp_Ym,     t(Yo_Yp),   Omega,     Omega,        Yp_PGSm, Yp_PGSm), #Yp1 Yp2
            cbind(Ym_Yp,     VY_Algebra,        t(Yo_Ym),   Ym_PGSp,   Ym_PGSp,      Omega,   Omega),   #Ym1 Ym2
            cbind(Yo_Yp,     Yo_Ym,     VY_Algebra,         thetaT,    thetaNT,      thetaT,  thetaNT), #Yo1 Yo2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaT),  k+gc,      gc,           gt,      gt),      #Tp1 Tp2
            cbind(t(Omega),  t(Ym_PGSp),t(thetaNT), gc,        k+gc,         gt,      gt),      #NTp1 NTp2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaT),  t(gt),     t(gt),        k+gc,    gc),      #Tm1 Tm2
            cbind(t(Yp_PGSm),t(Omega),  t(thetaNT), t(gt),     t(gt),        gc,      k+gc)),   #NTm1 NTm2
            dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

    # Expected means for all the variables:
        Means <- mxMatrix(type = "Full", nrow = 1, ncol = 14, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYp2", "meanYm1", "meanYm2", "meanYo1", "meanYo2", "meanTp1", "meanTp2", "meanNTp1", "meanNTp2", "meanTm1", "meanTm2", "meanNTm1", "meanNTm2"), 
        dimnames = list(NULL, c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2")),
        name = "expMeans")


    # put the mean and cov into a multivariate normal
        ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMeans", dimnames=c("Yp1", "Yp2", "Ym1", "Ym2", "Yo1", "Yo2", "Tp1", "Tp2", "NTp1", "NTp2", "Tm1", "Tm2", "NTm1", "NTm2"))
        
    # Convert data into a usable format for OpenMx:
        Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

    # Create fit function:
        FitFunctionML <- mxFitFunctionML()
        # Create fit function:
        #FitFunctionML <- mxFitFunctionWLS(type = "ULS", allContinuousMethod='marginals')
    # Specify what parameters we're going to be including in our model:
        Params <- list(
                    VY, 
                    #VF, 
                    VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                    VY_Algebra, 
                    VF_Algebra,  
                    Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                    VY_Constraint, 
                    #VF_Constraint, 
                    #VE_Constraint,
                    #Omega_Constraint, 
                    Gamma_Constraint, 
                    #adelta_Constraint,
                    j_constraint,
                    #gt_constraint, 
                    ht_constraint, 
                    #gc_constraint, 
                    hc_constraint, 
                    #gchc_constraint, 
                    itlo_constraint, 
                    itol_constraint, 
                    ic_constraint, 
                    v_constraint, 
                    w_constraint,
                    #wv_constraint,
                    thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                    CovMatrix, Means, ModelExpectations, FitFunctionML)
    # Create the model:
        options(warning.length = 8000)
        Model1 <- mxModel("BiSEM_PGS", Params, Example_Data_Mx)

        fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, showInits = F, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
        return(summary(fitModel1))

}

# The OpenMx function to fit the univariate model
fitUniSEMPGS <- function(data, max.cores = 2){
    library(OpenMx)
    library(data.table)
    library(stringr)

  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

#   options()$mxOptions$'Number of Threads'
#   options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  #mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Optimality tolerance","1e-7")
  #mxOption(NULL,"Analytic Gradients","No")

#   options()$mxOptions$'Feasibility tolerance'
#   #options()$mxOptions$'Analytic Gradients'
#   options()$mxOptions$'Gradient step size'  #1e-7
#   options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- data
  #dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel)]


  #Colnames
  colnames(dat_SEM) <- c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")
  nv=1

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS_lat)=.5 at t0
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.4,label="F",name="f",lbound=-.05)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.3,label="EnvPath",name="e",lbound=-.05)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.01,label="hap_PRS_cov",name="g")
  h <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.01,label="latent_hap_PRS_cov",name="h")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.3,label="obs_coef",name="delta",lbound=-.05)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.7,label="latent_coef",name="a",lbound=-.05)


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0.1,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.15,label="covAandF",name="w1")
  v1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.2,label="covAand_lF",name="v1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.05,label="copath",name="mu1")
  i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.02,label="covOL",name="i1")
  sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=1.5,label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.3,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.2, label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma1*(1 + sigma1*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega1*(1 + sigma1*mu1),name="w2")
  v2 <- mxAlgebra(2*f*Gamma1*(1 + sigma1*mu1),name="v2")
  i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(2*a*i1 + delta*k + 2*delta*g+.5*w1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra(2*a*h + a*j + 2*delta*i1 + .5*v1,name="Gamma2")
  #h <- mxAlgebra((g*a^2)/delta^2,name="h")
  h2 <- mxAlgebra(Gamma1^2*mu1,name="h2")


  # Equating nonlinear constraints and parameters
  xCon <- mxConstraint(x1==x2,name="xCon")
  wCon <- mxConstraint(w1==w2,name="wCon")
  vCon <- mxConstraint(v1==v2,name="vCon")
  iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  hCon <- mxConstraint(h==h2,name="hCon")
  sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma2 <- mxAlgebra(2*a*Gamma1 + 2*delta*Omega1 + a*v1 + delta*w1 + 2*f^2*sigma1*(1 + sigma1*mu1) + e^2,name="sigma2")
  ThetaNTM05 <- mxAlgebra((Omega1-delta*k),name="ThetaNTM05") #this is HALF theta_NTM in the math
  spsPRS <-mxAlgebra(sigma1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS
  PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigma1)*(1+mu1*sigma1),name="PO")
  spouse <-mxAlgebra(mu1*sigma1^2,name="spouse")

  v_Mean <- rep(0, 7)  
  names(v_Mean) <- c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Ym         Yp        Yo
    cbind(  k+g        ,g         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaNTM05   ),    #NTm
    cbind(  g          ,k+g       ,g          ,g        ,Omega1     ,spsPRS    ,Omega1        ),    #Tm
    cbind(  g          ,g         ,k+g        ,g        ,spsPRS     ,Omega1    ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,g          ,k+g      ,spsPRS     ,Omega1    ,Omega1     ),    #Tp
    cbind(  Omega1     ,Omega1    ,spsPRS     ,spsPRS   ,sigma1     ,spouse    ,PO         ),    #Ym
    cbind(  spsPRS     ,spsPRS    ,Omega1     ,Omega1   ,spouse     ,sigma1    ,PO         ),    #Yp
    cbind(  ThetaNTM05 ,Omega1    ,ThetaNTM05 ,Omega1   ,PO         ,PO        ,sigma1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE,  values= rep(0.1,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYm","meanYp","meanYo"),dimnames=list(NULL,c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm1","Tm1","NTp1","Tp1","Ym1","Yp1","Yo1"))

  params <- list(k ,j,  f , e , g , 
                 h, delta,a,
                 #x1 , 
                 w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                 #x2,
                 w2, v1,v2, i2, g2,
                 h2,
                 Omega2 , Gamma2, sigma2,
                 #xCon, 
                 wCon, 
                 iCon, 
                 gCon, 
                 OmegaCon, 
                 GammaCon, 
                 vCon, 
                 hCon, 
                 sigmaCon,
                 ThetaNTM05, 
                 spsPRS, PO, spouse,
                 CVmat, MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw",  means=v_Mean)
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=mxRun(modelAFE,intervals=TRUE,silent=FALSE)
  print("Model fitting done")
  # uncomment to print the openmx summary
  return(summary(AFE.Fit))

}