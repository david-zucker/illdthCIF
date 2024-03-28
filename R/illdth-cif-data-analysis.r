### Cumulative Incidence Function Estimation Based on Population-Based Biobank Data
### Malka Gorfine, David M. Zucker, and Shoval Shoham

### Code by David Zucker
### Version of 26 March 2024

### R code for applying to a dataset a new method
### for estimating the cumulative incidence function for disease
### in an illness-death model

library(etm)
library(Rcpp)

#SCROLL DOWN TO LINE 419 FOR MAIN FUNCTION

### FUNCTIONS ##################################################################

trans = function(u,t.opt) {
  if (t.opt == 1) ans = u
  if (t.opt == 2) ans = -log(1-u)
  if (t.opt == 3) ans = 0.5*pi - asin(sqrt(1-u))
  return(ans)
}

trans.deriv = function(u,t.opt) {
  didl = 1e-8
  if (t.opt == 1) ans = rep(1,length(u))
  if (t.opt == 2) ans = 1/(1-u)
  if (t.opt == 3) ans = 0.5 / (sqrt(u+didl)*sqrt(1-u+didl))
  return(ans)
}  

trans.inv = function(v,t.opt) {
  if (t.opt == 1) ans = v
  if (t.opt == 2) ans = 1-exp(-v)
  if (t.opt == 3) ans = 1 - (sin(0.5*pi-v)^2)
  return(ans)
}  

#KAPLAN-MEIER FUNCTION
KM = function(trecr,tfu,del,wt) {
  
  #IDENTIFY UNIQUE EVENT TIMES
  evtim = tfu[which(del==1)]
  time = unique(evtim)
  time = sort(time)
  ndist = length(time)

  #CALCULATIONS
  atrskvec = rep(0,ndist)
  km_surv = rep(0,ndist)
  skmcur = 1
  for (m in 1:ndist) {
    cur.tim = time[m] 
    d = sum(wt*(trecr < cur.tim)*(tfu==cur.tim)*del)
    atrsk = sum(wt*(trecr < cur.tim)*(tfu >= cur.tim))
    p = d/atrsk
    skmcur = skmcur * (1-p)
    atrskvec[m] = atrsk
    km_surv[m] = skmcur
  }  
  ans = list(time=time, surv=km_surv, atrsk=atrskvec)
  
  return(ans)

}
#END KAPLAN-MEIER FUNCTION

#C BACKEND FOR AJ ESTIMATOR
cppFunction('NumericVector x_aux_aj(
  int nsam, int n_aj, int nfst, NumericVector noncens, NumericVector prev, 
  NumericVector case_vec, NumericVector V1, NumericVector km_first_times,
  NumericVector ycal_circ, NumericVector ycal_circ_recip,
  NumericVector ycal_star, NumericVector ycal_star_recip,
  NumericVector sfac, NumericVector age_recr, int ntgrd, NumericVector tgrd,
  NumericVector age_first, NumericVector afun, NumericMatrix at_rsk_i,
  IntegerVector ixf_vec, int ncase, IntegerVector ix_case) {
	
NumericMatrix xa_mat(ntgrd,nsam);

double n_aj_d = n_aj;
double nsam_d = nsam;
double ppi = n_aj_d/nsam_d;

for (int it = 0; it < ntgrd; it++) {
  for (int ii = 0; ii < nsam; ii++) {
    xa_mat(it,ii) = 0;
  }
}  

for (int ic = 0; ic < ncase; ic++) {
  int ij = ix_case(ic) - 1;
  int ixfv_ij = ixf_vec(ij) - 1;
  double atrsk_cur = ycal_circ(ixfv_ij); 
	double atrsk_recip_cur = ycal_circ_recip(ixfv_ij);
  for (int ii = 0; ii < nsam; ii++) {
    if (prev(ii) == 0) {
      double trm1 = 0;
      double trm2 = 0;
      double trm3 = 0;
      double V3i = age_first(ii);
      if (noncens(ii) == 1) {
        int ixfv_ii = ixf_vec(ii) - 1;
        if (V3i < V1(ij)) {
           trm1 = sfac(ixfv_ii) * ycal_star_recip(ixfv_ii);
        }
      }
      for (int ixf = 0; ixf < nfst; ixf++) {
        double ftime = km_first_times[ixf];
        if ((ftime >= age_recr(ii)) && (ftime <= V3i) && (ftime < V1(ij))) {
          trm2 += sfac(ixf) * pow(ycal_star_recip(ixf),2);
        }
      }  
      trm2 = trm2 / n_aj;
      trm3 = atrsk_recip_cur*(at_rsk_i(ii,ixfv_ij) - atrsk_cur);
    	for (int it = 0; it < ntgrd; it++) {
      	if (V1(ij) <= tgrd(it)) {
	        xa_mat(it,ii) += -afun(ixfv_ij)*((trm1 - trm2)/ppi + trm3) / nsam;
      	}
    	}	
    }
  }    
}

return(xa_mat);

}')
#END C BACKEND FOR AJ ESTIMATOR
        
#FUNCTION FOR CIF COMPUTATION FOR AALEN-JOHANSEN ESTIMATOR
cifcmp.aj = function(age_recr, age_diag, age_death, age_end_fu, status_end, tgrd, wts) {
  
  #status_end
  #0 = alive without disease (i.e. censored)  
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease 
  
  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)
  
  #sample size  
  nsam = length(age_recr)

  #SETUPS
  noncens = (status_end != 0)
  ix.cen = which(status_end == 0)
  diseased = (status_end >= 2)
  prev = ((age_diag <= age_recr) & diseased)
  ix.prev = which(prev)
  case = (age_diag <= age_death) & (!prev) & (noncens)
  ix.case = which(case)
  ncase = length(ix.case)
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  cifest = rep(NA,ntgrd)
  cifest.sd1 = rep(NA,ntgrd)
  cifest.sd2 = rep(NA,ntgrd)
  age_first = pmin(age_diag, age_death, age_end_fu)
  
  #REMOVE PREVALENT CASES FOR COMPUTATION OF KM FOR FIRST TRANSITION
  nprev = length(ix.prev)
  n.aj = nsam - length(ix.prev)
  age_recr.aj = age_recr
  age_diag.aj = age_diag
  age_death.aj = age_death
  age_end_fu.aj = age_end_fu
  status_end.aj = status_end
  age_first.aj = age_first
  wts.aj = wts
  noncens.aj = noncens
  if (nprev > 0) {
    age_recr.aj = age_recr[-ix.prev]
    age_diag.aj = age_diag[-ix.prev]
    age_death.aj = age_death[-ix.prev]
    age_end_fu.aj = age_end_fu[-ix.prev]
    status_end.aj = status_end[-ix.prev]
    age_first.aj = age_first[-ix.prev]
    wts.aj = wts[-ix.prev]
    noncens.aj = noncens[-ix.prev]
  }  

  #KAPLAN-MEIER ESTIMATE FOR TIME TO FIRST TRANSITION
  #OMITTING PREVALENT CASES
  km_fst = KM(age_recr.aj, age_first.aj, noncens.aj, wts.aj)
  km_fst_times = km_fst$time
  km_fst_surv = km_fst$surv
  km_fst_surv1 = c(1,km_fst_surv)
  km_fst_at_rsk = km_fst$atrsk
  ycal.circ = km_fst_at_rsk/nsam
  ycal.star = km_fst_at_rsk/n.aj
  ycal.circ.recip = 1/ycal.circ
  ycal.star.recip = 1/ycal.star
  nfst = length(km_fst_times)
  sfac = km_fst_surv1[1:nfst]/km_fst_surv
  sfac[which(is.infinite(sfac))] = 0
  sfac[which(is.nan(sfac))] = 0
  afun = km_fst_surv1[1:nfst]/ycal.circ
  
  #PRELIMINARIES
  at_rsk_i = matrix(0,nsam,nfst)
  ixf_vec = rep(0,nsam)
  for (ix.f in 1:nfst) {
    icur = which((age_first==km_fst_times[ix.f]) & (!prev) & (noncens))
    ixf_vec[icur] = ix.f
    icur1 = which((age_recr < km_fst_times[ix.f]) &
      (age_first >= km_fst_times[ix.f]))
    at_rsk_i[icur1,ix.f] = 1
  }

  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam) 
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = afun[ixf_vec[ix.case]] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(x.main[ix.g,ix.case])/nsam 
    cifest.sd1[ix.g] = sqrt((sum(x.main[ix.g,ix.case]^2)/nsam - cifest[ix.g]^2)/nsam)  
  }
  
  #COMPUTATION OF TERMS GOING INTO IID REPRESENTATION FOR VARIANCE CALCULATION
  if (auxflg) {
    x.aux = x_aux_aj(nsam, n.aj, nfst, noncens, prev, case, age_diag, km_fst_times,
      ycal.circ, ycal.circ.recip, ycal.star, ycal.star.recip, sfac, age_recr, 
      ntgrd, tgrd, age_first, afun, at_rsk_i, ixf_vec, ncase, ix.case)
  }
  
  #WRAP UP IID REPRESENTATION AND VARIANCE CALCULATION
  eps = matrix(0,ntgrd,nsam)
  if (auxflg) {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] + x.aux[ix.g,] - cifest[ix.g]
      cifest.sd2[ix.g] = sqrt(mean(eps[ix.g,]^2)/nsam)
    }  
  } 
  else {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] - cifest[ix.g]
    }
    cifest.sd2 = cifest.sd1
  } 
 
  #COMPUTE ESTIMATOR USING etmCIF AND EXTRACT VARIANCE
  fstatus = rep(0,n.aj)
  fstatus[which(status_end.aj==1)] = 2
  fstatus[which(status_end.aj==2)] = 1
  fstatus[which(status_end.aj==3)] = 1
  mydat = data.frame(start=age_recr.aj,finish=age_first.aj,fstatus=fstatus)
  ajr = try(etmCIF(Surv(start,finish,fstatus != 0) ~ 1, data=mydat, etype=fstatus, failcode=1))
  ajr1 = summary(ajr, ci.fun='linear')
  ajr1 = ajr1[[1]]$'CIF 1'
  if (is.null(ajr1)) {
    aj_est_fin = cifest
    cifest.sd3=cifest.sd2
  }
  else {
  aj.time = ajr1$time
  aj.est = ajr1$P
  aj.sd = sqrt(ajr1$var)

  #ADAPT TO TIME GRID
  aj_est_fin = NULL
  cifest.sd3 = NULL
  ajcur = 0
  sdcur = 0
  for (ig in 1:ntgrd) {
    ixt = which(aj.time <= tgrd[ig])
    if (length(ixt)>0) {
      ixt1 = max(ixt)
      ajcur = aj.est[ixt1]
      sdcur = aj.sd[ixt1]
    }
    aj_est_fin = c(aj_est_fin, ajcur)
    cifest.sd3 = c(cifest.sd3, sdcur)
  }
  }
  
  #RETURN RESULT
  ans = list(cifest=cifest, cifest.sd1=cifest.sd1, cifest.sd2=cifest.sd2,
    cifest.sd3=cifest.sd3, eps=eps)
  return(ans)
  
}#END OF FUNCTION TO COMPUTE AJ ESTIMATOR

#C BACKEND FOR NEW ESTIMATOR
cppFunction('NumericVector x_aux_new(
  int nsam, int ndth, NumericVector died, NumericVector case_vec, NumericVector V1,
  NumericVector V2, NumericVector km_at_rsk_nrm, NumericVector km_at_rsk_nrm_recip,
  NumericVector sfac, NumericVector age_recr, int ntgrd, NumericVector tgrd, 
  NumericVector dth_vec, NumericVector afun, NumericMatrix at_rsk_i,
  IntegerVector ixd_vec, int ncase, IntegerVector ix_case) {
  
NumericMatrix xa_mat(ntgrd,nsam);

for (int ic = 0; ic < ncase; ic++) {
  int ij = ix_case(ic) - 1;
  int ixdv_ij = ixd_vec(ij) - 1;
  for (int ii = 0; ii < nsam; ii++) {
    double trm1 = 0;
    if ((died(ii) == 1) && (V2(ii) < V2(ij))) {
      int ixdv_ii = ixd_vec(ii) - 1;
      trm1 = sfac(ixdv_ii)*km_at_rsk_nrm_recip(ixdv_ii);
    }  
    double trm2 = 0;
    for (int ixd = 0; ixd < ndth; ixd++) {
      double dtime = dth_vec(ixd);
      if ((dtime >= age_recr(ii)) && (dtime <= V2(ii)) && (dtime < V2(ij))) {
        trm2 += sfac(ixd) * pow(km_at_rsk_nrm_recip(ixd),2);
      }       
    }
    trm2 = trm2 / nsam;
    double atrsk_cur = km_at_rsk_nrm(ixdv_ij); 
    double atrsk_recip_cur = km_at_rsk_nrm_recip(ixdv_ij);
    double trm3 = atrsk_recip_cur*(at_rsk_i(ii,ixdv_ij) - atrsk_cur);
    for (int it = 0; it < ntgrd; it++) {
      if (V1(ij) <= tgrd(it)) {
        xa_mat(it,ii) += -afun(ixdv_ij)*(trm1 - trm2 + trm3) / nsam;
      } 
    }         
  }    
}

return(xa_mat);
}')
#END C BACKEND FOR NEW ESTIMATOR 

#FUNCTION FOR CIF COMPUTATION FOR PROPOSED ESTIMATOR
cifcmp.new = function(age_recr, age_diag, age_death, age_end_fu, 
  status_end, tgrd, wts) {
  
  #status_end
  #0 = alive without disease (i.e. censored)  
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease 
  
  #sample size  
  nsam = length(age_recr)
  
  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)
  
  #SETUPS
  ix.cen = which(status_end == 0)
  died = (status_end == 1) | (status_end == 3)
  ix.died = which(died)
  case = (status_end == 3)
  ix.case = which(case)
  ncase = length(ix.case)
  cifest = rep(NA,ntgrd)
  cifest.sd1 = rep(NA,ntgrd)
  cifest.sd2 = rep(NA,ntgrd)
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  
  #KAPLAN-MEIER ESTIMATE OF DEATH TIME DISTN
  km_dth = KM(age_recr, age_end_fu, died, wts)
  km_dth_times = km_dth$time
  km_dth_surv = km_dth$surv
  km_dth_surv1 = c(1, km_dth_surv)
  km_at_rsk = km_dth$atrsk
  km_at_rsk_nrm = km_at_rsk/nsam   #estimate of script Y2 at death times
  km_at_rsk_nrm_recip = 1/km_at_rsk_nrm 
  ndth = length(km_dth_times)
  sfac = km_dth_surv1[1:ndth]/km_dth_surv
  sfac[which(is.infinite(sfac))] = 0
  sfac[which(is.nan(sfac))] = 0
  intwts.numer = -diff(km_dth_surv1)
  
  #PRELIMINARIES
  ixd_vec = rep(0,nsam)
  at_rsk_i = matrix(0,nsam,ndth)
  afun = rep(0,ndth)
  for (ix.d in 1:ndth) {
    icur = which(age_death == km_dth_times[ix.d])
    ixd_vec[icur] = ix.d
    ndth.cur.nrm = length(icur)/nsam
    afun[ix.d] = intwts.numer[ix.d] / ndth.cur.nrm 
    icur1 = which((age_recr < km_dth_times[ix.d]) &
      (age_death >= km_dth_times[ix.d]))
    at_rsk_i[icur1,ix.d] = 1
  } 
 
  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam) 
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = afun[ixd_vec[ix.case]] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(x.main[ix.g,ix.case])/nsam 
    cifest.sd1[ix.g] = sqrt((sum(x.main[ix.g,ix.case]^2)/nsam - cifest[ix.g]^2)/nsam)  
  }
  
  #COMPUTATION OF ADDITIONAL TERM GOING INTO IID REPRESENTATION FOR VARIANCE CALCULATION
  if (auxflg) {
    x.aux = x_aux_new(nsam, ndth, died, case, age_diag, age_death, km_at_rsk_nrm,
      km_at_rsk_nrm_recip, sfac, age_recr, ntgrd, tgrd, km_dth_times, afun, 
      at_rsk_i, ixd_vec, ncase, ix.case)
  }  
  
  #WRAP UP IID REPRESENTATION AND VARIANCE CALCULATION
  eps = matrix(0,ntgrd,nsam)
  if (auxflg) {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] + x.aux[ix.g,] - cifest[ix.g]
      cifest.sd2[ix.g] = sqrt(mean(eps[ix.g,]^2)/nsam)
    }  
  } 
  else {
    for (ix.g in 1:ntgrd) {
      eps[ix.g,] = x.main[ix.g,] - cifest[ix.g]
    }
    cifest.sd2 = cifest.sd1
  }  

  #RETURN RESULT
  ans = list(cifest=cifest, cifest.sd1=cifest.sd1, cifest.sd2=cifest.sd2, eps=eps)
  return(ans)
  
}
#END OF FUNCTION TO COMPUTE NEW ESTIMATOR

#FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS

#' Cumulative Incidence Function Estimation
#' 
#' @description
#' Function to implement the cumulative incidence function (CIF) estimator
#' of Gorfine, Zucker, and Shoham for the illness-death model
#' 
#' @param age_recr age at recruitment vector
#' @param age_diag age at disease diagnosis vector
#' @param age_death age at death vector
#' @param age_end_fu age at end of follow-up (death or censoring) vector (note that occurrence of disease does not end follow-up)
#' @param status_end status at end of follow-up vector (0=alive without disease, 1=died without disease, 2=alive with disease, 3=died with disease)
#' @param tgrd grid of ages at which the CIF will be computed (vector)
#' @param tgrd.cb grid of ages over which the simultaneous confidence band will be computed (vector)
#' @param auxflg logical flag for including (TRUE) or not including (FALSE) auxiliary terms
#' @param covpr desired confidence interval coverage probability
#' @param nresam number of bootstrap replications for confidence band
#' @param cmb.flg logical flag for including (TRUE) or not including (FALSE) combination estimator
#' @param adpwt logical flag for using (TRUE) or not using (FALSE) adaptive weights
#' @param cmb.wgt vector with length of length(tgrd) with weights for combination estimator (ignored if adpwt=TRUE)
#' @param t.opt option for transformation: 1=no transformation, 2=-log(1-u), 3=0.5*pi-arcsin(sqrt(1-u))
#' @param pt.ci option for pointwise CI's: 1=normal-theory, 2=bootstrap
#' @param bwts option for bootstrap weights: 1=normal, 2=exponential, 3=Poisson
#' 
#' @return a list containing CIF estimates, SD estimates, pointwise confidence limits, simultaneous confidence band limits
#' 
#' @export
cifcmp.full = function(age_recr, age_diag, age_death, age_end_fu, status_end, 
  tgrd, tgrd.cb, auxflg, covpr, nresam, cmb.flg, adpwt, cmb.wgt, t.opt, pt.ci, bwts) {
  
#ARGUMENTS
# age_recr = age at recruitment
# age_diag = age at diagnosis
# age_death = age at death
# age_end_fu = age at end of follow-up (death or censoring)
  #note that occurrence of disease does not end follow-up
# status_end = status at end of follow-up
  #0 = alive without disease (i.e. censored)  
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease
# tgrd = grid of ages at which the CIF will be computed
# tgrd.cb = grid of ages over which the simultaneous confidence band will be computed
# auxflg = logical flag for including (TRUE) or not including (FALSE) auxiliary terms
# covpr = desired confidence interval coverage probability
# nresam = number of bootstrap replications for confidence band
# cmb.flg = logical flag for including (TRUE) or not including (FALSE) combination estimator
# adpwt = logical flag for using (TRUE) or not using (FALSE) adaptive weights
# cmb.wgt = vector with length of length(tgrd) with weights for combination estimator 
  #ignored if adpwt=TRUE
#t.opt = option for transformation: 1=no transformation, 2=-log(1-u), 3=0.5*pi-arcsin(sqrt(1-u))
#pt.ci = option for pointwise CI's: 1=normal-theory, 2=bootstrap
#bwts = option for bootstrap weights: 1=normal, 2=exponential, 3=Poisson

#PRELIMINARIES
zcrit = qnorm(1-(1-covpr)/2)
didl = 1e-8 #a small number
n = length(age_recr)
sq.n = sqrt(n)
tgrd2 = tgrd.cb
ixcb2 = which(tgrd %in% tgrd2)

#AALEN-JOHANSEN ESTMATOR
print(noquote('Computing AJ estimator ...'))
cifest.aj = cifcmp.aj(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, rep(1,n))
aj.est = cifest.aj$cifest
aj.sd1 = cifest.aj$cifest.sd1
aj.sd2 = cifest.aj$cifest.sd2
aj.sd3 = cifest.aj$cifest.sd3
aj.est.tr = trans(aj.est,t.opt)
ajderv = trans.deriv(aj.est,t.opt)
aj.sd1.tr = aj.sd1*ajderv
aj.sd2.tr = aj.sd2*ajderv
aj.sd3.tr = aj.sd3*ajderv
eps.aj = cifest.aj$eps

#NEW ESTIMATOR
print(noquote('Computing new estimator ...'))
cifest.new = cifcmp.new(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, rep(1,n))
new.est = cifest.new$cifest
new.sd1 = cifest.new$cifest.sd1
new.sd2 = cifest.new$cifest.sd2
new.est.tr = trans(new.est,t.opt)
newderv = trans.deriv(new.est,t.opt)
new.sd1.tr = new.sd1*newderv
new.sd2.tr = new.sd2*newderv
eps.new = cifest.new$eps

#COMBINED ESTIMATOR
cmb.est = NULL
cmb.var = NULL
cmb.sd = NULL
cmb.ptwise.ci.width = NULL
cmb.ptwise.ci.lo = NULL
cmb.ptwise.ci.hi = NULL
cmb.wts = NULL
cmb.band.lo = NULL
cmb.band.hi = NULL
if (cmb.flg) {
  var.aj = aj.sd2^2
  var.new = new.sd2^2
  covar = rep(NA,ntgrd)
  for (ix.g in 1:ntgrd) {
    covar[ix.g] = cov(eps.aj[ix.g,],eps.new[ix.g,])/n
  }
  if (adpwt) {
    #adaptive weights
    wnumer = var.aj - covar
    wdenom = var.aj + var.new - 2*covar
    cmb.wts = ifelse(wdenom > 0, wnumer/wdenom, 1)
    cmb.wts = ifelse(cmb.wts < 1, cmb.wts, 1)
  }
  else {
    #fixed weights
    cmb.wts = cmb.wgt
  }
  wts1 = 1 - cmb.wts
  wts2 = cmb.wts
  cmb.est = wts1*aj.est + wts2*new.est
  cmb.var = (wts1^2)*var.aj + (wts2^2)*var.new + 2*wts1*wts2*covar
  cmb.sd = sqrt(pmax(cmb.var,1e-8))
  cmb.est.tr = trans(cmb.est,t.opt)
  cmbderv = trans.deriv(cmb.est,t.opt)
  cmb.sd.tr = cmb.sd*cmbderv
  eps.cmb = wts1*eps.aj + wts2*eps.new
}
  
#BOOTSTRAP RESAMPLING SCHEME
print(noquote('Bootstrap Reps for Confidence Bands ...'))
eps.mean.aj = apply(eps.aj,1,mean)
eps.mean.new = apply(eps.new,1,mean)
if (cmb.flg) {eps.mean.cmb = wts1*eps.mean.aj + wts2*eps.mean.new}  
if (pt.ci == 2) {
  aj.zstat.boot = matrix(NA,nresam,ntgrd)
  new.zstat.boot = matrix(NA,nresam,ntgrd)
  if (cmb.flg) {cmb.zstat.boot = matrix(NA,nresam,ntgrd)}
}
supvec.aj = rep(0,nresam)
supvec.new = rep(0,nresam)
if (cmb.flg) {supvec.cmb = rep(0,nresam)}
if (bwts == 1) {qqmat = matrix(rnorm(n*nresam),nresam,n)}
if (bwts == 2) {qqmat = matrix(rexp(n*nresam,1),nresam,n)}
if (bwts == 3) {qqmat = matrix(rpois(n*nresam,1),nresam,n)}
for (b in 1:nresam) {
  if (b %% 25 == 0) print(b)
  qq = qqmat[b,]
  qq = matrix(qq,1,n) %x% matrix(1,ntgrd,1)
  #AJ
  qe.aj = qq*eps.aj
  qe.mean.aj = apply(qe.aj,1,mean)
  dif.raw.aj = qe.mean.aj - eps.mean.aj
  aj.est.cur = aj.est + dif.raw.aj
  aj.est.cur = ifelse(aj.est.cur>0, aj.est.cur, 0)
  aj.est.cur = ifelse(aj.est.cur<1, aj.est.cur, 1)
  dif.aj = (trans(aj.est.cur,t.opt) - aj.est.tr) / 
    (aj.sd3*trans.deriv(aj.est.cur,t.opt)+didl)
  if (pt.ci == 2) {aj.zstat.boot[b,] = dif.aj}
  supvec.aj[b] = max(abs(dif.aj[ixcb2]))
  #NEW
  qe.new = qq*eps.new
  qe.mean.new = apply(qe.new,1,mean)
  dif.raw.new = qe.mean.new - eps.mean.new
  new.est.cur = new.est + dif.raw.new
  new.est.cur = ifelse(new.est.cur>0, new.est.cur, 0)
  new.est.cur = ifelse(new.est.cur<1, new.est.cur, 1)
  dif.new = (trans(new.est.cur,t.opt) - new.est.tr) / 
    (new.sd2*trans.deriv(new.est.cur,t.opt)+didl)
  if (pt.ci == 2) {new.zstat.boot[b,] = dif.new}
  supvec.new[b] = max(abs(dif.new[ixcb2]))
  #CMB
  if (cmb.flg) {
    qe.cmb = qq*eps.cmb
    qe.mean.cmb = apply(qe.cmb,1,mean)
    dif.raw.cmb = qe.mean.cmb - eps.mean.cmb
    cmb.est.cur = cmb.est + dif.raw.cmb
    cmb.est.cur = ifelse(cmb.est.cur>0, cmb.est.cur, 0)
    cmb.est.cur = ifelse(cmb.est.cur<1, cmb.est.cur, 1)
    dif.cmb = (trans(cmb.est.cur,t.opt) - cmb.est.tr) / 
      (cmb.sd*trans.deriv(cmb.est.cur,t.opt)+didl)
    if (pt.ci == 2) {cmb.zstat.boot[b,] = dif.cmb}
    supvec.cmb[b] = max(abs(dif.cmb[ixcb2]))
  }
}

#POINTWISE CONFIDENCE INTERVALS
if (pt.ci == 1) {
  aj.pt.crit = zcrit
  new.pt.crit = zcrit
  if (cmb.flg) {cmb.pt.crit = zcrit}
}
else {
  aj.pt.crit = apply(abs(aj.zstat.boot), 2, quantile, probs=covpr)
  new.pt.crit = apply(abs(aj.zstat.boot), 2, quantile, probs=covpr)
  if (cmb.flg) {cmb.pt.crit = apply(abs(cmb.zstat.boot), 2, quantile, probs=covpr)}
}
#AJ
aj.ptwise.ci.hwd.tr.1 = aj.pt.crit*aj.sd1.tr
aj.ptwise.ci.lo.1 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.hi.1 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.width.1 = aj.ptwise.ci.hi.1 - aj.ptwise.ci.lo.1
aj.ptwise.ci.hwd.tr.2 = aj.pt.crit*aj.sd2.tr
aj.ptwise.ci.lo.2 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.hi.2 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.1, t.opt)
aj.ptwise.ci.width.2 = aj.ptwise.ci.hi.2 - aj.ptwise.ci.lo.2
aj.ptwise.ci.hwd.tr.3 = aj.pt.crit*aj.sd3.tr
aj.ptwise.ci.lo.3 = trans.inv(aj.est.tr - aj.ptwise.ci.hwd.tr.3, t.opt)
aj.ptwise.ci.hi.3 = trans.inv(aj.est.tr + aj.ptwise.ci.hwd.tr.3, t.opt)
aj.ptwise.ci.width.3 = aj.ptwise.ci.hi.3 - aj.ptwise.ci.lo.3
#NEW
new.ptwise.ci.hwd.tr = new.pt.crit*new.sd2.tr
new.ptwise.ci.lo = trans.inv(new.est.tr - new.ptwise.ci.hwd.tr, t.opt)
new.ptwise.ci.hi = trans.inv(new.est.tr + new.ptwise.ci.hwd.tr, t.opt)
new.ptwise.ci.width = new.ptwise.ci.hi - new.ptwise.ci.lo
#COMBO
if (cmb.flg) {
  cmb.ptwise.ci.hwd.tr = cmb.pt.crit*cmb.sd.tr
  cmb.ptwise.ci.lo = trans.inv(cmb.est.tr - cmb.ptwise.ci.hwd.tr, t.opt)
  cmb.ptwise.ci.hi = trans.inv(cmb.est.tr + cmb.ptwise.ci.hwd.tr, t.opt)
  cmb.ptwise.ci.width = cmb.ptwise.ci.hi - cmb.ptwise.ci.lo
}

#CONFIDENCE BANDS
#AJ
aj.band.crit = quantile(supvec.aj,probs=covpr,na.rm=T)
aj.band.hwd.tr = aj.band.crit*aj.sd3[ixcb2]*ajderv[ixcb2]
aj.band.lo = trans.inv(aj.est.tr[ixcb2] - aj.band.hwd.tr, t.opt)
aj.band.hi = trans.inv(aj.est.tr[ixcb2] + aj.band.hwd.tr, t.opt)
aj.band.width = aj.band.hi - aj.band.lo
#NEW
new.band.crit = quantile(supvec.new,probs=covpr,na.rm=T)
new.band.hwd.tr = new.band.crit*new.sd2[ixcb2]*newderv[ixcb2]
new.band.lo = trans.inv(new.est.tr[ixcb2] - new.band.hwd.tr, t.opt)
new.band.hi = trans.inv(new.est.tr[ixcb2] + new.band.hwd.tr, t.opt)
new.band.width = new.band.hi - new.band.lo
#CMB
if (cmb.flg) {
  cmb.band.crit = quantile(supvec.cmb,probs=covpr,na.rm=T)
  cmb.band.hwd.tr = cmb.band.crit*cmb.sd[ixcb2]*cmbderv[ixcb2]
  cmb.band.lo = trans.inv(cmb.est.tr[ixcb2] - cmb.band.hwd.tr, t.opt)
  cmb.band.hi = trans.inv(cmb.est.tr[ixcb2] + cmb.band.hwd.tr, t.opt)
  cmb.band.width = cmb.band.hi - cmb.band.lo
}  
else {
  cmb.band.width = 0
  cmb.band.lo = 0
  cmb.band.hi = 0
}

ans = list(
  aj.est = aj.est,
  aj.sd1 = aj.sd1,
  aj.sd2 = aj.sd2,
  aj.sd3 = aj.sd3,
  aj.ptwise.ci.width.1 = aj.ptwise.ci.width.1,
  aj.ptwise.ci.lo.1 = aj.ptwise.ci.lo.1,
  aj.ptwise.ci.hi.1 = aj.ptwise.ci.hi.1, 
  aj.ptwise.ci.width.2 = aj.ptwise.ci.width.2,
  aj.ptwise.ci.lo.2 = aj.ptwise.ci.lo.2,
  aj.ptwise.ci.hi.2 = aj.ptwise.ci.hi.2, 
  aj.ptwise.ci.width.3 = aj.ptwise.ci.width.3,
  aj.ptwise.ci.lo.3 = aj.ptwise.ci.lo.3,
  aj.ptwise.ci.hi.3 = aj.ptwise.ci.hi.3, 
  new.est = new.est,
  new.sd1 = new.sd1,
  new.sd2 = new.sd2,
  new.ptwise.ci.width = new.ptwise.ci.width,
  new.ptwise.ci.lo = new.ptwise.ci.lo,
  new.ptwise.ci.hi = new.ptwise.ci.hi,
  cmb.est = cmb.est,
  cmb.sd1 = cmb.sd,
  cmb.sd2 = cmb.sd,
  cmb.wts = cmb.wts,
  cmb.ptwise.ci.width = cmb.ptwise.ci.width,
  cmb.ptwise.ci.lo = cmb.ptwise.ci.lo,
  cmb.ptwise.ci.hi = cmb.ptwise.ci.hi,
  aj.band.width = aj.band.width,
  aj.band.lo = aj.band.lo,
  aj.band.hi = aj.band.hi,
  new.band.width = new.band.width,
  new.band.lo = new.band.lo,
  new.band.hi = new.band.hi,
  cmb.band.width = cmb.band.width,
  cmb.band.lo = cmb.band.lo,
  cmb.band.hi = cmb.band.hi)
#ALTOGETHER 35 OUTPUT ITEMS

return(ans)

}
#END OF FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS
