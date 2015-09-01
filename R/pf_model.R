require(EpiModel)
#' 
#' @title Pair Formation Model - pf_model
#'
#' @details 
#' @description This file contains the function pf_model which implements a
#' pair-formation using the "DCM" custom model framework.  The latter
#' uses the deSolve package to solve compartmental epidemic models of 
#' infectious diseases.  In this case, the model allows for two groups and
#' three disease states: susceptible, in primary stage or chronic stage of
#  infections.  To model this population, it uses  six individual compartments
#' (male or female for each disease state) and nine pair-wise compartments
#' (covering all combinations of disease states).
#'
#' The model assumes purely heterogenous mixing (e.g., two sexes with purely heterosexual contact).
#' Pair-wise contacts are primarily monogamous, but the model allows for
#' additional acts outside of the partnership at rates specified by gender.
#' These are captured with the concur.rate and concur.rate.g2 parameters.
#'
#' The top level arguments as for all "DCM" models are as follows:
#' @param param an \code{EpiModel} object of class \code{\link{param.dcm}}.
#' @param init an \code{EpiModel} object of class \code{\link{init.dcm}}.
#' @param control an \code{EpiModel} object of class \code{\link{control.dcm}}.
#'
#' A special wrapper function, "init.pf", is provided for convenience.  It takes
#' the overall population counts by disease state and gender (s.num. s.num.g2,
#' etc.) and a "pairs.fract" statistic for the fraction of the population 
#' in paired relationships.  It then populates the 15 individual and pair-wise 
#' compartments for presentation #' to the "DCM" model as a init.dcm object.
#'
#' The parameters for init.pf are as follows:
#' @param pairs.fract initial fraction of the population in pairs of any type. 
#' @param s.num number of susceptible females in the population
#' @param i.num number of chronically infected females in the population. 
#' @param s.num.g2 number of susceptible males in the population
#' @param i.num.g2 number of chronically infected males in the population. 
#' @param p.num number of primary stage infected females in the population*
#' @param p.num.g2 number of primary stage infected males in the population*
#' (Note: The primary stage infection parameters are not required.)
#'
#' The wrapper is optional.  The user may also explicitly pass the following 
#' compartment arguments to dcm (in this order):
#' 
#' Single female compartments:s1.num, p1.num, i1.num,
#' Single male compartments: s1.num.g2, p1.num.g2, i1.num.g2,
#' Paired compartments: p_ss.num, p_smpf.num, p_smif.num, p_pmsf.num,
#'    p_pp.num, p_pmif.num, p_imsf.num, p_impf.num and p_ii.num
#' 
#' The param object shares some fields common to two-group, "SI" models, but
#' adds some additional fields to account for interactions within pairs and 
#' for the two disease sttes: primary and chronic.
#' The common parameters are:
#' @param trans.rate probability of transmission given an act or contact between
#'   a susceptible and an infected person in the population when the female is 
#'   the susceptible partner.
#' @param trans.rate.g2 probability of transmission given an act or contact between
#'   a susceptible and an infected person in the population when the male is 
#'   the susceptible partner.
#' @param act.rate average number of acts governing transmission per person  
#'   per unit time, between singles, regardless of disease status.
#'
#' The new parameters are:
#' @param formation.rate pairing rate for singles in the population 
#' @param separate.rate desolution rate for pairs in the population 
#' @param concur.rate rate of acts outside pair for females in pairs 
#' @param concur.rate.g2 rate of acts outside pair for males in pairs 
#' @param pfactor multiplier for increasing infectivity during primary stage
#' @param p_act.rate average number of acts per unit time within pairs.
#' @param pi.rate  rate of transition from primary to chronic infection 
#'
#' The vital parameters are:
#' @param b.rate or bs.rate birth/immigration rate for susceptible females into the population.
#' @param b.rate.g2 or bs.rate.g2 birth/immigration rate for susceptible males into the population.
#' @param bp.rate same for primary infected females into the population. (opt)
#' @param bp.rate.g2 same for primary infected males into the population. (opt)
#' @param bi.rate same for chronically infected females into the population.
#' @param bi.rate.g2 same for chronically infected males into the population.
#' @param ds.rate mortality rate for susceptible females
#' @param ds.rate.g2 mortality rate for susceptible males
#' @param dp.rate mortality rate for females during primary infection
#' @param dp.rate.g2 mortality rate for males during primary infection
#' @param di.rate mortality rate for females during chronic infection
#' @param di.rate.g2 mortality rate for males during chronic infection
#'
#'
#' @details 
#' ## Example: PF Model with special "init.pf" parameters
#'
#'source("pf_model.R")
#'
#'init<-init.pf(s.num=500000,
#'              s.num.g2=600000,
#'              i.num=37500,
#'              i.num.g2=60000,
#'              pairs.fract=0.4)
#'
#'param <-param.dcm(
#'              formation.rate=0.025,
#'              separate.rate=0.01,
#'              concur.rate=0.05,
#'              concur.rate.g2=0.2,
#'              b.rate=0.006,
#'              b.rate.g2=NA,
#'              ds.rate=0.0024,
#'              ds.rate.g2=0.0024,
#'              di.rate=0.012,
#'              di.rate.g2=0.012,
#'              trans.rate = 0.04,
#'              trans.rate.g2 = 0.01, 
#'	        balance="g1",
#'              p_act.rate = 4.0, 
#'              pi.rate = 0.25,
#'	        pfactor=5.,
#'              act.rate = 1.0)
#'
#'control<-control.dcm(new.mod=pf_model,
#'	      type = "SI",
#'            nsteps = 240,
#'	      dt=0.1)
#'
#'mod<-dcm(param,init,control)
#'
#' init.pf is the optional wrapper function for correctly initializing the
#' 15 compartments and packaging them as an init.dcm object.
#'
init.pf <- function(s.num, 
                     i.num, 
                     s.num.g2, 
                     i.num.g2, 
                     pairs.fract,
                     ...) {

      # derived totals required for initial conditions
      fs.num <- s.num 
      ms.num <- s.num.g2
      if (exists("p.num")) {
	fp.num <- p.num
      } else {
	fp.num <- 0
      }
      if (exists("p.num.g2")) {
	mp.num <- p.num.g2
      } else {
	mp.num <- 0
      }
      fi.num <- i.num
      mi.num <- i.num.g2

      m.num <- ms.num + mp.num + mi.num
      f.num <- fs.num + fp.num + fi.num
      num <- m.num + f.num

      pairs.num = num * pairs.fract / 2

      # pair compartments by disease state
      p_ss.num = pairs.num * (ms.num/m.num) * (fs.num/f.num)
      p_smpf.num = pairs.num * (ms.num/m.num) * (fp.num/f.num)
      p_smif.num = pairs.num * (ms.num/m.num) * (fi.num/f.num)
      p_pmsf.num = pairs.num * (mp.num/m.num) * (fs.num/f.num)
      p_pp.num = pairs.num * (mp.num/m.num) * (fp.num/f.num)
      p_pmif.num = pairs.num * (mp.num/m.num) * (fi.num/f.num)
      p_imsf.num = pairs.num * (mi.num/m.num) * (fs.num/f.num)
      p_impf.num = pairs.num * (mi.num/m.num) * (fp.num/f.num)
      p_ii.num = pairs.num * (mi.num/m.num) * (fi.num/f.num)

      msingles.num = m.num - pairs.num
      fsingles.num = f.num - pairs.num
      singles.num = msingles.num + fsingles.num

      # single compartments by disease state

      s1.num = fsingles.num * (fs.num/f.num)
      p1.num = fsingles.num * (fp.num/f.num)
      i1.num = fsingles.num * (fi.num/f.num)
      s1.num.g2 = msingles.num * (ms.num/m.num)
      p1.num.g2 = msingles.num * (mp.num/m.num)
      i1.num.g2 = msingles.num * (mi.num/m.num)

      # Initial conditions
      out <- as.list(c(s1.num = s1.num,
		p1.num = p1.num,
		i1.num = i1.num,
		s1.num.g2 = s1.num.g2,
		p1.num.g2 = p1.num.g2,
		i1.num.g2 = i1.num.g2,
		p_ss.num = p_ss.num,
		p_smpf.num = p_smpf.num,
		p_smif.num = p_smif.num,
		p_pmsf.num = p_pmsf.num,
		p_pp.num = p_pp.num,
		p_pmif.num = p_pmif.num,
		p_imsf.num = p_imsf.num,
		p_impf.num = p_impf.num,
		p_ii.num = p_ii.num))

      class(out) <- "init.dcm"
      return(out)
}
#'
#' pf_model - implements pair-formation model to be used with "DCM" framework
#'
pf_model <- function(t, t0, parms) {
  
  with(as.list(c(t0, parms)), {
    
      if (!exists("verbose")) {
	verbose=FALSE
      }

      # when not modeling primary infection, set pfactor to 1 (for SI and IS to II transitions)
      if (pi.rate < 0.001) {
	  pfactor = 1.0
      }
      # derived totals
      p_pairs.num = p_ss.num + p_smpf.num + p_smif.num + 
	      p_pmsf.num + p_pp.num + p_pmif.num + 
	      p_imsf.num + p_impf.num + p_ii.num

      fsingles.num = s1.num + p1.num + i1.num
      msingles.num = s1.num.g2 + p1.num.g2 + i1.num.g2
      singles.num = msingles.num + fsingles.num
      
      s2.num = p_ss.num + p_pmsf.num + p_imsf.num
      p2.num = p_pp.num + p_smpf.num + p_impf.num
      i2.num = p_ii.num + p_pmif.num + p_smif.num

      s2.num.g2 = p_ss.num + p_smpf.num + p_smif.num
      p2.num.g2 = p_pp.num + p_pmif.num + p_pmsf.num
      i2.num.g2 = p_ii.num + p_impf.num + p_imsf.num

      s.num=s1.num+s2.num
      s.num.g2=s1.num.g2+s2.num.g2
      p.num=p1.num+p2.num
      p.num.g2=p1.num.g2+p2.num.g2
      i.num=i1.num+i2.num
      i.num.g2=i1.num.g2+i2.num.g2

      num.g1 = s.num + p.num + i.num
      num.g2 = s.num.g2 + p.num.g2 + i.num.g2

      t_num = num.g1 + num.g2

      # contact rate for single females (excludes contacts w/ paired males)
      if (balance == "g1") {
	ct.g1 <- act.rate
	ct.g2 <- ct.g1 * fsingles.num / msingles.num

	pair.f <- formation.rate
	pair.m <- pair.f * fsingles.num / msingles.num
      }
      if (balance == "g2") {
	ct.g2 <- act.rate.g2
	ct.g1 <- ct.g2 * msingles.num / fsingles.num

	pair.m <- formation.rate.g2
	pair.f <- pair.m * msingles.num / fsingles.num
      }
      
      if (verbose) {
	print ("Single compartments")
	print(c(s1.num,s1.num.g2,p1.num,p1.num.g2,i1.num,i1.num.g2))
	print ("Pair compartments")
	print(c(p_ss.num,p_smpf.num,p_smif.num,p_pmsf.num,p_pp.num,p_pmif.num,p_imsf.num,p_impf.num,p_ii.num))
	print ("Totals")
	print(c(num.g1,num.g2,t_num))
      }
      
      # birth rates
      if (exists("b.rate")) {
	# if b.rate.g2 is "NA" just use female population for births
	if (is.na(b.rate.g2)) {
	  brs.f <- b.rate*num.g1
	  brs.m <- b.rate*num.g1
	  bs.rate=b.rate
	  bs.rate.g2=b.rate
	} else {
	  brs.f <- b.rate*num.g1
	  brs.m <- b.rate.g2*num.g2
	  bs.rate=b.rate
	  bs.rate.g2=b.rate.g2
	}
      } else {
	if (! exists("bs.rate")) {
	  brs.f <- 0
	  brs.m <- 0
	  bs.rate=0
	  bs.rate.g2=0
	}
	# if bs.rate.g2 is "NA" just use female population for births
	if (is.na(bs.rate.g2)) {
	  brs.f <- bs.rate*num.g1
	  brs.m <- bs.rate*num.g1
	  bs.rate.g2=bs.rate
	} else {
	  brs.f <- bs.rate*num.g1
	  brs.m <- bs.rate.g2*num.g2
	}
      }

      if (! exists ("bp.rate")) {
	bp.rate <- 0
      }
      brp.f <- bp.rate*num.g1

      if (! exists ("bp.rate.g2")) {
	bp.rate.g2 <- bp.rate
      }
      brp.m <- bp.rate.g2*num.g2

      if (! exists ("bi.rate")) {
	bi.rate <- 0
      }
      bri.f <- bi.rate*num.g1

      if (! exists ("bi.rate.g2")) {
	bi.rate.g2 <- bi.rate
      }
      bri.m <- bi.rate.g2*num.g2
      
      # death rates
      if (! exists ("ds.rate")) {
	ds.rate <- 0
      }
      if (! exists ("ds.rate.g2")) {
	ds.rate.g2 <- ds.rate
      }
      if (! exists ("dp.rate")) {
	dp.rate <- 0
      }
      if (! exists ("dp.rate.g2")) {
	dp.rate.g2 <- dp.rate
      }
      if (! exists ("di.rate")) {
	di.rate <- 0
      }
      if (! exists ("di.rate.g2")) {
	di.rate.g2 <- di.rate
      }
      if (verbose) {
	print ("Demographic rates")
	print(c(brs.f,brs.m,ds.rate,ds.rate.g2,dp.rate,dp.rate.g2,di.rate,di.rate.g2))
      }

      # ODEs for singles  

      # Demographic flows
      # birth, death and progression

      dem_dSm <- brs.m - ds.rate.g2*s1.num.g2 +
	      ds.rate*p_ss.num + dp.rate*p_smpf.num  + di.rate*p_smif.num 
      dem_dSf <- brs.f - ds.rate*s1.num +
	      ds.rate.g2*p_ss.num + dp.rate.g2*p_pmsf.num  + di.rate.g2*p_imsf.num 
      dem_dPm <- brp.m -dp.rate.g2*p1.num.g2 +
	      ds.rate*p_pp.num + ds.rate*p_pmsf.num  + ds.rate*p_pmif.num 
      dem_dPf <- brp.f -dp.rate*p1.num +
	      dp.rate.g2*p_pp.num + ds.rate.g2*p_smpf.num  + di.rate.g2*p_smpf.num 
      dem_dIm <- bri.m -di.rate.g2*i1.num.g2 +
	      di.rate*p_ii.num + ds.rate*p_imsf.num  + dp.rate*p_impf.num 
      dem_dIf <- bri.f -di.rate*i1.num +
	      di.rate.g2*p_ii.num + ds.rate.g2*p_smif.num  + dp.rate.g2*p_pmif.num 

      if (verbose) {
	print ("Single Demographic flows")
	print(c(dem_dSm,dem_dSf,dem_dPm,dem_dPf,dem_dIm,dem_dIf))
      }
      
      # Pair-wise Demographic flows
      dem_dSS <- -ds.rate*p_ss.num -ds.rate.g2*p_ss.num 
      dem_dSmPf <- -dp.rate*p_smpf.num - ds.rate.g2*p_smpf.num
      dem_dSmIf <- -di.rate*p_smif.num - ds.rate.g2*p_smif.num
      dem_dPmSf <- -ds.rate*p_pmsf.num - dp.rate.g2*p_pmsf.num
      dem_dPP   <- -dp.rate*p_pp.num - dp.rate.g2*p_pp.num
      dem_dPmIf <- -dp.rate.g2*p_pmif.num - di.rate*p_pmif.num
      dem_dImSf <- -di.rate.g2*p_imsf.num - ds.rate*p_imsf.num
      dem_dImPf <- -di.rate.g2*p_impf.num - dp.rate*p_impf.num
      dem_dII   <- -di.rate*p_ii.num - di.rate.g2*p_ii.num

      if (verbose) {
	print ("Pair Demographic flows")
	print(c(dem_dSS,dem_dSmPf,dem_dSmIf,dem_dPmSf,dem_dPP,
	  dem_dPmIf,dem_dImSf,dem_dImPf,dem_dII))
      }

      # Pair Formation flows
      #    formation (marriage) and separation (divorce)

      pf_dSm <- - pair.m * s1.num.g2 +
	      separate.rate * (p_ss.num + p_smif.num + p_smpf.num)
      pf_dSf <- - pair.f * s1.num +
	      separate.rate * (p_ss.num + p_imsf.num + p_pmsf.num)
      pf_dPm <- - pair.m * p1.num.g2 +
	      separate.rate * (p_pp.num + p_pmsf.num + p_pmif.num)
      pf_dPf <- - pair.f * p1.num + 
	      separate.rate * (p_pp.num + p_smpf.num + p_impf.num)
      pf_dIm <- - pair.m * i1.num.g2 +
	      separate.rate * (p_ii.num + p_impf.num + p_imsf.num)
      pf_dIf <- - pair.f * i1.num + 
	      separate.rate * (p_ii.num + p_pmif.num + p_smif.num)

      if (verbose) {
	print ("Single Pair formation flows")
	print(c(pf_dSm,pf_dSf,pf_dPm,pf_dPf,pf_dIm,pf_dIf))
      }
      
      # check to see that changes due to pairing for singles are balanced between sexes
      if(abs((pf_dSm+pf_dPm+pf_dIm)-(pf_dSf+pf_dPf+pf_dIf)) > 0.001){
	print ("Singles change from Pairing imbalance")
	print(pf_dSm+pf_dPm+pf_dIm)
	print(pf_dSf+pf_dPf+pf_dIf)
	exit(1)
      }

      # Pair Formation flows (to/from paired compartments)
      pf_dSS <- pair.f * s1.num*s1.num.g2/(msingles.num) - separate.rate * p_ss.num
      pf_dSmPf <- pair.f * p1.num*s1.num.g2/(msingles.num) - separate.rate * p_smpf.num
      pf_dSmIf <- pair.f * i1.num*s1.num.g2/(msingles.num) - separate.rate * p_smif.num
      pf_dPmSf <- pair.f * s1.num*p1.num.g2/(msingles.num) - separate.rate * p_pmsf.num
      pf_dPP   <- pair.f * p1.num*p1.num.g2/(msingles.num) - separate.rate * p_pp.num
      pf_dPmIf <- pair.f * i1.num*p1.num.g2/(msingles.num) - separate.rate * p_pmif.num
      pf_dImSf <- pair.f * s1.num*i1.num.g2/(msingles.num) - separate.rate * p_imsf.num
      pf_dImPf <- pair.f * p1.num*i1.num.g2/(msingles.num) - separate.rate * p_impf.num
      pf_dII   <- pair.f * i1.num*i1.num.g2/(msingles.num) - separate.rate * p_ii.num

      if (verbose) {
	print ("Pair Formation/Separation flows")
	print(c(pf_dSS,pf_dSmPf,pf_dSmIf,pf_dPmSf,pf_dPP,
	  pf_dPmIf,pf_dImSf,pf_dImPf,pf_dII))
      }

      # check to see that all changes pair formation / separation balance
      if(abs((pf_dSm+pf_dSf+pf_dPm+pf_dPf+pf_dIm+pf_dIf)+
	    2*(pf_dSS+pf_dSmPf+pf_dSmIf+pf_dPmSf+pf_dPP+
	    pf_dPmIf+pf_dImSf+pf_dImPf+pf_dII)) > 0.001) {
	print ("Pair imbalance")
	print(c(t,pf_dSS,pf_dSmPf,pf_dSmIf,pf_dPmSf,pf_dPP,
	  pf_dPmIf,pf_dImSf,pf_dImPf,pf_dII))
	print(pf_dSm+pf_dSf+pf_dPm+pf_dPf+pf_dIm+pf_dIf)
	print(2*(pf_dSS+pf_dSmPf+pf_dSmIf+pf_dPmSf+pf_dPP+
	    pf_dPmIf+pf_dImSf+pf_dImPf+pf_dII))
	exit(1)
      }


      # Singles Infection flows

      # Coital matchups of singles w/ potential for infection
      female1S.rate = trans.rate*((ct.g1*i1.num.g2/msingles.num)+pfactor*(ct.g1*p1.num.g2/msingles.num))
      if (p_pairs.num > 0)
	female1S.rate = female1S.rate + trans.rate*((concur.rate.g2*(p_imsf.num+p_impf.num+p_ii.num)/p_pairs.num)+pfactor*(concur.rate.g2*(p_pmif.num+p_pmsf.num+p_pp.num)/p_pairs.num))
      male1S.rate = trans.rate.g2*(((ct.g2*i1.num/fsingles.num))+pfactor*(ct.g2*p1.num/fsingles.num))
      if (p_pairs.num > 0)
	male1S.rate = male1S.rate + trans.rate.g2*((concur.rate*(p_smif.num+p_pmif.num+p_ii.num)/p_pairs.num)+pfactor*(concur.rate*(p_smpf.num+p_impf.num+p_pp.num)/p_pairs.num))

      inf_dSm <- -male1S.rate*s1.num.g2
      inf_dSf <- -female1S.rate*s1.num
      if (pi.rate > 0) {
	inf_dPm <- male1S.rate*s1.num.g2- pi.rate*p1.num.g2
	inf_dPf <- female1S.rate*s1.num- pi.rate*p1.num
	inf_dIm <- pi.rate*p1.num.g2
	inf_dIf <- pi.rate*p1.num
      } else {
	inf_dPm <- 0
	inf_dPf <- 0
	inf_dIm <- male1S.rate*s1.num.g2
	inf_dIf <- female1S.rate*s1.num
      }

      if (verbose) {
	print ("Single infectionflows")
	print(c(inf_dSm,inf_dSf,inf_dPm,inf_dPf,inf_dIm,inf_dIf))
      }
      
      # check to see that all changes due to infection state in singles add to 0
      if(abs(inf_dSm+inf_dSf+inf_dPm+inf_dPf+inf_dIm+inf_dIf) > 0.001) {
	print ("Singles infection imbalance")
	print (inf_dSm+inf_dSf+inf_dPm+inf_dPf+inf_dIm+inf_dIf)
	exit(1)
      }

      # group-specific foi (w/o counts)
      # rates when infection must comes from outside
      fSStoSP.rate = trans.rate*concur.rate*
	      (i1.num.g2/msingles.num + pfactor*p1.num.g2/msingles.num)
      mSStoPS.rate = trans.rate.g2*concur.rate.g2*
	      (i1.num/fsingles.num + pfactor*p1.num/fsingles.num)

      # rates when infection can comes from outside or within pair
      fSPtoPP.rate = fSStoSP.rate+trans.rate*p_act.rate*pfactor
      mSPtoPP.rate = mSStoPS.rate+trans.rate.g2*p_act.rate*pfactor

      fSItoPI.rate = fSStoSP.rate+trans.rate*p_act.rate
      mSItoPI.rate = mSStoPS.rate+trans.rate.g2*p_act.rate

      if (pi.rate > 0) {
	SS2SP.flow = fSStoSP.rate*(1 - mSStoPS.rate)*p_ss.num
	SS2PS.flow = mSStoPS.rate*(1 - fSStoSP.rate)*p_ss.num
	SS2PP.flow = mSStoPS.rate*fSStoSP.rate*p_ss.num
	SP2PP.flow = mSPtoPP.rate*(1-pi.rate)*p_smpf.num
	SP2PI.flow = mSPtoPP.rate*pi.rate*p_smpf.num
	SP2SI.flow = (1-mSPtoPP.rate)*pi.rate*p_smpf.num
	PS2PP.flow = fSStoSP.rate*(1-pi.rate)*p_pmsf.num
	PS2IS.flow = (1-fSStoSP.rate)*pi.rate*p_pmsf.num
	PS2IP.flow = fSStoSP.rate*pi.rate*p_pmsf.num
	PP2PI.flow = pi.rate*(1-pi.rate)*p_pp.num
	PP2IP.flow = pi.rate*(1-pi.rate)*p_pp.num
	PP2II.flow = pi.rate*pi.rate*p_pp.num
	IS2IP.flow = pi.rate*p_imsf.num
	SI2PI.flow = pi.rate*p_smif.num
	IP2II.flow = pi.rate*p_impf.num
	PI2II.flow = pi.rate*p_pmif.num
	inf_dSS <- -(SS2SP.flow + SS2PS.flow + SS2PP.flow)
	inf_dSmPf <- SS2SP.flow - SP2SI.flow - SP2PP.flow - SP2PI.flow
	inf_dSmIf <- SP2SI.flow - SI2PI.flow
	inf_dPmSf <- SS2PS.flow - PS2IS.flow - PS2PP.flow - PS2IP.flow
	inf_dPP   <- PS2PP.flow + SS2PP.flow + SP2PP.flow - PP2IP.flow - PP2PI.flow - PP2II.flow
	inf_dPmIf <- PP2PI.flow + SP2PI.flow + SI2PI.flow - PI2II.flow
	inf_dImSf <- PS2IS.flow - IS2IP.flow
	inf_dImPf <- IS2IP.flow + PS2IP.flow + PP2IP.flow - IP2II.flow
	inf_dII <- PI2II.flow + PP2II.flow + IP2II.flow
	msp.flow <- male1S.rate*s1.num.g2 + SS2PS.flow+SP2PP.flow+SI2PI.flow
	fsp.flow <- female1S.rate*s1.num + SS2SP.flow+PS2PP.flow+IS2IP.flow
	mpi.flow <- pi.rate*p1.num.g2 + PS2IS.flow+PP2IP.flow+PI2II.flow
	fpi.flow <- pi.rate*p1.num + SP2SI.flow+PP2PI.flow+IP2II.flow
      } else {
	SS2SI.flow = fSStoSP.rate*(1 - mSStoPS.rate)*p_ss.num
	SS2IS.flow = mSStoPS.rate*(1 - fSStoSP.rate)*p_ss.num
	SS2II.flow = mSStoPS.rate*fSStoSP.rate*p_ss.num
	SI2II.flow = mSPtoPP.rate*p_smif.num
	IS2II.flow = fSPtoPP.rate*p_imsf.num

	inf_dSS <- -SS2II.flow - SS2SI.flow - SS2IS.flow
	inf_dSmIf <- SS2SI.flow - SI2II.flow
	inf_dImSf <- SS2IS.flow - IS2II.flow
	inf_dII   <- SS2II.flow + SI2II.flow + IS2II.flow
	inf_dSmPf <- 0
	inf_dPmIf <- 0
	inf_dPmSf <- 0
	inf_dImPf <- 0
	inf_dPP <- 0
	msp.flow <- male1S.rate*s1.num.g2 + SS2IS.flow+SI2II.flow
	fsp.flow <- female1S.rate*s1.num + SS2SI.flow+IS2II.flow
	mpi.flow <- 0
	fpi.flow <- 0
      }


      if (verbose) {
	print ("Pair Infection flows")
	print(c(inf_dSS,inf_dSmPf,inf_dSmIf,inf_dPmSf,inf_dPP,
	  inf_dPmIf,inf_dImSf,inf_dImPf,inf_dII))
      }

      # check to see that all changes due to infection state in pairs add to 0
      if(abs(inf_dSS+inf_dSmPf+inf_dSmIf+inf_dPmSf+inf_dPP+
	  inf_dPmIf+inf_dImSf+inf_dImPf+inf_dII)>0.001) {
	print ("Pair infection imbalance")
	exit(1)
      }

      # Total singles flows
      dSm <- dem_dSm + pf_dSm + inf_dSm
      dPm <- dem_dPm + pf_dPm + inf_dPm
      dIm <- dem_dIm + pf_dIm + inf_dIm

      dSf <- dem_dSf + pf_dSf + inf_dSf
      dPf <- dem_dPf + pf_dPf + inf_dPf
      dIf <- dem_dIf + pf_dIf + inf_dIf

      if (verbose) {
	print ("Single flows")
	print(c(dSm,dSf,dPm,dPf,dIm,dIf))
      }

      # ODEs for pairs  


      # Total pair flows
      dSS <- dem_dSS + pf_dSS + inf_dSS
      dSmPf <- dem_dSmPf + pf_dSmPf + inf_dSmPf
      dSmIf <- dem_dSmIf + pf_dSmIf + inf_dSmIf
      dPmSf <- dem_dPmSf + pf_dPmSf + inf_dPmSf
      dPP <- dem_dPP + pf_dPP + inf_dPP
      dPmIf <- dem_dPmIf + pf_dPmIf + inf_dPmIf
      dImSf <- dem_dImSf + pf_dImSf + inf_dImSf
      dImPf <- dem_dImPf + pf_dImPf + inf_dImPf
      dII <- dem_dII + pf_dII + inf_dII

      if (verbose) {
	print ("Pair flows")
	print(c(dSS,dSmPf,dSmIf,dPmSf,dPP,
	  dPmIf,dImSf,dImPf,dII))
      }

      list(c(dSf, dPf, dIf,dSm, dPm, dIm, dSS,dSmPf,dSmIf,dPmSf,dPP,dPmIf,dImSf,dImPf,dII), 
	s1.num = s1.num,
	p1.num = p1.num,
	i1.num = i1.num,
	s1.num.g2 = s1.num.g2,
	p1.num.g2 = p1.num.g2,
	i1.num.g2 = i1.num.g2,
	p_ss.num = p_ss.num,
	p_smpf.num = p_smpf.num,
	p_smif.num = p_smif.num,
	p_pmsf.num = p_pmsf.num,
	p_pp.num = p_pp.num,
	p_pmif.num = p_pmif.num,
	p_imsf.num = p_imsf.num,
	p_impf.num = p_impf.num,
	p_ii.num = p_ii.num,
	p_sp.num = p_smpf.num+p_pmsf.num,
	p_si.num = p_smif.num+p_imsf.num,
	p_pi.num = p_pmif.num+p_impf.num,
	p_pairs.num = p_pairs.num,
	prev_p_ss.num = 100*p_ss.num/p_pairs.num,
	prev_p_smpf.num = 100*p_smpf.num/p_pairs.num,
	prev_p_smif.num = 100*p_smif.num/p_pairs.num,
	prev_p_pmsf.num = 100*p_pmsf.num/p_pairs.num,
	prev_p_pp.num = 100*p_pp.num/p_pairs.num,
	prev_p_pmif.num = 100*p_pmif.num/p_pairs.num,
	prev_p_imsf.num = 100*p_imsf.num/p_pairs.num,
	prev_p_impf.num = 100*p_impf.num/p_pairs.num,
	prev_p_ii.num = 100*p_ii.num/p_pairs.num,
	prev_p_sp.num = 100*(p_smpf.num+p_pmsf.num)/p_pairs.num,
	prev_p_si.num = 100*(p_smif.num+p_imsf.num)/p_pairs.num,
	prev_p_pi.num = 100*(p_pmif.num+p_impf.num)/p_pairs.num,
	sp.flow = msp.flow,
	sp.flow.g2 = fsp.flow,
	pi.flow = mpi.flow,
	pi.flow.g2 = fpi.flow,
	s2.num = s2.num,
	p2.num = p2.num,
	i2.num = i2.num,
	s2.num.g2 = s2.num.g2,
	p2.num.g2 = p2.num.g2,
	i2.num.g2 = i2.num.g2,
	s.num=s.num,
	s.num.g2=s.num.g2,
	t_s.num=s.num+s.num.g2,
	prev_s.num.g1=100*s.num/num.g1,
	prev_s.num.g2=100*s.num.g2/num.g2,
	prev_s.num=100*(s.num+s.num.g2)/t_num,
	prev_p.num.g1=100*p.num/num.g1,
	prev_p.num.g2=100*p.num.g2/num.g2,
	prev_p.num=100*(p.num+p.num.g2)/t_num,
	prev_i.num.g1=100*i.num/num.g1,
	prev_i.num.g2=100*i.num.g2/num.g2,
	prev_i.num=100*(i.num+i.num.g2)/t_num,
	p.num=p.num,
	p.num.g2=p.num.g2,
	t_p.num=p.num+p.num.g2,
	i.num=i.num,
	i.num.g2=i.num.g2,
	t_i.num=i.num+i.num.g2,
	num.g1=num.g1, num.g2=num.g2, t_num=t_num,
	b.flow = bs.rate + bp.rate + bi.rate,
	ds.flow = ds.rate * s.num,
	dp.flow = dp.rate * p.num,
	di.flow = di.rate * i.num,
	b.flow.g2 = bs.rate.g2 + bp.rate.g2 + bi.rate.g2,
	ds.flow.g2 = ds.rate.g2 * s.num.g2,
	dp.flow.g2 = dp.rate.g2 * p.num.g2,
	di.flow.g2 = di.rate.g2 * i.num.g2)
    })
}
