#' Salmon Mortality Framework Model
#' Neil Banas, Emma Tyldesley, Graeme Diack, Colin Bull
#' v0.8, March 2023:   replacing fry Ricker curve with a density-independent mortality
#'                     replacing parr Beverton-Holt with Ricker
#' v0.8.1, May 2023:   update to default parameters to improve model response to 2SW population NB

mortalityFramework <- function(p = list(
  N_initial = 2.4e6, # initial number of eggs
  
  # --- life history schedule ---
  
  yearday_eggDeposition = -60, # Matlab: datenum('1 Nov 0000') - 366
  baselineDuration_egg = 5, # in months
  yearday_endOfFry = 274, # Matlab: datenum('30 Sep 0000')
  baselineDuration_parr = 18,
  flexibleParrDuration = 0, # if 1, choose parr duration based on length at end of fry stage, as opposed to this being set by the user via parameter values
  baselineDuration_smolt = 1,
  baselineDuration_earlyPS = 3,
  baselineDuration_latePS = 5,
  baselineDuration_adultOc = 4,
  baselineDuration_adultCoastal = 2,
  baselineDuration_adultRiver = 4,
  
  # --- growth parameters ---
  W_egg = 1, # 1 gram
  Q10_dtegg = 6.5,   # Elliott and Hurley 1998, reciprocal of eq 1a:
  # 2.12x change over 4 degC
  exp_growth = 0.31, # Forseth et al. 2001
  # growth rates below are based on round numbers from previous versions, with a 
  # factor of 1.18 in freshwater and 0.98 in the ocean applied based on smolt and adult
  # length, R Bush 1980s
  gmaxFry = 0.0084, # 0.007 gives 7 cm at start of parr; can vary this +/- 30% to give 6-8 cm parr
  gmaxParr6 = 0.0181, # 0.0152 gives a 13 cm smolt for 6 mo parr
  gmaxParr18 = 0.0065, # 0.0054 gives a 13.5 cm smolt for 18 mo parr
  gmaxParr30 = 0.0040, # 0.0034 gives a 14 cm smolt for 30 mo parr
  # for comparison, c/100 in the Ratkowsky model used by Forseth et al. 2001
  # (growth at 1 g body weight and at optimal temp.) for the "mod. fast" category
  # would give 0.0205. These tuned gmax values should always be less than this,
  # to reflect the fact that food is actually quite seasonal, not year-round, and
  # that temperature can't ever be better than optimal
  gmaxOc1SW = 0.050, # 0.051 turns a 13-13.5 cm smolt into a 60 cm adult after 1SW
  gmaxOc2SW = 0.030, # 0.031 turns a 13-13.5 cm smolt into a 75 cm adult after 2SW
  ref_length_parr = 7,   # smaller than this at start of parr stage, add 12 mos to parr stage, if flexibleParrDuration = 1
  ref_length_earlyPS = 14, # just for scaling the equations, not tuning targets
  ref_length_adultRiver = 60,
  L3overW = 62^3 / 2500, # length in cm cubed over weight in g
  # calibrated using Bacon et al. 2009
  
  # --- mortality params ---
  m_egg = 0.1,
  m_fry = 0.95,
  m_smolt = 0.2, # 0.1 - 0.5
  parr_ricker_alpha = 0.9028,
  parr_ricker_beta = 9.244e-6,
    # smolt = alpha * parr * exp(-beta * parr)
    # these multipliers are there to make the egg-smolt Ricker relationship
    # come out to a fit to Bush data, Feb 2023:
    # egg-smolt alpha = 0.026 (0.015, 0.037)
    # egg-smolt beta = 4.16e-07 (2.25e-07, 6.06e-07)
    #	p.parr_ricker_alpha = 0.026 / (1-p.m_egg) / (1-p.m_fry) / (1-p.m_smolt)
    #     / (1 - p.mort_parr_annual) , if the egg-smolt data mainly reflects 18 mo parr 
    # p.parr_ricker_beta = 4.16e-07 / (1-p.m_egg) / (1-p.m_fry);
    #
    # note that there is a danger of alpha coming out greater than 1, which doesn't
	  # make biological sense. Compare (1-p.m_egg) * (1-p.m_fry) * (1-p.m_smolt) with
	  # 0.026 to keep this from happening
  mort_parr_annual = 0.2, # additional mortality if the parr take 18 mo instead of 6
  m_earlyPS_monthly = 0.37, # at ref_length_earlyPS; declines rapidly with size
  exp_sizeMort = -0.35, # dependence of daily mortality on weight
  rmort2SW = 1.09, # additional marine mortality (multiplier) for 2SW vs 1SW
    # chosen to produce a small but nonzero equilibrium population,
    # with marine survival around 1.7%
  m_adultOc_monthly = 0.03,
  m_adultCoastal = 0.1,
  m_adultRiver = 0.09, # R Bush 1980s has 0.4
 
  # --- fecundity parameters ---
  # Fecundity estimated as function of fork length (L_f) in cm:
  # log10(eggs)=m.log10(L_f)+c
  # Parameters from Hanson et al. (2019) by digitising results for fish with
  # smolt age 1-4 and sea winters 1-3 (figs 2 & 3).
  fecunditySlope = 2.9,
  fecundityIntercept = -1.52,
  sexRatio1SW = 0.5,
  sexRatioMSW = 0.7,
  
  # --- environmental scenario ---
  dTwinter = 0, # degrees relative to baseline, whatever baseline is;
  # applied to egg duration
  dgmaxFry = 1, # multiplier on gmaxFW during fry stage (first summer)
  dgmaxParr = 1, # multiplier on gmaxFW during parr stage
  dgmaxOc = 1), # multiplier on marine growth
  ...) { 
  
# --- stage structure ---
stages <- c("egg","fry","parr","smolt","earlyPS",
            "latePS","adultOc","adultCoastal",
            "adultRiver","adultSpawners","eggProduction")
stages_longnames <- c("egg","fry","parr","smolt","early post-smolt",
                      "late post-smolt","adult in ocean","adult on coast",
                      "adult in river","adult spawners","egg production")
nStages <- length(stages)
s <- data.frame( matrix(1:nStages,nrow=1))
colnames(s) <- stages # for more readable code: s$egg <- 1, s$juvSum1 <- 2, etc

# set proportion of spawners female
# 50:50 if 1SW returner; 70:30 if 2SW
# this is used to estimate egg production from returners
if (p$baselineDuration_adultOc > 12) {
propFemale = p$sexRatioMSW
} else {
  propFemale = p$sexRatio1SW
}

# -------------

# check p initialised
#str(p)

blank <- rep(NA,nStages)

# --- key state variables, defined at start of stage ---
res <- data.frame("N" = blank, # number of individuals
          "W" = blank,   # individual weight
          "t0" = blank,  # time (days)
          # --- associated quantities (- )could be back-calculated from the state variables
          # although in practice we do it the other way round)
          "L" = blank,   # length (from weight)
          "dt" = blank,  # stage durations
          "m" = blank    # mortality
          )

# initialise
res$N[1] <- p$N_initial
res$W[1] <- p$W_egg
res$L[1] <- (p$L3overW * p$W_egg) ^ (1/3)
res$t0[1] <- p$yearday_eggDeposition # Nov before "first year"

# -----------------------------
# --- Main loop over stages ---

for (i in 1:(nStages-2) ) {
  
  # --- stage duration ---
  
  if (i==s$egg) {
    # temperature-dependent egg duration
    dt_i_baseline <- p$baselineDuration_egg * (365/12)
    dt_i <- dt_i_baseline / p$Q10_dtegg ^ (p$dTwinter/10)
  } else if (i==s$fry) {
      # fry stage is defined as ending 30 Sep in first year
      dt_i <- p$yearday_endOfFry - res$t0[s$fry]
  } else if (i==s$parr) {
      dt_i <- p$baselineDuration_parr * (365/12)
      if (p$flexibleParrDuration==1) {
        if (res$L[s$parr] < p$ref_length_parr) {
          dt_i <- dt_i + 365
        }
      }
  } else {
    # all other stages
    dt_i <- p[[ paste("baselineDuration_",stages[i],sep="") ]] * (365/12)
  }
  
  # --- growth and weight ---
  
  r_size <- res$W[i] ^ -p$exp_growth # allometry
	g <- 0
	if (i == s$fry) {
	  g <- p$dgmaxFry * p$gmaxFry * r_size
	} else if (i == s$parr) {
	  if (dt_i > 365*2) {
	    gmax <- p$gmaxParr30
	  } else if (dt_i>365) {
	    gmax <- p$gmaxParr18
	  } else {
	    gmax <- p$gmaxParr6
	  }
	  g <- p$dgmaxParr * gmax * r_size
	} else if (i == s$smolt) {
	  g <- p$gmaxParr18 * r_size
	  # matters very little but have to put down something
	} else if (i >= s$earlyPS & i <= s$adultCoastal) {
	    if (p$baselineDuration_adultOc < 12) { # 1SW
	      gmax <- p$gmaxOc1SW
	    } else {
	      gmax <- p$gmaxOc2SW
	    }
	  g <- p$dgmaxOc * gmax * r_size
  }
	Wend_i <- res$W[i] * exp(g * dt_i)
	
	# --- mortality ---
	# density and duration dependence
	if (i == s$egg) {
	  daily_mort <- 1 - (1-p$m_egg) ^ (1/dt_i_baseline)
	  m_i <- 1 - (1-daily_mort) ^ (dt_i)
	} else if (i == s$fry) {
	  m_i <- p$m_fry
	} else if (i == s$parr) {
	  stock <- res$N[i]
	  recruits <- stock * p$parr_ricker_alpha * exp(-p$parr_ricker_beta * stock)
	  if (dt_i > 365*2) {
	    recruits <- recruits * (1 - p$mort_parr_annual) ^ 2
	    # additional penalty for 30 mo parr
	  } else if (dt_i > 365 ) {     	
	    recruits <- recruits * (1 - p$mort_parr_annual)
	    # additional penalty for 18 mo parr
	  }
	  m_i <- 1 - (recruits/stock) # total mortality over stage duration
	} else if (i == s$smolt) {
	    m_i <- p$m_smolt
	} else if (i == s$earlyPS || i == s$latePS) { # size-dependence during post-smolt only
	    refW <- p$ref_length_earlyPS^3 / p$L3overW
	    r_size <- (res$W[i]/refW) ^ p$exp_sizeMort
	    if (p$baselineDuration_adultOc > 12) { # additional mortality for 2SW
	      rlh <- p$rmort2SW
	    } else {
	      rlh <- 1
	    }
	    m_i <- 1 - max(0, 1 - p$m_earlyPS_monthly * r_size * rlh) ^ (dt_i/365*12)
	    if (p$baselineDuration_adultOc > 12) { # additional mortality for 2SW
	      m_i <- m_i * p$rmort2SW
	    }
	} else if (i == s$adultOc) {
	      if (p$baselineDuration_adultOc > 12) { # additional mortality for 2SW
	        rlh <- p$rmort2SW
	      } else {
	        rlh <- 1
	      }
	    m_i <- 1 - max(0, 1 - p$m_adultOc_monthly * rlh) ^ (dt_i/365*12)
	} else if (i == s$adultCoastal) {
	    m_i <- p$m_adultCoastal
	} else if (i == s$adultRiver) {
	    m_i = p$m_adultRiver # no adjustments
	}
	
	i
	m_i

	# --- bookkeepng ---
	res$m[i] <- m_i
	res$dt[i] <- dt_i
	res$N[i+1] <- res$N[i] * (1 - res$m[i])
	res$W[i+1] <- Wend_i
	res$L[i+1] <- (p$L3overW * Wend_i) ^ (1/3)
	res$t0[i+1] <- res$t0[i] + dt_i
 
} 

# --- Calculate egg production ---
# - 'nextGen' renamed to 'adultSpawners' to store survivors of adultRiver.
# - egg production stored in new stage 'eggProduction'.
spawners <- res$N[s$adultSpawners] * propFemale   # number of female spawners
spawner_L_f <- res$L[s$adultSpawners]             # mean spawner size (cm)
fecundity <- 10^( p$fecunditySlope*log10(spawner_L_f) + p$fecundityIntercept )
# eggs per female as function of body length

resultingEggs <- spawners * fecundity       # total egg production
res$m[s$adultSpawners] <- 1.0;               # zero survival of spawners
res$dt[s$adultSpawners] <- 0.0; 

# Note: in future we may want to model varying egg size or weight.
# For now, use same values as for initial eggs.
res$N[s$eggProduction] <- resultingEggs           # update nextGen numbers
res$W[s$eggProduction] <- p$W_egg                 # set W to egg W
res$L[s$eggProduction] <- (p$L3overW * p$W_egg) ^ (1/3)

res$t0[s$eggProduction] <- res$t0[s$adultSpawners] # no time passes in adultSpawner stage
# leave res.dt(s.eggProduction) and res.m(s.eggProduction) as NaN
# ---------------------------------
  

# --- end of main loop over stages ---
# ------------------------------------

res$stages <- stages
res$stages_longnames <- stages_longnames

# check res
#str(res)

# res contains results by stage. Now translate into daily and monthly values ---
# ET edit: should only do this up to adultSpawners stage

# Define daily -- data-frame again
nDays <- length( ceiling(res$t0[1]) : floor(res$t0[s$adultSpawners]) )
daily <- data.frame(matrix(ncol=5,nrow=nDays))
colnames(daily) <- c("t","mort","N","W","L")

daily$t <- ceiling(res$t0[1]) : floor(res$t0[s$adultSpawners])

# mortality shouldn't be interpolated from one stage to the next--we want a flat line
# within each stage
daily_mort <- 1 - (1-res$m[1:s$adultSpawners])^(1/res$dt[1:s$adultSpawners])
for (i in 1:(nStages-2)) {
	ff <- which( daily$t >= res$t0[i] & daily$t < res$t0[i+1])
	daily$mort[ff] <- daily_mort[i]
}
# N can be interpolated: a linear interp of log N matches the idea of constant daily
# loss rates
# Matlab's interp1 func replaced with R's approx
# But approx returns a list with xq and vq and we just want vq
# so need to use $y. Have checked they do the same thing.
daily$N <- exp( approx( res$t0[1:s$adultSpawners],log(res$N[1:s$adultSpawners]),daily$t,method="linear")$y )
# might as well treat W and L the same way
daily$W <- exp( approx(res$t0[1:s$adultSpawners],log(res$W[1:s$adultSpawners]),daily$t,method="linear")$y )
daily$L <- exp(approx(res$t0[1:s$adultSpawners],log(res$L[1:s$adultSpawners]),daily$t,method="linear")$y )

#str(daily)

# monthly output
# R doesn't like end of sequence not being a multiple of increment 365/12
# so have redefined
nMonths <- ceiling(res$t0[s$adultSpawners]/(365/12) )+1
monthly <- data.frame(matrix(ncol=5,nrow=nMonths))
colnames(monthly) <- c("t","mort","N","W","L")
 
monthly$t <- seq( from=ceiling(res$t0[1]), by=365/12, length.out=nMonths )
monthly_mort <- 1 - (1-res$m[1:s$adultSpawners])^(365/12/res$dt[1:s$adultSpawners])
 
for (i in 1:(nStages-1) ) {
   ff <- which(monthly$t >= res$t0[i] & monthly$t < res$t0[i+1])
   monthly$mort[ff] <- monthly_mort[i] 
 }
 
monthly$N <- exp( approx(res$t0[1:s$adultSpawners],log(res$N[1:s$adultSpawners]),monthly$t,method="linear")$y )
monthly$W <- exp( approx(res$t0[1:s$adultSpawners],log(res$W[1:s$adultSpawners]),monthly$t,method="linear")$y )
monthly$L <- exp( approx(res$t0[1:s$adultSpawners],log(res$L[1:s$adultSpawners]),monthly$t,method="linear")$y )

#str(monthly)
outputList <-  list(res,p,daily,monthly,lubridate::now())
names(outputList) <- c("res","p","daily","monthly","timestamp")
return(outputList)

}