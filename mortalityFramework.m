function [res,p,monthly,daily] = mortalityFramework(varargin);

% Salmon Mortality Framework Model
% v0.8, Jan 2023
% Neil Banas, Emma Tyldesley, Colin Bull

% replacing fry Ricker curve with a density-independent mortality
% replacing parr Beverton-Holt with Ricker


% --- stage structure ---
stages = {'egg','fry','parr','smolt','earlyPS','latePS',...
		  'adultOc','adultCoastal','adultRiver','adultSpawners','eggProduction'};
stages_longnames = ...
	{'egg','fry','parr','smolt',...
	'early post-smolt','late post-smolt',...
	'adult in ocean','adult on coast','adult in river','adult spawners','egg production'};
for i=1:length(stages)
	s.(stages{i}) = i; % for more readable code, s.egg = 1, s.fry = 2, ...
end

% you can get these variables by calling
% [stages,stages_longnames,s] = mortalityFramework('stages');
if nargin>0 && strcmpi(varargin{1},'stages')
	res = stages;
	p = stages_longnames;
	monthly = s;
	return;
end	

% --- model parameters ------------------------------------------------------------------
clear p
p.N_initial = 1e6;

% --- life history schedule ---

p.yearday_eggDeposition = datenum('1 Nov 0000') - 366;
p.baselineDuration_egg = 5; % in months
p.yearday_endOfFry = datenum('30 Sep 0000');
p.baselineDuration_parr = 18;
p.flexibleParrDuration = 0; % if 1, choose parr duration based on length at end
							% of fry stage, as opposed to this being set by the
							% user via parameter values
p.baselineDuration_smolt = 1;
p.baselineDuration_earlyPS = 3;
p.baselineDuration_latePS = 5;
p.baselineDuration_adultOc = 4;
p.baselineDuration_adultCoastal = 2;
p.baselineDuration_adultRiver = 4;

% --- growth params ---

p.W_egg = 1; % 1 gram
p.Q10_dtegg = 6.5;   % Elliott and Hurley 1998, reciprocal of eq 1a:
				     % 2.12x change over 4 degC
p.exp_growth = 0.31; % Forseth et al. 2001
p.gmaxFry = 0.007; % tuned to give 7 cm at start of parr
				   % vary this about +/- 30% to give 6-8 cm parr
p.gmaxParr6 = 0.0152; % tuned to give a 13 cm smolt for 6 mo parr
p.gmaxParr18 = 0.0054; % tuned to give a 13.5 cm smolt for 18 mo parr
p.gmaxParr30 = 0.0034; % tuned to give a 14 cm smolt for 30 mo parr
	% for comparison, c/100 in the Ratkowsky model used by Forseth et al. 2001
	% (growth at 1 g body weight and at optimal temp.) for the "mod. fast" category
	% would give 0.0205. These tuned gmax values should always be less than this,
	% to reflect the fact that food is actually quite seasonal, not year-round, and
	% that temperature can't ever be better than optimal
p.gmaxOc1SW = 0.051; % tuned to turn a 13-13.5 cm smolt into a 60 cm adult after 1SW
p.gmaxOc2SW = 0.031; % tuned to turn a 13-13.5 cm smolt into a 75 cm adult after 2SW
p.ref_length_parr = 7; % if flexibleParrDuration = 1, then if smaller than this at
					   % start of parr stage, add 12 mos to parr stage
					   
p.ref_length_earlyPS = 14; % just for scaling the equations, not tuning targets
p.ref_length_adultRiver = 60;
p.L3overW = 62^3 / 2500; % length in cm cubed over weight in g
						 % calibrated using Bacon et al. 2009
		
					 
% --- mortality params ---

p.m_egg = 0.1;
p.m_fry = 0.97;			% based (approximately) on the default in SalmonModeller
p.parr_ricker_slope = 0.33; % placeholder
p.parr_ricker_fryatmax = 50000; % placeholder
	% parr = fry * slope * exp(-fry/fryatmax)
p.mort_parr_annual = 0.2; % additional mortality if the parr take 18 mo instead of 6
p.m_smolt = 0.3; % 0.1 - 0.5
p.m_earlyPS_monthly = 0.40; % at ref_length_earlyPS; declines rapidly with size
p.exp_sizeMort = -0.37; % dependence of daily mortality on weight
p.rmort2SW = 1.07; % additional marine mortality (multiplier) for 2SW vs 1SW
p.m_adultOc_monthly = 0.03;
p.m_adultCoastal = 0.1;
p.m_adultRiver = 0.09;

% marine mortality parameters tuned based on 1SW, 2SW survival for the Bush,
% and so that a 25% change in smolt length (12-16 cm) has roughly a 2x effect
% on marine survival. Alternately, if we keep the exp_sizeMort consistent with
% Ricker 1976 -> Mangel 1994 -> IBASAM, and retune, we would have
% p.exp_sizeMort = -1.57;
% p.m_earlyPS_monthly = 0.55;
% p.rmort2SW = 1.04;
% wth same base-case marine survival but with really extreme variation as smolt length
% changes.


% --- fecundity params ---

% Fecundity estimated as function of fork length (L_f) in cm:
% log10(eggs)=m.log10(L_f)+c
% Parameters from Hanson et al. (2019) by digitising results for fish with
% smolt age 1-4 and sea winters 1-3 (figs 2 & 3).
p.fecunditySlope = 2.9;
p.fecundityIntercept = -1.52;
p.sexRatio1SW = 0.5; % proportion female
p.sexRatioMSW = 0.7;

% proportion of spawners female is 50:50 if 1SW returner; 70:30 if 2SW
% this is used to estimate egg production from returners
if p.baselineDuration_adultOc > 12
    propFemale = p.sexRatioMSW;
else 
    propFemale = p.sexRatio1SW;
end


% --- environmental scenario ---

p.dTwinter = 0; % degrees relative to baseline, whatever baseline is;
				% applied to egg duration
p.dgmaxFry = 1; % multiplier on gmaxFW during fry stage (first summer)
p.dgmaxParr = 1; % multiplier on gmaxFW during parr stage
	% vary these +/- 15% in combination to get a valid range of smolt sizes
p.dgmaxOc = 1; % multiplier on marine growth


% override defaults based on function inputs
if nargin==1 && isstruct(varargin{1})
	% input = parameter structure
	p = varargin{1};
else
	% input = name-value pairs
	for k=1:2:nargin
		if isfield(p,varargin{k})
			p.(varargin{k}) = varargin{k+1};
		else
			disp(['don''t recognise ''' varargin{k} '''; ignoring'])
		end
	end
end

nStages = length(stages);
blank = repmat(nan,[nStages 1]);
clear res; % structure to hold results
% --- key state variables, defined at start of stage ---
res.N = blank; % number of individuals
res.W = blank; % individual weight 
res.t0 = blank; % time (days)
% --- associated quantities (could be back-calculated from the state variables, although
%     in practice we do it the other way round)
res.L = blank; % length (from weight)
res.dt = blank; % stage duration
res.m = blank; % mortality

% initialise
res.N(1) = p.N_initial;
res.W(1) = p.W_egg;
res.L(1) = (p.L3overW .* p.W_egg) .^ (1/3);
res.t0(1) = p.yearday_eggDeposition;

% main loop over stages -----------------------------------------------------------------
for i = 1:nStages-2 % nStages-1 -> nStages-2

	% --- stage duration ---------------------------------
	
	if i==s.egg
		% temperature-dependent egg duration
		dt_i_baseline = p.baselineDuration_egg * (365/12);
		dt_i = dt_i_baseline ./ p.Q10_dtegg .^ (p.dTwinter./10);
	elseif i==s.fry
		% fry stage is defined as ending 30 Sep in first year
		dt_i = datenum('30 Sep 0000') - res.t0(s.fry);
	elseif i==s.parr
		dt_i = p.baselineDuration_parr * (365/12);
		if p.flexibleParrDuration
			if res.L(s.parr) < p.ref_length_parr
				dt_i = dt_i + 365;
			end
		end
	else % all other stages
		dt_i = p.(['baselineDuration_' stages{i}]) * (365/12);
	end


	% --- growth and weight -------------------------------
	
	r_size = res.W(i) .^ -p.exp_growth; % allometry
	g = 0;
	if i == s.fry
		g = p.dgmaxFry .* p.gmaxFry .* r_size;
	elseif i == s.parr
		if dt_i > 365*2
			gmax = p.gmaxParr30;
		elseif dt_i > 365
			gmax = p.gmaxParr18;
		else
			gmax = p.gmaxParr6;
		end
		g = p.dgmaxParr .* gmax .* r_size;
	elseif i == s.smolt
		g = p.dgmaxParr .* p.gmaxParr18 .* r_size; 
		% matters very little but have to put down something
	elseif i >= s.earlyPS & i <= s.adultCoastal
		if p.baselineDuration_adultOc < 12 % 1SW
			gmax = p.gmaxOc1SW;
		else
			gmax = p.gmaxOc2SW;
		end
		g = p.dgmaxOc .* gmax .* r_size;
	end
	Wend_i = res.W(i) .* exp(g .* dt_i);
	
	
	% --- mortality ----------------------------------------
	
	% density and duration dependence
	if i == s.egg
		daily_mort = 1 - (1-p.m_egg)    ^ (1/dt_i_baseline);
		m_i =        1 - (1-daily_mort) ^ (dt_i);
    elseif i == s.fry
		daily_mort = 1 - (1-p.m_fry)    ^ (1/dt_i_baseline);
		m_i =        1 - (1-daily_mort) ^ (dt_i);
    elseif i == s.parr
        stock = res.N(i);
        recruits = stock * p.parr_ricker_slope * exp(-stock/p.parr_ricker_fryatmax);
        if dt_i > 365*2
        	recruits = recruits .* (1 - p.mort_parr_annual) ^ 2;
        	% additional penalty for 30 mo parr
		elseif dt_i > 365        	
        	recruits = recruits .* (1 - p.mort_parr_annual);
        	% additional penalty for 18 mo parr
        end
        m_i = 1 - (recruits/stock); % total mortality over stage duration
	elseif i == s.smolt
		m_i =  p.m_smolt;
	elseif i == s.earlyPS || i == s.latePS % size-dependence during post-smolt only
		refW = p.ref_length_earlyPS^3 / p.L3overW;
		r_size = (res.W(i)/refW) ^ p.exp_sizeMort;
		if p.baselineDuration_adultOc > 12 % additional mortality for 2SW
			rlh = p.rmort2SW;
		else
			rlh = 1;
		end
		m_i = 1 - max(0, 1 - p.m_earlyPS_monthly * r_size * rlh) ^ (dt_i/365*12);
		if p.baselineDuration_adultOc > 12 % additional mortality for 2SW
			m_i = m_i * p.rmort2SW;
		end
	elseif i == s.adultOc
		if p.baselineDuration_adultOc > 12 % additional mortality for 2SW
			rlh = p.rmort2SW;
		else
			rlh = 1;
		end
		m_i = 1 - max(0, 1 - p.m_adultOc_monthly * rlh) ^ (dt_i/365*12);
	elseif i == s.adultCoastal
		m_i = p.m_adultCoastal;
	elseif i == s.adultRiver
		m_i = p.m_adultRiver; % no adjustments
	end

	
	% --- bookkeepng ---
	res.m(i) = m_i;
	res.dt(i) = dt_i;
	res.N(i+1) = res.N(i) * (1 - res.m(i));
	res.W(i+1) = Wend_i;
	res.L(i+1) = (p.L3overW .* Wend_i) .^ (1/3);
	res.t0(i+1) = res.t0(i) + dt_i;
    
end

% -- ET edit --
%
% --- Calculate egg production ---
% - 'nextGen' renamed to 'adultSpawners' to store survivors of adultRiver.
% - egg production stored in new stage 'eggProduction'.

spawners = res.N(s.adultSpawners) * propFemale;   % number of female spawners
spawner_L_f = res.L(s.adultSpawners);             % mean spawner size (cm)
fecundity = 10^( p.fecunditySlope*log10(spawner_L_f) + p.fecundityIntercept );
                                            % eggs per female as function of body length
resultingEggs = spawners * fecundity;       % total egg production
res.m(s.adultSpawners) = 1.0;               % zero survival of spawners
res.dt(s.adultSpawners) = 0.0;              

% Note: in future we may want to model varying egg size or weight.
% For now, use same values as for initial eggs.
res.N(s.eggProduction) = resultingEggs;     % update egg production numbers
res.W(s.eggProduction) = p.W_egg;           % set W to egg W
res.L(s.eggProduction) = (p.L3overW .* p.W_egg) .^ (1/3);
                                            % set L to egg L
res.t0(s.eggProduction) = res.t0(s.adultSpawners); % no time passes in adultSpawner stage
% leave res.dt(s.eggProduction) and res.m(s.eggProduction) as NaN

% ---------------------------------

res.stages = stages';
res.stages_longnames = stages_longnames';

% res contains results by stage. Now translate into daily and monthly values ------------

% ET edit: should only do this up to adultSpawners stage
daily.t = ceil(res.t0(1)) : floor(res.t0(s.adultSpawners));

% mortality shouldn't be interpolated from one stage to the next--we want a flat line
% within each stage
daily_mort = 1 - (1-res.m(1:s.adultSpawners)).^(1./res.dt(1:s.adultSpawners));
for i=1:nStages-2 % only do this up to adultSpawners stage
	ff = find(daily.t >= res.t0(i) & daily.t < res.t0(i+1));
	daily.mort(ff) = daily_mort(i);
end
% N can be interpolated: a linear interp of log N matches the idea of constant daily
% loss rates
daily.N = exp(interp1(res.t0(1:s.adultSpawners),log(res.N(1:s.adultSpawners)),daily.t));
% might as well treat W and L the same way
daily.W = exp(interp1(res.t0(1:s.adultSpawners),log(res.W(1:s.adultSpawners)),daily.t));
daily.L = exp(interp1(res.t0(1:s.adultSpawners),log(res.L(1:s.adultSpawners)),daily.t));

monthly.t = ceil(res.t0(1)) : 365/12 : floor(res.t0(s.adultSpawners));
monthly_mort = 1 - (1-res.m(1:s.adultSpawners)).^(365./12./res.dt(1:s.adultSpawners));

for i=1:nStages-2
	ff = find(monthly.t >= res.t0(i) & monthly.t < res.t0(i+1));
	monthly.mort(ff) = monthly_mort(i);
end
monthly.N = exp(interp1(res.t0(1:s.adultSpawners),log(res.N(1:s.adultSpawners)),monthly.t));
monthly.W = exp(interp1(res.t0(1:s.adultSpawners),log(res.W(1:s.adultSpawners)),monthly.t));
monthly.L = exp(interp1(res.t0(1:s.adultSpawners),log(res.L(1:s.adultSpawners)),monthly.t));



