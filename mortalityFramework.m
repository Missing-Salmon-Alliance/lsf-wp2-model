function [res,p,monthly,daily] = mortalityFramework(varargin);

% Salmon Mortality Framework Model
% v0.7, Aug 2022
% Neil Banas, Emma Tyldesley, Colin Bull



% --- stage structure ---
stages = {'egg','fry','parr','smolt','earlyPS','latePS',...
		  'adultOc','adultCoastal','adultRiver','nextGen'};
stages_longnames = ...
	{'egg','fry','parr','smolt',...
	'early post-smolt','late post-smolt',...
	'adult in ocean','adult on coast','adult in river','next generation'};
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
p.baselineDuration_parr = 6;
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
p.gmaxOc1SW = 0.051; % tuned to turn a 13-13.5 cm smolt into a 60 cm adult after 1SW
p.gmaxOc2SW = 0.031; % tuned to turn a 13-13.5 cm smolt into a 75 cm adult after 2SW
	% gmax is c/100 in the Ratkowsky model used by Forseth et al. 2001:
	% growth at 1 g body weight and at optimal temp. T = TM
	% base value is 2.05/100, which is the average over "mod. fast" category,
	% 5 rivers, Jonsson et al. 2001
	% but this needs to be reduced in FW to account for the fact that food is actually 
	% quite seasonal, not year-round
p.ref_length_parr = 7; % smaller than this at start of parr stage,
					   % add 12 mos to parr stage, if flexibleParrDuration = 1
					   
p.ref_length_earlyPS = 14; % just for scaling the equations, not tuning targets
p.ref_length_adultRiver = 60;
p.L3overW = 62^3 / 2500; % length in cm cubed over weight in g
						 % calibrated using Bacon et al. 2009

% temperature dependence: not used
p.TL_growth = 7.2;   % lower-bound, optimal, and upper-bound temperatures for growth:
p.TM_growth = 18.3;  % Jonsson et al. 2001
p.TU_growth = 24.5;
p.gR_growth = 0.175; % the g parameter in the Ratkowsky model
					 % average over "mod. fast" category, 5 rivers, Jonsson et al. 2001
		
					 
% --- mortality params ---

p.m_egg = 0.1;
p.maxParr = 50000;      % fry-stage carrying capacity
p.fryRicker = 0.08;     % hatching to parr Ricker stock-recruit parameter
p.maxSmolts = 27000;    % parr-stage carrying capacity
p.parrSmoltBH = 0.5;    % parr to smolt Beverton Holt stock-recruit parameter
    % the Ricker and BH parameters above all directly from Salmonmodeller
    % no references given
    % values for R. Bush data were maxParr 650000, fryRicker 0.259
    % BH curve based on 6 mo parr
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




% --- environmental scenario ---

p.dTwinter = 0; % degrees relative to baseline, whatever baseline is;
				% applied to egg duration
p.dgmaxFry = 1; % multiplier on gmaxFW during fry stage (first summer)
p.dgmaxParr = 1; % multiplier on gmaxFW during parr stage
	% vary these +/- 15% in combination to get a valid range of smolt sizes
p.dgmaxOc = 1; % placeholder; haven't evaluated how it behaves


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
for i = 1:nStages-1

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
	
	Teff = p.TM_growth; % assume optimal temperature for growth
	r_size = res.W(i) .^ -p.exp_growth; % allometry
	r_temp = (Teff - p.TL_growth) .* ...
			 (1 - exp(p.gR_growth .* (Teff - p.TU_growth))) ./ ...
			 (p.TM_growth - p.TL_growth) ./ ...
			 (1 - exp(p.gR_growth .* (p.TM_growth - p.TU_growth)));
	r_prey = 1; % haven't included any prey effects
		% could adjust r_prey for fry based on duration, using the idea that growth is 
		% concentrated in a 70 day window (Bacon et al. 2005), considering whether 
		% variation in duration (via dt_egg) falls during that window
	g = 0;
	if i == s.fry
		g = p.dgmaxFry .* p.gmaxFry .* r_size .* r_temp .* r_prey;
	elseif i == s.parr
		if dt_i > 365*2
			gmax = p.gmaxParr30;
		elseif dt_i > 365
			gmax = p.gmaxParr18;
		else
			gmax = p.gmaxParr6;
		end
		g = p.dgmaxParr .* gmax .* r_size .* r_temp .* r_prey;
	elseif i == s.smolt
		g = p.gmaxParr18 .* r_size .* r_temp .* r_prey; 
		% matters very little but have to put down something
	elseif i >= s.earlyPS & i <= s.adultCoastal
		if p.baselineDuration_adultOc < 12 % 1SW
			gmax = p.gmaxOc1SW;
		else
			gmax = p.gmaxOc2SW;
		end
		g = p.dgmaxOc .* gmax .* r_size .* r_temp .* r_prey;
	end
	Wend_i = res.W(i) .* exp(g .* dt_i);
	
	
	% --- mortality ----------------------------------------
	
	% density and duration dependence
	if i == s.egg
		daily_mort = 1 - (1-p.m_egg)    ^ (1/dt_i_baseline);
		m_i =        1 - (1-daily_mort) ^ (dt_i);
    elseif i == s.fry
        % apply Ricker density-dependent mortality (scramble competition)
        stock = res.N(i);
        recruits = p.fryRicker * stock * exp((-p.fryRicker/(exp(1)*p.maxParr)) * stock);
        m_i = 1 - (recruits/stock); % total mortality over stage duration
		% at the moment this is _not_ adjusted for stage duration, even though
		% fry duration changes in response to egg duration
    elseif i == s.parr
        % apply Beverton-Holt density-dependent mortality
        stock = res.N(i);
        recruits = (p.parrSmoltBH * stock) / (1 + (p.parrSmoltBH/p.maxSmolts) * stock);
            % numbers surviving over stage
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

res.stages = stages';
res.stages_longnames = stages_longnames';



% res contains results by stage. Now translate into daily and monthly values ------------


daily.t = ceil(res.t0(1)) : floor(res.t0(end));
% mortality shouldn't be interpolated from one stage to the next--we want a flat line
% within each stage
daily_mort = 1 - (1-res.m).^(1./res.dt);
for i=1:nStages-1
	ff = find(daily.t >= res.t0(i) & daily.t < res.t0(i+1));
	daily.mort(ff) = daily_mort(i);
end
% N can be interpolated: a linear interp of log N matches the idea of constant daily
% loss rates
daily.N = exp(interp1(res.t0,log(res.N),daily.t));
% might as well treat W and L the same way
daily.W = exp(interp1(res.t0,log(res.W),daily.t));
daily.L = exp(interp1(res.t0,log(res.L),daily.t));


monthly.t = ceil(res.t0(1)) : 365/12 : floor(res.t0(end));
monthly_mort = 1 - (1-res.m).^(365./12./res.dt);
for i=1:nStages-1
	ff = find(monthly.t >= res.t0(i) & monthly.t < res.t0(i+1));
	monthly.mort(ff) = monthly_mort(i);
end
monthly.N = exp(interp1(res.t0,log(res.N),monthly.t));
monthly.W = exp(interp1(res.t0,log(res.W),monthly.t));
monthly.L = exp(interp1(res.t0,log(res.L),monthly.t));



