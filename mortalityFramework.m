function [res,p,monthly,daily] = mortalityFramework(varargin);

% Salmon Mortality Framework Model
% v0.5.1, July 2022
% Neil Banas, Emma Tyldesley, Colin Bull

% to do -----------------------------
%
% adult duration as a function of size, and affecting marine mortality
%  for the Bush, 6 mo at sea = 0.95 marine mortality,
%               18 mo at sea = 0.98 marine mortality 
%
% actual cycles of prey and temperature in FW, to enable
%     T effects beyond egg
%     season-specific management scenarios


% -----------------------------------


% --- stage structure ---
stages = {'egg','juvSum1','juvLater','smolt','earlyPS','latePS',...
		  'adultOc','adultCoastal','adultRiver','nextGen'};
stages_longnames = ...
	{'egg','juvenile first summer','juvenile after first summer','smolt',...
	'early post-smolt','late post-smolt',...
	'adult in ocean','adult on coast','adult in river','next generation'};
for i=1:length(stages)
	s.(stages{i}) = i; % for more readable code, s.egg = 1, s.juvSum1 = 2, ...
end


% --- model parameters ------------------------------------------------------------------
clear p
p.N_initial = 1e6;
p.baselineDuration_months = [5 6 6 1 3 5 4 2 4];

% --- growth params ---

p.W_egg = 1; % 1 gram
p.Q10_dtegg = 6.5;   % Elliott and Hurley 1998, reciprocal of eq 1a:
				     % 2.12x change over 4 degC
p.exp_growth = 0.31; % Forseth et al. 2001
p.gmaxFry = 0.007; % tuned to give 7 cm at start of parr 
p.gmaxParr6 = 0.013; % tuned to give a 13 cm smolt for 6 mo parr
p.gmaxParr18 = 0.0052; % tuned to give a 13 cm smolt for 18 mo parr
p.gmaxOc = 0.051; % tuned to turn a 13 cm smolt into a 60 cm adult		   	     
	% gmax is c/100 in the Ratkowsky model used by Forseth et al. 2001:
	% growth at 1 g body weight and at optimal temp. T = TM
	% base value is 2.05/100, which is the average over "mod. fast" category,
	% 5 rivers, Jonsson et al. 2001
	% but this needs to be reduced in FW to account for the fact that food is actually 
	% quite seasonal, not year-round
p.ref_length_parr = 7; % smaller than this at start of parr stage,
					   % add 12 mos to parr stage
p.ref_length_earlyPS = 13;
p.ref_length_adultRiver = 60;
p.L3overW = 62^3 / 2500; % length in cm cubed over weight in g
						 % calibrated using Bacon et al. 2009
	     
p.TL_growth = 7.2;   % lower-bound, optimal, and upper-bound temperatures for growth:
p.TM_growth = 18.3;  % Jonsson et al. 2001
p.TU_growth = 24.5;
p.gR_growth = 0.175; % the g parameter in the Ratkowsky model
					 % average over "mod. fast" category, 5 rivers, Jonsson et al. 2001
		
					 
% --- mortality params ---

p.m_egg = 0.1;
p.m_smolt = 0.3; % 0.1 - 0.5
p.m_earlyPS = 0.5;
p.mort_marine_monthly = 0.66; % applied to early/latePS and adult stages.
						% defined at start of earlyPS; declines rapidly with size
						% tuned to make overall marine survival 5%
p.m_adultRiver = 0.09;
p.exp_sizeMort = -1.57; % beta from Ricker 1976 -> Mangel 1994 -> IBASAM
						% (dependence of daily mortality on weight)
p.Lcutoff_sizeMort = 60; % beyond this stops applying size-dependent mortality.
						% makes very little difference whether this is 25 or 60

p.maxParr = 50000;      % juv first summer carrying capacity
p.fryRicker = 0.08;     % hatching to parr Ricker stock-recruit parameter
p.maxSmolts = 27000;    % juv later carrying capacity
p.parrSmoltBH = 0.5;    % parr to smolt Beverton Holt stock-recruit parameter
    % the Ricker and BH parameters above all directly from Salmonmodeller
    % no references given
    % values for R. Bush data were maxParr 650000, fryRicker 0.259
    % BH curve based on 6 mo parr; have to think about how to adjust mortality for
    % 18 mo parr
p.mort_parr_annual = 0.2; % additional mortality if the parr take 18 mo instead of 6


% --- environmental scenario ---

p.dTwinter = 0; % degrees relative to baseline, whatever baseline is;
				% applied to egg duration
p.dgmaxFry = 0.9; % multiplier on gmaxFW during juvSum1.
	% Vary this +/- 20% to get 6-8 cm at start of parr
	% default value is < 1 to tip the base case more clearly into the 18 mo parr scenario
p.dgmaxParr = 1; % multiplier on gmaxFW during juvLater
	% Vary this about 20% to get 11-15 cm smolts
% m_smolt is the 4th control knob 



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
res.t0(1) = datenum('1 Nov 0000') - 366; % Nov before "first year"



% main loop over stages -----------------------------------------------------------------
for i = 1:nStages-1

	% --- stage duration ---------------------------------
	
	dt_i_baseline = p.baselineDuration_months(i) * (365/12);
	if i==s.egg
		% temperature-dependent egg duration
		dt_i = dt_i_baseline ./ p.Q10_dtegg .^ (p.dTwinter./10);
	elseif i==s.juvSum1
		% juvSum1 (as opposed to juvLater) is defined as ending 30 Sep in first year
		dt_i = datenum('30 Sep 0000') - res.t0(s.juvSum1);
	elseif i==s.juvLater
		if res.L(s.juvLater) < p.ref_length_parr
			dt_i = dt_i + 365;
		end
	else
		% all other stages
		dt_i = dt_i_baseline; % placeholder
	end


	% --- growth and weight -------------------------------
	
	Teff = p.TM_growth; % assume optimal temperature for growth.
		% Don't show this to the user at this point--it's too ambiguous.
		% Would the adjustment be around optimal temp or on the rising slope?
	r_size = res.W(i) .^ -p.exp_growth; % allometry
	r_temp = (Teff - p.TL_growth) .* ...
			 (1 - exp(p.gR_growth .* (Teff - p.TU_growth))) ./ ...
			 (p.TM_growth - p.TL_growth) ./ ...
			 (1 - exp(p.gR_growth .* (p.TM_growth - p.TU_growth)));
	r_prey = 1; % haven't included any prey effects
		% could adjust r_prey for fry based on duration, using the idea that growth is 
		% concentrated in a 70 day window (Bacon et al. 2005), considering whether 
		% variation in duration (via dt_egg) falls during that window. Without an 
		% adjustment for this, will underestimate the effect of small variations in
		% dt_egg on growth
	g = 0;
	if i == s.juvSum1
		g = p.dgmaxFry .* p.gmaxFry .* r_size .* r_temp .* r_prey;
	elseif i == s.juvLater
		if dt_i > 365
			gmax = p.gmaxParr18;
		else
			gmax = p.gmaxParr6;
		end
		g = p.dgmaxParr .* gmax .* r_size .* r_temp .* r_prey;
	elseif i == s.smolt
		g = p.gmaxParr18 .* r_size .* r_temp .* r_prey; 
		% matters very little but have to put down something
	elseif i >= s.earlyPS & i <= s.adultCoastal
		g = p.gmaxOc .* r_size .* r_temp .* r_prey;
	end
	Wend_i = res.W(i) .* exp(g .* dt_i);
	
	
	% --- mortality ----------------------------------------
	
	% density and duration dependence
	if i == s.egg
		daily_mort = 1 - (1-p.m_egg)    ^ (1/dt_i_baseline);
		m_i =        1 - (1-daily_mort) ^ (dt_i);
    elseif i == s.juvSum1
        % apply Ricker density-dependent mortality (scramble competition)
        stock = res.N(i);
        recruits = p.fryRicker * stock * exp((-p.fryRicker/(exp(1)*p.maxParr)) * stock);
        m_i = 1 - (recruits/stock); % total mortality over stage duration
		% at the moment this is _not_ adjusted for stage duration, even though
		% fry duration changes in response to egg duration
    elseif i == s.juvLater
        % apply Beverton-Holt density-dependent mortality
        stock = res.N(i);
        recruits = (p.parrSmoltBH * stock) / (1 + (p.parrSmoltBH/p.maxSmolts) * stock);
            %  numbers surviving over stage
        if dt_i > 365
        	recruits = recruits .* (1 - p.mort_parr_annual);
        	% the B-H curve is parameterised for 6 mo parr. This is an additional
        	% penalty for 18 mo parr
        end
        m_i = 1 - (recruits/stock); % total mortality over stage duration
	elseif i == s.smolt
		daily_mort = 1 - (1-p.m_smolt)  ^ (1/dt_i_baseline);
		m_i =        1 - (1-daily_mort) ^ (dt_i);
	elseif i >= s.earlyPS && i <= s.adultCoastal
		refW = p.ref_length_earlyPS^3 / p.L3overW;
		Wcutoff = p.Lcutoff_sizeMort^3 / p.L3overW;
		r_size = (min(res.W(i), Wcutoff)/refW) ^ p.exp_sizeMort;
		m_i = 1 - (1 - p.mort_marine_monthly * r_size) ^ (dt_i/365*12);
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



