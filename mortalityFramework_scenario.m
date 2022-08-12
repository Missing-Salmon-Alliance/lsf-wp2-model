function [res,p,monthly,daily] = ...
	mortalityFramework_scenario(smoltAge,seaWinters,varargin);
	
% wrapper for calling mortalityFramework.m in terms of a high-level scenario,
% with the translation into specific parameter values handled internally.
%
% works with v0.7 of mortalityFramework.m.
%
% currently can capture 1,2,3 year old smolts (parr stage of 6, 18, 30 mo)
% 	's1'	's2'	's3'
% and 1,2 sea winters (adult-in-ocean stage of 4 or 16 mo)
%	'1SW'	'2SW'
%
% pass other parameter settings after these, e.g.
%	res = mortalityFramework_scenario('s1','1SW','dgmaxParr',1.2);

paramList = {};

% right now the adjustment to growth rate is handled within mortalityFramework.m:
% we simply pass in a parr duration and the model figures out which parr growth
% rate to use. That complexity could potentially be moved here.
% likewise, we currently do not assume that differences in parr duration grow out of
% differences in fry growth (e.g. s2,s3 fish start out as smaller parr): that could be
% added here.
if strcmpi(smoltAge,'s1') | smoltAge==1
	paramList = cat(2,paramList,'baselineDuration_parr',6,'flexibleParrDuration',0);
elseif strcmpi(smoltAge,'s2') | smoltAge==2
	paramList = cat(2,paramList,'baselineDuration_parr',18,'flexibleParrDuration',0);
elseif strcmpi(smoltAge,'s3') | smoltAge==3
	paramList = cat(2,paramList,'baselineDuration_parr',30,'flexibleParrDuration',0);
else
	warning('unrecognised smolt age');
end

if strcmpi(seaWinters,'1SW') | seaWinters==1
	paramList = cat(2,paramList,'baselineDuration_adultOc',4);
elseif strcmpi(seaWinters,'2SW') | seaWinters==2
	paramList = cat(2,paramList,'baselineDuration_adultOc',16);
else
	warning('unrecognised number of sea winters');
end

if nargin > 2
	paramList = cat(2,paramList,varargin);
end

[res,p,monthly,daily] = mortalityFramework(paramList{:});
