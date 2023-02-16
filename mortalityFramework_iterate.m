function [res,p] = mortalityFramework_iterate(Ngen,varargin);

% res = mortalityFramework_iterate(Ngen,...);
%
%    iterates the mortality framework model over _Ngen_ generations.
%    additional arguments are scenario/model parameters. Values can be given as constants 
%    or as vectors of length Ngen.
%
%    for example
%       gen = mortalityFramework_iterate(5,2,1,'dgmaxfry',1.1);
%       runs 7 generations of 2yo smolts, 1SW, rapid fry growth
%

[stages,stages_longnames,s] = mortalityFramework('stages');

% figure out what N_initial to use for the first generation. Since this might be set
% by the defaults, the scenario, or explicitly as additional parameters, running the
% first generation an extra time is the safest way to do this
args = varargin;
for k = 2:2:length(args)
	if length(args{k})==Ngen, args{k} = args{k}(1); end
end
[~,p] = mortalityFramework_scenario(args{:});
Negg = p.N_initial;

% set up output structure
clear res
res.stages = stages;
res.stages_longnames = stages_longnames;
fields = {'N','W','t0','L','dt','m'}; % output fields to save

% iterate
for n = 1:Ngen
	% param values for generation n
	args = varargin;
	for k = 2:2:length(args)
		if length(args{k})==Ngen
			args{k} = args{k}(n);
		end
	end
	resn = mortalityFramework_scenario(args{:},'N_initial',Negg);
	% save output
	for k = 1:length(fields)
		res.(fields{k})(:,n) = resn.(fields{k});
	end
	% set up next generation
	Negg = res.N(s.eggProduction,n);
end

% fix timebase
for n = 2:Ngen
	res.t0(:,n) = res.t0(:,n) + res.t0(end,n-1);
end
% in the saved output, there's an offset of 60 d between end of generation n and
% start of generation n+1 that needs cleaned up

