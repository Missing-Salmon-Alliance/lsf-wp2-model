% sensitivity experiment for mortalityFramework v0.5, July 2022
%
% for a set of parameters representing things that management action can influence,
% and another set representing environmental circumstances,
%
% - effect on # of returners of a standardised amount of change in each param
% - how much change in each is required to produce a 10% change in returners?
% - what range on parameters is safe for a user who is varying them in all combinations? 


params = {	'm_egg',...
			'dgmaxFry','maxParr','fryRicker',...
			'dgmaxParr','maxSmolts','parrSmoltBH',...
			'm_smolt',...
			'dTwinter','gmaxOc','m_earlyPS','mort_marine_monthly'};

[res0,p0] = mortalityFramework; % p0 = baseline parameters

dp = 0.2; % +/- 20% change unless otherwise specified
clear ranges
for k=1:length(params)
	ranges(k,:) = p0.(params{k}) .* [1-dp 1+dp];
end
% can put overrides here
ranges(strmatch('m_smolt',params),:)= [0.1 0.5];
ranges(strmatch('dTwinter',params),:) = [-2 2];


stages = {'egg','juvSum1','juvLater','smolt','earlyPS','latePS',...
		  'adultOc','adultCoastal','adultRiver','nextGen'};
for i=1:length(stages)
	s.(stages{i}) = i; % for more readable code, s.egg = 1, s.juvSum1 = 2, ...
end



% one-at-a-time sensitivity
clear out_single
for k=1:length(params)
	p=p0;
	p.(params{k}) = ranges(k,1); % param k, low
	res = mortalityFramework(p);
	out_single(k,1) = max(0, res.N(s.adultRiver));
	p=p0;
	p.(params{k}) = ranges(k,2); % param k, high
	res = mortalityFramework(p);
	out_single(k,2) = max(0, res.N(s.adultRiver));
end
figure
out0 = res0.N(s.adultRiver);
plot(1:length(params),out_single(:,1),'bo',1:length(params),out_single(:,2),'ro');
hold on
plot(xlim,out0.*[1 1],'b:');
set(gca,'xtick',1:length(params),'xticklabel',...
	params,'xticklabelrotation',90,'ticklabelinterpreter','none');
ylabel('Number of returning adults');



% all-at-once sensitivity
N = 1e4;
r0 = rand(length(params),N);
r0 = 0.5 + (r0-0.5) .* repmat((1:N)/N,[length(params) 1]);
vals = repmat(ranges(:,1),[1 N]) + ...
	   r0 .* repmat(ranges(:,2)-ranges(:,1),[1 N]);
clear out
out.N = nan.*ones(1,N);
out.Ladult = nan.*ones(1,N);
out.Lsmolt = nan.*ones(1,N);
out.survFW = nan.*ones(1,N);
out.survOc = nan.*ones(1,N);
for i=1:N
	p = p0;
	for k=1:length(params)
		p.(params{k}) = vals(k,i);
	end
	res = mortalityFramework(p);
	out.N(i) = max(0,res.N(s.adultRiver));
	out.survFW(i) = res.N(s.earlyPS)./res.N(s.egg);
	out.survOc(i) = res.N(s.adultRiver)./res.N(s.earlyPS);
	out.Ladult(i) = res.L(s.adultRiver);
	out.Lsmolt(i) = res.L(s.smolt);
end

% how much variation is too much?
r=r0;
r(strmatch('m_smolt',params)) = 0; % omit these hand-selected ranges from the analysis
r(strmatch('dTwinter',params)) = 0;
figure
subplot 221
plot(max(abs(r-0.5)),out.survFW,'o')
title('FW survival');
subplot 222
plot(max(abs(r-0.5)),out.survOc,'o')
title('Marine survival');
subplot 223
plot(max(abs(r-0.5)),out.Ladult,'o')
title('Ladult')
subplot 224
plot(max(abs(r-0.5)),out.Lsmolt,'o')
title('Lsmolt')




% nonlinear response to dgmaxFry
dg = linspace(0.7,1.3,20);
clear N Lsmolt
for i=1:length(dg)
	res = mortalityFramework('dgmaxFry',dg(i));
	N(i) = res.N(s.adultRiver);
	Lsmolt(i) = res.L(s.smolt);
end
figure
subplot 221
plot(dg,N,'o-');
xlabel('dgmaxFry');
ylabel('adults in river');
subplot 222
plot(dg,Lsmolt,'o-');
xlabel('dgmaxFry');
ylabel('smolt length');
subplot 223
plot(Lsmolt,N,'o')
xlabel('smolt length');
ylabel('adults in river');

