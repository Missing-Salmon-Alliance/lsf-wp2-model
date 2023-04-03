% model experiments to support hand-tuning of mortalityFramework v0.7.
% each experiment has to be incorporated into the default parameter set
% in mortalityFramework.m for the next experiment to give the final result.
%
% Neil Banas aug 2022


[stages,stages_longnames,s] = mortalityFramework('stages');


% translation between dgmaxfry and parr length

clear fryExpt
fryExpt.dg = 0.5 : 0.02 : 1.5;
for i=1:length(fryExpt.dg)
	[res,p] = mortalityFramework('dgmaxFry',fryExpt.dg(i));
	fryExpt.Lparr(i)  = res.L(s.parr);
end
fryExpt.gmaxFry = p.gmaxFry;
figure
plot(fryExpt.dg,fryExpt.Lparr,'o');
grid;
xlabel('dgmaxFry');
ylabel('Lparr');
disp('gmaxFry that gives 6,7,8 cm parr:');
disp(fryExpt.gmaxFry .* interp1(fryExpt.Lparr,fryExpt.dg,[6 7 8]));
gmaxFry_best = fryExpt.gmaxFry .* interp1(fryExpt.Lparr,fryExpt.dg,7);


% parr-length - dgmaxParr combinations that hit the smolt length targets for 6,18,30 mo

clear parrExpt
[parrExpt.gmaxFry, parrExpt.gmaxParr, parrExpt.dur] = ...
	ndgrid(0.0045 : 0.0005 : 0.009,  0.002 : 0.00025 : 0.018,  [6 18 30]);
for i=1:size(parrExpt.gmaxFry,1)
	for j=1:size(parrExpt.gmaxParr,2)
		for k=1:size(parrExpt.dur,3)
			res = mortalityFramework('flexibleParrDuration',0,...
									 'dgmaxFry',1, 'dgmaxParr',1,...
									 'gmaxFry',parrExpt.gmaxFry(i,j,k),...
									 'gmaxParr6',parrExpt.gmaxParr(i,j,k),...
									 'gmaxParr18',parrExpt.gmaxParr(i,j,k),...
									 'gmaxParr30',parrExpt.gmaxParr(i,j,k),...
									 'baselineDuration_parr',parrExpt.dur(i,j,k));
			parrExpt.Lparr(i,j,k) = res.L(s.parr);
			parrExpt.Lsmolt(i,j,k) = res.L(s.smolt);
			parrExpt.Ladult(i,j,k) = res.L(s.adultRiver);
		end
	end
end
ii = find(abs(parrExpt.gmaxFry(:,1)-gmaxFry_best) == ...
		  min(abs(parrExpt.gmaxFry(:,1)-gmaxFry_best)));
figure
plot(squeeze(parrExpt.gmaxParr(ii,:,:)), squeeze(parrExpt.Lsmolt(ii,:,:)), 'o-');
xlabel(['gmaxParr, at gmaxFry = ' num2str(parrExpt.gmaxFry(ii))]);
ylabel('Lsmolt');
legend('6 mo','18 mo','30 mo');
ylim([6 20]);
grid
disp(['at gmaxFry = ' num2str(parrExpt.gmaxFry(ii)) ',']);
disp('gmaxParr that gives 13 cm, 1-year smolts:');
disp(interp1(parrExpt.Lsmolt(ii,:,1),parrExpt.gmaxParr(ii,:,1),13));
disp('gmaxParr that gives 13.5 cm, 2-year smolts:');
disp(interp1(parrExpt.Lsmolt(ii,:,2),parrExpt.gmaxParr(ii,:,2),13.5));
disp('gmaxParr that gives 14 cm, 3-year smolts:');
disp(interp1(parrExpt.Lsmolt(ii,:,3),parrExpt.gmaxParr(ii,:,3),14));
gmaxParr_best = [interp1(parrExpt.Lsmolt(ii,:,1),parrExpt.gmaxParr(ii,:,1),13);
				 interp1(parrExpt.Lsmolt(ii,:,2),parrExpt.gmaxParr(ii,:,2),13.5);
				 interp1(parrExpt.Lsmolt(ii,:,3),parrExpt.gmaxParr(ii,:,3),14)];


% fry - parr growth combinations that span realistic smolt length ranges

gmaxParr_best3 = repmat(reshape(gmaxParr_best,[1 1 3]),size(parrExpt.dur(:,:,1)));
parrExpt.valid = (parrExpt.dur==6 & parrExpt.Lsmolt >= 12 & parrExpt.Lsmolt < 14) ...
			   | (parrExpt.dur==18 & parrExpt.Lsmolt >= 12 & parrExpt.Lsmolt < 14.5) ...
			   | (parrExpt.dur==30 & parrExpt.Lsmolt >= 12 & parrExpt.Lsmolt < 15);
figure
scatter(parrExpt.gmaxFry(parrExpt.valid) ./ gmaxFry_best, ...
		parrExpt.gmaxParr(parrExpt.valid) ./ gmaxParr_best3(parrExpt.valid), ...
		20 + 2*parrExpt.dur(parrExpt.valid), parrExpt.dur(parrExpt.valid), 'o');
colormap([146 54 38; 191 146 54; 34 65 106]./255);
xlabel('dgmaxFry');
ylabel('dgmaxParr');
title('fry-parr growth combinations giving valid smolt sizes')


% adult size for 1SW and 2SW returners.

clear gmaxOcExpt
gmaxOcExpt.g = 0.02 : 0.005 : 0.06;
for i=1:length(gmaxOcExpt.g)
	res = mortalityFramework('gmaxOc1SW',gmaxOcExpt.g(i), 'baselineDuration_adultOc',4);
	gmaxOcExpt.Ladult1SW(i) = res.L(s.adultRiver);
	gmaxOcExpt.survOc1SW(i) = res.N(s.adultRiver) / res.N(s.earlyPS);
	res = mortalityFramework('gmaxOc2SW',gmaxOcExpt.g(i), 'baselineDuration_adultOc',16);
	gmaxOcExpt.Ladult2SW(i) = res.L(s.adultRiver);
	gmaxOcExpt.survOc2SW(i) = res.N(s.adultRiver) / res.N(s.earlyPS);	
end
figure % tuning gmaxOc based on adult size
plot(gmaxOcExpt.g,gmaxOcExpt.Ladult1SW,'o-',gmaxOcExpt.g,gmaxOcExpt.Ladult2SW,'o-')
xlabel('gmaxOc');
ylabel('L adult');
legend('1SW','2SW');
disp('gmaxOc that gives a 60 cm 1SW adult:')
disp(interp1(gmaxOcExpt.Ladult1SW,gmaxOcExpt.g,60));
disp('gmaxOc that gives a 75 cm 2SW adult:')
disp(interp1(gmaxOcExpt.Ladult2SW,gmaxOcExpt.g,75));


% marine mortality

% tuning to three reference numbers:
% - 140 cm smolt, 1SW, marine survival = 4.9%
% - 140 cm smolt, 2SW, marine survival = 0.75%
% - 120 vs 160 cm smolt, 1SW, marine survival varies by 2x
% tuning three parameters:
% - m_earlyPS_monthly
% - rmort2SW
% - exp_sizeMort
clear ocExpt
[ocExpt.m_earlyPS_monthly, ocExpt.rmort2SW, ocExpt.exp_sizeMort] = ...
	ndgrid(0.35 : 0.01 : 0.6,    1.02 : 0.005 : 1.08,    -1.57 : 0.1 : 0);
for i=1:size(ocExpt.rmort2SW,1)
	for j=1:size(ocExpt.rmort2SW,2)
		for k=1:size(ocExpt.rmort2SW,3)
			% basic 1SW and 2SW experiments to tune earlyPS mortality and the 2SW penalty
			res = mortalityFramework_scenario('s2','1SW',...
				'm_earlyPS_monthly',ocExpt.m_earlyPS_monthly(i,j,k),...
				'rmort2SW',ocExpt.rmort2SW(i,j,k),...
				'exp_sizeMort',ocExpt.exp_sizeMort(i,j,k));
			ocExpt.Lsmolt1SW(i,j,k) = res.L(s.smolt);
			ocExpt.Ladult1SW(i,j,k) = res.L(s.adultRiver);
			ocExpt.survOc1SW(i,j,k) = res.N(s.adultRiver) / res.N(s.earlyPS);
			res = mortalityFramework_scenario('s2','2SW',...
				'm_earlyPS_monthly',ocExpt.m_earlyPS_monthly(i,j,k),...
				'rmort2SW',ocExpt.rmort2SW(i,j,k),...
				'exp_sizeMort',ocExpt.exp_sizeMort(i,j,k));
			ocExpt.Lsmolt2SW(i,j,k) = res.L(s.smolt);
			ocExpt.Ladult2SW(i,j,k) = res.L(s.adultRiver);
			ocExpt.survOc2SW(i,j,k) = res.N(s.adultRiver) / res.N(s.earlyPS);
			
			% now turn FW growth up 10% to test how smolt size and
			% mortality are related (and what exp_sizeMort should be)
			Lsmoltlo = ocExpt.Lsmolt1SW(i,j,k);
			survlo = ocExpt.survOc1SW(i,j,k);
			res = mortalityFramework_scenario('s2','1SW',...
				'm_earlyPS_monthly',ocExpt.m_earlyPS_monthly(i,j,k),...
				'rmort2SW',ocExpt.rmort2SW(i,j,k),...
				'exp_sizeMort',ocExpt.exp_sizeMort(i,j,k),...
				'dgmaxFry',1.1, 'dgmaxParr', 1.1);
			Lsmolthi = res.L(s.smolt);
			survhi = res.N(s.adultRiver) / res.N(s.earlyPS);
			ocExpt.Lsmolt_ratio(i,j,k) = Lsmolthi/Lsmoltlo;
			ocExpt.survRatio(i,j,k) = survhi/survlo;
		end
	end
end	
f = find(ocExpt.survOc1SW > 0.049 * 0.9 & ocExpt.survOc1SW < 0.049 * 1.1 & ...
		 ocExpt.survOc2SW > 0.0075 * 0.9 & ocExpt.survOc2SW < 0.0075 * 1.1);
figure
subplot 221
plot(ocExpt.m_earlyPS_monthly(:),ocExpt.survOc1SW(:),'.',...
	 ocExpt.m_earlyPS_monthly(f),ocExpt.survOc1SW(f),'*');
xlabel('m earlyPS monthly');
ylabel('1SW marine survival');
subplot 222
plot(ocExpt.rmort2SW(:),ocExpt.survOc2SW(:),'.',...
	 ocExpt.rmort2SW(f),ocExpt.survOc2SW(f),'*');
xlabel('rmort2SW');
ylabel('2SW marine survival');
subplot 223
plot(ocExpt.exp_sizeMort(:),ocExpt.survRatio(:),'.',...
	 ocExpt.exp_sizeMort(f),ocExpt.survRatio(f),'*');
xlabel('exp sizeMort');
ylabel('survival ratio, 14.7 vs 13.6 cm smolts');
subplot 224
plot(ocExpt.m_earlyPS_monthly(:), ocExpt.exp_sizeMort(:),'.',...
	 ocExpt.m_earlyPS_monthly(f), ocExpt.exp_sizeMort(f),'*');
xlabel('m earlyPS monthly');
ylabel('exp sizeMort');
% at default exp_sizeMort of -1.57, m_earlyPS_monthly = 0.55, rmort2SW = 1.04,
% but smolt size sensitivity is extreme.
% at exp_sizeMort = -0.37, size sensitivity is 1.26x mortality for 8% change in length,
% consistent with 2x mortality for 25% change (12->16cm),
% and m_earlyPS_monthly = 0.40, rmort2SW = 1.07 give a good fit to 1SW, 2SW survival.


% vary FW and marine growth continuously; how does marine survival vs length look
% for each of these parameterisations?
clear LsurvExpt
[LsurvExpt.dgfw, LsurvExpt.dgoc] = meshgrid(0.8 : 0.02 : 1.2, 0.7 : 0.02 : 1.1);
for i=1:size(LsurvExpt.dgfw,1)
	for j=1:size(LsurvExpt.dgfw,2)
		res = mortalityFramework_scenario('s2','1SW',...
			'dgmaxFry',LsurvExpt.dgfw(i,j),'dgmaxParr',LsurvExpt.dgfw(i,j),...
			'dgmaxOc',LsurvExpt.dgoc(i,j));
		LsurvExpt.Lsmolt(i,j) = res.L(s.smolt);
		LsurvExpt.Ladult(i,j) = res.L(s.adultRiver);
		LsurvExpt.survOc(i,j) = res.N(s.adultRiver) / res.N(s.earlyPS);

		res = mortalityFramework_scenario('s2','1SW',...
			'dgmaxFry',LsurvExpt.dgfw(i,j),'dgmaxParr',LsurvExpt.dgfw(i,j),...
			'dgmaxOc',LsurvExpt.dgoc(i,j),...
			'm_earlyPS_monthly',0.40, 'rmort2SW', 1.07, 'exp_sizeMort', -0.37);
		LsurvExpt.Lsmolt_alt(i,j) = res.L(s.smolt);
		LsurvExpt.Ladult_alt(i,j) = res.L(s.adultRiver);
		LsurvExpt.survOc_alt(i,j) = res.N(s.adultRiver) / res.N(s.earlyPS);
	end
end
figure
subplot 221
plot(LsurvExpt.Lsmolt(:),LsurvExpt.survOc(:),'o');
xlabel('smolt length');
ylabel('marine survival, 1SW');
title('(default size dependence)')
subplot 222
plot(LsurvExpt.Ladult(:),LsurvExpt.survOc(:),'o');
xlabel('adult length');
ylabel('marine survival, 1SW');
subplot 223
plot(LsurvExpt.Lsmolt_alt(:),LsurvExpt.survOc_alt(:),'o');
xlabel('smolt length');
ylabel('marine survival, 1SW');
title('(weaker size dependence)')
subplot 224
plot(LsurvExpt.Ladult_alt(:),LsurvExpt.survOc_alt(:),'o');
xlabel('adult length');
ylabel('marine survival, 1SW');
% on the basis of this comparison, the "alternate" mortality parameters 
% (turning size sensitivity way down) have been inserted as the defaults
% in v0.6.1

