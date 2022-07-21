% model experiments to support hand-tuning of mortalityFramework v0.6.
% each experiment has to be incorporated into the default parameter set
% in mortalityFramework.m for the next experiment to give the final result.
%
% Neil Banas july 2022


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


% adult size for 1SW and 2SW returners

clear gmaxOcExpt
gmaxOcExpt.g = 0.02 : 0.005 : 0.06;
for i=1:length(gmaxOcExpt.g)
	res = mortalityFramework('gmaxOc',gmaxOcExpt.g(i), 'baselineDuration_adultOc',4);
	gmaxOcExpt.Ladult1SW(i) = res.L(s.adultRiver);
	res = mortalityFramework('gmaxOc',gmaxOcExpt.g(i), 'baselineDuration_adultOc',4+12);
	gmaxOcExpt.Ladult2SW(i) = res.L(s.adultRiver);
end
figure
plot(gmaxOcExpt.g,gmaxOcExpt.Ladult1SW,'o-',gmaxOcExpt.g,gmaxOcExpt.Ladult2SW,'o-')
xlabel('gmaxOc');
ylabel('L adult');
legend('1SW','2SW');
disp('gmaxOc that gives a 60 cm 1SW adult:')
disp(interp1(gmaxOcExpt.Ladult1SW,gmaxOcExpt.g,60));
disp('gmaxOc that gives a 75 cm 2SW adult:')
disp(interp1(gmaxOcExpt.Ladult2SW,gmaxOcExpt.g,75));



