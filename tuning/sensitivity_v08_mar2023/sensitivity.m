[stages,stages_longnames,s] = mortalityFramework('stages');
Ni = linspace(0,5e6,30);


% read list of experiments
[~,~,raw] = xlsread('sensitivityExperiments.xlsx');
clear expt
for k=1:size(raw,1)-1 % one header row
	expt.name{k} = raw{k+1,1};
	p = {raw{k+1,2:end}};
	pairs = length(find(~isnan([p{2:2:end}])));
	expt.params{k} = {p{1:pairs*2}};
end

baseCase = strmatch('2010s-base',expt.name);

for k=1:length(expt.name)
	expt.curve{k}.egg = Ni;
	disp(['experiment ' num2str(k) ' ------']);
	disp(strvcat(expt.params{k}{:}));
	for i=1:length(Ni)
		[res,p] = mortalityFramework_scenario(2,1,...
			expt.params{baseCase}{:}, expt.params{k}{:}, 'N_initial',Ni(i));
		expt.curve{k}.smolt(i) = res.N(s.earlyPS);
	end
	curve = expt.curve{k};
	expt.mrep(k) = res.N(s.earlyPS) / res.N(s.eggProduction); % constant across Ni cases
	expt.smolt(k) = max(curve.smolt);
	expt.eggatmax(k) = curve.egg(curve.smolt==max(curve.smolt));
	expt.smoltreq(k) = expt.mrep(k) .* expt.eggatmax(k);
end


figure
cmap = pairedCatColours;
cmap = cat(1,cmap(1,:),[0 0 0],cmap(3:end,:),cmap(2,:));
for k=1:length(expt.name)
 	plot(expt.curve{k}.egg, expt.curve{k}.smolt + 50*k,...
 		'k-','color',cmap(k,:),'linewidth',1);
	hold on
	yl = ylim;
	plot(xlim,xlim.*expt.mrep(k) + 50*k,'k--','color',cmap(k,:),'linewidth',1);
	ylim([0 5e4]);
	plot(expt.eggatmax(k),expt.smolt(k),'k*','markeredgecolor',cmap(k,:));
	plot(expt.eggatmax(k),expt.smoltreq(k),'k^','markeredgecolor',cmap(k,:));
	grid;
	xlabel('Eggs');
	ylabel('Smolts');
end



figure
subplot 311
for k=1:length(expt.name)
	h = bar(k, expt.eggatmax(k));
	set(h,'FaceColor',cmap(k,:));
	hold on
end
xlim([0 length(expt.name)+1]);
set(gca,'xtick',1:length(expt.name),'xticklabel',expt.name,'xticklabelrotation',-45);
ylabel('Eggs at max smolt production');

subplot 312
for k=1:length(expt.name)
	h = bar(k, expt.smolt(k));
	set(h,'FaceColor',cmap(k,:));
	hold on
end
set(gca,'xtick',1:length(expt.name),'xticklabel',expt.name,'xticklabelrotation',-45);
ylabel('Max smolt production');

subplot 313
for k=1:length(expt.name)
	h = bar(k, expt.smolt(k) - expt.smoltreq(k));
	set(h,'FaceColor',cmap(k,:));
	hold on
end
set(gca,'xtick',1:length(expt.name),'xticklabel',expt.name,'xticklabelrotation',-45);
ylabel('Surplus smolt production');
