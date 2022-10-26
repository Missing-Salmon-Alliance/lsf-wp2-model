%% plot mortalityframework output

nStages = length(res.stages);

close all
n.mos = length(monthly.t)
mos = datetime(0000,1,1)+calmonths((1:n.mos)-3);
mos.Format='MMM';
figure(1), clf
subplot(4,1,1)
semilogy(mos,monthly.N,'linewidth',1.5)
ylabel('number of fish')
grid on
xlim([min(mos) max(mos)+calmonths(1)])
xticks(mos)
xticklabels(cellstr(mos))
title('Number of fish')
% add number of eggs
hold on
plot(mos(end),res.N(end),'ro')
% add stage labels
subplot(4,1,2)
stairs(mos,monthly.mort,'linewidth',1.5)
title('Mortality')
ylabel('mortality rate (per month)')
grid on
xlim([min(mos) max(mos)+calmonths(1)])
xticks(mos)
xticklabels(cellstr(mos))
subplot(4,1,3)
plot(mos,monthly.W,'linewidth',1.5)
title('Weight')
ylabel('weight of individual (g)')
grid on
xlim([min(mos) max(mos)+calmonths(1)])
xlim([min(mos) max(mos)+calmonths(1)])
xticks(mos)
xticklabels(cellstr(mos))
subplot(4,1,4)
plot(mos,monthly.L,'linewidth',1.5)
title('Length')
ylabel('length of individual (cm)')
grid on
xlim([min(mos) max(mos)+calmonths(1)])
xticks(mos)
xticklabels(cellstr(mos))
xlabel('Time')
sgtitle('mortality framework output by month')
%%
figure(2), clf
semilogy(res.t0(1:(end-1)),res.N(1:(end-1)),'-','linewidth',1.5)
xlabel('time (days relative to 1 Jan in hatching year)')
ylabel('number of fish')
grid on
hold on
semilogy(res.t0(end),res.N(end),'ro')

for i=1:nStages-1
   semilogy([res.t0(i) res.t0(i)],[100 max(res.N)],'k')
   text(res.t0(i)+10,110,res.stages_longnames(i),'Rotation',90,'FontSize',12)
end
text(res.t0(end)+10,res.N(end),res.stages_longnames(end),'FontSize',12)
title('mortality framework native output')