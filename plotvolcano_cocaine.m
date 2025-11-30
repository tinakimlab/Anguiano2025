load('cocaine_hits_filtered.mat');

% "cocaine_hits" is the final filtered cocaine hits, with column1=log2fc,
% and column2=P-value

Qthresh=0.727;

figure;
plot(allproteins(:,1),-log10(allproteins(:,2)),'.k');
hold on;
plot(cocaine_hits(:,1),-log10(cocaine_hits(:,2)),'.r');
hold on
plot([-6.5 6.5],[-log10(Qthresh) -log10(Qthresh)],'--k');
plot([0.585 0.585],[0 8],'--k')
plot([-0.585 -0.585],[0 8],'--k')
ylim([0 8])
xlim([-6.5 6.5])