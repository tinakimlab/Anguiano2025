% Plot Nrxn-TurboID NAc or PAG volcano plots for scRNAseq datasets
% Note: proteins have been filtered by those with >=2 peptides, & "_MOUSE"
% proteome
% First column = log2FC. Second column = adj. pval.
close all
clear all
load('projections_turbo.mat');

%% PLOT NAc

thresh=0.05;

%Remove NAc proteins that were not found in NAc RNA dataset (NA values)
[NAc, NAc_delidx]=rmmissing(NAc);
NAc_sizes=NAcproteins_foldchange;
NAc_sizes(NAc_delidx)=[];


% Plot datapoints with significance below cutoff on graph (adj. pvalue insignificant)
NAc_Pbelow = NAc(NAc(:,2)>=thresh,:);
NAc_Pbelow_idx=find(NAc(:,2)>=thresh);

% Create graph
figure;
hold on

for a=1:length(NAc_Pbelow)
    bubblechart(NAc_Pbelow(a,1),-log10(NAc_Pbelow(a,2)),NAc_sizes(NAc_Pbelow_idx(a)),'MarkerFaceColor','k','MarkerEdgeColor','none')
end

bubblelim([min(NAc_sizes) max(NAc_sizes)]);
bubblesize([4 12]);

% Plot datapoints with adj. p-value above threshold on graph (adj. pvalue < 0.05)
NAc_pabove = NAc(NAc(:,2)<thresh,:);
NAc_Pabove_idx=find(NAc(:,2)<thresh);


for a=1:length(NAc_pabove)
    if abs(NAc_pabove(a,1))>=0.585
        bubblechart(NAc_pabove(a,1),-log10(NAc_pabove(a,2)),NAc_sizes(NAc_Pabove_idx(a)),'MarkerFaceColor','red','MarkerEdgeColor','none')
    else
        bubblechart(NAc_pabove(a,1),-log10(NAc_pabove(a,2)),NAc_sizes(NAc_Pabove_idx(a)),'MarkerFaceColor','k','MarkerEdgeColor','none')
    end
end

plot([-10 10],[-log10(thresh) -log10(thresh)],'--k');
plot([-0.585 -0.585],[0 18],'--k');
plot([0.585 0.585],[0 18],'--k');
ylim([-0.5 18])
xlim([-4 4])


NAcproteins_NAcRNAhits=length(find(NAc_pabove(:,1)>=0.585));
NAcproteins_PAGRNAhits=length(find(NAc_pabove(:,1)<=-0.585));


%%

figure;
hold on;

%Remove PAG proteins that were not found in PAG RNA dataset (NA values)
[PAG, PAG_delidx]=rmmissing(PAG);
PAG_sizes=PAGproteins_foldchange;
PAG_sizes(PAG_delidx)=[];


% Plot datapoints with significance below cutoff on graph (adj. pvalue insignificant)
PAG_Pbelow = PAG(PAG(:,2)>=thresh,:);
PAG_Pbelow_idx=find(PAG(:,2)>=thresh);

for a=1:length(PAG_Pbelow)
    bubblechart(PAG_Pbelow(a,1),-log10(PAG_Pbelow(a,2)),PAG_sizes(PAG_Pbelow_idx(a)),'MarkerFaceColor','k','MarkerEdgeColor','none')
end

bubblelim([min(PAG_sizes) max(PAG_sizes)]);
bubblesize([4 12]);

% Plot datapoints with ad. p-value above threshold on graph (adj. pvalue < 0.05)
PAG_pabove = PAG(PAG(:,2)<thresh,:);
PAG_Pabove_idx=find(PAG(:,2)<thresh);

for a=1:length(PAG_pabove)
    if abs(PAG_pabove(a,1))>=0.585
        bubblechart(PAG_pabove(a,1),-log10(PAG_pabove(a,2)),PAG_sizes(PAG_Pabove_idx(a)),'MarkerFaceColor','g','MarkerEdgeColor','none')
    else
        bubblechart(PAG_pabove(a,1),-log10(PAG_pabove(a,2)),PAG_sizes(PAG_Pabove_idx(a)),'MarkerFaceColor','k','MarkerEdgeColor','none')
    end
end


plot([-10 10],[-log10(thresh) -log10(thresh)],'--k');
plot([-0.585 -0.585],[0 18],'--k');
plot([0.585 0.585],[0 18],'--k');
ylim([-0.5 18])
xlim([-4 4])


PAGproteins_NAcRNAhits=length(find(PAG_pabove(:,1)>=0.585));
PAGproteins_PAGRNAhits=length(find(PAG_pabove(:,1)<=-0.585));


%%
%NAc_above = NAc_pabove(abs(NAc_pabove(:,1))>=0.585,1);
%PAG_above = PAG_pabove(abs(PAG_pabove(:,1))>=0.585,1);

% plot only RNA NAc/PAG fold changes that are significant
NAc_above = NAc_pabove(:,1);
PAG_above = PAG_pabove(:,1);

%NAc_above=NAc(:,1);
%PAG_above=PAG(:,1);


figure; histogram(NAc_above(:,1),'normalization','pdf','binwidth',0.2925)
hold on
histogram(PAG_above(:,1),'normalization','pdf','binwidth',0.2925)
plot([median(NAc_above) median(NAc_above)],[0 1.5],'-k')
plot([median(PAG_above) median(PAG_above)],[0 1.5],'-k')
xlim([-4 4])

[p,h]=ranksum(NAc_above,PAG_above);


%% Plot piecharts

%plot matched, mismatched, and below foldchange threshold

NAc_match=length(find(NAc_pabove(:,1)>=0.585))
NAc_mismatch=length(find(NAc_pabove(:,1)<=-0.585))
NAc_neutral=length(NAc_above)-NAc_match-NAc_mismatch

PAG_match = length(find(PAG_pabove(:,1)>=0.585))
PAG_mismatch = length(find(PAG_pabove(:,1)<=-0.585))
PAG_neutral = length(PAG_above)-PAG_match-PAG_mismatch

figure;
piechart([NAc_match,NAc_mismatch,NAc_neutral],["Match","Mismatch","Neutral"]);

figure;
piechart([PAG_match,PAG_mismatch,PAG_neutral],["Match","Mismatch","Neutral"]);
