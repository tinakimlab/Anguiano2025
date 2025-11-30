% Plot Nrxn-TurboID vs TurboID-CAAX exclusive proteins

% data is stored with 1st column = L2fc, 2nd column = QVAL, and 3rd
% column = # peptides found.

% "nrxnhits_caaxvals" = proteins that were identified as true candidates in
% Nrxn-TurboID PFC +biotin vs -biotin dataset, and then looked for in the
% TurboID-CAAX PFC +biotin vs -biotin dataset. Vice versa for
% "caaxhits_nrxnvals". These only include proteins that were detected in
% the other dataset.

% Nrxn-TurboID: 1,513 exclusive proteins (not found in filtered
% TurboID-CAAX candidate list). 938 of those proteins were not identified
% at all in the TurboID-CAAX +biotin vs -biotin dataset, and are not
% included here.

% TurboID-CAAX: 83 exclusive proteins (not found in filtered Nrxn-TurboID
% candidate list). 42 of those proteins were not identifeid at all in the
% Nrxn-TurboID +biotin vs -biotin dataset, and are not included here.

close all
clear all

load('nrxn_caax_exclusive_comparison.mat')

%% Nrxn exclusive in CAAX dataset

nrxn_trueEXC=883; % not found anywhere in CAAX datasets
nrxn_sig_subthreshL2FC=length(find(nrxnhits_caaxvals(:,2)<0.05 & nrxnhits_caaxvals(:,1)>0 & nrxnhits_caaxvals(:,1)<0.585));
nrxn_sig_negL2FC=length(find(nrxnhits_caaxvals(:,2)<0.05 & nrxnhits_caaxvals(:,1)<0));
nrxn_notsig=length(find(nrxnhits_caaxvals(:,2)>=0.05));
nrxn_sig_L2FC_1pep=length(find(nrxnhits_caaxvals(:,2)<0.05 & nrxnhits_caaxvals(:,1)>=0.585 & nrxnhits_caaxvals(:,3)<2));

figure;
piechart([nrxn_trueEXC, nrxn_notsig, nrxn_sig_negL2FC, nrxn_sig_subthreshL2FC,]);


figure;
plot(nrxnhits_caaxvals(:,1),-log10(nrxnhits_caaxvals(:,2)),'.k');
hold on
plot([-5 5],[-log10(0.05) -log10(0.05)],'--k')
plot([0.585 0.585],[0 4.5],'--k')
plot([-0.585 -0.585],[0 4.5],'--k')
ylim([0.5 4.5])
xlim([-5 5])

%% CAAX exclusive in Nrxn dataset
CAAX_trueEXC=43; % not found anywhere in Nrxn datasets
CAAX_sig_subthreshL2FC=length(find(CAAXhits_nrxnvals(:,2)<0.05 & CAAXhits_nrxnvals(:,1)>0 & CAAXhits_nrxnvals(:,1)<0.585));
CAAX_sig_negL2FC=length(find(CAAXhits_nrxnvals(:,2)<0.05 & CAAXhits_nrxnvals(:,1)<0));
CAAX_notsig=length(find(CAAXhits_nrxnvals(:,2)>=0.05));
CAAX_sig_L2FC_1pep=length(find(CAAXhits_nrxnvals(:,2)<0.05 & CAAXhits_nrxnvals(:,1)>=0.585 & CAAXhits_nrxnvals(:,3)<2));

figure;
piechart([CAAX_trueEXC, CAAX_sig_negL2FC, CAAX_notsig, CAAX_sig_subthreshL2FC, CAAX_sig_L2FC_1pep]);

figure;
plot(CAAXhits_nrxnvals(:,1),-log10(CAAXhits_nrxnvals(:,2)),'.k');
hold on
plot([-8 8],[-log10(0.05) -log10(0.05)],'--k')
plot([0.585 0.585],[0 3.5],'--k')
plot([-0.585 -0.585],[0 3.5],'--k')
ylim([1 3.5])
xlim([-8 8])
