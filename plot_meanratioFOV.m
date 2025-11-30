close all
clear all

load('alldata_meanratio_FOV.mat');
% data matrices is stored in the following order by column:
% TurboID, Nrxn-Turbo, Nrxn-CIBN-Turbo, CD4-CIBN-Turbo, CAAX1, CAAX2

figure;
hold on;
for a=1:6
    currdata=[alldata_nobiotin(:,a) alldata_biotin(:,a)];
    currdata=[currdata(~isnan(currdata(:,1)),1) currdata(~isnan(currdata(:,2)),2)];
 
    bar([a*2-1 a*2], mean(currdata));
    errorbar([a*2-1 a*2], mean(currdata), std(currdata)./sqrt(length(currdata)));
    scatter(ones(length(currdata),1)*(a*2-1),currdata(:,1),'jitter',0.1);
    scatter(ones(length(currdata),1)*a*2,currdata(:,2),'jitter',0.1);

end