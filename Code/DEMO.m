% Causal Centrality Demo

% MIT License
% Copyright (c) 2023 Tahereh S. Zarghami (tszarghami{at}gmail.com)
% Date: 14 April 2023
% The code has dependencies on SPM: https://www.fil.ion.ucl.ac.uk/spm/

%% Import DCMs from Data folder
load('CSDs.mat')

%% PEB + Optimize with exploratory BMR

PEB = spm_dcm_peb(CSDs(:));
PEB_BMR = spm_dcm_bmr_all(PEB);
rE = PEB_BMR.Ep;
rC = PEB_BMR.Cp;
PEB_BMR = spm_dcm_reduce(PEB,rE,rC);

%% Causal Centrality

[CCent_PEB] = CausalCent(CSDs,PEB_BMR);

%% Plot

figure;
n = size(CCent_PEB,1); % #nodes
subplot(1,2,1)
Ap = reshape(PEB_BMR.Ep,n,n);
imagesc(Ap);axis square
title('Causal graph')
colormap(bluewhitered);colorbar

subplot(1,2,2)
b = bar(CCent_PEB,'k');
xticks(1:n);axis square;xlabel('Node')
title('Causal Centrality')
