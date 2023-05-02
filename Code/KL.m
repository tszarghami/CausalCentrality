function [dKL] = KL(p,q)

% Kullback Leibler divergence (KL) KL(p,q)
% Input: 
%     p: p.mu, p.sigma
%     q: q.mu, q.sigma
% Output: 
%     dKL: the KL divergence between p and q

% Copyright (c) 2014, Meizhu Liu (meizhu.liu{at}yahoo.com)
% All rights reserved.
% See LICENSE file
% Source: https://www.mathworks.com/matlabcentral/fileexchange/46090-kl-divergence-between-gaussian-distributions

% I have slightly modified Meizhu's original code --> pinv and regularization

mu0 = p.mu;
mu1 = q.mu;
k = length(mu0);
sigma0 = p.sigma;
sigma1 = q.sigma;

% tmp = inv(sigma1)*sigma0;
tmp = pinv(sigma1)*sigma0;

% dKL = 0.5*(trace(tmp)+(mu1-mu0)'*inv(sigma1)*(mu1-mu0)-k-log(det(tmp)));
dKL = 0.5*(trace(tmp)+(mu1-mu0)'*pinv(sigma1)*(mu1-mu0)-k-log(det(tmp)+1e-128));% modification: pinv & regularize
