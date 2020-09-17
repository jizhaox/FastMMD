function [d1, d2] = MMDFastfood(xPos, xNeg, allSgm, nBasis)
% Approximate the MMD by Fastfood
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] Ji Zhao, Deyu Meng. FastMMD: Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
%  nBasis: number of basis for approximating p(w), see our paper.
% Output
%  d1: estimate of biased MMD
%  d2: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 12/19/2013

MAX_SIZE = 1e7;

nSgm = numel(allSgm);
k0 = 1; % K(0,0)=1
nDim = size(xPos, 2);
d1 = zeros(nSgm, 1);
d2 = zeros(nSgm, 1);
nPos = size(xPos, 1);
nNeg = size(xNeg, 1);

para = FastfoodPara(nDim, nBasis);

n = numel(para.S);
bsz = max(ceil(MAX_SIZE/n/2), 1);

%%
nBlock = ceil(nPos/bsz);
phiPos = zeros(n*2, nSgm);
thtPos = Fastfood(xPos', para);
for ii = 1:nBlock
    i1 = (ii-1)*bsz + 1;
    i2 = min(ii*bsz, nPos);
    t = thtPos(:, i1:i2);
    for jj = 1:nSgm
        sgm = allSgm(jj);
        t1 = t/sgm;
        t2 = [cos(t1); sin(t1)];
        phiPos(:, jj) = phiPos(:, jj) + sum(t2,2);
    end
end
phiPos = phiPos/nPos * n^(-1/2);
clear thtPos;

%%
nBlock = ceil(nNeg/bsz);
phiNeg = zeros(n*2, nSgm);
thtNeg = Fastfood(xNeg', para);
for ii = 1:nBlock
    i1 = (ii-1)*bsz + 1;
    i2 = min(ii*bsz, nPos);
    t = thtNeg(:, i1:i2);
    for jj = 1:nSgm
        sgm = allSgm(jj);
        t1 = t/sgm;
        t2 = [cos(t1); sin(t1)];
        phiNeg(:, jj) = phiNeg(:, jj) + sum(t2,2);
    end
end
phiNeg = phiNeg/nNeg * n^(-1/2);

%%
for ii = 1:nSgm
    d1(ii) = norm(phiPos(:, ii)-phiNeg(:, ii), 2);
    t = d1(ii)^2 + norm(phiPos(:, ii),2)^2/(nPos-1) + norm(phiNeg(:, ii),2)^2/(nNeg-1) - k0*(nPos+nNeg-2)/(nPos-1)/(nNeg-1);
    d2(ii) = (max(t,0))^(1/2);
end
