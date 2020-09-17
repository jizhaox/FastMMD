function [d1, d2] = MMDFourierFeature(xPos, xNeg, allSgm, nBasis)
% Approximate the MMD by Random Fourier Features (Random Kitchen Sinks)
%  with block process
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] Ji Zhao, Deyu Meng. Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
%      NIPS Workshop on Randomized Methods for Machine Learning (RMML2013), 2013.
%  [2] Ji Zhao, Deyu Meng. FastMMD: Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
%      Neural Computation, 2015.
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
%  nBasis: number of basis for approximating p(w), see our paper.
% Output
%  d1: estimate of biased MMD
%  d2: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 02/18/2014

MAX_SIZE = 1e7;

if 0
    rng('default');
end

k0 = 1; % K(0,0)=1
nDim = size(xPos, 2);
%W = randn(nBasis, nDim);
d1 = zeros(numel(allSgm), 1);
d2 = zeros(numel(allSgm), 1);

nPos = size(xPos, 1);
nNeg = size(xNeg, 1);

bsz = max(ceil(MAX_SIZE/nDim), 1);

%%
nBlock1 = ceil(nBasis/bsz); % for W
nBlock2 = ceil(nPos/bsz); % for dataset 1
nBlock3 = ceil(nNeg/bsz); % for dataset 2
phiPos = zeros(nBasis*2, numel(allSgm));
phiNeg = zeros(nBasis*2, numel(allSgm));

for ii = 1:nBlock1
    i1 = (ii-1)*bsz + 1;
    i2 = min(ii*bsz, nBasis);
    i1 = i1*2-1; % need double space to store [cos sin]
    i2 = i2*2;
    tmp = min(bsz, nBasis-(ii-1)*bsz);
    W = randn(tmp, nDim);
    for jj = 1:nBlock2
        j1 = (jj-1)*bsz + 1;
        j2 = min(jj*bsz, nPos);
        t = W*xPos(j1:j2, :)';
        for kk = 1:numel(allSgm)
            sgm = allSgm(kk);
            t1 = t/sgm;
            t2 = [cos(t1); sin(t1)];
            phiPos(i1:i2, kk) = phiPos(i1:i2, kk) + sum(t2,2);
        end
    end
    for jj = 1:nBlock3
        j1 = (jj-1)*bsz + 1;
        j2 = min(jj*bsz, nNeg);
        t = W*xNeg(j1:j2, :)';
        for kk = 1:numel(allSgm)
            sgm = allSgm(kk);
            t1 = t/sgm;
            t2 = [cos(t1); sin(t1)];
            phiNeg(i1:i2, kk) = phiNeg(i1:i2, kk) + sum(t2,2);
        end
    end
end
phiPos = phiPos/nPos * nBasis^(-1/2);
phiNeg = phiNeg/nNeg * nBasis^(-1/2);

%%
for ii = 1:numel(allSgm)
    d1(ii) = norm(phiPos(:, ii)-phiNeg(:, ii), 2);
    t = d1(ii)^2 + norm(phiPos(:, ii),2)^2/(nPos-1) + norm(phiNeg(:, ii),2)^2/(nNeg-1) - k0*(nPos+nNeg-2)/(nPos-1)/(nNeg-1);
    d2(ii) = (max(t,0))^(1/2);
end

