function [d1, d2, ds1, ds2, ds3] = MMD3(xPos, xNeg, allSgm)
% maximum mean discrepancy (MMD) via pdist and pdist2 and block process
% vectorized version, fast and need limit memory
% See Eqn. (3-5) in [Gretton et al. 2012 JMLR]
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] A. Gretton, K. M. Borgwardt, M. J. Rasch, B. Scholkopf, and A. Smola. A kernel two-sample test. 
%      JMLR, 13(3):723–773, 2012.
%
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
% Output
%  d1: biased estimate of MMD
%  d2: unbiased estimate of MMD
%  ds1: square of biased MMD without positive constraint
%  ds2: square of unbiased MMD without positive constraint
%  ds3: square of unbiased MMD without positive constraint when (nPos == nNeg)
%       used for B-test

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 01/09/2014

MAX_SIZE = 1e6;

nDim = size(xPos, 2);
bsz = ceil(sqrt(MAX_SIZE/nDim)); % block size
nSgm = numel(allSgm);
nPos = size(xPos, 1);
nNeg = size(xNeg, 1);
if nPos==0 || nNeg==0
    d1 = [];
    d2 = [];
    return;
end
allGamma = 1 ./ (2*allSgm.^2);

s1 = zeros(nSgm, 1);
s2 = zeros(nSgm, 1);
s3 = zeros(nSgm, 1);

%%
nBlkPos = ceil(nPos/bsz);
for ii = 1:nBlkPos
    idx1 = (ii-1)*bsz+1;
    idx2 = min(ii*bsz, nPos);
    tmp = pdist(xPos(idx1:idx2, :), 'euclidean').^2;
    for kk = 1:nSgm
        gm = allGamma(kk);
        tmp2 = exp(-tmp*gm);
        s1(kk) = s1(kk) + sum(tmp2(:));
    end
end
for ii = 1:nBlkPos
    for jj = 1:(ii-1)
        idx1 = (ii-1)*bsz+1;
        idx2 = min(ii*bsz, nPos);
        idx3 = (jj-1)*bsz+1;
        idx4 = min(jj*bsz, nPos);
        tmp = pdist2(xPos(idx1:idx2,:), xPos(idx3:idx4,:), 'euclidean').^2;
        for kk = 1:nSgm
            gm = allGamma(kk);
            tmp2 = exp(-tmp*gm);
            s1(kk) = s1(kk) + sum(tmp2(:));
        end
    end
end

%%
nBlkNeg = ceil(nNeg/bsz);
for ii = 1:nBlkNeg
    idx1 = (ii-1)*bsz+1;
    idx2 = min(ii*bsz, nNeg);
    tmp = pdist(xNeg(idx1:idx2, :), 'euclidean').^2;
    for kk = 1:nSgm
        gm = allGamma(kk);
        tmp2 = exp(-tmp*gm);
        s2(kk) = s2(kk) + sum(tmp2(:));
    end
end
for ii = 1:nBlkNeg
    for jj = 1:(ii-1)
        idx1 = (ii-1)*bsz+1;
        idx2 = min(ii*bsz, nNeg);
        idx3 = (jj-1)*bsz+1;
        idx4 = min(jj*bsz, nNeg);
        tmp = pdist2(xNeg(idx1:idx2,:), xNeg(idx3:idx4,:), 'euclidean').^2;
        for kk = 1:nSgm
            gm = allGamma(kk);
            tmp2 = exp(-tmp*gm);
            s2(kk) = s2(kk) + sum(tmp2(:));
        end
    end
end

%%
for ii = 1:nBlkPos
    for jj = 1:nBlkNeg
        idx1 = (ii-1)*bsz+1;
        idx2 = min(ii*bsz, nPos);
        idx3 = (jj-1)*bsz+1;
        idx4 = min(jj*bsz, nNeg);
        tmp = pdist2(xPos(idx1:idx2,:), xNeg(idx3:idx4,:), 'euclidean').^2;
        for kk = 1:nSgm
            gm = allGamma(kk);
            tmp2 = exp(-tmp*gm);
            s3(kk) = s3(kk)+sum(tmp2(:));
        end
    end
end

%%
ds1 = (s1*2+nPos)/nPos^2 + (s2*2+nNeg)/nNeg^2 - s3*2/nPos/nNeg;
d1 = sqrt( max(ds1,0) );

ds2 = (s1/nPos/(nPos-1) + s2/nNeg/(nNeg-1) - s3/nPos/nNeg)*2;
d2 = sqrt( max(ds2,0) );

%%
if nargout==5
    if (nPos~=nNeg)
        error('The sample number of two sets should be equal!');
    end
    s4 = zeros(nSgm, 1);
    tmp = sum((xPos-xNeg).^2, 2);
    for kk = 1:nSgm
        gm = allGamma(kk);
        tmp2 = exp(-tmp*gm);
        s4(kk) = s3(kk)-sum(tmp2(:));
    end
    ds3 = (s1+s2-s4)*2/nPos/(nPos-1);
end
