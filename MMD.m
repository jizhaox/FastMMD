function [d, s] = MMD(xPos, xNeg, allSgm, opt)
% Naive implementation of maximum mean discrepancy (MMD)
% this is very slow and we suggest function MMD3
% See Eqn. (3-5) in [Gretton et al. 2012 JMLR]
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
%
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
%  opt: biased, unbiased, btest
% Output
%  d: biased estimate of MMD
%  ds: square of MMD without positive constraint

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 10/26/2013

if nargin < 4
    opt = 'biased';
end

nSgm = numel(allSgm);
nPos = size(xPos, 1);
nNeg = size(xNeg, 1);
if nPos==0 || nNeg==0
    d = [];
    return;
end
allGamma = 1 ./ (2*allSgm.^2);

s = zeros(nSgm, 1);
if strcmp(opt, 'biased')
for kk = 1:nSgm
    s1 = 0;
    s2 = 0;
    s3 = 0;
    gm = allGamma(kk);
    for ii = 1:nPos
        for jj = 1:nPos
            tmp = xPos(ii, :) - xPos(jj, :);
            s1 = s1 + exp( -sum(tmp(:).^2) * gm );
        end
    end
    for ii = 1:nNeg
        for jj = 1:nNeg
            tmp = xNeg(ii, :) - xNeg(jj, :);
            s2 = s2 + exp( -sum(tmp(:).^2) * gm );
        end
    end
    for ii = 1:nPos
        for jj = 1:nNeg
            tmp = xPos(ii, :) - xNeg(jj, :);
            s3 = s3 + exp( -sum(tmp(:).^2) * gm );
        end
    end
    s(kk) = s1/nPos^2 + s2/nNeg^2 - s3*2/nPos/nNeg;
end
    d = sqrt( max(s,0) );
end

if strcmp(opt, 'unbiased')
for kk = 1:nSgm
    s1 = 0;
    s2 = 0;
    s3 = 0;
    gm = allGamma(kk);
    for ii = 1:nPos
        for jj = 1:nPos
            if jj ~= ii
                tmp = xPos(ii, :) - xPos(jj, :);
                s1 = s1 + exp( -sum(tmp(:).^2) * gm );
            end
        end
    end
    for ii = 1:nNeg
        for jj = 1:nNeg
            if jj ~= ii
                tmp = xNeg(ii, :) - xNeg(jj, :);
                s2 = s2 + exp( -sum(tmp(:).^2) * gm );
            end
        end
    end
    for ii = 1:nPos
        for jj = 1:nNeg
            tmp = xPos(ii, :) - xNeg(jj, :);
            s3 = s3 + exp( -sum(tmp(:).^2) * gm );
        end
    end
    s(kk) = s1/nPos/(nPos-1) + s2/nNeg/(nNeg-1) - s3*2/nPos/nNeg;
end
    d = sqrt( max(s,0) );
end

if strcmp(opt, 'btest')
if nPos ~= nNeg % nPos = nNeg
    error('The sample number of two sets should be equal!');
end
for kk = 1:nSgm
    s1 = 0;
    s2 = 0;
    s3 = 0;
    gm = allGamma(kk);
    for ii = 1:nPos
        for jj = 1:nPos
            if jj ~= ii
                tmp = xPos(ii, :) - xPos(jj, :);
                s1 = s1 + exp( -sum(tmp(:).^2) * gm );
            end
        end
    end
    for ii = 1:nNeg
        for jj = 1:nNeg
            if jj ~= ii
                tmp = xNeg(ii, :) - xNeg(jj, :);
                s2 = s2 + exp( -sum(tmp(:).^2) * gm );
            end
        end
    end
    for ii = 1:nPos
        for jj = 1:nNeg
            if jj~=ii %% notice here
                tmp = xPos(ii, :) - xNeg(jj, :);
                s3 = s3 + exp( -sum(tmp(:).^2) * gm );
            end
        end
    end
    s(kk) = s1/nPos/(nPos-1) + s2/nNeg/(nNeg-1) - s3*2/nPos/(nNeg-1);
end
    d = sqrt( max(s,0) );
end
