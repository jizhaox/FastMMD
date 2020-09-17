function d = MMDBlockTest(xPos, xNeg, allSgm, blocksize)
% MMD by B-tests
% based on the code of Wojciech Zaremba
% https://github.com/wojzaremba/btest/
% Reference: 
% [1] W. Zaremba, A. Gretton and M. Blaschko.
%     B-tests: Low Variance Kernel Two-Sample Tests. NIPS, 2013.
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
% Output
%  d: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 12/23/2013

nSgm = numel(allSgm);
m = min( size(xPos, 1), size(xNeg, 1));
if nargin<4
    blocksize = floor(sqrt(m));
end
m2 = floor(m / blocksize); % number of blocks
hh = zeros(m2, nSgm);
allGamma = 1 ./ (2*allSgm.^2);

for x = 1 : blocksize
    for y = 1 : blocksize
        if (x ~= y)
            idx1 = ((m2 * (x - 1)) + 1) : (m2 * x);
            idx2 = ((m2 * (y - 1)) + 1) : (m2 * y);
            dists = sum( (xPos(idx1, :)-xPos(idx2, :)).^2, 2);
            for k = 1:nSgm
                gm = allGamma(k);
                H = exp(-dists*gm);
                hh(:, k) = hh(:, k) + H;
            end
            dists = sum( (xNeg(idx1, :)-xNeg(idx2, :)).^2, 2);
            for k = 1:nSgm
                gm = allGamma(k);
                H = exp(-dists*gm);
                hh(:, k) = hh(:, k) + H;
            end
            dists = sum( (xPos(idx1, :)-xNeg(idx2, :)).^2, 2);
            for k = 1:nSgm
                gm = allGamma(k);
                H = exp(-dists*gm);
                hh(:, k) = hh(:, k) - H;
            end
            dists = sum( (xNeg(idx1, :)-xPos(idx2, :)).^2, 2);
            for k = 1:nSgm
                gm = allGamma(k);
                H = exp(-dists*gm);
                hh(:, k) = hh(:, k) - H;
            end
        end
    end
end
d = mean(hh, 1);
d = max(d, 0);
d = sqrt(d/blocksize/(blocksize-1));
d = d(:);
