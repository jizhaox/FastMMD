function d = MMDlinear(xPos, xNeg, allSgm)
% Approximate the MMD by linear time
% See Lemma 14 of [Gretton et al. 2012 JMLR]
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] A. Gretton, K. M. Borgwardt, M. J. Rasch, B. Scholkopf, and A. Smola. A kernel two-sample test. 
%      JMLR, 13(3):723–773, 2012.
%
% Input:
%  xPos, xNeg: two sample sets
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
% Output
%  d: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 11/03/2013

nSgm = numel(allSgm);
N = min(size(xPos, 1), size(xNeg, 1));
N = floor(N/2);
allGamma = 1 ./ (2*allSgm.^2);
d = zeros(nSgm, 1);
for kk = 1:nSgm
    gm = allGamma(kk);
    s = 0;
    for ii = 1:N
        t = xPos(ii*2-1, :) - xPos(ii*2, :);
        s = s + exp( -sum(t(:).^2) * gm );
        t = xNeg(ii*2-1, :) - xNeg(ii*2, :);
        s = s + exp( -sum(t(:).^2) * gm );
        t = xPos(ii*2-1, :) - xNeg(ii*2, :);
        s = s - exp( -sum(t(:).^2) * gm );
        t = xPos(ii*2, :) - xNeg(ii*2-1, :);
        s = s - exp( -sum(t(:).^2) * gm );
    end
    s = max(s, 0);
    d(kk) = sqrt(s/N);
end
