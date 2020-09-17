function [D1, D2, aPos, thPos, aNeg, thNeg, nPos, nNeg] = CalcuAmpOnline(aPos, thPos, aNeg, thNeg, nPos, nNeg, X, Y, W, allSgm)
% Online calculate the amplitude of the combined sine function with same frequency
% Approximate the MMD by feature maps and sampling
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] Ji Zhao, Deyu Meng. Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
%      NIPS Workshop on Randomized Methods for Machine Learning (RMML2013), 2013.
%  [2] Ji Zhao, Deyu Meng. FastMMD: Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
%      Neural Computation, 2015.
%
% Output:
%  D1: estimate of biased MMD
%  D2: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 11/07/2013

nSgm = numel(allSgm);
D1 = zeros(nSgm, 1);
D2 = zeros(nSgm, 1);
if Y == 1
    nPos = nPos+1;
elseif Y == -1
    nNeg = nNeg+1;
end
for ii = 1:nSgm
    % from original space to Fourier space
    w = W/allSgm(ii);
    th = w * X';
    if Y == 1
        tmp1 = sqrt(aPos{ii}.^2 + 1 + 2*aPos{ii}.*cos(thPos{ii}-th));
        tmp2 = atan2(aPos{ii}.*sin(thPos{ii})+sin(th), aPos{ii}.*cos(thPos{ii})+cos(th));
        aPos{ii} = tmp1;
        thPos{ii} = tmp2;
    elseif Y == -1
        tmp1 = sqrt(aNeg{ii}.^2 + 1 + 2*aNeg{ii}.*cos(thNeg{ii}-th));
        tmp2 = atan2(aNeg{ii}.*sin(thNeg{ii})+sin(th), aNeg{ii}.*cos(thNeg{ii})+cos(th));
        aNeg{ii} = tmp1;
        thNeg{ii} = tmp2;
    end
    if nPos > 0 && nNeg > 0
        s = aPos{ii} / nPos;
        t = -aNeg{ii} / nNeg;
        D = s.^2 + t.^2 + 2*s.*t.*cos(thPos{ii}-thNeg{ii});
        D1(ii) = sqrt( mean(D) );
    end
    if nPos > 1 && nNeg >1
        %% unbiased estimation of MMD
        bsAll = mean(D); % squared biased MMD for all samples
        bsPos = mean(s.^2); % squared biased MMD for all positive samples
        bsNeg = mean(t.^2); % squared biased MMD for all negative samples
        k0 = 1; % K(0) = 1;
        tmp =  bsAll + (bsPos-k0)/(nPos-1) + (bsNeg-k0)/(nNeg-1);
        D2(ii) = sqrt(max(tmp,0));
        
    end
end
if nPos==0 || nNeg==0
    D1 = [];
end
if nPos<=1 || nNeg<=1
    D2 = [];
end
