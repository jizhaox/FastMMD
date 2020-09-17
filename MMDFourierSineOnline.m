function [d1, d2] = MMDFourierSineOnline(X, Y, allSgm, nBasis, isPlot)
% Approximate the MMD by feature maps and sampling
% This is an online version
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
%  [1] Ji Zhao, and Deyu Meng. FastMMD: Ensemble of Circular Discrepancy for Efficient Two-Sample Test.
%      Neural Computation, 27(6): 1345 - 1372, 2015.
%
% Input:
%  X: samples, each row is a sample
%  Y: labels, +1 or -1
%  allSgm: bandwidth parameter for Gaussian kernel, scale or vector
%  nBasis: number of basis for approximating p(w), see our paper.
%  isPlot: flag for visualization
% Output
%  d1: estimate of biased MMD
%  d2: estimate of unbiased MMD

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 11/07/2013

if nargin<5
    isPlot = 0;
end
if 0
    rng('default');
end

N = numel(Y);
nSgm = numel(allSgm);
nDim = size(X, 2);
W = randn(nBasis, nDim);

aPos = cell(nSgm,1);
thPos = cell(nSgm,1);
aNeg = cell(nSgm,1);
thNeg = cell(nSgm,1);
for ii = 1:nSgm
    aPos{ii} = zeros(nBasis,1);
    thPos{ii} = zeros(nBasis,1);
    aNeg{ii} = zeros(nBasis,1);
    thNeg{ii} = zeros(nBasis,1);
end
nPos = 0;
nNeg = 0;
if isPlot
    figure;
    cnt = 0;
end
for jj = 1:N
    [d1, d2, aPos, thPos, aNeg, thNeg, nPos, nNeg] = ...
        CalcuAmpOnline(aPos, thPos, aNeg, thNeg, nPos, nNeg, X(jj,:), Y(jj), W, allSgm);

    if isPlot && mod(jj, 20)==0
        cnt = cnt + 1;
        semilogx(allSgm, d1, allSgm, d2); %hold on
        xlabel('\sigma in Gaussian Kernel')
        ylabel('Maximum Mean Discrepancy')
        legend('MMD-biased', 'MMD-unbiased')
        axis([min(allSgm(:)) max(allSgm(:)) 0 0.5])
        title(['sample number: ', num2str(jj)]);
        drawnow;
        %filename = 'onlineMMD.gif';
        %frame = getframe(1);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        %if cnt == 1;
        %    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        %else
        %    imwrite(imind,cm,filename,'gif','WriteMode','append');
        %end
    end
end




