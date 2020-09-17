function [THT, PHI] = Fastfood(x, para, sgm)
% FastFood for apprimating Gaussian kernels
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
% [1] Q. V. Le, T. Sarlos, and A. J. Smola. 
%     Fastfood - Approximating Kernel Expansions in Loglinear Time. ICML, 2013.
% Input:
%  x: samples, each row is a sample
%  para: parameters for Fastfood.
%  sgm: bandwidth parameter for Gaussian kernel
% Output
%  THT: angles used for feature mapping
%  PHI: feature mapping for input x
%
% Example
% d = 128; % dimension of data
% n = d*10; % number of basis for approximation
% x1 = rand(d, 1); % input pattern
% x2 = rand(d, 1);
% sgm = 5;
% k1 = exp( -norm(x1-x2,2)^2/(2*sgm^2) ) % exact kernel
% para = FastfoodPara(d, n);
% [~, PHI1] = Fastfood(x1, para, sgm);
% [~, PHI2] = Fastfood(x2, para, sgm);
% k2 = PHI1'*PHI2 % approximation of kernel by fastfood

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 12/19/2013

if nargin < 3
    sgm = 1;
end

%% preprocessing parameter
[d1, m]= size(x);
l = ceil(log2(d1));
d = 2^l;
if d == d1
    xx = x;
else
    xx = zeros(d, m);
    xx(1:d1, :) = x; 
end
    
k = numel(para.B);
n = d * k;
THT = zeros(n, m);

try
    % test whether we can use fwht_spiral
    fwht_spiral([1; 1]);
    for ii = 1:k
        B = para.B{ii};
        G = para.G{ii};
        II = para.II{ii};
        xx = bsxfun(@times, xx, B);
        t = fwht_spiral(xx);
        t = t(II, :);
        t = bsxfun(@times, t, G);
        THT(((ii-1)*d+1):(ii*d), :) = fwht_spiral(t);
    end
    S = para.S;
    THT = bsxfun(@times, THT, S/d^(1/2));
catch
    fprintf(1, 'Cannot perform Walsh–Hadamard transform using Spiral.\n');
    fprintf(1, 'Use Matlab function fwht instead. It is much slower.\n');
    for ii = 1:k
        B = para.B{ii};
        G = para.G{ii};
        II = para.II{ii};
        xx = bsxfun(@times, xx, B);
        t = fwht(xx, d, 'hadamard');
        t = t(II, :);
        t = bsxfun(@times, t, G*d);
        THT(((ii-1)*d+1):(ii*d), :) = fwht(t, d, 'hadamard');
    end
    S = para.S;
    THT = bsxfun(@times, THT, S*d^(1/2));
end

if nargout==2
    t = THT/sgm;
    PHI = [cos(t); sin(t)] * n^(-1/2);
end
