function para = FastfoodPara(d, n)
% Generate Fastfood parameters for apprimating Gaussian kernels
% The kernel is Gaussian: k(x,y) = exp(-||x-y||_2^2 / (2*sigma^2)).
% Reference:
% [1] Q. V. Le, T. Sarlos, and A. J. Smola. 
%     Fastfood - Approximating Kernel Expansions in Loglinear Time. ICML, 2013.
% Input:
%  d: dimension of data
%  n: number of basis for approximation
% Output:
%  para: parameters for Fastfood.

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 12/19/2013

if nargin<2
    n = d*10;
end
%% preprocessing parameter
% pad the vectors with zeros until d = 2^l holds.
l = ceil(log2(d));
d = 2^l;
k = ceil(n/d);
n = d * k;

B = cell(k, 1);
G = cell(k, 1);
II = cell(k, 1);
S = cell(k, 1);
for ii = 1:k
    B{ii} = randsrc(d, 1, [1 -1]);
    G{ii} = randn(d, 1);
    t = randperm(d);
    II{ii} = t(:);
    t = gammaincinv(rand(d,1), d/2, 'lower');
    t = (t*2).^(1/2);
    S{ii} = t * norm(G{ii}, 'fro')^(-1);
end
S1 = zeros(n, 1);
for ii = 1:k
    S1(((ii-1)*d+1):(ii*d)) = S{ii};
end

para.B = B;
para.G = G;
para.II = II;
para.S = S1;
