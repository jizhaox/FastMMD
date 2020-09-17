function [X, Y] = GenSamp(is_rand)
% Output:
%  X: samples, each row is a sample
%  Y: labels, +1 or -1

% Ji Zhao@CMU
% zhaoji84@gmail.com
% 2011

%%
if ~is_rand
    rng('default');
end
N = 320; % sample number

X = rand(N, 2);
X = (X - 0.5) * 10; %% [-5, 5]

Y = -ones(N, 1);
t = sqrt(sum(X.^2, 2));
Y( (t>1 & t<4) ) = 1;

% two sets have the same number of sampes
idx1 = find(Y==1);
idx2 = find(Y==-1);
N1 = numel(idx1);
N2 = numel(idx2);

if N1>N2
    X(idx1(N2+1:end), :) = [];
    Y(idx1(N2+1:end), :) = [];
elseif N2>N1
    X(idx2(N1+1:end), :) = [];
    Y(idx2(N1+1:end), :) = [];
end
