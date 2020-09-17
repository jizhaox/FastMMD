% FastMMD
% Ji Zhao@CMU
% zhaoji84@gmail.com
% 10/26/2013

clear;

allSgm = 10.^(-2:0.2:2);
allSgm = allSgm(:);
nBasis = 2^10;

%% generate data
[X, Y] = GenSamp(1); % 1 -- is_rand
xPos = X(Y == 1, :);
xNeg = X(Y == -1, :);

fprintf(1, '-------------beginning-------------\n')
%% MMD
%tic, d1 = MMD(xPos, xNeg, allSgm, 'biased'); f1 = MMD(xPos, xNeg, allSgm, 'unbiased'); toc
tic, [d1, f1] = MMD3(xPos, xNeg, allSgm); toc

%% FastMMD via Random Fourier Feature
tic, [d2, f2] = MMDFourierFeature(xPos, xNeg, allSgm, nBasis); toc

%% FastMMD via Fastfood
tic, [d3, f3] = MMDFastfood(xPos, xNeg, allSgm, nBasis); toc

%% block test (B-test)
tic, f4 = MMDBlockTest(xPos, xNeg, allSgm); toc

%% MMD-linear
tic, f5 = MMDlinear(xPos, xNeg, allSgm); toc

%% FastMMD by calculating the amplitude of sinusoids
%[d6, f6] = MMDFourierSineOnline(X, Y, allSgm, nBasis, true);

%%
figure, semilogx(allSgm, d1, allSgm, d2, allSgm, d3);
legend('MMD-biased', 'FastMMD-Fourier', 'FastMMD-Fastfood')
xlabel('\sigma in Gaussian Kernel')
ylabel('Maximum Mean Discrepancy')
axis([min(allSgm(:)) max(allSgm(:)) 0 0.5])

figure, semilogx(allSgm, f1, allSgm, f2, allSgm, f3, allSgm, f4, allSgm, f5);
legend('MMD-unbiased', 'FastMMD-Fourier', 'FastMMD-Fastfood', 'B-test', 'MMD-linear');
xlabel('\sigma in Gaussian Kernel')
ylabel('Maximum Mean Discrepancy')
axis([min(allSgm(:)) max(allSgm(:)) 0 0.5])