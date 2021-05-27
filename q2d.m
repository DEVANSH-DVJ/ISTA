clc;
clear;
close all;

% Setting the seed
rng(100);

% Reading
orig = zeros(100,1);
ind = randi(100, [10 1]);
orig(ind) = randi(10, [10 1]);

A = zeros(100, 100);
A(1:4, 1) = [10 3 2 1];
A(1:5, 2) = [6 4 3 2 1];
A(1:6, 3) = [3 3 4 3 2 1];
for i=4:97
     A(i-3:i+3,i) = [1 2 3 4 3 2 1];
end
A(95:100, 98) = [1 2 3 4 3 3];
A(96:100, 99) = [1 2 3 4 6];
A(97:100, 100) = [1 2 3 10];
A = A/16;

alpha = floor(eigs(A'*A,1)) + 1;
iter = 10000;
lambda = 0.01;

y = A*orig;
recon = ista(y, A, lambda, alpha, iter);
fprintf('RMSE : %f\n', norm(recon - orig)/norm(orig));
