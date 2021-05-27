clc;
clear;
close all;

%% Set seed
rng(0);

%% Read the image
orig = cast(imread("data/barbara256.png"),'double');

%% Constants
% Set the Height and Width of the image
H = size(orig, 1);
W = size(orig, 2);
% Set Patch size
ps = 8;

% Define the sensing matrix as iid Gaussian Matrix
phi = randn(ps*ps/2, ps*ps);

%% Reconstruction of the original image
% Define the orthonormal matrix in which the patches are sparse - here, 2D-DCT
psi = kron(dctmtx(ps)', dctmtx(ps)');

% Define the sensing matrix w.r.t. DCT coefficients
A = phi * psi;

% Set alpha, lambda and number of iterations for ISTA
alpha = floor(eigs(A'*A, 1)) + 2;
lambda = 3;
iter = 100;

% Initialize reconstructed image and averaging matrix
recon = zeros(H, W, 'double');
avg_mat = zeros(H, W, 'double');

tic;
% For every (overlapping) patch
for i=1:H-ps+1
    for j=1:W-ps+1
        % Get the compressed measurement of the patch
        y = phi * reshape(orig(i:i+ps-1,j:j+ps-1), [8*8 1]);

        % Use ISTA to obtain the DCT coefficients
        theta = ista(y, A, lambda, alpha, iter);

        % Update the reconstructed patch from the coefficients
        recon(i:i+ps-1,j:j+ps-1) = recon(i:i+ps-1,j:j+ps-1) + reshape(psi * theta, [8 8]);
        avg_mat(i:i+ps-1,j:j+ps-1) = avg_mat(i:i+ps-1,j:j+ps-1) + ones(8,8);

        % Print the co-ordinates of the patch, to check for speed and debugging
        fprintf('(%i, %i)\n', i, j);
    end
end

%% Normalize the reconstructed frames
recon = 2 * recon ./ avg_mat;
recon = cast(recon, 'uint8');
orig = cast(orig, 'uint8');

%% Save the result and Compute RMSE (Relative Mean Squared Error)
% Display and Save the reconstructed frame
figure; imshow([orig, recon]);
imwrite(orig, 'results/q2b_orig.png');
imwrite(recon, 'results/q2b_recon.png');
% RMSE of the reconstructed image
fprintf('RMSE : %f\n', (norm(double(orig) - double(recon), 'fro')^2 / norm(double(orig), 'fro')^2));

% Evaluate the time taken
toc;
