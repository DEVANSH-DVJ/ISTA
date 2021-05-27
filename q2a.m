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
% Set standard deviation of Gaussian Noise
noise_std = 2;
% Set Patch size
ps = 8;

%% Get the noisy image using additive gaussian noise
noisy = orig + noise_std*randn(H, W);

%% Reconstruction of the original image
% Define the orthonormal matrix in which the patches are sparse - here, 2D-DCT
psi = kron(dctmtx(ps)', dctmtx(ps)');

% Define the sensing matrix as Identity
phi = diag(ones(ps*ps, 1));

% Define the sensing matrix w.r.t. DCT coefficients
A = phi * psi;

% Set alpha, lambda and number of iterations for ISTA
alpha = floor(eigs(A'*A, 1)) + 2;
lambda = 3 * noise_std;
iter = 100;

% Initializing reconstructed image and averaging matrix
recon = zeros(H, W, 'double');
avg_mat = zeros(H, W, 'double');

tic;
% For every (overlapping) patch
for i=1:H-ps+1
    for j=1:W-ps+1
        % Get the noisy patch
        y = phi * reshape(noisy(i:i+ps-1,j:j+ps-1), [ps*ps 1]);

        % Use ISTA to obtain the DCT coefficients
        theta = ista(y, A, lambda, alpha, iter);

        % Update the reconstructed patch from the coefficients
        recon(i:i+ps-1, j:j+ps-1) = recon(i:i+ps-1, j:j+ps-1) + reshape(psi * theta, [ps ps]);
        avg_mat(i:i+ps-1, j:j+ps-1) = avg_mat(i:i+ps-1, j:j+ps-1) + ones(ps, ps);

        % Print the co-ordinates of the patch, to check for speed and debugging
        fprintf('(%i, %i)\n', i, j);
    end
end

%% Normalize the reconstructed frames
recon = recon ./ avg_mat;
recon = cast(recon, 'uint8');
orig = cast(orig, 'uint8');

%% Save the result and Compute RMSE (Relative Mean Squared Error)
% Display and Save the reconstructed frame
figure; imshow([orig, recon]);
imwrite(orig, 'results/q2a_orig.png');
imwrite(recon, 'results/q2a_recon.png');
% RMSE of the reconstructed image
fprintf('RMSE : %f\n', (norm(double(recon - orig), 'fro')^2 / norm(double(orig), 'fro')^2));

% Evaluate the time taken
toc;
