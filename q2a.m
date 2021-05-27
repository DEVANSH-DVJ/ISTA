clc;
clear;
close all;

% Setting seed
rng(0);

% Reading
orig = cast(imread("data/barbara256.png"),'double');
H = size(orig, 1);
W = size(orig, 2);
% figure; imshow(cast(orig, 'uint8'));

% Adding Gaussian Noise
noise_img = orig + 2*randn(H, W);
% figure; imshow(cast(noise_img, 'uint8'));

% Calculating phi, psi and thus, A.
phi = diag(ones(64,1));
psi = kron(dctmtx(8)', dctmtx(8)');
A = phi*psi;

% Setting alpha, lambda and number of iterations
alpha = floor(eigs(A'*A,1)) + 1;
lambda = 1;
iter = 100;

% Initializing reconstructed image and averaging matrix
recon_img = zeros(H, W, 'double');
avg_mat = zeros(H, W, 'double');

tic; % Timer start

% Iterating over all possible 8x8 patches in the image
for i=1:H-7
    for j=1:W-7
        y = phi * reshape(noise_img(i:i+7,j:j+7), [8*8 1]);
        theta = ista(y, A, lambda, alpha, iter);
        recon_img(i:i+7,j:j+7) = recon_img(i:i+7,j:j+7) + reshape(psi * theta, [8 8]);
        avg_mat(i:i+7,j:j+7) = avg_mat(i:i+7,j:j+7) + ones(8,8);
        i, j % Prints the coordinates, to check for speed and debugging
    end
end

% Normalize the reconstructed image
recon_img(:,:) = recon_img(:,:)./avg_mat(:,:);
recon_img(recon_img < 0) = 0;
recon_img(recon_img > 255) = 255;

% Save the image and calculate RMSE
figure; imshow(cast([recon_img(:,:), orig(:,:)], 'uint8'));
imwrite(cast([recon_img(:,:), orig(:,:)], 'uint8'), 'results/q2a.png');
fprintf('RMSE : %f\n', norm(recon_img(:,:) - orig(:,:), 'fro')/norm(orig(:,:), 'fro'));


toc; % Timer end
