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

% Calculating phi, psi and thus, A.
phi = randn(32, 64);
psi = haarmtx(64);
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
        y = phi * reshape(orig(i:i+7,j:j+7), [8*8 1]);
        theta = ista(y, A, lambda, alpha, iter);
        recon_img(i:i+7,j:j+7) = recon_img(i:i+7,j:j+7) + reshape(psi * theta, [8 8]);
        avg_mat(i:i+7,j:j+7) = avg_mat(i:i+7,j:j+7) + ones(8,8);
        i, j % Prints the coordinates, to check for speed and debugging
    end
end

% Normalize the reconstructed image
recon_img(:,:) = 2*recon_img(:,:)./avg_mat(:,:);
recon_img(recon_img < 0) = 0;
recon_img(recon_img > 255) = 255;

% Save the image and calculate RMSE
figure; imshow(cast([recon_img(:,:), orig(:,:)], 'uint8'));
imwrite(cast([recon_img(:,:), orig(:,:)], 'uint8'), 'results/q2c.png');
fprintf('RMSE : %f\n', norm(recon_img(:,:) - orig(:,:), 'fro')/norm(orig(:,:), 'fro'));

toc; % Timer end


function H = haarmtx(n)
    % HAARMTX Compute Haar orthogonal transform matrix.
    %
    %   H = HAARMTX(N) returns the N-by-N HAAR transform matrix.  H*A
    %   is the HAAR transformation of the columns of A and H'*A is the inverse
    %   transformation of the columns of A (when A is N-by-N).
    %   If A is square, the two-dimensional Haar transformation of A can be computed
    %   as H*A*H'. This computation is sometimes faster than using
    %   DCT2, especially if you are computing large number of small
    %   Haar transformation, because H needs to be determined only once.
    %
    %   Class Support
    %   -------------
    %   N is an integer scalar of class double. H is returned
    %   as a matrix of class double.
    %
    %
    %   I/O Spec
    %   N - input must be double
    %   D - output DCT transform matrix is double
    %
    %   Author : Frederic Chanal (f.j.chanal@student.tue.nl) - 2004
    a=1/sqrt(n);
    for i=1:n
        H(1,i)=a;
    end
    for k=1:n-1
        p=fix(log2(k));
        q=k-2^p+1;
        t1=n/2^p;
        sup=fix(q*t1);
        mid=fix(sup-t1/2);
        inf=fix(sup-t1);
        t2=2^(p/2)*a;
        for j=1:inf
            H(k+1,j)=0;
        end
        for j=inf+1:mid
            H(k+1,j)=t2;
        end
        for j=mid+1:sup
            H(k+1,j)=-t2;
        end
        for j=sup+1:n
            H(k+1,j)=0;
        end
    end
end