
function theta = ista(y, A, lambda, alpha, iter)
    % Input:
    %   y : Signal
    %   A : Overcomplete dictionary
    %   lambda : Weight of l1 norm
    %   alpha : Step size
    %   iter : Number of iterations
    % Output:
    %   theta : Sparse coefficients
    % Brief:
    %   Iterative Shrinkage and Thresholding Algorithm (ISTA) for solving
    %    cost function J(theta) = ||y - A*theta||_2^2 + lambda * ||theta|_1

    % Initialize theta
    theta = zeros(size(A, 2), 1);

    % Set soft thershold
    thres = lambda/(2*alpha);

    % Iteratively converge
    for i=1:iter
        theta = soft(theta + (A' * (y - A*theta))/alpha, thres);
    end
end

function y = soft(x,T)
    % Input:
    %   x : Input value
    %   T : Soft thresholding value (> 0)
    % Output:
    %   y : Output value
    % Brief:
    %   if x < -T, y = x + T
    %   if -T < x < T, y = 0
    %   if x > T, y = x - T

    y = sign(x).*(max(0, abs(x)-T));
end
