function theta = ista(y, A, lambda, alpha, iter)
    theta = zeros(size(A, 2), 1);
    thres = lambda/(2*alpha);
    for i=1:iter
        theta = soft(theta + (A'*(y - A*theta))/alpha, thres);
    end
end

function y = soft(x,T)
    y = sign(x).*(max(0, abs(x)-T));
end
