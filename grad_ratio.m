function theta = grad_ratio(x)
% x: N x 1 vector
% theta: size = size(x) - 2
% calculate theta = grad_prev/grad_next
d = diff(x);

theta = d(1:end-1)./d(2:end);

theta(isnan(theta) | isinf(theta)) = 0;
end

