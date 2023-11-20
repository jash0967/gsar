function [f, grad] = baseline(beta,y,X,G)
% Compute the loglikelihood value using the nested fixed point
% Each input is cell

L = 0;

grad = zeros(1,size(beta,2));

for g=1:G
    X_g = X{g};
    y_g = y{g};

    u0 = exp(X_g*beta');
    F = u0./(1+u0);
    F(F==1)=1-1e-10; % If u0 is too big, F can be exact unity, which leads to a log of zero in the following.

    L_g = -(y_g.*log(F)+(1-y_g).*log(1-F));
    L = L+sum(L_g);

    % Obtain a gradient group-wise using the vector form and summurize over
    % individuals to get a 1X(Q+K) gradient vector.

    eta = X_g*beta';
    l_deriv_beta = y_g .* X_g - exp(eta)./(1+exp(eta)) .* X_g;
    grad_group = l_deriv_beta;
    
    % Due to the negative likelihood, the sum is subtracted.
    grad = grad - sum(grad_group,1);

end

f = L;

end

