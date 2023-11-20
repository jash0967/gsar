function [f, grad] = twostep_obj(lam,beta,y,X,W,G)
% Compute the loglikelihood value using the nested fixed point
% Each input is cell

L = 0;

grad = zeros(1,1);

for g=1:G
    W_g = W{g};
    X_g = X{g};
    y_g = y{g};

    N = size(X_g,1);

    W_aggregated=lam.*W_g;

    pstar = fxp_p_alt(beta, X_g, W_aggregated);

    u0 = exp(W_aggregated*pstar + X_g*beta');
    F = u0./(1+u0);
    F(F==1)=1-1e-10; % If u0 is too big, F can be exact unity, which leads to a log of zero in the following.

    L_g = -(y_g.*log(F)+(1-y_g).*log(1-F));
    L = L+sum(L_g);

    % Obtain a gradient group-wise using the vector form and summarize over
    % individuals to get a 1X(Q+K) gradient vector.

    grad_group = zeros(N,1);

    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2);
    phi_diag = diag(phi);

    p_deriv_rho = (eye(N) - phi_diag * W_aggregated) \ phi_diag *W_g * pstar;
    l_deriv_rho = y_g .* (W_aggregated*p_deriv_rho + W_g * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_rho + W_g * pstar);
    grad_group(:,1) = l_deriv_rho;

    % Due to the negative likelihood, the sum is subtracted.
    grad = grad - sum(grad_group,1);

end

f = L;

end
