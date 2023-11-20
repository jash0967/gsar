function [f, grad] = nfxp_flex(coef,y,X,W,Q,G)
% Compute the loglikelihood value using the nested fixed point
% Each input is cell

lam = coef(1:Q);
sig = coef(Q+1);
beta = coef(Q+2:end);

L = 0;

grad = zeros(1,size(coef,2));

for g=1:G
    W_g = W{g};
    X_g = X{g};
    y_g = y{g};

    N = size(X_g,1);

    W_aggregated=zeros(N);
    D = zeros(N,N,Q);

    if Q>1
        for q=1:Q
            degree = diag(sum(W_g(:,:,q)));
            D(:,:,q) = degree;
            W_aggregated=W_aggregated+lam(q).*(W_g(:,:,q)-sig.*degree);
        end

    elseif Q==1
        degree = diag(sum(W_g));
        D(:,:,1) = degree;
        W_aggregated=lam(1).*(W_g-sig.*degree);

    end

    pstar = fxp_p_alt(beta, X_g, W_aggregated);

    u0 = exp(W_aggregated*pstar + X_g*beta');
    F = u0./(1+u0);
    F(F==1)=1-1e-10; % If u0 is too big, F can be exact unity, which leads to a log of zero in the following.

    L_g = -(y_g.*log(F)+(1-y_g).*log(1-F));
    L = L+sum(L_g);

    % Obtain a gradient group-wise using the vector form and summarize over
    % individuals to get a 1X(Q+K) gradient vector.

    grad_group = zeros(N,size(coef,2));

    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2);
    phi_diag = diag(phi);
    p_deriv_beta = (eye(N) - phi_diag * W_aggregated) \ phi_diag * X_g;

    l_deriv_beta = y_g .* (W_aggregated*p_deriv_beta + X_g) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_beta + X_g);
    grad_group(:,Q+2:end) = l_deriv_beta;

    sig_grad_temp = zeros(N,1);
    D_sum = zeros(N,N);

    % for q=1:Q
    %     p_deriv_rho = (eye(N) - phi_diag * W_aggregated) \ phi_diag * W_g(:,:,q) * pstar;
    %     l_deriv_rho = y_g .* (W_aggregated*p_deriv_rho + (W_g(:,:,q)-sig.*D(:,:,q)) * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_rho + (W_g(:,:,q)-sig.*D(:,:,q)) * pstar);
    %     grad_group(:,q) = l_deriv_rho;
    %
    %     p_deriv_sig = (eye(N) - phi_diag * W_aggregated) \ phi_diag * lam(q) * (-D(:,:,q)) * pstar;
    %     sig_grad_temp = sig_grad_temp + p_deriv_sig;
    %     D_sum = D_sum + lam(q) .* D(:,:,q);
    % end
    %
    % grad_group(:,Q+1) = y_g .* (W_aggregated*sig_grad_temp - D_sum * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*sig_grad_temp - D_sum * pstar);
    %
    % % Due to the negative likelihood, the sum is subtracted.
    % grad = grad - sum(grad_group,1);

    for q=1:Q
        p_deriv_rho = (eye(N) - phi_diag * W_aggregated) \ phi_diag * W_g(:,:,q) * pstar;
        l_deriv_rho = y_g .* (W_aggregated*p_deriv_rho + (W_g(:,:,q)-sig.*D(:,:,q)) * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_rho + (W_g(:,:,q)-sig.*D(:,:,q)) * pstar);
        grad_group(:,q) = l_deriv_rho;

        p_deriv_sig = (eye(N) - phi_diag * W_aggregated) \ phi_diag * lam(q) * (-D(:,:,q)) * pstar;
        l_deriv_sig = y_g .* (W_aggregated*p_deriv_sig - lam(q) .* D(:,:,q) * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_sig - lam(q) .* D(:,:,q) * pstar);
        grad_group(:,Q+1) = grad_group(:,Q+1) + l_deriv_sig;
    end

    % Due to the negative likelihood, the sum is subtracted.
    grad = grad - sum(grad_group,1);

end

f = L;

end

