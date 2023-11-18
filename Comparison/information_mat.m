function mat = information_mat(coef,y,X,W,Q,G)
% Compute the loglikelihood value using the nested fixed point
% Each input is cell

if Q==0
    % Information matrix for the standard binary choice model
    mat = zeros(1,size(coef,2));

    for g=1:G
        X_g = X{g};
        y_g = y{g};

        % Obtain a gradient group-wise using the vector form and summarize over
        % individuals to get a 1X(Q+K) gradient vector.

        eta = X_g*coef';
        l_deriv_beta = y_g .* X_g - exp(eta)./(1+exp(eta)) .* X_g;
        grad_group = l_deriv_beta;

        mat = mat + grad_group' * grad_group;
    end

else
    % Information matrix for the SAR model
    lam = coef(1:Q);
    beta = coef(Q+1:end);

    mat = zeros(1,size(coef,2));

    for g=1:G
        W_g = W{g};
        X_g = X{g};
        y_g = y{g};

        N = size(X_g,1);

        W_aggregated=zeros(N);
        if Q>1
            for q=1:Q
                W_aggregated=W_aggregated+lam(q).*W_g(:,:,q);
            end
        end
        if Q==1
            W_aggregated=lam(1).*W_g;
        end

        % Obtain a gradient group-wise using the vector form and summarize over
        % individuals to get a 1X(Q+K) gradient vector.

        grad_group = zeros(N,size(coef,2));

        pstar = fxp_p_alt(beta, X_g, W_aggregated);
        eta = W_aggregated*pstar + X_g*beta';
        phi = exp(-eta)./ ((1+exp(-eta)).^2);
        phi_diag = diag(phi);
        p_deriv_beta = (eye(N) - phi_diag * W_aggregated) \ phi_diag * X_g;

        l_deriv_beta = y_g .* (W_aggregated*p_deriv_beta + X_g) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_beta + X_g);
        grad_group(:,Q+1:end) = l_deriv_beta;

        for q=1:Q
            p_deriv_rho = (eye(N) - phi_diag * W_aggregated) \ phi_diag *W_g(:,:,q) * pstar;
            l_deriv_rho = y_g .* (W_aggregated*p_deriv_rho + W_g(:,:,q) * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_rho + W_g(:,:,q) * pstar);
            grad_group(:,q) = l_deriv_rho;
        end
        mat = mat + grad_group' * grad_group;
    end
end

end

