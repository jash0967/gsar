function output = smtp(coef,y,X,W,Q,G)

lam = coef(1:Q);
beta = coef(Q+1:end);
smtp = zeros(1,G);

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

    pstar = fxp_p_alt(beta, X_g, W_aggregated);
    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2);
    phi_diag = diag(phi);

    smtp(g) = 1/det(eye(N)-phi_diag*W_aggregated);

end
output = mean(smtp);

end

