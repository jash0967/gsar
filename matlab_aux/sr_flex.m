function [rad, networks] = sr_flex(coef, W, X, Q, G)
%SR Summary of this function goes here
%   Detailed explanation goes here
lam = coef(1:Q);
beta = coef(Q+2:end);
sig = coef(Q+1);

rad = zeros(G,1);
networks = cell(1,G);

for g=1:G
    W_g = W{g};
    X_g = X{g};

    N = size(W_g,1);

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
    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2);
    phi_diag = diag(phi);

    rad(g) = eigs(phi_diag*W_aggregated, 1);
    % rad(g) = eigs(W_aggregated, 1);
    networks{g} = W_aggregated;
end

end

