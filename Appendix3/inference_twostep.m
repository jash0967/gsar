function [stars, stderror] = inference_twostep(lam,beta,y,X,W,G)
%INFERENCE Summary of this function goes here
%   Detailed explanation goes here

grad = information_mat_twostep(lam,beta, y, X, W, G);

cov_matrix = inv(grad);
intervals_HO1 = zeros(1,2);
intervals_HO2 = zeros(1,2);
intervals_HO3 = zeros(1,2);
intervals_HO1(:,1) = -sqrt(diag(cov_matrix)) *1.64 + lam;
intervals_HO1(:,2) = sqrt(diag(cov_matrix)) *1.64+ lam;
intervals_HO2(:,1) = -sqrt(diag(cov_matrix)) *1.96 + lam;
intervals_HO2(:,2) = sqrt(diag(cov_matrix)) *1.96+ lam;
intervals_HO3(:,1) = -sqrt(diag(cov_matrix)) *2.57 + lam;
intervals_HO3(:,2) = sqrt(diag(cov_matrix)) *2.57+ lam;

stderror = sqrt(diag(cov_matrix));

stars = zeros(1,1);
stars = string(stars);
for i=1:size(stars)
    if 0 < intervals_HO1(i,1) || 0 > intervals_HO1(i,2)
        stars(i) = "*";
    end

    if 0 < intervals_HO2(i,1) || 0 > intervals_HO2(i,2)
        stars(i) = "**";
    end

    if 0 < intervals_HO3(i,1) || 0 > intervals_HO3(i,2)
        stars(i) = "***";
    end

    if 0 > intervals_HO1(i,1) && 0 < intervals_HO1(i,2)
        stars(i) = "";
    end

end

end


function mat = information_mat_twostep(lam,beta,y,X,W,G)
% Compute the loglikelihood value using the nested fixed point
% Each input is cell

% Information matrix for the SAR model
mat = 0;

for g=1:G
    W_g = W{g};
    X_g = X{g};
    y_g = y{g};

    N = size(X_g,1);

    W_aggregated=lam.*W_g;

    % Obtain a gradient group-wise using the vector form and summarize over
    % individuals to get a 1X(Q+K) gradient vector.

    grad_group = zeros(N,1);

    pstar = fxp_p_alt(beta, X_g, W_aggregated);
    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2);
    phi_diag = diag(phi);

    p_deriv_rho = (eye(N) - phi_diag * W_aggregated) \ phi_diag *W_g * pstar;
    l_deriv_rho = y_g .* (W_aggregated*p_deriv_rho + W_g * pstar) - exp(eta)./(1+exp(eta)) .* (W_aggregated*p_deriv_rho + W_g * pstar);
    grad_group(:,1) = l_deriv_rho;

    mat = mat + grad_group' * grad_group;
end

end
