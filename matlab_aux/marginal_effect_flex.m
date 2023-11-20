function [direct, indirect, naive] = marginal_effect_flex(coef,y,X,W,Q,G)

lam = coef(1:Q);
sig = coef(Q+1);
beta = coef(Q+2:end);

K=size(beta,2);

naive_k=zeros(1,K);
naive_q=zeros(1,Q);
direct=zeros(1,K);
indirect=zeros(1,K);

population = 0;

disp('Evaluating marginal effects...')
for g=1:G
    disp(['Village no.',num2str(g)])

    W_g = W{g};
    X_g = X{g};
    y_g = y{g};

    N = size(X_g,1);
    population = population + N;

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

    naive_temp_k=zeros(N,K);
    naive_temp_q=zeros(N,Q);
    direct_temp=zeros(N,K);
    indirect_temp=zeros(N,K);

    pstar = fxp_p_alt(beta, X_g, W_aggregated);
    eta = W_aggregated*pstar + X_g*beta';
    phi = exp(-eta)./ ((1+exp(-eta)).^2); % P[y=1] * (1-P[y=1])
    phi_diag = diag(phi);

    % For the equations, see Lee et al.(2014) p.411, Chomsisengphet et al.(2018) p.75
    % The "naive" marginal effect
    parfor i=1:N
        for k=1:K
            if k==2||k==3
                naive_temp_k(i,k) = phi_diag(i,i) .* beta(k);
            else
                cov_temp_1 = X_g(i,:);
                cov_temp_0 = X_g(i,:);
                cov_temp_1(k) = 1;
                cov_temp_0(k) = 0;

                eta_temp_1 = W_aggregated(i,:) * pstar + cov_temp_1*beta' ;
                eta_temp_0 = W_aggregated(i,:) * pstar + cov_temp_0*beta' ;

                p1 = exp(eta_temp_1)./(1+exp(eta_temp_1));
                p0 = exp(eta_temp_0)./(1+exp(eta_temp_0));

                naive_temp_k(i,k) = p1-p0;
            end
        end

        for q=1:Q
            naive_temp_q(i,q) = phi_diag(i,i) .* lam(q);
        end

    end
    naive_k(1,:) = naive_k(1,:) + sum(naive_temp_k,1);
    naive_q(1,:) = naive_q(1,:) + sum(naive_temp_q,1);

    % The "sophisticated" marginal effect
    parfor i=1:N
        parallel_temp=zeros(N,K);
        for k=1:K
            if k==2||k==3
                deriv_X=zeros(N,K);
                deriv_X(i,k)=1;

                M_partial = (eye(N) - phi_diag * W_aggregated) \ phi_diag * deriv_X * beta';
                direct_k = phi_diag(i,i) * (beta(k) + W_aggregated(i,:) * M_partial);

                for j=1:N
                    if j~=i
                        indirect_k = phi_diag(j,j) * (W_aggregated(j,:) * M_partial);
                        parallel_temp(j,k) = indirect_k;
                    end
                end

            else
                cov_temp_1 = X_g;
                cov_temp_0 = X_g;
                cov_temp_1(i,k) = 1;
                cov_temp_0(i,k) = 0;

                pstar_1 = fxp_p_alt(beta, cov_temp_1, W_aggregated);
                pstar_0 = fxp_p_alt(beta, cov_temp_0, W_aggregated);

                direct_k = pstar_1(i)-pstar_0(i);

                for j=1:N
                    if j~=i
                        indirect_k = pstar_1(j)-pstar_0(j);
                        parallel_temp(j,k) = indirect_k;
                    end
                end
            end
            direct_temp(i,k) = direct_k;
        end
        indirect_temp(i,:) = sum(parallel_temp,1)./(N-1);
    end
    direct(1,:) = direct(1,:) + sum(direct_temp,1);
    indirect(1,:) = indirect(1,:) + sum(indirect_temp,1);
end

% Take an average with the computed population
direct=direct'./population;
indirect=indirect'./population;
naive_k=naive_k./population;
naive_q=naive_q./population;
naive = [naive_k, naive_q]';
end

