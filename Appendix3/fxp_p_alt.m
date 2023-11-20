function p=fxp_p_alt(beta,X,W)

% This function solves the nonlinear equation and find the fixed point, p*.
% In this alt. ver., lambda is absorbed in W, which is provided as a
% combined higher-order weight matrix.

N = size(X,1);
p = zeros(N,1);

d=ones(N,1);
p1=ones(N,1);

counter = 1;

while d'*d > 1.0e-5
    p0=p1;

    u0=exp(W*p0 + X(:,:)*beta');
    p1 = u0./(1+u0);
    d=p1-p0;

    counter = counter + 1;
    if counter>1.0e+4
        break;
    end
    
end

p(:,1)=p1;
end