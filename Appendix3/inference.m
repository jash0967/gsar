function [stars, stderror] = inference(coef,y,X,W,Q,G)
%INFERENCE Summary of this function goes here
%   Detailed explanation goes here

if Q == 0
    K = size(X{1,1},2);

    grad = information_mat(coef, y, X, W, Q, G);

    cov_matrix = inv(grad);
    intervals_HO1 = zeros(Q+K,2);
    intervals_HO2 = zeros(Q+K,2);
    intervals_HO3 = zeros(Q+K,2);
    intervals_HO1(:,1) = -sqrt(diag(cov_matrix)) *1.64 + coef';
    intervals_HO1(:,2) = sqrt(diag(cov_matrix)) *1.64+ coef';
    intervals_HO2(:,1) = -sqrt(diag(cov_matrix)) *1.96 + coef';
    intervals_HO2(:,2) = sqrt(diag(cov_matrix)) *1.96+ coef';
    intervals_HO3(:,1) = -sqrt(diag(cov_matrix)) *2.57 + coef';
    intervals_HO3(:,2) = sqrt(diag(cov_matrix)) *2.57+ coef';

    stderror = sqrt(diag(cov_matrix));

    stars = zeros(Q+K,1);
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

else
    K = size(X{1,1},2);

    grad = information_mat(coef, y, X, W, Q, G);

    cov_matrix = inv(grad);
    intervals_HO1 = zeros(Q+K,2);
    intervals_HO2 = zeros(Q+K,2);
    intervals_HO3 = zeros(Q+K,2);
    intervals_HO1(:,1) = -sqrt(diag(cov_matrix)) *1.64 + coef';
    intervals_HO1(:,2) = sqrt(diag(cov_matrix)) *1.64+ coef';
    intervals_HO2(:,1) = -sqrt(diag(cov_matrix)) *1.96 + coef';
    intervals_HO2(:,2) = sqrt(diag(cov_matrix)) *1.96+ coef';
    intervals_HO3(:,1) = -sqrt(diag(cov_matrix)) *2.57 + coef';
    intervals_HO3(:,2) = sqrt(diag(cov_matrix)) *2.57+ coef';

    stderror = sqrt(diag(cov_matrix));

    stars = zeros(Q+K,1);
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

end

