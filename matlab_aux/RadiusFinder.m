function [radius_bound, radius] = RadiusFinder(network, G, Q)
%RADIUSFINDER Summary of this function goes here
%   Detailed explanation goes here

radius = zeros(1,G);

for g=1:G
    temp = network{1,g};
    sum = zeros(size(temp(:,:,1)));

    for q=1:Q

        temp_q = temp(:,:,q);
        sum = sum + temp_q;
    end

    radius_g = eigs(sum,1);
    radius(1,g) = abs(radius_g);
end

radius_bound = max(radius);

end

