function [x_new, P_new] = merge(Gaussian_mix)
    x_new = sum( Gaussian_mix.weight .* Gaussian_mix.mean , 2);
    P_new = zeros(size(Gaussian_mix.covariance(:,:,1)));
    %P_new = sum( Gaussian_mix.weight .* Gaussian_mix.covariance , 3);
    for i = 1:length(Gaussian_mix.weight)
        P_new =P_new + Gaussian_mix.weight(i) .* Gaussian_mix.covariance(:,:,i) + Gaussian_mix.weight(i) * (x_new - Gaussian_mix.mean(:,i)) * (x_new - Gaussian_mix.mean(:,i))';
    end
end