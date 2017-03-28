function Z = plotGMM_1D(x, mixtures, p_i)
%
% Plots 1D gaussian mixture probablilty density function
%
Z = zeros(length(x), 1);

h = length(mixtures);

for i =1:h
   gm = mixtures(i);
   z = gausspdf(x, gm.mu, gm.sigma)*p_i(i);
   Z = Z + z;
end


