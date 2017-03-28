function [old_mixtures, old_p_i, old_o, old_chisqr] = gmmem(X, mixtures, histogram, bin_centers, tol, sf)

dim = size(X, 1);
N = size(X, 2);
M = length(mixtures);
mixtures.mu;

% Data histogram fro chi-square test
histX = histogram;
xmax = max(max(X));
xmin = min(min(X));
nbins = length(histX);
dm = (xmax-xmin)/ nbins;

x_chi = bin_centers;
normZ = sum(histX);

% Initializations
p_i = ones(1, M)/M;
minEigenvalue = 2;

chisqr = 1e99;

iters =  3200;
o = 1;

% Chi-squared stopping criterion
   x_s = x_chi - dm/(2*sf);
   x_s_max = x_chi(length(x_chi)) + dm/(2*sf);
   x_s = [x_s, x_s_max];


do
   old_mixtures = mixtures;
   old_chisqr = chisqr;
   old_p_i = p_i;
   old_o = o;
   for i=1:M
      gm = mixtures(i);

      V = X - gm.mu * ones(1, N);
      Q(:, i) = gausspdf(X, gm.mu, gm.sigma) * p_i(i); % p(x_d|i)p(i)
   end

   su = sum(Q, 2);
   Q = Q./repmat(su, [1 M]);

   % log likelyhood
   logl = sum(log(su));

   for i=1:M
      gm = mixtures(i);
      sum_p_i_x = sum(Q(:, i)); %sum from d=1 to N p(i|x_d)

      p_i(i) = (1/size(Q, 1)) * sum_p_i_x;
      gm.mu = X*Q(:, i)./sum_p_i_x;

      sigma = zeros(dim);

      for j = 1:N
         V = X(:, j) - gm.mu;
         sigma = sigma + Q(j, i)*(V*V');
      end

      old_sigma = gm.sigma;

      gm.sigma = sigma/sum_p_i_x;

       [P, L] = eig(gm.sigma);
       if any(diag(L)<minEigenvalue)
       disp("Failed eigenvalue test")
          gm.sigma = old_sigma;
       end
      mixtures(i) = gm;

      % For chi-squared test
      tmixtures(i).mu = gm.mu/sf;
      tmixtures(i).sigma = gm.sigma/sf^2;
   end
   o = o + 1;

   
   
   tZ = plotGMM_1D(x_s, tmixtures, p_i); % Model histogram

   nu = nbins - 2 * M - 1; % Number of degress of freedom

   F = 0;
   for i=1:M
    for j = 1 : nbins
      xs_bin = linspace(x_s(j), x_s(j+1), 50);
      y = plotGMM_1D(xs_bin, tmixtures(i), 1);
      F_comp(j) = trapz(xs_bin', y)*p_i(i);
    end  
      %sumFcomp = sum(F_comp)
      F = F + F_comp;
   end
   
  
   %chisqr = (1/nu)*sum((histX' .- tZ*dm*(normZ/sf)).^2 ./ tZ*dm(normZ/sf));
   %chisqr = (1/nu)*sum((histX .- F*dm*(normZ/sf)).^2 ./ F*dm*(normZ/sf));
   
   d = histX .- F*normZ;
   chisqr = sum(d.*d)/N;
   diff = chisqr .- old_chisqr;
 
until ((diff > -2e-6) || (o == iters));

return;
