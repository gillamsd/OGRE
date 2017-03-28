function gausspdf = gausspdf(x, mu, sigma)
% Gaussian Mixture Model
% By Andrew L. Hanna, Nov 3, 2006 at www.allbookez.com/pdf/494jto/
%
% Inputs:
%     x = data vector
%     mu = A-priori vector of means
%     sigma = A-priori covriance matrix
%     mu, sigma resprecnt are the parameters of a multiple gaussian model (a mixture of gaussians)
%
% Outputs:
%     gausspdf = a vector containing the gaussain mixture probablity density function
%
%
dim = size(x, 1);
d = length(mu);
c = (2*pi)^(d/2);
c = c*(det(sigma)^(1/2));

for i = 1:size(x, 2)
      gausspdf(i, :) = exp((-1/2)*((x(:, i)-mu)'*inv(sigma)*(x(:, i)-mu)));
end
gausspdf = gausspdf ./c;
return
