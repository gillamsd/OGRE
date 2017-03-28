function [hess_array, n_cell] = hess_gen_p(data, bluetxt, redtxt, mincolor, ...
  maxcolor, minmag, maxmag, nclrs, nmags, nstar_min, ctol)
%
% Create a CMD density (Hess) array from a CM array or CM3 array
% INPUTS:
%   data = Color-magnitude or color-color array in file cma_namne.
%   bluetxt = phtometric errors for blue filter
%   redtxt = TXT file containg phtometric errors for red filter
%   mincolor and maxcolor are the lower and upper limits on colors that will be plotted
%   minmag and maxmag are the lower and upper limits on mags. that will be plotted
%   nclrs = number of color bins
%   nmags = number of magnitude bins 
%   nstar_min = minimum number of stars for statistics
%   ctol =  color cut
%
% OUTPUTS:
%   hess_array = 2-D Array containing number of stars in each color-magnitude bin
%   n_cell = number of stars in the data
%
% CMA entries will be summed over non-overlapping boxes dcolor x dmag in area.
%
% Save current working directory.

currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

% Star coordinates from  .TXT files
blue_RAs = bluetxt(:,3);
blue_Decs = bluetxt(:, 4);
red_RAs = redtxt(:,3);
red_Decs = redtxt(:, 4);

% Find the intersections of the cma and the TXT files
[cleaned, ia, ib] = intersect([blue_RAs, blue_Decs], [data(:,1), data(:,2)], 'rows');
blue_sigmas = bluetxt(ia, 9);
[blue_rows, blue_cols] = size(blue_sigmas);

[cleaned, ia, ib] = intersect([red_RAs, red_Decs], [data(:,1), data(:,2)], 'rows');
red_sigmas = redtxt(ia, 9);
[red_rows, red_cols] = size(red_sigmas);

% Add sigmas to data array
data = data(ib, :);

% Color errors
color_sigmas = sqrt(blue_sigmas .* blue_sigmas .+ red_sigmas .* red_sigmas);

% Phtotometric cut
indc = find(color_sigmas <= ctol);
data = data(indc, :);
red_sigmas = red_sigmas(indc);
blue_sigmas = blue_sigmas(indc);
color_sigmas = color_sigmas(indc);

% Number of observations
nrows = rows(data);

clr_locs = linspace(mincolor, maxcolor, nclrs + 1);
mag_locs = linspace(minmag, maxmag, nmags + 1);

n_cell = 0;
hess_array = [];

% Distribute stars' red mags over the mag bins and loactions in 
% the magnitude dimension
if nrows >= nstar_min;
  p_m = hist_fractions_p(data(:, 4), red_sigmas, mag_locs);

  % Distribute stars' colors over the color bins and location in the color
  % dimesion.
  p_c = hist_fractions_p(data(:, 3), color_sigmas, clr_locs);

%  [p_c_rows, p_c_cols] = size(p_c);
%  [p_m_rows, p_m_cols] = size(p_m);

%hess_array = repmat(p_c(1, :), nmags, 1).* repmat(p_m(1, :)', 1, nclrs);

  a = p_c(1, :);
  b = p_m(1, :)';
  a_inds = [1:nclrs].*ones(nmags,1);  
  b_inds = [1:nmags]'.*ones(1,nclrs);
  %hess_array = a([1:size(a,2)].*ones(nmags,1)).*b([1:size(b,1)]'.*ones(1,nclrs));
  hess_array = a(a_inds).*b(b_inds);

  for k = 2:nrows
    a = p_c(k, :);
    b = p_m(k, :)';
    % hess_array = hess_array .+ repmat(p_c(k, :), nmags, 1).* repmat(p_m(k, :)', 1, nclrs);
    %hess_array = hess_array .+ a([1:size(a,2)].*ones(nmags,1)).*b([1:size(b,1)]'.*ones(1,nclrs));
    hess_array = hess_array .+ a(a_inds).*b(b_inds);
  end
% Number of stars in (RA, Dec) cell
  n_cell = sum(sum(hess_array));

end

clear p_c;
clear p_m;

cd(currentdir);

return;
