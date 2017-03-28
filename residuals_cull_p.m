function cleaned = residuals_cull_p(residuals_file)
%
% Remove data from a residuals file created by geotran_p.m
%
% INPUTS
%   FILE: residuals_file
%
% OUTPUTS
%   FILE: Replaces input file
%
% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.
%
% Load the reiduals file
if exist(residuals_file, 'file')
  residuals = load(residuals_file);
 else
  disp("That file does not exist");
  exit()
end
[resrows, rescols] = size(residuals);

% Plot Histogram of residuals
rng = abs(range(residuals(:, 7:8)));
mresid = std(abs(residuals(:, 7:8)));
nbins_RA = rng(:, 1)/(.3*mresid(:, 1));
nbins_dec = rng(:, 2)/(.3*mresid(:, 2));

figure(31);
fig = gcf ();
ax = gca();
hist(residuals(:,7),nbins_RA)
title (cstrcat('Right Ascemsion Residuals'), 'fontsize', 12)
xlabel("Observed RA - Estimated RA (Degrees)", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axes

figure(32);
fig = gcf ();
ax = gca();
hist(residuals(:,8),nbins_dec)
title (cstrcat('DeclinationResiduals'), 'fontsize', 12)
xlabel("Observed Dec - Estimated Dec (Degrees)", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axes

% Select data ranges to keep.
figure(31);
[x1, y1, button] = ginput(2);
figure(32);
[x2, y2, button] = ginput(2);
min_RA_res  =  min(x1);
max_RA_res  =  max(x1);
min_dec_res  =  min(x2);
max_dec_res  =  max(x2);

% Cull the outliers
RAindex = find((residuals(:,7) >= min_RA_res) & (residuals(:,7) <= max_RA_res) );
[br, bc] = size(residuals);
residuals = residuals(RAindex, :);

decindex = find((residuals(:,8) >= min_dec_res) & residuals(:,8) <= max_dec_res);
residuals = residuals(decindex, :);
[ar, ac] = size(residuals);
num_culled = br - ar;
disp(cstrcat(int2str(num_culled), " observations removed from the residuals file"));


if exist('temp.resids', 'file')
   delete('temp.resids')
end

dlmwrite ('temp.resids', residuals,  "precision", "%10.7f", "delimiter", " ");
[status, msg, msgid] = copyfile ('temp.resids', residuals_file);

delete('temp.resids');
figure(31);
fig = gcf ();
ax = gca();
hist(residuals(:,7),nbins_RA)
title (cstrcat('Right Ascemsion Residuals'), 'fontsize', 12)
xlabel("Observed RA - Estimated RA (Degrees)", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axes

figure(32);
fig = gcf ();
ax = gca();
hist(residuals(:,8),nbins_dec)
title (cstrcat('DeclinationResiduals'), 'fontsize', 12)
xlabel("Observed Dec - Estimated Dec (Degrees)", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axes

cd(currentdir);