function transformed_mags = magtran_p(residuals_file, positions_list, cma_in, cma_out)
%
% Transform a set of magnitudes in an observed system to a reference syste.
% Get the observed and reference positions from the geotran_p.m residuals file.
% Get the observed magnitudes from the Positions file used by geotran_p.m
% Get the reference magnitudes from the CMA file.
%
% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.
%
% Load the OctaveForge Optimization package
% This assumes that it is installed.
pkg load optim
%
% Load CMA file
cma = load(cma_in);
[cmar, cmac] = size(cma);
  
% Load Already measured positions
% Note: residuals = [ref_coords', obs_coords', f', ref_coords'.-f'];% Save the residuals 
if exist(residuals_file, 'file')
  residuals = load(residuals_file);
else
  residuals = [];
end
[resrows, rescols] = size(residuals);

 % Read Sandquist and Hess File
%[RA, Dec, type, m1, sig_m1, m2, sig_m2] = textread(positions_list, '%10.6f %9.6f  %s %6.3f %6.3f %6.3f %6.3f');

% Read Stetson File
[RA, Dec, type, m1, sig_m1, m2, sig_m2] = textread(positions_list, '%16.12f %15.12f  %s %6.3f %6.3f %6.3f %6.3f');

sandh = [RA, Dec];
nsandh = length(RA);


sandh = [RA, Dec];
nsandh = length(RA);

% Find the CMA and positions list entries that match the ones in the residuals file
slindices = [];
clindices = [];
for i = 1: resrows
  delta = cma(:, 1:2).- repmat(residuals(i, 1:2), cmar, 1);
  rad = (delta(:, 1).^2 .+ delta(:, 2).^2).^5;
  [cl, clindex] = min(rad);
  
  deltas = sandh(:, 1:2).- repmat(residuals(i, 3:4),nsandh, 1);
  rads = (deltas(:, 1).^2 .+ deltas(:, 2).^2).^5;
  [sl, slindex] = min(rads);
  
  slindices = [slindices, slindex];
  clindices = [clindices, clindex];
end

szclindex = length(clindices);
szslindex = length(slindices);

obs_mags1 = m1(slindices);
obs_mags2 = m2(slindices);
sig_obs1 = sig_m1(slindices);
sig_obs2 = sig_m2(slindices);
obs_clrs = obs_mags1 - obs_mags2;
sig_obs_clrs = (sig_obs1.^2 .+ sig_obs2.^2).^0.5;
ref_mags2 = cma(clindices, 4);
ref_clrs = cma(clindices, 3);
ref_mags1 = ref_mags2 + ref_clrs;
mag_diffs1 = ref_mags1 - obs_mags1;
mag_diffs2 = ref_mags2 - obs_mags2;
ndiffs = length(mag_diffs1);

% Quality Information
ref_flags = cma(clindices, 5:8);
ref_CIs = cma(clindices, 10:13);
empty_col  = zeros(ndiffs, 1);

% Plot labels & fontsize
xlabel_string = "F606W-F814W";
ylabel_string = "I- F814W";
fsize = 14;

% Plots
h = figure(33);
% Set toolkit to use for this figure (= handle)
graphics_toolkit(h, 'fltk')
% Clear this plot
clf(h)
fig = gcf ();
ax = gca();
set (fig, 'visible', 'on')
set(ax, 'position', [0.12,0.15,0.75,0.78]);

% For Screen viewing
set(gcf,'position', [200,400,950,950]);
set(ax, 'xminortick', 'on');
set(ax, 'yminortick', 'on');

%plot(ref_clrs, ref_clrs-obs_clrs, "k+", "markersize", 10.0)
plot(ref_clrs, obs_mags2-ref_mags2, "k+", "markersize", 10.0)
hold on

%limits= [min(ref_clrs), max(ref_clrs), min(ref_clrs-obs_clrs), max(ref_clrs-obs_clrs)];
limits= [min(ref_clrs), max(ref_clrs), min(obs_mags2-ref_mags2), max(obs_mags2-ref_mags2)];

axis (limits);
%xlabel("WFC3 Colors", 'fontsize', 12);
%ylabel("WFPC3 Colors - Sandquist & Hess Colors", 'fontsize', 12);
xlabel(xlabel_string, 'fontsize', fsize);
ylabel(ylabel_string, 'fontsize', fsize);
legend(cstrcat("N = ", int2str(ndiffs)))
set(gca,'fontsize', fsize); % sets font of numbers on axes

select = 1;
pindices = [];
while (select == 1)

  % Select data ranges to keep.
  figure(33);
  [x, y, button] = ginput(2);
  x1  =  min(x);
  x2  =  max(x);
  y1  =  min(y);
  y2  =  max(y);

  % Cull the outliers
 % pindex = find((ref_clrs >= x1) & (ref_clrs <= x2) &...
 %               (obs_clrs >= y1) & (ref_clrs-obs_clrs <= y2 ));
 
  pindex = find((ref_clrs >= x1) & (ref_clrs <= x2) &...
                (obs_mags2 .- ref_mags2 >= y1) & (obs_mags2 .- ref_mags2 <= y2 ));
  [br, bc] = size(obs_mags1);

  % Mark rejected obs.
  plot(ref_clrs(pindex), obs_mags2(pindex)-ref_mags2(pindex), "ro", "markersize", 10.0)
  %plot(ref_clrs(pindex), ref_clrs(pindex)-obs_clrs(pindex), "ro", "markersize", 10.0)
  
  select = yes_or_no("Continue removing observations? ");
  
  % Accumulat the list of culled obs. for later.
  pindices = [pindices; pindex];
end
npind = length(pindices);

obsref1 = [obs_mags1, ref_mags1, obs_clrs, ref_clrs, sig_obs_clrs, ref_flags, empty_col, ref_CIs];
obsref2 = [obs_mags2, ref_mags2, obs_clrs, ref_clrs, sig_obs_clrs, ref_flags, empty_col, ref_CIs];

% Keep these
if npind > 0
  obsref1 = setdiff(obsref1, obsref1(pindices, :), "rows");
  obsref2 = setdiff(obsref2, obsref2(pindices, :), "rows");
end

obs_mags1 = obsref1(:, 1);
ref_mags1 = obsref1(:, 2);
obs_mags2 = obsref2(:, 1);
ref_mags2 = obsref2(:, 2);

% If using first S&H magnitude
%obs_clrs = obsref1(:, 3);
%ref_clrs = obsref1(:, 4);
%sig_obs_clrs = obsref1(:, 5);
%ref_flags = obsref1(:, 6:9);
%empty_col = obsref1(:, 10);
%ref_CIs = obsref1(:, 11:14);
% If using second Stetson magnitude
obs_clrs = obsref2(:, 3);
ref_clrs = obsref2(:, 4);
sig_obs_clrs = obsref2(:, 5);
ref_flags = obsref2(:, 6:9);
empty_col = obsref2(:, 10);
ref_CIs = obsref2(:, 11:14);

[ar, ac] = size(obs_mags1);
num_culled = br - ar;
disp(cstrcat(int2str(num_culled), " observations discarded"));
  
mag_diffs1 = ref_mags1 - obs_mags1;
mag_diffs2 = ref_mags2 - obs_mags2;

clr_diffs = ref_clrs - obs_clrs;
ndiffs = length(mag_diffs2);


clc
% Initial guesses for first mag.
%c0_init1 = mean(mag_diffs1);
%sig_c0_init1 = std(mag_diffs1);
%model_params1 = [c0_init1, 0, 0];

% Initial guesses for second mag.
c0_init1 = mean(mag_diffs2);
sig_c0_init1 = std(mag_diffs2);
model_params1 = [c0_init1, 0, 0];

% Compute transformation coeeficinets 
sig_obs_clrs = 1.*sig_obs_clrs;
%data_wts = 1./sig_obs_clrs;
%data_wts = ones(ndiffs, 1);
data_wts = 1./abs(obs_clrs).^.5;

% For the first Sandquist magnitude (m1)
%[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]= leasqr(obs_clrs, mag_diffs1, model_params1, ...
%                            'Magnitude_Transformation_Model_p', 1e-9, 10000, data_wts);

%obs_trans1 = obs_mags1 + f;
%mag_residuals = ref_mags1 .- obs_trans1; 
%RMS_mag_resids = (mean(mag_residuals.^2)).^.5
%nmresids = rows(mag_residuals);
%nmbins = (max(mag_residuals) - min(mag_residuals))/(0.3*std(mag_residuals));
                            
                            
% For the second setson magnitudes (m2)                           
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]= leasqr(obs_clrs, mag_diffs2, model_params1, ...
                            'Magnitude_Transformation_Model_p', 1e-9, 10000, data_wts);                           
                                              
obs_trans2 = obs_mags2 + f;
mag_residuals = ref_mags2 .- obs_trans2; 
RMS_mag_resids = (mean(mag_residuals.^2)).^.5
nmresids = rows(mag_residuals);
nmbins = (max(mag_residuals) - min(mag_residuals))/(0.3*std(mag_residuals));


hold off
% Mark final selection

% Plots
h = figure(33);
% Set toolkit to use for this figure (= handle)
graphics_toolkit(h, 'fltk')
% Clear this plot
clf(h)
fig = gcf ();
ax = gca();
set (fig, 'visible', 'on')
set(ax, 'position', [0.12,0.15,0.75,0.78]);

% For Screen viewing
set(gcf,'position', [200,400,950,950]);
set(ax, 'xminortick', 'on');
set(ax, 'yminortick', 'on');

plot(ref_clrs, obs_mags2-ref_mags2, "b+", "markersize", 10.0)

legend(cstrcat("N = ", int2str(ndiffs)))
%limits= [min(ref_clrs), max(ref_clrs), min(ref_clrs-obs_clrs), max(ref_clrs-obs_clrs)];
limits= [min(ref_clrs), max(ref_clrs), min(obs_mags2-ref_mags2), max(obs_mags2-ref_mags2)];
axis (limits);
%xlabel("WFC3 Instrumental Colors", 'fontsize', 12);
%ylabel("Comparison - WFC3 Magnitudes", 'fontsize', 12);
xlabel(xlabel_string, 'fontsize', fsize);
ylabel(ylabel_string, 'fontsize', fsize);
legend(cstrcat("N = ", int2str(ndiffs), "  RMS=", num2str(RMS_mag_resids, 1)))
set(gca,'fontsize', fsize); % sets font of numbers on axes

% Display Results

c0_1 = p(1);
c1_1 = p(2);
c2_1 = p(3);

sig_c0_1 = covp(1,1)^.5;
sig_c1_1 = covp(2,2)^.5;
sig_c2_1 = covp(3,3)^.5;

% Compute transformation coefficients color (m1-m2) --> WFC3 color
% Initial guesses
c0_init2 = mean(mag_diffs1 .- mag_diffs2);
sig_c0_init2 = std(c0_init2);
model_params2 = [c0_init2, 0, 0];
[f,p,kvg,iter,corp,covp,covr,stresid,Z,r2]= leasqr(obs_clrs, clr_diffs, model_params2, ...
                            'Magnitude_Transformation_Model_p', 1e-9, 10000, data_wts);
                                              
% If Transformed mag1 is estimated - construct transformed mag2
%obs_trans2 = obs_mags2 + f;
%mag_residuals2 = ref_mags2 .- obs_trans2; 
%chi_trans_mags2 = (mean(mag_residuals2.^2)).^.5
%nmresids2 = rows(mag_residuals2);
%nmbins = (max(mag_residuals) - min(mag_residuals))/(0.3*std(mag_residuals));

% If Transformed mag2 is estimated - cocstruct transformed mag1
obs_trans1 = obs_mags2 - f;
mag_residuals1 = ref_mags1 .- obs_trans1; 
RMS_color_resids = (mean(mag_residuals1.^2)).^.5
nmresids1 = rows(mag_residuals1);
nmbins = (max(mag_residuals) - min(mag_residuals))/(0.3*std(mag_residuals));



% Display Results

c0_2 = p(1);
c1_2 = p(2);
c2_2 = p(3);

sig_c0_2 = covp(1,1)^.5;
sig_c1_2 = covp(2,2)^.5;
sig_c2_2 = covp(3,3)^.5;

disp(" Transformation: TMAG = SMAG + c0 + c1 TCOL + c2 TCOL^2");
disp(" Where, TMAG = magnitude in WFP3 system")
disp("        SMAG = magnitude in WFPC2 system")
disp("        TCOL = color in WFPC2 system")
disp("  ")
disp ("Initial model for first filter magnitude")
disp(cstrcat("c0_1 = ", num2str(model_params1(1), "%7.4f")))
disp(cstrcat("c1_1 = ", num2str(model_params1(2), "%7.4f")))
disp(cstrcat("c2_1 = ", num2str(model_params1(3), "%7.4f")))
disp("  ")
disp ("Solution for first filter magnitude")
disp(cstrcat("c0_1 = ", num2str(c0_1, "%7.4f"), " +/- ", num2str(sig_c0_1, "%7.4f")))
disp(cstrcat("c1_1 = ", num2str(c1_1, "%7.4f"), " +/- ", num2str(sig_c1_1, "%7.4f")))
disp(cstrcat("c2_1 = ", num2str(c2_1, "%7.4f"), " +/- ", num2str(sig_c2_1, "%7.4f")))
disp("  ")
disp ("Initial model for color transformation")
disp(cstrcat("c0_2 = ", num2str(model_params2(1), "%7.4f")))
disp(cstrcat("c1_2 = ", num2str(model_params2(2), "%7.4f")))
disp(cstrcat("c2_2 = ", num2str(model_params2(3), "%7.4f")))
disp("  ")
disp ("Solution for color transformation")
disp(cstrcat("c0_2 = ", num2str(c0_2, "%7.4f"), " +/- ", num2str(sig_c0_2, "%7.4f")))
disp(cstrcat("c1_2 = ", num2str(c1_2, "%7.4f"), " +/- ", num2str(sig_c1_2, "%7.4f")))
disp(cstrcat("c2_2 = ", num2str(c2_2, "%7.4f"), " +/- ", num2str(sig_c2_2, "%7.4f")))
disp("  ")
disp ("Solution for second filter magnitude")
sig0 = (sig_c0_1^2 + sig_c0_2^2)^0.5;
sig1 = (sig_c1_1^2 + sig_c1_2^2)^0.5;
sig2 = (sig_c2_1^2 + sig_c2_2^2)^0.5;
disp(cstrcat("k0 = ", num2str(c0_1-c0_2, "%7.4f"), " +/- ", num2str(sig0, "%7.4f")))
disp(cstrcat("k1 = ", num2str(c1_1-c1_2, "%7.4f"), " +/- ", num2str(sig1, "%7.4f")))
disp(cstrcat("k2 = ", num2str(c2_1-c2_2, "%7.4f"), " +/- ", num2str(sig2, "%7.4f")))



% Plot results
tmags = obs_mags1 .+ c0_1 .+ (c1_1.* obs_clrs) .+ (c2_1 .* obs_clrs.^2);
tcols = obs_clrs .+ c0_2 .+ (c1_2.* obs_clrs) .+ (c2_2 .* obs_clrs.^2);
tmags2 = tmags .- tcols;
%plot (obs_mags1, tmags, "mo")

% Save residuals with rejected observations removed.
if exist('tempm.resids', 'file')
    delete('tempm.resids')
end

% Keep the residuals corresponding to the obs. selscted from the mag. transformation 
residuals = setdiff(residuals, residuals(pindices, :), "rows");
dlmwrite ('tempm.resids', residuals,  "precision", "%10.7f", "delimiter", " ");
[status, msg, msgid] = copyfile ('tempm.resids', residuals_file);

delete('tempm.resids');

% Make up a CMA file of the cross-identified stars
if exist('tempm.cma', 'file')
    delete('tempm.cma')
end
cma_array = [residuals(:, 1:2), ref_clrs, ref_mags1, ref_flags, empty_col, ref_CIs]; 
[r,c] = size(cma_array)
dlmwrite ('tempm.cma', cma_array,  "precision", "%10.7f", "delimiter", " ");
[status, msg, msgid] = copyfile ('tempm.cma', cma_out);
delete('tempm.cma');select

% Display magnitude residuals
figure(34);
fig = gcf ();
ax = gca();
hist(mag_residuals,nmbins, "facecolor", "g")
title (cstrcat('Magnitude Residuals'), 'fontsize', 12)
xlabel("Observed Magnitude - Estimated Magnitude", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axe s
legend(cstrcat("N = ", int2str(nmresids)))  


% Return to starting directory
% Unload Optimization package
pkg unload optim
cd(currentdir);
