function geotran_p(cma_in, positions_list, dRA, ddec, RAl, RAu, decl, decu)
#
% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.
%
% Load the OctaveForge Optimization package
% This assumes that it is installed.
pkg load optim
% 
% Load command parameters
param_file = "Param_Files\\geotran.par";
if exist(param_file, 'file')
   [params, npar] = read_param_file(param_file);

   cma_in = params.CMA_FILE;
   positions_list = params.POSITIONS_LIST;
   residuals_file = params.RESIDUALS_FILE;
   stellar_type = params.STELLAR_TYPE;
   model = params.MODEL;
   footprint1 = params.FOOTPRINT1;
   footprint2 = params.FOOTPRINT2;
   footprint3 = params.FOOTPRINT3;
   footprint4 = params.FOOTPRINT4;
   footprint5 = params.FOOTPRINT5;
   footprint6 = params.FOOTPRINT6;
   cma_cl = str2num(params.CMA_COLOR_LOWER_LIMIT);
   cma_cu = str2num(params.CMA_COLOR_UPPER_LIMIT);
   cma_ml = str2num(params.CMA_MAGNITUDE_LOWER_LIMIT);
   cma_mu = str2num(params.CMA_MAGNITUDE_UPPER_LIMIT);
   snh_cl = str2num(params.SNH_COLOR_LOWER_LIMIT);
   snh_cu = str2num(params.SNH_COLOR_UPPER_LIMIT);
   snh_ml = str2num(params.SNH_MAGNITUDE_LOWER_LIMIT);
   snh_mu = str2num(params.SNH_MAGNITUDE_UPPER_LIMIT);
   dRA = str2num(params.RA_SHIFT_TO_LIST_POSNS);
   ddec = str2num(params.DEC_SHIFT_TO_LIST_POSNS);
   RAl = str2num(params.DISPLAY_MIN_RA);
   RAu = str2num(params.DISPLAY_MAX_RA);
   decl = str2num(params.DISPLAY_MIN_DEC);
   decu = str2num(params.DISPLAY_MAX_DEC);
   use_external_model = params.USE_MODEL_PARAMS
   shift_RA = str2num(params.RA_OFFSET);
   shift_dec = str2num(params.DEC_OFFSET);
   theta0 = str2num(params.THETA0);
   RA_scale = str2num(params.RA_SCALE_FACTOR);
   dec_scale = str2num(params.DEC_SCALE_FACTOR);
   non_par_factor = str2num(params.NON_PARALLEL_FACTOR);
   c0 = str2num(params.TMAG_C0);
   c1 = str2num(params.TMAG_C1);
   c2 = str2num(params.TMAG_C2);
else
   disp ("No parameter file available. Using command line arguments")
   npar = 0;
end

% Load CMA file
cma = load(cma_in);
%sandh = load(positions_list);
[cmar, cmac] = size(cma);
 
% Load Already measured positions
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
manual_shift = [dRA, ddec];
snh_clr = m1 - m2;

% Defaults
[co_min_def] = min(cma)(1:2);
[co_max_def] = max(cma)(1:2);

if  nargin < 8 && npar < 8
   RAl = co_min_def(1);
   decl = co_min_def(2);
   RAu = co_max_def(1);
   decu = co_max_def(2);
end

% Apply selection criteria to data
index = find(cma(:, 3) >= cma_cl & cma(:, 3) <= cma_cu & cma(:, 4) >= cma_ml & cma(:, 4) <= cma_mu & ...
            cma(:, 1) >= RAl & cma(:, 1) <= RAu & cma(:, 2) >= decl & cma(:, 2) <= decu);

% For comparison with V filter in S&H
%shindex = find( snh_clr >= snh_cl & snh_clr <= snh_cu & m1 >= snh_ml & m1 <= snh_mu);          
            
% For comparison with I filter in Stetson
shindex = find( snh_clr >= snh_cl & snh_clr <= snh_cu & m2 >= snh_ml & m2 <= snh_mu);           
   
% Stellar type amd mag/color selections from Sandquist data
sandh = sandh(shindex, :);
type = type(shindex);
m1 = m1(shindex);
m2 = m2(shindex);
tindex = find(char(type) == stellar_type);
sandh = sandh(tindex, :);
type = type(tindex);
m1 = m1(tindex);
m2 = m2(tindex);
nsandh = length(tindex);

clc
disp(cstrcat("Number of observations from ", positions_list, " displayed = ",  int2str(nsandh)));

% Plot chosen section of S&H CMD
%h = figure(37);
% Set toolkit to use for this figure (= handle)
%graphics_toolkit(h, 'fltk')
% Clear this plot
%clf(h)
%fig = gcf ();
%ax = gca();
%set (fig, 'visible', 'on')
%set(ax, 'position', [0.12,0.15,0.75,0.78]);

% For Screen viewing
%set(gcf,'position', [200,400,800,800]);
%set(ax, 'xminortick', 'on');
%set(ax, 'yminortick', 'on');

%m1 = m1(tindex);
%m2 = m2(tindex);
%plot(m1-m2, m1, "ko", "markersize", 5.0)

%limits= [min(m1-m2), max(m1-m2), min(m1), max(m1)];
%axis (limits);
%xlabel("(m1-m2)", 'fontsize', 12);
%ylabel("m1", "fontsize",  12);
%legend(cstrcat("N = ", int2str(nsandh)))
%set(gca,'fontsize',12); % sets font of numbers on axes
%set (gca, 'ydir', 'reverse')

% Adding a user specified offset to RA and Dec
%[sr, sc] = size(sandh);
%sandh = sandh .+ repmat(manual_shift, sr, 1);

% Plotting map         
h = figure(30);
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

if resrows > 0
  plot(cma(:, 1)(index), cma(:, 2)(index), "b+", "markersize", 10.0, sandh(:, 1), sandh(:, 2), "ro", "markersize", 10.0,... 
       residuals(:,1), residuals(:,2), "r*", "markersize", 10.0, residuals(:,3), residuals(:,4), "b*", "markersize", 20.0);
       text(sandh(:,1),sandh(:,2), type, 'horizontalalignment', 'left', 'verticalalignment', 'bottom')
       text(cma(:,1)(index), cma(:,2)(index), int2str(index), 'horizontalalignment', 'left', 'verticalalignment', 'bottom')
       text(cma(:,1)(index), cma(:,2)(index), int2str(cma(:, 5)(index)), ...
                            'horizontalalignment', 'right', 'verticalalignment', 'bottom')
       text(cma(:,1)(index), cma(:,2)(index), int2str(cma(:, 6)(index)), ...
                            'horizontalalignment', 'right', 'verticalalignment', 'top')
      hold on
else 
  plot(cma(:, 1)(index), cma(:, 2)(index), "b+", "markersize", 10.0, sandh(:, 1), sandh(:, 2), "ro", "markersize", 10.0);
  text(sandh(:,1),sandh(:,2), type, 'horizontalalignment', 'left', 'verticalalignment', 'bottom')
  text(cma(:,1)(index), cma(:,2)(index), int2str(index), 'horizontalalignment', 'left', 'verticalalignment', 'bottom')
  text(cma(:,5)(index), cma(:,2)(index), int2str(index), 'horizontalalignment', 'right', 'verticalalignment', 'top')
 % text(cma(:,6)(index), cma(:,2)(index), int2str(index), 'horizontalalignment', 'right', 'verticalalignment', 'bottom')
  hold on 
 end
 
limits= [RAl, RAu, decl, decu];
axis (limits);
xlabel("Right Ascension (Degrees)", 'fontsize', 12);
ylabel("Declination (Degrees)", 'fontsize', 12);
legend(cstrcat("N = ", int2str(nsandh)))
set(gca,'fontsize',12); % sets font of numbers on axes

fsize = 14;

% Draw The WFC3 Footprints (black) this work 
% and the WFPC2 footprints used by S&H (red) on the map
draw_primitive_p (footprint1, 'k', 1, fsize);
hold on
draw_primitive_p (footprint2, 'k', 4, fsize);
hold on
draw_primitive_p (footprint3, 'g', 8, fsize);
hold on
draw_primitive_p (footprint4, 'm', 7, fsize);
hold on
draw_primitive_p (footprint5, 'c', 2, fsize);
hold on
draw_primitive_p (footprint6, 'c', 2, fsize);
hold on

cont = 1;
while (cont == 1)

  % Collect position pairs for coordinate transformation
  [x, y, button] = ginput();
  bindex = find (button==1);
  coords = [x(bindex), y(bindex)];
  [nr, nc] = size(coords);
  oindex = [1:2:nr];
  refindex = [2:2:nr];
  obs_coords = coords(oindex, :);
  ref_coords = coords(refindex, :);
  
  % Identify reference objects selected on the map in the CMA file.
  % Identify observations selected on the map in the Sandquist and Hess data file.
  
 
  disp("               ID        WFC3 Mag    Color")
 slindices = [];
 clindices = [];
 for i = 1:nr/2 
    deltas = sandh(:, 1:2).- repmat(obs_coords(i, :), nsandh, 1);
    rads = (deltas(:, 1).^2 .+ deltas(:, 2).^2).^5;
    [sl, slindex] = min(rads);
  
    delta = cma(:, 1:2).- repmat(ref_coords(i, :), cmar, 1);
    rad = (delta(:, 1).^2 .+ delta(:, 2).^2).^5;
    [cl, clindex] = min(rad);
    
    % Transform S&H magnitudes (WFPC2) to Gillam magnitudes (WFPC3)
    %tcol = cma(clindex,3);
    tcol = m1(slindex) - m2(slindex);
    refmag = cma(clindex,4);
    
    % For S&H V  
    %tmag = m1(slindex) .+ c0 .+ c1 .* tcol + c2 .* tcol.^2;
    
    % For Stetson I
    tmag = m2(slindex) .+ c0 .+ c1 .* tcol + c2 .* tcol.^2;
    
    disp(cstrcat("Observed       ", strjust(char(type(slindex))), "      (", num2str(tmag, "%5.3f"), ")      ", num2str(tcol, "%5.3f") ))                                               
    
    
    disp(cstrcat("Reference     ", num2str(clindex, "%6i"), "        ", num2str(refmag, "%5.3f"), "      ", num2str(cma(clindex,3), "%5.3f")))
                                                 
    disp("  ")

    slindices = [slindices, slindex];
    clindices = [clindices, clindex];
  end
  
  % Temporary copies for predicting locations
  tobs = sandh(slindices, 1:2);
  tref = cma(clindices, 1:2);
  
    
  % Use the catalog coordinates instead of the measure ones 
  savethese = yes_or_no("Add these stars to the data? ");
  
  if (savethese == 1) 
      obs_coords = tobs';
      ref_coords = tref';
  end
  
  if resrows > 0 
      if (savethese == 0) 
        ref_coords = [residuals(:, 1:2)'];
        obs_coords = [residuals(:, 3:4)'];
      else
  
      % Add to previous measurements
        ref_coords = [residuals(:, 1:2)', ref_coords];
        obs_coords = [residuals(:, 3:4)', obs_coords];
      end
  end
  
  if (rows(tobs)*savethese == 0) && (resrows == 0)
      disp ("Select at least five matches to start the residuals file")
      return;
  end
  
  
  % Fit Scaled rotation model to measured coordinates
  % Starting guess at parameters
  % RA_shift, dec_shift - the shift of the center of the object positions to the reference positions
  % rotation angle, theta, = 0;
  % RA and dec Scale factors both 1.
  %
  % Initial guesses
  RA_shift = mean(ref_coords(1,:) - obs_coords(1, :));
  dec_shift = mean(ref_coords(2,:) - obs_coords(2, :));
  model_params = [RA_shift, dec_shift, 0, 1, 1];

  if strcmp(use_external_model, "Yes")
      model_params = [shift_RA, shift_dec, theta0, RA_scale, dec_scale, non_par_factor];
  end
  delta_RA_init = model_params(1);
  delta_dec_init = model_params(2);
  Theta_init = model_params(3);
  S_RA_init = model_params(4);
  S_dec_init = model_params(5);
  np_init = model_params(6);

  data_wts = 1./obs_coords.^0.5;
  [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]= leasqr(obs_coords, ref_coords, model_params, 'Scaled_Rotation_Model_p', 1e-9, 10000, data_wts);

  f = reshape(f, 2, length(f)/2);

  chi_RA = (mean(mean((f(1,:) .- ref_coords(1,:)).^2))).^.5
  chi_Dec = (mean(mean((f(2,:) .- ref_coords(2,:)).^2))).^.5

  if exist('temp.resids', 'file')
    delete('temp.resids')
  end

  % Plot current solution
  figure(30)
  hold on
  predictions = Scaled_Rotation_Model_p(tobs', p);
  plot(predictions(1, :),predictions(2, :),"rs", "markersize", 10.0)
  plot(f(1, :), f(2, :), "bs", "markersize", 10.0)
  hold on
  [pr, pc] = size(predictions);
  
  % Print predictions to the screen
  blanks = repmat("    ", pc, 1);
  disp("Predicted reference positions for")
  disp("   RA  (J2000)   Dec         ID");
  disp(cstrcat(num2str(predictions', "%10.7f   "), blanks, strjust(char(type(slindices)))))

  % Save the residuals 

  residuals = [ref_coords', obs_coords', f', ref_coords'.-f'];
  [resrows, rescols] = size(residuals);
 
  dlmwrite ('temp.resids', residuals,  "precision", "%10.7f", "delimiter", " ");
  [status, msg, msgid] = copyfile ('temp.resids', residuals_file);

  delete('temp.resids');
  
  cont = yes_or_no("Continue collecting data? ");
  figure(30)
                
end
hold off

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
set(gca,'fontsize',14); % sets font of numbers on axe s
legend(cstrcat("N = ", int2str(resrows)))

figure(32);
fig = gcf ();
ax = gca();
hist(residuals(:,8),nbins_dec)
title (cstrcat('DeclinationResiduals'), 'fontsize', 12)
xlabel("Observed Dec - Estimated Dec (Degrees)", 'fontsize', 12);
ylabel("Frequency", 'fontsize', 12);
set(gca,'fontsize',14); % sets font of numbers on axes
legend(cstrcat("N. = ", int2str(resrows)))

% Display Results

delta_RA = p(1);
delta_dec = p(2);
Theta = p(3);
S_RA = p(4);
S_dec = p(5);
np = p(6);

sig_delta_RA = covp(1,1)^.5;
sig_delta_dec = covp(2,2)^.5;
sig_Theta = covp(3,3)^.5;
sig_S_RA = covp(4,4)^.5;
sig_S_dec = covp(5,5)^.5;
sig_np = covp(6,6)^.5;

disp ("Initial Guess")
disp(cstrcat("Center shift: RA (J2000), Dec (J2000) = ", num2str(delta_RA_init, "%10.7f"), " deg ", num2str(delta_dec_init, "%10.7f"), " deg"))
disp(cstrcat("Rotation: Theta= ", num2str(Theta_init, "%10.7f"), " radians"))
disp(cstrcat("Scale factors: RA= ", num2str(S_RA_init, "%9.4f"), " Dec = ", num2str(S_dec_init, "%9.4f"), ...
                                 " Non-Orthogonality = ", num2str(np_init, "%9.4f")))

disp ("Solution")
disp(cstrcat("Center shift: RA (J2000), Dec (J2000) = ", num2str(delta_RA, "%10.7f"), " +/- ", num2str(sig_delta_RA, "%10.7f"), " deg ", ...
                                                   num2str(delta_dec,"%10.7f"), " +/- ", num2str(sig_delta_dec, "%10.7f"), " deg"))

disp(cstrcat("Rotation: Theta= ", num2str(Theta, "%10.7f"), " +/- ", num2str(sig_Theta, "%10.7f"), " radians" ))
disp(cstrcat("Scale factors: RA= ", num2str(S_RA, "%9.4f"), " +/- ", num2str(sig_S_RA, "%9.4f"), " Dec = ", ...
                                    num2str(S_RA, "%9.4f"), " +/- ", num2str(sig_S_RA, "%9.4f"), " Non-Orthogonality = ", ...
                                    num2str(np, "%9.4f"), " +/- ", num2str(sig_np, "%9.4f")))
  
  
% Return to starting directory
% Unload Optimization package
pkg unload optim
cd(currentdir);


