function color_dist_p(infile, nbmag, nbcolor, ng, maglower, magupper, clower, cupper, rcl, rcu, sfin)
%
% Divide a CMD sequence into a number of magnitude bins and plot the distributions of color within them
%
% Inputs:
%     infile = Color-magnitude array
%     nbmag = number of magnitude bins
%     nbcolor = number of color bins
%     ng = number of Gaassians in the mixture model
%     maglower = lower magnitude limit of CMA sample
%     magupper = upper magnitude limit of CMA sample
%     clower = lower mag. limit of CMD sample
%     cupper = upper mag. limit of CMD sample
%     rcl = lower rectified color limit of sample used to model color distribution
%     rcu = upper rectified color limit of sample used to model color distribution
%     sf = Scale factor on Gaussians to ensure that variance is greater than unity.
%     (Note: if the computed variance is less than unity gmmem.m fails.)
%
%
% Outputs
%     Output file containing an array of histogram points at various magnitudes.
%     nbmag histograms to the screen
%     color_dist.JPEG - A JPEG file containing the histograms.
%
% Scipts called
%     mag2newcolor = fit_ridge_line_p(infile, data, use, dcl, dma,  clower, cupper, maglower, magupper)
%     [mixtures, p_i, o, chisqr] = gmmem(X, mixtures, f, cbin, .1, sf)
%     Z = plotGMM_1D(x, mixtures, p_i)
%
% Use gnuplot because the TeX interpreter does not work with fltk
graphics_toolkit("gnuplot")

% Output file base name
out_file_base = strtok (infile, ".");
ofb1 = strsplit (out_file_base, "\\");
[rofb1, cofb1] = size(ofb1);
out_file_base = cstrcat(char(ofb1(cofb1)));
outfile = cstrcat("vrtcma\\",out_file_base,".vcma");

% Extract filters from input file name
%front = strtok(infile, ".");
%front = strsplit(front, "_");
%[rf, cf] = size(front);
%filter1 = char(front(2));
%filter2 = char(front(3));

%fcolor = cstrcat("[",filter1,"-",filter2,"}");

% Choose first or second filter depending on if the
% filename contains the string "ALT" or not.

%mlabel = cstrcat("m_{",filter1'"}");

% Make up CMD axis lables from file name
parts = strsplit(infile, ".");
front = char(parts(1));
ext = char(parts(2));
front = strsplit(front, "_");
[rf, cf] = size(front);
filter1 = char(front(2));
filter2 = char(front(3));

[c, r] = size(strsplit(char(front(3)), "-"));

if r == 1
   clabel = cstrcat("m_{",filter1, "}", " - ", "m_{", filter2,"}")
end
if r == 2
   fil2 = strsplit(filter2, "-");
   f1 = char(fil2(1));
   f2 = char(fil2(2));
   clabel = cstrcat("m_{",f1, "}"," - ", "m_{",f2, "}");
end

mlabel = cstrcat("m_{", filter1, "}");

for i=3:cf
   altind = char(front(i));
   if strcmp(altind, "ALT")
      mlabel = cstrcat("m_{",filter2,"}");
   end
end

cdlabel = "Color Distribution";
if (nbmag > 1)
   cdlabel = cstrcat(cdlabel,"s");
end

dclabel = cstrcat('{\Delta}', "(", clabel, ")");

% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

% Load the CMA file
if exist(infile, 'file')
   data = load(infile);
else
   disp(cstrcat("File ", infile, " Dose not exist"));
   return
end

% Total number of observations in CMA
total = rows(data);
ncols = columns(data);

% Get max. and min. magnitudes anc colors
maxmag = max(data(:, 4));
minmag = min(data(:, 4));

% User specified magnitude limits
if (nargin > 5 )
   maxmag = min(magupper, maxmag);
   minmag = max(maglower, minmag);
end

index = find(data(:, 4) >= minmag & data(:, 4) < maxmag);

dmag = (maxmag - minmag) / nbmag;

mincolor = min(data(:, 3)(index));
maxcolor = max(data(:, 3)(index));

% User specified color limits
if (nargin >= 8)
   maxcolor = cupper;
   mincolor = clower;
end

% User specified rectified color limits
if (nargin >= 10)
   maxrcolor = rcu;
   minrcolor = rcl;
else
   % Defaults that do nothing.
   maxrcolor = 10000;
   minrcolor = -10000;
end

% Apply the color-magnitude selection criteria
index = find(data(:, 4) >= minmag & data(:, 4) < maxmag & data(:, 3) >= mincolor & data(:, 3) < maxcolor);

dcolor =  (maxcolor - mincolor) / nbcolor;
pcbin = linspace(mincolor+dcolor/2, maxcolor-dcolor/2, nbcolor);

colrs = data(:, 3)(index);
mags = data(:, 4)(index);

% Number of observations in sample
nrows = rows(data(:, 3)(index));

drcolor = (maxrcolor - minrcolor)/nbcolor;
rcbin = linspace(minrcolor+drcolor/2, maxrcolor-drcolor/2, nbcolor);

%clf

% Putting the label in the top-right corner of the plot(s), appriately scaled.
tfontsize = [13,   13,   12,   12,   11,   11,   10,   10,   9,    9   ];
xploc =     [0.67, 0.65, 0.67, 0.67, 0.71, 0.71, 0.71, 0.71, 0.74, 0.75] -.32;
yploc =     [0.98, 0.96, 0.95, 0.93, 0.91, 0.90, 0.89, 0.87, 0.83, 0.81];

% Load ridgeline
polyfile= cstrcat("HessPly\\",out_file_base, ".hess.ply");
if exist(polyfile, 'file')
   c = load(polyfile); % the ridge line, upper and lower bounds of the ROI
%else
%   clc
%   disp(cstrcat("To continue: Make a ridgeline file from ", infile));
%   return
end

% Find the ridge-line of the most populous group and make a spline fit to it.
dcl = 0.025;
dma = 0.025;
use = 0;
nc = 300;
mag2newcolor = fit_ridge_line_p(infile, data, use, dcl, dma,  clower, cupper, maglower, magupper);
fidmags = linspace(maglower, magupper, nc);
new_colors =  ppval(mag2newcolor, fidmags);
magl = max(maglower, min(data(:, 4)));
magu = min(magupper, max(data(:, 4)));
mag_upper = max(magu, max(fidmags));
mag_lower = min(magl, min(fidmags));
dm = (mag_upper-mag_lower)/nc;

% Set up figure
histlist = [];
h = figure(16);
clf(h);
fig = gcf ();
ax = gca();
set(ax, 'xminortick', 'on');
set(ax, 'yminortick', 'on');
set (fig, 'visible', 'on')
set(gcf,'position', [200,25,1600,900]);
set(gca,'interpreter', "tex");

% Display Raw CMA
left_edge = 0.07;
col_spacing1 = 0.30;
col_spacing2 = 0.30;
ledge2 = left_edge + col_spacing1;
ledge3 = ledge2 + col_spacing2;

bottom = 0.07;
height = 0.88;
top = bottom + height;
width = 0.25;

subplot("position", [left_edge, bottom, width, height]);
hold off
plot(colrs,mags,"s","markersize", 1.0)
hold 

% PLot ridge-line of most populous group
plot(new_colors, fidmags, "r", "linewidth", 2.5);

set (gca, 'ydir', 'reverse')
%title (cstrcat(infile), 'fontsize', 18)
axis([mincolor, maxcolor, magl, magu]);
limits = axis;
xlabel(clabel, 'fontsize', 14);
ylabel(mlabel, 'fontsize', 14, 'rotation', 90);
%legend(cstrcat("# Obs. = ", int2str(nrows)));
hold off

% Subtract rige-line colors from raw data
colrs = colrs - ppval(mag2newcolor, mags);  % Rectified colors for display
colvrt = data(:,3) - ppval(mag2newcolor, data(:, 4)); % For the output verticalized CMA file.
vrtcma = [data(:, 1:2), colvrt, data(:, 4)];
if ncols > 4
   vrtcam = [vrtcma, data(:, 5:ncols)];
end
% Save list of selected objects
dlmwrite ('temp.vrt',  vrtcma, "precision", "%3.8f", "delimiter", " ");
[status, msg, msgid] = copyfile ('temp.vrt', outfile);
delete('temp.vrt');

% Display Color distibutions
% Find nax. and min. colors of entire sample.
maxcol = max(colrs);
mincol = min(colrs);
m = 1; % Dimension of dataset. i.e 1,2,3....
sf = sfin;
p_is = [];
mus = [];
sigs = [];
ndata = [];
dclabels = [];
mups = [];
mls = [];

for i = 1:nbmag
   % Select obs. in the magnitude bin
   ml = minmag +(i-1)*dmag;
   mup = minmag + i*dmag;

   % For the statistics output
   mups = [mups, mup];
   mls = [mls, ml];

   % Select a magnitude slice of the data
   cindex = find(mags >= ml & mags < mup);
   colors = colrs(cindex);
   
   % Select observations based on their rectified colors.
   rindex = find(colors > minrcolor & colors <= maxrcolor);
   colors = colors(rindex);

   %subplot(nbmag, 3, 3*i, 'align');
   %subplot(nbmag, 3, 3*i, "position", [0.7, 1-i/nbmag, 0.4, 1/nbmag]);
   subplot("position", [ledge3, top-height*i/nbmag, width, 0.999*height/nbmag]);
      
   % Scale data so as to make varaiance greater than one.
   % EM routine fails if Var(X) < 1.
   X = sf*colors';

   %xmax = max(max(X)); % This is where the problem is!
   %xmin = min(min(X));
   xmin = rcl*sf;
   xmax = rcu*sf;
   %dm = (xmax-xmin)/nbcolor;
   dm = (rcu-rcl)*sf/nbcolor;
   
   size(X);
   
   % Data histogram
   %[f, cbin] = hist (X/sf, nbcolor);
   [f, cbin] = hist (X/sf, nbcolor);
      
   normZ = sum(f); % Total number of stars in the sample.
   
   % Initial mixture model means and variances

   for k = 1:ng
      % Contrived guestimate for mean1 and sigma1
      meanX = mean(X, 2)+randn(m,1)*2;
      mixtures(k).mu = meanX;
      mixtures(k).sigma = eye(m,m)*std(X)^2;
   end
  mixtures(1).mu = 0.0;
  mixtures(2).mu = -0.04;
  %mixtures(2).mu = -0.07;
  %mixtures(1).sigma = 0.02;
  %mixtures(2).sigma = 0.02;
  %mixtures(2).sigma = 0.02;
  
   % EM Stopping tolerance (third argument) set to delta-log-likelyhood= 1e-5.
   [mixtures, p_i, o, chisqr] = gmmem(X, mixtures, f, cbin, .1, sf);
   
   p_is = [p_is; p_i];
   iters(i) = o;
   chis(i) = chisqr;
   ndata = [ndata, sum(f)];

   % Means and variances
   means = [];
   sigmas = [];
   for k = 1:ng
      mu = mixtures(k).mu/sf;
      sig = mixtures(k).sigma/sf^2;

      % For output stats.
      means = [means, mu];
      sigmas = [sigmas, sig];
   end
   mus = [mus; means];
   sigs = [sigs; sigmas];

   % For plotting
   x = [xmin:dm/50:xmax]; % Points on the total prob. distribution
   Z = plotGMM_1D(x, mixtures, p_i);
    
   %hist(X/sf, nbcolor, "facecolor", "g")
   hist(X/sf, rcbin, "facecolor", "g")

   % Plot the estimated population PDF
   hold on
   
   
   plot(x/sf, Z*dm*normZ, "r", "linewidth", 2.5);
   for ii = 1:ng
      z = plotGMM_1D(x, mixtures(ii), 1);
      normz(ii) = p_is(ii)*sum(f)*dm;
      plot(x/sf, z*normz(ii), "b", "linewidth", 1.75);
   end
   
   axis([minrcolor, maxrcolor])
   set(gca,'interpreter', "tex");
   limits = axis;
   Xp = limits(1) + 0.97 * (limits(2) - limits(1));
   Yp = limits(3) + yploc(nbmag) * (limits(4) - limits(3));

   axis("tic", "labely");
   
   text(Xp, 0.92*Yp, sprintf('%3.2f - %3.2f', ml, mup), ...
                        "horizontalalignment", "right", ....
                        "verticalalignment", "middle","fontsize", 10);
   text(Xp, 0.8*Yp, sprintf("{\\chi}^{2}=%3.2f", chis(i)), ...
                        "horizontalalignment", "right", ....
                        "verticalalignment", "middle", "fontsize", 10)
  
   if (i == nbmag)
      xlabel (dclabel, 'fontsize', 14);
      axis("tic", "label[xy]");
   end
   if (i == int32(nbmag/2))
      ylabel ("N", 'fontsize', 14, 'rotation', 0);
   end
   if (i == 1)
 %     title(cdlabel, 'fontsize', 16);
   end

   % Save for statistics.
   dclabels = [dclabels, dclabel];

   hold off
end

% Display rectified CMA
subplot("position", [ledge2, bottom, width, height]);

hold on
plot(colrs, mags,"s","markersize", 0.25)
plot([0,0], [magl, magu], "r", "linewidth", 2.5);

set(gca, 'ydir', 'reverse')
set(gca,'interpreter', "tex")
axis("tic", "labelx");
axis([minrcolor, maxrcolor, magl, magu]);
limits = axis;
%title (cstrcat(infile), 'fontsize', 18)
xlabel(dclabel, 'fontsize', 14);
%ylabel(mlabel, 'fontsize', 12);
%legend(cstrcat("# Obs. = ", int2str(nrows)));
hold off

% Save the plot to a JPEG file.
%print -deps color_dist


% Statistics

% Save processing history
log = cstrcat("logs\\color_dist_p_", int2str(int64(time)),".log");
diary (log)
processing_time = strftime("%Y-%d-%m %T %z\n", localtime(time()));
disp ("          ");
disp ('-------------------------------------------------------------------------');
disp ("                           color_dist_p.m Statistics       ");
disp ('-------------------------------------------------------------------------');
disp (cstrcat("Input CMA: ", infile));
disp (cstrcat("Input ridge line: ", polyfile));
disp (cstrcat("Processed: ", processing_time));
disp ("          ");
disp (cstrcat("Color: ", dclabel));
disp (cstrcat("Color range: ",num2str(rcl), " to ",num2str(rcu)));
disp (cstrcat("Number of color bins           = ",int2str(nbcolor)));
disp (cstrcat("Number of magnitude bins       = ",int2str(nbmag)));
disp ("  ");
for i = 1:nbmag
   disp (cstrcat("Magnitude bin                  = ", num2str(mls(i)),"<",  mlabel, "<=", num2str(mups(i))));
   disp (cstrcat("Number of Gaussians in model   = ",int2str(ng)));
   disp (cstrcat("Proportion of each Gaussian    = ",num2str(p_is(i,:), " %3.3e")))
   disp ( cstrcat("Means                          = ", num2str(mus(i, :), " %3.3e")) );
   disp ( cstrcat("Variances                      = ", num2str(sigs(i, :), " %3.3e")) );
   disp ( cstrcat("N                              = ", num2str(ndata(i), " %7i")) );
   disp ( cstrcat("Chi-square                     = ", num2str(chis(i),"% 3.3f") )  );
   disp (cstrcat("Number of iterations           = ",int2str(iters(i))));
   if i < nbmag
      disp ("  ");
   end
end
disp ('-------------------------------------------------------------------------');
disp ("   ");

disp(cstrcat("Screen statistics has been saved to ", log))
diary off

% Return to original directory
cd(currentdir);

return;

