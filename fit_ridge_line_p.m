function mag2newcolor = fit_ridge_line_p(infile, data, use, dcolor, dmag, cl, cu, ml, mu)
%
% Fit a ridge line to a CMD 
%
% infile = the input file name
% data = Color-magnitude array
% dmag = size of magnitude side of statistics box
% dcolor = size of color side of statistics box
% cl and cu are the lower and upper limits on colors that will be plotted
% ml and mu are the lower and upper limits on mags. that will be plotted
% use = flag for whether to use ridge line from previous run.
%
% CMA entries will be summed over non-overlapping boxes dcolor x dmag in area.
%
% Make up the output file names based on the input file name
%
if nargin == 0
   clc
   disp ("This script creates a color-magnitude density plot (Hess diagram).");
   disp ("It imposes a grid on the color-magnitude diagram and counts the number of");
   disp ("points inside each grid rectangle.");
   disp ("  ");
   disp ("You must provide at least the name of the .cma or .cm3 file");
   disp ("  ");
   disp ("Syntax: hess_plot_p(arg1, arg2, arg3, arg4, arg5)");
   disp ("Where,");
   disp ("     arg1 = the .cma or .cm3 file name");
   disp ("     arg2 = size of color dimension of summing box");
   disp ("     arg3 = size of magnitude dimension of summing box");
   disp ("     arg4 and arg5 are the lower and upper limits on mags. that will be plotted");
   disp ("  ");
   disp ("If these restrictions are violated, the relevant axis will display bin numbers instead of magnitudes.");
   disp ("  ");
   disp ("If less than four arguments are given, default magnitude limits wil be used.");
   disp ("If less than two arguments are given, default mag. bines and magnitude limits wil be used.");
   disp ("  ");
   return
end

% Number of observations
 nrows = rows(data);

 % Make up input ridge line file name.
out_file_base = strtok (infile, ".");
ofb = strsplit (out_file_base, "\\");
[rofb, cofb] = size(ofb);
out_file_base = cstrcat("HessPly\\", char(ofb(cofb)));
outfile = strcat(char(out_file_base(1, :)), '.hess');
polyfile = strcat(outfile, '.ply');
 
% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

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
   clabel = cstrcat(filter1, " - ", filter2);
end
if r == 2
   clabel = cstrcat(filter2);
end

mlabel = cstrcat(filter1);


% Choose first or second filter depending on if the
% filename contains the string "ALT" or not.

mlabel = cstrcat(filter1);
for i=3:cf
   altind = char(front(i));
   if strcmp(altind, "ALT")
      mlabel = cstrcat(filter2);
   end
end

h = figure(500);
% Set toolkit to use for this figure (= handle)
graphics_toolkit(h, 'fltk')
% Clear this plot
clf(h)
fig = gcf ();
ax = gca();
set (fig, 'visible', 'on')
set(ax, 'position', [0.15,0.08,0.81,0.90]);
%set(gcf,'position', [200,75,375,715]);
   
set(ax, 'xminortick', 'on');
set(ax, 'yminortick', 'on');
plot(data(:,3), data(:,4), 'bs', 'markersize', 1.5, "markerfacecolor", "auto")
hold on
axis ([cl, cu, ml, mu]);
set (gca, 'ydir', 'reverse')
xlabel(clabel, 'fontsize', 18);
ylabel(mlabel, 'fontsize', 18);

set(gca,'fontsize',18); % sets font of numbers on axes
set(ax, 'xminortick', 'on');
set(ax, 'yminortick', 'on');

% Draw a curve on the CMD
% Save ridge points

% Load existing ridge-line
if exist(polyfile, 'file')
   saved_ridge = load(polyfile);
else
   use = 0;
end

if use == 1
  mhats = saved_ridge(:,1);
  nc_smooth = saved_ridge(:,2);
  
  % Fit a spline to the saved ridge points
  color2mag = spline_fit_ISOs_p(3, 3, nc_smooth, mhats);
  mag2color = spline_fit_ISOs_p(3, 3, mhats, nc_smooth);
else 
  more off
  disp("Select ridge points on the plot by left-mouse-clcks")
  more on
  [xr, yr, button] = ginput();
  binds = find(button == 1);
  xr = xr(binds);
  yr = yr(binds);
  
  

  % Fit a spline to the ridge points
  color2mag = spline_fit_ISOs_p(2, 3, xr, yr);
  mag2color = spline_fit_ISOs_p(2, 3, yr, xr);

  % Plot the ridge-line
  % Set toolkit to use for this figure (= handle)
  h=figure(550);
  graphics_toolkit(h, 'fltk')
  % Clear this plot
  clf(h)
  fig = gcf ();
  ax = gca();
  set (fig, 'visible', 'on')
  set(ax, 'position', [0.15,0.08,0.81,0.90]);
  %set(gcf,'position', [200,75,375,715]);
   
  set(ax, 'xminortick', 'on');
  set(ax, 'yminortick', 'on');
  plot(data(:,3), data(:,4), 'bs', 'markersize', 1.5, "markerfacecolor", "auto")
  hold on
  axis ([cl, cu, ml, mu]);
  set (gca, 'ydir', 'reverse')
  xlabel(clabel, 'fontsize', 18);
  ylabel(mlabel, 'fontsize', 18);

  set(gca,'fontsize',18); % sets font of numbers on axes
  set(ax, 'xminortick', 'on');
  set(ax, 'yminortick', 'on');
  plot(xr, yr, "b+", "markersize", 3.5)
  plot(xr, ppval(color2mag, xr), "ro", "markersize", 3.5, "markerfacecolor", "auto")
end

if nargin == 9
   minmag = ml;
   maxmag = mu;
end

% Auto-select the number CMD boxes
if nargin >= 7 && nargin < 9
   minmag = min(data(:,4));
   maxmag = max(data(:,4));
end

if nargin == 3
  minmag = min(data(:,4));
  maxmag = max(data(:,4));
  mincol = min(data(:,3));
  maxcol = max(data(:,3));  
  dcolor = 0.02;
  dmag = 0.05;
end

nmags = int32(ceil((maxmag - minmag) / dmag) + 1);

mhats = linspace(minmag, maxmag, nmags);
chats =  ppval(mag2color, mhats);
clower = chats - dcolor/2;
cupper = chats + dcolor/2;

for j = 1:nmags
  maglower = minmag + (j-1)*dmag;
  magupper = maglower + dmag;

  ind = find( data(:,3) > clower(j) & data(:,3) <= cupper(j));
  ninds = length(ind);
  if ninds == 0
     median_color = chats(j);
   else
     median_color = median(data(ind, 3));
   end
   new_colors(j) = median_color;
end
   
% Fit a spline to the ridge points
mag2newcolor = spline_fit_ISOs_p(3, 3, mhats, new_colors);

dnc_smooth =  ppval(mag2newcolor, mhats);  
plot(dnc_smooth, mhats, "g-", "linewidth", 3.5)
hold off
% Axes label control
set (gca, 'ydir', 'reverse')
%title (title_str, 'fontsize', 18)
xlabel(clabel, 'fontsize', 16);
ylabel(mlabel, 'fontsize', 16);
set(gca,'fontsize',15); % sets font of numbers on axes


% Statistics
% Save processing history
log = cstrcat("hess_plot_p_", int2str(int64(time)),".log");
diary (log);
processing_time = strftime("%Y-%d-%m %T %z\n", localtime(time()));
disp ("          ");
disp ('------------------------------------------------------');
disp ("           fit_ridge_line_p.m Statistics       ");
disp ('------------------------------------------------------');
disp (cstrcat(infile, ' ------> ', outfile));
disp (cstrcat(infile, ' ------> ', polyfile));
disp (cstrcat("Processed: ", processing_time));
disp ("  ");
disp (cstrcat("Number of observations  = ",int2str(nrows)));
disp ("  ");
disp ("Statistics box spec.");
disp (cstrcat("Color interval          = ",num2str(dcolor), " mag"));
disp (cstrcat("Magnitude interval      = ",num2str(dmag), " mag"));
disp ('------------------------------------------------------');
disp ("          ");

disp(cstrcat("Screen statistics has been saved to ", log))
diary off;

cd(currentdir);

return;
