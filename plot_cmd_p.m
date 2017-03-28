function plot_cmd_p(infile, plot_type, rcl, rcu, rml, rmu, rRAl, rRAu, rDecl, rDecu)
%
% Plot color-magnitude and color_color diagrams from two HST Daophot Photometry Lists
% comma-separated text file
% This primitive is intended for use in other functions and
% programs
%
% Usage: plot_cmd_p(stars_XXXXX_YYYYY.cma, type)
%   stars_XXXXX_YYYYY.cma = input color-magnitude array
% where XXXXX is something like F555W
%       YYYYY is something like F814W
%
% infile = input color-magnitude array. (e.g. 'w3w3_F555W_F814W.cm3')
%           or
%          input color-color array (e.g. 'w3w3w3w3_F555W-F775W_F555W-F814W.cc3')
%
% plot_type = "cmd" - make a color magnitude diagram (or a color-color diagram) only.
%             "mono" - make a CMD with monochrome (blue) symbols
%             "montage" - arrays several mono plots
%             "overplot" - like montage, except it does not clear the figure first. 
%                          You must set the marker color by setting the PLOT_COLOR environment variable
%             "magtran" - apply a third-order polynomial transformation before plotting
%                         Reads the coefficeints from the geotran parameter file
%             "shape" - Plot a shape over the CMD. (Currently one shape is hard-coded).
%             "map" - make a map only
%             "both" - make a cmd or color-color diagarm, and a map
%             "none" - no plots
%
% RAl, RAu = Loweer and upper Right Ascensions to plot.
% Decl, Decu = Loweer and upper Declinationss to plot.
% rcl, rcu =  lower and upper color limits of a CMD box
% rml, rmu =  lower and upper magnitude limits of a CMD box
clc

% Set graphics tookit to fltk
% This is used because there are number of bugs in GNUplot that prevent scatter plots with
% marker colors that vary with data point from working.

graphics_toolkit fltk

% Output file base names
out_file_base = strtok (infile, ".");
ofb1 = strsplit (out_file_base, "\\");
[rofb1, cofb1] = size(ofb1);
out_file_base = cstrcat(char(ofb1(cofb1)));
db_out = strcat(out_file_base, '.db');

% Save current working directory
thisdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

if  nargin < 2
   clc
   disp("The color magnitude (or color-color) file name and plot type must be supplied:");
   disp("   Argument 1 - Color-magnitude array. (e.g. 'w3w3_F555W_F814W.cm3')");
   disp('   Argument 2 - Plot type - "map", "cmd", "both", or "none"');
   disp(' The "cmd" option: If you supply a color-color diagram (.cc3) file instead of a CMA file,');
   disp(" the script will plot a color-color diagram instead.");
   disp("  ");
   disp("Optional arguments are:");
   disp("   Argument 3 - Lower color limit");
   disp("   Argument 4 - Upper color limit");
   disp("   Argument 5 - Lower magnitude limit");
   disp("   Argument 6 - Upper magnitude limit");
   disp("   Argument 7 - Lower Right Ascension limit");
   disp("   Argument 8 - Upper Right Ascension limit");
   disp("   Argument 9 - Lower Declination limit");
   disp("   Argument 10 - Upper Declinatoinlimit");
   disp(" If you supply an odd number of args, the program will use all the complete pairs and use the");
   disp(" defualts on remaining parameters.");
   disp("  ");
   disp("Examples:");
   disp("   plot_cmd_p('w3w3_F555W_F814W.cm3', 'cmd', 1, 2, 23, 26) - will plot a CMD restricting the");
   disp("   color range to 1 - 2 amgs and the magnitude range to 23 - 26 mags");
   disp(" ");
   return
end

% Load four-column CMD array.
if exist(infile, 'file')
   data = load (infile);
else
   disp(cstrcat(infile, " does not exist."));
   return
end

% Plot window and selection defaults
cl = -20;
cu = 20;
ml = -20;
mu = 50;
RAl = 0;
RAu = 360;
Decl = -90;
Decu = 90;

% User color limits
if nargin >=4
   cl = rcl;
   cu = rcu;
end

% User magnitude limits
if nargin >=6
   ml = rml;
   mu = rmu;
end

% User RA limits
if nargin >=8
   RAl = rRAl;
   RAu = rRAu;
end

% User Dec limits
if nargin == 10
   Decl = rDecl;
   Decu = rDecu;
end

% Apply selection criteria
index = find(data(:, 3) >= cl & data(:, 3) <= cu & data(:, 4) >= ml & data(:, 4) <= mu & ...
                data(:, 1) >= RAl & data(:, 1) <= RAu & data(:, 2) >= Decl & data(:, 2) <= Decu );

t1 = time;

RAs = data(:, 1)(index);
Decs = data(:, 2)(index);
colors = data(:, 3)(index);
magnitudes = data(:, 4)(index);
[rdata, cdata] = size(data);
if cdata >= 13
   CI1s = data(:, 12)(index);
   CI2s = data(:, 13)(index);
end
[nrows, ncols] = size(RAs);

% Sort by magnitudes
[sorted_mags, index] = sort(magnitudes);
sorted_colors = colors(index);
sorted_RAs = RAs(index);
sorted_Decs = Decs(index);
if cdata >= 13
   sorted_CI1s = CI1s(index);
   sorted_CI2s = CI2s(index);
end

% Output files go to thw work directory
%Index the output
nsorted = length(sorted_Decs);
sindex = [1:1:nsorted];

% Sort output by CI radius
%CIs = ((sorted_CI1s - 2.3).^2 + (sorted_CI2s - 2.5).^2).^.5;
%[sorted_CIs, index] = sort(CIs, "ascend");

% Save the Coordinates of the CMD selection
%dlmwrite ('temp.db',  [sorted_RAs(index), sorted_Decs(index), sorted_CI1s(index), sorted_CI2s(index), sindex'], %"precision", "%3.7f", "delimiter", " ");

dlmwrite ('temp.db',  [sorted_RAs(index), sorted_Decs(index)], "precision", "%3.7f", "delimiter", " ");
[status, msg, msgid] = copyfile ('temp.db', db_out);

delete('temp.db');

% Return to current working directory
cd(thisdir)


% Add labels to color-magnitude diagram
% Marker sizes represent brightness
mag_max = max(sorted_mags);
mag_min = min(sorted_mags);


% Make up CMD axis lables from file name

parts = strsplit(infile, ".");
front = char(parts(1));
ext = char(parts(2));
front = strsplit(front, "_");
[rf, cf] = size(front);
filter1 = char(front(2));
filter2 = char(front(3));
mlabel = cstrcat(filter1);

% Get size of second filter string in filnename
% If it holds one filter name the plot will be a color-magnitude diagram.
% If it holds two, it will be color-color diagram.
collab = strsplit(char(front(3)), "-");
[c, r] = size(collab);

% Plot CMD or color-color diagram as indicated by the input file name.
if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
   if r == 1
      col_label = cstrcat(filter1, " - ", filter2);
      if strcmp(plot_type, "ABSMAG")
         col_label = cstrcat("M_", filter1, "- M_", filter2);
      end
      if strcmp(filter2(1), "C")
         col_label = filter2;
      end
   end
   if r == 2
      col_label = cstrcat(filter2);
   end
   if r == 4
      col_label = cstrcat("[", char(collab(1)), "-", char(collab(2)), "] - [", char(collab(3)), "-", char(collab(4)),"]");
   end
   title_str = cstrcat('Color-Magnitude Diagram of ', getenv('TARGET_CLUSTER'));
else
   col_label = cstrcat(filter2);
   title_str = cstrcat('Color-Color Diagram of ', getenv('TARGET_CLUSTER'));
end

% Label the vetical axis as the first or second filter mag. depending on the file_name
for i=3:cf
   altind = char(front(i));
   if strcmp(altind, "ALT")
      mlabel = cstrcat(filter2);
      if strcmp(plot_type, "ABSMAG")
         mlabel = cstrcat("M_", filter2);
      end
   end
end


% Colormap
colormap(jet(1024));
symbols = ['bs'; 'ro'; 'gs'; 'md'; 'b>'; 'r<'; 'gp'; 'mo'; 'bs'; 'rd';];

% Plot Mono Magnitude Diagrams of Putative Stars
plot_type = toupper(plot_type);

if strcmp(plot_type, "ABSMAG")
   % Get figure number from environment variable
   figno = str2num(getenv('FIGURE_NUMBER'))

   % Reest every ten plots to avoid number getting too big
   if figno > 5
      setenv('FIGURE_NUMBER', "1");
      figno = 1;
   end

   xpos = 380*(figno-1);
   h = figure(figno)
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk');
   % Clear this plot
   clf(h)
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.15,0.18,0.75,0.73]);
   set(gcf,'position', [xpos,25,400,400]);
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');

   plot(sorted_colors, sorted_mags,  'k+', 'markersize', 2.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 12)
   xlabel(col_label, 'fontsize', 10);
   ylabel(mlabel, 'fontsize', 10);
   %legend(cstrcat("# Obs. = ", int2str(nrows)));
   set(gca,'fontsize',8); % sets font of numbers on axes

   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_mono.jpeg

   % Increment figure number ernvironment variable
   setenv('FIGURE_NUMBER', num2str(figno + 1));
end

if strcmp(plot_type, "MONTAGE")
   % Get figure number from environment variable
   figno = str2num(getenv('FIGURE_NUMBER'))

   % Reest every ten plots to avoid number getting too big
   if figno > 25
      setenv('FIGURE_NUMBER', "1");
      figno = 1;
   end

   xpos = 380*(figno-1);
   h = figure(figno)
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk');
   % Clear this plot
   clf(h)
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.15,0.18,0.75,0.73]);
   set(gcf,'position', [xpos,25,400,400]);
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');

   plot(sorted_colors, sorted_mags,  'k+', 'markersize', 2.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 12)
   xlabel(col_label, 'fontsize', 10);
   ylabel(mlabel, 'fontsize', 10);
   %legend(cstrcat("# Obs. = ", int2str(nrows)));
   set(gca,'fontsize',8); % sets font of numbers on axes

   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_mono.jpeg

   % Increment figure number ernvironment variable
   setenv('FIGURE_NUMBER', num2str(figno + 1));
end

% Plot Mono Magnitude Diagrams of Putative Stars
plot_type = toupper(plot_type);

if strcmp(plot_type, "MONO")
   h = figure(2);
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk')
   % Clear this plot
   clf(h)
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.15,0.08,0.81,0.90]);
   % Publishing size
   % set(gcf,'position', [200,25,300,450]);
   % For Scrreen viewing
   %
   %set(gcf,'position', [200,75,800,800]);
   %
   % Bellazzini 5:11 aspect ratio
  % set(gcf,'position', [200,75,375,715]);
   
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   plot(sorted_colors, sorted_mags,  'k+', 'markersize', 2.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 18)
   xlabel(col_label, 'fontsize', 18);
   ylabel(mlabel, 'fontsize', 18);
   legend(cstrcat("N = ", int2str(nrows)))
   set(gca,'fontsize',18); % sets font of numbers on axes
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_mono.jpeg
end

if strcmp(plot_type, "SHAPE")
   h = figure(2);
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'gnuplot')
   % Clear this plot
   clf(h)
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.15,0.08,0.81,0.90]);
   % Publishing size
   % set(gcf,'position', [200,25,300,450]);
   % For Scrreen viewing
   %
   %set(gcf,'position', [200,75,800,800]);
   %
   % Bellazzini 5:11 aspect ratio
   set(gcf,'position', [200,75,375,715]);
   
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   plot(sorted_colors, sorted_mags,  'b+', 'markersize', 0.25)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 18)
   xlabel(cstrcat("(", col_label, ")_0"), 'fontsize', 14);
   ylabel(cstrcat("M_{", mlabel, "}"), 'fontsize', 14);
   legend(cstrcat("N = ", int2str(nrows)))
   set(gca,'fontsize',14); % sets font of numbers on axes
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
  hold on
  
  % Plot The bellazzini GB Selection Box
  pb = [0.77, -1.99; 0.89,-1.97; 0.67, 3.05; 0.75, 3.53; ...
        0.98,  4.40; 1.86, 6.25; 1.62, 6.80; 1.66, 6.57; ...
        0.75,  7.55; 0.35, 7.50; 0.20, 7.75; -0.14, 7.75; ...
        0.04,  4.70; 0.12, 3.30; 0.20, 2.50; 0.29, 2.22; ...
        0.35,  2.00; 0.41, 1.95; 0.77, -1.99;]; 
        
  [rpb, cpb] = size(pb);
 
  for i=1:rpb-1
    plot([pb(i,1), pb(i+1, 1)], [pb(i, 2), pb(i+1, 2)], "k", "linewidth", 3.0)
  end
  
  hold off
  
  % Save the plot to an jpg file.
  print -Ggswin64c.exe -f2 "-S300,550" -color CMD_shape.jpg
end
  
  
  
  
  
  
% Overplot (red) Magnitude Diagrams of Putative Stars
plot_type = toupper(plot_type);

if strcmp(plot_type, "MAGTRAN")
   color = getenv("PLOT_COLOR");
   if strcmp(color, "")
      color = "r+";
   end 
   
 % Load command parameters
param_file = "Param_Files\\geotran.par";
if exist(param_file, 'file')
  [params, npar] = read_param_file(param_file);
  c0 = str2num(params.TMAG_C0);
  c1 = str2num(params.TMAG_C1);
  c2 = str2num(params.TMAG_C2);
  c3 = str2num(params.TMAG_C3);
  c4 = str2num(params.TMAG_C4);
  c5 = str2num(params.TMAG_C5);
  % Apply magnitude transformation
  tmags1 = sorted_mags .+ c0 .+ (c1.* sorted_colors) .+ (c2 .* sorted_colors.^2);
  tcols = sorted_colors .+ c3 .+ (c4.* sorted_colors) .+ (c5 .* sorted_colors.^2);
 
else
  disp ("No parameter file available. Using command line arguments")
  npar = 0;
end
   
   h = figure(2);
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk')
   % Clear this plot
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.09,0.09,0.83,0.83]);
   % Publishing size
   % set(gcf,'position', [200,25,300,450]);
   % For Scrreen viewing
   set(gcf,'position', [200,75,800,800]);

   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   plot(tcols, tmags1,  color, 'markersize', 2.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 18)
   xlabel(col_label, 'fontsize', 10);
   ylabel(mlabel, 'fontsize', 10);
   %legend(cstrcat("# Obs. = ", int2str(nrows)))
   set(gca,'fontsize',10); % sets font of numbers on axes
  
   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_mono.jpeg
end


% Overplot (red) Magnitude Diagrams of Putative Stars
plot_type = toupper(plot_type);

if strcmp(plot_type, "OVERPLOT")
   color = getenv("PLOT_COLOR");
   if strcmp(color, "")
      color = "r+";
   end 
   h = figure(2);
   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk')
   % Clear this plot
   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'position', [0.09,0.09,0.83,0.83]);
   % Publishing size
   % set(gcf,'position', [200,25,300,450]);
   % For Scrreen viewing
   set(gcf,'position', [200,75,800,800]);

   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   plot(sorted_colors, sorted_mags,  color, 'markersize', 4.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 18)
   xlabel(col_label, 'fontsize', 10);
   ylabel(mlabel, 'fontsize', 10);
   %legend(cstrcat("# Obs. = ", int2str(nrows)))
   set(gca,'fontsize',10); % sets font of numbers on axes
  
   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_mono.jpeg
end


% Plot colorized Color Magnitude Diagrams of Putative Stars
plot_type = toupper(plot_type);

if strcmp(plot_type, "CMD") || strcmp(plot_type, "BOTH")

   star_colors = sorted_colors;

   % Set graphics handle
   h = figure(2);

   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk')
   % Clear this plot
   clf(h)

   fig = gcf ();
   ax = gca();
   set (fig, 'visible', 'on')
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   set(ax, 'position', [0.15,0.15,0.75,0.70]);
   set(gcf,'position', [200,25,600,900]);
   scatter(sorted_colors, sorted_mags, 2.0, star_colors, "filled")

   %plot(sorted_colors, sorted_mags,  's', 'markersize', 3.0)
   axis ([cl, cu, ml, mu]);
   if cstrcat(ext, ".cma") || cstrcat(ext, ".cm3")
      set (gca, 'ydir', 'reverse')
   end
   %title (title_str, 'fontsize', 18)
   xlabel(col_label, 'fontsize', 16);
   ylabel(mlabel, 'fontsize', 16);
   legend(cstrcat("# Obs. = ", int2str(nrows)));
   set(gca,'fontsize',15); % sets font of numbers on axes

   % Save the plot to a JPEG file.
   % print -Ggswin64c.exe CMD_colord.jpeg
end


% Plot Cluster Map
if strcmp(plot_type, "MAP") || strcmp(plot_type, "BOTH")

   star_sizes = 10.*10.^(-0.7.*(sorted_mags-mag_min)/(mag_max-mag_min));
   star_colors = sorted_colors;

   % Set graphics handle
   h = figure(1);

   % Set toolkit to use for this figure (= handle)
   graphics_toolkit(h, 'fltk')
   % Clear this plot
   clf(h)

   figure(3);
   fig = gcf ();
   ax = gca();
   set(ax, 'xminortick', 'on');
   set(ax, 'yminortick', 'on');
   set (fig, 'visible', 'on')
   set(ax,'position', [0.15,0.15,0.75,0.75]);
   set(gcf,'position', [200,25,600,600]);
   scatter(sorted_RAs, sorted_Decs, star_sizes, star_colors, 'filled')
   axis ([RAl, RAu, Decl, Decu]);
   set (gca, 'xdir', 'reverse')
   %title (cstrcat('Map of Star Locations in ', getenv('TARGET_CLUSTER')),'fontsize', 18);
   xlabel('Right Ascension (degrees)','fontsize', 16 );
   ylabel('Declination (degrees)','fontsize', 16);
   legend(cstrcat("# Obs. = ", int2str(nrows)));
   set(gca,'fontsize',15); % sets font of numbers on axes

   % Save the plot to a JPEG file.
   %print -Ggswin64c.exe Map_colord.jpeg
end

t2 = time;

Time_taken = [num2str((t2 - t1)/60), " ", "minutes"]

% Counting stars

if !strcmp(plot_type, "OVERPLOT")
  cstars = 0;
  cstars = yes_or_no("Count stars? ");

  nstars = 0;
  while cstars == 1
    [x,y,button] = ginput();
    nboxes = rows(x)-1;
    if nboxes > 0
      for i = 1: nboxes
        x1  =  min([x(i), x(i+1)]);
        x2  =  max([x(i), x(i+1)]);
        y1  =  min([y(i), y(i+1)]);
        y2  =  max([y(i), y(i+1)]);
    
        c1 = [x1, y1];
        c2 = [x2, y2];

        pindex = find((colors >= x1) & (colors < x2) &...
                    (magnitudes >= y1) & (magnitudes < y2 ));
        [pr, pc] = size(pindex);
        nstars = nstars + pr;
        pr = pr
        figure(2)
        hold on
        plot([x1, x1], [y1, y2], "r", "linewidth", 1.0, [x2, x2], [y1, y2], "r", "linewidth", 1.0, ...
             [x1, x2], [y1, y1], "r", "linewidth", 1.0, [x1, x2], [y2, y2], "r", "linewidth", 1.0)
      end
    end
    cstars = yes_or_no("Count in another region? ");
    disp(cstrcat("Number of stars = ", num2str(nstars, "%7i")))
  end
 end