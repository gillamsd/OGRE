function plot_hess_p(cma_name, hess_array, mincolor, maxcolor, minmag, maxmag, xal, yal, cb)
%
% Draw a Hes Plot
%
% INPUTS:
%   hess_array = array containing numbers of stars in color-magnitude bins.
%   cb = display color bar (yes or no).
%   xal = display x-axis labels (yes or no).
%   yal = display y-axis labels (yes or no).
%
% Make up CMD axis lables from file name and the .TXT file names (They contain
% the photometric errors

parts = strsplit(cma_name, ".");
front = char(parts(1));
ext = char(parts(2));
front = strsplit(front, "_");
[rf, cf] = size(front);
filter1 = char(front(2));
filter2 = char(front(3));
corr = "";
if cf >= 3
  corr = char(front(4));
end
[c, r] = size(strsplit(char(front(3)), "-"));

if strcmp(corr, "CORR")
    magtype = "M";  
else
    magtype = "m";
end

if r == 1
   clabel = cstrcat(magtype, "_{", filter1,"}", " - ", "m_{", filter2, "}");
end
if r == 2
   fil2 = strsplit(filter2, "-");
   f1 = char(fil2(1));
   f2 = char(fil2(2));
   clabel = cstrcat(magtype, "_{",f1, "}"," - ", magtype, "_{",f2, "}");
end

% Choose first or second filter depending on if the
% filename contains the string "ALT" or not.

mlabel = cstrcat(magtype, "_{", filter1, "}");
for i=3:cf
   altind = char(front(i));
   if strcmp(altind, "ALT")
      mlabel = cstrcat(magtype, "_{", filter2, "}");
   end
end

% Image display
xrange = [mincolor:0.2:maxcolor];
yrange = [minmag:0.5:maxmag];
cma = colormap(jet(2048));
h = image(xrange, yrange, hess_array);
set(h, 'cdatamapping', 'scaled');
set(gca,'interpreter', "tex");

% Axes label control
  set (gca, 'ydir', 'reverse')
  %title (title_str, 'fontsize', 18)
  
if strcmp(xal, "yes")
  xlabel(clabel, 'fontsize', 16);
end

if strcmp(yal, "yes")
  ylabel(mlabel, 'fontsize', 16, 'rotation', 90);
end
set(gca,'fontsize',15); % sets font of numbers on axes

 % Color-bar control
if strcmp(cb, "yes")
  caxis([min(min(hess_array)), max(max(hess_array))])
  c = colorbar('East');
  labels = {};
  for v=get(c, 'ytick'), labels{end+1} = sprintf('%1.1e ',v); end
  set(c, 'yticklabel', labels);
end