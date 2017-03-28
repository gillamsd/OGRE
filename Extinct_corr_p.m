function Extinct_corr_p(cmafile, R, excess, mM0)
%
% Apply xxtinction Corrections, per Cardelli, Clayton and Mathis Ap> J. v345 p245-256, 189 October 1,
% to the magnitudes in CMA file.
% Also subtract the distance modulus to make CMD comparable with theory
%
% The CMA file will be in the sense [M1, C]
% Where, C = M2-M3
% and,   M1 = HLA magnitude in filter 1
%        M2 = HLA magnitude in filter 2
%        M3 = HLA magnitude in filter 3
%
% New M1 = M1 + A1
% New M2 = M2 + A2
% New C = New M2 - New M3
%
%
% Inputs:
% cmafile = A color-magnitude array in a CMA file.
% R = ratio of selective to total extinction                % Guessed
% excess = E(B-V)                                           % measured
% mM0 = distance modulus                                    % Measured
%
% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

% Make up the output and temporary file names, based on the input file name
out_file_base = substr (cmafile, 1, length(cmafile)-4);

cma_out = strcat(out_file_base, '_EXCOR.cma');

% Get the source filters from the CMA file name
front = strtok(cmafile, ".");
base = front;
front = strsplit(front, "_");
[rf, cf] = size(front);
[linst, nchars] = size(char(front(1)));
instrum1 = char(front(1))(nchars-3:nchars-2);
instrum2 = char(front(1))(nchars-1:nchars);
filt1 = char(front(2));    %1st filter of FXXXW_FYYYW sytle CMA file
lf2 = length(char(front(3)));
if lf2 == 5
   filt2 = char(front(2));  %2nd filter of FXXXW_FYYYW sytle CMA file
   filt3 = char(front(3));  %3rd filter of FXXXW_FYYYW sytle CMA file
else
   filt2 = char(front(3))(1:5);  %2nd filter of FXXXW_FYYYW-FZZZW sytle CMA file
   filt3 = char(front(3))(7:11); %3rd filter of FXXXW_FYYYW-FZZZW sytle CMA file
end

% For the FXXXW_FYYYW_ALT style of CMA file name
for i=3:cf
   altind = char(front(i));
   if strcmp(altind, "ALT")
      filt1 = char(front(3));    %1st filter of FXXXW_FYYYW_ALT sytle CMA file
      filt2 = char(front(2));    %2nd filter of FXXXW_FYYYW_ALT sytle CMA file
      filt3 = char(front(3));    %3rd filter of FXXXW_FYYYW_ALT sytle CMA file
   end
end

filt1 = toupper(filt1);
filt2 = toupper(filt2);
filt3 = toupper(filt3);

% Load color-magnitude array (CMA) one
cma = load(cmafile);

% Pivot Wavelengths for HST filters from "Wide Field Camera 3 Instrument Mini-Handbook for Cycle 13
wfc3_filters = ["F225W"; "F275W"; "F300X"; "F336W"; "F390W"; "F438W"; "F555W"; "F606W"; "F625W"; "F814W"]; % Filter name
pivot_lambda = [ 2284.3;  2742.4;  2819.6;  3361.1;  3904.6;  4318.7;  5309.8;  5932.3;  6254.0;  8304.7]; % center wavelength in microns

% Temporary fix for WFPC2 F555W and F814W filters
pivot_lambda = [ 2284.3;  2742.4;  2819.6;  3361.1;  3904.6;  4318.7;  5439.0;  5932.3;  6254.0;  8012.0]; % center wavelength in

% Extinction corrections structure for wfc3 only
nrows_wfc3_filters = rows(wfc3_filters);

for i = 1: nrows_wfc3_filters
   wfc3_lambda.(deblank(wfc3_filters(i,:))) = pivot_lambda(i);
end

lambdas = [wfc3_lambda.(filt1), wfc3_lambda.(filt2),  wfc3_lambda.(filt3)] / 1e4;


Av = R*excess; % Number of magnitudes of extinction in the V band.

% Extinction Corrections per Cardelli, Clayton and Mathis Ap> J. v345 p245-256, 189 October 1
% Get extinction ratio for each filter
ext_ratios = extinction_p(lambdas, R, excess);  % Extinction ratios at wavelengths of  M1 and M2 filters

% Calculate E(lamda2 - lambda3)
A_lambda = ext_ratios(2, :)*Av; % Magnitudes of extinction in the M1,  M2 and M3 filters
E = A_lambda(2) - A_lambda(3);  % E(filter2 - filter3)

% Number of observations in the output files
[rows, cols] = size(cma);

% Apply extinction corrections to cols. 3 and 4 of the CMA
corrections = [0, 0, -E, -A_lambda(1)- mM0, 0, 0, 0, 0, 0, 0, 0, 0, 0](:, 1:cols); % Corrections

cma = cma + corrections; % Deredenned filter1 and color

if exist('temp.ec', 'file')
   delete('temp.ec')
end

% Save the output files
dlmwrite ('temp.ec',  cma, "precision", "%3.7f", "delimiter", " ");
[status, msg, msgid] = copyfile ('temp.ec', cma_out);

delete('temp.ec');


% Return to original directory
cd (currentdir);


% Statistics
disp ("          ");
disp ('-----------------------------------------------------------');
disp ("                Extinct_corr_p.m Statistics       ");
disp ('-----------------------------------------------------------');
disp (cstrcat("INPUT FILE:  ", cmafile));
disp (cstrcat("OUTPUT FILE: ", cma_out));
disp ("          ");
disp (cstrcat('Color [', filt2, '-', filt3,']'));
disp (cstrcat('Filter1 ', filt1, " (", num2str(ext_ratios(1,1)), " Angstroms)"));
disp (cstrcat('Filter2 ', filt2, " (", num2str(ext_ratios(1,2)), " Angstroms)"));
disp (cstrcat('Filter3 ', filt3, " (", num2str(ext_ratios(1,3)), " Angstroms)"));
disp ("          ");
disp ("Assumed parameters");
disp (cstrcat("(m=M)o = ",num2str(mM0, "%3.3f"), ' mag'));
disp (cstrcat("R = ",num2str(R, "%3.2f")));
disp (cstrcat("E(B-V) = ", num2str(excess, "%3.3f"), ' mag'));
disp (cstrcat("A(V) = R*E(B-V)= ", num2str(Av, "%3.3f"), ' mag'));
disp ("          ");
disp (cstrcat("Extinction corrections applied:"));
disp (cstrcat("To ", filt1," = ", num2str(-A_lambda(1),"%3.3f"), " mag"));
disp (cstrcat("To ", filt2," = ", num2str(-A_lambda(2),"%3.3f"), " mag"));
disp (cstrcat("To ", filt3," = ", num2str(-A_lambda(3),"%3.3f"), " mag"));
disp (cstrcat('To [',filt2, '-', filt3,'] = ', num2str(-E, "%3.3f"), ' mag'));
disp ('---------------------------------------------------------');
disp ("          ");


