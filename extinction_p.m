function ratios = extinction_p(lambda, R)
%
% Calcultes the ratio A(lambda)/A(V) per Cardelli, Clayton and Mathis Ap> J. v345 p245-256, 189 October 1
%
% Inputs:
% lambda = array of effective filter wavelengths in micrometers
% R is formally  A(V)/E(B-V)
%
% Outputs:
% An array containing the input wavelengths (top row) and corresponding extinction ratios (bottom row).

% Convert wavelengths to inverse wavelengths
xin = 1./lambda;  % Units 1/micrometers

% Sort input array into IR, Optical/NIR and UV/FarUV subsets

ind = find(xin<1.1);
i = find(xin(ind)>=0.29);
indIR = ind(i);           % Indices of IR subset

ind = find(xin<3.3);
i =  find(xin(ind)>=1.1);
indOpt = ind(i);          % Indices of optical/NIR subset

ind = find(xin<8.0);
i = find(xin(ind)>=3.3);
indUV = ind(i);            % Indices of UV subset

ind = find(xin<=10.0);
i = find(xin(ind)>=8.0);
indFUV = ind(i);            % Indices of Far-UV subset

ind1 = find(xin<0.29);        % Indices of item outside range of extinction functions
ind2 = find(xin>10.0);

[rIR, cIR] = size(indIR);
[rOpt, cOpt] = size(indOpt);
[rUV, cUV] = size(indUV);
[rFUV, cFUV] = size(indFUV);

[r1, c1] = size(ind1);
[r2, c2] = size(ind2);

% Basic relation

% Coefficients a(x) and b(x)

a = [];  % For collecting the a(x) coeffs. in the order that they are calc'ed.
b = [];  % For collecting th b(x) coeffs. in the order that they are calc'ed.
lindex = []; % For collecting indices of wavelengths in the order that they are calc'ed.

% Infrared (0.29 um-1 to 1.1um-1, 3.333 um to .0909 um)
if cIR > 0
   x = xin(indIR);
   a = [a, 0.574*x.^1.61];
   b = [b, -0.527*x.^1.61];
   lindex = [lindex, indIR];
end

% Optical/NIR (1.1u um-1 to 3.3 um-1, 0.0909 um to 0.30303 um)
if cOpt > 0
   x = xin(indOpt);
   y = x - 1.82;
   a = [a, 1 + 0.17699.*y  - 0.50447.*y.^2 - 0.02427.*y.^3  + 0.72085.*y.^4 + 0.01979.*y.^5  - 0.77530.*y.^6 + 0.32999.*y.^7];
   b = [b, 1.41338.*y + 2.28305.*y.^2  + 1.07233.*y.^3  - 5.38434.*y.^4 - 0.62251.*y.^5  + 5.30260.*y.^6   - 2.09002.*y.^7];
   lindex = [lindex, indOpt];
end

% UV (3.3 um-1 to 8 um-1, 0.30303 um to 0.125 um)
if cUV > 0
   x = xin(indUV);
   Fa = -0.04473.*(x - 5.9).^2 - 0.009779.*(x-5.9).^3;
   Fb = 0.2130.*(x-5.9).^2 + 0.1207.*(x-5.9).^3;
   next_a = 1.752- 0.316.*x- 0.104./((x- 4.67).^2 + 0.341) .+ Fa;
   a = [a, next_a];
   next_b = -3.090 + 1.825.*x + 1.206./((x - 4.62).^2  + 0.263) .+  Fb;
   b = [b, next_b];
   lindex = [lindex, indUV];
end

% FUV (8.0 to 10.0 um-1)
if cFUV > 0
   x = xin(indFUV);
   a = -1.073 - 0.628.*(x-8) + 0.137.*(x-8).^2 - 0.070.*(x-8).^3;
   b = 13.670 + 4.257.*(x-8) - 0.420.*(x-8).^2 + 0.374.*(x-8).^3;
   lindex = [lindex, indFUV];
end

if c1 > 0
   for j = 1:c1
      a = [a, -999];
      b = [b, 0];
   end
   lindex = [lindex, ind1];
end

if c2 > 0
   for j = 1: c2
      a = [a, 999];
      b = [b, 0];
   end
%   wl = [wl, lambda(ind2)];
   lindex = [lindex, ind2];
end


% Sort extinction ratios in order of wavelength
[s, i] = sort(lindex);
k = a + b./R;  % a + b./R = <A(lambda)/A(V)>

ratios = [lambda; k(i)];


