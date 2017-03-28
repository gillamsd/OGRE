function cmd_gen_p(hst_file1, hst_file2, tol, RA_shift, Dec_shift)
%
% Read the files containing the sections of the input data
% and make one array for each.
%

% Save current working directory
currentdir = pwd;
cd (getenv('datdir'))  % Photometry lists are in 'datdir'.

% Some constants
tol = tol / 3600;
zp = 25; % Default IRAF zero point
pix_size = 0.04; % arc sec
ApRad = 0.3; # Photometry apperture radius in Arc sec.
phi = pi*ApRad^2;  # Photometry aperture area
deg2rad = pi/180;

% Default RA and Dec coordinate shifts.
if nargin >=3 && nargin < 5
   RA_shift = 0;
   Dec_shift = 0;
end

% Make header file and split data files from raw phtometry files
[files1, info1] = read_HST_dao_file_p(hst_file1);
[files2, info2] = read_HST_dao_file_p(hst_file2);

% Get Header Information
[Instrument1, Start_time1, Exp_Len1, pRA1, pDec1, Gain1] =  get_header_info_p(info1);
[Instrument2, Start_time2, Exp_Len2, pRA2, pDec2, Gain2] =  get_header_info_p(info2);

% Zeroeth order pointing offset between pictures from camera pointing.
delRA0 = pRA1 - pRA2;             %degrees around the equator
delDec0 = pDec1 - pDec2;          %degress around the equator


% Degress of arc relative to RA=0, Dec=0
% i.e. correct for spherical geometry
pRA1_arc = pRA1*cos(pDec1*pi/180);
pRA2_arc = pRA2*cos(pDec2*pi/180);
delRA0_arc = pRA1_arc - pRA2_arc;

% Add filters to output file names
filter1_1 = info1.Filter_1;
filter1_2 = info1.Filter_2;
filter2_1 = info2.Filter_1;
filter2_2 = info2.Filter_2;

filter1 = filter1_1;
if isspace(filter1_1)
   filter1 = filter1_2;
end

filter2 = filter2_1;
if isspace(filter2_1)
   filter1 = filter2_2;
end

cam = cstrcat("w",Instrument1(length(Instrument1): length(Instrument1)), "w", Instrument2(length(Instrument2): length(Instrument2)));
outFileBase = cstrcat(cam, "_", filter1, "_", filter2);
cmd_file_name = cstrcat(outFileBase, ".cma");
db_file_name =  cstrcat(outFileBase, ".db");

% Number of Julian Days between firt and second epoch
delJD = Start_time1 - Start_time2;

% Pointing difference
RA_offset = delRA0_arc*3600;
Dec_offset = delDec0*3600;


% Make an operations list.
[nfiles1, nchars] = size(files1);
[nfiles2, nchars] = size(files2);

% Delete output CMD and SIMBAD database query files, and open new ones for appending
if exist(cmd_file_name, 'file')
      delete(cmd_file_name);
end
star_cmd = fopen(cmd_file_name, 'a');

if exist(db_file_name, 'file')
      delete(db_file_name)
end
star_db = fopen(db_file_name, 'a');

% Match all data files with all other data files
size_of_cmd_array = 0;

t1 = time;

% Initialize histogram accumulator
histogram_array = [];

for i = 1:nfiles1

%file1 = char(files1(i,:))

   data1 = load_data_p( files1(i,:) );

   [ID1, RA1, Dec1, TotMag1, CIs1] = get_cmd_data_p(data1, "n", Instrument1);
   nr1 = length(RA1);

   RA_arc1 = acos( ( sin(Dec1*deg2rad) ).^2 .+ ( cos(RA1*deg2rad) .* ( cos(Dec1*deg2rad) ).^2 ) ) / deg2rad;

   for j = 1:nfiles2

 %  file2 = files2(j,:);

      data2 = load_data_p( files2(j,:) );
      [ID2, RA2, Dec2, TotMag2, CIs2] = get_cmd_data_p(data2, "n", Instrument2);
      nr2 = length(RA2);

      RA_arc2 = acos( ( sin(Dec2*deg2rad) ).^2 .+ ( cos(RA2*deg2rad) .* ( cos(Dec2*deg2rad) ).^2 ) ) / deg2rad;

      rep_RA_arc1 = repmat(RA_arc1', nr2, 1);
      rep_RA1 = repmat(RA1', nr2, 1);
      rep_Dec1 = repmat(Dec1', nr2, 1);
      rep_TM1 = repmat(TotMag1', nr2, 1);
      rep_CIs1 = repmat(CIs1', nr2, 1);

      rep_RA_arc2 = repmat(RA_arc2, 1, nr1);
      rep_RA2 = repmat(RA2, 1, nr1);
      rep_Dec2 = repmat(Dec2, 1, nr1);
      rep_TM2 = repmat(TotMag2, 1, nr1);
      rep_CIs2 = repmat(CIs2, 1, nr1);

      repColor = rep_TM1 - rep_TM2;

      % Find closest postional matches within acircle radius tol seconds of arc
      % RA_arc = RA * cos(dec)

      % Differences in Right Ascensions (in degrees of right asenscion)
      delta_RAs = rep_RA2 - rep_RA1;
      sgns = sign(delta_RAs);


      %                                                  + (RA2, Dec2)
      %                                                  |
      %                                                  | dtheta(Dec) = dDec
      %                    dTheta(RA)                    |
      %   +----------------------------------------------o   <------ Small arc of a lattitude circle
      % (RA1, Dec1)
      %

      % Convert to the arc-length (degrees of arc) on the lattitude circle at declinations delta1

      cdelRAs =  (sin(rep_Dec1*deg2rad)).^2 .+  ( cos(delta_RAs*deg2rad) .* ( cos(rep_Dec1*deg2rad) ).^2 );
      cdelRAs = min(cdelRAs, 1);
      cdelRAs = max(cdelRAs, -1);
      del_RAs =  acos (cdelRAs) / deg2rad;



      % Add estimate of shift, degree of arc,  from Epoch 1 to Epoch 2
      del_RAs =  sgns.*del_RAs .- RA_shift/3600 ;

      %del_RAs = (rep_RA_arc2 - rep_RA_arc1) .- RA_shift/3600 ;
      del_Decs = (rep_Dec2 - rep_Dec1) .- Dec_shift/3600;

      % Every cross-matched detection
      index = find((del_RAs.^2 + del_Decs.^2) < tol^2);

      % Outputs
      % Keepers - objects marked with flags = 0
      stars_mag = rep_TM1(index);
      stars_color = repColor(index);
      stars_RA1 = rep_RA1(index);
      stars_Dec1 = rep_Dec1(index);
      stars_d_RA = del_RAs(index);
      stars_d_Dec = del_Decs(index);
      stars_RA2 = rep_RA2(index);
      stars_Dec2 = rep_Dec2(index);
      stars_CIs1 = rep_CIs1(index);
      stars_CIs2 = rep_CIs2(index);

      % Output CMD files
      stars_CMD = [stars_RA1, stars_Dec1, stars_color, stars_mag, (stars_d_RA.^2 + stars_d_Dec.^2).^.5, ...
      stars_d_RA, stars_d_Dec, stars_RA2, stars_Dec2 ];

      % Check: Second epoch coordinates adjusted to first epoch system.
      RAsh_deg = RA_shift/3600;
      sgn2 = sign(RAsh_deg);

      cdRA = (cos(RAsh_deg*deg2rad).- (sin(stars_Dec1*deg2rad).^2) ) ./ (cos(stars_Dec1*deg2rad)).^2;
      % Deal with numberical error leading to argument of acos > 1
      cdRA = min(cdRA, 1);
      cdRA = max(cdRA, -1);
      delta_RA = sgn2*acos(cdRA)/deg2rad;

      %delta_RA = acos( (cos(RAsh_deg*deg2rad).- (sin(stars_Dec1*deg2rad).^2) ) ./ (cos(stars_Dec1*deg2rad)).^2 ) / deg2rad

      stars_RA2 = rep_RA2(index) .- delta_RA;
      stars_Dec2 = stars_Dec2 .- Dec_shift/3600;

      stars_CMD = [stars_CMD, stars_RA2, stars_Dec2, stars_CIs1, stars_CIs2];

      % Accumulate shifts in RA and Dec.

      histogram_array = [histogram_array; [stars_d_RA*3600, stars_d_Dec*3600]];

      fprintf (star_cmd, "%3.7f %3.7f %3.3f %3.3f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.3f %3.3f\n", stars_CMD');
      fprintf (star_db, "%3.7f %3.7f \n", [stars_RA1, stars_Dec1]');

      % Accumulate matching statistics
      size_of_cmd_array =  size_of_cmd_array + length(rep_TM1(index));

      % Clear the large 2D arrays out of memory
      clear rep_RA_arc2;
      clear rep_RA2;
      clear rep_Dec2;
      clear rep_TM2;
      clear rep_CIs2;
      clear rep_RA_arc1;
      clear rep_RA1;
      clear rep_Dec1;
      clear rep_TM1;
      clear rep_CIs1;
      clear repColor;
      clear del_Decs;
      clear del_RAs;
      clear data1;
      clear data2;
      clear ID2;
      clear RA2;
      clear Dec2;
      clear RA_arc2;
      clear TotMag2;
      clear CIs2;
      clear star_CIs2;
   end
   % Clear arrays
   clear ID1;
   clear RA1;
   clear Dec1;
   clear RA_arc1;
   clear TotMag1;
   clear CIs1;

   clear stars*

   % delete intermediate data files
   delete( files1(i,:));

end

[RADec_shift] = median( histogram_array);
[RADec_errors] = std( histogram_array);

delete(files2);

fclose(star_cmd);
fclose(star_db);

%whos;

% Return to current working directory
cd(currentdir);

shifts_p(histogram_array, 32);
clc
% Statistics
disp ("          ");
disp ('------------------------------------------------');
disp ("        cmd_gen_p.m Statistics          ");
disp ('------------------------------------------------');
disp (cstrcat(hst_file1, ' ---> ', cmd_file_name));
disp (cstrcat(hst_file2, ' ---> ', cmd_file_name));
disp (cstrcat("FILTER1: ", filter1, "  FILTER2: ", filter2));
disp (cstrcat("Epoch 1: ", num2str(Start_time1, "%7.6f"), "  Epoch 2: ", num2str(Start_time2, "%7.6f"), " MJD"));
disp (" ");
disp (cstrcat("Tolerance circle radius = ", num2str(tol*3600, "%1.3f"), " arc. sec."));
disp ("  ");
disp ("Differences from epoch1 to epoch2.");
disp ("e.g [RA(Epoch2) - RA(Epoch1)]");
disp ( cstrcat( "Target difference:"));
disp ( cstrcat( "RA difference  = ", num2str(RA_offset, '%4.3f'), " arc. sec.") );
disp ( cstrcat( "Dec differsnce = ", num2str(Dec_offset, '%4.3f'), " arc. sec.") );
disp ("  ");
disp ( cstrcat( "Coordinate system adjustments:"));
disp ( cstrcat( "RA  =            ", num2str(RA_shift, '%4.7f'), " arc. sec.") );
disp ( cstrcat( "Dec  =           ", num2str(Dec_shift, '%4.7f'), " arc. sec.") );
disp ("  ");
disp ( cstrcat( "Adjusted coordinate differences:"));
disp ( cstrcat( "RA difference  = ", num2str(RADec_shift(1), '%4.7f'), " arc. sec.") );
disp ( cstrcat( "Dec difference = ", num2str(RADec_shift(2), '%4.7f'), " arc. sec.") );
disp ("  ");
disp ( cstrcat( "Apply the following adjustments on the"));
disp ( cstrcat( "next iteration:"));
disp ( cstrcat( "RA adjustment  = ", num2str(RA_shift + RADec_shift(1), '%4.7f'), " arc. sec.") );
disp ( cstrcat( "Dec adjustment = ", num2str(Dec_shift + RADec_shift(2), '%4.7f')," arc. sec.") );
disp ("  ");
disp ( cstrcat( "Number of observations in common = ", int2str(size_of_cmd_array) ) );
disp ('------------------------------------------------');
disp ("          ");

t2 = time;

Time_taken = [num2str((t2 - t1)/60), " ", "minutes"]

end







