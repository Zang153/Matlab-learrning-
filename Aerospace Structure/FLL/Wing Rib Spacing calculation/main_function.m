%% Script is to calculate the rib spacing as stringer X-section area and thickness vary. 
% *** ATTENTION: CLOSE ALL OTHER PROGRAMS BEFORE RUNNING ***
% ***       ENSURE ON POWER SUPPLY BEFORE RUNNING        ***

% This assumes all stringers are of equal thickness

% LE HENG LAURENCE LU
% 28.02.2023

% Only doing for top surface as compression drives buckling

%% Section 1: Determine stations where ribs MUST to be placed
% [engine 1, engine 2, engine 3]
key_rib_points = [3.548, 3.548+3.22, 3.548+3.22+3.22]; % in [m]
key_rib_points = key_rib_points * 1000; % in [mm]
key_rib_points_idx = 1;

%% Section 2: Wing information (Entered SI units)
unit_convr = 10^3; % convert to [mm]

S_half_wing = (0.5 * 70.8) *unit_convr * unit_convr;
half_b = 31/2 * unit_convr;
root_c = 2.85 * unit_convr;
taper_ratio = 0.53;
sweep_angle = 7; % in [DEG]
taper_distance_start = 3.548 * unit_convr;
% in %-section chord
front_spar = 0.15;
rear_spar = 0.65;
thickness = 0.12;

skin_mat_E = 73850; % in [MPa]
%skin_mat_E = 70*1000;

rib_mat_E = 72800; % in [MPa]
rib_mat_compr_yield = 657.5; % in [MPa]
    
stringer_area = 180:10:400; % range of stringer X-section area (180 - 400 mm^2)
% below stringer_area range --> unable to interpolate from
% FARRAR diagram. Based on SolidWorks 400 mm^2 is max stringer area.
tst2 = 0.6:0.1:1.4;


%% Section 3: Get skin thickness (using thickest portion, from wing root, throughout wing) [CHECKED W/ EXCEL]

desired_pitch = 80; % desired stringer pitch in [mm]
% [box width, actual possible pitch, box height, Bending Moment]
% At wing root
[c, b1, b2, BM] = wing_prop(0, desired_pitch, half_b ,root_c, taper_ratio, taper_distance_start, front_spar, rear_spar, thickness); % <-- sub-fn
t2 = ( ((BM / (c*b2))*b1^2) / (3.62*skin_mat_E) )^(1/3);

%% Section 4: Get rib spacing as fn of stringer X-sect area and thickness ratio [CHECKED W/ EXCEL] and min. rib sizes

% row --> ts/t2. col --> rib location from root. page --> stringer area
max_length = 500;
L = zeros(length(tst2), max_length + 1, length(stringer_area));
R = zeros(length(tst2), max_length, length(stringer_area)); % Arr storing rib thickness in [mm]

% Set up for loops
idxL = 1;
tempL = 0;
arrL = zeros(1, max_length);
dist_root = 0; % rib distance from wing root in [mm]
%crit_dist = 500; % in [mm] crit distance a calculated rib is from a mandatory rib

idxR = 1;
tempR = 0;
arrR = zeros(1, max_length);

% Find rib spacing for each given stringer area and stringer / skin
% thickness ratio. All in [mm].
% Thereafter, find the min. rib thickness at that rib location to prevent
% buckling / yielding.
for idxA = 1:length(stringer_area)
    for idxT = 1:length(tst2)
        while (dist_root < half_b)
            % Get station information
            [c, b1, b2, BM] = wing_prop(dist_root, desired_pitch, half_b, root_c, ...
                taper_ratio, taper_distance_start, front_spar, rear_spar, ...
                thickness);
            % Get sigma0
            sigma0 = (BM / (c*b2)) / t2;
            % Get sigma_ratio from Catchpole diagram
            sigma_ratio = Call_C_ratio((stringer_area(idxA) / (b1*t2)), tst2(idxT)); % sub-fn to get Catchpole diagram values
            % Get Sigma_cr
            sigma_cr = sigma_ratio * sigma0;
            % Get Farrar efficiency factor
            F = Call_F_ratio((stringer_area(idxA) / (b1*t2)), tst2(idxT)); % sub-fn to get FARRAR diagram values
            % Get rib distance from point of analysis
            tempL = ((F / sigma_cr)^2) * (skin_mat_E * (BM / (c * b2))); % rib spacing
            dist_root = dist_root + tempL;
            
            % Check if dist_root is NaN (Not a Number due to issues calling
            % Farrar efficiency factor) --> if so, ignore
            if (isnan(dist_root))
                arrL = zeros(1, max_length);
                
                disp([num2str(tst2(idxT)), ' / ', num2str(tst2(end)), ' and ', ...
                num2str(stringer_area(idxA)), ' mm^2/ ', ...
                num2str(stringer_area(end)), ' mm^2 NaN detected.'])
                
                break
            end

            % Check if dist_root just past to mandatory rib location --> 
            % if so, overwrite 
            % (i.e if dist_root just past engine 1/2/3)
            if (key_rib_points_idx < (length(key_rib_points) + 1))
                if (dist_root > key_rib_points(key_rib_points_idx))
                    dist_root = key_rib_points(key_rib_points_idx);
                    key_rib_points_idx = key_rib_points_idx + 1;
                end
            end
            
            % Check if dist_root exceeded wing span
            if (dist_root >= half_b)
                break
            end

            % Store into temp array
            arrL(idxL) = dist_root;
            
            % Determine min. rib thickness
            eqv_panel_thickness = t2 + stringer_area(idxA) / b1;
            if (idxL == 1)
                rib_spacing = dist_root; % for theoretical rib at root
            else
                rib_spacing = arrL(idxL) - arrL(idxL - 1); % for other ribs
            end
            I = (c * eqv_panel_thickness^3)/12 + (c * eqv_panel_thickness) * (b2 * 0.5)^2; % Second moment of area
            CrushLoad = ((BM^2) * rib_spacing * b2 * eqv_panel_thickness * c) / (2 * rib_mat_E * I^2);
            BucklingCrit_thickness = ( (CrushLoad * b2^2) / (3.62*rib_mat_E*c) )^(1/3); % min. thickness --> buckling
            YieldCrit_thickness = CrushLoad / (rib_mat_compr_yield * c); % min. thickness --> yielding
            [MaxRibThickness, ~] = max([BucklingCrit_thickness, YieldCrit_thickness]);

            arrR(idxR) = MaxRibThickness;
            
            % Set up for next itr
            idxL = idxL + 1;
            idxR = idxR + 1;
        end
        % Determine min. rib thickness for the last rib
        eqv_panel_thickness = t2 + stringer_area(idxA) / b1;
        rib_spacing = half_b - arrL(idxL); % dist from wing tip to last rib
        I = (c * eqv_panel_thickness^3)/12 + (c * eqv_panel_thickness) * (b2 * 0.5)^2; % Second moment of area
        CrushLoad = ((BM^2) * rib_spacing * b2 * eqv_panel_thickness * c) / (2 * rib_mat_E * I^2);
        BucklingCrit_thickness = ( (CrushLoad * b2^2) / (3.62*rib_mat_E*c) )^(1/3); % min. thickness --> buckling
        YieldCrit_thickness = CrushLoad / (rib_mat_compr_yield * c);
        [MaxRibThickness, ~] = max([BucklingCrit_thickness, YieldCrit_thickness]);
        arrR(idxR) = MaxRibThickness;

        % Save to master array
        L(idxT, :, idxA) = [0, arrL]; % incld theoretical wing root rib
        R(idxT, :, idxA) = arrR;
        
        % Set up for next itr
        idxL = 1;
        tempL = 0;
        arrL = zeros(1, max_length);
        dist_root = 0;
        idxR = 1;
        tempR = 0;
        arrR = zeros(1, max_length);
        key_rib_points_idx = 1;
        
        % Display to know progress
        disp([num2str(tst2(idxT)), ' / ', num2str(tst2(end)), ' and ', ...
            num2str(stringer_area(idxA)), ' mm^2/ ', ...
            num2str(stringer_area(end)), ' mm^2 completed.'])
    end
end

%% Section 5: Plot the graphs
    
    multi_plot(stringer_area, tst2, L, R, key_rib_points, half_b) % sub-fn

    disp(['The skin thickness is: ', num2str(t2), ' mm.'])
