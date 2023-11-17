%% Script is to calculate the rib spacing as stringer X-section area and thickness vary. 
% This assumes all stringers are of equal thickness

% LE HENG LAURENCE LU
% 24.02.2023

% Only doing for top surface as compression drives buckling
% For debugging, using values from Dr Yang's Excel sheet (non-Copy ver.)
%% Section 1: Determine stations where ribs MUST to be placed
% [engine 1, engine 2, engine 3]
key_rib_points = [3.548, 3.548+3.22, 3.548+3.22+3.22]; % in [m]
key_rib_points = key_rib_points * 1000; % in [mm]
key_rib_points_idx = 1;

%% Section 2: Wing information (Entered SI units)
unit_convr = 10^3; % convert to [mm]

half_b = 31/2 * unit_convr;
root_c = 2.85 * unit_convr;
taper_ratio = 0.53;
sweep_angle = 7; % in [DEG]
taper_distance_start = 3.548 * unit_convr;
% in %-section chord
front_spar = 0.15;
rear_spar = 0.65;
thickness = 0.12;

%skin_mat_E = 73850; % in [MPa]
skin_mat_E = 70*1000;

% stringer_area = 10:10:150; % range of stringer X-section area (10 - 150 mm^2)
% tst2 = 0.6:0.1:1.4;
stringer_area = 57.6;
tst2 = 1/1.0473;


%% Section 3: Get skin thickness (using thickest portion, from wing root, throughout wing) [CHECKED W/ EXCEL]

desired_pitch = 80; % desired stringer pitch in [mm]
% [box width, actual possible pitch, box height, Bending Moment]
% At wing root
% [c, b1, b2, BM] = wing_prop(0, desired_pitch, root_c, sweep_angle, taper_distance_start, front_spar, rear_spar, thickness);

c = 1000;
b1 = 50;
b2 = 0.2 * 1000;
BM = 23286.8 * 1000;

t2 = ( ((BM / (c*b2))*b1^2) / (3.62*skin_mat_E) )^(1/3);
%% Section 4: Get rib spacing as fn of stringer X-sect area and thickness ratio

% row --> ts/t2. col --> rib location from root. page --> stringer area
L = zeros(length(tst2), 100, length(stringer_area));

idxL = 1;
tempL = 0;
arrL = zeros(1, 100);
dist_root = 0; % rib distance from wing root in [mm]
crit_dist = 500; % in [mm] crit distance a calculated rib is from a mandatory rib

for idxA = 1:length(stringer_area)
    for idxT = 1:length(tst2)
        while (dist_root < 200)
            % Get station information
%             [c, b1, b2, BM] = wing_prop(dist_root, desired_pitch, root_c, ...
%                 sweep_angle, taper_distance_start, front_spar, rear_spar, ...
%                 thickness);
            % Get sigma0
            sigma0 = (BM / (c*b2)) / t2;
            % Get sigma_ratio from Catchpole diagram
            sigma_ratio = Call_C_ratio((stringer_area(idxA) / (b1*t2)), tst2(idxT));
            % Get Sigma_cr
            sigma_cr = sigma_ratio * sigma0;
            % Get Farrar efficiency factor
            F = Call_F_ratio((stringer_area(idxA) / (b1*t2)), tst2(idxT));
            % Get rib distance from point of analysis
            tempL = ((F / sigma_cr)^2) * (skin_mat_E * (BM / (c * b2)));
            dist_root = dist_root + tempL;
            
            % Check if dist_root close to mandatory rib location --> 
            % if so, overwrite 
            % (i.e if dist_root just slightly past engine 1/2/3)
            if (key_rib_points_idx < 4)
                if (dist_root - key_rib_points(key_rib_points_idx) <= crit_dist ...
                        && dist_root - key_rib_points(key_rib_points_idx) > 0)
                
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
            
            % Set up for next itr
            idxL = idxL + 1;
        end
        % Save to master array
        L(idxT, :, idxA) = arrL;
        % Set up for next itr
        idxL = 1;
        tempL = 0;
        arrL = zeros(1,100);
        dist_root = 0;
        key_rib_points_idx = 1;
    end
end

%% Section 5: Plot the graph
L