%% Script is to calculate the optimal stringer panal number. This assumes equally spaced stringers
% LE HENG LAURENCE LU
% 21.02.2023

% Only doing for top surface as compression drives buckling

%% Section 1: Determine stations to perform analysis
% [Wing root, engine 1, engine 2, engine 3]
key_stn_points = [0, 3.548, 3.548+3.22, 3.548+3.22+3.22]; % in [m]

analysis_point = input("Input station of analysis. 1 - Wing Root. 2 - Engine 1. 3 - Engine 2. 4 - Engine 3. ");
stn_analyse = key_stn_points(analysis_point);

%% Section 2A: Wing information (all in SI units)
half_b = 31/2;
root_c = 2.85;
taper_ratio = 0.53;
sweep_angle = 7; % in [DEG]
taper_distance_start = 3.548;
front_spar = 0.15;
rear_spar = 0.65;
thickness = 0.12;

skin_mat_E = 73850 * 10^6; % in [Pa]

%% Section 2B: Wing information at analysis point (all in SI units)
% section chord
if (stn_analyse <= taper_distance_start)
    section_chord = root_c;
else
    section_chord = root_c - (stn_analyse - taper_distance_start) * tan(sweep_angle*pi/180);
end

% Wing box width
c = rear_spar * section_chord - front_spar * section_chord;
%c = 1; % For debugging purposes

% Wing box height
b2 = thickness * section_chord;


%% Section 2C: Get wing loading data
% loading data = [SF, BM, T]
loading_data = get_wing_loading(stn_analyse); % <-- calling sub-fn

BM = loading_data(2);

%% Section 3A: Run analysis - Control, no stringers [CHECKED W/ EXCEL]
% Compressive load
N0 = BM / (c*b2);
  
% Skin thickness without stringers
t2_0 = ((N0*c^2)/(3.62*skin_mat_E))^(1/3);

% Skin X-section area
A_skin0 = c*t2_0;

%% Section 3B: Run analysis - With Stringers [CHECKED W/ EXCEL]
% Weight of stringer + skin only depends on X-sectional area of stringer.
%
% Catchpole diagram depends on stringer thickness but this can be
% determined later once stringer spacing value is obtained.
% Noting that stringer spacing value has to be constant throughout wing
% span therefore, number of panels changes from root to tip
% 
% Since stringer width and flange length is linked by a ratio, with
% determined stringer thickness and found stringer X-sectional area -->
% able to get overall dimensions of stringer such that 
% 2 * stringer flange < stringer pitch distance

% number of panels = n. number of stringers = (n - 1)
n = 2:24;
% assumed stringer thickness
stringer_thickness = 1 * 10^-3; % in [m]
% stringer flange to web ratio
w = 0.3; % Assumed
% min X-section area of stringer (stringer height must be > 1/w --> about 4mm --> d = 1.2mm)
min_stringer_area = (1.2*10^-3)*stringer_thickness*(2 + 1/w);
% resolution
delta_stringer_area = 1 * 10^-6; % delta = 1 mm^2

% stringer pitch distance
b = zeros(1, length(n));
% skin thickness
t2 = zeros(1, length(n));
% sigma_o
sigma_o = zeros(1, length(n));
% effective panel area
max_col_skin = 1000;
A_stringer = zeros(length(n), max_col_skin);
A_skin_efft = zeros(length(n), max_col_skin);
% weight saving
weight_saving = zeros(length(n), max_col_skin);


for idx = 1:length(n) % looping through by number of stringers
    
    % determine stringer pitch distance
    b(idx) = c/n(idx);

    % determine max X-section area of stringer
    max_stringer_area = (b(idx)/2)*stringer_thickness*(2 + 1/w);

    % determine skin panel thickness
    t2(idx) = ((N0*b(idx)^2)/(3.62*skin_mat_E))^(1/3);

    % determine sigma_o
    sigma_o(idx) = N0 / t2(idx);

    % determine effective panel area --> determine weight saving
    % ( weight proportional to volume proportional
    % to X-sectional area )
    A_stringer_range = min_stringer_area:delta_stringer_area:max_stringer_area;
    A_stringer(idx, [1:length(A_stringer_range)]) = A_stringer_range;
    for idx2 = 1:length(A_stringer_range)
        A_skin_efft(idx, idx2) = n(idx) * b(idx) * (t2(idx) + A_stringer_range(idx2) / b(idx));
        weight_saving(idx, idx2) = abs(A_skin0 - A_skin_efft(idx, idx2) ) / A_skin0 * 100;
    end

    % Display highest weight saving
    [max_weight_saved, idx2] = max(weight_saving(idx, :));
    disp(['For ', num2str(n(idx)), [' panels, ' ...
        'max weight saving of '], num2str(max_weight_saved), ...
        ' % , achieved with ', num2str(t2(idx)*10^6), ' mm^2 skin area, ', ...
        num2str(A_stringer_range(idx2)*10^6), ' mm^2 stringer area.\n'])    

end

%% Section 4A: Plot Weight saving - Stringer area graph

figure(1)
hold on
for idx = 1:length(n)
    % count number of non-zero elements in the row
    non_zero_ele = nnz(A_stringer(idx,:));

    % Plot graph
    plot(A_stringer(idx,[1:non_zero_ele])*10^6, weight_saving(idx, [1:non_zero_ele]))
    legends{idx} = sprintf('%d panels', n(idx));

    pause(0.2)

end

legend(legends)
xlabel('Stringer Cross-Section area (mm^2)')
ylabel('Weight Saving (%)')
grid on
hold off

%% Section 4B: Plot Weight saving - Effective panel area graph

figure(2)
hold on
for idx = 1:length(n)
    % count number of non-zero elements in the row
    non_zero_ele = nnz(A_stringer(idx,:));

    % Plot graph
    plot(A_skin_efft(idx,[1:non_zero_ele])*10^6, weight_saving(idx, [1:non_zero_ele]))
    legends{idx} = sprintf('%d panels', n(idx));

    pause(0.2)

end

legend(legends)
xlabel('Effective Panel Cross-Section area (mm^2)')
ylabel('Weight Saving (%)')
grid on
hold off

%% Section 5: Output to table
A_stringer_ouput = A_stringer .* 10^6; % in [mm^2]
A_skin_efft_output = A_skin_efft .* 10^6; % in [mm^2]
weight_saving_output = weight_saving;
% stringer_pitch_output = b' * 10^3; % in [mm]
Thickness_skin_output = t2' * 10^3; % in [mm]