%% This code just gets the important wing properties

%% Input data

wing_type = input('Enter wing type. 1 - Main wing; 2 - Horizontal Tail; 3 - Vertical Tail: '); % 1 - main wing; 2 - horizontal tail; 3 - vertical tail

span_dist = input('Spanwise distance from wing root (mm): '); % in [mm]
%% Input data (manual change)
switch wing_type
    case 1 % main wing
        desired_stringer_pitch = 80;
    case 2 % horizontal tail
        desired_stringer_pitch = 100;
    case 3 % vertical tail
        desired_stringer_pitch = 100;
end

if (wing_type > 3 || wing_type < 0)
    disp('[ERROR] Wrong wing type entered.')
    clear all
    return
end

%% Wing parameters
unit_convr = 10^3; % convert to [mm]
switch wing_type
    case 1 % main wing
        key_rib_points = [3.548, 3.548+3.22, 3.548+3.22+3.22]; % in [m]
        key_rib_points = key_rib_points * 1000; % in [mm]

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

    case 2 % horizontal tail
        S_half_wing = (0.5 * 20.20) *unit_convr * unit_convr;
        half_b = 10.7/2 * unit_convr;
        root_c = 2.36 * unit_convr;
        taper_ratio = 0.6;
        sweep_angle = 10; % in [DEG]
        taper_distance_start = 0 * unit_convr;
        % in %-section chord
        front_spar = 0.15;
        rear_spar = 0.7;
        thickness = 0.09;

    case 3 % vertical tail
        S_half_wing = 12.29 *unit_convr * unit_convr;
        half_b = 3.5 * unit_convr;
        root_c = 4.62 * unit_convr;
        taper_ratio = 0.52;
        sweep_angle = 45; % in [DEG]
        taper_distance_start = 0 * unit_convr;
        % in %-section chord
        front_spar = 0.15;
        rear_spar = 0.7;
        thickness = 0.09;
        tail_weight = 712 * 9.81;

end

if (span_dist > half_b)
    disp('[ERROR] Span distnace exceeded wing half span')
    clear all
    return
end

%% Get wing properties
% Distance from wing root
disp(['Distance from root: ', num2str(span_dist), ' mm'])

% Section chord
    if (span_dist <= taper_distance_start)
        section_chord = root_c;
    else % Assume as a regular trapezium
        temp_b = half_b - taper_distance_start;
        temp_dist = span_dist - taper_distance_start;
        section_chord = (root_c * taper_ratio) + (temp_b - temp_dist) * (root_c - root_c * taper_ratio) / temp_b;
    end

    disp(['Section chord: ', num2str(section_chord), ' mm'])

% Wing box length
c = rear_spar * section_chord - front_spar * section_chord;

disp(['Section wing box length (c): ', num2str(c), ' mm'])

% Wing box height
b2 = thickness * section_chord;

disp(['Section wing box height (b2): ', num2str(b2), ' mm'])

% Stringer pitch
panel_number = c / desired_stringer_pitch;
panel_number = panel_number - mod(panel_number, 1); % round down to nearest int
b1 = c / panel_number;

disp(['Desired stringer pitch: ', num2str(desired_stringer_pitch), ' mm'])
disp(['Required number of panels (n): ', num2str(panel_number), ' '])
disp(['Actual stringer pitch (b1 or b): ', num2str(b1), ' mm'])

% Loading data = [SF, BM, T]]
loading_data = get_wing_loading2(span_dist/1000, wing_type); % Stn in wingloading table in [m] but dist in [mm]. <-- sub-fn
disp(['Shear Force: ', num2str(loading_data(1)), ' N'])
disp(['Bending Moment: ', num2str(loading_data(2)), ' Nm'])
disp(['Torque: ', num2str(loading_data(3)), ' Nm'])

disp('***END OF PROGRAM***')

% Clear memory
clear all
