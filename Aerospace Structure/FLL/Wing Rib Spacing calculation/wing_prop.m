%% Sub-function is to calculate the wing section properties
% LE HENG LAURENCE LU
% 23.02.2023

function [c, b1, b2, BM] = wing_prop(dist, desired_stringer_pitch, half_b, root_chord, taper_ratio, taper_dist_start, front_spar, rear_spar, thickness)
    
    % get wing section wing box dimensions
    if (dist <= taper_dist_start)
        section_chord = root_chord;
    else % Assume as a regular trapezium
        temp_b = half_b - taper_dist_start;
        temp_dist = dist - taper_dist_start;
        section_chord = (root_chord * taper_ratio) + (temp_b - temp_dist) * (root_chord - root_chord * taper_ratio) / temp_b;
    end
    
    % Wing box length
    c = rear_spar * section_chord - front_spar * section_chord;
    
    % Wing box height
    b2 = thickness * section_chord;

    % Stringer pitch
    panel_number = c / desired_stringer_pitch;
    panel_number = panel_number - mod(panel_number, 1); % round down to nearest int
    b1 = c / panel_number;

    % get loading data
    % loading_data = [SF, BM, T]]
    loading_data = get_wing_loading(dist/1000); % Stn in wingloading table in [m] but dist in [mm]. <-- sub-fn
    BM = loading_data(2) * 1000; % in [Nmm]. wing loading data in [Nm]

end