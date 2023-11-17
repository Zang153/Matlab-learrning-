%% Function is to get required wing loading for given spanwise station (in meters)
% LE HENG LAURENCE LU
% 19.02.2023

function loading = get_wing_loading2(stn_analysis, tail_type)
    %% Section 1: Load wing loading information
    switch tail_type
        case 1 % Main Wing
            load_type = input('1 - Positive lift max loading. 2 - Negtaive lift max loading. Enter: ');
            switch load_type
                case 1
                    load('wing loading 3.75G2.mat')
                    wing_load = wingloading375G;

                case 2
                    load('wing loading -1.5G.mat')
                    wing_load = WingLS;

            end
            
        case 2 % Horizontal Tail
            load('HTcase2.mat')
            wing_load = HTcase2;

        case 3 % Vertical Tail
            load('VTcase2.mat')
            wing_load = VTcase2;
    end
    %% Section 2: Find required row
    idx = 1;
    check = 1;
    [max_row, ~] = size(wing_load);
    while (check == 1)
        if(idx > max_row) % fail safe
            fprintf('Error. Analysis station exceed wing loading index bounds\n')
            return
        end

        if (wing_load(idx, 1) == stn_analysis)
            check = -1;
            break
        end

        if (wing_load(idx, 1) > stn_analysis)
            check = 0;
            break
        end
        
        idx = idx+1;
    end

    %% Section 3: Interpolate
    if (check == -1) % excact value
        loading_temp = wing_load(idx, :);
    elseif (check == 0) % interpolate values
        loading_temp = ones(1, 4);
        loading_temp(1,1) = stn_analysis;
        for (i = 2:4)
            loading_temp(1,i) = interpolater(wing_load(idx-1, i), wing_load(idx,i), wing_load(idx-1, 1), wing_load(idx,1),stn_analysis);
        end
    end

    %% Return values to caller
    % [SF, BM, T]
    loading = loading_temp(2:4);
end