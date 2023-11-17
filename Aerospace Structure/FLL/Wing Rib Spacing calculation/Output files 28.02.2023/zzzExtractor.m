%% This code snippet is to extract the 3D matrix, and output it into a 2D matrix

% LE HENG LAURENCE LU
% 28.02.2023

disp(['The 3D matrix raw data of rib spacings / min. rib thickness saved ' ...
    'from the main function ' ...
    'is a series of 2D matrices. ' ...
    'Each 2D matrix (a page) represents a stringer cross sectional area, ' ...
    'storing the stringer thickness / skin thickness ratio against ' ...
    'rib spacing / min. rib thickness values for that given stringer area.' ...
    'Use this tool to extract the required 2D matrix from the 3D matrix' ...
    'for easy Excel data entry'])
    
disp(' ') % paragraph break

data_type = input('Which data to load? 1 - Rib spacings. 2 - Min. rib thickness: ');

switch  data_type
    case 1
        % Load rib spacings
        load('Output_file_raw_data_wo_NaNL.mat')
        mtx_choosen = L;
    case 2
        % Load min. rib thickness
        load('Output_file_raw_data_wo_NaNR.mat')
        mtx_choosen = R;
end

% Check line 34 of main function for spercific stringer area processed
stringer_area = 180:10:400;
disp(['Stringer area varies from ', num2str(stringer_area(1)), ...
    ' mm^2 to ', num2str(stringer_area(end)), ' mm^2 at 10mm^2 intervals'])

choose_stringer_area = input(['Enter stringer area "page" ' ...
    'you wish to extract (mm^2): '] );


for idx = 1: length(stringer_area)
    if(choose_stringer_area == stringer_area(idx))
        break
    end
end
out_mtx = mtx_choosen(:,:,idx);

disp(['Extraction completed! Double click "out_mtx" in the MATLAB Workspace' ...
    ' to access the data'])

