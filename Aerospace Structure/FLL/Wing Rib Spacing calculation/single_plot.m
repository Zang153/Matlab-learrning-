%% This function plots all the rib positions into one single graph
% LE HENG LAURENCE LU
% 24.02.2023

function single_plot(stringer_area, tst2, L, key_rib_points, half_b)


figure(1)
hold on
for idxA = 1:length(stringer_area)
    % Different markers for different stringer area
    switch idxA
        case 1
            marking = 'ob';
        case 2
            marking = '+b';
        case 3
            marking = '*b';
        case 4
            marking = '.b';
        case 5
            marking = 'xb';
        case 6
            marking = 'squareb';
        case 7
            marking = 'diamondb';
        case 8
            marking = '^b';
        case 9
            marking = 'pentagramb';
        case 10
            marking = 'hexagramb';
        case 11
            marking = 'om';
        case 12
            marking = '+m';
        case 13
            marking = '*m';
        case 14
            marking = '.m';
        case 15
            marking = 'xm';
    end
    
    for idxT = 1:length(tst2)
        plot(L(idxT, :, idxA), tst2(idxT)*ones(length(L(idxT, :, idxA))), marking);
    end
end
xline([key_rib_points, half_b], 'r')

xlabel('Spanwise Station (mm)')
ylabel('stringer thickness/skin thickness ratio')

hold off

end