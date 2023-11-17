%% This function plots all the rib positions into multiple graphs
% One figure = one given stringer cross-sectional area

% LE HENG LAURENCE LU
% 28.02.2023

function multi_plot(stringer_area, tst2, L, R, key_rib_points, half_b)
x_limit = ( half_b - mod(half_b, 1) ) + 1;
for idxA = 1:length(stringer_area)
    figure(idxA)
    title_text = ['Stringer Cross Section Area: ', num2str(stringer_area(idxA)), ' mm^2'];

    for idxT = 1:length(tst2)
        plot(L(idxT, :, idxA), tst2(idxT)*ones(length(L(idxT, :, idxA))), '*b');
        for idxR = 1:nnz(R(idxT, :, idxA)) % indicates min. rib thickness (in mm) for a particular rib
            text(L(idxT, idxR, idxA), tst2(idxT), ['\uparrow', num2str(R(idxT, idxR, idxA))], 'FontSize', 7, 'VerticalAlignment','top')
        end
        % indicates total number of ribs
        text(L(idxT, idxR, idxA), tst2(idxT), ['   ', num2str(idxR), ' ribs'], 'FontSize', 8, 'VerticalAlignment', 'bottom')
        hold on
    end

    xline(key_rib_points, '-r')
    xline(half_b, '-r')
    xlim([0, x_limit])
    ylim([tst2(1) - 0.1, tst2(end) + 0.1])
    xlabel('Spanwise Station (mm)')
    ylabel('stringer thickness/skin thickness ratio')
    title(title_text)

    hold off

end