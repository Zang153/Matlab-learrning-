x_limit = ( half_b - mod(half_b, 1) ) + 1;
for idxA = 1:length(stringer_area)
    
    figure(idxA)
   
    title_text = ['Stringer Cross Section Area: ', num2str(stringer_area(idxA)), ' mm^2'];

    for idxT = 1:length(tst2)
        plot(L(idxT, :, idxA), tst2(idxT)*ones(length(L(idxT, :, idxA))), '*b');
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