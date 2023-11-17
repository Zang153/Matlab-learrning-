%% Test rib spacing, if 1mm stringer thickness --> min rib spacing [all units in mm or MPa]
t2 = 1.0473;
ts = 1;
As = 57.6;
b = 50;
%N = 116434/1000;
N = 100000/1000;
sigma0 = N / t2;
E = 70*10^3;

ratio_t = ts/t2;
ratio_b = As / (b*t2);

Catchpole_ratio = Call_C_ratio(ratio_b, ratio_t);
F_Factor = Call_F_ratio(ratio_b, ratio_t);

sigma_cr = sigma0 * Catchpole_ratio;

L = 100;
delta_L = 0.1;
crit = 1;
check = 0;
idx = 1;

while (check == 0)
    sigma_F = F_Factor * sqrt((N*E)/L);
    
    diff = sigma_cr - sigma_F;
    if (diff > 0)
        check = 1;
        break
    end
    if (idx >= 10^5)
        check = 1;
        disp('Unable to converge')
        break
    end
    L = L + delta_L;
    idx = idx + 1;
end

disp(['Final est. rib spacing is: ', num2str(L), ' mm'])

%% Test rib spacing, if 1mm stringer thickness --> min rib spacing [all units in mm or MPa]
t2 = 1.0473;
ts = 1;
As = 57.6;
b = 50;
N = 116434/1000;
sigma0 = N / t2;
E = 70*10^3;

% ratio_t = ts/t2;
ratio_t = 0.6:0.1:1.4;
ratio_b = As / (b*t2);

final_L = zeros(1, length(ratio_t));


for i = 1:length(ratio_t)
    L = 100;
    delta_L = 0.1;
    crit = 1;

    check = 0;
    idx = 1;

    Catchpole_ratio = Call_C_ratio(ratio_b, ratio_t(i));
    F_Factor = Call_F_ratio(ratio_b, ratio_t(i));

    sigma_cr = sigma0 * Catchpole_ratio;

    while (check == 0)
        sigma_F = F_Factor * sqrt((N*E)/L);
        diff = sigma_cr - sigma_F;
        if (diff > 0)
            check = 1;
            break
        end
        if (idx >= 10^5)
            check = 1;
            disp('Unable to converge')
            break
        end
        L = L + delta_L;
        idx = idx + 1;
    end

    final_L(i) = L;
end

plot(ratio_t, final_L)
xlabel('ts/t ratio')
ylabel('Rib spacing (mm)')
