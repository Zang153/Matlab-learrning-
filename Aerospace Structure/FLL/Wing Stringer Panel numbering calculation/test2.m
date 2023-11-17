%% Variable in disp --> output to console
x = 10;
y = 11;

disp(['X is ', num2str(x), ' Y is ', num2str(y)])

%% Plot figure in each loop iter and adding a legend
figure
hold on
x = 1:5;
for idx = 1:length(x)
    y = idx + 2*x;

    plot(x, y)
    legends{idx} = sprintf('Plot: %d', idx);
    pause(2);
end
legend(legends)
grid on
hold off

%% Test
A = zeros(2, 10);
B = 1:5;
A(1, 1:length(B)) = B