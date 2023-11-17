x = 1:10;

y = ones(1, length(x));

z = 100:120;


plot(x, y, '*b')
for idx = 1:length(x)
    text(x(idx), y(idx), ['\uparrow', num2str(z(idx))], 'FontSize', 10, 'VerticalAlignment','top')
end
