%% This function is to call FARRAR Factor Diagram
function [Factor_ratio] = Call_F_ratio(X,Y)
%   X for A_s/bt 
%   Y for t_s/t

x = [0.4 0.8 1.2 1.6 2.0];
y = 0.4:0.2:1.8;
[a,b] = meshgrid(x,y);

c = [0.75 0.75 0.6 0.5 0.45;0.75 0.85 00.82 0.76 0.72;...
    0.69 0.87 0.88 0.83 0.79;0.6 0.77 0.92 0.92 0.87;...
    0.58 0.77 0.87 0.93 0.94;0.53 0.73 0.82 0.88 0.91;...
    0.5 0.68 0.78 0.83 0.86; 0.49 0.65 0.74 0.78 0.83];

Factor_ratio = interp2(a,b,c,X,Y);


end


