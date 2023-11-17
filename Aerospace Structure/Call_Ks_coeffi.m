function [Ks] = Call_Ks_coeffi(X)
%   X for a/b 

x = 1:0.5:6;
y = [15 9.8 9.2 8.8 8.6 8.4 8.38 8.35 8.3 8.28 8.25];

Ks = interp1(x, y, X);


end

