clc;
clear all;

%%
load("WingGeo.mat");
WingGeo = table2array(WingGeo);

%%  Material properties
%   Al 2024
%   Young's Modulus in Pa 
%   Posion ratio
%   Shear Modulus in Pa
%   Bulk Modulus in Pa
E = 70*10^9;
mu = 0.3;
G = 27*10^9;
K = 58*10^9;

%%
half_b = 31/2;
root_c = 2.85;
taper_ratio = 0.53;
sweep_angle = 7; % in [DEG]
taper_distance_start = 3.548;
front_spar = 0.15;
rear_spar = 0.65;
thickness = 0.12;

figure
line((0:half_b),1,"*")

%%  Essential Definition for wing
%   Web height vary from root to tip
WebHeight = zeros(size(WingGeo(:,1)));
WebHeightRatio = 0.09;                      %   Change later
WebHeight = WingGeo(:,4).*WebHeightRatio; 