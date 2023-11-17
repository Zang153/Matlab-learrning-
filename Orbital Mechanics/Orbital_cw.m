% Coursework for AERO70016/ AERO97079 Orbital Mechanics 
% Date: March 17th, 2023
% Programmed by: YX Zang
% Script: Orbital_cw.m
% This script mainly defines two functions to achieve Cowell's method and
% Encke's method.
% Reference book, "Orbital Mechanics for Engineering Students 3nd Edition"
% by Howard D. Curtis
 
clc
clear;
close all;
% format long;

% Define the parameters
year = 365.25;                                                  % in days
yearinsecond = 86400 * year;                        % in seconds per year
au = 149597870700 / 1000;                                         % in km
t0 = 0;
tf = 20*yearinsecond;

r_init = [-1.05 * au, 0, 0]';                                     % in km
v_init = [0, -6.1316 * au / yearinsecond, 0]';                  % in km/s

mu = 39.4769 * au^3 / yearinsecond^2;           % gravitational parameter

at_init = 1 / 3 * 10^(-7);                                    % in km/s^2

% Part 1 
%-------------------------------------------------------------------------%
x_init = [r_init; v_init]';

options = [odeset('MaxStep',tf/1000)];

[~, un] = ode45(@(t, x) unpertubed(x, mu), [t0 tf], x_init, options);

[t2, x] = ode45(@(t,x) cowell(x, mu, au, at_init), [t0 tf], x_init, ...
    options);

Cowell_final = x(end,:);
display(Cowell_final);

plotQ1(un, x);

% Part 2
%-------------------------------------------------------------------------%
R0 = r_init;
V0 = v_init;

del_t = tf/500;
options = odeset('MaxStep',del_t);

t1 = t0;
tsave = t0;
y = [R0; V0]';
del_y0 = zeros(6,1);

t1 = t1 + del_t;

while t1 < tf+del_t/2
        
        % calculate delta
        [dum, z] = ode45(@(t, x) encke(t, x, R0, V0, t0, at_init, au, ...
            mu), [t0 t1], del_y0, options);
        
        % calculate osculating orbit
        [Rosc, Vosc] = rv_from_r0v0_dt(R0, V0, t1 - t0, mu);
        
        % Update
        R0 = Rosc + z(end, 1:3)';
        V0 = Vosc + z(end, 4:6)';
        t0 = t1;

        tsave = [tsave;t1]; %#ok
        y = [y; [R0;V0]']; %#ok
        t1 = t1 + del_t;
        
end

Encke_final = y(end,:);
display(Encke_final);
plotQ2(y, x);

% Part 3
%-------------------------------------------------------------------------%
analysis_a_e(y, mu);

m = [y(10,:);y(125,:);y(250,:);y(375,:);y(500,:)];

draw_osculating_orbit(m, mu);
hold on;
plot(x(:,1), x(:,2), 'LineWidth', 2, 'Color', 'magenta');

% P4 Functions
%-------------------------------------------------------------------------%
function dxdt = unpertubed(x, mu)
    
    r = x(1:3);
    v = x(4:6);

    rMag = norm(r);
    
    a = -mu*r/rMag^3;

    dxdt = [v; a];

end

function dxdt = cowell(x, mu, au, at_init)
    
    r = norm([x(1) x(2) x(3)]);
    v = norm([x(4) x(5) x(6)]);
   
    vx = x(4);
    vy = x(5);
    vz = x(6);
    ax = -mu * x(1) / r^3 + at_init * (au / r)^2 * x(4) / v;
    ay = -mu * x(2) / r^3 + at_init * (au / r)^2 * x(5) / v;
    az = -mu * x(3) / r^3 + at_init * (au / r)^2 * x(6) / v;

    dxdt = transpose([vx vy vz ax ay az]);
    
end

function dxdt = encke(t, x, r0, v0, t0, at_init, au, mu)
    
    del_r = x(1:3);
    del_v = x(4:6);
    
    dt = t - t0;
    [Rosc, Vosc] = rv_from_r0v0_dt(r0, v0, dt, mu);

    Rpp = Rosc + del_r;
    Vpp = Vosc + del_v;

    rosc = norm(Rosc);
    rpp = norm(Rpp);
       
    % F = 1 - (rosc/rpp)^3;
    q = dot(del_r, del_r - 2.*Rpp)/rpp^2;
    F = - q*(3 + 3*q +q^2)/(1 + (1 + q)^1.5);
    ap = at_init.*(au./rpp).^2.*(Vpp/norm(Vpp));
    del_a = -mu/rosc^3*(del_r - F*Rpp) + ap;
    dxdt = [del_v; del_a];

end

function [r, v] = rv_from_r0v0_dt(r_init, v_init, dt, mu)
%{
  This function computes the state vector (r, v) from the initial
  state vector (r0, v0) and the change in time.

  mu  - gravitational parameter (km^3/s^2)
  dt  - change in time (s)
  r0  - position vector at time t0 (km)
  v0  - velocity vector at time t0 (km/s)
  r   - final position vector (km)
  v   - final velocity vector (km/s)
%}   
    if dt ~= 0
        coe_init = coe_from_sv(r_init, v_init, mu);

        TA_init = coe_init(6);
        e_init = coe_init(2);
        a_init = coe_init(7);
        n_init = sqrt(mu/a_init^3);
    
        dTA = cal(n_init, dt, TA_init, e_init);
        [r, v] = rv_from_r0v0_ta(r_init, v_init, dTA, mu);
    else
        r = r_init;
        v = v_init;
    end
end

function [coe, E] = coe_from_sv(R, V, mu)
%{
  This function computes teh classical orbital elements (coe)
  from the state vector (R, V).

  mu  - gravitational parameter (km^3/s^2)
  R   - position vector in the geocentric equatorial frame (km)
  V   - velocity vector in the geocentric equatorial frame (km)
  r,v - the magnitudes of R and V
  vr  - radial velocity component (km/s)
  H   - the angular momentum vector (km^2/s)
  h   - the magnitude of H (km^2/s)
  incl - inclination of the orbit (rad)
  N   - the node line vector (km^2/s)
  n   - the magnitude of N
  cp  - cross product of N and R
  RA  - right ascension of the ascending node (rad)
  E   - eccentricity vector 
  e   - eccentricity (magnitude of E)
  eps - a small number below which the eccentricity is considered to be
        zero
  w   - argument of perigee (rad)
  TA  - true anomaly (rad)
  a   - semimajor axis (km)
  pi  - 3.1415926...
  coe - vector of orbital elements [h e RA incl w TA a]
%}

eps = 1.e-10;

r = norm(R);
v = norm(V);

vr = dot(R,V)/r;

H = cross(R,V);
h = norm(H);

incl = acos(H(3)/h);

N = cross([0 0 1],H);
n = norm(N);

if n ~= 0
    RA = acos(N(1)/n);
    if N(2) < 0 
        RA = 2*pi - RA;
    end
else
    RA = 0;
end
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

if n ~= 0
    if e > eps 
        w = acos(dot(N,E)/n/e);
        if E(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end

if e > eps
    TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0 
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end


a = h^2/mu/(1 - e^2);

coe = [h e RA incl w TA a E(1) E(2) E(3)];

end

function dTA = cal(n, dt, TA_init, e)
%{
  This function convert the delta time to the delta angle between two
  position along the orbit.

  mu  - gravitational parameter (km^3/s^2)
  dt  - change in time (s)
  n   - sqrt(mu/a^3) (1/s)
  TA_init - TA at r0 position (degrees)
  e   - eccentricity
%}
   
    E = 2*atan(sqrt((1 - e)/(1 + e))*tan(TA_init/2));

    M_init = E - e*sin(E);

    M = M_init + n*dt;
    
    E = kepler_E(e, M);

    TA_final = 2*atan(sqrt((1 + e)/(1 - e))*tan(E/2));

    dTA = TA_final - TA_init;
end

function E = kepler_E(e, M)
%{
  This function uses Newton's Method to solve Kepler's equation 
  E - e*sin(E) = M for the eccentric anomaly, given the eccentricity and
  the mean anomaly.

  E  - eccentric anomaly (radians)
  e  - eccentricity, passed from the calling program
  M  - mean anomaly (radians). passed fro the calling program
  pi - 3.1415926...
%}

error = 1.e-10;

if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

end

function [r, v] = rv_from_r0v0_ta(r0, v0, dt, mu)
%{
  This function computes the state vector (r, v) from the initial
  state vector (r0, v0) and the change in true anomaly.

  mu  - gravitational parameter (km^3/s^2)
  dt  - change in true anomaly (degrees)
  r0  - position vector at time t0 (km)
  v0  - velocity vector at time t0 (km/s)
  r   - final position vector (km)
  v   - final velocity vector (km/s)
%}

%global mu

[f, g] = f_and_g_ta(r0, v0, dt, mu);
[fdot, gdot] = fDot_and_gDot_ta(r0, v0, dt, mu);
r = f*r0 + g*v0;
v = fdot*r0 + gdot*v0;

end

function [f, g] = f_and_g_ta(r0, v0, dt, mu)
%{
  This funciton calculates the Lagrange f and g coefficients form the
  change in true anomaly since time t0.

  mu  - gravitational parameter (km^3/s^2)
  dt  - change in true anomaly (degrees)
  r0  - position vector at time t0 (km)
  v0  - velocity vector at time t0 (km/s)
  h   - angular momentum (km^2/s)
  vr0 - radial component of v0 (km/s)
  r   - radial position after the change in anomaly
  f   - the Lagrange f coefficient (dimensionless)
  g   - the Lagrange g coefficient (s)
%}

h   = norm(cross(r0, v0));
vr0 = dot(v0, r0)/norm(r0);
r0  = norm(r0);
s   = sin(dt);
c   = cos(dt);

r   = h^2/mu/(1 + (h^2/mu/r0 - 1)*c - h*vr0*s/mu);

f   = 1 - mu*r*(1 - c)/h^2;

g   = r*r0*s/h;
end

function [fdot, gdot] = fDot_and_gDot_ta(r0, v0, dt, mu)
%{
  This funciton calculates the time derivatives of the Lagrange f and g 
  coefficients form the change in true anomaly since time t0.

  mu    - gravitational parameter (km^3/s^2)
  dt    - change in true anomaly (degrees)
  r0    - position vector at time t0 (km)
  v0    - velocity vector at time t0 (km/s)
  h     - angular momentum (km^2/s)
  vr0   - radial component of v0 (km/s)
  r     - radial position after the change in anomaly
  fdot  - the time derivatives of the Lagrange f coefficient(1/s)
  gdot  - the time derivatives of the Lagrange g coefficient(dimensionless)
%}

h   = norm(cross(r0, v0));

vr0 = dot(v0, r0)/norm(r0);

r0  = norm(r0);

s   = sin(dt);
c   = cos(dt);

fdot = mu/h*(vr0/h*(1 - c) - s/r0);
gdot = 1 - mu*r0/h^2*(1 - c);

end

function analysis_a_e(statevector, mu)
%{
  This function calculate the changes of a and e for the osculating orbit.
%}
    
    OE = zeros(length(statevector)-1, 10);
    r = statevector(:, 1:3);
    v = statevector(:, 4:6);

    for i = 1:length(statevector)-1
   
        OE(i,:) =  coe_from_sv(r(i, :), v(i, :), mu);
       
    end
    
    e = OE(:, 2);
    a = OE(:, 7);
    plotQ3(e, a);
end

function draw_osculating_orbit(statevector, mu)
%{
  This function is used to calculate the osculating orbits
%}
    
    figure
    OE = zeros(length(statevector)-1, 10);
    r = statevector(:, 1:3);
    v = statevector(:, 4:6);
    
    for i = 1:length(statevector)-1
   
        OE(i,:) =  coe_from_sv(r(i, :), v(i, :), mu);
        
        TA = OE(i,6);

        e = OE(i,2);

        ellipse(r(i, :), v(i, :), TA, e, mu);
        legends{i} = sprintf('Time %d', i); %#ok
        pause(0.2);
        
    end
    
    ax = gca;
    ax.FontWeight = "bold";
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.FontSize = 15;

    xlabel('x - axis  in km','FontSize',13);
    ylabel('y - axis  in km','FontSize',13);
    
    xlim([-3.5e8, 4.e8]);
    ylim([-3.5e8, 3.5e8]);
    title('Osculating Orbit VS perturbed trajectory','FontSize',16, ...
        'FontWeight','bold');
    % legend(legends);
    grid on;
    axis equal;
    fprintf('\nPlot finished!');
end

function ellipse(r, v, TA, e, mu)
%{
  This function is used for draw the ellipitical orbit.
%}
    
    [r0, ~] = rv_from_r0v0_ta(r, v, -TA, mu);
    
    [r1, ~] = rv_from_r0v0_ta(r, v, pi - TA, mu);
    
    draw_ellipse(r0(1), r0(2), r1(1), r1(2), e)

end

function plotQ1(un, cowell)
%{
  This function is used for draw the result of Cowell's method.

  un     - state vector of the unperturbed trajectory.
  cowell - state vector of the perturbed trajectory.
%}

    figure
    unx = un(:, 1);
    uny = un(:, 2);

    cox = cowell(:, 1);
    coy = cowell(:, 2);

    plot(unx, uny, 'LineStyle', '--', 'LineWidth', 1.5);
    hold on;
    plot(cox, coy, 'LineWidth', 2);
    hold on;   
    plot([0, 0], 'r*', 'MarkerSize', 20, 'LineWidth', 2);

    ax = gca;
    ax.FontWeight = "bold";
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.FontSize = 15;

    xlabel('x - axis  in km','FontSize',13);
    ylabel('y - axis  in km','FontSize',13);
    
    xlim([-3.5e8, 4.e8]);
    ylim([-3.5e8, 3.5e8]);
    title('Cowell method','FontSize',16, 'FontWeight','bold');
    legend('Unperturbed', 'Cowell', 'Location', 'northeast');
    grid on;
    axis equal;
    fprintf('Plot finished!\n');
end

function plotQ2(encke, cowell)
%{
  This function is used for compare the results from two different method.

  encke - state vector of the Encke's method.
  cowell - state vector of the Cowell's method.
%}
    figure
    enx = encke(:, 1);
    eny = encke(:, 2);

    cox = cowell(:, 1);
    coy = cowell(:, 2);

    plot(enx, eny, 'ok');
    hold on;
    plot(cox, coy, 'LineWidth', 1.5);
    hold on;   
    plot([0, 0], 'r*', 'MarkerSize', 20, 'LineWidth', 2);

    ax = gca;
    ax.FontWeight = "bold";
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.FontSize = 15;

    xlabel('x - axis  in km','FontSize',13);
    ylabel('y - axis  in km','FontSize',13);
    
    xlim([-3.5e8, 4.e8]);
    ylim([-3.5e8, 3.5e8]);
    title('Encke method compared with cowell','FontSize',16, ...
        'FontWeight','bold');
    legend('Encke', 'Cowell', 'Location', 'northeast');
    grid on;
    axis equal;
    fprintf('\nPlot finished!');
end

function plotQ3(e, a)

    figure
    ax = gca;
    ax.FontWeight = "bold";
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.FontSize = 15;

    plot(e,'LineWidth',1.5);
    xlabel('Time  in year','FontSize',13);
    ylabel('Eccentricity','FontSize',13);
    xlim([-10, 510]);
    ylim([-0.002, 0.035]);
    xticks(0:125:500)
    xticklabels({'0y','5y','10y','15y','20y'})
    title('Variation of encentricity','FontSize',16, ...
        'FontWeight','bold');
    % legend('Encke', 'Cowell', 'Location', 'northeast');
    grid on;

    figure
    ax = gca;
    ax.FontWeight = "bold";
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.FontSize = 15;

    plot(a,'LineWidth', 1.5);
    xlabel('Time  in year','FontSize',13);
    ylabel('semi-major axis in km','FontSize',13);
    xlim([-10, 510]);
    ylim("auto");
    xticks(0:125:500)
    xticklabels({'0y','5y','10y','15y','20y'})
    title('Variation of semi-major axis','FontSize',16, ...
        'FontWeight','bold');
    % legend('Encke', 'Cowell', 'Location', 'northeast');
    grid on;
      
    fprintf('\nPlot finished!');
end

function draw_ellipse(x1, y1, x2, y2, e)
%{
  This function is used for drawing a ellipse.
%}

% x1, y1, x2, y2 are the coordinates of the two vertices of the ellipse
% e is the eccentricity of the ellipse

% calculate the major and minor axes
a = sqrt((x2-x1)^2 + (y2-y1)^2) / 2;
b = a * sqrt(1 - e^2);

% calculate the center of the ellipse
xc = (x1 + x2) / 2;
yc = (y1 + y2) / 2;

% calculate the angle between the major axis and the x-axis
theta = atan2(y2-y1, x2-x1);

% generate the parameter values
t = linspace(0, 2*pi);

% calculate the x and y coordinates of the ellipse
x = xc + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
y = yc + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);

distance1 = norm([x1, y1]);
distance2 = norm([x2, y2]);
% plot the ellipse
plot(x, y, 'LineWidth', 2, 'LineStyle',':')
hold on;
plot([0, 0], 'r*', 'MarkerSize', 20, 'LineWidth', 2);
hold on;
if (distance1 < distance2)
    plot(x1, y1, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
else
    plot(x2, y2, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
axis equal
end
