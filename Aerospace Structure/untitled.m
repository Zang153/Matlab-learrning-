clc;
close all;
clear all;

load('wing loading 3.75G.mat');
load('wing loading -1.5G.mat');

stn = wingloading375G(:,1);

SFP = wingloading375G(:,2);
SFN = -1*WingLS(:,2);
BMP = wingloading375G(:,3);
BMN = -1*WingLS(:,3);
TP = wingloading375G(:,4);
TN = -1*WingLS(:,4);

figure
figure_FontSize=12;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(stn, SFP, stn, SFN);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Shear force (N)', 'FontSize', figure_FontSize,'FontWeight','bold');
legend('3.75G', '-1.5G');
grid on
%------------------------------------------------------------------------%
figure
plot(stn, BMP, stn, BMN);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Bending moment (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');
legend('3.75G', '-1.5G');
grid on
%------------------------------------------------------------------------%
figure
plot(stn, TP, stn, TN);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Torque (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');
legend('3.75G', '-1.5G');
grid on
%------------------------------------------------------------------------%

%% 

close all;

figure_FontSize=10;

stn = wingaeroload(:,1);

inertial_force = -1*winginertialload(:,5);

aero_force = wingaeroload(:,3);

inertial_SF = -1*winginertialload(:,6);

inertial_BM = -1*winginertialload(:,8);

aero_SF = wingaeroload(:,4);

aero_BM = wingaeroload(:,6);

all_SF = aeroplusinertial(:,1);

all_BM = aeroplusinertial(:,2);

ET_SF = bendingduetoengine(:,3);

ET_BM = bendingduetoengine(:,5);

%--------------------------------------------------------------------------%
% figure
% title('1G flight');
% plot(stn, inertial_force, stn, aero_force);
% xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
% ylabel('load (N)', 'FontSize', figure_FontSize,'FontWeight','bold');
% 
% legend('Weight','Lift')
% grid on;
% %--------------------------------------------------------------------------%
% figure
% title('Shear force');
% plot(stn, inertial_SF, stn, aero_SF, stn, all_SF);
% xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
% ylabel('Shear force (N)', 'FontSize', figure_FontSize,'FontWeight','bold');
% 
% legend('inertial SF', 'aero SF','all SF');
% grid on;
% %--------------------------------------------------------------------------%
% figure
% title('Shear force');
% plot(stn, inertial_BM, stn, aero_BM, stn, all_BM);
% xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
% ylabel('Bending Moment (N/m)', 'FontSize', figure_FontSize,'FontWeight','bold');
% 
% legend('inertial BM', 'aero Bm', 'all BM');
% grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, ET_BM);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Bending Moment (N/m)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('Thrust BM');
grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, ET_SF);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Shear force (N)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('Thrust SF');
grid on;
%%
close all;

load('HT diagram data.mat');

stn = tailcase1(:,1);

c1_SF = -1*tailcase1(:,2);

c1_BM = -1*tailcase1(:,3);

c1_Tq = -1*tailcase1(:,4);

c2_SF = tailcase2(:,2);

c2_BM = tailcase2(:,3);

c2_Tq = tailcase2(:,4);

figure_FontSize = 12;
figure
plot(stn, c1_SF, stn, c2_SF);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Shear force (N)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case1','case2')
grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, c1_Tq, stn, c2_Tq);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Torque (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case1','case2')
grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, c1_BM , stn, c2_BM);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Bending Moment (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case1','case2')
grid on;



%%

close all;

load('HT diagram data.mat');

stn = VTcase2(:,1);

c2_SF = VTcase2(:,2);

c2_BM = VTcase2(:,3);

c2_Tq = VTcase2(:,4);

figure_FontSize = 12;
figure
plot(stn, c2_SF);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Shear force (N)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case2')
grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, c2_Tq);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Torque (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case2')
grid on;
%--------------------------------------------------------------------------%
figure
plot(stn, c2_BM);
xlabel('Station points (per station)', 'FontSize', figure_FontSize,'FontWeight','bold','Units','centimeters');
ylabel('Bending Moment (N m)', 'FontSize', figure_FontSize,'FontWeight','bold');

legend('case2')
grid on;




















































