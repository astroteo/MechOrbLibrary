%% PERTURBATION PLOTS COMPARISON
clear all
close all

%% DATA
mu_Eu=3.1888e3;    %[km^3/s^2]
R_Eu=1560.8;       %[km]
h_p=120;           %[km]
rp_f=h_p+R_Eu;     %[km]
h_a=200;           %[km]
ra_f=h_a+R_Eu;     %[km]
a_f=(ra_f+rp_f)/2; %[km]
e_f=1-(rp_f/a_f);  %[-]
p_f=a_f*(1-e_f^2); %[km]
i_f=78;     %[deg]
OM_f=250;   %[deg]
om_f=80;    %[deg]
theta_f=0;  %[deg]
T_f=sqrt((a_f^3)/mu_Eu)*2*pi; %[s]

KP_in=[a_f;e_f;i_f;OM_f;om_f;theta_f];

options=odeset('RelTol',1e-12,'AbsTol',1e-24);

%% SHORT PERIOD
dt_short=(0:T_f);
[t_pert_NUM,KP_pert]= ode113('kp_mass',dt_short,KP_in,options);
[t_pert_SRP,solar_wind_pert]=ode113('kp_solar_wind',dt_short,KP_in,options);
[t_pert_3RDB,third_body_pert]=ode113('kp_3rdb',dt_short,KP_in,options);
[t_pert,total_pert]= ode113('kp_total',dt_short,KP_in,options);

i=0;
r=[2000 2000 2000];
while i<length(dt_short) && norm(r)>R_Eu
    i=i+1;
    [r]=OrbPar2RV(KP_pert(i,1),KP_pert(i,2),KP_pert(i,3),KP_pert(i,4),KP_pert(i,5),KP_pert(i,6),mu_Eu);
    IMPACT_TIME=dt_short(i);
end

if i~=length(dt_short)
    disp('SATELLITE IMPACTS ON SOIL DUE TO PERTURBATIONS AFTER')
    [YEARS,MONTHS,DAYS,HOURS,MINUTES,SECONDS]=time2impact(IMPACT_TIME)
end

% NON UNIFORM MASS PLOTS
figure(1)
plot(t_pert,KP_pert(:,1),'--c')
title('PERTURBATION PLOTS COMPARISON - "a" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('a [Km]')
grid on
hold on

figure(2)
plot(t_pert,KP_pert(:,2),'--c')
title('PERTURBATION PLOTS COMPARISON - "e" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('e [-]')
grid on
hold on

figure(3)
plot(t_pert,KP_pert(:,3),'--c')
title('PERTURBATION PLOTS COMPARISON - "i" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('i [deg]')
grid on
hold on

figure(4)
plot(t_pert,KP_pert(:,4),'--c')
title('PERTURBATION PLOTS COMPARISON - "OM" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('OM [deg]')
grid on
hold on

figure(5)
plot(t_pert,KP_pert(:,5),'--c')
title('PERTURBATION PLOTS COMPARISON - "om" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('om [deg]')
grid on
hold on

% SOLAR RADIATION PRESSURE PLOTS
figure(1)
plot(t_pert,solar_wind_pert(:,1),'--g')

figure(2)
plot(t_pert,solar_wind_pert(:,2),'--g')

figure(3)
plot(t_pert,solar_wind_pert(:,3),'--g')

figure(4)
plot(t_pert,solar_wind_pert(:,4),'--g')

figure(5)
plot(t_pert,solar_wind_pert(:,5),'--g')

% THIRD BODY PLOTS
figure(1)
plot(t_pert,third_body_pert(:,1),'--b')

figure(2)
plot(t_pert,third_body_pert(:,2),'--b')

figure(3)
plot(t_pert,third_body_pert(:,3),'--b')

figure(4)
plot(t_pert,third_body_pert(:,4),'--b')

figure(5)
plot(t_pert,third_body_pert(:,5),'--b')

% TOTAL PLOTS
figure(1)
plot(t_pert,total_pert(:,1),'r')
legend('Non uniform mass','Solar radiation pressure','Third body presence','Global')

figure(2)
plot(t_pert,total_pert(:,2),'r')
legend('Non uniform mass','Solar radiation pressure','Third body presence','Global')

figure(3)
plot(t_pert,total_pert(:,3),'r')
legend('Non uniform mass','Solar radiation pressure','Third body presence','Global')

figure(4)
plot(t_pert,total_pert(:,4),'r')
legend('Non uniform mass','Solar radiation pressure','Third body presence','Global')

figure(5)
plot(t_pert,total_pert(:,5),'r')
legend('Non uniform mass','Solar radiation pressure','Third body presence','Global')

