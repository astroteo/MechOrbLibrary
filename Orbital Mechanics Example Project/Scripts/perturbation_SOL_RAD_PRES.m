%% PERTURBATION DUE TO SOLAR RADIATION PRESSURE
clear all
close all

%% DATA
% Final orbit data
mu_Eu=3.1888e3;    %[km^3/s^2]
r_e=1560.8;       %[km]
h_p=120;           %[km]
rp_f=h_p+r_e;     %[km]
h_a=200;           %[km]
ra_f=h_a+r_e;     %[km]
a_f=(ra_f+rp_f)/2; %[km]
e_f=1-(rp_f/a_f);  %[-]
p_f=a_f*(1-e_f^2); %[km]
i_f=78;     %[deg]
OM_f=250;   %[deg]
om_f=80;    %[deg]
theta_f=0; %[deg]
T_f=sqrt((a_f^3)/mu_Eu)*2*pi; %[s]

KP_in=[a_f;e_f;i_f;OM_f;om_f;theta_f];

options=odeset('RelTol',1e-12,'AbsTol',1e-24);


%% SHORT PERIOD
dt_short=(0:T_f);
%[t_pert,solar_wind_pert]=ode113('pert_solar_wind',(0:1:T_f),KP_in);
[t_pert,solar_wind_pert]=ode113('kp_solar_wind',dt_short,KP_in,options);

i=0;
r=[2000 2000 2000];
while i<length(dt_short) && norm(r)>r_e
    i=i+1;
    [r]=OrbPar2RV(solar_wind_pert(i,1),solar_wind_pert(i,2),solar_wind_pert(i,3),solar_wind_pert(i,4),solar_wind_pert(i,5),solar_wind_pert(i,6),mu_Eu);
    IMPACT_TIME=dt_short(i);
end

if i~=length(dt_short)
    disp('SATELLITE IMPACTS ON SOIL DUE TO PERTURBATIONS AFTER')
    [YEARS,MONTHS,DAYS,HOURS,MINUTES,SECONDS]=time2impact(IMPACT_TIME)
end

% Short period perturbation plots
figure(1)
plot(t_pert,solar_wind_pert(:,1))
title('SOLAR RADIATION PRESSURE - "a" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('a [Km]')
grid on

figure(2)
plot(t_pert,solar_wind_pert(:,2))
title('SOLAR RADIATION PRESSURE - "e" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('e [-]')
grid on

figure(3)
plot(t_pert,solar_wind_pert(:,3))
title('SOLAR RADIATION PRESSURE - "i" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('i [deg]')
grid on

figure(4)
plot(t_pert,solar_wind_pert(:,4))
title('SOLAR RADIATION PRESSURE - "OM" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('OM [deg]')
grid on

figure(5)
plot(t_pert,solar_wind_pert(:,5))
title('SOLAR RADIATION PRESSURE - "om" Perturbations over one period','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('om [deg]')
grid on



%% LONG PERIOD
AU=149597870.700 ;
mu_sun=1.372e+11;
T_year=2*pi*sqrt(AU^3/mu_sun);

% lp means long period
step=1e3;
dt_long=(0:step:T_year);
%[t_pert_lp,solar_wind_pert_lp]=ode113('pert_solar_wind',(0:step:T_year),KP_in);
[t_pert_lp,solar_wind_pert_lp]=ode113('kp_solar_wind',dt_long,KP_in,options);

i=0;
r=[2000 2000 2000];
while i<length(dt_long) && norm(r)>r_e
    i=i+1;
    [r]=OrbPar2RV(solar_wind_pert_lp(i,1),solar_wind_pert_lp(i,2),solar_wind_pert_lp(i,3),solar_wind_pert_lp(i,4),solar_wind_pert_lp(i,5),solar_wind_pert_lp(i,6),mu_Eu);
    IMPACT_TIME_lp=dt_long(i);
    t=i;
end

if i~=length(dt_long)
    disp('SATELLITE IMPACTS ON SOIL DUE TO PERTURBATIONS AFTER')
    [YEARS,MONTHS,DAYS,HOURS,MINUTES,SECONDS]=time2impact(IMPACT_TIME_lp)
    t_pert_lp=t_pert_lp(1:t);
    solar_wind_pert_lp=solar_wind_pert_lp(1:t,:);
end

% Long period perturbation plots
figure(6)
plot(t_pert_lp,solar_wind_pert_lp(:,1))
title('SOLAR RADIATION PRESSURE - "a" Perturbations over a solar year','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('a [Km]')
grid on
hold on

figure(7)
plot(t_pert_lp,solar_wind_pert_lp(:,2))
title('SOLAR RADIATION PRESSURE - "e" Perturbations over a solar year','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('e [-]')
grid on
hold on

figure(8)
plot(t_pert_lp,solar_wind_pert_lp(:,3))
title('SOLAR RADIATION PRESSURE - "i" Perturbations over a solar year','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('i [deg]')
grid on
hold on

figure(9)
plot(t_pert_lp,solar_wind_pert_lp(:,4))
title('SOLAR RADIATION PRESSURE - "OM" Perturbations over a solar year','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('OM [deg]')
grid on
hold on

figure(10)
plot(t_pert_lp,solar_wind_pert_lp(:,5))
title('SOLAR RADIATION PRESSURE - "om" Perturbations over a solar year','FontSize',13,'FontWeight','bold')
xlabel('time [s]')
ylabel('om [deg]')
grid on
hold on

sprintf('\t\t\t ORBITAL PARAMETERS AFTER 1 YEAR OR IMPACT \n\n \t\t\t\t    a \t\t   e \t\t   i \t\t   OM \t\t   om \n Initial operating orbit: \t %4.2f [km] \t %1.4f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg] \n Final operating orbit: \t %4.2f [km] \t %1.4f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg]',a_f,e_f,i_f,OM_f,om_f,solar_wind_pert_lp(i,1),solar_wind_pert_lp(i,2),solar_wind_pert_lp(i,3),solar_wind_pert_lp(i,4),solar_wind_pert_lp(i,5))

%% FILTERING LONG PERIOD
%Mean values evaluated over NN periods
NN=50;
[N,p]=size(solar_wind_pert_lp);
dif=abs(t_pert_lp-NN*T_f);
for i=1:length(dif)
    if dif(i)==min(dif)
        n=i;
    end
end
M=floor(N/n);
solar_wind_pert_fp=zeros(M,p);
time_pert_fp=zeros(M,1);
for i=1:6
    for j=1:M
        solar_wind_pert_fp(j,i)=sum(solar_wind_pert_lp(((j-1)*n+1:j*n),i))/n;
        time_pert_fp(j)=t_pert_lp(n*(j-1)+round(n/2));
    end
end
time_pert_fp=[0;time_pert_fp;t_pert_lp(end)];
solar_wind_pert_fp=[sum(solar_wind_pert_lp(1:round(n/3),:))/round(n/3);solar_wind_pert_fp;sum(solar_wind_pert_lp(end-round(n/3):end,:))/(round(n/3)+1)];

figure(6)
plot(time_pert_fp,solar_wind_pert_fp(:,1),'r')
legend('a_p_e_r_t','a_F_I_L_T_ _p_e_r_t')

figure(7)
plot(time_pert_fp,solar_wind_pert_fp(:,2),'r')
legend('e_p_e_r_t','e_F_I_L_T_ _p_e_r_t')

figure(8)
plot(time_pert_fp,solar_wind_pert_fp(:,3),'r')
legend('i_p_e_r_t','i_F_I_L_T_ _p_e_r_t')

figure(9)
plot(time_pert_fp,solar_wind_pert_fp(:,4),'r')
legend('OM_p_e_r_t','OM_F_I_L_T_ _p_e_r_t')

figure(10)
plot(time_pert_fp,solar_wind_pert_fp(:,5),'r')
legend('om_p_e_r_t','om_F_I_L_T_ _p_e_r_t')

if t~=length(dt_long)
    figure(6)
    plot(IMPACT_TIME_lp,solar_wind_pert_lp(t,1),'*g','LineWidth',2)
    legend('a_p_e_r_t','a_F_I_L_T_ _p_e_r_t','Impact on soil')
    figure(7)
    plot(IMPACT_TIME_lp,solar_wind_pert_lp(t,2),'*g','LineWidth',2)
    legend('e_p_e_r_t','e_F_I_L_T_ _p_e_r_t','Impact on soil')
    figure(8)
    plot(IMPACT_TIME_lp,solar_wind_pert_lp(t,3),'*g','LineWidth',2)
    legend('i_p_e_r_t','i_F_I_L_T_ _p_e_r_t','Impact on soil')
    figure(9)
    plot(IMPACT_TIME_lp,solar_wind_pert_lp(t,4),'*g','LineWidth',2)
    legend('OM_p_e_r_t','OM_F_I_L_T_ _p_e_r_t','Impact on soil')
    figure(10)
    plot(IMPACT_TIME_lp,solar_wind_pert_lp(t,5),'*g','LineWidth',2)
    legend('om_p_e_r_t','om_F_I_L_T_ _p_e_r_t','Impact on soil')
end