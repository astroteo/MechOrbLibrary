function solar_wind_pert= kp_solar_wind(t,KP)
%% DATA
% Initial orbital parameters SATELLITE
a=KP(1);            %[Km]
e=KP(2);            %[-]
i=KP(3);            %[deg]
Omega=KP(4);        %[deg]
omega=KP(5);        %[deg]
theta0=KP(6);       %[deg]

% Gravity constants
%mu_J=126711995;   %[km^3/s^2]
mu_Eu=3.1888e3;   %[km^3/s^2]

% Sun-Satellite data
A_sat=3;          %[m^2]
m_sat=780;        %[kg]
eps=0.34;         %[-]
solar_flux=50.57; %[W/m^2]
c=299792458;      %[m/s]

% Initial time data
%year0=2018;
load('OPTIMUM_data.mat')
t_tot=t0+tof_1+tof_2; %[s]

% Data EUROPA
load('HORIZON_data.mat','R_Eu_data','R_Ju_data')

% Plantes radius
r_j=69911;   %[km]
r_eu=1560.8; %[km]


[r_imp]=OrbPar2RV(a,e,i,Omega,omega,theta0,mu_Eu);

if norm(r_imp)>r_eu
    %% JUPITER-EUROPA-SATELLITE ABSOLUT POSITIONS [I,J,K]
    [R_J]=interp_data(R_Ju_data,2,t+t_tot);                     %[km]
    % DAY=365*(year0-2000)+t_tot/(3600*24)+t/(3600*24);
    % [E,mu_s]=uplanet_mick(DAY,5);
    % a_J=E(1);     %[km]
    % e_J=E(2);     %[-]
    % i_J=E(3);     %[deg]
    % OM_J=E(4);    %[deg]
    % om_J=E(5);    %[deg]
    % theta_J=E(6); %[deg]
    % [R_J]=OrbPar2RV(a_J,e_J,i_J,OM_J,om_J,theta_J,mu_s);         %[km]
    [R_Eu]=interp_data(R_Eu_data,2,t+t_tot);                     %[km]
    % [theta_Eu]=time2theta(a_Eu,e_Eu,theta0_Eu,t+t_tot,mu_J);     %[deg]
    % [R_Eu]=OrbPar2RV(a_Eu,e_Eu,i_Eu,OM_Eu,om_Eu,theta_Eu,mu_J);  %[km]
    [R_sat]=OrbPar2RV(a,e,i,Omega,omega,theta0,mu_Eu);           %[km]
    vers_s=(R_J+R_Eu+R_sat)/norm(R_J+R_Eu+R_sat);


    %% ECLIPSE FACTOR
    R_Eu_sun=-R_J-R_Eu;        %[km]
    r_Eu_shade=r_j;           %[km]
    r_sat_shade=r_eu;         %[km]
    beta_Eu_shade=acos(dot(R_Eu,-R_J)/(norm(R_Eu)*norm(R_J)));               %[rad]
    beta_sat_shade=acos(dot(R_sat,R_Eu_sun)/(norm(R_sat)*norm(R_Eu_sun)));  %[rad]
    if beta_Eu_shade>pi/2 && norm(R_Eu)*sin(beta_Eu_shade)<r_Eu_shade
        v=0;
    else
        if beta_sat_shade>pi/2 && norm(R_sat)*sin(beta_sat_shade)<r_sat_shade
            v=0;
        else
            v=1;
        end
    end


    %% SOLAR WIND ACCELERATION
    F_sw=((v*solar_flux*A_sat*(1+eps))/(c*m_sat*1e3))*vers_s;  %[km/s^2]

    % Acceleration in [r,theta,h]
    i=i*(pi/180);
    Omega=Omega*(pi/180);
    omega=omega*(pi/180);
    theta=theta0*(pi/180);

    ROT=[cos(omega)*cos(Omega)-sin(omega)*cos(i)*sin(Omega) cos(omega)*sin(Omega)+sin(omega)*cos(i)*cos(Omega) sin(omega)*sin(i)
        -sin(omega)*cos(Omega)-cos(omega)*cos(i)*sin(Omega) -sin(omega)*sin(Omega)+cos(omega)*cos(i)*cos(Omega) cos(omega)*sin(i)
        sin(i)*sin(Omega) -sin(i)*cos(Omega) cos(i)];

    rot=[cos(theta),sin(theta),0
        -sin(theta),cos(theta),0
        0,0,1];

    F_sw_rth=rot*ROT*F_sw; %[km/s^2]

    Fr=F_sw_rth(1);        %[Km/s^2]
    Ftheta=F_sw_rth(2);    %[Km/s^2]
    Fh=F_sw_rth(3);        %[Km/s^2]

    r = a*(1-e^2)/(1+e*cos(theta));     %[Km]
    p=a*(1-e^2);                        %[km]
    h=sqrt(mu_Eu*p);                    %[km]
    E= 2*atan(tan(theta/2)*sqrt((1-e)/(1+e)));         %[rad]


    %% VARIATION OF KEPLERIAN PARAMETERS IN [r,theta,h]
    a_dot = 2*sqrt(a^3/(mu_Eu*(1-e^2)))*(e*sin(theta)*Fr+(1+e*cos(theta))*Ftheta);
    e_dot = sqrt(a*(1-e^2)/mu_Eu)*(sin(theta)*Fr+(cos(theta)+cos(E))*Ftheta);
    i_dot = sqrt(a*(1-e^2)/mu_Eu)*(cos(theta+omega)/(1+e*cos(theta)))*Fh;
    Omega_dot = sqrt(a*(1-e^2)/mu_Eu)*(sin(theta+omega)/((1+e*cos(theta))*sin(i)))*Fh;
    omega_dot = sqrt(a*(1-e^2)/mu_Eu)*(-(cos(theta)/e)*Fr+sin(theta)*((2+e*cos(theta))/(e*(1+e*cos(theta))))*Ftheta-((sin(theta+omega)*cot(i))/(1+e*cos(theta)))*Fh);
    theta_dot=h/(r^2)+(p*cos(theta)*Fr-(p+r)*sin(theta)*Ftheta)/(e*h);
else
    a_dot=0;
    e_dot=0;
    i_dot=0;
    Omega_dot=0;
    omega_dot=0;
    theta_dot=0;
end

solar_wind_pert= [a_dot;e_dot;i_dot*(180/pi);Omega_dot*(180/pi);omega_dot*(180/pi);theta_dot*(180/pi)];
