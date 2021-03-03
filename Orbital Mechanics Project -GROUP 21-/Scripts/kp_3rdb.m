function third_body_pert= kp_3rdb(t,KP)
%% DATA
% Initial orbital parameters SATELLITE
a=KP(1);            %[Km]
e=KP(2);            %[-]
i=KP(3);            %[deg]
Omega=KP(4);        %[deg]
omega=KP(5);        %[deg]
theta0=KP(6);       %[deg]

% Gravity constants
mu_J=126711995;   %[km^3/s^2]
mu_Eu=3.1888e3;   %[km^3/s^2]

% Initial time data
load('OPTIMUM_data.mat')
t_tot=t0+tof_1+tof_2; %[s]

% Data EUROPA
load('HORIZON_data.mat','R_Eu_data')

r_e=1560.8;       %[km]
[r_imp]=OrbPar2RV(a,e,i,Omega,omega,theta0,mu_Eu);

if norm(r_imp)>r_e
    %% EUROPA-SATELLITE ABSOLUT POSITIONS [I,J,K]
    [R_Eu]=interp_data(R_Eu_data,2,t+t_tot);                     %[km]
    [R_sat]=OrbPar2RV(a,e,i,Omega,omega,theta0,mu_Eu);           %[km]


    %% THIRD BODY ACCELERATION
    R_sat_J=R_Eu+R_sat;
    F_3rdb=mu_J*((R_sat_J-R_sat)/(norm(R_Eu)^3)-R_sat_J/(norm(R_sat_J)^3));  %[km/s^2]

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

    F_3rdb_rth=rot*ROT*F_3rdb; %[km/s^2]

    Fr=F_3rdb_rth(1);        %[Km/s^2]
    Ftheta=F_3rdb_rth(2);    %[Km/s^2]
    Fh=F_3rdb_rth(3);        %[Km/s^2]

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

third_body_pert= [a_dot;e_dot;i_dot*(180/pi);Omega_dot*(180/pi);omega_dot*(180/pi);theta_dot*(180/pi)];
