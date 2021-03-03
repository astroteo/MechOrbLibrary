function KP_d= kp_mass(t,KP)
%% DATA
a=KP(1);            %[Km]
e=KP(2);            %[-]
i=KP(3)*(pi/180);            %[rad]
omega=KP(5)*(pi/180);        %[rad]
theta=KP(6)*(pi/180);        %[rad]

mu_Eu=3.1888e+03;
J2=435.5e-6;
r_eu=1560.8;

r = a*(1-e^2)/(1+e*cos(theta));     %[Km]
p=a*(1-e^2);                        %[km]
h=sqrt(mu_Eu*p);                    %[km]
E= 2*atan(tan(theta/2)*sqrt((1-e)/(1+e)));         %[rad]

[r_imp]=OrbPar2RV(a,e,i,Omega,omega,theta0,mu_Eu);

if norm(r_imp)>r_eu
    %% ACCELERATION IN [r,tan,h]
    Fr = (-(mu_Eu*J2*(r_eu^2)*3)/(2*r^4))*(1-3*(sin(i))^2*(sin(theta+omega))^2);            %[Km/s^2]
    Ftheta = (-(mu_Eu*J2*(r_eu^2)*3)/(r^4))*(sin(i))^2*sin(theta+omega)*cos(theta+omega);     %[Km/s^2]
    Fh =(-(mu_Eu*J2*(r_eu^2)*3)/(r^4))*cos(i)*sin(i)*sin(theta+omega);                      %[Km/s^2]

    %% VARIATIONS OF KEPLERIAN PARAMETERS IN [r,tan,h]
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

KP_d= [a_dot; e_dot; i_dot*(180/pi); Omega_dot*(180/pi); omega_dot*(180/pi); theta_dot*(180/pi)];