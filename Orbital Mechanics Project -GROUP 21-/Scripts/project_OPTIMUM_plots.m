%% OPTIMUM PLOTS
clear all
close all
%% DATA
t0=input('\nDEFINE THE DEPARTURE DATE \n t0[s]:  ');
tof_1=input('\nDEFINE THE 1^(st) TIME OF FLIGHT \n tof_1[s]:  ');
tof_2=input('\nDEFINE THE 2^(nd) TIME OF FLIGHT \n tof_2[s]:  ');

% GANYMEDE-EUROPA-IO data
load('HORIZON_data.mat')

% Reference axis
I=[1;0;0];
J=[0;1;0];
K=[0;0;1];

% Planetary gravity constants
mu_J=126711995;
mu_Io=((8.9319E+22)/(5.9736E+24))*398600;
mu_Eu=3.1888e3;

% Plantes radius
r_j=69911;    %[km]
r_g=5262/2;   %[km]
r_i=1821;     %[km]
r_e=3121/2;   %[km]

[x,y,z]=sphere;

% Orbits plot
r_Ga_t0=R_Ga_data(t0/7.2e3+1,:);
r_Eu_t0=R_Eu_data(t0/7.2e3+1,:);
r_Io_t0=R_Io_data(t0/7.2e3+1,:);
v_Ga_t0=V_Ga_data(t0/7.2e3+1,:);
v_Eu_t0=V_Eu_data(t0/7.2e3+1,:);
v_Io_t0=V_Io_data(t0/7.2e3+1,:);
[a_Ga,e_Ga,i_Ga,OM_Ga,om_Ga,theta0_Ga]=RV2OrbPar(r_Ga_t0,v_Ga_t0,mu_J);
[a_Io,e_Io,i_Io,OM_Io,om_Io]=RV2OrbPar(r_Io_t0,v_Io_t0,mu_J);
[a_Eu,e_Eu,i_Eu,OM_Eu,om_Eu]=RV2OrbPar(r_Eu_t0,v_Eu_t0,mu_J); 
r_Ga=ones(3,360);
v_Ga=ones(3,360);
r_Io=ones(3,360);
v_Io=ones(3,360);
r_Eu=ones(3,360);
v_Eu=ones(3,360);
for k=1:360
[r_Ga(:,k),v_Ga(:,k)]=OrbPar2RV(a_Ga,e_Ga,i_Ga,OM_Ga,om_Ga,k,mu_J);
[r_Io(:,k),v_Io(:,k)]=OrbPar2RV(a_Io,e_Io,i_Io,OM_Io,om_Io,k,mu_J);
[r_Eu(:,k),v_Eu(:,k)]=OrbPar2RV(a_Eu,e_Eu,i_Eu,OM_Eu,om_Eu,k,mu_J);
end
figure()
plot3(r_Ga(1,:),r_Ga(2,:),r_Ga(3,:),'--b',r_Io(1,:),r_Io(2,:),r_Io(3,:),'--r',r_Eu(1,:),r_Eu(2,:),r_Eu(3,:),'g',0,0,0,'ok')
legend('Ganymede','Io','Europa','Jupiter')
grid on

figure()
plot3(r_Ga(1,:),r_Ga(2,:),r_Ga(3,:),'--b',r_Io(1,:),r_Io(2,:),r_Io(3,:),'--r',r_Eu(1,:),r_Eu(2,:),r_Eu(3,:),'g')
hold on
surf(x*r_j,y*r_j,z*r_j)
colormap('summer')
axis([-1.2e6 1.2e6 -1.2e6 1.2e6 -1.2e6 1.2e6])
legend('Ganymede','Io','Europa','Jupiter')
grid on



%% LAMBERT TRANSFERS
% First transfer arc.
r1=R_Ga_data(t0/7.2e3+1,:);
v1=V_Ga_data(t0/7.2e3+1,:);
t_1=t0+tof_1;
r2=R_Io_data(t_1/7.2e3+1,:);
v2=V_Io_data(t_1/7.2e3+1,:);
[a_t1,p_L1,e_t1,error1,v1_L1,v2_L1,theta1]=lambert_mick(r1,r2,tof_1,mu_J);
[a_L1,e_L1,i_L1,OM_L1,om_L1,theta1_L1]=RV2OrbPar(r1,v1_L1',mu_J);
[theta2_L1]=time2theta(a_L1,e_L1,theta1_L1,tof_1,mu_J);
if theta1_L1>theta2_L1
    theta2_L1=theta2_L1+360;
end
L1_it=(theta1_L1:1:theta2_L1);
r_L1=ones(3,length(L1_it));
v_L1=ones(3,length(L1_it));
for k=1:length(L1_it)
[r_L1(:,k),v_L1(:,k)]=OrbPar2RV(a_L1,e_L1,i_L1,OM_L1,om_L1,L1_it(k),mu_J);
end

cost_L=norm(v1_L1-v1);

% Second transfer arc
t_2=t0+tof_1+tof_2;
r3=R_Eu_data(t_2/7.2e3+1,:);
v3=V_Eu_data(t_2/7.2e3+1,:);
[a_t2,p_L2,e_t2,error2,v1_L2,v2_L2,theta2]=lambert_mick(r2,r3,tof_2,mu_J);
[a_L2,e_L2,i_L2,OM_L2,om_L2,theta1_L2]=RV2OrbPar(r2,v1_L2',mu_J);
[theta2_L2]=time2theta(a_L2,e_L2,theta1_L2,tof_2,mu_J);
if theta1_L2>theta2_L2
    theta2_L2=theta2_L2+360;
end
L2_it=(theta1_L2:1:theta2_L2);
r_L2=ones(3,length(L2_it));
v_L2=ones(3,length(L2_it));
for k=1:length(L2_it)
[r_L2(:,k),v_L2(:,k)]=OrbPar2RV(a_L2,e_L2,i_L2,OM_L2,om_L2,L2_it(k),mu_J);
end


% Lambert's transfers plot
[a_Io,e_Io,i_Io,OM_Io,om_Io,theta0_Io]=RV2OrbPar(r2,v2,mu_J);
[a_Eu,e_Eu,i_Eu,OM_Eu,om_Eu,theta0_Eu]=RV2OrbPar(r3,v3,mu_J); 
for k=1:360
[r_Io(:,k),v_Io(:,k)]=OrbPar2RV(a_Io,e_Io,i_Io,OM_Io,om_Io,k,mu_J);
[r_Eu(:,k),v_Eu(:,k)]=OrbPar2RV(a_Eu,e_Eu,i_Eu,OM_Eu,om_Eu,k,mu_J);
end

figure(3)
plot3(r_Ga(1,:),r_Ga(2,:),r_Ga(3,:),'--b',r_Io(1,:),r_Io(2,:),r_Io(3,:),'--g',r_Eu(1,:),r_Eu(2,:),r_Eu(3,:),'r',r_L1(1,:),r_L1(2,:),r_L1(3,:),'c',r_L2(1,:),r_L2(2,:),r_L2(3,:),'--c',r1(1),r1(2),r1(3),'*b',r2(1),r2(2),r2(3),'*g',r3(1),r3(2),r3(3),'*r',0,0,0,'ok')
legend('Ganymede','Io','Europa','Lambert_1','Lambert_2','start','GA','target','Jupiter')
grid on
hold on

figure()
plot3(r_Ga(1,:),r_Ga(2,:),r_Ga(3,:),'--b',r_Io(1,:),r_Io(2,:),r_Io(3,:),'--g',r_Eu(1,:),r_Eu(2,:),r_Eu(3,:),'r',r_L1(1,:),r_L1(2,:),r_L1(3,:),'c',r_L2(1,:),r_L2(2,:),r_L2(3,:),'--c')
hold on
surf(x*r_j,y*r_j,z*r_j)
colormap('summer')
surf(x*r_g+r1(1),y*r_g+r1(2),z*r_g+r1(3))
surf(x*r_i+r2(1),y*r_i+r2(2),z*r_i+r2(3))
surf(x*r_e+r3(1),y*r_e+r3(2),z*r_e+r3(3))
axis([-1.2e6 1.2e6 -1.2e6 1.2e6 -1.2e6 1.2e6])
legend('Ganymede','Io','Europa','Lambert_1','Lambert_2','Jupiter','Ganymede','Io','Europa')
grid on



%% POWERED GRAVITY ASSIST
v_pl=v2;
v_inf_i=v2_L1-v_pl;
v_inf_f=v1_L2-v_pl;
a_h1=-mu_Io/((norm(v_inf_i)^2));
a_h2=-mu_Io/((norm(v_inf_f)^2));
DELTA=acos(dot(v_inf_i,v_inf_f)/(norm(v_inf_i)*norm(v_inf_f))); %[rad]        DELTA=delta1+delta2
R_imp=r_i+5;
f=@(d1) 1-(a_h2/a_h1)*(1-1./sin(DELTA-d1))-(1./sin(d1));
[del,it]=bisez(0,DELTA,1e-8,f);
d1_ga=del(end);
d2_ga=DELTA-d1_ga;
e_h1=1/sin(d1_ga);
r_p_ga=a_h1*(1-e_h1);
if r_p_ga<R_imp;
    disp('POWERED GRAVITY ASSIST MANUEVER IMPACTS ON IO SOIL')
end
e_h2=1-(r_p_ga/a_h2);
p_h1=a_h1*(1-(e_h1^2));
p_h2=a_h2*(1-(e_h2^2));
vp_h1=sqrt(mu_Io/p_h1)*(1+e_h1);
vp_h2=sqrt(mu_Io/p_h2)*(1+e_h2);

cost_powered_ga=abs(vp_h2-vp_h1);
DV_GA=norm(v_inf_f-v_inf_i)-cost_powered_ga;

% Powered gravity assist plot
h_vers=cross(v_inf_i,v_inf_f)/norm(cross(v_inf_i,v_inf_f));
i_h=acos(dot(K,h_vers))*(180/pi); %[deg]
n_vers=cross(K,h_vers);
if dot(n_vers,J)>0;
   OM_h=acos(dot(n_vers,I))*(180/pi); %[deg]
else 
    OM_h=(2*pi-acos(dot(n_vers,I)))*(180/pi); %[deg]
end
h_vect=(r_p_ga*vp_h1)*h_vers;
e_vect=(cross(v_inf_i,h_vect)/mu_Io)+v_inf_i/norm(v_inf_i);
if dot(e_vect,K)>0
    om_h=acos(dot(n_vers,e_vect)/norm(e_vect))*(180/pi); %[deg]
else 
    om_h=(2*pi-acos(dot(n_vers,e_vect)/norm(e_vect)))*(180/pi); %[deg]
end
theta_inf=round(((3*pi)/2-d1_ga)*(180/pi))+1;
theta_INF=round((pi/2+d2_ga)*(180/pi))-1;
r_h1=ones(3,length(theta_inf:360));
v_h1=ones(3,length(theta_inf:360));
for k=theta_inf:360
[r_h1(:,k+1-theta_inf),v_h1(:,k+1-theta_inf)]=OrbPar2RV(a_h1,e_h1,i_h,OM_h,om_h,k,mu_Io);
end
r_h2=ones(3,length(1:theta_INF));
v_h2=ones(3,length(1:theta_INF));
for k=1:theta_INF+1
[r_h2(:,k),v_h2(:,k)]=OrbPar2RV(a_h2,e_h2,i_h,OM_h,om_h,k-1,mu_Io);
end
figure()
plot3(r_h1(1,:),r_h1(2,:),r_h1(3,:),'--r',r_h2(1,:),r_h2(2,:),r_h2(3,:),'--b',r_h2(1,1),r_h2(2,1),r_h2(3,1),'*k',0,0,0,'ok')
legend('Hyp_I_N','Hyp_O_U_T','x_i_m_p_u_l_s_e','Io')
grid on

figure()
plot3(r_h1(1,:),r_h1(2,:),r_h1(3,:),'--r',r_h2(1,:),r_h2(2,:),r_h2(3,:),'--b',r_h2(1,1),r_h2(2,1),r_h2(3,1),'*k')
hold on
surf(x*r_i,y*r_i,z*r_i)
colormap('spring')
axis([-1.2e4 1.2e4 -1.2e4 1.2e4 -1.2e4 1.2e4])
legend('Hyp_I_N','Hyp_O_U_T','x_i_m_p_u_l_s_e','Io')
grid on




%% TARGET TRAJECTORY APPROACH
% Final orbit data.
h_p=120;
h_a=200;
rp_f=h_p+r_e;
ra_f=h_a+r_e;
a_f=(ra_f+rp_f)/2;
e_f=1-(rp_f/a_f);
p_f=a_f*(1-e_f^2);
i_f=78;
OM_f=250;
om_f=80;
[RP_F,VP_F]=OrbPar2RV(a_f,e_f,i_f,OM_f,om_f,0,mu_Eu);
T_f=2*pi*sqrt(a_f^3/mu_Eu);
h_vers_f=cross(RP_F,VP_F)/(norm(cross(RP_F,VP_F)));
n_vers_f=cross(K,h_vers_f);

% 1st impulse: Circularization of the hyperbola
rp_in=rp_f;
vc_1=sqrt(mu_Eu/rp_f);
v_inf=v2_L2'-v3';
a_in=-(mu_Eu)/(norm(v_inf)^2);
e_in=1-(rp_in/a_in);
p_in=a_in*(1-e_in^2);
vp_in=sqrt(mu_Eu/p_in)*(1+e_in);
h_vers_in=cross(n_vers_f,v_inf)/(norm(cross(n_vers_f,v_inf)));
h_vect_in=(rp_in*vp_in)*h_vers_in;
i_in=acos(dot(h_vers_in,K))*(180/pi);
e_vect_in=(cross(v_inf,h_vect_in)/mu_Eu)+v_inf/(norm(v_inf));
e_vers_in=e_vect_in/(norm(e_vect_in));
if dot(e_vect_in,K)>0
    om_in=acos(dot(n_vers_f,e_vect_in)/norm(e_vect_in))*(180/pi); %[deg]
else 
    om_in=(2*pi-acos(dot(n_vers_f,e_vect_in)/norm(e_vect_in)))*(180/pi); %[deg]
end
om_IN=acos(dot(n_vers_f,e_vers_in))*(180/pi);
OM_in=OM_f;
a_c1=rp_f;
e_c1=0;
i_c1=i_in;
OM_c1=OM_in;
om_c1=0;

cost_c1=abs(vp_in-vc_1);

% 2nd impulse: Changin plane
a_c2=a_c1;
e_c2=0;
i_c2=i_f;
OM_c2=OM_c1;
om_c2=0;

cost_c2=abs(2*vc_1*sin((i_f-i_in)/2*(pi/180)));

% 3rd impulse: Immission into final orbit
cost_f=abs(norm(VP_F)-vc_1);

% Immision final orbit plot
theta_in=round(((3*pi)/2-asin(1/e_in))*(180/pi))+10;
r_in=ones(3,length(theta_in:360));
for k=theta_in:360
[r_in(:,k+1-theta_in)]=OrbPar2RV(a_in,e_in,i_in,OM_in,om_in,k,mu_Eu);
end

theta_c1=om_in;          %[deg]
s_2nd_imp=sqrt(a_c1^3/mu_Eu)*((360-theta_c1)*pi/180);    %[s]
theta_c2=om_f;   %[deg]
s_3rd_imp=sqrt(a_c1^3/mu_Eu)*(theta_c2*pi/180);  %[s]
step=100;
it_c1=linspace(theta_c1,360,step);
it_c2=linspace(0,theta_c2,step);
it_f=linspace(0,359,step);
r_c1=ones(3,step);
r_c2=ones(3,step);
r_f=ones(3,step);
for k=1:step
[r_c1(:,k)]=OrbPar2RV(a_c1,e_c1,i_c1,OM_c1,om_c1,it_c1(k),mu_Eu);
[r_c2(:,k)]=OrbPar2RV(a_c2,e_c2,i_c2,OM_c2,om_c2,it_c2(k),mu_Eu);
[r_f(:,k)]=OrbPar2RV(a_f,e_f,i_f,OM_f,om_f,it_f(k),mu_Eu);
end
    
figure()
plot3(r_in(1,:),r_in(2,:),r_in(3,:),'--m',r_c1(1,:),r_c1(2,:),r_c1(3,:),'r',r_c2(1,:),r_c2(2,:),r_c2(3,:),'b',r_f(1,:),r_f(2,:),r_f(3,:),'g',0,0,0,'ok',r_c1(1,1),r_c1(2,1),r_c1(3,1),'*m',r_c1(1,end),r_c1(2,end),r_c1(3,end),'*r',r_f(1,1),r_f(2,1),r_f(3,1),'*b')
legend('Hyp_I_N','Circ_o_r_b_1','Circ_o_r_b_2','Operational_o_r_b','Europa','3rd impulse','4th impulse','5th impulse')
grid on

figure()
plot3(r_in(1,:),r_in(2,:),r_in(3,:),'--m',r_c1(1,:),r_c1(2,:),r_c1(3,:),'r',r_c2(1,:),r_c2(2,:),r_c2(3,:),'b',r_f(1,:),r_f(2,:),r_f(3,:),'g',r_c1(1,1),r_c1(2,1),r_c1(3,1),'*m',r_c1(1,end),r_c1(2,end),r_c1(3,end),'*r',r_f(1,1),r_f(2,1),r_f(3,1),'*b')
hold on
surf(x*r_e,y*r_e,z*r_e)
colormap('cool')
axis([-0.3e4 0.3e4 -0.3e4 0.3e4 -0.3e4 0.3e4])
legend('Hyp_I_N','Circ_o_r_b_1','Circ_o_r_b_2','Operational_o_r_b','3rd impulse','4th impulse','5th impulse','Europa')
grid on

TOTAL_COST=cost_L+cost_powered_ga+cost_c1+cost_c2+cost_f;

%% RESULTS PRINT
sprintf('\t\t\t TRANSFERT ORBITS PARAMETERS \n\n \t\t\t\t    a \t\t   e \t\t   i \t\t   OM \t\t   om \n Lambert transfert GA-IO: \t %6.1f [km] \t %1.4f [-] \t %2.3f [deg] \t %3.2f [deg] \t %3.2f [deg] \n Incoming hyperbola on Io: \t %3.2f [km] \t %2.4f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg] \n Outcoming hyperbola on Io: \t %3.2f [km] \t %2.4f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg] \n Lambert transfert IO-EU: \t %6.1f [km] \t %1.4f [-] \t %2.3f [deg] \t %3.2f [deg] \t %3.2f [deg] \n Incoming hyperbola on Europa: \t %3.2f [km] \t %1.4f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg] \n 1st Circular orbit on Europa: \t %4.2f [km] \t %1.1f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg] \n 2nd Circular orbit on Europa: \t %4.2f [km] \t %1.1f [-] \t %3.2f [deg] \t %3.2f [deg] \t %3.2f [deg]',a_L1,e_L1,i_L1,OM_L1,om_L1,a_h1,e_h1,i_h,OM_h,om_h,a_h2,e_h2,i_h,OM_h,om_h,a_L2,e_L2,i_L2,OM_L2,om_L2,a_in,e_in,i_in,OM_in,om_in,a_c1,e_c1,i_c1,OM_c1,om_c1,a_c2,e_c2,i_c2,OM_c2,om_c2)

sprintf('\t\t\t TRANSFERT COSTS \n \t\t\t on board propulsion \n\n AROUND GANYMEDE: \n Departure cost \t\t %3.3f km/s \n\n AROUND IO: \n Powered gravity assist cost \t %3.3f km/s \n\n AROUND EUROPA: \n Breaking cost \t\t\t %3.3f km/s \n Change of plane cost \t\t %3.3f km/s \n Immision to final orbit cost \t %3.3f km/s \n\n TOTAL COST \t\t\t %3.3f km/s \n\n \t\t\t gravity assisted contributions \n AROUND IO: \n gravity assist contribution \t %3.3f km/s',cost_L,cost_powered_ga,cost_c1,cost_c2,cost_f,TOTAL_COST,DV_GA)

[YEAR_0,MONTH_0,DAY_0,HOUR_0,MINUTE_0,SECOND_0]=departure_date(t0);
t_ga=t0+tof_1;
[YEAR_ga,MONTH_ga,DAY_ga,HOUR_ga,MINUTE_ga,SECOND_ga]=departure_date(t_ga);
t_1st_imp=t0+tof_1+tof_2;
[YEAR_1st,MONTH_1st,DAY_1st,HOUR_1st,MINUTE_1st,SECOND_1st]=departure_date(t_1st_imp);
t_2nd_imp=t_1st_imp+s_2nd_imp;
[YEAR_2nd,MONTH_2nd,DAY_2nd,HOUR_2nd,MINUTE_2nd,SECOND_2nd]=departure_date(t_2nd_imp);
t_3rd_imp=t_2nd_imp+s_3rd_imp;
[YEAR_3rd,MONTH_3rd,DAY_3rd,HOUR_3rd,MINUTE_3rd,SECOND_3rd]=departure_date(t_3rd_imp);

sprintf('\t\t\t MANEUVERING DATES \n\n DEPARTURE DATE \t\t\t %d:%d:%d   %d/%d/%d \n\n Powered gravity assist maneuver date \t %d:%d:%d   %d/%d/%d \n Breaking maneuver date \t\t %d:%d:%d   %d/%d/%d \n Change of plane maneuver date \t\t %d:%d:%d %d/%d/%d \n\n ARRIVAL DATE \t\t\t\t %d:%d:%d  %d/%d/%d',HOUR_0,MINUTE_0,SECOND_0,DAY_0,MONTH_0,YEAR_0,HOUR_ga,MINUTE_ga,SECOND_ga,DAY_ga,MONTH_ga,YEAR_ga,HOUR_1st,MINUTE_1st,SECOND_1st,DAY_1st,MONTH_1st,YEAR_1st,HOUR_2nd,MINUTE_2nd,SECOND_2nd,DAY_2nd,MONTH_2nd,YEAR_2nd,HOUR_3rd,MINUTE_3rd,SECOND_3rd,DAY_3rd,MONTH_3rd,YEAR_3rd)


 save('OPTIMUM_data.mat','t0','tof_1','tof_2');