%% OPTIMUM RESEARCH
clear all
close all


%% DATA
% GANYMEDE-EUROPA-IO data
load('HORIZON_data.mat')

% Reference axis
I=[1 0 0];
J=[0 1 0];
K=[0 0 1];

% Planetary gravity constants
mu_J=126711995;
mu_Io=((8.9319E+22)/(5.9736E+24))*398600;
mu_Eu=3.1888e3;

% Final orbit data
r_Eu=1560.8;
h_p=120;
h_a=200;
rp_f=h_p+r_Eu;
ra_f=h_a+r_Eu;
a_f=(ra_f+rp_f)/2;
e_f=1-(rp_f/a_f);
p_f=a_f*(1-e_f^2);
i_f=78;
OM_f=250;
om_f=80;
[RP_F,VP_F]=OrbPar2RV(a_f,e_f,i_f,OM_f,om_f,0,mu_Eu);
T_f=2*pi*sqrt(a_f^3/mu_Eu);
h_vers_f=cross(RP_F,VP_F)./(norm(cross(RP_F,VP_F)));
n_vers_f=cross(K,h_vers_f);

% Iteration data
N_h_t0=input('\nDEPARTURE DATE TIME STEP [h]:\n(must be a multiple of two hours)\n(2<step<52560)\n');                                                               % N_h = cycle step in terms of hours (must be pair)
N_h_tof=input('\nTIME OF FLIGHT TIME STEP [h]:\n(must be a multiple of two hours)\n(2<step<96)\n');
t0_it=(0:3.6e3*N_h_t0:5*365*24*3600+24*3600);
tof_it=(57.6e3:3.6e3*N_h_tof:4*24*3600);
TOTAL_COST=zeros(length(t0_it),length(tof_it),length(tof_it));

%% RESEARCH OPTIMAL MANUEVER

for i=1:length(t0_it)
    t0=t0_it(i);
    for j=1:length(tof_it)
        tof_1=tof_it(j);
        for k=1:length(tof_it)
            tof_2=tof_it(k);
            
%% LAMBERT TRANSFERS
% First transfer arc.
r1=R_Ga_data(t0/7.2e3+1,:);
v1=V_Ga_data(t0/7.2e3+1,:);
t_1=t0+tof_1;
r2=R_Io_data(t_1/7.2e3+1,:);
v2=V_Io_data(t_1/7.2e3+1,:);
[a_t1,p_t1,e_t1,error1,v1_t1,v2_t1,theta1]=lambert_mick(r1,r2,tof_1,mu_J);

cost_L=norm(v1_t1-v1);

% Second transfer arc
t_2=t0+tof_1+tof_2;
r3=R_Eu_data(t_2/7.2e3+1,:);
v3=V_Eu_data(t_2/7.2e3+1,:);
[a_t2,p_t2,e_t2,error2,v1_t2,v2_t2,theta2]=lambert_mick(r2,r3,tof_2,mu_J);




%% POWERED GRAVITY ASSIST

vinf_i=v2_t1-v2;
vinf_f=v1_t2-v2;
a_h1=-mu_Io/((norm(vinf_i)^2));
a_h2=-mu_Io/((norm(vinf_f)^2));
DELTA=acos(dot(vinf_i,vinf_f)/(norm(vinf_i)*norm(vinf_f))); %[rad]        DELTA=delta1+delta2
r_Io=(3642.6/2);
R_imp=r_Io+0.5;
f=@(d1) 1-(a_h2/a_h1)*(1-1./sin(DELTA-d1))-(1./sin(d1));
[x,it]=bisez(0,DELTA,1e-8,f);
d1_ga=x(end);
d2_ga=DELTA-d1_ga;
e_h1=1/sin(d1_ga);
r_p_ga=a_h1*(1-e_h1);
e_h2=1-(r_p_ga/a_h2);
p_h1=a_h1*(1-(e_h1^2));
p_h2=a_h2*(1-(e_h2^2));
vp_h1=sqrt(mu_Io/p_h1)*(1+e_h1);
vp_h2=sqrt(mu_Io/p_h2)*(1+e_h2);
if r_p_ga<R_imp;
    cost_powered_ga=inf;
else
    cost_powered_ga=abs(vp_h2-vp_h1);
end



%% TARGET TRAJECTORY APPROACH

% 1st impulse: Circularization of the hyperbola
rp_in=rp_f;
vc_1=sqrt(mu_Eu/rp_f);
v_inf=v2_t2-v3;
a_in=-(mu_Eu)/(norm(v_inf)^2);
e_in=1-(rp_in/a_in);
p_in=a_in*(1-e_in^2);
vp_in=sqrt(mu_Eu/p_in)*(1+e_in);
h_vers_in=cross(n_vers_f,v_inf)/(norm(cross(n_vers_f,v_inf)));
i_in=acos(dot(h_vers_in,K))*(180/pi);

cost_c1=abs(vp_in-vc_1);

% 2nd impulse: Changin plane
di=(i_f-i_in)*pi/180;
cost_c2=abs(2*vc_1*sin(di/2));

% 3rd impulse: Immission into final orbit
cost_f=abs(norm(VP_F)-vc_1);


TOTAL_COST(i,j,k)=cost_L+cost_powered_ga+cost_c1+cost_c2+cost_f;
        end
    end
end


%% DATE OPTIMUM VALUE
opt_cost=min(min(min(TOTAL_COST)))
for u=1:1:length(t0_it)
for t=1:1:length(tof_it)
for r=1:1:length(tof_it)
if TOTAL_COST(u,t,r)==opt_cost
i_opt=u;
j_opt=t;
k_opt=r;
end
end
end
end
t0=t0_it(i_opt)
tof_1=tof_it(j_opt)
tof_2=tof_it(k_opt)