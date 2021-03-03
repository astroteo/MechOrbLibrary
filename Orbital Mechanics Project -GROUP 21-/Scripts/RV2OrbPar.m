function [a,em,i,OM,om,theta]=RV2OrbPar(r,v,Mass_or_mu)
mu_t=398600;
M_t=5.9736e24;
if nargin==2
    Mass_or_mu=mu_t;
end
if  Mass_or_mu>1e13
    mu=(Mass_or_mu/M_t)*mu_t;
else
    mu=Mass_or_mu;
end

rm=norm(r);
vm=norm(v);
a=mu/(((2*mu)/rm)-(vm^2));
h=cross(r,v);
hm=norm(h);
l=cross(v,h);
e=l/mu-r/rm;
em=norm(e);
i=acos(h(3)/hm);
k=[0;0;1];
j=[0;1;0];
iv=[1;0;0];
n=cross(k,h)/norm(cross(k,h));
if dot(n,j)>0
   OM=acos(dot(n,iv));
else
    OM=2*pi-acos(dot(n,iv));
end
if dot(e,k)>0
    om=acos(dot(n,e)/em);
else
    om=2*pi-acos(dot(n,e)/em);
end
if dot(r,v)>0
    theta=acos(dot(r,e)/(rm*em));
else
    theta=2*pi-acos(dot(r,e)/(rm*em));
end
i=(180*i)/pi;
OM=(180*OM)/pi;
om=(180*om)/pi;
theta=(180*theta)/pi;