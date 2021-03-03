function [data_at_t]=interp_data(data,N_h,t)

%R,V represents the matrices containig position and velocity of the body at
%all the time instants computed by Horizon.
%N_h represents the step choosen  to compute the state in hours
%!!!MUST BE THE SAME TIME STEP IMPOSED IN HORIZON's FILE.TXT!!!

bottom_h=floor(t/(N_h*3600));
top_h=ceil(t/(N_h*3600));

tiempo=t-(bottom_h*N_h*3600);


data_at_t=data(bottom_h+1,:)+((data(top_h+1,:)-data(bottom_h+1,:))/(N_h*3600))*tiempo;
data_at_t=data_at_t';