%% GENERATION FILE DATA
clear all
close all


% Import data
[A_Ju]=importdata('data_Jupiter.txt');
[A_Ga]=importdata('data_Ganymede.txt');
[A_Eu]=importdata('data_Europa.txt');
[A_Io]=importdata('data_Io.txt');
j=1;
k=1;
l=length(A_Ga);
R_Ju_data=ones(length(2:3:l-1),3);
V_Ju_data=ones(length(2:3:l-1),3);
R_Ga_data=ones(length(2:3:l-1),3);
V_Ga_data=ones(length(2:3:l-1),3);
R_Eu_data=ones(length(2:3:l-1),3);
V_Eu_data=ones(length(2:3:l-1),3);
R_Io_data=ones(length(2:3:l-1),3);
V_Io_data=ones(length(2:3:l-1),3);

% define data matrix
for i=2:3:l-1
  R_Ju_data(j,:)=str2num(cell2mat(A_Ju(i)));
  R_Ga_data(j,:)=str2num(cell2mat(A_Ga(i)));
  R_Eu_data(j,:)=str2num(cell2mat(A_Eu(i)));
  R_Io_data(j,:)=str2num(cell2mat(A_Io(i)));
  j=j+1;
 end


 
 for i=3:3:l
  V_Ju_data(k,:)=str2num(cell2mat(A_Ju(i)));
  V_Ga_data(k,:)=str2num(cell2mat(A_Ga(i)));
  V_Eu_data(k,:)=str2num(cell2mat(A_Eu(i)));
  V_Io_data(k,:)=str2num(cell2mat(A_Io(i)));
  k=k+1;
 end

 % generation file
 save('HORIZON_data.mat','R_Ju_data','V_Ju_data','R_Ga_data','V_Ga_data','R_Eu_data','V_Eu_data','R_Io_data','V_Io_data');
