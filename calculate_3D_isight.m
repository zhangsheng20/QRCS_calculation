
%QRCS_program_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\QRCS_calculation';
%ISIGHT_output_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\isight_missile_qrcs_optimization\output';
%child_dirname=1;
if(exist('child_dirname','var')==0)
   child_dirname=1;
end
if(exist('QRCS_program_path','var')==0)
   QRCS_program_path='C:\Users\Administrator\Desktop\quantum stealth\mygit'; 
end
if(exist('ISIGHT_output_path','var')==0)
   ISIGHT_output_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\isight_missile_qrcs_optimization\output';
end

addpath(QRCS_program_path);
my_str=num2str(child_dirname);
dat_filename=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'.dat');
step=0.2;
theta_begin=pi/2;
theta_end=pi/2;
phi_begin=0;
phi_end=pi/6;
% fileformat='catia';
fileformat='gridpro';
%filename='C:\Users\Administrator\Desktop\quantum stealth\mygit\dom-1.dat'
calculate_3D(dat_filename,step,theta_begin,theta_end,phi_begin,phi_end,fileformat,0)

output_mat_path=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'.mat');
load(output_mat_path);


sigma_Q_1=10.^(l_sigma_Q_1/10);
value_mean1=mean(sigma_Q_1);
mean_QRCS=10*log10(abs(value_mean1));















