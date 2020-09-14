
%QRCS_program_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\QRCS_calculation';
%ISIGHT_output_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\isight_missile_qrcs_optimization\output';
%child_dirname=1;
addpath(QRCS_program_path);
my_str=num2str(child_dirname);
output_path=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'.dat');
step=0.1;
theta_begin=pi/2;
theta_end=pi/2;
phi_begin=-pi/12;
phi_end=pi/12;
% fileformat='catia';
fileformat='gridpro';
%filename='C:\Users\Administrator\Desktop\quantum stealth\mygit\dom-1.dat'
calculate_3D(output_path,step,theta_begin,theta_end,phi_begin,phi_end,fileformat)