
%QRCS_program_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\QRCS_calculation';
%ISIGHT_output_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\isight_missile_qrcs_optimization\output';
%child_dirname=1;
clear;
grid_size=0.01;
% fileformat='gridpro';
fileformat='catia';
if(exist('child_dirname','var')==0)
   child_dirname='sphere-r-0.5';
end
if(exist('QRCS_program_path','var')==0)
   QRCS_program_path='C:\Users\Administrator\Desktop\quantum stealth\mygit'; 
end
if(exist('ISIGHT_output_path','var')==0)
   %ISIGHT_output_path='C:\Users\Administrator\Desktop\QRCS_optimization_ws\isight_missile_qrcs_optimization\output';
    ISIGHT_output_path='C:\Users\Administrator\Desktop\quantum stealth\mygit\model';
end

addpath(QRCS_program_path);
if( ~isa(child_dirname,'string')||~isa(child_dirname,'char') )
    my_str=num2str(child_dirname);
end

if(~exist('grid_size','var'))
    dat_filename=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'.dat');
else
    if(strcmp(fileformat,'gridpro'))
        dat_filename=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'-grid-',num2str(grid_size),'-gridpro','.dat');
    else
        dat_filename=strcat(ISIGHT_output_path,'\',my_str,'\',my_str,'-grid-',num2str(grid_size),'.dat')
    end
end
step=0.1;
theta_begin=0;
theta_end=pi/2;
phi_begin=0;
phi_end=0;
% fileformat='catia';
%fileformat='gridpro';
%filename='C:\Users\Administrator\Desktop\quantum stealth\mygit\dom-1.dat'
calculate_3D(dat_filename,step,theta_begin,theta_end,phi_begin,phi_end,fileformat,1)

output_mat_path = strrep(dat_filename,'.dat','.mat');
load(output_mat_path);


sigma_Q_1=10.^(l_sigma_Q_1/10);
value_mean1=mean(sigma_Q_1);
mean_QRCS=10*log10(abs(value_mean1));















