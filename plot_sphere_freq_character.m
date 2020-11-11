
clear;

R=0.1;
c=3.0e+8;

ii=0;
l_ka=0:1:10;
l_NIs=zeros(size(l_ka));
frequency_G=zeros(size(l_ka));
for ka=l_ka
    ii=ii+1;
    k=ka/R;
    lambda=2*pi/k;
    frequency_G(ii)=c/lambda/10^9;
    
    allsum=pi*R*1.0i/k*(1-exp(1.0i*2*k*R));
%     allsum=sin(k*R);
    NIs=(abs(allsum))^2/4*pi;
    l_NIs(ii)=10*log10(NIs);
   
end
plot(l_ka,l_NIs,'-')

 
grid_size=0.001;
fileformat='gridpro';
% fileformat='catia';
child_dirname='sphere-r-0.1';

QRCS_program_path=''; 
ISIGHT_output_path='model';

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

step=0.5;
theta_begin=0;
theta_end=0;
phi_begin=0;
phi_end=0;
l_G=calculate_3D(dat_filename,frequency_G,step,...
     theta_begin,theta_end,phi_begin,phi_end,fileformat,0);

l_G=10*log10(l_G);
plot(l_ka,l_NIs,'-')
hold on;
plot(l_ka,l_G,'-')
hold on;
plot(l_ka,l_NIs-l_G,'-')
xlabel('ka');
ylabel('(dB/m^2)');  
legend('NIs','G','NIs/G','Location','northeast')
grid on;












