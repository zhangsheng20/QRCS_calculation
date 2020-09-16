function [ sigma_Q ] = calc_sigma_Q(theta,shape)

%���ۼ������Ŀ�����sigma_Q


c=3.0e+8;
lambda=0.25;
omega=2*pi*c/lambda;
k=2*pi/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'rect')
%����ƽ������sigma_Q
a=1;
b=1;
sigma_Q=4*pi.*((a.*b).^2)/(lambda^2).* abs(cos(theta)).*( sin(k*a.*sin(theta))./(k*a.*sin(theta))   ).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'circle')
%Բ�ΰ�
R=1;
A=pi*R^2;
sigma_Q=16*pi*A^2/(lambda^2)*abs(cos(theta))*(besselj(1,2*k*R*sin(theta))/(2*k*R*sin(theta)))^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(shape,'triangle')
%�����ΰ�
h=1.0;
b=1.0;

sigma_Q=pi*b^2*h^2/lambda^2*abs(cos(theta))*(sinc(k*b/pi/2*sin(theta)))^4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

