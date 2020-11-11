function [sigama_Q] = calculate_sphere_G2(k,R)


if(nargin<2)
    k=8*pi;
    R=0.1;
end
a=calculate_integral(k,R,0)

l_theta_p=0:0.01:pi;
l_out=zeros(size(l_theta_p));
l_NIs=zeros(size(l_theta_p));
for ii=1:length(l_theta_p)
    theta_p=l_theta_p(ii);
    l_NIs(ii)=calculate_integral(k,R,theta_p);  
    l_out(ii)=l_NIs(ii)*sin(theta_p);
end
  plot(l_theta_p,10*log10(l_NIs),'-')
  grid on;

G=sum(l_out);
% sigama_Q=4*pi *pi*R^2 *(2*pi*R/k*sin(k*R))^2/G;

end

function [NIs] = calculate_integral(k,R,theta_p)
fun = @(theta,phi) exp(1i*k*R*(-sin(theta_p).*sin(theta).*cos(phi)+(-1-cos(theta_p))*cos(theta)))...
     .*abs(sin(theta))*R^2;

% fun = @(theta,phi) R^2*abs(sin(theta))+0*theta.*phi
    thetamin = @(phi) -atan(1./cos(phi)/tan(theta_p));
    thetamax =  pi/2;
    phimin = -pi/2;
    phimax = pi/2;

q= integral2(fun,phimin,phimax,thetamin,thetamax');
NIs=(abs(q))^2;
end

