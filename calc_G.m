function [ G ] = calc_G(theta_s,phi_s,shape,d_theta)

%Çó½â·ÖÄ¸
if nargin==3
    d_theta=0.1;
    d_phi=0.1; % if this value larger than 0.03  the result will be influented
end


G=0;
cnt=0;
if strcmp(shape,'dihedral')
   for theta_d=pi/4:d_theta:3*pi/4
       for phi_d=-pi/2:d_phi:pi/2
        cnt=cnt+1;
%          fprintf('calculate G: %d \\ %d \n',cnt,fix(pi/2/d_theta*pi*2/d_phi));
        G=G+sin(theta_d)*calc_N_I_s( theta_s,phi_s,shape,theta_d,phi_d);  
       end
   end 
else   
   for theta_d=0:d_theta:pi/2
       for phi_d=0:d_phi:pi*2
        cnt=cnt+1;
%          fprintf('calculate G: %d \\ %d \n',cnt,fix(pi/2/d_theta*pi*2/d_phi));
        G=G+sin(theta_d)*calc_N_I_s( theta_s,phi_s,shape,theta_d,phi_d);  
       end
   end
end 

G=G*d_phi*d_theta;
% fprintf('theta_s: %f  G: %d \f \n',theta_s,G);

end


