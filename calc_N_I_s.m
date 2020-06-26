function [ N_I_s ] = calc_N_I_s(  theta_s,phi_s,shape,theta_d,phi_d)
%一个一个原子计算
if nargin==0
  theta_s=0;
  phi_s=0;
  theta_d=0;
  phi_d=0; 

elseif nargin==3 
  theta_d=theta_s;
  phi_d=phi_s;
    else
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%公共量
c=3.0e+8;
lambda=0.25;
omega=2*pi*c/lambda;
delta=0.1*lambda;
r=100000000*lambda;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'rect')
%矩形板

a=1;
b=1;

%r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)];
%r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
allsum=0;
for xi=-a/2:delta:a/2
    for yi=-b/2:delta:b/2
       if nargin==3
            Delta_Ri=2*sqrt(((xi-r.*sin(theta_s).*cos(phi_s)).^2+(yi-r.*sin(theta_s).*sin(phi_s)).^2+(r.*cos(theta_s)).^2));
           %Delta_Ri=2*norm(r_d-r_i); 
       else
           %r_i=[xi,yi,0];
           %Delta_Ri=norm(r_s-r_i)+norm(r_d-r_i);  %norm算得太慢
           Delta_Ri=sqrt((r.*sin(theta_s).*cos(phi_s)-xi).^2+(r.*sin(theta_s).*sin(phi_s)-yi).^2+(r.*cos(theta_s)).^2)...
                     +sqrt((r*sin(theta_d).*cos(phi_d)-xi).^2+(r.*sin(theta_d).*sin(phi_d)-yi).^2+(r.*cos(theta_d)).^2);
       end
       allsum=allsum+exp(1.000i*omega*Delta_Ri/c);
    end
end
N_I_s=(abs(allsum)).^2;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'circle')
%圆形板

% %% 新的
% R=1;
% ii=0;
% coordinate=zeros(fix(4*R^2/delta/delta),2);
% for xi=-R:delta:R
%     for yi=-sqrt(R^2-xi^2):delta:sqrt(R^2-xi^2)
%        ii=ii+1;
%        coordinate(ii,:)=[xi yi];
%     end 
% end
% coordinate((ii+1:length(coordinate)),:)=[];
% Mat_Delta_Ri=2*sqrt(((coordinate(:,1)-r.*sin(theta_s).*cos(phi_s)).^2+(coordinate(:,2)-r.*sin(theta_s).*sin(phi_s)).^2+(r.*cos(theta_s)).^2));
% d=exp(1.00000i*omega*Mat_Delta_Ri/c);
% allsum= sum(d);
% N_I_s=(abs(allsum)).^2;


% 旧的
R=1;
%r_d=[r.*sin(theta).*cos(phi),r.*sin(theta).*sin(phi),r.*cos(theta)];
allsum=0;
for xi=-R:delta:R
    for yi=-sqrt(R^2-xi^2):delta:sqrt(R^2-xi^2)
         if nargin==3
            Delta_Ri=2*sqrt(((xi-r.*sin(theta_s).*cos(phi_s)).^2+(yi-r.*sin(theta_s).*sin(phi_s)).^2+(r.*cos(theta_s)).^2));    
         else
           Delta_Ri=sqrt((r.*sin(theta_s).*cos(phi_s)-xi).^2+(r.*sin(theta_s).*sin(phi_s)-yi).^2+(r.*cos(theta_s)).^2)...
                     +sqrt((r*sin(theta_d).*cos(phi_d)-xi).^2+(r.*sin(theta_d).*sin(phi_d)-yi).^2+(r.*cos(theta_d)).^2);
         end
       allsum=allsum+exp(1.000i*omega*Delta_Ri/c);
    end
end
N_I_s=(abs(allsum)).^2;

 
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'triangle')
%三角形板

h=1;
b=1;

%r_d=[r.*sin(theta).*cos(phi),r.*sin(theta).*sin(phi),r.*cos(theta)];
allsum=0;
for yi=0:delta:h
    for xi=b/2/h*(yi-h):delta:-b/2/h*(yi-h)   
         if nargin==3
            Delta_Ri=2*sqrt(((xi-r.*sin(theta_s).*cos(phi_s)).^2+(yi-r.*sin(theta_s).*sin(phi_s)).^2+(r.*cos(theta_s)).^2));    
         else
           Delta_Ri=sqrt((r.*sin(theta_s).*cos(phi_s)-xi).^2+(r.*sin(theta_s).*sin(phi_s)-yi).^2+(r.*cos(theta_s)).^2)...
                     +sqrt((r*sin(theta_d).*cos(phi_d)-xi).^2+(r.*sin(theta_d).*sin(phi_d)-yi).^2+(r.*cos(theta_d)).^2);
         end
       allsum=allsum+exp(1.000i*omega*Delta_Ri/c);
    end
end

N_I_s=(abs(allsum)).^2;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'dihedral')
%二面角
%公共量
c=3.0e+8;
lambda=0.25;
omega=2*pi*c/lambda;
delta=0.1*lambda;
r=100000*lambda;


a=0.5;
b=0.5;
l_x=delta:delta:a;
l_y=-b/2:delta:b/2;
Mat_x=ones(length(l_y),1)*l_x;
Mat_y=(ones(length(l_x),1)*l_y)';
MatA=zeros(length(l_x)*length(l_y),3);
MatB=MatA;

MatA(:,1)=Mat_x(:);
MatA(:,2)=Mat_y(:);
MatB(:,3)=Mat_x(:);
MatB(:,2)=Mat_y(:);

ra=pi/4;
MatRotate=[cos(ra) 0 sin(ra);0 1 0;-sin(ra) 0 cos(ra)];
MatA=MatA*MatRotate';
MatB=MatB*MatRotate';
r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)];
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
allsum=0;

for ii=1:length(MatA)
    for kk=1:length(MatB)  
        r_i=MatA(ii,:);
        r_k=MatB(kk,:);
        Delta_Ri=norm(r_s-r_i)+norm(r_k-r_d)+norm(r_i-r_k);
      
%         Delta_Ri=sqrt((r.*sin(theta_s).*cos(phi_s)-MatA(ii,1)).^2+(r.*sin(theta_s).*sin(phi_s)-MatA(ii,2)).^2+(r.*cos(theta_s)).^2)...
%                  +sqrt((r*sin(theta_d)*cos(phi_d)).^2+(r*sin(theta_d)*sin(phi_d)-MatB(kk,2)).^2+(r*cos(theta_d)-MatB(kk,3)).^2)...
%                  +sqrt((MatA(ii,1)-MatB(kk,1)).^2+(MatA(ii,2)-MatB(kk,2)).^2+(MatA(ii,3)-MatB(kk,3)).^2);
     
    dd=exp(1.000i*omega*Delta_Ri/c);
    allsum=allsum+dd;
    end   
end

N_I_s=(abs(allsum)).^2;

end
end

