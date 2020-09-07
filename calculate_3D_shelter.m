function calculate_3D_shelter

route={'.\irregular_1m_1m_0.5m_0.03m.dat'};
data_process(route);


l_theta=-pi/2:0.1:pi/2;
l_sigma_Q_1=zeros(size(l_theta));
l_sigma_Q_2=zeros(size(l_theta));

ii=0;
for theta= l_theta
    ii=ii+1;
    fprintf('progress: %d \\ %d \n',ii,length(l_theta));
    %A=calc_A_3D(theta,0);
    %l_sigma_Q_1(ii)=(calc_G_3D(theta,0));
    l_sigma_Q_1(ii)=10*log10(calc_3D_shelter_N_I_s(pi/6,0,theta,0));
    l_sigma_Q_2(ii)=10*log10(calc_3D_N_I_s(pi/6,0,theta,0));
    
    
end
figure(2);
plot(l_theta,l_sigma_Q_1,l_theta,l_sigma_Q_2)
legend('have shelter','no shelter','Location','northeast')
xlabel('\theta(rad)');
ylabel('\sigma_Q(dB/m^2)');  


end




%% 
function [ output ] = calc_3D_I_s( theta_s,phi_s,theta_d,phi_d )
%公共量


c=3.0e+8;
lambda=0.25;
omega=2*pi*c/lambda;
r=100000000*lambda;
 
sum=0;
global nf;
global fr0;
global fn;
N=0;
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
if(nargin==4)
  r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)]; 
else
    r_d=r_s;
end
for ii=1:1:nf
    if (dot(r_s,fn(ii,:))/r>0.0001) &&  (dot(r_d,fn(ii,:))/r>0.0001) %被照射到的光子
    %handles.normalline=line([fr0(ii,1),fr0(ii,1)+fn(ii,1)]',[fr0(ii,2),fr0(ii,2)+fn(ii,2)]',[fr0(ii,3),fr0(ii,3)+fn(ii,3)]','color','r','linewidth',2);
    %if fr0(ii,3)>0.999
        ri_x=fr0(ii,1);
        ri_y=fr0(ii,2);
        ri_z=fr0(ii,3);%fr0 表示面元的中心坐标
        if(nargin==2)
            Delta_Ri=2*sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2));     
        else
            Delta_Ri=sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2))...
                      +sqrt(((ri_x-r.*sin(theta_d).*cos(phi_d)).^2+(ri_y-r.*sin(theta_d).*sin(phi_d)).^2+(ri_z-r.*cos(theta_d)).^2))
        end
        sum=sum+exp(1.00000i*omega*Delta_Ri/c);
        N=N+1;
    end    
end
output=(abs(sum)).^2/N;
fprintf('     N: %d \n',N);
end


%% 
function [ output ] = calc_3D_shelter_N_I_s( theta_s,phi_s,theta_d,phi_d)
%公共量


remark=shelter_prcess(theta_s,phi_s);
remark=shelter_prcess(theta_d,phi_d,remark);

global c;
global lambda;
omega=2*pi*c/lambda;
r=100000000*lambda;
 
sum=0;
global nf;
global fr0;
global fn;
N=0;
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
if(nargin==4)
  r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)]; 
else
    r_d=r_s;
end

for ii=1:1:nf
     if remark(ii)~=0
    %if fr0(ii,3)>0.999
        ri_x=fr0(ii,1);
        ri_y=fr0(ii,2);
        ri_z=fr0(ii,3);%fr0 表示面元的中心坐标
        if(nargin==2)
            Delta_Ri=2*sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2));     
        else
            Delta_Ri=sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2))...
                      +sqrt(((ri_x-r.*sin(theta_d).*cos(phi_d)).^2+(ri_y-r.*sin(theta_d).*sin(phi_d)).^2+(ri_z-r.*cos(theta_d)).^2));
        end
        sum=sum+exp(1.00000i*omega*Delta_Ri/c);
        N=N+1;
    end    
end
output=(abs(sum)).^2;
fprintf('shelter     N: %d  theta_s:%f  theta_d:%f  \n',N,theta_s/pi*180,theta_d/pi*180);
end


%% 
function [ output ] = calc_3D_N_I_s( theta_s,phi_s,theta_d,phi_d )
%公共量

global c;
global lambda;
omega=2*pi*c/lambda;
r=100000000*lambda;
 
sum=0;
global nf;
global fr0;
global fn;
N=0;
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
if(nargin==4)
  r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)]; 
else
    r_d=r_s;
end

for ii=1:1:nf
     if (dot(r_s,fn(ii,:))/r>0.0001) &&  (dot(r_d,fn(ii,:))/r>0.0001) %被照射到的光子
    %if fr0(ii,3)>0.999
        ri_x=fr0(ii,1);
        ri_y=fr0(ii,2);
        ri_z=fr0(ii,3);%fr0 表示面元的中心坐标
        if(nargin==2)
            Delta_Ri=2*sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2));     
        else
            Delta_Ri=sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2))...
                      +sqrt(((ri_x-r.*sin(theta_d).*cos(phi_d)).^2+(ri_y-r.*sin(theta_d).*sin(phi_d)).^2+(ri_z-r.*cos(theta_d)).^2));
        end
        sum=sum+exp(1.00000i*omega*Delta_Ri/c);
        N=N+1;
    end    
end
output=(abs(sum)).^2;
fprintf('no shelter    N: %d  theta_s:%f  theta_d:%f  \n',N,theta_s/pi*180,theta_d/pi*180);
end

%% 
function [ G ] = calc_G_3D(theta_s,phi_s)
%求解分母

d_theta=0.1;
d_phi=0.1; % if this value larger than 0.03  the result will be influented

G=0;
cnt=0;
for theta_d=0:d_theta:pi
    for phi_d=0:d_phi:2*pi
        cnt=cnt+1
        G=G+sin(theta_d)*calc_3D_N_I_s( theta_s,phi_s,theta_d,phi_d);  
    end
end
G=G*d_phi*d_theta;

end

%% 
function [ A ] = calc_A_3D(theta_s,phi_s)
%求解照射正交面积
global lambda;
r=100000000*lambda;
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
global fn;
global fA;
global nf;
sum_A_i=0;
for ii=1:1:nf
    if (dot(r_s,fn(ii,:))/r>0.0001)%被照射到的光子
        sum_A_i=sum_A_i+fA(ii)*dot(r_s,fn(ii,:))/r;
    end  
end
A=sum_A_i;

end

%% 遮挡点剔除
function [ new_remark ] = shelter_prcess(theta_s,phi_s,remark)
%remark 迭代值
global vertex;
global vertex_s;
global nf;
global fL_s;
global fr0_s;
global fn_s;
global fA_s;
global modle_max_length;

% if nargin==2
%     new_remark=zeros(1,nf);
% end
new_remark=zeros(1,nf);

modle_max_length=2;%投影图形的最大长度，决定Z-buff的大小
d_l=0.01;

m=fix(modle_max_length/2/d_l); %n正半轴或者负半轴的网格数目
buff_I=zeros(2*m+10,2*m+10);%存投影的编号
buff_Z=ones(2*m+10,2*m+10)-modle_max_length; %没有被照射到的深度设置为0


%旋转模型
z_j=[sin(theta_s)*cos(phi_s) sin(theta_s)*sin(phi_s) cos(theta_s)]';
x_j=[cos(theta_s)*cos(phi_s) cos(theta_s)*sin(phi_s) -sin(theta_s)]';
y_j=cross(z_j,x_j);
L=[x_j,y_j,z_j];
vertex_s(:,:,1)=vertex(:,:,1)*L';
vertex_s(:,:,2)=vertex(:,:,2)*L';%point (面元编号，x/y/z，第几个点) 
vertex_s(:,:,3)=vertex(:,:,3)*L';

%将模型平移便于编号
vertex_s(:,1,:)=vertex_s(:,1,:)+(m+1)*d_l; 
vertex_s(:,2,:)=vertex_s(:,2,:)+(m+1)*d_l;


fL_s=circshift(vertex_s,[0,0,-1])-vertex_s;%edge
fr0_s=mean(vertex_s,3);%center, reference point
frc_s=(circshift(vertex_s,[0,0,-1])+vertex_s)/2;% center point of edge
fn_s=cross(fL_s(:,:,1),fL_s(:,:,2),2);%cross product
fA_s=sqrt(dot(fn_s,fn_s,2))/2;%area of every facet
fn_s=bsxfun(@rdivide,fn_s,2*fA_s);%normal
 
%
for ii =1:1: nf
    if(fn_s(ii,3)>0.001) %面元法向量朝上
          for x_i=floor(min(vertex_s(ii,1,:))/d_l ):1:ceil(max(vertex_s(ii,1,:))/d_l )
              for y_i=floor(min(vertex_s(ii,2,:))/d_l ):1:ceil(max(vertex_s(ii,2,:))/d_l )
                  A=[vertex_s(ii,1,1),vertex_s(ii,2,1)];
                  B=[vertex_s(ii,1,2),vertex_s(ii,2,2)];
                  C=[vertex_s(ii,1,3),vertex_s(ii,2,3)];
                  P=[x_i*d_l,y_i*d_l];
                  if(is_in_Triangle(A,B,C,P)==1)
                      if(buff_Z(x_i,y_i)<fr0_s(ii,3))
                          buff_I(x_i,y_i)=ii;    
                          buff_Z(x_i,y_i)=fr0_s(ii,3);
                      end
                  end
              end
          end 
    end
end

[a,b]=size(buff_I);
for ii=1:a
    for jj=1:b
        if(buff_I(ii,jj)~=0)
            new_remark(buff_I(ii,jj))=1;       
        end
    end
end

if nargin==3
    new_remark=new_remark&remark;  
end

end
%% 判断一个点是否在三角形内
function [  result ] = is_in_Triangle(A,B,C,P)

v0 = C - A ;
v1 = B - A ;
v2 = P - A ;

dot00 = dot(v0,v0) ;
dot01 = dot(v0,v1) ;
dot02 = dot(v0,v2) ;
dot11 = dot(v1,v1) ;
dot12 = dot(v1,v2) ;

inverDeno = 1 / (dot00 * dot11 - dot01 * dot01) ;

u = (dot11 * dot02 - dot01 * dot12) * inverDeno ;
if (u < 0 || u > 1) %if u out of range, return directly
    result= 0 ;
    return;
end

v = (dot00 * dot12 - dot01 * dot02) * inverDeno ;
if (v < 0 || v > 1) % if v out of range, return directly
    result=0 ;
    return;
end
    
if u + v <= 1 
    result=1;
else
    result=0;
end

end


%% 数据读取与预处理
function data_process(route)

%cope for code.m
partfile=route;

%% 有几个变量其他函数计算的时候要用，声明为全局变量
global nf;
global fr0;
global fn;
global fA;
global c;
global lambda;
global vertex;
global mNumber;
c=3.0e+8;
lambda=0.25;

%% 
np=1;
ratio=1;
ha=[1;0;0];
ua=[0;0;1];
la=cross(ua,ha);
if abs(la)<.5
    errordlg('坐标轴设置错误');
    return;
end

vTmp=cell(1,np);
fTmp=cell(1,np);
nvTmp=zeros(1,np); %v代表节点
nfTmp=zeros(1,np); %f代表面元
for iPart=1:np
fid=fopen(partfile{iPart});%file id
    nvTmp(iPart)=0;%vertice number
    nfTmp(iPart)=0;%facets number
    while(~feof(fid))%end of file
        str=fscanf(fid,'%s',1);%scan string
        if(strcmp(str,'GRID*'))%if match GRID*
            nvTmp(iPart)=nvTmp(iPart)+1;%vertices number plus
            fscanf(fid,'%s',7);%scan useless string
        elseif(strcmp(str,'CTRIA3'))%if match CATIA Triangle,3 vertice
            nfTmp(iPart)=nfTmp(iPart)+1;%facets number plus
            fscanf(fid,'%s',5);%useless string
        end
    end
    vTmp{iPart}=zeros(nvTmp(iPart),3);
    fTmp{iPart}=zeros(nfTmp(iPart),3);
    ivT=0;
    ifT=0;
    frewind(fid);
    while(~feof(fid))%end of file
        str=fscanf(fid,'%s',1);%scan string
        if(strcmp(str,'GRID*'))%if match GRID*
            ivT=ivT+1;%vertices number plus
            fscanf(fid,'%s',1);%scan useless string
            vTmp{iPart}(ivT,[1,2])=fscanf(fid,'%f',2);%scan coordinate x,y
            fscanf(fid,'%s',3);%useless string
            vTmp{iPart}(ivT,3)=fscanf(fid,'%f',1);%scan coordinate z
        elseif(strcmp(str,'CTRIA3'))%if match CATIA Triangle,3 vertice
            ifT=ifT+1;%facets number plus
            fscanf(fid,'%s',2);%useless string
            fTmp{iPart}(ifT,:)=fscanf(fid,'%d',3);%indice of facet
        end
    end
    fclose(fid);%close file
end

vertex=vTmp{1};
facet=fTmp{1}; %每个面元的节点号
mNumber=ones(sum(nfTmp),1);
for iPart=2:np%deliminate same vertex
    fTmp{iPart}=fTmp{iPart}+size(vertex,1);
    for p2=nvTmp(iPart):-1:1
        tmp=bsxfun(@minus,vTmp{iPart}(p2,:),vertex);
        p1=find(dot(tmp,tmp,2)<myeps,1);
        if p1
            fTmp{iPart}(fTmp{iPart}==p2+size(vertex,1))=p1;
            fTmp{iPart}(fTmp{iPart}>p2+size(vertex,1))=fTmp{iPart}(fTmp{iPart}>p2+size(vertex,1))-1;
            vTmp{iPart}(p2,:)=[];
            nvTmp(iPart)=nvTmp(iPart)-1;
        end
    end
    vertex=[vertex;vTmp{iPart}];
    facet=[facet;fTmp{iPart}];
    mNumber(sum(nfTmp(1:iPart-1))+1:sum(nfTmp(1:iPart)))=iPart;
end

vertex=vertex*[ha,la,ua]*ratio; 
nf=size(facet,1);
nv=size(vertex,1);
vertex=cat(3,vertex(facet(:,1),:),vertex(facet(:,2),:),vertex(facet(:,3),:));%point (面元编号，x/y/z，第几个点) 
clear vTmp;
clear fTmp;
clear nXTmp;
clear nFTmp;

%% preprocess for PO
fL=circshift(vertex,[0,0,-1])-vertex;%edge
fr0=mean(vertex,3);%center, reference point
frc=(circshift(vertex,[0,0,-1])+vertex)/2;% center point of edge
fn=cross(fL(:,:,1),fL(:,:,2),2);%cross product
fA=sqrt(dot(fn,fn,2))/2;%area of every facet
fn=bsxfun(@rdivide,fn,2*fA);%normal
%axes(handles.axes_mesh);
fill3(permute(vertex(:,1,:),[3,1,2]),permute(vertex(:,2,:),[3,1,2]),permute(vertex(:,3,:),[3,1,2]),mNumber');
%handles.normalline=line([fr0(:,1),fr0(:,1)+fn(:,1)*0.1]',[fr0(:,2),fr0(:,2)+fn(:,2)*0.1]',[fr0(:,3),fr0(:,3)+fn(:,3)*0.1]','color','r','linewidth',1);
axis equal tight;
xlabel('x');
ylabel('y');
zlabel('z');


end


function [  ] = draw(remark)
global fL_s;
global fr0_s;
global fn_s;
global fA_s;
global nf;

global vertex_s;
global mNumber;
fill3(permute(vertex_s(:,1,:),[3,1,2]),permute(vertex_s(:,2,:),[3,1,2]),permute(vertex_s(:,3,:),[3,1,2]),mNumber');
xlabel('x');
ylabel('y');
zlabel('z');
for ii=1:nf
    if remark(ii)~=0
        line([fr0_s(ii,1),fr0_s(ii,1)+fn_s(ii,1)*sqrt(fA_s(ii))]',[fr0_s(ii,2),fr0_s(ii,2)+fn_s(ii,2)*sqrt(fA_s(ii))]',...
        [fr0_s(ii,3),fr0_s(ii,3)+fr0_s(ii,3).*sqrt(fA_s(ii))]','color','r','linewidth',2);
    end
end

end





