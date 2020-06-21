


function calculate_3D
%cope for code.m
%filename='.\cube_1m_1m_0.03m_.dat';
%filename='Cylinder-r-0.5-l-0.1-grid-0.025.dat';
%  filename='sphere-r-0.5-grid-0.01.dat';
% filename='sphere-r-1-grid-0.1.dat';
% filename='cube-l-1-w-1-grid-0.03.dat';
%filename='Cylinder-r-0.3-l-0.5-grid-0.03.dat';
% filename='sphere-r-0.5-grid-0.025.dat'
% filename='TriangularPrism-d-1-h-1-l-0.2-grid-0.025.dat';
filename='Cylinder-r-1-l-0.1-grid-0.03.dat';
partfile={filename};

%% 有几个变量其他函数计算的时候要用，声明为全局变量
global nf;
global fr0;
global fn;
global fA;
global c;
global lambda;
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
%handles.normalline=line([fr0(:,1),fr0(:,1)+fn(:,1).*sqrt(fA)]',[fr0(:,2),fr0(:,2)+fn(:,2).*sqrt(fA)]',[fr0(:,3),fr0(:,3)+fn(:,3).*sqrt(fA)]','color','r','linewidth',2);
axis equal tight;
xlabel('x');
ylabel('y');
zlabel('z');

%% draw picture

l_theta=-pi/2:0.01:pi/2;
l_sigma_Q_1=zeros(size(l_theta));

ii=0;

for theta= l_theta
    ii=ii+1;
    A=calc_A_3D(theta,0);

    G=calc_G_3D(theta,0);
    fprintf('progress: %d \\ %d  A:%f   G:%f  G/A^2:%f  G/A:%f  \n', ...
                      ii,length(l_theta),A,G,G/A^2,G/A);
    N_I_s= calc_3D_N_I_s_2(theta,0,theta,0);        
    l_sigma_Q_1(ii)=20*log10(4*pi*A*N_I_s/G);
    fprintf('N_I_s:%f  \n',N_I_s)
    %l_sigma_Q_1(ii)=20*log10(calc_3D_N_I_s_2(theta,0,theta,0));
end
figure(2);
plot(l_theta,l_sigma_Q_1,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)');  
grid on
filename=strrep(filename,'.dat','.mat');
save (filename,'l_theta','l_sigma_Q_1');


end


%%
function [ output ] = calc_3D_N_I_s_2( theta_s,phi_s,theta_d,phi_d )
%公共量

global c;
global lambda;
omega=2*pi*c/lambda;
r=100000*lambda;
 
global fr0;
global fn;
global fA;

r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
if(nargin==4)
  r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)]; 
else
    r_d=r_s;
end

if(nargin==2)
    Mat_Delta_Ri=2*sqrt(((fr0(:,1)-r*sin(theta_s).*cos(phi_s)).^2+(fr0(:,2)-r*sin(theta_s).*sin(phi_s)).^2+(fr0(:,3)-r*cos(theta_s)).^2));     
else
    Mat_Delta_Ri=sqrt(((fr0(:,1)-r*sin(theta_s).*cos(phi_s)).^2+(fr0(:,2)-r*sin(theta_s).*sin(phi_s)).^2+(fr0(:,3)-r*cos(theta_s)).^2))...
              +sqrt(((fr0(:,1)-r*sin(theta_d).*cos(phi_d)).^2+(fr0(:,2)-r*sin(theta_d).*sin(phi_d)).^2+(fr0(:,3)-r*cos(theta_d)).^2));
end

a=fn*r_s'/r;
b=fn*r_d'/r;
mask=(a>0.08) &  (b>0.08) & fn(:,3)>0.5;
d=exp(1.00000i*omega*Mat_Delta_Ri/c).*fA.*mask;
allsum=sum(d);
output=(abs(allsum)).^2;

% N=sum(mask);
% fprintf('     N: %d \n',N);
end




%% 
function [ G ] = calc_G_3D(theta_s,phi_s)
%求解分母

d_theta=0.03;
d_phi=0.2; % if this value larger than 0.03  the result will be influented

G=0;
cnt=0;
for theta_d=0:d_theta:pi
    for phi_d=0:d_phi:2*pi
        cnt=cnt+1;
        G=G+sin(theta_d)*calc_3D_N_I_s_2( theta_s,phi_s,theta_d,phi_d);  
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

a=fn*r_s';
mask=a>0;
Mat_sum=mask.*fA.*a/r;
A=sum(Mat_sum);


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
    if (dot(r_s,fn(ii,:))/r>0) &&  (dot(r_d,fn(ii,:))/r>0) %被照射到的光子
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
N
end

%% 
function [ output ] = calc_3D_N_I_s( theta_s,phi_s,theta_d,phi_d )
%公共量

global c;
global lambda;
omega=2*pi*c/lambda;
r=100000*lambda;
 
sum=0;
global nf;
global fr0;
global fn;
global fA;
N=0;
r_s=[r*sin(theta_s)*cos(phi_s),r*sin(theta_s)*sin(phi_s),r*cos(theta_s)];
if(nargin==4)
  r_d=[r*sin(theta_d)*cos(phi_d),r*sin(theta_d)*sin(phi_d),r*cos(theta_d)]; 
else
    r_d=r_s;
end

for ii=1:1:nf 
 
  a=dot(r_s,fn(ii,:))/r;
  b=dot(r_d,fn(ii,:))/r;
  %fprintf(' %d:a=%f,b=%f \n',ii,a,b); 
   if (a>0.1) &&  (b>0.1)  && fn(ii,3)>0.5         %被照射到的原子
     %line([fr0(ii,1),fr0(ii,1)+fn(ii,1)*0.1]',[fr0(ii,2),fr0(ii,2)+fn(ii,2)*0.1]',[fr0(ii,3),fr0(ii,3)+fn(ii,3)*0.1]','color','r','linewidth',1);  
        ri_x=fr0(ii,1);
        ri_y=fr0(ii,2);
        ri_z=fr0(ii,3);%fr0 表示面元的中心坐标
        if(nargin==2)
            Delta_Ri=2*sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2));     
        else
            Delta_Ri=sqrt(((ri_x-r.*sin(theta_s).*cos(phi_s)).^2+(ri_y-r.*sin(theta_s).*sin(phi_s)).^2+(ri_z-r.*cos(theta_s)).^2))...
                      +sqrt(((ri_x-r.*sin(theta_d).*cos(phi_d)).^2+(ri_y-r.*sin(theta_d).*sin(phi_d)).^2+(ri_z-r.*cos(theta_d)).^2));
        end
        sum=sum+exp(1.00000i*omega*Delta_Ri/c)*fA(ii);
        N=N+1;
    end    
end
output=(abs(sum)).^2;
fprintf('     N: %d \n',N);
end








