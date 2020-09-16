
addpath('mexcgns3_1');
startup_mexcgns

[xs, elems, typestr, var_nodes] = readcgns('dom-4.cgns');

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

facet=elems;
vertex =xs;

vertex=vertex*[ha,la,ua]*ratio; 
nf=size(facet,1);
nv=size(vertex,1);
vertex=cat(3,vertex(facet(:,1),:),vertex(facet(:,2),:),vertex(facet(:,3),:));%point (面元编号，x/y/z，第几个点) 


%% preprocess for PO
fL=circshift(vertex,[0,0,-1])-vertex;%edge
fr0=mean(vertex,3);%center, reference point
frc=(circshift(vertex,[0,0,-1])+vertex)/2;% center point of edge
fn=cross(fL(:,:,1),fL(:,:,2),2);%cross product
fA=sqrt(dot(fn,fn,2))/2;%area of every facet
fn=bsxfun(@rdivide,fn,2*fA);%normal
%axes(handles.axes_mesh);
fill3(permute(vertex(:,1,:),[3,1,2]),permute(vertex(:,2,:),[3,1,2]),permute(vertex(:,3,:),[3,1,2]),3');
%normalline=line([fr0(:,1),fr0(:,1)+fn(:,1).*sqrt(fA)]',[fr0(:,2),fr0(:,2)+fn(:,2).*sqrt(fA)]',[fr0(:,3),fr0(:,3)+fn(:,3).*sqrt(fA)]','color','r','linewidth',2);
axis equal tight;
xlabel('x');
ylabel('y');
zlabel('z');



