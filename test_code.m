
function test_code
test008();
end


%%  二面角数值解
function test008

if ~exist('dihedral.mat')
    l_theta=0:0.1:pi/2;
    l_sigma_Q_1=zeros(size(l_theta));

    ii=1;
    for theta= l_theta
        A=1*sqrt(2)*cos(theta-pi/4);
        G=calc_G(theta,0,'dihedral');
        l_sigma_Q_1(ii)=20*log10(4*pi*A*calc_N_I_s(theta,0,'dihedral')/G);
        fprintf('%d\\%d   G:%f\n',ii,length(l_theta),1);
        ii=ii+1;
    end
else
   load('dihedral.mat') 
end

plot(l_theta,l_sigma_Q_1,'-*');
% ylim([-5, 5]);
grid on
save('dihedral','l_theta','l_sigma_Q_1');
end


%%  理论解和数值解图像
function test003
test004();
hold on
test005('circle');
legend('数值解','理论解','Location','northeast')
hold on
% ylim([-5, 5]);
grid on
end

%% 求均值 (球体)
function test007
   load('sphere-r-0.5-grid-0.025.mat')
% load('sphere-r-1-grid-0.1.mat')
%  load('sphere-r-0.5-grid-0.01.mat')
l_sigma_Q_2=zeros(size(l_theta));



ii=1;
for theta=l_theta
     sigma_Q_1(ii)=10^(l_sigma_Q_1(ii)/20);
     sigma_Q_2(ii)=pi*0.25;
      ii=ii+1; 
end

ii=1;
for theta=l_theta
    if isnan(sigma_Q_1(ii)) || isnan(sigma_Q_2(ii))
        sigma_Q_1(ii)=[];
        sigma_Q_2(ii)=[];
        l_theta(ii)=[]; 
    else
        ii=ii+1; 
    end    
end


 plot(l_theta,20*log10(sigma_Q_1),'-*',l_theta,20*log10(sigma_Q_2),'-*');
 value_mean1=mean(sigma_Q_1);
 fprintf('数值解均值:%f  = %f dB  \n',value_mean1,  20*log10(abs(value_mean1))  );
 value_mean2=mean(sigma_Q_2);
 fprintf('理论解解均值:%f  = %f dB  \n',value_mean2,  20*log10(abs(value_mean2))  );
 fprintf('均值差：%f dB \n',20*log10(abs(value_mean1))-20*log10(abs(value_mean2)));

end

%% 求均值
function test006
%  load('cube-l-1-w-1-grid-0.05.mat')
% load('Cylinder-r-1-l-0.1-grid-0.03.mat')
load('TriangularPrism-d-1-h-1-l-0.2-grid-0.025.mat')

sigma_Q_1=zeros(size(l_theta));
sigma_Q_2=zeros(size(l_theta));

ii=1;
for theta=l_theta
     sigma_Q_1(ii)=10^(l_sigma_Q_1(ii)/20);
     sigma_Q_2(ii)=calc_sigma_Q(theta,'triangle');
      ii=ii+1; 
end

ii=1;
for theta=l_theta
    if isnan(sigma_Q_1(ii)) || isnan(sigma_Q_2(ii))
        sigma_Q_1(ii)=[];
        sigma_Q_2(ii)=[];
        l_theta(ii)=[]; 
    else
        ii=ii+1; 
    end    
end
 plot(l_theta,20*log10(sigma_Q_1),'-*',l_theta,20*log10(sigma_Q_2),'-*');
 value_mean1=mean(sigma_Q_1);
 fprintf('数值解均值:%f  = %f dB  \n',value_mean1,  20*log10(abs(value_mean1))  );
 
 value_mean2=mean(sigma_Q_2);
 fprintf('理论解解均值:%f  = %f dB  \n',value_mean2,  20*log10(abs(value_mean2))  );
 fprintf('均值差：%f dB \n',20*log10(abs(value_mean1))-20*log10(abs(value_mean2)));

end




%% 绘制理论解图像
function test005(shape)
l_theta=-pi/2+0.1:0.01:pi/2-0.1;
l_sigma_Q_1=zeros(size(l_theta));

ii=0;
for theta= l_theta
    ii=ii+1;
    l_sigma_Q_1(ii)=20*log10(calc_sigma_Q(theta,shape));
end
plot(l_theta,l_sigma_Q_1,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)');  
legend('理论解','Location','northeast')
end

%% 绘制 数值解
function test004
% load('cube-l-1-w-1-grid-0.05.mat');
% load('Cylinder-r-1-l-0.5-grid-0.1.mat');
% load('sphere-r-0.5-grid-0.025.mat');
% load('TriangularPrism-d-1-h-1-l-0.2-grid-0.1.mat');
% load('cube-l-1-w-1-grid-0.05.mat');
% load('sphere-r-0.5-grid-0.025.mat');
%  load('TriangularPrism-d-1-h-1-l-0.2-grid-0.025.mat')
load('Cylinder-r-1-l-0.1-grid-0.03.mat')
plot(l_theta,l_sigma_Q_1,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)'); 
grid on;
%save test_code;
end



%% 球体理论解图像
function test001
l_theta=0:0.01:pi/2;
l_sigma_Q_2=zeros(size(l_theta));

ii=0;
for theta= l_theta
    ii=ii+1;
    l_sigma_Q_2(ii)=20*log10(pi*0.5*0.5);
end
plot(l_theta,l_sigma_Q_2,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)');  
%ylim([-1, 2]);
end

