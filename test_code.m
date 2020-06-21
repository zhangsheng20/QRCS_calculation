
function test_code
test003();

test005('triangle');
hold on;
test004();
legend('理论解','数值解','Location','northeast')
grid on
xlim([-1.2,1.2])
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
figure(2);
plot(l_theta,l_sigma_Q_1,'-*')
xlabel('\theta(rad)');
ylabel('\sigma_Q(dB/m^2)');  
legend('理论解','Location','northeast')
end

%% 绘制圆盘数值解
function test004
%load('cube_1m_1m_0.03m_.mat')
% load('Cylinder-r-1-l-0.5-grid-0.1.mat')
load('sphere-r-0.5-grid-0.025.mat')
% load('TriangularPrism-d-1-h-1-l-0.2-grid-0.1.mat')
% load('cube-l-1-w-1-grid-0.05.mat')
%load('sphere-r-0.5-grid-0.025.mat')
plot(l_theta,l_sigma_Q_1,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)');  
%save test_code;
end


%% 球体理论解和数值解图像
function test003
test004();
hold on
test001();
legend('数值解','理论解','Location','northeast')
hold on
ylim([-5, 5]);
grid on
end

%% 球体理论解图像
function test001
l_theta=0:0.01:pi/2;
l_sigma_Q_2=zeros(size(l_theta));

ii=0;
for theta= l_theta
    ii=ii+1;
    l_sigma_Q_2(ii)=pi*0.5*0.5;
end
plot(l_theta,l_sigma_Q_2,'-*')
xlabel('\theta(rad)');
ylabel('QRCS(dB/m^2)');  
%ylim([-1, 2]);
end

