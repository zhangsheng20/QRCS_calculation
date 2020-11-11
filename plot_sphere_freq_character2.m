
clear;

R=0.1;
c=3.0e+8;

ii=0;
l_ka=0:1:20;
l_NIs=zeros(size(l_ka));
l_sigama_Q=zeros(size(l_ka));
frequency_G=zeros(size(l_ka));
for ka=l_ka
    ii=ii+1
    k=ka/R;
    lambda=2*pi/k;
    frequency_G(ii)=c/lambda/10^9;
    l_sigama_Q(ii) = calculate_sphere_G2(k,R);
end

plot(l_ka,10*log10(l_sigama_Q),'-')

 