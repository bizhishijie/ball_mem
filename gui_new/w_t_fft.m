load('w_t.mat');
load('p_t.mat');

fs1=1000;
T1=1/fs1;
L1=length(z_1);
t1=(0:L1-1)*T1;
N1=2^nextpow2(L1);
y1=fft(z_1,N1)/L1;
f1=fs1/2*linspace(0,1,N1/2);
figure(1)
plot(f1,2*abs(y1(1:N1/2)));

fs2=1000;
T2=1/fs2;
L2=length(z_2);
t2=(0:L2-1)*T2;
N2=2^nextpow2(L2);
y2=fft(z_2,N2)/L2;
f2=fs2/2*linspace(0,1,N2/2);
figure(2)
plot(f2,2*abs(y2(1:N2/2)));
