close all
noisepath = 'sound/noise/';
[noise3,fs] = audioread([noisepath,'“ÙπÏ-4.wav']);
noise0 = audioread([noisepath,'“ÙπÏ.wav']);
noise6 = audioread([noisepath,'“ÙπÏ-7.wav']);
noise4 = audioread([noisepath,'“ÙπÏ-5.wav']);
noise1 = audioread([noisepath,'“ÙπÏ-2.wav']);
noise5 = audioread([noisepath,'“ÙπÏ-6.wav']);
noise2 = audioread([noisepath,'“ÙπÏ-3.wav']);
noise = [noise0,noise1,noise2,noise3,noise4,noise5,noise6];

P11 = cpsd(noise0,noise0,hanning(256),128,256);
P22 = cpsd(noise1,noise1,hanning(256),128,256);
P12 = cpsd(noise0,noise1,hanning(256),128,256);

d = 0.042;
k = 2.8;
% coherence
msc12 = P12./sqrt((P11.*P22));

% find optimal k
E0 = 100;
k_optimal = 0;
for k = 0.5:0.05:3.5
    T = sin(2*pi*f*d*k/c)./(2*pi*f*d*k/c);T(1) = 0.998;
    E = sum((real(msc12)'-T).^2);
    if(E<E0)
        k_optimal = k;
    end
    E0 = E;
end

T = sin(2*pi*f*d*k_optimal/c)./(2*pi*f*d*k_optimal/c);T(1) = 0.998;

% error
E = sum((real(msc12)'-T).^2)

figure,plot(real(msc12))
hold on,plot(T)