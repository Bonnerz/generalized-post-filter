%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence  
%
%  Author: wangwei
%  Data  : 6/15/2017
%  env   : test in matlab2012b
%  refer in:
%  "Microphone Array Post-Filter Based on Noise Field Coherence" IEEE 2003
%  "MICROPHONE ARRAY POST-FILTER FOR DIFFUSE NOISE FIELD"  IEEE 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all;
c = 340; % speed of sound
N = 7;

%% load recorded office noise audio
noisepath = '../sound/noise/';
[noise3,fs] = audioread([noisepath,'“ÙπÏ-4.wav']);
noise0 = audioread([noisepath,'“ÙπÏ.wav']);
noise6 = audioread([noisepath,'“ÙπÏ-7.wav']);
noise4 = audioread([noisepath,'“ÙπÏ-5.wav']);
noise1 = audioread([noisepath,'“ÙπÏ-2.wav']);
noise5 = audioread([noisepath,'“ÙπÏ-6.wav']);
noise2 = audioread([noisepath,'“ÙπÏ-3.wav']);
noise = [noise0,noise1,noise2,noise3,noise4,noise5,noise6];

%% 

pathname = '../sound/';
pathname = '../sound/num2_MIC5/';
% pathname = '../sound/meetingroom_MIC3/';
pathname = '../sound/“Ù∆µ≤…—˘Ωµ‘Î0417/01/';
[wav3,fs] = audioread([pathname,'“ÙπÏ-4.wav']);
wav0 = audioread([pathname,'“ÙπÏ.wav']);
wav6 = audioread([pathname,'“ÙπÏ-7.wav']);
wav4 = audioread([pathname,'“ÙπÏ-5.wav']);
wav1 = audioread([pathname,'“ÙπÏ-2.wav']);
wav5 = audioread([pathname,'“ÙπÏ-6.wav']);
wav2 = audioread([pathname,'“ÙπÏ-3.wav']);

%% 
% w0 = resample(wav0,10,1);[d1, cv, pro] = delayesttm(w1,w0,fs);delay(1) = (d1*fs)
% w1 = resample(wav1,10,1);
% w2 = resample(wav2,10,1);
% w3 = resample(wav3,10,1);
% w4 = resample(wav4,10,1);
% w5 = resample(wav5,10,1);
% w6 = resample(wav6,10,1);


delay = zeros(N,1);
anchor = wav1;
s = struct('interpord',0);
[d1, cv, pro] = delayesttm(anchor,wav0,fs,s);delay(1) = (d1*fs);
[d2, cv, pro] = delayesttm(anchor,wav1,fs,s);delay(2) = (d2*fs);
[d3, cv, pro] = delayesttm(anchor,wav2,fs,s);delay(3) = (d3*fs);
[d4, cv, pro] = delayesttm(anchor,wav3,fs,s);delay(4) = (d4*fs);
[d5, cv, pro] = delayesttm(anchor,wav4,fs,s);delay(5) = (d5*fs);
[d6, cv, pro] = delayesttm(anchor,wav5,fs,s);delay(6) = (d6*fs);
[d7, cv, pro] = delayesttm(anchor,wav6,fs,s);delay(7) = (d7*fs);
delay
d = -1*delay;
md = max(d);
x = [wav0(d(1)+1:end-md+d(1)),wav1(d(2)+1:end-md+d(2)),wav2(d(3)+1:end-md+d(3)),wav3(d(4)+1:end-md+d(4)),wav4(d(5)+1:end-md+d(5)),wav5(d(6)+1:end-md+d(6)),wav6(d(7)+1:end-md+d(7))];
DS = sum(x,2)/N;
% audiowrite('wav/msc/meetingroom_MIC3/TestMIC1_DS.wav',DS,fs);
%% 
L = 1024;                   %16KHz,32 ms

% window = hamming(L);
window = triang(L);
% window = ones(L,1);
% window = repmat(ham,1,N)';

P_len = 129;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);

Pss = zeros(1,P_len);

Pss_e = zeros((N*N-N)/2,P_len);

z = zeros(1,length(x(:,1)));

t = 1;
%%
tic

%the distance between two adjacent microphone,XMOS microphone array board
%was used
r = 0.045; 
r = 0.042;
% the distance between two microphone
dij = [r,r,r,r,r,r,...
       r,r*sqrt(3),r*2,r*sqrt(3),r,...
       r,r*sqrt(3),r*2,r*sqrt(3),...
       r,r*sqrt(3),r*2,...
       r,r*sqrt(3),...
       r];
f = 0:fs/256:fs/2;
c = 346;    %velocity of sound wave spread in 25 celsius degree air
    
Inc = L/2;   
hwt = waitbar(0,'general poster filter');
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-L*2
    for i = 1:N
        % estimate auto  power spectral density
        Pxii(i,:) = cpsd(x(p:p+L-1,i),x(p:p+L-1,i),hanning(256),128); 
    end
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            % estimates the cross powerspectral density,
            % Welch's average method,the window width can make a trade-off between variance and
            % resolution
            Pxij(t,:) = cpsd(x(p:p+L-1,i),x(p:p+L-1,j),hanning(256),128);
            
            %diffuse noise coherence function,can be modeled as sin(x)/(x) 
            T = sin(2*pi*f*dij(t)*2/c)./(2*pi*f*dij(t)*2/c);T(1) = 0.998;%T(2) = 0.996;
           
            % estimate source signal's PSD
            Pss_e(t,:) = (real(Pxij(t,:)) - 0.5*real(T).*(Pxii(i,:)+Pxii(j,:)))...
                         ./...
                         (ones(1,P_len) - real(T));
             t = t+1;
        end
    end
    t = 1;
    % take the average of multichanel signal to improve robustness
    Pss = sum(Pss_e)*2/(N*N-N); 
    
    % handle the indeterminite soulution when MSC°÷1
    Pss(Pss<0) = 0;
    
    % calculate the frequency domain filter coefficient
    W_e = real(Pss)./Pssnn;
   
    % interpolate filter to match the length of the signal
    W_e_i = interp1(1:length(W_e),W_e,1:0.25:length(W_e));

    % get the final filter coefficient,matlab CPSD function we used above only get the
    % oneside PSD,but fft(signal) is twoside,so we add the conjugate side
    W = [W_e_i,conj(fliplr(W_e_i(2:512)))];
    % transfor the signal to frequency domain
    Xds = fft([DS(p:p+L-1)'].*window');
    
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:L);
    
    % keep the signal
    z(p:p+L-1) = z(p:p+L-1) + s_est;
    waitbar(p/(length(x(:,1))));
end
close(hwt);
toc
% audiowrite('wav/msc/meetingroom_MIC3/TestMIC1_inc512_dijt042X2_0998_hanning256128_Pss0.wav',z,fs);

