%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence  
%
%  Author: wangwei
%  Data  : 6/15/2017
%  refer in:
%  "Microphone Array Post-Filter Based on Noise Field Coherence" IEEE 2003
%  "MICROPHONE ARRAY POST-FILTER FOR DIFFUSE NOISE FIELD"  IEEE 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all;
c = 340; % speed of sound


%% load recorded office noise audio
noisepath = '../sound/rec1/';
[noise1,fs] = audioread([noisepath,'音轨-2.wav']);
noise2 = audioread([noisepath,'音轨-3.wav']);
noise3 = audioread([noisepath,'音轨-4.wav']);
noise4 = audioread([noisepath,'音轨-5.wav']);
noise = [noise1,noise2,noise3,noise4];

%% 

pathname = '../sound/';
N = size(noise,2);%通道数
%use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'speech.wav']);
close all
%% use RIR to generate reverberation output
source_pos = [1.94 3.5 2];              % Source position [x y z] (m)

% [s,fs] = audioread('an103-mtms-senn4.wav');

scale = 10;

h = RIR_generator( source_pos,0.2);
h = h*scale;

x1 = conv(speech,h(1,:));
x2 = conv(speech,h(2,:));
x3 = conv(speech,h(3,:));
x4 = conv(speech,h(4,:));
x = [x1,x2,x3,x4];


% noise = noise(1:length(x(:,1)),:);
noise0 = noise(:,1);
%scale source signal to obtain 0 dB input SNR    
% speech = speech(1:length(noise1))/2;   

        

%% 模拟阵列接收信号,添加噪声

fs = 16000;

% use a clean audio plus a real office noise to simulate as a perfect array
% time-aligned output
noise = noise(1:length(x(:,1)),:);
x = x+noise;

% DS output
DS = sum(x,2)/N;  


% compute input SNR
speech = x1; % desired signal
noiseDS = sum(noise,2);
SNR_input = 10*log10(sum(speech.^2)/sum(noiseDS.^2))


%% 
L = 512;                   %16KHz,32 ms

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
znoise = zeros(1,length(x(:,1)));
zspeech = zeros(1,length(x(:,1)));
t = 1;
%%
tic

%the distance between two adjacent microphone,XMOS microphone array board
%was used
r = 0.032; 

% the distance between two microphone

dij = [r*sqrt(2),2*r,r*sqrt(2),r*sqrt(2),r*2,r*sqrt(2)];
f = 0:fs/256:fs/2;
c = 343;
    
Inc = 128;   
k_optimal = 2;
hwt = waitbar(0,'general poster filter');
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-L
    for i = 1:N
        % estimate auto  power spectral density
        Pxii(i,:) = cpsd(x(p:p+L-1,i),x(p:p+L-1,i),hanning(128),64); 
    end
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            % estimates the cross powerspectral density,
            % Welch's average method,the window width can make a trade-off between variance and
            % resolution
            Pxij(t,:) = cpsd(x(p:p+L-1,i),x(p:p+L-1,j),hanning(128),64);
            
            %diffuse noise coherence function,can be modeled as sin(x)/(x) 
            T = sin(k_optimal*pi*f*dij(t)*1/c)./(k_optimal*pi*f*dij(t)*1/c);T(1) = 0.998;%T(2) = 0.996;
           
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
    
    % handle the indeterminite soulution when MSC≈1
    Pss(Pss<0) = 0;
    
    % calculate the frequency domain filter coefficient
    W_e = real(Pss)./Pssnn;
   
    % interpolate filter to match the length of the signal
    W_e_i = interp1(1:length(W_e),W_e,1:.5:length(W_e));

    % get the final filter coefficient,matlab CPSD function we used above only get the
    % oneside PSD,but fft(signal) is twoside,so we add the conjugate side
    W = [W_e_i,conj(fliplr(W_e_i(2:256)))];
    % transfor the signal to frequency domain
    Xds = fft([DS(p:p+L-1)'].*window');
    
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:L);
    
    % keep the signal
    z(p:p+L-1) = z(p:p+L-1) + s_est;
    
    
    % perform the same operation to noise and desired signal in order to calculate overall
    % SNR_output,if on need to observe SNR,comment this section
    XnoiseDS = fft(noiseDS(p:p+L-1)'.*window');
    noiseDS_filtered = W.*(XnoiseDS);
    iNoise = ifft(noiseDS_filtered);
    znoise(p:p+L-1) = znoise(p:p+L-1) + iNoise(1:L);
    
    Xspeech = fft(speech(p:p+L-1)'.*window');
    speech_filtered = W.*(Xspeech);
    iSpeech = ifft(speech_filtered);
    zspeech(p:p+L-1) = zspeech(p:p+L-1) + iSpeech(1:L);

    waitbar(p/(length(x(:,1))));
end
close(hwt);

% compute output SNR
SNR_output = 10*log10(sum(zspeech.^2)/sum(znoise.^2))
toc


