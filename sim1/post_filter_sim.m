%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence  
%
%  Author: wangwei
%  Data  : 6/15/2017
%  refer to:
%  "Microphone Array Post-Filter Based on Noise Field Coherence" IEEE 2003
%  "MICROPHONE ARRAY POST-FILTER FOR DIFFUSE NOISE FIELD"  IEEE 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all;
c = 340; % speed of sound


%% generate diffuse noise
[ noise,Pos ] = noise_gen;

noise = noise'/20;
N = size(noise,2);        %Í¨µÀÊý
%% 

pathname = '../sound/';

%use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'an101-mtms-senn3.wav']);
%scale source signal to obtain 0 dB input SNR     
noise = noise(1:length(speech),:);
% [speech ,fs] = audioread([pathname,'1KHz_16000.wav']);
s = repmat(speech,1,N);

pathname = '../sound/speech/';
               

%% 

fs = 16000;

% use a clean audio plus a real office noise to simulate as a perfect array
% time-aligned output
x = s+noise;

% DS output
% perfect time alignment
DS = sum(x,2)/N;  

% compute input SNR
noiseDS = sum(noise,2);
SNR_input = 10*log10(sum(speech.^2)/sum(noiseDS.^2))


%% 
L = 2048;                   %16KHz,32 ms
N_FFT = 256;

% window = hamming(L);
window = hamming(N_FFT);
% window = ones(L,1);
% window = repmat(ham,1,N)';

P_len = N_FFT/2+1;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);
Pxii_pre = ones(N,P_len);
Pssnn_pre = ones(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);
Pxij_pre = ones((N*N-N)/2,P_len);

Pss = zeros(1,P_len);

Pss_e = zeros((N*N-N)/2,P_len);

z = zeros(1,length(x(:,1)));
znoise = zeros(1,length(x(:,1)));
zspeech = zeros(1,length(x(:,1)));
t = 1;
%%
tic

%respeaker 4MIC array V2.0
r = 0.032; 

% the distance between two microphone
dij = [r*sqrt(2),2*r,r*sqrt(2),r*sqrt(2),r*2,r*sqrt(2)];

f = 0:fs/256:fs/2;
c = 343;

alpha = 0.95;

w = 2*pi*fs*(0:N_FFT/2)/N_FFT;
Inc = 128; 
k_optimal = 1;
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
%         Xi = fft(x(p:p+N_FFT-1,i).*window);
%         Xi2 = abs(Xi).^2;
%         Pxii_curr = Xi2;
%         Pxii(i,:) = alpha*Pxii_pre(i,:)+(1-alpha)*Pxii_curr(1:N_FFT/2+1).';      
    end
    Pxii_pre = Pxii;
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            % estimates the cross powerspectral density,
            % Welch's average method,the window width can make a trade-off between variance and
            % resolution
%             Xi = fft(x(p:p+N_FFT-1,i).*window).';
%             Xj = fft(x(p:p+N_FFT-1,j).*window).';
%             Pxij(t,:) = alpha*Pxij_pre(t,:)+(1-alpha)*Xi.*conj(Xj);
            
            Pxij(t,:) = cpsd(x(p:p+L-1,i),x(p:p+L-1,j),hanning(128),64);
            
            %diffuse noise coherence function,can be modeled as sin(x)/(x)
            T = sin(2*pi*f*dij(t)*1/c)./(2*pi*f*dij(t)*1/c);T(1) = 0.998;%T(2) = 0.996;
           
            % estimate source signal's PSD
            Pss_e(t,:) = (real(Pxij(t,:)) - 0.5*real(T).*(Pxii(i,:)+Pxii(j,:)))...
                         ./...
                         (ones(1,P_len) - real(T));
             t = t+1;
        end
    end
    Pxij_pre = Pxij;
    t = 1;
    % take the average of multichanel signal to improve robustness
    Pss = sum(Pss_e)*2/(N*N-N); 
    
    % handle the indeterminite soulution when MSC¡Ö1
    Pss(Pss<0) = 0;
    
    % calculate the frequency domain filter coefficient
    W_e = real(Pss)./Pssnn;
   
    % interpolate filter to match the length of the signal
%     W_e_i = interp1(1:length(W_e),W_e,1:.5:length(W_e));

    % get the final filter coefficient,matlab CPSD function we used above only get the
    % oneside PSD,but fft(signal) is twoside,so we add the conjugate side
    W = [W_e,conj(fliplr(W_e(2:128)))];
    % transfor the signal to frequency domain
    Xds = fft([DS(p:p+N_FFT-1)'].*window');
    
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:N_FFT);
    
    % keep the signal
    z(p:p+N_FFT-1) = z(p:p+N_FFT-1) + s_est;
    
    
    % perform the same operation to noise and desired signal in order to calculate overall
    % SNR_output,if on need to observe SNR,comment this section
    XnoiseDS = fft(noiseDS(p:p+N_FFT-1)'.*window');
    noiseDS_filtered = W.*(XnoiseDS);
    iNoise = ifft(noiseDS_filtered);
    znoise(p:p+N_FFT-1) = znoise(p:p+N_FFT-1) + iNoise(1:N_FFT);
    
    Xspeech = fft(speech(p:p+N_FFT-1)'.*window');
    speech_filtered = W.*(Xspeech);
    iSpeech = ifft(speech_filtered);
    zspeech(p:p+N_FFT-1) = zspeech(p:p+N_FFT-1) + iSpeech(1:N_FFT);

    waitbar(p/(length(x(:,1))));
end
close(hwt);

% compute output SNR
SNR_output = 10*log10(sum(zspeech.^2)/sum(znoise.^2))
toc


