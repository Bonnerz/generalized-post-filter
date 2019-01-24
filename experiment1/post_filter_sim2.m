%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence 
%  simulation configuration:
%    use RIR method to generate input signal 
%    additive noise type:diffuse noise
%    array type:5 mic ULA,5cm
%  dependencies:
%    RIR-Generator
%    Signal-Generator
%
%  Author: wangwei
%  Data  : 6/15/2017
%  refer to:
%  "Microphone Array Post-Filter Based on Noise Field Coherence" IEEE 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
c = 340; % speed of sound

%% generate diffuse noise
[ noise,Pos,dm] = noise_gen_ULA;

noise = noise'/20;
N = size(noise,2);        %Channels
%% 

pathname = '../sound/';

%use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'an101-mtms-senn3.wav']);

pathname = '../sound/speech/';
               
%% use RIR to generate reverberation output
d = 0.05;
angle = [90 0]/180*pi;          % source direction [0,180]
[x1,y1,z1]=sph2cart(angle(1),angle(2),0.7);    % source position 1
source_pos = [x1,y1,z1];              % Source position [x y z] (m)

%source_pos = [2-d*2 1.5+0.7 2];              % Source position [x y z] (m)

scale = 10;

h = RIR_generator( source_pos,0.42);
h = h*scale;
h = h(:,1:1024);

x1 = conv(speech,h(1,:));
x2 = conv(speech,h(2,:));
x3 = conv(speech,h(3,:));
x4 = conv(speech,h(4,:));
x5 = conv(speech,h(5,:));
s = [x1,x2,x3,x4,x5];
s = s(1:end-1024,:);

% signal+difuse noise
if(size(s,1)>size(noise,1))
    x = s(1:size(noise,1),:)+noise;
else
    x = s+noise(1:size(s,1),:);
end

%% 
% DS output
% perfect time alignment
DS = sum(x,2)/N;  

% compute input SNR
noiseDS = sum(noise,2);
SNR_input = 10*log10(sum(speech.^2)/sum(noiseDS.^2))

%% 
N_FFT = 256;

window = hamming(N_FFT);


P_len = N_FFT/2+1;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);
Pxii_pre = ones(N,P_len);
Pssnn_pre = ones(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);
Pxij_pre = ones((N*N-N)/2,P_len);
Pxij_curr = ones((N*N-N)/2,P_len);

Pss = zeros(1,P_len);

Pss_e = zeros((N*N-N)/2,P_len);

z = zeros(1,length(x(:,1)));
znoise = zeros(1,length(x(:,1)));
zspeech = zeros(1,length(x(:,1)));
t = 1;
%%
tic

% array spacing
r = 0.05; 

f = 0:fs/256:fs/2;
w = 2*pi*fs*(0:N_FFT/2)/N_FFT;

M = N;

alpha = 0.3;

% fixed wideband MVDR using pre-defined noise cross-spectral matrix
[ MVDR_out, x1,H,DI,Fvv] = superdirectiveMVDR_ULA(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle(1)-90/180*pi);
DS = MVDR_out;
% [sc,F]=mycohere(x1(:,1),x1(:,2),256,fs,hanning(256),0.75*256);
Inc = 128; 
k_optimal = 1;
hwt = waitbar(0,'general poster filter');
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-N_FFT
    for i = 1:N
        Xi = fft(x(p:p+N_FFT-1,i).*window);
        Pxii_curr = abs(Xi).^2;
        % eq.11
        Pxii(i,:) = alpha*Pxii_pre(i,:)+(1-alpha)*Pxii_curr(1:N_FFT/2+1).';      
    end
    Pxii_pre = Pxii;
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            
            Xi = fft(x(p:p+N_FFT-1,i).*window).';
            Xj = fft(x(p:p+N_FFT-1,j).*window).';
            % cross-spectral
            Pxij_temp = Xi.*conj(Xj);
            % half bin
            Pxij_curr(t,:) = Pxij_temp(1:N_FFT/2+1);
            % average
            Pxij(t,:) = alpha*Pxij_pre(t,:)+(1-alpha)*Pxij_curr(t,:);
                 
            % eq.22 estimate source signal's PSD
            Pss_e(t,:) = (real(Pxij(t,:)) - 0.5*real(Fvv(:,i,j)').*(Pxii(i,:)+Pxii(j,:)))...
                         ./...
                         (ones(1,P_len) - real(Fvv(:,i,j)'));
             t = t+1;
        end
    end
    Pxij_pre = Pxij;
    t = 1;
    % eq.23 
    % take the average of multichanel signal to improve robustness
    Pss = sum(Pss_e)*2/(N*N-N); 
    
    % handle the indeterminite soulution when MSC¡Ö1
    % Pss(Pss<0) = 1e-3;
    
    % eq.23 
    % calculate the frequency domain filter coefficient
    W_e = real(Pss)./Pssnn;

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
    
    Xspeech = fft(DS(p:p+N_FFT-1)'.*window');
    speech_filtered = W.*(Xspeech);
    iSpeech = ifft(speech_filtered);
    zspeech(p:p+N_FFT-1) = zspeech(p:p+N_FFT-1) + iSpeech(1:N_FFT);

    waitbar(p/(length(x(:,1))));
end
close(hwt);

% compute output SNR
SNR_output = 10*log10(sum(zspeech.^2)/sum(znoise.^2))
toc


