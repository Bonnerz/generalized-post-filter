%--------------------------------------------------------------------------
%
% Example script for cinf_1D.m and cinf_3D.m
%
% Author        : dr.ir. Emanuel A.P. Habets
% Date          : 16-09-2010
%
% Related paper : E.A.P. Habets and S. Gannot, 'Generating sensor signals
%                 in isotropic noise fields', Submitted to the Journal
%                 of the Acoustical Society of America, May, 2007.
%
% Comment       :
%
%--------------------------------------------------------------------------

clear;
close all;

%% Initialization
fs = 8000;                       % Sample frequency
NFFT = 256;                      % Number of frequency bins (for analysis)
w = 2*pi*fs*(0:NFFT/2)/NFFT; 
c = 340;                         % Speed of sound
L = 2^18;                        % Data length
r = 0.032;
[x1,y1,z1]=sph2cart(0,0,r);    % Sensor position 1
[x2,y2,z2]=sph2cart(90/180*pi,0,r);    % Sensor position 2
[x3,y3,z3]=sph2cart(180/180*pi,0,r);    % Sensor position 1
[x4,y4,z4]=sph2cart(270/180*pi,0,r);    % Sensor position 2
P = [x1 x2 x3 x4;y1 y2 y3 y4;z1 z2 z3 z4]; % Construct position matrix
M = size(P,2);                           % Number of sensors

% Calculate sensor distances w.r.t. sensor 1
d = zeros(1,M);
for m = 2:M
    d(m) = norm(P(:,m)-P(:,1),2);
end

%% Generate sensor signals
params.c = c;
params.fs = fs;
params.N_phi = 64;

% 1D example
z = cinf_1D(d,L,params); 

% 3D example
% z = cinf_3D(P,L,params); 

%% Calculate spatial coherences
sc_sim = zeros(M-2,NFFT/2+1);
sc_theory = zeros(M-2,NFFT/2+1);
for m = 1:M-1
    [sc,F]=mycohere(z(1,:)',z(m+1,:)',NFFT,fs,hanning(NFFT),0.75*NFFT);
    sc_sim(m,:) = real(sc');

    sc_theory(m,:) = besselj(0,w*d(m+1)/c);
end

%% Plot results
% Sensor pair 1-2
figure(1); 
m=1;                                   
plot(F/1000,sc_sim(m,:),'k')
hold on;
plot(F/1000,sc_theory(m,:),'--k')
hold off;
xlabel('Frequency [kHz]');
ylabel('Spatial Coherence');
title(sprintf('Distance %1.2f m',d(m+1)));
set(gca,'DataAspectRatio',[1 0.75 1]);
legend('Simulation','Theory');
grid on;

% Sensor pair 1-3
figure(2);
m=2; 
plot(F/1000,sc_sim(m,:),'k')
hold on;
plot(F/1000,sc_theory(m,:),'--k')
hold off;
xlabel('Frequency [kHz]');
ylabel('Spatial Coherence');
title(sprintf('Distance %1.2f m',d(m+1)));
set(gca,'DataAspectRatio',[1 0.75 1]);
legend('Simulation','Theory');
grid on;