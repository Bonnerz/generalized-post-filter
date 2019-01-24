function [ DS, x1,H,DI,Fvv] = superdirectiveMVDR_ULA( x,fs,N,frameLength,inc,r,angle,Fvv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% superdirective beamformer using uniform linear array
%   
%      input :
%          x : input signal ,samples * channel
%          fs: sample rate
%          N : fft length,frequency bin number
%frameLength : frame length,usually same as N
%        inc : step increment
%          d : array element spacing
%      angle : incident angle,[0,180]
%
%     output :
%         DS : delay-sum output
%         x1 : presteered signal,same size as x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 256;
% inc = 32;
% frameLength = 256;
c = 340;
Nele = size(x,2);
omega = zeros(frameLength,1);
H = ones(N/2+1,Nele)';

theta = angle;%*pi/180; % 方位角 0 < angle <180,attention for pi

tao = r*cos(theta)*[0:Nele-1]/c;     %

yds = zeros(length(x(:,1)),1);
x1 = zeros(size(x));
DI = zeros(N/2+1,1);

M = Nele;
N_FFT = N;
f = 0:fs/256:fs/2;
Fvv = zeros(N_FFT/2+1,M,M);
k_optimal = 1;%1.7;
for i = 1:M
    for j = 1:M   
        if i == j
            Fvv(:,i,j) = ones(N_FFT/2+1,1);
        else
            dij = r*abs(i-j);
            Fvv(:,i,j) = sin(2*pi*f*dij*k_optimal/c)./(2*pi*f*dij*k_optimal/c);Fvv(1,i,j) = 0.998;%T(2) = 0.996;
        end
    end
end


for k = 1:N/2+1
%         inv_Fvv = inv(squeeze(Fvv(k,:,:)));
        omega(k) = 2*pi*(k-1)*fs/N;    
        
        % 方向向量
        d = exp(-1j*omega(k)*tao');
%         H(k,:) = [1;exp(-1j*omega(k)*tao);exp(-1j*omega(k)*2*tao)];
        %对齐向量，以第一个阵元为参考，
        %例如若第一个阵元最慢(theata>0),则将第2、3、....个阵元分别延迟exp(-j*w*m*tao)
%         H(k,:) = [1;exp(-1j*omega(k)*tao);exp(-1j*omega(k)*2*tao);];

        % MVDR soulution
%         Fvvk = diag(ones(1,Nele));%squeeze(Fvv(k,:,:));
        Fvv_k = (squeeze(Fvv(k,:,:))+1e-3*eye(Nele));
        if(1)%k~=31&&k~=20&&k~=16)
        H(:,k) =    Fvv_k\d ...
                 ./(d'/Fvv_k*d);
%        H(:,k) =   d'/M; % DS weights
        DI(k) = (abs(H(:,k)'*d))^2 ...
                /(H(:,k)'*squeeze(Fvv(k,:,:))*H(:,k));
        else
         H(:,k) = d;
        end
end
for i = 1:inc:length(x(:,1))-frameLength

%     d = zeros(N,Nele);
%     d(33,:) = 1;
    d = fft(x(i:i+frameLength-1,:).*hamming(frameLength));
%     x_fft = dot(H.',d(1:129,:),2);
    x_fft = H.'.*d(1:N/2+1,:);
    yf = sum(x_fft,2);
    Cf = [yf;conj(flipud(yf(2:N/2)))];
    
    % 恢复延时累加的信号
    yds(i:i+frameLength-1) = yds(i:i+frameLength-1)+(ifft(Cf));
    
    
    % 恢复各路对齐后的信号
    xf  = [x_fft;conj(flipud(x_fft(2:N/2,:)))];
    x1(i:i+frameLength-1,:) = x1(i:i+frameLength-1,:)+(ifft(xf));
end
DS = real(yds); 
x1 = real(x1);
DI = pow2db(abs(DI));

end

