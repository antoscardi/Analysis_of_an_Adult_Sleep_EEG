clear variables; close all; clc;

%% Plot the signal EEG to analize the different waves
load record.mat;
record=record-mean(record); % we subtract the mean in order to have a zero mean signal
L=length(record); 
fs_original=512;

% the filtfilt function requires the input variables to be class double, 
% the double function didn't work so we are changing every variable from
% single to double.
 s = whos;
 for i = 1:length(s)
      if strcmp(s(i).class,'single')
         name = s(i).name;
         assignin('base', name, double(evalin('base', name)));
      end
 end

%Anti-aliasing FIR low-pass filter
Fs = 512;
Hd = fdesign.lowpass('Fp,Fst,Ap,Ast',60,70,0.1,100,Fs);
LP = design(Hd,'equiripple');
fvtool(LP,'Analysis','freq','Fs',512);

LP_filtered= filter(LP , record); %we filter the signal using filtfilt function

% We downsample the signal to have less compuational complexity
fs_down=128; % The characteristic frequency of the EEG signals is between 0.5-60 Hz, so we select a downsampling frequency of 128 Hz > 60*2=120 Hz
record_down=downsample(LP_filtered, fs_original/fs_down);
L_down=length(record_down);

% Definition of the time vectors
t=0:1/fs_down:L_down/fs_down-1/fs_down;
t_o=0:1/fs_original:L/fs_original-1/fs_original;

% EEG graph
figure(1), 
plot(t_o/60, record); 
xlabel('Time [min]');xlim([0 540]);ylabel('Amplitude [microV]');ylim([-.3e03 +.3e03]);xticks(0:90:540);title('EEG');
legend('Normal');
hold on;
plot(t/60, record_down);
legend('Original','Downsampled');
fs=fs_down; % For the sake of semplicity we rename the downsampling frequency as 'fs' and the length as 'L'
L=L_down;


% NOISE CANCELLATION with a bandpass filter at 0.5-35Hz
Fn = fs/2;                                      % Nyquist Frequency
Wp = [0.5   34.5]/Fn;                           % Passband Frequency (Normalised)
Ws = [0.4   40.5]/Fn;                           % Stopband Frequency (Normalised)
Rp =   1;                                       % Passband Ripple (dB)
Rs = 150;                                       % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                 % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                      % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                    % Convert To Second-Order-Section For Stability

h=fvtool(sosbp,'Analysis','freq','Fs',128);     % Plotting filter response



EEG_filtered = filtfilt(sosbp, gbp, record_down);
figure(3),
plot(t/60,EEG_filtered, t/60,record_down);
xlabel('Time [min]');xlim([0 540]);ylabel('Amplitude [microV]');ylim([-280 280]);xticks(0:90:540);title('Noise filtered cheb');
legend('Noise Filtered','Downsampled');    

% Spectrogram of the whole signal, with windows of 20 min 
figure(4)
spectrogram(EEG_filtered, 307200, 'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis


% BANDPASS FILTERING
% Designing different bandpass filters to extract alpha,beta,theta,delta waves

% ALPHA waves extraction (8-13 Hz), filters: 
Fn = fs/2;                                      % Nyquist Frequency
Wp = [8.5   12.5]/Fn;                           % Passband Frequency (Normalised)
Ws = [7.5   13.5]/Fn;                           % Stopband Frequency (Normalised)
Rp =   1;                                       % Passband Ripple (dB)
Rs = 150;                                       % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                 % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                      % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                    % Convert To Second-Order-Section For Stability
fvtool(sosbp,'Analysis','freq');                % Plotting filter response

% Filtering and plotting the waves
alpha_filtered = filtfilt(sosbp, gbp, EEG_filtered);
figure(5), 
plot(t/60,EEG_filtered,t/60,alpha_filtered);
xlabel('Time [min]');xlim([180 180.3]);ylabel('Amplitude [microV]');ylim([-40 40]);xticks(180:0.05:180.3);title('Alpha filtered cheb');
legend('Noise filtered','Alpha wave');         

%spectrogram of the alpha waves
figure(6),
spectrogram(alpha_filtered,307200,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('alpha');

%PSD of alpha waves
y=fft(EEG_filtered); %Fourier transform of the EEG
Px=y.*conj(y)/L;
freq=[0:fs/L:fs-1/L]';
y1=fft(alpha_filtered);  %Fourier transform of the alpha waves
Px1=y1.*conj(y1)/L;
figure(7),
subplot(2,1,1); plot(freq(1:L/2),abs(y(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('EEG Spectrum');legend('Noise Filtered','alpha');
subplot(2,1,2); plot(freq(1:L/2),abs(y1(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('Alpha Spectrum');legend('Noise Filtered','alpha');
figure(8),
loglog(freq(1:L/2),(Px(1:L/2)).^10);hold on; loglog(freq(1:L/2),(Px1(1:L/2)).^10);
legend('Noise Filtered','alpha');xlabel('Frequency [Hz]');xlim([0 100]);ylabel('PSD (dB)');title('Alpha PSD');legend('Noise Filtered','alpha');



% BETA waves extraction (14-35 Hz):                               
Wp = [14.5   34.5]/Fn;                          % Passband Frequency (Normalised)
Ws = [13.5   35.5]/Fn;                          % Stopband Frequency (Normalised)
Rp =  1;                                        % Passband Ripple (dB)
Rs = 150;                                       % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                 % Filter Order
[z,b,k] = cheby2(n,Rs,Ws);                      % Filter Design
[sosbp,gbp] = zp2sos(z,b,k);                    % Convert To Second-Order-Section For Stability
fvtool(sosbp,'Analysis','freq');                % Plotting filter response

% Filtering and plotting the waves
beta_filtered = filtfilt(sosbp, gbp, EEG_filtered);
figure(9), 
plot(t/60,EEG_filtered,t/60,beta_filtered);
xlabel('Time [min]');xlim([180 180.3]);ylabel('Amplitude [microV]');ylim([-40 40]);xticks(180:0.05:180.3);title('Beta filtered cheb');
legend('Noise filtered','Beta wave'); 

%Plotting the spectrogram of beta waves
figure(10)
spectrogram(beta_filtered,307200,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('beta'); 

%PSD of beta waves
y1=fft(beta_filtered);
Px1=y1.*conj(y1)/L;
figure(11),
subplot(2,1,1); plot(freq(1:L/2),abs(y(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('EEG Spectrum');legend('Noise Filtered','beta');
subplot(2,1,2); plot(freq(1:L/2),abs(y1(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('Beta Spectrum');legend('Noise Filtered','beta');
figure(12),
loglog(freq(1:L/2),(Px(1:L/2)).^10);hold on; loglog(freq(1:L/2),(Px1(1:L/2)).^10);
legend('Noise Filtered','beta');xlabel('Frequency [Hz]');xlim([0 100]);ylabel('PSD (dB)');title('Beta PSD');legend('Noise Filtered','beta');



% THETA waves extraction (4-7 Hz)                               
Wp = [4.5   6.5]/Fn;                            % Passband Frequency (Normalised)
Ws = [3.5   7.5]/Fn;                            % Stopband Frequency (Normalised)
Rpb =  1;                                       % Passband Ripple (dB)
Rs= 150;                                        % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rpb,Rs);                % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                      % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                    % Convert To Second-Order-Section For Stability
fvtool(sosbp,'Analysis','freq');                % Plotting filter response

% Filtering and plotting the waves
theta_filtered = filtfilt(sosbp, gbp, EEG_filtered);
figure(13),
plot(t/60,EEG_filtered,t/60,theta_filtered);
xlabel('Time [min]');xlim([180 180.3]);ylabel('Amplitude [microV]');ylim([-40 40]);xticks(180:0.05:180.3);title('Theta filtered cheb');
legend('Noise filtered','Theta wave'); 

%Spectrogram of theta waves 
figure(14),
spectrogram(theta_filtered,307200,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('theta');

%PSD of theta waves
y1=fft(theta_filtered);
Px1=y1.*conj(y1)/L;
figure(15),
subplot(2,1,1); plot(freq(1:L/2),abs(y(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('EEG Spectrum');legend('Noise Filtered','theta');
subplot(2,1,2); plot(freq(1:L/2),abs(y1(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 50]);ylabel('Amplitude [microV]');xticks(0:2:50);title('Theta Spectrum');legend('Noise Filtered','theta');
figure(16),
loglog(freq(1:L/2),(Px(1:L/2)).^10);hold on; loglog(freq(1:L/2),(Px1(1:L/2)).^10);
legend('Noise Filtered','theta');xlabel('Frequency [Hz]');xlim([0 100]);ylabel('PSD (dB)');title('Theta PSD');legend('Noise Filtered','theta');



% DELTA waves extraction (0.5-4 Hz)
%Chebychef                               
Wp= [0.6   3.8]/Fn;                             % Passband Frequency (Normalised)
Ws = [0.4  4.2]/Fn;                             % Stopband Frequency (Normalised)
Rpb = 1;                                        % Passband Ripple (dB)
Rs = 150;                                       % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rpb,Rs);                % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                      % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                    % Convert To Second-Order-Section For Stability
fvtool(sosbp,'Analysis','freq');                % Plotting filter response

% Filtering and plotting the waves
delta_filtered = filtfilt(sosbp, gbp, EEG_filtered);
figure(17),
plot(t/60,EEG_filtered,t/60, delta_filtered);
xlabel('Time [min]');xlim([180 180.3]);ylabel('Amplitude [microV]');ylim([-40 40]);xticks(180:0.05:180.3);title('Delta filtered cheb');
legend('Noise filtered','Delta wave');         
            
%Spectrogram of theta waves
figure(18),
spectrogram(delta_filtered,307200,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('delta');

%PSD of delta waves
y1=fft(delta_filtered);
Px1=y1.*conj(y1)/L;
figure(19),
subplot(2,1,1); plot(freq(1:L/2),abs(y(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 10]);ylabel('Amplitude [microV]');xticks(0:0.5:10);title('EEG Spectrum');legend('Noise Filtered','delta');
subplot(2,1,2); plot(freq(1:L/2),abs(y1(1:L/2)));
xlabel('Frequency [Hz]');xlim([0 10]);ylabel('Amplitude [microV]');xticks(0:0.5:10);title('Delta Spectrum');legend('Noise Filtered','delta');
figure(20),
loglog(freq(1:L/2),(Px(1:L/2)).^10);hold on; loglog(freq(1:L/2),(Px1(1:L/2)).^10);
legend('Noise filtered','delta');xlabel('Frequency [Hz]');xlim([0 100]);ylabel('PSD (dB)');title('Delta PSD');legend('Noise Filtered','delta');

%SLEEP STAGES ANALYSIS
%The time windows of the spectrogram have been choosen according to the 
%usual time for every sleep stage, present in literature.  

%AWAKE
%Plotting spectrogram of both alpha and beta waves in order to see when the patient is AWAKE  
AWAKE = alpha_filtered + beta_filtered;  
figure(21)
spectrogram(AWAKE,153600,'yaxis');
title('AWAKE');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
%it seems that the patient fell asleep after 20 min the recording started 

%STAGE 1
%Plotting spectrogram of theta waves in order to see when the patient is in STAGE 1  
S1 = theta_filtered;  
figure(22),
spectrogram(S1,76800,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('STAGE 1');

%STAGE 3 and 4  
%Plotting spectrogram of delta waves in order to see when the patient is in STAGES 3/4  
S3_4 = delta_filtered;  
figure(23)
spectrogram(S3_4,230400,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('STAGE 3/4');

%REM
%Plotting spectrogram of both beta and theta waves in order to distinguish
%REM phases from the others  (clim = [-5,15]
REM = beta_filtered + theta_filtered;
figure(24),
spectrogram(REM,153600,'yaxis');
hfig=gcf; hfig.CurrentAxes.CLim = [0 15]; %scaling the intensity in dBm of the spectrogram for better analysis
title('REM');
%first rem after 70 min he fell asleep lasted 10 min or 40 min 


