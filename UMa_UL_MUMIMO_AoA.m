clc;
clear all;
close all;
digits(100)
warning('off')

%% === Configure CSI-RS Resources  ===
carrier = nrCarrierConfig;
carrier.NSizeGrid = 1; % Bandwidth in RB
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'Normal';
ofdmInfo = nrOFDMInfo(carrier);


% === Configure Transmit and Receive Antenna Arrays ===
simParams.fc = 6e9; % freqRange = 'FR1';
simParams.c = physconst('LightSpeed');

simParams.lambda = simParams.c/simParams.fc;
simParams.NumTx = 1; 
simParams.NumRx = 8;
n = 0:simParams.NumRx-1;eta=pi;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% numFrames = 1; % 10 ms frames -> 1 subframe - 10 slots - 14 symbols
K = carrier.NSizeGrid * 12;
L = carrier.SymbolsPerSlot;
simParams.NUser = 1;
simParams.NPaths = 7;
SNRmin = 0; SNRmax = 20; step=1; nSNRLevels = length(SNRmin:step:SNRmax);% dB
nQPSKSym = 12;
SqErr_snr_AoA_MP_individual = zeros(simParams.NPaths,nSNRLevels);
SqErr_snr_AoA_DFT_individual = zeros(simParams.NPaths,nSNRLevels);
SqErr_snr_Gamma_MP = zeros(1,nSNRLevels);
SqErr_snr_Gamma_DFT = zeros(1,nSNRLevels);
SqErr_snr_Gamma_LinearMMSE = zeros(1,nSNRLevels);
SqErr_snr_Gamma_True = zeros(1,nSNRLevels); 
Sq_snr_Hmk_true = zeros(1,nSNRLevels);

NMSE_snr_Gamma_MP = zeros(1,nSNRLevels);
NMSE_snr_Gamma_DFT = zeros(1,nSNRLevels);
NMSE_snr_Gamma_MMSE = zeros(1,nSNRLevels);
NMSE_snr_Gamma_True = zeros(1,nSNRLevels);

snrIdx = 1;
totalMCTrials = 1; % Monte-Carlo trials
totalNumScatMIMOChannels = 100;
qpskRef = pskmod([0 1 2 3],4,pi/4);
tolFloating = 1e-2;
dist_CU_User = zeros(simParams.NUser,simParams.NPaths); 

BER_MP = zeros(nSNRLevels,1); 
BER_DFT = zeros(nSNRLevels,1);
BER_MMSE = zeros(nSNRLevels,1);
BER_theoretical = zeros(nSNRLevels,1);
EbN0 = zeros(nSNRLevels,1);
Beta = 4.1290e-10; % Empirical Large-scale fading coefficient
%% >>>>>>>>>>>>>>>> MAIN SIMULATION LOOP -- Pilot Training <<<<<<<<<<<<<<<<<<<
for snrdB=SNRmin:step:SNRmax
    snrPower = 10^(snrdB/10); % power
    
    fprintf(">>>>>>>>>>>>>>>>>>> Current SNR %d >>>>>>>>>>>>>>>>>>>>> \n", snrdB)
    for iMCTrial = 1:totalMCTrials
        ichannel = randi([1 totalNumScatMIMOChannels]);
        for iSlot = 0:carrier.SlotsPerFrame-1
            carrier.NSlot = iSlot; % current slot number
            [txInfo,channel,chanInfo,rxInfo]=genChannelOutput(ichannel,ofdmInfo,simParams,carrier,snrdB);
            txGrid = txInfo.txGrid;
            AoAs_True = txInfo.AoAs_True;
            Amk_True = txInfo.Amk_True;
            txAmp = txInfo.txAmp;
            offset = chanInfo.offset;
            pathGains = chanInfo.pathGains;
            pathDelays = chanInfo.pathDelays;
            maxChDelay = chanInfo.maxChDelay;
            fadWave = rxInfo.fadWave;
            rxArrayStv = rxInfo.rxArrayStv;
            % Correct timing offset
            fadWave = fadWave(1+offset:end,:);
            % Perform OFDM demodulation
            sigGrid = nrOFDMDemodulate(carrier,fadWave);
            sigPow = computeResourceGridPower(sigGrid,ofdmInfo.Nfft);
            txPow = computeResourceGridPower(txGrid,ofdmInfo.Nfft);
            noiseGrid = 1/sqrt(2).*complex(randn(size(sigGrid)),randn(size(sigGrid)));
            noisePow = computeResourceGridPower(noiseGrid,ofdmInfo.Nfft);
            rxGrid = sigGrid + noiseGrid;
            freqIdx = randi([1 K]); % in case of freq-selective fading
    
            %% ========= Perform Angle-Domain Channel Estimation Algo. ===============
            Xpilot = txGrid(freqIdx,1:nQPSKSym);
            Ypilot = reshape(rxGrid(freqIdx,1:nQPSKSym,:),[],nQPSKSym,simParams.NumRx);     % pilot recv: NPilotSymb X NRx
            Ypilot = permute(Ypilot,[3,2,1]); % NumRx X NumPilotSym
            Npilot = reshape(noiseGrid(freqIdx,1:nQPSKSym,:),[],nQPSKSym,simParams.NumRx);  % pilot recv: NPilotSymb X NRx
            Npilot = permute(Npilot,[3,2,1]); % NumRx X NumPilotSym
            pilotNorm = norm(Xpilot);
            
            AOAs_matpencil = zeros(1,simParams.NPaths);
            AOAs_dft = zeros(1,simParams.NPaths);
            for iSample = 1:nQPSKSym
                ySampled = Ypilot(:,iSample)*Xpilot(iSample)'./norm(Xpilot(iSample))^2;
                AOAs_matpencil = AOAs_matpencil + sort(matpencil_aoa(ySampled,simParams.NPaths),'ascend');
                AOAs_dft = AOAs_dft + reshape(sort(dft_aoa(ySampled,simParams.NumRx,simParams.NPaths),'ascend'),1,simParams.NPaths);
            end
            AOAs_matpencil = AOAs_matpencil ./ nQPSKSym;
            AOAs_dft = AOAs_dft ./ nQPSKSym;
            %%  Computing Angle Error Metrics 
	        for path=1:simParams.NPaths
		        SqErr_snr_AoA_MP_individual(path,snrIdx) = SqErr_snr_AoA_MP_individual(path,snrIdx) + abs(AoAs_True(1,path) - AOAs_matpencil(path)).^2;
		        SqErr_snr_AoA_DFT_individual(path,snrIdx) = SqErr_snr_AoA_DFT_individual(path,snrIdx) + abs(AoAs_True(1,path) - AOAs_dft(path)).^2;
            end
    
            %%  Channel Amplitude Estimation
	        Amk_hat_matpencil = rxArrayStv(simParams.fc,[AOAs_matpencil;zeros(1,simParams.NPaths)]);
	        Amk_hat_dft = rxArrayStv(simParams.fc,[AOAs_dft;zeros(1,simParams.NPaths)]);
            Gmk_hat_mp = Amk_hat_matpencil* (Amk_hat_matpencil\ (Ypilot*Xpilot'/pilotNorm^2)) ;
            Gmk_hat_dft = Amk_hat_dft * (Amk_hat_dft \ (Ypilot*Xpilot'/pilotNorm^2));
            Gmk_true = Amk_True * (Amk_True \ (Ypilot*Xpilot'/pilotNorm^2));
    
            %% Perform MMSE Channel Estimation
            h_LinMMSE = h_MMSE_CE(Ypilot,Xpilot,Npilot,Beta);
    
            %% ========= Perform True Channel Estimation Algo. ===============
            refInd = [freqIdx:12:freqIdx+nQPSKSym*12-1];
            refSym = txGrid(refInd);
            [hEstGrid,noiseEst] = nrChannelEstimate(rxGrid,txGrid); % 12x14x8
            hmk_true = reshape(hEstGrid(freqIdx,1,:),simParams.NumRx,[]);

            %% ========= Perform Error Estimation ===============
            sqErr_MP = mean(norm(hmk_true - Gmk_hat_mp)^2);
            sqErr_DFT = mean(norm(hmk_true - Gmk_hat_dft)^2);
            sqErr_MMSE = mean(norm(hmk_true - h_LinMMSE)^2);
            SqErr_snr_Gamma_MP(snrIdx) = SqErr_snr_Gamma_MP(snrIdx) + sqErr_MP;
	        SqErr_snr_Gamma_DFT(snrIdx) = SqErr_snr_Gamma_DFT(snrIdx) + sqErr_DFT;
	        SqErr_snr_Gamma_LinearMMSE(snrIdx) = SqErr_snr_Gamma_LinearMMSE(snrIdx) + sqErr_MMSE;

            
            sqNorm_Gamma_true = mean(norm(hmk_true)^2);
            nmse_MP = sqErr_MP / sqNorm_Gamma_true;
            nmse_DFT = sqErr_DFT / sqNorm_Gamma_true;
            nmse_MMSE = sqErr_MMSE / sqNorm_Gamma_true;
            NMSE_snr_Gamma_MP(snrIdx) = NMSE_snr_Gamma_MP(snrIdx) + nmse_MP;
            NMSE_snr_Gamma_DFT(snrIdx) = NMSE_snr_Gamma_DFT(snrIdx) + nmse_DFT;
            NMSE_snr_Gamma_MMSE(snrIdx) = NMSE_snr_Gamma_MMSE(snrIdx) + nmse_MMSE;
            
	        %% <<<<<<<<<<<<<<<< Channel Estimation End <<<<<<<<<<<<<<<<<<<

            %% >>>>>>>>>> SER Estimation Start >>>>>>>>>>>>
            symEnc = txGrid(freqIdx,1:nQPSKSym) ./ txAmp;
            symDecMP = Gmk_hat_mp'*Ypilot;
            symDecDFT = Gmk_hat_dft'*Ypilot;
            symDecMMSE = h_LinMMSE'*Ypilot;
            [symRx,symHest] = nrExtractResources(refInd,rxGrid,hEstGrid);
            symDecTrue = nrEqualizeMMSE(symRx,symHest,noiseEst);
        
            symEncQPSK = nrSymbolDemodulate(symEnc.','QPSK','DecisionType','Hard');
            symDemodMP = nrSymbolDemodulate(symDecMP.','QPSK','DecisionType','Hard');
            symDemodDFT = nrSymbolDemodulate(symDecDFT.','QPSK','DecisionType','Hard');
            symDemodMMSE = nrSymbolDemodulate(symDecMMSE.','QPSK','DecisionType','Hard');
            symDemodTrue = nrSymbolDemodulate(symDecTrue,'QPSK','DecisionType','Hard');
        
            berMP = biterr(symEncQPSK,symDemodMP);
            berDFT = biterr(symEncQPSK,symDemodDFT);
            berLMMSE = biterr(symEncQPSK,symDemodMMSE);
            berTrue = biterr(symEncQPSK,symDemodTrue);
            BER_MP(snrIdx) = BER_MP(snrIdx) + berMP;
            BER_DFT(snrIdx) = BER_DFT(snrIdx) + berDFT;
            BER_MMSE(snrIdx) = BER_MMSE(snrIdx) + berLMMSE;
            BER_theoretical(snrIdx) = BER_theoretical(snrIdx) + berTrue;
        end
    end
    

    BER_MP(snrIdx) = BER_MP(snrIdx) ./ (2*nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    BER_DFT(snrIdx) = BER_DFT(snrIdx) ./ (2*nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    BER_MMSE(snrIdx) = BER_MMSE(snrIdx) ./ (1*nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    BER_theoretical(snrIdx) = BER_theoretical(snrIdx) ./ (2*nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);

    Eb = mean(rms(sigGrid(freqIdx,1:nQPSKSym,:)),3)^2;
    N0 = mean(rms(noiseGrid(freqIdx,1:nQPSKSym,:)),3)^2;
    % <<<<<<<<<< BER Estimation End <<<<<<<<<<<<<<
    EbN0(snrIdx) = Eb/N0;
    
    NMSE_snr_Gamma_MP(snrIdx) = NMSE_snr_Gamma_MP(snrIdx) ./ (nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    NMSE_snr_Gamma_DFT(snrIdx) = NMSE_snr_Gamma_DFT(snrIdx) ./ (nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    NMSE_snr_Gamma_MMSE(snrIdx) = NMSE_snr_Gamma_MMSE(snrIdx) ./ (nQPSKSym*totalMCTrials*carrier.SlotsPerFrame);
    
    snrIdx = snrIdx + 1;
end

%% >>>>>>>>>>>>>>>>>>>>> Plotting Start >>>>>>>>>>>>>>>>>>>>
RMSE_snr_AoA_MP_individual = zeros(simParams.NPaths,nSNRLevels);
RMSE_snr_AoA_DFT_individual = zeros(simParams.NPaths,nSNRLevels);
RMSE_snr_Gamma_MP = sqrt(SqErr_snr_Gamma_MP / totalMCTrials);
RMSE_snr_Gamma_DFT = sqrt(vpa(SqErr_snr_Gamma_DFT) / totalMCTrials);
RMSE_snr_Hmk_LinearMMSE = sqrt(vpa(SqErr_snr_Gamma_LinearMMSE) / totalMCTrials);
RMSE_snr_Gamma_True = sqrt(vpa(SqErr_snr_Gamma_True) / totalMCTrials);
job=string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss"));

fig1=figure;
hold on; grid on
plot(SNRmin:step:SNRmax,20*log10(RMSE_snr_Gamma_MP),'Color','#0099ff','DisplayName','Hmk Matrix Pencil');
plot(SNRmin:step:SNRmax,20*log10(RMSE_snr_Gamma_DFT),'Color','red','DisplayName','Hmk DFT');
plot(SNRmin:step:SNRmax,20*log10(RMSE_snr_Hmk_LinearMMSE),'DisplayName','Hmk Linear-MMSE');
plot(SNRmin:step:SNRmax,20*log10(RMSE_snr_Gamma_True),'DisplayName','Hmk True-AoA');
xlabel('SNR (dB)')
ylabel('RMSE of Channel Response Estimation (dB scale)')
title(sprintf('RMSE of Channel Response Estimation: NumPilot=%d NAnt=%d',nQPSKSym,simParams.NumRx));
legend show
hold off
pngfile=sprintf('RMSE_snr_CSI_Est_j%s',job);
print(fig1,pngfile,'-dpng')

fig2=figure;
hold on; grid on
semilogy(SNRmin:step:SNRmax,NMSE_snr_Gamma_MP,'Color','#0099ff','DisplayName','Hmk Matrix Pencil');
semilogy(SNRmin:step:SNRmax,NMSE_snr_Gamma_DFT,'Color','red','DisplayName','Hmk DFT');
semilogy(SNRmin:step:SNRmax,NMSE_snr_Gamma_MMSE,'DisplayName','Hmk MMSE');
xlabel('SNR (dB)')
ylabel('NMSE of Channel Response Estimation')
title(sprintf('NMSE of Channel Response Estimation: NumPilot=%d NAnt=%d',nQPSKSym,simParams.NumRx));
legend show
hold off
pngfile=sprintf('NMSE_snr_CSI_Est_j%s',job);
print(fig2,pngfile,'-dpng')

fig3=figure;
hold on
grid on
semilogy(SNRmin:step:SNRmax,BER_MP,'o','DisplayName','BER Matrix-Pencil');
semilogy(SNRmin:step:SNRmax,BER_DFT,'.','DisplayName','BER DFT');
semilogy(SNRmin:step:SNRmax,BER_MMSE,'*','DisplayName','BER Linear-MMSE');
semilogy(SNRmin:step:SNRmax,BER_theoretical,'square','DisplayName','BER theoretical');
xlabel('Eb/N0 (dB)')
ylabel('BER')
title(sprintf('BER vs. Transmitter Power: NumPilot=%d NAnt=%d',nQPSKSym,simParams.NumRx));
legend show
hold off
pngfile=sprintf('BER_TxPower%s',job);
print(fig3,pngfile,'-dpng')

%% <<<<<<<<<<<<<<<<<< Plotting End <<<<<<<<<<<<<<<<<<<<<<<<

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%======== Channel Estimation ========
function h_LinMMSE = h_MMSE_CE(y,x,n,Beta)
% y                 = Frequency-domain received signal
%                   NumRx X NumSym
% d                 = Large-scale fading parameter
%                   NumRx X NumPaths
% x                 = pilot symbol
% noise
[~,nQPSKSym] = size(y);
yExtracted = y*x'/(norm(x)^2);
zeta = Beta / Beta + (norm(n)/nQPSKSym)^2;
W_MMSE = 1/(1+1/zeta);
h_LinMMSE = W_MMSE*yExtracted;
end



%% AoA Estimation Related
%% ===== Extended DFT + Angle-Rotation =====
function [AoA_Est_DFT] = dft_aoa(ymk_Sampled,N,L)
% Input:    yt, NumRx, N_MultiPath
% Output:   AOA_estimated, beta_estimated
FN = dftmtx(N)/sqrt(N);
Ndft_points = 100; %% can choose whatever # you want

AoA_Est_DFT = zeros(L,1);
hDFT = FN * ymk_Sampled;

% Coarse Peak Finding
% -- Find the central point (qInit) of each path
[qInits,isNeg] = findInitDFTBin(hDFT,N,L);

for l=1:L
    fNq = FN(:,qInits(l));
    ymk_DFT = fNq .* ymk_Sampled;
    
    angles_in_phi = [-Ndft_points/2: Ndft_points /2]*pi/ Ndft_points; %% Ndft_points in the phi domain
    st_vec_mtx = exp(1i* [0:N-1]' * angles_in_phi);  %% N \times Ndft_points matrix of Ndft_points steering vectors

    % Now if x is the data vector
    angle_info = abs(st_vec_mtx' * ymk_DFT);
    [~, max_angle_location] = max(angle_info);
    phi_location = angles_in_phi(max_angle_location);
    theta_init = 2*qInits(l)/N;
    
    if isNeg(l)
        theta = -theta_init + phi_location/pi; %% since \phi = kd\sin\theta
    else
        theta = theta_init - phi_location/pi;
    end

    if abs(theta) > 1
        theta = findNextPeak(theta,qInits(l),angle_info,angles_in_phi,N);
    end

    if isNeg(l)
        AoA_Est_DFT(l) = -1*real(acosd(theta));
    else
        AoA_Est_DFT(l) = real(acosd(theta));
    end
end
AoA_Est_DFT = sort(AoA_Est_DFT);
end


function [Q,isNeg] = findInitDFTBin(hDFT,N,L)
    [~,I] = sort(abs(hDFT),'descend');
    threshold = floor(N/2);
    Q = zeros(1,L);
    isNeg = zeros(1,L);
    pl = 1;
    
    for l=1:N
        if I(l) >= (threshold + 1)
            Q(pl) = I(l)-2;
            pl = pl+1;
            isNeg(l) = 1;
        else
            Q(pl) = I(l);
            pl = pl+1;
            isNeg(l) = 0;
        end
        if pl > L
            break;
        end
    end
end

function new_theta = findNextPeak(prev_theta,curr_qInit,angle_info,angles_in_phi,N)
    [ang_pks,ang_loc] = findpeaks(angle_info);
    [~,sorted_ang_loc_ind] = sort(ang_pks,'descend');
    ang_loc_sorted = ang_loc(sorted_ang_loc_ind);
    
    ang_locs_L = ang_loc_sorted(2:end);
    isNeg = 0; 
    new_theta = 2; % dummy init value
    if sign(prev_theta) < 1
        isNeg = 1;
    end
    idx=1; NPks = length(ang_locs_L);
    while (abs(new_theta) > 1 && idx <= NPks)
        curr_max_angle_loc = ang_locs_L(idx);
        curr_phi_loc = angles_in_phi(curr_max_angle_loc);
        if isNeg
            new_theta = -2*curr_qInit/N + curr_phi_loc/pi;
        else
            new_theta = 2*curr_qInit/N - curr_phi_loc/pi;
        end
        idx = idx+1;
    end
end


%% ====== AoA estimation using Matrix Pencil ======
function AOAs = matpencil_aoa(ymk_Sampled,L)
% ymk_Sampled = Nx1
% y = As + n
% N: Number of array elements
% M: Number of paths
% L: Matrix Pencil parameter

[N,~] = size(ymk_Sampled);
P = L; % size of window
ymk_Sampled = resample(ymk_Sampled,4,2);
x = ymk_Sampled(:,1); % Kx1=>only one time sample
Y1 = zeros(N-P,P);
Y2 = zeros(N-P,P);
for p=1:P
    Y1(:,p) = x(p:N-P-1+p,1);
    Y2(:,p) = x(p+1:N-P+p,1);
end

Y1_pinv = (Y1'*Y1)\Y1';
z_hat = sort(eig(Y1_pinv*Y2),'descend');

% AOAs = sort(acosd(imag(log(z_hat(1:L))) ./ pi));
% AOAs = AOAs';

% find the upper bound (need |e^i pi cos(theta)| <=1 )
[~,i] = mink(abs(abs(z_hat)-1),L);
z_hat_cap = z_hat(i);
AOAs = zeros(1,L);

for pl=1:L
    z_tmp = z_hat_cap(pl);

    theta_tmp = acosd(imag(log(z_tmp)) ./ pi);
    theta_tmpReal = asind(imag(log(z_tmp)) ./ pi);
    if sign(theta_tmpReal) > 0
        AOAs(pl) = theta_tmp;
    else
        AOAs(pl) = -theta_tmp;
    end
end
% AOAs = sort(AOAs);
end


function [txInfo,channel,chanInfo,rxInfo]=genChannelOutput(ichannel,ofdmInfo,simParams,carrier,snrdB)
% Configure transmitter/receiver/Scatterers positions
refax = [[1;0;0] [0;1;0] [0;0;0]];
chanFileName = fullfile("ScatMIMOChanModels",sprintf("chanModel%d.mat",ichannel));
file = java.io.File(chanFileName);
fullpath = char(file.getAbsolutePath());
try
    chanModel = load(fullpath);
catch
    chanFileName = fullfile("ScatMIMOChanModels",sprintf("chanModel%d.mat",1));
    file = java.io.File(chanFileName);
    fullpath = char(file.getAbsolutePath());
    chanModel = load(fullpath);
end
channel = chanModel.channel;
simParams.scatPos = channel.ScattererPosition;
simParams.posTx = channel.TransmitArrayPosition;
simParams.posRx = channel.ReceiveArrayPosition;
[~,~,tau] = channel(complex(randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx), ...
    randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx)));
maxChDelay = ceil(max(tau)*ofdmInfo.SampleRate);

% Transmitter Setup
txAmp = 10000*10^(snrdB/20);
txGrid = nrResourceGrid(carrier,simParams.NumTx);
qamSymbols = nrSymbolModulate(randi([0,1],numel(txGrid)*2,1),'QPSK');
txGrid(:) = txAmp.*qamSymbols;


%% === Send the Waveform through the Channel ===
% OFDM Modulation
[txWaveform,~] = nrOFDMModulate(carrier,txGrid);
% Append zeros to the transmit waveform to account for channel delay
txWaveform = [txWaveform; zeros(maxChDelay,simParams.NumTx)];
% Pass the waveform through the channel
[fadWave,pathGains,pathDelays] = channel(txWaveform);

% Estimate timing offset
%% >>>>>>>>>>>>>>>>>> Channel Estimation Start >>>>>>>>>>>>>>>>>>>>
offset = nrTimingEstimate(carrier,fadWave,txGrid);
if offset > maxChDelay
    offset = 0;
end

% Receiver Setup
rxArrayStv = phased.SteeringVector('SensorArray',channel.ReceiveArray,'PropagationSpeed',simParams.c);
% Compute True AoA based on ScatPos and posRx
[~,AoAs_True] = rangeangle(simParams.scatPos,simParams.posRx,refax);
AoAs_True(1,:) = sort(AoAs_True(1,:),'ascend');
Amk_True = rxArrayStv(simParams.fc,AoAs_True(1,:));

%% >>>>>>>>>> Output Args >>>>>>>>>>>>
txInfo.txGrid = txGrid;
txInfo.AoAs_True = AoAs_True;
txInfo.Amk_True = Amk_True;
txInfo.txAmp = txAmp;
chanInfo.pathGains = pathGains;
chanInfo.pathDelays = pathDelays;
chanInfo.maxChDelay = maxChDelay;
chanInfo.offset = offset;
rxInfo.fadWave = fadWave;
rxInfo.rxArrayStv = rxArrayStv;


%% Plot the Scattering MIMO Scenario
% Configure the MIMO scene parameters
[~,scatRxAng] = rangeangle(simParams.scatPos(:,1),simParams.posRx,refax);
azRxBeamWidth = 30; % In degrees
elRxBeamWidth = 30; % In degrees
rxAng = getInitialBeamDir(scatRxAng,azRxBeamWidth,elRxBeamWidth);
wR = rxArrayStv(simParams.fc,rxAng);
% sceneParams.TxArray = channel.TransmitArray;
% sceneParams.RxArray = channel.ReceiveArray;
% sceneParams.TxArrayPos = simParams.posTx;
% sceneParams.RxArrayPos = simParams.posRx;
% sceneParams.ScatterersPos = simParams.scatPos;    
% sceneParams.Lambda = simParams.lambda;
% sceneParams.ArrayScaling = 100;   % To enlarge antenna arrays in the plot
% sceneParams.MaxTxBeamLength = 10; % Maximum length of transmit beams in the plot
% sceneParams.MaxRxBeamLength = 25; % Maximum length of receive beam in the plot
% wT = 1;
% hPlotSpatialMIMOScene(sceneParams,wT,wR);
% view(2)

end


function beamDir = getInitialBeamDir(scatAng,azBeamWidth,elBeamWidth)
%   getInitialBeamDir returns the initial beam direction BEAMDIR, given the
%   angle of scatterer position with respect to transmit or receive antenna
%   array SCATANG, beamwidth of transmit or receive beam in azimuth plane
%   AZBEAMWIDTH, and beamwidth of transmit or receive beam in elevation
%   plane ELBEAMWIDTH.

    % Azimuth angle boundaries of all transmit/receive beams
    azSSBSweep = -180:azBeamWidth:180;
    % Elevation angle boundaries of all transmit/receive beams
    elSSBSweep = -90:elBeamWidth:90;
    
    % Get the azimuth angle of transmit/receive beam
    azIdx1 = find(azSSBSweep <= scatAng(1),1,'last');
    azIdx2 = find(azSSBSweep >= scatAng(1),1,'first');
    azAng = (azSSBSweep(azIdx1) + azSSBSweep(azIdx2))/2;
    
    % Get the elevation angle of transmit/receive beam
    elIdx1 = find(elSSBSweep <= scatAng(2),1,'last');
    elIdx2 = find(elSSBSweep >= scatAng(2),1,'first');
    elAng = (elSSBSweep(elIdx1) + elSSBSweep(elIdx2))/2;
    
    % Form the azimuth and elevation angle pair (in the form of [az;el])
    % for transmit/receive beam
    beamDir = [azAng;elAng];
end
