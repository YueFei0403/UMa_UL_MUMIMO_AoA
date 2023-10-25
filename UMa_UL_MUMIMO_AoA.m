clc;
clear all;
close all;


%% === Configure CSI-RS Resources  ===
carrier = nrCarrierConfig;
carrier.NSizeGrid = 52; % Bandwidth in RB
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'Normal';
ofdmInfo = nrOFDMInfo(carrier);


% === Configure Transmit and Receive Antenna Arrays ===
simParams.fc = 6e9; % freqRange = 'FR1';
simParams.c = physconst('LightSpeed');

simParams.lambda = simParams.c/simParams.fc;
simParams.NumTx = 1; 
simParams.NumRx = 8;
% Configure the transmit and receive antenna elements
txAntenna = phased.ShortDipoleAntennaElement;                  % To avoid transmission beyond +/- 90
                                                                 % degrees from the broadside, baffle
                                                                 % the back of the transmit antenna
                                                                 % element by setting the BackBaffled
                                                                 % property to true
rxAntenna = phased.IsotropicAntennaElement('BackBaffled',false); % To receive the signal from 360 degrees,
                                                                 % set the BackBaffled property to false

simParams.txArray = phased.NRRectangularPanelArray('Size',[1, 1, 1, 1],'ElementSet', {txAntenna},...
            'Spacing',[0.5*simParams.lambda,0.5*simParams.lambda,3*simParams.lambda,3*simParams.lambda]);
simParams.rxArray = phased.ULA('Element',rxAntenna, ...
    'NumElements',simParams.NumRx,'ElementSpacing',0.5*simParams.lambda,'ArrayAxis','x');

% Configure transmitter/receiver/Scatterers positions
% Configure transmitter/receiver/Scatterers positions
refax = [[1;0;0] [0;1;0] [0;0;0]];

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% numFrames = 1; % 10 ms frames -> 1 subframe - 10 slots - 14 symbols
K = carrier.NSizeGrid * 12;
L = carrier.SymbolsPerSlot;
carrier.NFrame = 10;

% Configure SRS
srs = nrSRSConfig;
srs.NumSRSSymbols = 4;
srs.SymbolStart = 0;
srs.NumSRSPorts = simParams.NumTx;
srs.FrequencyStart = 0;
srs.CSRS = 14; % use the whole bandwidth
srs.KTC = 2;
srs.Repetition = 2;
srs.SRSPeriod = [1 0]; % present in all slots
srs.ResourceType = 'periodic';
simParams.NUser = 1;
simParams.NPaths = 7;
SNRmin = 1; SNRmax = 30; step=1; nSNRLevels = length(SNRmin:step:SNRmax);% dB
firstSymb = 1; lastSymb = srs.NumSRSSymbols; % 3GPP Specs
freqSCRange = 1:srs.KTC:K;    % subcarrier range in freq range
NumFreqSymb = length(freqSCRange);
NumPilotSymb = length(firstSymb:1:lastSymb); 
SqErr_snr_AoA_MP_individual = zeros(simParams.NPaths,nSNRLevels);
SqErr_snr_AoA_DFT_individual = zeros(simParams.NPaths,nSNRLevels);
SqErr_snr_Gamma_MP = zeros(1,nSNRLevels);
SqErr_snr_Gamma_DFT = zeros(1,nSNRLevels);
SqErr_snr_Hmk_LinearMMSE = zeros(1,nSNRLevels);
snrIdx = 1;
totalMCTrials = 5000; % Monte-Carlo trials
totalNumScatMIMOChannels = 100;
qpskRef = pskmod([0 1 2 3],4,pi/4);
tolFloating = 1e-2;
dist_CU_User = zeros(simParams.NUser,simParams.NPaths); 
rng('default')

BER_MP = zeros(nSNRLevels,1); 
BER_DFT = zeros(nSNRLevels,1);
BER_MMSE = zeros(nSNRLevels,1);
BER_theoretical = zeros(nSNRLevels,1);
nQPSKSym = 14;
%% >>>>>>>>>>>>>>>> MAIN SIMULATION LOOP <<<<<<<<<<<<<<<<<<<
for snrdB=SNRmin:step:SNRmax
    snrPower = 10^(snrdB/10); % power
    txSNRAmp = 10^(snrdB/20);
    tmp = 0;
    fprintf(">>>>>>>>>>>>>>>>>>> Current SNR %d >>>>>>>>>>>>>>>>>>>>> \n", snrdB)
    for iMCTrial = 1:totalMCTrials
        % Get maximum channel delay by transmitting random signal
        toShuffle = ~mod(iMCTrial-1,1000);
        if toShuffle 
            ichannel = randi([1 3]);
            chanFileName = fullfile("ScatMIMOChanModels",sprintf("chanModel%d.mat",ichannel));
            chanModel = load(chanFileName);
            channel = chanModel.channel;
            simParams.scatPos = channel.ScattererPosition;
            simParams.posRx = channel.ReceiveArrayPosition;
            [~,~,tau] = channel(complex(randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx), ...
                randn(ofdmInfo.SampleRate*1e-3,simParams.NumTx)));
            maxChDelay = ceil(max(tau)*ofdmInfo.SampleRate);
            % Compute True AoA based on ScatPos and posRx
            [~,AoA_True] = rangeangle(simParams.scatPos,simParams.posRx,refax);
            AoA_True(1,:) = sort(AoA_True(1,:),'ascend');
            n = 0:simParams.NumRx-1;eta=pi;
            Amk_true = zeros(simParams.NumRx,simParams.NPaths);
            for npath=1:simParams.NPaths
                Amk_true(:,npath) = exp(-1i.*n'*eta*cosd(AoA_True(1,npath)))./sqrt(simParams.NumRx);
            end
            
            nSlot = 1;
            carrier.NSlot = nSlot; % current slot number
            txGrid = nrResourceGrid(carrier,simParams.NumTx);
            qamSymbols = nrSymbolModulate(randi([0,1],numel(txGrid)*2,1),'QPSK');
            txAmp = txSNRAmp;
            txGrid(:) = txAmp*qamSymbols;
            
            
            %% === Send the Waveform through the Channel ===
            % OFDM Modulation
            [txWaveform,~] = nrOFDMModulate(carrier,txGrid);
            % Append zeros to the transmit waveform to account for channel delay
            txWaveform = [txWaveform; zeros(maxChDelay,simParams.NumTx)];
            % Pass the waveform through the channel
            [fadWave,pathGains,~] = channel(txWaveform);
            pg = permute(pathGains,[2,3,1]);
            [~,pgSortedInd] = sort(mean(abs(pg),1),2);
            pgSorted = pg(:,pgSortedInd);
            Dmk_true = diag(mean(pgSorted,1));
            hmk_true = sum(Amk_true * Dmk_true,2);
            pathLoss = pgSorted(:,1);
            % Estimate timing offset
            %% >>>>>>>>>>>>>>>>>> Channel Estimation Start >>>>>>>>>>>>>>>>>>>>
            offset = nrTimingEstimate(carrier,fadWave,txGrid);
            if offset > maxChDelay
                offset = 0;
            end
        end
        dist_CU_User = computeDist(channel);

        % Correct timing offset
        syncTdWaveform = fadWave(1+offset:end,:);

        % % Perform OFDM demodulation
        sigGrid = nrOFDMDemodulate(carrier,syncTdWaveform);
        N0freq = sqrt(1/(simParams.NumRx*ofdmInfo.Nfft*2));
        noiseGrid = zeros(size(sigGrid));
        noiseGrid(:) = N0freq*(randn(size(noiseGrid)) + 1i*randn(size(noiseGrid)));

        sigPower = rms(sigGrid(:)).^2;
        noisePower = rms(noiseGrid(:)).^2;
        rxGrid = sigGrid + noiseGrid;

        freqIdx = 1;
        Eb = norm(sigGrid(1,1:nQPSKSym))^2 / nQPSKSym;
        N0 = norm(noiseGrid(1,1:nQPSKSym,1))^2 / nQPSKSym;       
        Xpilot = txGrid(freqIdx,firstSymb:lastSymb);
        Npilot = reshape(noiseGrid(freqIdx,firstSymb:lastSymb,:),[],NumPilotSymb,simParams.NumRx); % pilot recv: NPilotSymb X NRx
        Ypilot = reshape(rxGrid(freqIdx,firstSymb:lastSymb,:),[],NumPilotSymb,simParams.NumRx); % pilot recv: NPilotSymb X NRx
        Ypilot = permute(Ypilot,[3,2,1]); % NumRx X NumPilotSym
        Npilot = permute(Npilot,[3,2,1]); % NumRx X NumPilotSym
        h_LinMMSE = h_MMSE_CE(Ypilot,Xpilot,Npilot);


        %% ========= Perform other Channel Estimation Algo. ===============
        ySampled = Ypilot * Xpilot'/norm(Xpilot)^2;
        AOAs_matpencil = matpencil_aoa(ySampled,simParams.NPaths);
        AOAs_dft = dft_aoa(ySampled,simParams.NumRx,simParams.NPaths);
        
		for path=1:simParams.NPaths
			SqErr_snr_AoA_MP_individual(path,snrIdx) = SqErr_snr_AoA_MP_individual(path,snrIdx) + abs(AoA_True(1,path) - AOAs_matpencil(path)).^2;
			SqErr_snr_AoA_DFT_individual(path,snrIdx) = SqErr_snr_AoA_DFT_individual(path,snrIdx) + abs(AoA_True(1,path) - AOAs_dft(path)).^2;
		end
		
		%% Channel Amplitude Estimation
		Amk_hat_matpencil = zeros(simParams.NumRx,simParams.NPaths);
		Amk_hat_dft = zeros(simParams.NumRx,simParams.NPaths);
		for path=1:simParams.NPaths
		   Amk_hat_matpencil(:,path) = exp(-1i.*n*eta*cosd(AOAs_matpencil(path)))./sqrt(simParams.NumRx);
		   Amk_hat_dft(:,path) = exp(-1i.*n*eta*cosd(AOAs_dft(path))) ./ sqrt(simParams.NumRx);
		end
		Dmk_hat_mp = pinv(Amk_hat_matpencil'*Amk_hat_matpencil)*Amk_hat_matpencil'*ySampled;
		Rd_hat_mp = (Dmk_hat_mp)*(Dmk_hat_mp');
		Bmk_hat_mp = sort(diag(sqrt(Rd_hat_mp)));
		Gmk_hat_mp = Amk_hat_matpencil * Bmk_hat_mp / sqrt(simParams.NPaths);
		
		Dmk_hat_dft = pinv(Amk_hat_dft'*Amk_hat_dft)*Amk_hat_dft'*ySampled;
		Rd_hat_dft = (Dmk_hat_dft)*(Dmk_hat_dft');
		Bmk_hat_dft = sort(diag(sqrt(Rd_hat_dft)));
		Gmk_hat_dft = Amk_hat_dft * Bmk_hat_dft / sqrt(simParams.NPaths);

        SqErr_snr_Gamma_MP(snrIdx) = SqErr_snr_Gamma_MP(snrIdx) + sum(mean(abs(hmk_true - Gmk_hat_mp).^2,1));
		SqErr_snr_Gamma_DFT(snrIdx) = SqErr_snr_Gamma_DFT(snrIdx) + sum(mean(abs(hmk_true - Gmk_hat_dft).^2,1));
		SqErr_snr_Hmk_LinearMMSE(snrIdx) = SqErr_snr_Hmk_LinearMMSE(snrIdx) + sum(mean(abs(hmk_true - h_LinMMSE).^2));
		%% <<<<<<<<<<<<<<<< Channel Estimation End <<<<<<<<<<<<<<<<<<<

        %% >>>>>>>>>> BER Estimation Start >>>>>>>>>>>>
        qpskSymRecv = reshape(rxGrid(freqIdx,1:nQPSKSym,:),[],nQPSKSym,simParams.NumRx); % pilot recv: NPilotSymb X NRx
        qpskSymRecv = permute(qpskSymRecv,[3,2,1]); % NumRx X NumPilotSym
        symEnc = txGrid(freqIdx,1:nQPSKSym) ./ txAmp;
        symDecMP = qpskDemap(pinv(Gmk_hat_mp)*qpskSymRecv,qpskRef);
        symDecDFT = qpskDemap(pinv(Gmk_hat_dft)*qpskSymRecv,qpskRef);
        symDecLMMSE = qpskDemap(pinv(h_LinMMSE)*qpskSymRecv,qpskRef);
        
        berMP = nnz(~(abs(symEnc - symDecMP) < tolFloating));
        BER_MP(snrIdx) = BER_MP(snrIdx) + berMP; 
        berDFT = nnz(~(abs(symEnc - symDecDFT) < tolFloating));
        BER_DFT(snrIdx) = BER_DFT(snrIdx) + berDFT; 
        berMMSE = nnz(~(abs(symEnc - symDecLMMSE) < tolFloating));
        BER_MMSE(snrIdx) = BER_MMSE(snrIdx) + berMMSE;

        
        BER_theoretical(snrIdx) = BER_theoretical(snrIdx) + computeIdealBER(Eb,N0);
        %% <<<<<<<<<< BER Estimation End <<<<<<<<<<<<<<

    end
    snrIdx = snrIdx + 1;
end


RMSE_snr_AoA_MP_individual = zeros(simParams.NPaths,nSNRLevels);
RMSE_snr_AoA_DFT_individual = zeros(simParams.NPaths,nSNRLevels);
RMSE_snr_Gamma_MP = sqrt(SqErr_snr_Gamma_MP / totalMCTrials);
RMSE_snr_Gamma_DFT = sqrt(vpa(SqErr_snr_Gamma_DFT) / totalMCTrials);
RMSE_snr_Hmk_LinearMMSE = sqrt(vpa(SqErr_snr_Hmk_LinearMMSE) / totalMCTrials);
for path=1:simParams.NPaths
	RMSE_snr_AoA_MP_individual(path,:) = SqErr_snr_AoA_MP_individual(path,:) / totalMCTrials;
	RMSE_snr_AoA_DFT_individual(path,:) = SqErr_snr_AoA_DFT_individual(path,:) / totalMCTrials;
end
RMSE_snr_AoA_MP_individual = sqrt(RMSE_snr_AoA_MP_individual / simParams.NPaths);
RMSE_snr_AoA_DFT_individual = sqrt(RMSE_snr_AoA_DFT_individual / simParams.NPaths);
job=string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss"));
% job=getenv('SLURM_JOB_ID');
% for np = 1:simParams.NPaths
%     fig1=figure;
%     hold on
%     grid on
%     plot(SNRmin:step:SNRmax,log10(RMSE_snr_AoA_MP_individual(np,:)),'Color','#0099ff','DisplayName','AoA Matrix-Pencil');
%     plot(SNRmin:step:SNRmax,log10(RMSE_snr_AoA_DFT_individual(np,:)),'Color','red','DisplayName','AoA DFT');
%     xlabel('SNR (dB)')
%     ylabel('RMSE of AoA Estimation (dB-Scale)')
%     title(sprintf('RMSE of AoA Estimation -- No.%d Path NAnt=%d',np,simParams.NumRx))
%     legend show
%     hold off
%     pngfile=sprintf('RMSE_snr_AoA_path_%d_j%s',np,job);
%     print(fig1,pngfile,'-dpng')
% end

fig2=figure;
hold on; grid on
plot(SNRmin:step:SNRmax,log10(RMSE_snr_Gamma_MP),'Color','#0099ff','DisplayName','Hmk Matrix Pencil');
plot(SNRmin:step:SNRmax,log10(RMSE_snr_Gamma_DFT),'Color','red','DisplayName','Hmk DFT');
plot(SNRmin:step:SNRmax,log10(RMSE_snr_Hmk_LinearMMSE),'DisplayName','Hmk Linear-MMSE');
xlabel('SNR (dB)')
ylabel('RMSE of Channel Response Estimation (dB scale)')
title(sprintf('RMSE of Channel Response Estimation: NumPilot=%d NAnt=%d',NumPilotSymb,simParams.NumRx));
legend show
hold off
pngfile=sprintf('RMSE_snr_CSI_Est_j%s',job);
print(fig2,pngfile,'-dpng')

BER_MP = BER_MP ./ totalMCTrials; 
BER_DFT = BER_DFT ./ totalMCTrials;
BER_MMSE = BER_MMSE ./ totalMCTrials;
BER_theoretical = BER_theoretical ./ totalMCTrials;
fig3=figure;
hold on; grid on
plot(SNRmin:step:SNRmax,log10(BER_MP),'Color','#0099ff','DisplayName','BER Matrix Pencil');
plot(SNRmin:step:SNRmax,log10(BER_DFT),'Color','red','DisplayName','BER DFT');
plot(SNRmin:step:SNRmax,log10(BER_MMSE),'DisplayName','BER Linear-MMSE');
plot(SNRmin:step:SNRmax,log10(BER_theoretical),'DisplayName','BER theoretical');
xlabel('SNR (dB)')
ylabel('BER (dB scale)')
title(sprintf('BER vs. Transmitter Power: NumPilot=%d NAnt=%d',NumPilotSymb,simParams.NumRx));
legend show
hold off
pngfile=sprintf('BER_TxPower%s',job);
print(fig3,pngfile,'-dpng')

%% Plot the Scattering MIMO Scenario
% rxArrayStv = phased.SteeringVector('SensorArray',simParams.rxArray,'PropagationSpeed',simParams.c);
% [~,scatRxAng] = rangeangle(simParams.scatPos(:,1),simParams.posRx,refax);
% azRxBeamWidth = 30; % In degrees
% elRxBeamWidth = 30; % In degrees
% rxAng = getInitialBeamDir(scatRxAng,azRxBeamWidth,elRxBeamWidth);
% wR = rxArrayStv(simParams.fc,rxAng);
% Configure the MIMO scene parameters
% sceneParams.TxArray = simParams.txArray;
% sceneParams.RxArray = simParams.rxArray;
% sceneParams.TxArrayPos = simParams.posTx;
% sceneParams.RxArrayPos = simParams.posRx;
% sceneParams.ScatterersPos = simParams.scatPos;
% sceneParams.Lambda = simParams.lambda;
% sceneParams.ArrayScaling = 100;   % To enlarge antenna arrays in the plot
% sceneParams.MaxTxBeamLength = 10; % Maximum length of transmit beams in the plot
% sceneParams.MaxRxBeamLength = 25; % Maximum length of receive beam in the plot
% wT=1;
% hPlotSpatialMIMOScene(sceneParams,wT,wR);
% view(2)


% hSRSGrid(carrier,srs,1,true); % Create and display a single-slot resource grid containing SRS
% title(['Resource Grid Containing SRS. NRRC = ' num2str(srs.NRRC)]);
% hSRSAnnotations(carrier,srs);


%%======== Channel Estimation ========
function h_LinMMSE = h_MMSE_CE(y,x,n)
% y                 = Frequency-domain received signal
%                   NumRx X NumSym
% d                 = Large-scale fading parameter
%                   NumRx X NumPaths
% x                 = pilot symbol
% noise
[NumRx,NumPilotSymb] = size(y);
yExtracted = y * x'  /norm(x)^2;
noiseSigVar = norm(n)^2 /norm(x)^2.* eye(NumRx);
Dmk = y * y' / NumPilotSymb;
Dmk_LS = yExtracted*yExtracted';
h_LinMMSE = Dmk * pinv(Dmk +noiseSigVar) * yExtracted;
% fprintf("==== y and yExtracted ====\n")
% disp(norm(y))
disp(norm(Dmk_LS * pinv(Dmk + noiseSigVar)))
disp(norm(yExtracted))
end



%% AoA Estimation Related
%% ===== Extended DFT + Angle-Rotation =====
function [AoA_Est_DFT] = dft_aoa(ymk_Sampled,N,L)
% Input:    yt, NumRx, N_MultiPath
% Output:   AOA_estimated, beta_estimated
FN = dftmtx(N)/sqrt(N);
Ndft_points = 5000; %% can choose whatever # you want

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

    AoA_Est_DFT(l) = acosd(theta);
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
        if I(l) >= (threshold + 2)
            Q(pl) = I(l)-1;
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
    zR = real(z_tmp);

    theta_tmp = acosd(imag(log(z_tmp)) ./ pi);
    eR = real(exp(1i*pi*cosd(theta_tmp)));
    if sign(eR) == sign(zR)
        AOAs(pl) = theta_tmp;
    end
end
AOAs = sort(AOAs);
end



    

    
    