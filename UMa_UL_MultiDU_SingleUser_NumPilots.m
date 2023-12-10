clc;
clear all;
close all;
digits(100)
warning('off')

%% === Configure Base-station transmitters  ===
carrier = nrCarrierConfig;
carrier.NSizeGrid = 1; % Bandwidth in RB
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'Normal';
ofdmInfo = nrOFDMInfo(carrier);
K = carrier.NSizeGrid * 12;

% === Configure Transmit and Receive Antenna Arrays ===
simParams.fc = 6e9; % freqRange = 'FR1';
simParams.c = physconst('LightSpeed');

simParams.lambda = simParams.c/simParams.fc;
simParams.NumTx = 1; 
simParams.NumRx = 8;
simParams.NPaths = 7;
simParams.NumDU = 4; % intensity of DUs -- per circle of radius 500m
simParams.NumRxMultiDU = simParams.NumRx*simParams.NumDU; % Number of receivers at the user end
simParams.posRx = [0;0;0];
% Configure Scatterers
simParams.refax = [[1;0;0] [0;1;0] [0;0;0]];

simParams.serveRadius = [20 30 40 50 100 200 500];
simParams.numServeRadius = length(simParams.serveRadius);
simParams.NChannelModel = 30;
simParams.folderName = "MultiDUChannelModels";

% Configure the transmit and receive antenna elements for each pair of single DU
% and single user
simParams.txAntenna = phased.IsotropicAntennaElement;            % To avoid transmission beyond +/- 90
                                                                 % degrees from the broadside, baffle
                                                                 % the back of the transmit antenna
                                                                 % element by setting the BackBaffled
                                                                 % property
                                                                 % to true.
simParams.rxAntenna = phased.IsotropicAntennaElement('BackBaffled',false); % To receive the signal from 360 degrees,
                                                                 % set the BackBaffled property to false

simParams.txArray = phased.NRRectangularPanelArray('Size',[1, 1, 1, 1],'ElementSet', {simParams.txAntenna},...
            'Spacing',[0.5*simParams.lambda,0.5*simParams.lambda,3*simParams.lambda,3*simParams.lambda]);
simParams.rxArray = phased.ULA('Element',simParams.rxAntenna, ...
    'NumElements',simParams.NumRx,'ElementSpacing',0.5*simParams.lambda,'ArrayAxis','x');
rxArrayStv = phased.SteeringVector('SensorArray',simParams.rxArray,'PropagationSpeed',simParams.c);

n = 0:simParams.NumRx-1;eta=pi;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% numFrames = 1; % 10 ms frames -> 1 subframe - 10 slots - 14 symbols
K = carrier.NSizeGrid * 12;
L = carrier.SymbolsPerSlot;
simParams.NUser = 1;
simParams.NPaths = 7;
pilotLengths = [1 3 6 9 12];% dB
nPilotLengths = length(pilotLengths);

%% Hold data separately for each DU in $$_SingleDU$$ data structures
SqErr_snr_Gamma_MP_SingleDU = zeros(simParams.NumDU,nPilotLengths);
SqErr_snr_Gamma_DFT_SingleDU = zeros(simParams.NumDU,nPilotLengths);
SqErr_snr_Gamma_LinearMMSE_SingleDU = zeros(simParams.NumDU,nPilotLengths);
SqErr_snr_Gamma_True_MultiDU = zeros(simParams.NumDU,nPilotLengths); 

NMSE_snr_Gamma_MP_SingleDU = zeros(simParams.NumDU,nPilotLengths);
NMSE_snr_Gamma_DFT_SingleDU = zeros(simParams.NumDU,nPilotLengths);
NMSE_snr_Gamma_MMSE_SingleDU = zeros(simParams.NumDU,nPilotLengths);

BER_MP_SingleDU = zeros(simParams.NumDU,nPilotLengths);
BER_DFT_SingleDU = zeros(simParams.NumDU,nPilotLengths);
BER_MMSE_SingleDU = zeros(simParams.NumDU,nPilotLengths);
BER_theoretical_SingleDU = zeros(simParams.NumDU,nPilotLengths);

totalMCTrials = 1000; % Monte-Carlo trials
tolFloating = 1e-2;
Beta = 4.1290e-10; % Empirical Large-scale fading coefficient
nQPSKSymb=12;
%% >>>>>>>>>>>>>>>> MAIN SIMULATION LOOP -- Pilot Training <<<<<<<<<<<<<<<<<<<
snrdB = 10;snrPower = 10^(snrdB/10); % power
serveRadius = 100;
EbN0 = zeros(nPilotLengths,1);
for pilotIdx=1:nPilotLengths
    currPilotL = pilotLengths(pilotIdx);
    fprintf(">>>>>>>>>>>>>>>>>>> Current NumPilots=%d >>>>>>>>>>>>>>>>>>>>> \n", currPilotL)
    for iMCTrial = 1:totalMCTrials
        ichannel = mod(iMCTrial,30)+1;
        for iSlot = 0:carrier.SlotsPerFrame-1
            carrier.NSlot = iSlot; % current slot number
            %% Hold concatenated channels (8x4=32 receiver antennas) in $$_MultiDU$$ data structures
            Hmk_snr_MP_MultiDU = zeros(simParams.NumRxMultiDU,1);
            Hmk_snr_DFT_MultiDU = zeros(simParams.NumRxMultiDU,1);
            Hmk_snr_LinearMMSE_MultiDU = zeros(simParams.NumRxMultiDU,1);
            Hmk_snr_True_MultiDU = zeros(simParams.NumRxMultiDU,1);

            [AoAs_True,Amk_True,sigGridMultiDU,noiseGridMultiDU,txGrid,txAmp] = genMultiDUChannelOutput(ofdmInfo,simParams,carrier,snrdB,serveRadius,ichannel);
            freqIdx = randi([1 K]); % in case of freq-selective fading-> choose a random row
            
            %txGrid:           NumSubcarrier X NumOFDMSym X
            %                       NumTransmitterAntenna 
            %sigGridMultiDU:   NumSubcarrier X NumOFDMSym X
            %                                   (NumReceiverAntennaPerDU X
            %                                   NumDU)
            for DUIdx = 1:simParams.NumDU
                antennaRange = (DUIdx-1)*8+1:DUIdx*8;
                sigGrid = sigGridMultiDU(:,:,antennaRange);
                noiseGrid = noiseGridMultiDU(:,:,antennaRange);
                rxGrid = sigGrid + noiseGrid;
                %% ========= Perform Angle-Domain Channel Estimation Algo. ===============
                Xpilot = txGrid(freqIdx,1:currPilotL);
                Ypilot = reshape(rxGrid(freqIdx,1:currPilotL,:),[],currPilotL,simParams.NumRx);     % pilot recv: NPilotSymb X NRx
                Ypilot = permute(Ypilot,[3,2,1]); % NumRx X NumPilotSym
                Npilot = reshape(noiseGrid(freqIdx,1:currPilotL,:),[],currPilotL,simParams.NumRx);  % pilot recv: NPilotSymb X NRx
                Npilot = permute(Npilot,[3,2,1]); % NumRx X NumPilotSym
                pilotNorm = norm(Xpilot);
                
                AOAs_matpencil = zeros(1,simParams.NPaths);
                AOAs_dft = zeros(1,simParams.NPaths);
                for iSample = 1:currPilotL
                    ySampled = Ypilot(:,iSample)*Xpilot(iSample)'./norm(Xpilot(iSample))^2;
                    AOAs_matpencil = AOAs_matpencil + sort(matpencil_aoa(ySampled,simParams.NPaths),'ascend');
                    AOAs_dft = AOAs_dft + reshape(sort(dft_aoa(ySampled,simParams.NumRx,simParams.NPaths),'ascend'),1,simParams.NPaths);
                end
                AOAs_matpencil = AOAs_matpencil ./ currPilotL;
                AOAs_dft = AOAs_dft ./ currPilotL;
                
                %  Channel Amplitude Estimation
	            Amk_hat_matpencil = rxArrayStv(simParams.fc,[AOAs_matpencil;zeros(1,simParams.NPaths)]);
	            Amk_hat_dft = rxArrayStv(simParams.fc,[AOAs_dft;zeros(1,simParams.NPaths)]);
                Gmk_hat_mp = Amk_hat_matpencil* (Amk_hat_matpencil\ (Ypilot*Xpilot'/pilotNorm^2)) ;
                Gmk_hat_dft = Amk_hat_dft * (Amk_hat_dft \ (Ypilot*Xpilot'/pilotNorm^2));
                Gmk_true = Amk_True * (Amk_True \ (Ypilot*Xpilot'/pilotNorm^2));
        
                %% ========= Perform MMSE Channel Estimation    ===============
                h_LinMMSE = h_MMSE_CE(Ypilot,Xpilot,Npilot,Beta);
        
                %% ========= Perform True Channel Estimation Algo. ===============
                refInd = [freqIdx:12:freqIdx+currPilotL*12-1];
                refSym = txGrid(refInd);
                [hEstGrid,noiseEst] = nrChannelEstimate(rxGrid,txGrid); % rxGrid=12x14x8; txGrid=12x14x1
                hmk_true = reshape(hEstGrid(freqIdx,1,:),simParams.NumRx,[]);
    
                %% ========= Perform Error Estimation ===============
                sqErr_MP = mean(norm(hmk_true - Gmk_hat_mp)^2);
                sqErr_DFT = mean(norm(hmk_true - Gmk_hat_dft)^2);
                sqErr_MMSE = mean(norm(hmk_true - h_LinMMSE)^2);
                SqErr_snr_Gamma_MP_SingleDU(DUIdx,pilotIdx) = SqErr_snr_Gamma_MP_SingleDU(DUIdx,pilotIdx) + sqErr_MP;
	            SqErr_snr_Gamma_DFT_SingleDU(DUIdx,pilotIdx) = SqErr_snr_Gamma_DFT_SingleDU(DUIdx,pilotIdx) + sqErr_DFT;
	            SqErr_snr_Gamma_LinearMMSE_SingleDU(DUIdx,pilotIdx) = SqErr_snr_Gamma_LinearMMSE_SingleDU(DUIdx,pilotIdx) + sqErr_MMSE;
    
                sqNorm_Gamma_true = mean(norm(hmk_true)^2);
                nmse_MP = sqErr_MP / sqNorm_Gamma_true;
                nmse_DFT = sqErr_DFT / sqNorm_Gamma_true;
                nmse_MMSE = sqErr_MMSE / sqNorm_Gamma_true;
                NMSE_snr_Gamma_MP_SingleDU(DUIdx,pilotIdx) = NMSE_snr_Gamma_MP_SingleDU(DUIdx,pilotIdx) + nmse_MP;
                NMSE_snr_Gamma_DFT_SingleDU(DUIdx,pilotIdx) = NMSE_snr_Gamma_DFT_SingleDU(DUIdx,pilotIdx) + nmse_DFT;
                NMSE_snr_Gamma_MMSE_SingleDU(DUIdx,pilotIdx) = NMSE_snr_Gamma_MMSE_SingleDU(DUIdx,pilotIdx) + nmse_MMSE;
                
	            % <<<<<<<<<<<<<<<< Channel Estimation End <<<<<<<<<<<<<<<<<<<
            
                %% >>>>>>>>>> BER Estimation Start at Single DU >>>>>>>>>>>>
                symQPSKInd =  [freqIdx:12:freqIdx+currPilotL*12-1];
                yQPSK = reshape(rxGrid(freqIdx,1:currPilotL,:),[],currPilotL,simParams.NumRx);     % pilot recv: NPilotSymb X NRx
                yQPSK = permute(yQPSK,[3,2,1]); % NumRx X NumPilotSym
                symEnc = txGrid(freqIdx,1:currPilotL) ./ txAmp;
                symDecMP = Gmk_hat_mp'*yQPSK;
                symDecDFT = Gmk_hat_dft'*yQPSK;
                symDecMMSE = h_LinMMSE'*yQPSK;
                symDecTrue = hmk_true'*yQPSK;
            
                symEncQPSK = nrSymbolDemodulate(symEnc.','QPSK','DecisionType','Hard');
                symDemodMP = nrSymbolDemodulate(symDecMP.','QPSK','DecisionType','Hard');
                symDemodDFT = nrSymbolDemodulate(symDecDFT.','QPSK','DecisionType','Hard');
                symDemodMMSE = nrSymbolDemodulate(symDecMMSE.','QPSK','DecisionType','Hard');
                symDemodTrue = nrSymbolDemodulate(symDecTrue.','QPSK','DecisionType','Hard');
            
                berMP = biterr(symEncQPSK,symDemodMP);
                berDFT = biterr(symEncQPSK,symDemodDFT);
                berLMMSE = biterr(symEncQPSK,symDemodMMSE);
                berTrue = biterr(symEncQPSK,symDemodTrue);
                BER_MP_SingleDU(DUIdx,pilotIdx) = BER_MP_SingleDU(DUIdx,pilotIdx) + berMP/(2*currPilotL);
                BER_DFT_SingleDU(DUIdx,pilotIdx) = BER_DFT_SingleDU(DUIdx,pilotIdx) + berDFT/(2*currPilotL);
                BER_MMSE_SingleDU(DUIdx,pilotIdx) = BER_MMSE_SingleDU(DUIdx,pilotIdx) + berLMMSE/(2*currPilotL);
                BER_theoretical_SingleDU(DUIdx,pilotIdx) = BER_theoretical_SingleDU(DUIdx,pilotIdx) + berTrue/(2*currPilotL);
            end
            EbN0(pilotIdx) = EbN0(pilotIdx) + norm(Ypilot)^2 / norm(Npilot)^2; % we only count EbN0 for single DU
        end
    end
    BER_MP_SingleDU(:,pilotIdx) = BER_MP_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    BER_DFT_SingleDU(:,pilotIdx) = BER_DFT_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    BER_MMSE_SingleDU(:,pilotIdx) = BER_MMSE_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    BER_theoretical_SingleDU(:,pilotIdx) = BER_theoretical_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);

    % <<<<<<<<<< BER Estimation End <<<<<<<<<<<<<<
    NMSE_snr_Gamma_MP_SingleDU(:,pilotIdx) = NMSE_snr_Gamma_MP_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    NMSE_snr_Gamma_DFT_SingleDU(:,pilotIdx) = NMSE_snr_Gamma_DFT_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    NMSE_snr_Gamma_MMSE_SingleDU(:,pilotIdx) = NMSE_snr_Gamma_MMSE_SingleDU(:,pilotIdx) ./ (totalMCTrials*carrier.SlotsPerFrame);
    
    EbN0(pilotIdx) = EbN0(pilotIdx)  ./ (totalMCTrials*carrier.SlotsPerFrame);
end

%% >>>>>>>>>>>>>>>>>>>>> Plotting Start >>>>>>>>>>>>>>>>>>>>
job=string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss"));

for DUIdx=1:simParams.NumDU
    fig1=figure;
    hold on; 
    semilogy(pilotLengths,NMSE_snr_Gamma_MP_SingleDU(DUIdx,:),'-o','Color','#0099ff','DisplayName','Hmk Matrix Pencil');
    semilogy(pilotLengths,NMSE_snr_Gamma_DFT_SingleDU(DUIdx,:),'-^','Color','red','DisplayName','Hmk DFT');
    semilogy(pilotLengths,NMSE_snr_Gamma_MMSE_SingleDU(DUIdx,:),'-*','DisplayName','Hmk MMSE');
    set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
    grid on
    xlabel('Number of Pilots')
    ylabel('NMSE of Channel Response Estimation')
    title(sprintf('NMSE CSI: NAnt=%d ServeRadius=%d DUIdx=%d',simParams.NumRx,serveRadius,DUIdx));
    legend show
    hold off
    pngfile=sprintf('NMSE_pilot_CSI_DUIdx%dj%s',DUIdx,job);
    print(fig1,pngfile,'-dpng')
end

for DUIdx=1:simParams.NumDU
    fig2=figure;
    hold on
    semilogy(pilotLengths,BER_MP_SingleDU(DUIdx,:),'-o','DisplayName','BER Matrix-Pencil');
    semilogy(pilotLengths,BER_DFT_SingleDU(DUIdx,:),'-^','DisplayName','BER DFT');
    semilogy(pilotLengths,BER_MMSE_SingleDU(DUIdx,:),'-*','DisplayName','BER Linear-MMSE');
    semilogy(pilotLengths,BER_theoretical_SingleDU(DUIdx,:),'-square','DisplayName','BER theoretical');
    set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
    grid on
    xlabel('Number of Pilots')
    ylabel('BER')
    title(sprintf('BER vs Number of Pilots ServeRadius=%d DUIdx=%d',serveRadius,DUIdx));
    legend show     
    hold off
    pngfile=sprintf('BER_Pilot_DUIdx%dj%s',DUIdx,job);
    print(fig2,pngfile,'-dpng')
end


fig3=figure;
hold on
semilogy(pilotLengths,EbN0,'-o','DisplayName','EbN0');
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
grid on
xlabel('Number of Pilots')
ylabel('EbN0')
title(sprintf('Num. of Pilotsvs EbN0 NumDU=%d ServeRadius=%d',simParams.NumDU,serveRadius));
legend show
hold off
pngfile=sprintf('EbN0_Pilot_MultiDU%s',job);
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
[~,currPilotL] = size(y);
yExtracted = y*x'/(norm(x)^2);
zeta = Beta / Beta + ((norm(n)^2/norm(y)^2) / currPilotL);
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