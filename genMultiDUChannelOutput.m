function [AoAs_True_MultiDU,Amk_True_MultiDU,sigGridMultiDU,noiseGridMultiDU,txGrid,txAmp] = genMultiDUChannelOutput( ...
    ofdmInfo,simParams,carrier,snrdB,serveRadius.ichannel)
refax = [[1;0;0] [0;1;0] [0;0;0]];

% Transmitter Setup
txAmp = 2000*10^(snrdB/20); % 30dBm
txGrid = nrResourceGrid(carrier,simParams.NumTx);
qamSymbols = nrSymbolModulate(randi([0,1],numel(txGrid)*2,1),'QPSK');
txGrid(:) = txAmp.*qamSymbols;

% Multi-DU Setup
sigGridMultiDU = [];
AoAs_True_MultiDU = [];
Amk_True_MultiDU = [];
rxArrayStv = phased.SteeringVector('SensorArray',simParams.rxArray,'PropagationSpeed',simParams.c);
for DUIdx = 1:simParams.NumDU
        chanFileName = fullfile("MultiDUChannelModels",sprintf("R%d-Chan%d-DU%d.mat", ...
                        serveRadius,ichannel,DUIdx));
        file = java.io.File(chanFileName);
        fullpath = char(file.getAbsolutePath());
        try
            chanModel = load(fullpath);
        catch
            chanFileName = fullfile("MultiDUChannelModels",sprintf("R%d-Chan%d-DU%d.mat", ...
                        serveRadius,ichannel,DUIdx));
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
        
        
        %% === Send the Waveform through the Channel ===
        % OFDM Modulation
        [txWaveform,~] = nrOFDMModulate(carrier,txGrid);
        % Append zeros to the transmit waveform to account for channel delay
        txWaveform = [txWaveform; zeros(maxChDelay,simParams.NumTx)];
        % Pass the waveform through the channel
        [fadWave,~,~] = channel(txWaveform);
        
        % Estimate timing offset
        %% >>>>>>>>>>>>>>>>>> Channel Estimation Start >>>>>>>>>>>>>>>>>>>>
        offset = nrTimingEstimate(carrier,fadWave,txGrid);
        if offset > maxChDelay
            offset = 0;
        end
        
        % Receiver Setup
        
        % Compute True AoA based on ScatPos and posRx
        [~,AoAs_True] = rangeangle(simParams.scatPos,simParams.posRx,refax);
        AoAs_True(1,:) = sort(AoAs_True(1,:),'ascend');
        Amk_True = rxArrayStv(simParams.fc,AoAs_True(1,:));
        
        
        % Correct timing offset
        fadWave = fadWave(1+offset:end,:);
        % Perform OFDM demodulation
        sigGrid = nrOFDMDemodulate(carrier,fadWave);           
       
        AoAs_True_MultiDU = cat(2,AoAs_True_MultiDU,AoAs_True);
        Amk_True_MultiDU = cat(2,Amk_True_MultiDU,Amk_True);
        sigGridMultiDU = cat(3,sigGridMultiDU,sigGrid);
end
noiseGridMultiDU = 1/sqrt(2).*complex(randn(size(sigGridMultiDU)),randn(size(sigGridMultiDU)));
end
