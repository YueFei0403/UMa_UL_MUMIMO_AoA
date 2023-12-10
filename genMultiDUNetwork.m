clc;
clear all;
close all;


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

simParams.numDU = 4; % intensity of DUs -- per circle of radius 500m
simParams.numRxMIMO = simParams.NumRx*simParams.numDU; % Number of receivers at the user end
simParams.posRx = [0;0;0];
% Configure Scatterers
simParams.refax = [[1;0;0] [0;1;0] [0;0;0]];

simParams.serveRadius = [20 30 40 50 100 200 500];
simParams.numServeRadius = length(simParams.serveRadius);

% Configure the transmit and receive antenna elements
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

simParams.NChannelModel = 30;
simParams.folderName = "MultiDUChannelModels";
if ~exist(simParams.folderName,'dir')
    fprintf("MultiDU-Channel Model Folder Doesn't Exist\n" + ...
        "Creating New Folder >> \n");
    mkdir(simParams.folderName);
    genMultiDUChan(simParams,ofdmInfo);
else
    fprintf("== Channel Model Folder Exists ==\n");
    rmdir(simParams.folderName,'s');
    fprintf("Old Folder deleted << \n");
    mkdir(simParams.folderName);
    fprintf("New Folder created >> \n");
    genMultiDUChan(simParams,ofdmInfo)
end


function genMultiDUChan(simParams,ofdmInfo)
x0=0;y0=0;
for radiusIdx = 1:simParams.numServeRadius
    xcart = simParams.serveRadius(radiusIdx) * rand(simParams.numDU,1);
    ycart = simParams.serveRadius(radiusIdx) * rand(simParams.numDU,1);
    xcart = xcart + x0;
    ycart = ycart + y0;
    distX1 = ceil(min(xcart)); distY1 = ceil(min(ycart));
    distX2 = floor(max(xcart)); distY2 = floor(max(ycart));
    xcoord = randi([distX1,distX2],simParams.NPaths);
    ycoord = randi([distY1,distY2],simParams.NPaths);

    for ichannel = 1:simParams.NChannelModel
        for DUIdx = 1:simParams.numDU
            txPosX = xcart(DUIdx);
            txPosY = ycart(DUIdx);
            simParams.posTx = [txPosX;txPosY;0];
            for npath = 1:simParams.NPaths
                xScat = xcoord(npath); yScat=ycoord(npath);
                simParams.scatPos(:,npath) =[xScat;yScat;0];
            end
            channel = genScatteringMIMOChannel(simParams,ofdmInfo);
            channelFile = fullfile(simParams.folderName,sprintf("R%d-Chan%d-DU%d.mat", ...
                simParams.serveRadius(radiusIdx),ichannel,DUIdx));
            save(channelFile,"channel");
        end
    end
end
end



function channel = genScatteringMIMOChannel(simParams,ofdmInfo)
% === Configure the Channel ===
channel = phased.ScatteringMIMOChannel;
channel.PropagationSpeed = simParams.c;
channel.CarrierFrequency = simParams.fc;
channel.SimulateDirectPath = false;
channel.Polarization = 'None';
channel.SpecifyAtmosphere = false;
channel.ChannelResponseOutputPort = true;
channel.SampleRate = ofdmInfo.SampleRate;
channel.MaximumDelaySource = 'Property';
channel.MaximumDelay = 3e-6; % empirical data

% Configure Scatterers
scatCoeff = ones(1,simParams.NPaths);
channel.ScattererSpecificationSource = 'Property'; %  All scatterer velocities are zero.
channel.ScattererPosition = simParams.scatPos; 
channel.ScattererCoefficient = scatCoeff;

% Configure transmit array parameters
channel.TransmitArrayMotionSource = 'Property';
channel.TransmitArray = simParams.txArray;
channel.TransmitArrayPosition = simParams.posTx;
channel.TransmitArrayOrientationAxes = [[1;0;0] [0;1;0] [0;0;-1]];

% Configure receive array parameters
channel.ReceiveArrayMotionSource = 'Property';
channel.ReceiveArray = simParams.rxArray;
channel.ReceiveArrayPosition = simParams.posRx;
channel.ReceiveArrayOrientationAxes = [[1;0;0] [0;1;0] [0;0;-1]];
end



%% Simulate the Locations of DU following Poisson Distribution
% % Suppose user at origin
% % 
% numServeRadius = length(serveRadius);
% areaRef = pi*500^2;
% x0=0;y0=0;
% figure
% hold on
% for radiusIdx = 1:numServeRadius
%     currRadius = serveRadius(radiusIdx);
%     currArea = pi*currRadius^2 / areaRef;
%     currNumDU = floor(poissrnd(currArea*meanNumDU)); %round down
%     fprintf("\n--- Current number of DUs=%d ----\n",currNumDU)
%     theta = 2*pi*(rand(currNumDU,1)); % angular coordinates
%     rho = currRadius*sqrt(rand(currNumDU,1));  % radial coordinates
%     [xcart,ycart] = pol2cart(theta,rho);
%     xcart = xcart + x0;
%     ycart = ycart + y0;
%     %Plotting
% 
%     [nR,dR] = rat(currRadius);
%     scatter(xcart,ycart,'DisplayName',sprintf("r=%d/%d",nR,dR));
% 
% end
% xlabel('x');ylabel('y');
% axis square;
% legend show
% hold off









