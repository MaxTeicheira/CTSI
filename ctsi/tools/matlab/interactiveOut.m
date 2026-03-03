% This script parses interactiveOut.txt
clear; clc;

data = load('../output/interactiveOut.txt');
baseInd = 0;
timeLen = data(1);
baseInd = baseInd+1;

% vector<double> &timeVec
timeVec = data(baseInd+1:baseInd+timeLen);
baseInd = baseInd+timeLen;

% vector<int> anCollectorID
numAnCollectedCharge = data(baseInd+1);
baseInd = baseInd+1;
anCollectorID = data(baseInd+1:baseInd+numAnCollectedCharge);
baseInd = baseInd+numAnCollectedCharge;

% vector<int> caCollectorID
numCaCollectedCharge = data(baseInd+1);
baseInd = baseInd+1;
caCollectorID = data(baseInd+1:baseInd+numCaCollectedCharge);
baseInd = baseInd+numCaCollectedCharge;

%==========================================================================
% Anode section
disp(['Importing:']);
numAnTrig = data(baseInd+1);
baseInd = baseInd+1;
% vector<int> &anodeID
disp(['anodeID']);
anodeID = data(baseInd+1:baseInd+numAnTrig);
baseInd = baseInd+numAnTrig;
% vector<double> &anEnergy
disp(['anEnergy']);
anEnergy = data(baseInd+1:baseInd+numAnTrig);
baseInd = baseInd+numAnTrig;
% vector<double> &anTriggerTime
disp(['anTriggerTime']);
anTriggerTime = data(baseInd+1:baseInd+numAnTrig);
baseInd = baseInd+numAnTrig;
% vector<double> &noisyAnEnergy
disp(['noisyAnEnergy']);
noisyAnEnergy = data(baseInd+1:baseInd+numAnTrig);
baseInd = baseInd+numAnTrig;
% vector<double> &noisyAnTriggerTime
disp(['noisyAnTriggerTime']);
noisyAnTriggerTime = data(baseInd+1:baseInd+numAnTrig);
baseInd = baseInd+numAnTrig;

% vector<vector<double> > &anVsTime
disp(['anVsTime']);
numAnCollectedCharge = data(baseInd+1);
baseInd = baseInd+1;
anVsTime = cell(1, numAnCollectedCharge);
for i = 1:numAnCollectedCharge
    anVsTime{i} = data(baseInd+1:baseInd+timeLen);
    baseInd = baseInd+timeLen;
end

%==========================================================================
% Cathode section
numCaTrig = data(baseInd+1);
baseInd = baseInd+1;
% vector<int> &cathodeID
disp(['cathodeID']);
cathodeID = data(baseInd+1:baseInd+numCaTrig);
baseInd = baseInd+numCaTrig;
% vector<double> &caEnergy
disp(['caEnergy']);
caEnergy = data(baseInd+1:baseInd+numCaTrig);
baseInd = baseInd+numCaTrig;
% vector<double> &caTriggerTime
disp(['caTriggerTime']);
caTriggerTime = data(baseInd+1:baseInd+numCaTrig);
baseInd = baseInd+numCaTrig;
% vector<double> &noisyCaEnergy
noisyCaEnergy = data(baseInd+1:baseInd+numCaTrig);
baseInd = baseInd+numCaTrig;
% vector<double> &noisyCaTriggerTime
noisyCaTriggerTime = data(baseInd+1:baseInd+numCaTrig);
baseInd = baseInd+numCaTrig;

% vector<vector<double> > &caVsTime
disp(['caVsTime']);
numCaCollectedCharge = data(baseInd+1);
baseInd = baseInd+1;
caVsTime = cell(1, numCaCollectedCharge);
for i = 1:numCaCollectedCharge
    caVsTime{i} = data(baseInd+1:baseInd+timeLen);
    baseInd = baseInd+timeLen;
end

%==========================================================================
% vector<vector<struct trailLog> > &eventTrails
disp(['eventTrails']);
numFieldsInTrail = data(baseInd+1);
baseInd = baseInd+1;
numIntrxn = data(baseInd+1);
baseInd = baseInd+1;
eventTrails = cell(1, numIntrxn);
% Loop over all interactions
for i = 1:numIntrxn
    disp(['Interaction ' num2str(i)]);
    % Loop over all charge elements
    numChargeElems = data(baseInd+1);
    baseInd = baseInd+1;
    for j = 1:numChargeElems
		% xPosE
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.xPosE = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% yPosE
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.yPosE = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% zPosE
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.zPosE = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% qETrapped
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.qETrapped = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% qEMobile
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.qEMobile = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% xPosH
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.xPosH = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% yPosH
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.yPosH = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% zPosH
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.zPosH = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% qHTrapped
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.qHTrapped = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
		% qHMobile
        len = data(baseInd+1);
        baseInd = baseInd+1;
        eventTrails{i}{j}.qHMobile = data(baseInd+1:baseInd+len);
        baseInd = baseInd+len;
    end
end

% vector<vector<vector<struct phiLog> > > &eventAnPhi
disp(['eventAnPhi']);
numAnodes = data(baseInd+1);
baseInd = baseInd+1;
% Loop over anodes that collected charge
for i = 1:numAnodes
    % Loop over interactions
    numIntrxn = data(baseInd+1);
    baseInd = baseInd+1;
    for j = 1:numIntrxn
        % Loop over charge elements
        numChargeElems = data(baseInd+1);
        baseInd = baseInd+1;
        for h = 1:numChargeElems
            % EPhi
            len = data(baseInd+1);
            baseInd = baseInd+1;
            eventAnPhi{i}{j}{h}.EPhi = data(baseInd+1:baseInd+len);
            baseInd = baseInd+len;
            % HPhi
            len = data(baseInd+1);
            baseInd = baseInd+1;
            eventAnPhi{i}{j}{h}.HPhi = data(baseInd+1:baseInd+len);
            baseInd = baseInd+len;
        end
    end
end

% vector<vector<vector<struct phiLog> > > &eventCaPhiIn
disp(['eventCaPhi']);
numCathodes = data(baseInd+1);
baseInd = baseInd+1;
% Loop over cathodes that collected charge
for i = 1:numCathodes
    % Loop over interactions
    numIntrxn = data(baseInd+1);
    baseInd = baseInd+1;
    for j = 1:numIntrxn
        % Loop over charge elements
        numChargeElems = data(baseInd+1);
        baseInd = baseInd+1;
        for h = 1:numChargeElems
            % EPhi
            len = data(baseInd+1);
            baseInd = baseInd+1;
            eventCaPhi{i}{j}{h}.EPhi = data(baseInd+1:baseInd+len);
            baseInd = baseInd+len;
            % HPhi
            len = data(baseInd+1);
            baseInd = baseInd+1;
            eventCaPhi{i}{j}{h}.HPhi = data(baseInd+1:baseInd+len);
            baseInd = baseInd+len;
        end
    end
end

% vector<vector<int> > &trailSizeE
disp(['trailSizeE']);
numIntrxn = data(baseInd+1);
baseInd = baseInd+1;
for i = 1:numIntrxn
    % Loop over interactions
    numChargeElem = data(baseInd+1);
    baseInd = baseInd+1;
    for j = 1:numChargeElem
        trailSizeE{i}{j} = data(baseInd+1);
        baseInd = baseInd+1;
    end
end

% vector<vector<int> > &trailSizeH
disp(['trailSizeH']);
numIntrxn = data(baseInd+1);
baseInd = baseInd+1;
for i = 1:numIntrxn
    % Loop over interactions
    numChargeElem = data(baseInd+1);
    baseInd = baseInd+1;
    for j = 1:numChargeElem
        trailSizeH{i}{j} = data(baseInd+1);
        baseInd = baseInd+1;
    end
end

% vector<double> &timeVecPreamp
disp(['timeVecPreamp']);
preampLen = data(baseInd+1);
baseInd = baseInd+1;
timeVecPreamp = data(baseInd+1:baseInd+preampLen);
baseInd = baseInd+preampLen;

% vector<vector<double> > &anVsTimeReg
% Loop over all triggered anodes
disp(['anVsTimeReg']);
for i = 1:numAnCollectedCharge
	anVsTimeReg{i} = data(baseInd+1:baseInd+preampLen);
	baseInd = baseInd+preampLen;
end

% vector<vector<double> > &caVsTimeReg
% Loop over all triggered cathodes
disp(['caVsTimeReg']);
for i = 1:numCaCollectedCharge
    caVsTimeReg{i} = data(baseInd+1:baseInd+preampLen);
	baseInd = baseInd+preampLen;
end

% vector<vector<double> > &noisyAnVsTimeReg
disp(['noisyAnVsTimeReg']);
for i = 1:numAnCollectedCharge
	noisyAnVsTimeReg{i} = data(baseInd+1:baseInd+preampLen);
	baseInd = baseInd+preampLen;
end

% vector<vector<double> > &noisyCaVsTimeReg
disp(['noisyCaVsTimeReg']);
for i = 1:numCaCollectedCharge
	noisyCaVsTimeReg{i} = data(baseInd+1:baseInd+preampLen);
	baseInd = baseInd+preampLen;
end

% vector<vector<double> > &anVsTimePreamp
disp(['anVsTimePreamp']);
for i = 1:numAnCollectedCharge
    anVsTimePreamp{i} = data(baseInd+1:baseInd+preampLen);
    baseInd = baseInd+preampLen;
end
			
% vector<vector<double> > &caVsTimePreamp
disp(['caVsTimePreamp']);
for i = 1:numCaCollectedCharge
    caVsTimePreamp{i} = data(baseInd+1:baseInd+preampLen);
    baseInd = baseInd+preampLen;
end

% vector<vector<double> > &noisyAnVsTimePreamp
disp(['noisyAnVsTimePreamp']);
for i = 1:numAnCollectedCharge
    noisyAnVsTimePreamp{i} = data(baseInd+1:baseInd+preampLen);
    baseInd = baseInd+preampLen;
end

% vector<vector<double> > &noisyCaVsTimePreamp
disp(['noisyCaVsTimePreamp']);
for i = 1:numCaCollectedCharge
    noisyCaVsTimePreamp{i} = data(baseInd+1:baseInd+preampLen);
    baseInd = baseInd+preampLen;
end

clear i j h len;


%% plot anode current vs time

figure
plot(timeVec,anVsTime{1},'b')
xlim([0,2000e-9])

% figure
% plot(timeVecPreamp,anVsTimeReg{1},'b')
% xlim([0,2000e-9])
% 
% figure
% plot(timeVecPreamp,anVsTimePreamp{1},'b')
% xlim([0,2000e-9])


figure
plot(timeVec,caVsTime{1})
xlim([0,2000e-9])

%% Subplot
subplot(2, 1, 1);
plot(timeVec, anVsTime{1}, 'b');
% xlim([0, 2000e-9]);
title('Anode waveform', 'FontSize', 20);
xlabel('Time (s)', 'FontSize', 20);
ylabel('Current', 'FontSize', 20);
set(gca, 'FontSize', 20);

subplot(2, 1, 2);
plot(timeVec, caVsTime{1});
% xlim([0, 2000e-9]);
title('Cathode waveform', 'FontSize', 20);
xlabel('Time (s)', 'FontSize', 20);
ylabel('Current', 'FontSize', 20);
set(gca, 'FontSize', 20);