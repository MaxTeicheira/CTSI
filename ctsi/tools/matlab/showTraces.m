clear; clc;
close all;
interactiveOut

fig = figure;
for i = 1:numIntrxn
    for j = 1:numChargeElems
        plot3(eventTrails{i}{j}.xPosE, eventTrails{i}{j}.yPosE, eventTrails{i}{j}.zPosE, 'b-');
        hold on;
        plot3(eventTrails{i}{j}.xPosH, eventTrails{i}{j}.yPosH, eventTrails{i}{j}.zPosH, 'r.-');
    end
end
grid on;
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
title(['Blue: Electrons, Red: Holes']);

% axis tight; set(gca,'nextplot','replacechildren');
% v = VideoWriter('rotatingPlot2deg.mp4');
% open(v);
% for az = 1:2:360
%     view(az, 30);
%     frame = getframe(fig);
%     writeVideo(v,frame);
% end
% close(v);


figure;
for i = 1:numIntrxn
    for j = 1:numChargeElems
        plot(eventTrails{i}{j}.xPosE, eventTrails{i}{j}.zPosE, '.-');
        hold on;
        plot(eventTrails{i}{j}.xPosH, eventTrails{i}{j}.zPosH, 'r.-');
    end
end
grid on;
xlabel('x (cm)');
zlabel('z (cm)');
title(['Blue: electron, red: holes']);

fig = figure;
for i = 1:numIntrxn
    for j = 1:numChargeElems
        plot(eventTrails{i}{j}.yPosE, eventTrails{i}{j}.zPosE, '.-');
        hold on;
        plot(eventTrails{i}{j}.yPosH, eventTrails{i}{j}.zPosH, 'r.-');
    end
end
grid on;
xlabel('y (cm)');
zlabel('z (cm)');
title(['Blue: electron, red: holes']);


for i = 1:numAnCollectedCharge
	figure;
	plot(timeVec, anVsTime{i}, '.-');
	title(['Anode ' num2str(anCollectorID(i))]);
end

for i = 1:numCaCollectedCharge
	figure;
	plot(timeVec, caVsTime{i}, '.-');
	title(['Cathode ' num2str(caCollectorID(i))]);
end

for i = 1:numAnCollectedCharge
	figure;
	plot(timeVecPreamp, anVsTimePreamp{i}, '.-');
	title(['CTSI preamplified anode ' num2str(anCollectorID(i))]);
end

for i = 1:numAnCollectedCharge
	figure;
	plot(timeVecPreamp, noisyAnVsTimeReg{i}, '-');
	title(['Noisy anode ' num2str(anCollectorID(i))]);
end

for i = 1:numAnCollectedCharge
	figure;
	plot(timeVecPreamp, noisyAnVsTimePreamp{i}, '-');
	title(['Noisy CTSI preamplified anode ' num2str(anCollectorID(i))]);
end

for i = 1:numCaCollectedCharge
	figure;
	plot(timeVecPreamp, caVsTimePreamp{i}, '.-');
	title(['CTSI preamplified cathode ' num2str(caCollectorID(i))]);
end

for i = 1:numCaCollectedCharge
	figure;
	plot(timeVecPreamp, noisyCaVsTimeReg{i}, '-');
	title(['Noisy anode ' num2str(caCollectorID(i))]);
end

for i = 1:numCaCollectedCharge
	figure;
	plot(timeVecPreamp, noisyCaVsTimePreamp{i}, '-');
	title(['Noisy CTSI preamplified cathode ' num2str(caCollectorID(i))]);
end

anodeID
anEnergy
anTriggerTime
noisyAnEnergy
noisyAnTriggerTime

cathodeID
caEnergy
caTriggerTime
noisyCaEnergy
noisyCaTriggerTime

% for h = 1:numAnTrig
% 	for i = 1:numIntrxn
%         for j = 1:numChargeElem
%             len = length(eventAnPhi{anodeID(h) + 1}{i}{j}.EPhi);
%             figure;
%             plot(time(1:len), eventAnPhi{anodeID(h) + 1}{i}{j}.EPhi, '.-');
%             xlabel(['Time (second)']);
%             ylabel(['Weighting potential']);
%             title(['Anode ' num2str(anodeID(h) + 1) ', intrxn ' num2str(i) ', electron element ' num2str(j)]);
%         end
% 	end
% end
% 
% for h = 1:numAnTrig
% 	for i = 1:numIntrxn
%         for j = 1:numChargeElem
%             len = length(eventAnPhi{anodeID(h) + 1}{i}{j}.HPhi);
%             figure;
%             plot(time(1:len), eventAnPhi{anodeID(h) + 1}{i}{j}.HPhi, '.-');
%             xlabel(['Time (second)']);
%             ylabel(['Weighting potential']);
%             title(['Cathode ' num2str(anodeID(h) + 1) ', intrxn ' num2str(i) ', hole element ' num2str(j)]);
%         end
% 	end
% end