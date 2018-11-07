close all;
clear;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
turnRate = 5;
ParPath = Struct_Param_Path_Init(true,turnRate,90,[6000,4000],400,100);
ParGen = Struct_Param_Gen_Init(1/100000); % clutterdensity inside the parenthesis

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ParPath,devTimeVec]= Path_Target_With_Models(ParPath,ParGen);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Timings] = Timing_Model(ParPath,ParGen);
if length(Timings.impInstVec)< 30
   while(length(Timings.impInstVec)< 30)
       [ParPath,devTimeVec,~]= Path_Target(ParPath,ParGen);
       [Timings] = Timing_Model(ParPath,ParGen);
   end
end                                   
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParStdDev = Struct_Param_Std_Dev([40,10],deg2rad(1),.5);
[Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev);                                                            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ClutterStruct = Clutter_Generation(ParPath.TarMaxRho,ParGen.clutterDensity,Timings.impInstVec);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transWave.true = Planner_True(ParPath,Timings.impInstVec);
impLength = length(Timings.impInstVec);
Threshold.cw = 250;
Threshold.fm = 250;
figureN = 1;
transWaveAll = string(zeros(impLength,5));
caseLoopIndx = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for caseLoop = ["C","F","I","T","M"]
    [FilterEstimate,transWave.estIMM,Measurements] = IMM_Filter(ParPath,ParGen,Timings,Measurements,...
                                                            ClutterStruct,caseLoop,ParStdDev,transWave.true);
    caseLoopIndx = caseLoopIndx+1;
    transWaveAll(:,caseLoopIndx) = transWave.estIMM;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diffBtwTrueExp = (ParPath.TarPathMat(Timings.impInstVec,:) - (FilterEstimate(1:2,:))').^2;
    diffBtwTrueExp = sqrt(diffBtwTrueExp(:,1) + diffBtwTrueExp(:,2)); 
    CwIndx = find(transWave.estIMM == "C");
    FmIndx = find(transWave.estIMM == "F");
    trackPointsUKF = sum(diffBtwTrueExp(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExp(FmIndx) <= Threshold.fm);
    UKFTP = (trackPointsUKF/impLength)*100;

    trackCommInd = [FmIndx(diffBtwTrueExp(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp(CwIndx) <= Threshold.cw)];
    UKFRMS       = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)/length(trackCommInd));
    disp(['Track Probability in ',char(caseLoop),' is ', num2str(UKFTP)]);
    disp(['RMS in ',char(caseLoop),' is ', num2str(UKFRMS)]);
     
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    figureN = figureN + 1;
    figure(figureN)
    plot(FilterEstimate(1,:),FilterEstimate(2,:),'r-*');
    hold on;        
    plot(ParPath.TarPathMat(Timings.impInstVec,1),ParPath.TarPathMat(Timings.impInstVec,2),'-o');
    title('Filtered and True ');
    hold off;  

    figureN = figureN + 1;
    figure(figureN)
    plot(diffBtwTrueExp,'-*');
    title('diff true exp');

    figureN = figureN + 1;
    figure(figureN)
    EstPosCwIndx = find(transWave.estIMM == "C");
    EstPosFmIndx = find(transWave.estIMM == "F");
    plot(FilterEstimate(1,:),FilterEstimate(2,:));
    hold on;
    plot(FilterEstimate(1,EstPosCwIndx),FilterEstimate(2,EstPosCwIndx),'ro');
    plot(FilterEstimate(1,EstPosFmIndx),FilterEstimate(2,EstPosFmIndx),'bd');
    hold off;
    title("red-circle: CW blue-diamond: FM")

%     figure(5)  
%     subplot(2,2,1)
%     plot(Timings.impInstVec,ParPath.TrueVelocities(:,1));
%     title('true x vel')
% 
%     subplot(2,2,3)
%     plot(Timings.impInstVec,ParPath.TrueVelocities(:,2));
%     title('true y vel')
% 
%     subplot(2,2,2)
%     plot(Timings.impInstVec,FilterEstimate(3,:));
%     title('Estimated x vel UKF')
% 
%     subplot(2,2,4)
%     plot(Timings.impInstVec,FilterEstimate(4,:));
%     title('Estimated y vel UKF')  


%     figure(6)
%     polarplot(ParPath.TarPathMatInPol(:,1),ParPath.TarPathMatInPol(:,2));
%     hold on;
%     [thetaEst,rhoEst] = cart2pol(FilterEstimate(1,:),FilterEstimate(2,:));
%     polarplot(thetaEst,rhoEst);
%     polarscatter(thetaEst(EstPosCwIndx),rhoEst(EstPosCwIndx),'red','o','filled');
%     polarscatter(thetaEst(EstPosFmIndx),rhoEst(EstPosFmIndx),'blue','d','filled');
%     hold off;
%     title("red-circle: CW blue-diamond: FM")      
end
