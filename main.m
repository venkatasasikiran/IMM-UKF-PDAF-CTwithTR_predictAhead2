close all;
clear;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
turnRate = 5.5;
ParPath = Struct_Param_Path_Init(true,turnRate,90,[6000,4000],400,100);
ParGen = Struct_Param_Gen_Init(1/100000); % clutterdensity inside the parenthesis

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ParPath,devTimeVec]= Path_Target_With_Models(ParPath,ParGen);
% ParPath= Path_Straight_line(ParPath,ParGen);

% figure(1)
% plot(ParPath.TarPathMat(:,1),ParPath.TarPathMat(:,2));
figure(2)
polarplot(ParPath.TarPathMatInPol(:,1),ParPath.TarPathMatInPol(:,2))
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Timings] = Timing_Model(ParPath,ParGen);
if length(Timings.impInstVec)< 30
   while(length(Timings.impInstVec)< 30)
       [ParPath,devTimeVec,~]= Path_Target(ParPath,ParGen);
       [Timings] = Timing_Model(ParPath,ParGen);
   end
end                                   
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParStdDev = Struct_Param_Std_Dev([120,25],deg2rad(1),1);
[Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev);                                                            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ClutterStruct = Clutter_Generation(ParPath.TarMaxRho,ParGen.clutterDensity,Timings.impInstVec);
ClutterStruct = Clutter_Generation(ParPath.TarMaxRho+150,ParGen.clutterDensity,Timings.impInstVec);
%use the second one when we are switching to clutter adjustment case
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transWave.true = Planner_True(ParPath,Timings.impInstVec);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% [FilterEstimate,transWave.estUKF,Measurements] = IMM_Filter(ParPath,ParGen,Timings,Measurements,...
%                                                             ClutterStruct,"C",ParStdDev,transWave.true);
[FilterEstimate,transWave.estIMM,Measurements] = IMM_Filter2(ParPath,ParGen,Timings,Measurements,...
                                                                        ClutterStruct,ParStdDev);
                                                        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
impLength = length(Timings.impInstVec);
Threshold.cw = 250;
Threshold.fm = 250;    

diffBtwTrueExp1 = (ParPath.TarPathMat(Timings.impInstVec(1:end-1),:) - (FilterEstimate.Update(1:2,:))').^2;
diffBtwTrueExp1 = sqrt(diffBtwTrueExp1(:,1) + diffBtwTrueExp1(:,2)); 
CwIndx = find(transWave.estIMM(1:end-1) == "C");
FmIndx = find(transWave.estIMM(1:end-1) == "F");
trackPoints1 = sum(diffBtwTrueExp1(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExp1(FmIndx) <= Threshold.fm);
trackCommInd1 = [FmIndx(diffBtwTrueExp1(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp1(CwIndx) <= Threshold.cw)];
UKFTP1        = (trackPoints1/(impLength-1))*100;
UKFRMS1       = sqrt(sum(diffBtwTrueExp1(trackCommInd1).^2)/length(trackCommInd1));
display(UKFTP1);
display(UKFRMS1);


diffBtwTrueExp2 = (ParPath.TarPathMat(Timings.impInstVec,:) - (FilterEstimate.Predict(1:2,:))').^2;
diffBtwTrueExp2 = sqrt(diffBtwTrueExp2(:,1) + diffBtwTrueExp2(:,2)); 
CwIndx = find(transWave.estIMM == "C");
FmIndx = find(transWave.estIMM == "F");
trackPoints2 = sum(diffBtwTrueExp2(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExp2(FmIndx) <= Threshold.fm);
trackCommInd2 = [FmIndx(diffBtwTrueExp2(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp2(CwIndx)<= Threshold.cw)];
UKFTP2        = (trackPoints2/(impLength))*100;
UKFRMS2       = sqrt(sum(diffBtwTrueExp2(trackCommInd2).^2)/length(trackCommInd2));
display(UKFTP2);
display(UKFRMS2);


     
%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% figure(2)
% plot(FilterEstimate(1,:),FilterEstimate(2,:),'r-*');
% hold on;        
% plot(ParPath.TarPathMat(Timings.impInstVec(1:end-2),1),ParPath.TarPathMat(Timings.impInstVec(1:end-2),2),'-o');
% title('Filtered and True UKF');
% hold off;  
            
figure(3);
plot(diffBtwTrueExp1,'-*');
title('update');

figure(4);
plot(diffBtwTrueExp2,'-*');
title('predict')



% figure(4);
% EstPosCwIndx = find(transWave.estUKF(1:end-2) == "C");
% EstPosFmIndx = find(transWave.estUKF(1:end-2) == "F");
% plot(FilterEstimate(1,:),FilterEstimate(2,:));
% hold on;
% plot(FilterEstimate(1,EstPosCwIndx),FilterEstimate(2,EstPosCwIndx),'ro');
% plot(FilterEstimate(1,EstPosFmIndx),FilterEstimate(2,EstPosFmIndx),'bd');
% hold off;<= Threshold.cw)]
% subplot(2,2,2)
% plot(Timings.impInstVec(1:end-2),FilterEstimate(3,:));
% title('Estimated x vel UKF')
% 
% subplot(2,2,4)
% plot(Timings.impInstVec(1:end-2),FilterEstimate(4,:));
% title('Estimated y vel UKF')  
% 
% 
% figure(6)
% polarplot(ParPath.TarPathMatInPol(:,1),ParPath.TarPathMatInPol(:,2));
% hold on;
% [thetaEst,rhoEst] = cart2pol(FilterEstimate(1,:),FilterEstimate(2,:));
% polarplot(thetaEst,rhoEst);
% polarscatter(thetaEst(EstPosCwIndx),rhoEst(EstPosCwIndx),'red','o','filled');
% polarscatter(thetaEst(EstPosFmIndx),rhoEst(EstPosFmIndx),'blue','d','filled');
% hold off;
% title("red-circle: CW blue-diamond: FM")      

