close all;
clear;
load('MaxAcc.mat');
clutterVersion = [100000];
for clutterLoop = 1:length(clutterVersion)
    yawDevRateLoopMin = 3;
    yawDevRateLoopJump = .5;
    yawDevRateLoopMax = 5.5;
    nYaw = (yawDevRateLoopMax/yawDevRateLoopJump)+1;
    nYawStart = yawDevRateLoopMin/yawDevRateLoopJump;
    if clutterLoop == 1
        delete performance.mat
        delete loopParam.mat                   
        IMM_TP.C = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.C = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.F = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.F = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.I = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.I = zeros(length(clutterVersion),nYaw-nYawStart);
%         IMM_TP.T = zeros(length(clutterVersion),nYaw-nYawStart);
%         IMM_RMS.T = zeros(length(clutterVersion),nYaw-nYawStart);
        
        UKF_TP.C = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS.C = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_TP.F = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS.F = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_TP.I = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS.I = zeros(length(clutterVersion),nYaw-nYawStart);
%         UKF_TP.T = zeros(length(clutterVersion),nYaw-nYawStart);
%         UKF_RMS.T = zeros(length(clutterVersion),nYaw-nYawStart);
%         tempValues(nYaw-nYawStart) = struct;
        save('performance','IMM_TP','IMM_RMS','UKF_TP','UKF_RMS');        
    end
    for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
    totalSim = 500;        
    tempTP.C = zeros(totalSim,1);
    tempRMS.C = zeros(totalSim,1);    
    tempTP.F = zeros(totalSim,1);
    tempRMS.F = zeros(totalSim,1);
    tempTP.I = zeros(totalSim,1);
    tempRMS.I = zeros(totalSim,1);
%     tempTP.T = zeros(totalSim,1);
%     tempRMS.T = zeros(totalSim,1);
    
    tempUTP.C = zeros(totalSim,1);
    tempURMS.C = zeros(totalSim,1);    
    tempUTP.F = zeros(totalSim,1);
    tempURMS.F = zeros(totalSim,1);
    tempUTP.I = zeros(totalSim,1);
    tempURMS.I = zeros(totalSim,1);
%     tempUTP.T = zeros(totalSim,1);
%     tempURMS.T = zeros(totalSim,1);
        for nSim = 1:totalSim
            if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                save('loopParam','clutterVersion','totalSim')                   
            else                
                load('loopParam');
            end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
            ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],400,100);
            ParGen = Struct_Param_Gen_Init(1/clutterVersion(clutterLoop));% clutterdensity inside the parenthesis

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
%             ClutterStruct = Clutter_Generation(ParPath.TarMaxRho,ParGen.clutterDensity,Timings.impInstVec);
            ClutterStruct = Clutter_Generation(ParPath.TarMaxRho+150,ParGen.clutterDensity,Timings.impInstVec);
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            transWave.true = Planner_True(ParPath,Timings.impInstVec);
            impLength = length(Timings.impInstVec);
            Threshold.cw = 250;
            Threshold.fm = 250;
            acc = MaxAcc(2,(MaxAcc(1,:) == yawDevRateLoop));
            noiseAcc = [acc;acc];
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            for caseLoop = ["C","F","I"]
                [FilterEstimate,transWave.estIMM,~] = IMM_Filter(ParPath,ParGen,Timings,Measurements,...
                                                                        ClutterStruct,caseLoop,ParStdDev,transWave.true);
                [FilterEstimateUKF,transWave.estUKF,~] = UKF_PDAF(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,caseLoop,ParStdDev,transWave.true);           
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                diffBtwTrueExp = (ParPath.TarPathMat(Timings.impInstVec,:) - (FilterEstimate(1:2,:))').^2;
                diffBtwTrueExp = sqrt(diffBtwTrueExp(:,1) + diffBtwTrueExp(:,2)); 
                CwIndxIMM = find(transWave.estIMM == "C");
                FmIndxIMM = find(transWave.estIMM == "F");
                trackPoints = sum(diffBtwTrueExp(CwIndxIMM) <= Threshold.cw) + sum(diffBtwTrueExp(FmIndxIMM) <= Threshold.fm);
                trackCommInd = [FmIndxIMM(diffBtwTrueExp(FmIndxIMM) <= Threshold.fm);CwIndxIMM(diffBtwTrueExp(CwIndxIMM) <= Threshold.cw)];
                
                XPlusVecUKF = FilterEstimateUKF.XkPlusVec(:,2:end)';
                diffBtwTrueExpUKF = (ParPath.TarPathMat(Timings.impInstVec,:) - XPlusVecUKF(:,1:2)).^2;
                diffBtwTrueExpUKF = sqrt(diffBtwTrueExpUKF(:,1) + diffBtwTrueExpUKF(:,2)); 
                CwIndx = find(transWave.estUKF == "C");
                FmIndx = find(transWave.estUKF == "F");
                trackPointsUKF = sum(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);
                trackCommIndUKF = [FmIndx(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw)];
                
                trackComm = intersect(trackCommInd,trackCommIndUKF);
                switch caseLoop
                    case "C" 
                        tempTP.C(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.C(nSim) = sqrt(sum(diffBtwTrueExp(trackComm).^2)...
                                                /length(trackComm));
                        tempUTP.C(nSim) = (trackPointsUKF/impLength)*100;        
                        tempURMS.C(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackComm).^2)...
                                                /length(trackComm));
                    case "F"
                        tempTP.F(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.F(nSim) = sqrt(sum(diffBtwTrueExp(trackComm).^2)...
                                                /length(trackComm));
                        tempUTP.F(nSim) = (trackPointsUKF/impLength)*100;        
                        tempURMS.F(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackComm).^2)...
                                                /length(trackComm));
                    case "I"
                        tempTP.I(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.I(nSim) = sqrt(sum(diffBtwTrueExp(trackComm).^2)...
                                                /length(trackComm));
                        tempUTP.I(nSim) = (trackPointsUKF/impLength)*100;        
                        tempURMS.I(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackComm).^2)...
                                                /length(trackComm));
%                     case "T"
%                         tempTP.T(nSim) = (trackPoints/impLength)*100;        
%                         tempRMS.T(nSim) = sqrt(sum(diffBtwTrueExp(trackComm).^2)...
%                                                 /length(trackComm));               
%                         tempUTP.T(nSim) = (trackPointsUKF/impLength)*100;        
%                         tempURMS.T(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackComm).^2)...
%                                                 /length(trackComm));
                end
                if nSim == totalSim
                    load('performance');
                    IMM_TP.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.C);
                    IMM_RMS.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.C);
                    IMM_TP.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.F);
                    IMM_RMS.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.F);
                    IMM_TP.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.I);
                    IMM_RMS.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.I);
%                     IMM_TP.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.T);
%                     IMM_RMS.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.T);
%                     IMM_TP.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.T);
%                     IMM_RMS.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.T);                    
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).TP = tempTP;
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).RMS = tempRMS;
                    
                    UKF_TP.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempUTP.C);
                    UKF_RMS.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempURMS.C);
                    UKF_TP.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempUTP.F);
                    UKF_RMS.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempURMS.F);
                    UKF_TP.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempUTP.I);
                    UKF_RMS.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempURMS.I);
%                     UKF_TP.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempUTP.T);
%                     UKF_RMS.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempURMS.T);
%                     UKF_TP.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempUTP.T);
%                     UKF_RMS.T(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempURMS.T);                    
                    tempValU(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).TP = tempUTP;
                    tempValU(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).RMS = tempURMS;
                    
                end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            end
            save('loopParam','clutterVersion','totalSim');
            if nSim == totalSim
                save('performance','IMM_TP','IMM_RMS','UKF_TP','UKF_RMS','tempVal','tempValU');
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart MaxAcc;
            else
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart tempTP tempRMS tempUTP tempURMS MaxAcc;
            end 
        end
        disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
    end
    disp(['cluttterLoop:',num2str(clutterLoop)]);    
    clearvars -except MaxAcc;
end
