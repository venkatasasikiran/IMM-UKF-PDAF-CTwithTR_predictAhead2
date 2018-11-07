close all;
clear;
clutterVersion = [100000];
for clutterLoop = 1:length(clutterVersion)
    yawDevRateLoopMin = 4.5;
    yawDevRateLoopJump = .5;
    yawDevRateLoopMax = 6;
    nYaw = (yawDevRateLoopMax/yawDevRateLoopJump)+1;
    nYawStart = yawDevRateLoopMin/yawDevRateLoopJump;
    if clutterLoop == 1
        delete performance.mat
        delete loopParam.mat                   
        IMM_TP.C1 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.C1 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.F1 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.F1 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.I1 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.I1 = zeros(length(clutterVersion),nYaw-nYawStart);
        
        IMM_TP.C2 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.C2 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.F2 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.F2 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_TP.I2 = zeros(length(clutterVersion),nYaw-nYawStart);
        IMM_RMS.I2 = zeros(length(clutterVersion),nYaw-nYawStart); 
        save('performance','IMM_TP','IMM_RMS');        
    end
    for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
    totalSim = 1000;        
    tempTP.C = zeros(totalSim,3);
    tempRMS.C = zeros(totalSim,3);    
    tempTP.F = zeros(totalSim,3);
    tempRMS.F = zeros(totalSim,3);
    tempTP.I = zeros(totalSim,3);
    tempRMS.I = zeros(totalSim,3);    
        for nSim = 1:totalSim
            if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                save('loopParam','clutterVersion','totalSim')                   
            else                
                load('loopParam');
                
            end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
            ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[5000,3000],300,100);
            ParGen = Struct_Param_Gen_Init(1/clutterVersion(clutterLoop));% clutterdensity inside the parenthesis

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ParPath,devTimeVec]= Path_Target_With_Models(ParPath,ParGen);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Timings] = Timing_Model(ParPath,ParGen);                                           
    %% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ParStdDev = Struct_Param_Std_Dev([120,25],deg2rad(1),1);
            [Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev);                                                            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ClutterStruct = Clutter_Generation(ParPath.TarMaxRho+150,ParGen.clutterDensity,Timings.impInstVec);
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            transWave.true = Planner_True(ParPath,Timings.impInstVec);
            impLength = length(Timings.impInstVec);
            Threshold.cw = 250;
            Threshold.fm = 250;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            for caseLoop = ["C","F","I"]
                [FilterEstimate,transWave.estIMM,Measurements] = IMM_Filter(ParPath,ParGen,Timings,Measurements,...
                                                                        ClutterStruct,caseLoop,ParStdDev,transWave.true);

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                diffBtwTrueExp = (ParPath.TarPathMat(Timings.impInstVec,:) - (FilterEstimate(1:2,:))').^2;
                diffBtwTrueExp = sqrt(diffBtwTrueExp(:,1) + diffBtwTrueExp(:,2)); 
                CwIndx = find(transWave.estIMM == "C");
                FmIndx = find(transWave.estIMM == "F");
                trackPoints = sum(diffBtwTrueExp(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExp(FmIndx) <= Threshold.fm);
                
                switch caseLoop
                    case "C" 
                        tempTP.C(nSim,2) = trackPoints;
                        tempTP.C(nSim,3) = impLength;
                        tempTP.C(nSim,1) = (trackPoints/impLength)*100;        
                        
                        
                        diffVec.C  = diffBtwTrueExp;
                        trackInd.C = [FmIndx(diffBtwTrueExp(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp(CwIndx) <= Threshold.cw)];                    
                        tempRMS.C(nSim,1) = sqrt(sum(diffBtwTrueExp(trackInd.C).^2)...
                                                /length(trackInd.C));
                    case "F"
                        tempTP.F(nSim,2) = trackPoints;
                        tempTP.F(nSim,3) = impLength;
                        tempTP.F(nSim,1) = (trackPoints/impLength)*100;
                        tempTP.F(nSim,1) = (trackPoints/impLength)*100;        
                        
                        
                        diffVec.F  = diffBtwTrueExp;
                        trackInd.F = [FmIndx(diffBtwTrueExp(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp(CwIndx) <= Threshold.cw)];                    
                        tempRMS.F(nSim,1) = sqrt(sum(diffBtwTrueExp(trackInd.F).^2)...
                                                /length(trackInd.F));
                    case "I"
                        tempTP.I(nSim,2) = trackPoints;
                        tempTP.I(nSim,3) = impLength;
                        tempTP.I(nSim,1) = (trackPoints/impLength)*100;
                        tempTP.I(nSim,1) = (trackPoints/impLength)*100;        
                        
                        
                        
                        diffVec.I  = diffBtwTrueExp;
                        trackInd.I = [FmIndx(diffBtwTrueExp(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp(CwIndx) <= Threshold.cw)];                    
                        tempRMS.I(nSim,1) = sqrt(sum(diffBtwTrueExp(trackInd.I).^2)...
                                                /length(trackInd.I));
                        
                        trackCommInd = intersect(intersect(trackInd.C,trackInd.F),trackInd.I);
                        tempRMS.C(nSim,2) = sum(diffVec.C(trackCommInd).^2);
                        tempRMS.C(nSim,3) = length(trackCommInd);
                        tempRMS.F(nSim,2) = sum(diffVec.F(trackCommInd).^2);
                        tempRMS.F(nSim,3) = length(trackCommInd);
                        tempRMS.I(nSim,2) = sum(diffVec.I(trackCommInd).^2);
                        tempRMS.I(nSim,3) = length(trackCommInd);
                end
                if nSim == totalSim
                    load('performance');
                    IMM_TP.C1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.C(:,1));
                    IMM_RMS.C1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.C(:,1));
                    IMM_TP.F1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.F(:,1));
                    IMM_RMS.F1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.F(:,1));
                    IMM_TP.I1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.I(:,1));
                    IMM_RMS.I1(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.I(:,1));
                    
                    IMM_TP.C2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = ...
                                                                         (sum(tempTP.C(:,2))/sum(tempTP.C(:,3)))*100;
                    IMM_RMS.C2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) =...
                                                                        sqrt(sum(tempRMS.C(:,2))/sum(tempRMS.C(:,3)));
                    IMM_TP.F2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = ...
                                                                        (sum(tempTP.F(:,2))/sum(tempTP.F(:,3)))*100;
                    IMM_RMS.F2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = ...
                                                                        sqrt(sum(tempRMS.F(:,2))/sum(tempRMS.F(:,3)));
                    IMM_TP.I2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = ...
                                                                        (sum(tempTP.I(:,2))/sum(tempTP.I(:,3)))*100;
                    IMM_RMS.I2(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = ...
                                                                        sqrt(sum(tempRMS.I(:,2))/sum(tempRMS.I(:,3)));
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).TP = tempTP;
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).RMS = tempRMS;
                end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            end
            save('loopParam','clutterVersion','totalSim');
            if nSim == totalSim
                save('performance','IMM_TP','IMM_RMS','tempVal');
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart;
            else
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart tempTP tempRMS;
            end 
        end
        disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
    end
    disp(['cluttterLoop:',num2str(clutterLoop)]);    
    clear;
end
