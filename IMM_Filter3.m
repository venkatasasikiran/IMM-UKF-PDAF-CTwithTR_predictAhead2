function [FilterEstimate,transWave,Measurements] = IMM_Filter3(ParPath,ParGen,Timings,Measurements,...
                                                        ClutterStruct,ParStdDev)
    
    n = size(Timings.impInstVec,1);
    extImpInstVec = [1;Timings.impInstVec];
    stateSizeCV = 4;
    stateSizeCT = 4;
    transWave = string(zeros(n,1));
    FilterEstimate = zeros(min(stateSizeCV,stateSizeCT),n);
    
    FilterEstimateCV = Struct_Param_Filter_Init(ParPath,stateSizeCV,n,"CV");
    FilterEstimateCTpos = Struct_Param_Filter_Init(ParPath,stateSizeCT,n,"CT");
    FilterEstimateCTneg = Struct_Param_Filter_Init(ParPath,stateSizeCT,n,"CT");
    ParPDAF = Struct_Param_PDAF_init;
    
    Measurements.MeasFinal.Pos = zeros(n,2);
    Measurements.MeasFinal.rangeRate = nan(n,1);
    
    modeTransProb = [0.90 0.05  0.05;
                     0.19 0.80  0.01;
                     0.19 0.01  0.80];
    modeProb      = [0.90;0.05;0.05];
    
    for i = 1:n
        if i == 1
            T = (extImpInstVec(i+1) - extImpInstVec(i))/ParGen.f;
            [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,cNorm] = Mixing_Interaction(FilterEstimateCV,FilterEstimateCTpos,...
                                                                            FilterEstimateCTneg,modeTransProb,modeProb);
                                                                  
            FilterMat = Struct_Mat_Filter(T,ParStdDev,ParPath.yawDevRate);                                    
            [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg] = Estimating_States(FilterEstimateCV,FilterEstimateCTpos,...
                                                                                        FilterEstimateCTneg,FilterMat,i);
            FilterPredict              = cNorm(1)*FilterEstimateCV.XkMinus + cNorm(2)*FilterEstimateCTpos.XkMinus...
                                        + cNorm(3)*FilterEstimateCTneg.XkMinus;          
            transWave(i,:)             = Planner_Est(FilterPredict,ParPath,"I");
        end
        [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,z,LambdaAll] = Update_States_Clutter_Adjust(FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,...
                                                                                         FilterMat,Measurements,ClutterStruct,transWave(i),ParPDAF,i);
                     
        modeProb = Update_Mode_Prob(LambdaAll,cNorm);
        FilterEstimate(:,i) = FilterEstimateCV.XkPlus*modeProb(1) + FilterEstimateCTpos.XkPlus*modeProb(2)...
                                        + FilterEstimateCTneg.XkPlus*modeProb(3);
        switch transWave(i)
            case "C"
                Measurements.MeasFinal.Pos(i,:) = z(1:2)';
                Measurements.MeasFinal.rangeRate(i) = z(3);
            case "F"
                Measurements.MeasFinal.Pos(i,:) = z';
        end
        if i < n
            T2 = (extImpInstVec(i+2) - extImpInstVec(i+1))/ParGen.f;
            [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,cNorm] = Mixing_Interaction(FilterEstimateCV,FilterEstimateCTpos,...
                                                                            FilterEstimateCTneg,modeTransProb,modeProb);
                                                                  
            FilterMat = Struct_Mat_Filter(T2,ParStdDev,ParPath.yawDevRate);                                    
            [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg] = Estimating_States(FilterEstimateCV,FilterEstimateCTpos,...
                                                                                        FilterEstimateCTneg,FilterMat,i+1);
            if i ==1            
                FilterPredict              = cNorm(1)*FilterEstimateCV.XkMinus + cNorm(2)*FilterEstimateCTpos.XkMinus...
                                            + cNorm(3)*FilterEstimateCTneg.XkMinus;          
                transWave(i+1,:)             = Planner_Est(FilterPredict,ParPath,"I");
            end
        end
        if i < n-1
            T3 = (extImpInstVec(i+3) - extImpInstVec(i+2))/ParGen.f;
            [FilterEstimateCV2,FilterEstimateCTpos2,FilterEstimateCTneg2,cNorm2] = Mixing_Interaction2(FilterEstimateCV,FilterEstimateCTpos,...
                                                                                FilterEstimateCTneg,modeTransProb,cNorm);
            FilterMat2 = Struct_Mat_Filter(T3,ParStdDev,ParPath.yawDevRate);
            [FilterEstimateCV2,FilterEstimateCTpos2,FilterEstimateCTneg2] = Estimating_States(FilterEstimateCV2,FilterEstimateCTpos2,...
                                                                                FilterEstimateCTneg2,FilterMat2,i+2);
            FilterPredict                   = cNorm2(1)*FilterEstimateCV2.XkMinus + cNorm2(2)*FilterEstimateCTpos2.XkMinus...
                                            + cNorm2(3)*FilterEstimateCTneg2.XkMinus;
            transWave(i+2,:)                = Planner_Est(FilterPredict,ParPath,"I");
        end
    end
        
       


end

function [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,cNorm] = Mixing_Interaction(FilterEstimateCV,...
                                                                            FilterEstimateCTpos,FilterEstimateCTneg,modeTransProb,modeProb)
    nModes = length(modeProb);
    cNorm = zeros(nModes,1);
    MixingProb = zeros(nModes);  
    

    for i1=1:nModes
        cNorm(i1) = (modeTransProb(:,i1))'*modeProb;        
    end

    for i=1:nModes
        for j = 1:nModes
           MixingProb(i,j) = modeTransProb(i,j)*modeProb(i)/cNorm(j);
        end
    end
    
    
    XPlusVec = [FilterEstimateCV.XkPlus FilterEstimateCTpos.XkPlus FilterEstimateCTneg.XkPlus];
    PPlus(1).Mat = FilterEstimateCV.PkPlus;
    PPlus(2).Mat = FilterEstimateCTpos.PkPlus;
    PPlus(3).Mat = FilterEstimateCTneg.PkPlus;
    
    FilterEstimateCV.XkPlusMix = XPlusVec(:,1)*MixingProb(1,1) + XPlusVec(:,2)*MixingProb(2,1) + XPlusVec(:,3)*MixingProb(3,1);
    FilterEstimateCTpos.XkPlusMix = XPlusVec(:,1)*MixingProb(1,2) + XPlusVec(:,2)*MixingProb(2,2) + XPlusVec(:,3)*MixingProb(3,2);
    FilterEstimateCTneg.XkPlusMix = XPlusVec(:,1)*MixingProb(1,3) + XPlusVec(:,2)*MixingProb(2,3) + XPlusVec(:,3)*MixingProb(3,3);
    
    FilterEstimateCV.PkPlusMix = zeros(4);
    FilterEstimateCTpos.PkPlusMix = zeros(4);
    FilterEstimateCTneg.PkPlusMix = zeros(4);
    for i2 = 1:nModes
        FilterEstimateCV.PkPlusMix = FilterEstimateCV.PkPlusMix + MixingProb(i2,1)*(PPlus(i2).Mat ...
                                   + (XPlusVec(:,i2)- FilterEstimateCV.XkPlusMix)*(XPlusVec(:,i2)- FilterEstimateCV.XkPlusMix)');
    end
    for i3 = 1:nModes
        FilterEstimateCTpos.PkPlusMix = FilterEstimateCTpos.PkPlusMix + MixingProb(i3,2)*(PPlus(i3).Mat ...
                                      + (XPlusVec(:,i3)- FilterEstimateCTpos.XkPlusMix)*(XPlusVec(:,i3)- FilterEstimateCTpos.XkPlusMix)');
    end
    for i4 = 1:nModes
        FilterEstimateCTneg.PkPlusMix = FilterEstimateCTneg.PkPlusMix + MixingProb(i4,3)*(PPlus(i4).Mat ...
                                      + (XPlusVec(:,i4)- FilterEstimateCTneg.XkPlusMix)*(XPlusVec(:,i4)- FilterEstimateCTneg.XkPlusMix)');
    end
end
function [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,cNorm2] = Mixing_Interaction2(FilterEstimateCV,...
                                                                            FilterEstimateCTpos,FilterEstimateCTneg,modeTransProb,cNorm)
    nModes = length(cNorm);
    cNorm2 = zeros(nModes,1);
    MixingProb = zeros(nModes);  
    

    for i1=1:nModes
        cNorm2(i1) = (modeTransProb(:,i1))'*cNorm;        
    end

    for i=1:nModes
        for j = 1:nModes
           MixingProb(i,j) = modeTransProb(i,j)*cNorm(i)/cNorm2(j);
        end
    end
    
    
    XPlusVec = [FilterEstimateCV.XkMinus FilterEstimateCTpos.XkMinus FilterEstimateCTneg.XkMinus];
    PPlus(1).Mat = FilterEstimateCV.PkMinus;
    PPlus(2).Mat = FilterEstimateCTpos.PkMinus;
    PPlus(3).Mat = FilterEstimateCTneg.PkMinus;
    
    FilterEstimateCV.XkPlusMix = XPlusVec(:,1)*MixingProb(1,1) + XPlusVec(:,2)*MixingProb(2,1) + XPlusVec(:,3)*MixingProb(3,1);
    FilterEstimateCTpos.XkPlusMix = XPlusVec(:,1)*MixingProb(1,2) + XPlusVec(:,2)*MixingProb(2,2) + XPlusVec(:,3)*MixingProb(3,2);
    FilterEstimateCTneg.XkPlusMix = XPlusVec(:,1)*MixingProb(1,3) + XPlusVec(:,2)*MixingProb(2,3) + XPlusVec(:,3)*MixingProb(3,3);
    
    FilterEstimateCV.PkPlusMix = zeros(4);
    FilterEstimateCTpos.PkPlusMix = zeros(4);
    FilterEstimateCTneg.PkPlusMix = zeros(4);
    for i2 = 1:nModes
        FilterEstimateCV.PkPlusMix = FilterEstimateCV.PkPlusMix + MixingProb(i2,1)*(PPlus(i2).Mat ...
                                   + (XPlusVec(:,i2)- FilterEstimateCV.XkPlusMix)*(XPlusVec(:,i2)- FilterEstimateCV.XkPlusMix)');
    end
    for i3 = 1:nModes
        FilterEstimateCTpos.PkPlusMix = FilterEstimateCTpos.PkPlusMix + MixingProb(i3,2)*(PPlus(i3).Mat ...
                                      + (XPlusVec(:,i3)- FilterEstimateCTpos.XkPlusMix)*(XPlusVec(:,i3)- FilterEstimateCTpos.XkPlusMix)');
    end
    for i4 = 1:nModes
        FilterEstimateCTneg.PkPlusMix = FilterEstimateCTneg.PkPlusMix + MixingProb(i4,3)*(PPlus(i4).Mat ...
                                      + (XPlusVec(:,i4)- FilterEstimateCTneg.XkPlusMix)*(XPlusVec(:,i4)- FilterEstimateCTneg.XkPlusMix)');
    end
end


function FilterMat = Struct_Mat_Filter(T,ParStdDev,yawDevRate)

    FilterMat.T = T;
    FilterMat.CV.F = [1 0 T 0;...
                      0 1 0 T;...
                      0 0 1 0;...
                      0 0 0 1];
    FilterMat.CT.Fpos = [1         0         sin(yawDevRate*T)/yawDevRate         (cos(yawDevRate*T)- 1)/yawDevRate;...
                         0         1      (1 - cos(yawDevRate*T))/yawDevRate       sin(yawDevRate*T)/yawDevRate;...
                         0         0          cos(yawDevRate*T)           -sin(yawDevRate*T);...
                         0         0          sin(yawDevRate*T)            cos(yawDevRate*T)];
    FilterMat.CT.Fneg = [1         0         sin(-yawDevRate*T)/-yawDevRate         (cos(-yawDevRate*T)- 1)/-yawDevRate;...
                         0         1      (1 - cos(-yawDevRate*T))/-yawDevRate       sin(-yawDevRate*T)/-yawDevRate;...
                         0         0          cos(-yawDevRate*T)           -sin(-yawDevRate*T);...
                         0         0          sin(-yawDevRate*T)            cos(-yawDevRate*T)];
    FilterMat.CV.L = [0.5*(T)^2 0                 ;...     
                         0 0.5*(T)^2 ;...
                         T 0                 ;...
                         0 T]; 
    FilterMat.CT.L = [0.5*(T)^2 0                 ;...     
                         0 0.5*(T)^2 ;...
                         T 0                 ;...
                         0 T];                   
    qX = .1;%noiseAcc(1)*1;
    qY = .1;%noiseAcc(2)*1;
    Qxy = [qX^2 0;
           0    qY^2];
    FilterMat.CV.Q = FilterMat.CV.L*Qxy*FilterMat.CV.L';
    FilterMat.CT.Q = FilterMat.CT.L*Qxy*FilterMat.CT.L';
    
    FilterMat.R.Cw = [ParStdDev.Range.cw^2             0               0;  
                               0               ParStdDev.theta^2 0;
                               0                           0               ParStdDev.RangeRate.cw^2];
    FilterMat.R.Fm = [ParStdDev.Range.fm^2             0;  
                               0              ParStdDev.theta^2];
end


function ParPDAF = Struct_Param_PDAF_init()
    ParPDAF.gateGamma = 16;  
    ParPDAF.Pd = 1;
    ParPDAF.Pg.fm = 0.9997;
    ParPDAF.Pg.cw = 0.9989;
    ParPDAF.clutterDensity = 1/100000;
end

function newModeProb = Update_Mode_Prob(Lambda,cNorm)
    if Lambda(1) <= 10^-25 && Lambda(2) <= 10^-25 && Lambda(3) <= 10^-25
        newModeProb = cNorm;
%         newModeProb = [1/3;1/3;1/3];
    else
        newModeProb = zeros(length(cNorm),1);
        c = Lambda*cNorm;
        for i = 1:length(cNorm)
            newModeProb(i) = Lambda(i)*cNorm(i)/c;
        end
    end
%     NrNaN = sum(isnan(newModeProb(:)));
%     if NrNaN > 0
%         newModeProb = prevModeProb;
%     end
end