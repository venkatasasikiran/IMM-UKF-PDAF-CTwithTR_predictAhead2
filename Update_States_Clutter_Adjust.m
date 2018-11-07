function [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,z,LambdaAll] = Update_States_Clutter_Adjust(FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg,FilterMat,Measurements,ClutterStruct,transWave,ParPDAF,loopNum)
    switch transWave
        case "C"
            z = [Measurements.MeasRange.cw(loopNum);Measurements.MeasTheta(loopNum);...
                Measurements.MeasRangeRate.cw(loopNum)];
            
        case "F"
            z = [Measurements.MeasRange.fm(loopNum);Measurements.MeasTheta(loopNum)];            
    end        
[zHatCV,ZPointsCV] = Measurement_Predict(FilterEstimateCV,transWave);
[zHatCTp,ZPointsCTp] = Measurement_Predict(FilterEstimateCTpos,transWave);
[zHatCTn,ZPointsCTn] = Measurement_Predict(FilterEstimateCTneg,transWave);
maxEstRange = max([zHatCV(1),zHatCTn(1),zHatCTp(1)]) + 100;
if maxEstRange > max(ClutterStruct(loopNum).rhoClutter)
    ClutterStruct = Clutter_Adjustment(maxEstRange,ParPDAF.clutterDensity,ClutterStruct,loopNum);
end
[FilterEstimateCV,LambdaCV,~] = UKF_Update_CV(FilterEstimateCV,FilterMat,z,loopNum,transWave,ParPDAF,ClutterStruct,zHatCV,ZPointsCV);
[FilterEstimateCTpos,LambdaCTpos,~] = UKF_Update_CT(FilterEstimateCTpos,FilterMat,z,loopNum,transWave,ParPDAF,ClutterStruct,zHatCTp,ZPointsCTp);
[FilterEstimateCTneg,LambdaCTneg,~] = UKF_Update_CT(FilterEstimateCTneg,FilterMat,z,loopNum,transWave,ParPDAF,ClutterStruct,zHatCTn,ZPointsCTn);
LambdaAll = [LambdaCV,LambdaCTpos,LambdaCTneg];

end

function [zHat,ZPoints] = Measurement_Predict(FilterEstimate,transWave)
    XMinusPts = FilterEstimate.XMinusPts;
    XkMinus   = FilterEstimate.XkMinus;    
    n = length(XkMinus);
    alpha = 10^-3;
    kappa = 0;    
    lamda = alpha^2*(n+kappa) - n;    

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for i=1:2*n
        wPts(i+1)       = 1/(2*(n+lamda));
    end
    
    ZPoints             = Trans_H(XMinusPts,transWave);
    zHat                = ZPoints*wPts';
end

function  ClutterStruct = Clutter_Adjustment(maxRho,clutterDensity,ClutterStruct,loopNum)
    areaSquare = (2*maxRho)^2;
    ptsGen = ceil(areaSquare*clutterDensity);    
    xClutter = -maxRho + 2*maxRho*rand(ptsGen,1);          %This is uniform distributiom
    yClutter = -maxRho + 2*maxRho*rand(ptsGen,1);          %This is uniform distributiom
    [thetaClutter,rhoClutter] = cart2pol(xClutter,yClutter);
    indx = find(rhoClutter<=maxRho);
    ClutterStruct(loopNum).thetaClutter = thetaClutter(indx);
    ClutterStruct(loopNum).rhoClutter = rhoClutter(indx);
    noClutterPts = length(indx);        
    ClutterStruct(loopNum).rangeRateClutter = 5*randn(noClutterPts,1);
end

function [FilterEstimateCV,LambdaCV,M]= UKF_Update_CV(FilterEstimateCV,FilterMat,z,loopNum,transWave,ParPDAF,ClutterStruct,zHat,ZPoints)
    XMinusPts = FilterEstimateCV.XMinusPts;
    XkMinus   = FilterEstimateCV.XkMinus;
    PkMinus   = FilterEstimateCV.PkMinus;
    n = length(XkMinus);
    alpha = 10^-3;
    kappa = 0;
    beta = 2;
    lamda = alpha^2*(n+kappa) - n;    

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for i=1:2*n
        wPts(i+1)       = 1/(2*(n+lamda));
    end
    
%     ZPoints             = Trans_H(XMinusPts,transWave);
%     zHat                = ZPoints*wPts';
    if transWave == "C"
        estRangeRate        = zHat(end);
        FilterEstimateCV.estRangeRate(loopNum) = estRangeRate;
    end    
    m1 = size(XMinusPts,1);
    m2 = size(zHat,1);
    Pzz = zeros(m2);
    for i = 0:2*n
        if i == 0
            Pzz = Pzz...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(ZPoints(:,i+1)- zHat)*(ZPoints(:,i+1)- zHat)';
        else
            Pzz = Pzz + wPts(i+1)...
                      *(ZPoints(:,i+1)- zHat)*(ZPoints(:,i+1)- zHat)';
        end
        
    end
    
    Pxz = zeros(m1,m2);
    for i = 0:2*n
        if i == 0
            Pxz = Pxz...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(XMinusPts(:,i+1)-XkMinus)*(ZPoints(:,i+1)- zHat)';
        else
            Pxz = Pxz + wPts(i+1)...
                      *(XMinusPts(:,i+1)- XkMinus)*(ZPoints(:,i+1)- zHat)';
        end        
    end    
    
    switch transWave
        case "C"
            S_k = Pzz + FilterMat.R.Cw;
            K_k = Pxz/(S_k);
            nZ = 3;
            cNZ = pi^(nZ/2)/gamma(nZ/2 + 1);
            ClutterMat = [ClutterStruct(loopNum).rhoClutter';ClutterStruct(loopNum).thetaClutter';ClutterStruct(loopNum).rangeRateClutter'];
            ClutterMat = [ClutterMat z];
        case "F"
            S_k = Pzz + FilterMat.R.Fm;
            K_k = Pxz/(S_k);
            nZ = 2;
            cNZ = pi^(nZ/2)/gamma(nZ/2 + 1);
            ClutterMat = [ClutterStruct(loopNum).rhoClutter';ClutterStruct(loopNum).thetaClutter'];
            ClutterMat = [ClutterMat z];
    end 
    [zGateVec,M] = Gating(ClutterMat,zHat,S_k,ParPDAF);
    vVec = zGateVec - repmat(zHat,1,size(zGateVec,2));
    betaVec = Data_Assoc(vVec,M,S_k,cNZ,ParPDAF,nZ,transWave);   
    if M == 0
        PkPlus = PkMinus;
        XkPlus = XkMinus;
        vVec = zeros(length(z),1);
    else
        vCom = vVec*betaVec(2:end);
        PkPlus = betaVec(1)*PkMinus + (1-betaVec(1))*(PkMinus - K_k*S_k*K_k')...
                 + P_Bar(betaVec(2:end),vVec,vCom,K_k,M);
        XkPlus = XkMinus + K_k*(vCom);
    end
%     PkPlus = PkMinus - K_k*S_k*K_k';
%     XkPlus = XkMinus + K_k*(z-zHat);
    FilterEstimateCV.XkPlus = XkPlus;
    FilterEstimateCV.PkPlus = PkPlus;
    FilterEstimateCV.XkPlusVec(:,loopNum+1) = XkPlus;
    FilterEstimateCV.PkVec(loopNum+1).Plus = PkPlus;
    LambdaCV = Calc_LH(vVec,M,S_k,cNZ,ParPDAF,nZ,transWave); % calcultaiton of likelyhood
end

function [FilterEstimateCT,LambdaCT,M] = UKF_Update_CT(FilterEstimateCT,FilterMat,z,loopNum,transWave,ParPDAF,ClutterStruct,zHat,ZPoints)
    XMinusPts = FilterEstimateCT.XMinusPts;
    XkMinus   = FilterEstimateCT.XkMinus;
    PkMinus   = FilterEstimateCT.PkMinus;
    n = length(XkMinus);
    alpha = 10^-3;
    kappa = 0;
    beta = 2;
    lamda = alpha^2*(n+kappa) - n;    

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for i=1:2*n
        wPts(i+1)       = 1/(2*(n+lamda));
    end
    
%     ZPoints             = Trans_H(XMinusPts,transWave);
%     zHat                = ZPoints*wPts';
    if transWave == "C"
        estRangeRate        = zHat(end);
        FilterEstimateCT.estRangeRate(loopNum) = estRangeRate;
    end    
    m1 = size(XMinusPts,1);
    m2 = size(zHat,1);
    Pzz = zeros(m2);
    for i = 0:2*n
        if i == 0
            Pzz = Pzz...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(ZPoints(:,i+1)- zHat)*(ZPoints(:,i+1)- zHat)';
        else
            Pzz = Pzz + wPts(i+1)...
                      *(ZPoints(:,i+1)- zHat)*(ZPoints(:,i+1)- zHat)';
        end
        
    end
    
    Pxz = zeros(m1,m2);
    for i = 0:2*n
        if i == 0
            Pxz = Pxz...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(XMinusPts(:,i+1)-XkMinus)*(ZPoints(:,i+1)- zHat)';
        else
            Pxz = Pxz + wPts(i+1)...
                      *(XMinusPts(:,i+1)- XkMinus)*(ZPoints(:,i+1)- zHat)';
        end        
    end    
    
    switch transWave
        case "C"
            S_k = Pzz + FilterMat.R.Cw;
            K_k = Pxz/(S_k);
            nZ = 3;
            cNZ = pi^(nZ/2)/gamma(nZ/2 + 1);
            ClutterMat = [ClutterStruct(loopNum).rhoClutter';ClutterStruct(loopNum).thetaClutter';ClutterStruct(loopNum).rangeRateClutter'];
            ClutterMat = [ClutterMat z];
        case "F"
            S_k = Pzz + FilterMat.R.Fm;
            K_k = Pxz/(S_k);
            nZ = 2;
            cNZ = pi^(nZ/2)/gamma(nZ/2 + 1);
            ClutterMat = [ClutterStruct(loopNum).rhoClutter';ClutterStruct(loopNum).thetaClutter'];
            ClutterMat = [ClutterMat z];
    end 
    [zGateVec,M] = Gating(ClutterMat,zHat,S_k,ParPDAF);
    vVec = zGateVec - repmat(zHat,1,size(zGateVec,2));
    betaVec = Data_Assoc(vVec,M,S_k,cNZ,ParPDAF,nZ,transWave);   
    if M == 0
        PkPlus = PkMinus;
        XkPlus = XkMinus;
        vVec = zeros(length(z),1);%vVec*betaVec(2:end);
    else
        vCom = vVec*betaVec(2:end);
        PkPlus = betaVec(1)*PkMinus + (1-betaVec(1))*(PkMinus - K_k*S_k*K_k')...
                 + P_Bar(betaVec(2:end),vVec,vCom,K_k,M);
        XkPlus = XkMinus + K_k*(vCom);
    end 
%     PkPlus = PkMinus - K_k*S_k*K_k';
%     XkPlus = XkMinus + K_k*(z-zHat);
    FilterEstimateCT.XkPlus = XkPlus;
    FilterEstimateCT.PkPlus = PkPlus;
    FilterEstimateCT.XkPlusVec(:,loopNum+1) = XkPlus;
    FilterEstimateCT.PkVec(loopNum+1).Plus = PkPlus;
    LambdaCT = Calc_LH(vVec,M,S_k,cNZ,ParPDAF,nZ,transWave);

end

function output = Trans_H(input,transWave)

    
    [~,col] = size(input);    
    switch transWave        
        case "C"
           output = zeros(3,col);
           for i = 1: col              
               output(1,i) = sqrt(input(1,i)^2 + input(2,i)^2);
               output(2,i) = atan2(input(2,i),input(1,i));
               output(3,i) = (input(1,i)*input(3,i) + input(2,i)*input(4,i))...
                             /sqrt(input(1,i)^2 + input(2,i)^2);
           
           end
        case "F"
            output = zeros(2,col);
            for i = 1: col              
               output(1,i) = sqrt(input(1,i)^2 + input(2,i)^2);
               output(2,i) = atan2(input(2,i),input(1,i));     
           end
    end
    
end

function Lambda = Calc_LH(vVec,M,S_k,cNZ,ParPDAF,nZ,transWave)
    switch transWave
        case "C"
           Pg = ParPDAF.Pg.cw; 
        case "F"
           Pg = ParPDAF.Pg.fm;
    end    
    Pd = ParPDAF.Pd;
%     if M == 0
%         e =((1/(((2*pi)^(nZ/2)) * sqrt(det(S_k)))) * exp(-0.5*vVec'*(S_k \ vVec)))/Pg;
%     else
        e = zeros(1,M);
        for i = 1:M
            e(i) = ((1/(((2*pi)^(nZ/2)) * sqrt(det(S_k)))) * exp(-0.5*vVec(:,i)'*(S_k \ vVec(:,i))))/Pg;          % 3.4.3-12
        end
%     end
    
    gateGamma = ParPDAF.gateGamma;
    gateVolume = cNZ*sqrt(gateGamma^nZ * det(S_k));
    b = M*(1- Pd*Pg)/(Pd*Pg*gateVolume);
    if M == 0
%         Lambda = (b + sum(e))* (Pd*Pg/ParPDAF.clutterDensity) * gateVolume^(-M);
%         Lambda = (1/(((2*pi)^(nZ/2)) * sqrt(det(S_k)))) * exp(-0.5*vVec'*(S_k \ vVec));
        Lambda = (1- Pg)* ParPDAF.clutterDensity;
%         Lambda = (1- Pg);
    else
%         Lambda = (b + sum(e))* (Pd*Pg/M) * gateVolume^(-M+1);
        Lambda = (1- Pg)* ParPDAF.clutterDensity + sum(e);
%           Lambda = (1- Pg)* gateVolume^(-M) + gateVolume^(-M+1)*(Pg/M)* sum(e);
    end
end

function [zGateVec,M] = Gating(ClutterMat,zHat,S_k,ParPDAF)
    zGateVec = zeros(size(ClutterMat));
    M = 0;   
    gateGamma = ParPDAF.gateGamma;
    for i = 1:size(ClutterMat,2)       
        error = (ClutterMat(:,i)-zHat)'*(S_k \ (ClutterMat(:,i)-zHat));
        if error <= gateGamma
            M = M+1;
            zGateVec(:,M) = ClutterMat(:,i);            
        end
    end
    zGateVec = zGateVec(:,1:M);
 
end
% 
function betaVec = Data_Assoc(vVec,M,S_k,cNZ,ParPDAF,nZ,transWaveEstUKF)

    e = zeros(1,M);
    for i = 1:M
        e(i) = exp(-0.5*vVec(:,i)'*(S_k \ vVec(:,i)));          % 3.4.3-12
    end
    gateGamma = ParPDAF.gateGamma;
    switch transWaveEstUKF
        case "C"
           Pg = ParPDAF.Pg.cw; 
        case "F"
           Pg = ParPDAF.Pg.fm;
    end
    b = (2*pi/gateGamma)^(nZ/2)*M/cNZ*(1-ParPDAF.Pd*Pg)/ParPDAF.Pd;      % 3.4.3-14

    betaZero = b/(b + sum(e));                        %3.4.3-11b
    betaRest = e ./ (b+sum(e));
    betaVec = [betaZero;betaRest'];
end

function pBarMat = P_Bar(partBetaVec,vVec,vCom,K_k,M)

    sum = zeros(size(vCom,1));
    for i = 1:M
        sum = sum + partBetaVec(i)*(vVec(:,i)* vVec(:,i)');
    end
    pBarMat = K_k*(sum - vCom*vCom')*K_k';
end

% e(i) = exp(-0.5*vVec(:,i)'*(S_k \ vVec(:,i))); 
% b = (2*pi/gateGamma)^(nZ/2)*M/cNZ*(1-ParPDAF.Pd*Pg)/ParPDAF.Pd;