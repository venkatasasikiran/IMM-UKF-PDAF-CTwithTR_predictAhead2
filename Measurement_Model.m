function [Measurements,ParPath] = Measurement_Model(ParPath,impInstVec,ParStdDev)

    TruePosMat = ParPath.TarPathMat(impInstVec,:);
    TruePosMatPol = ParPath.TarPathMatInPol(impInstVec,:);
    TrueVelocities = True_Vel(ParPath,impInstVec);
    trueRangeRate = (TruePosMat(:,1).*TrueVelocities(:,1)...
                    +TruePosMat(:,2).*TrueVelocities(:,2))...
                    ./sqrt(TruePosMat(:,1).^2+TruePosMat(:,2).^2);
    
    MeasRange.cw = TruePosMatPol(:,2) + ParStdDev.Range.cw*randn(size(TruePosMat,1),1);
    MeasRange.fm = TruePosMatPol(:,2) + ParStdDev.Range.fm*randn(size(TruePosMat,1),1);  

    MeasTheta = TruePosMatPol(:,1) + ParStdDev.theta*randn(size(TruePosMat,1),1);
    
    MeasRangeRate.cw = trueRangeRate + ParStdDev.RangeRate.cw*randn(size(trueRangeRate));
    
    
    Measurements.MeasRange = MeasRange;
    Measurements.MeasTheta = MeasTheta;
    Measurements.MeasRangeRate = MeasRangeRate;
    
    ParPath.TrueVelocities = TrueVelocities;
    
%     Measurements.maxStdDevRange = max(ParStdDev.Range.cw,ParStdDev.Range.fm);
%     Measurements.maxStdDevTheta = max(ParStdDev.theta);
%     Measurements.maxStdDevRangeRate = max(ParStdDev.RangeRate.cw);
end


function TrueVelocities = True_Vel(ParPath,impInstVec)
    if impInstVec(end) == size(ParPath.TarPathMat,1)
        impInstVec(end) = impInstVec(end) - 1;
    end
    xVec = [ParPath.TarPathMat(impInstVec,1)';ParPath.TarPathMat(impInstVec+1,1)'];
    yVec = [ParPath.TarPathMat(impInstVec,2)';ParPath.TarPathMat(impInstVec+1,2)'];
    velDir = atan2(diff(yVec),diff(xVec));
    TrueVelocities = [ParPath.initTarV*cos(velDir') ParPath.initTarV*sin(velDir')];
end