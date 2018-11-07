function transWave = Planner_True(ParPath,impInstVec)

    trueTargetDirTheta = True_Target_Direction(ParPath,impInstVec);
    srcWaves = ['F','C']; 
    transWave = zeros(size(impInstVec,1),1);
    transWave = string(transWave);
    ChanceNodeMat = zeros(size(impInstVec,1),size(srcWaves,2));
    for i = 1: length(srcWaves)
        ChanceNodeMat(:,i) = True_Eval(trueTargetDirTheta',srcWaves(i));
    end
    transWave(ChanceNodeMat(:,2) >= ChanceNodeMat(:,1)) = "C";
    transWave(ChanceNodeMat(:,1) > ChanceNodeMat(:,2)) = "F";    
end

function trueTargetDirTheta = True_Target_Direction(ParPath,impInstVec)
    TarPathMat = ParPath.TarPathMat;
    if impInstVec(end) == size(TarPathMat,1)
        impInstVec(end) = impInstVec(end) - 1;
    end
    xVec                = [TarPathMat(impInstVec,1)';TarPathMat(impInstVec+1,1)'];
    yVec                = [TarPathMat(impInstVec,2)';TarPathMat(impInstVec+1,2)'];
    velDirSlope         = diff(yVec)./diff(xVec);
    slopeWithOrigin     = (yVec(1,:)- ParPath.selfY)./(xVec(1,:)- ParPath.selfX);
    trueTargetDirTheta  = atan(abs((velDirSlope-slopeWithOrigin)...
                         ./(1+ (velDirSlope.*slopeWithOrigin))));
end
    
function scoreVec = True_Eval(trueTargetDirTheta,srcWave)

    switch srcWave
        case 'F'
            fm = 1;
            cw = 0;
        case 'C'
            fm = 0;
            cw = 1;
    end

    scoreVec = (((cos(trueTargetDirTheta)).^2) * cw...
                + ((sin(trueTargetDirTheta)).^2)* fm);
end