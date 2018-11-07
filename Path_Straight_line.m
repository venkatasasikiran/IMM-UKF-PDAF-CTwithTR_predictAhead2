function [varargout] = Path_Straight_line(ParPath,ParGen)

    f               = ParGen.f;   
    [~,tarR]        = cart2pol(ParPath.initTarX,ParPath.initTarY);
    tarTheta        = deg2rad(360*rand);
    [tarX,tarY]     = pol2cart(tarTheta,tarR);
    T               = 0:1/f:ParPath.tLimit;    
    TarPathMat      = zeros(size(T,2),2);
    TarPathMat(1,1) = tarX;
    TarPathMat(1,2) = tarY;
    
    destR           = tarR;
    destTheta       = deg2rad(360*rand);
    [destX, destY]  = pol2cart(destTheta,destR);
    dirTheta      = atan2(destY - tarY, destX - tarX);
    ParPath.initTarX = tarX;
    ParPath.initTarY = tarY;
    ParPath.initVelDir = dirTheta;
    
    isChangeSpeed = false;
    isChangeTest = true;
    changeTime   = 90*f;
    changeTimer = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initState     = [tarX;tarY;ParPath.initTarV*cos(dirTheta);...
                    ParPath.initTarV*sin(dirTheta)];
    TransCV       = [1    0    1/f    0;...
                     0    1     0    1/f;...
                     0    0     1     0;...
                     0    0     0     1];
               
    
    %% declaration of required variables
    for i = 2:length(T)
        if round(T(i)*f) == changeTime && isChangeTest 
            isChangeSpeed     = true;
            isChangeTest      = false;
            acceleration = rand*4;
        end
        if isChangeSpeed
            changeTimer = changeTimer + 1;
            speed       = sqrt(initState(3)^2 + initState(4)^2) + acceleration*1/f;
            initState   = [initState(1:2);speed*cos(dirTheta);speed*sin(dirTheta)];
            tarState    = TransCV*initState;          
            tarX                 = tarState(1);
            tarY                 = tarState(2);
            TarPathMat(i,1)      = tarX;
            TarPathMat(i,2)      = tarY;
            initState            = tarState;
            
            if changeTimer >= 10*f
                changeTimer = 0;
                isChangeTest = true;
                isChangeSpeed = false;
                changeTime = i + 90*f;
            end
            
        else        
            tarState             = TransCV*initState;
            tarX                 = tarState(1);
            tarY                 = tarState(2);
            TarPathMat(i,1)      = tarX;
            TarPathMat(i,2)      = tarY;
            initState            = tarState;
        end
    end
    [tarThetaVec,tarRhoVec] = cart2pol(TarPathMat(:,1),TarPathMat(:,2));
    TarPathMatInPol = [tarThetaVec tarRhoVec];    
    ParPath.TarPathMat = TarPathMat;
    ParPath.TarPathMatInPol = TarPathMatInPol;
    ParPath.TarMaxRho = max(TarPathMatInPol(:,2));
    
    varargout{1} = ParPath;    
end

