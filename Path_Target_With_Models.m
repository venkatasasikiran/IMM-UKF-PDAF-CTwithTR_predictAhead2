function [varargout] = Path_Target_With_Models(ParPath,ParGen)

    f               = ParGen.f;
    selfX           = ParPath.selfX;
    selfY           = ParPath.selfY;
    tarX            = ParPath.initTarX;
    tarY            = ParPath.initTarY;
    yawDevRate      = ParPath.yawDevRate;
    T               = 0:1/f:ParPath.tLimit;
    devVecIndxLimit = ceil(ParPath.tLimit/5);
    TarPathMat      = zeros(size(T,2),2);
    TarPathMat(1,1) = ParPath.initTarX;
    TarPathMat(1,2) = ParPath.initTarY;
    
    
    %% declaration of required variables
%     devTime       = Rand_Gen;
    devTime       = Gen_Dev_Time(ParPath);      
    devTimeVec    = zeros(devVecIndxLimit,2);
    if ParPath.isManeuver
        devTimeVec(1,1) = devTime*f +1;
    end
    devVecIndx    = 2;
    hapnDev       = 0;
    dirTheta      = atan2(selfY - tarY, selfX - tarX);
    isDevChoose   = true;                                 %Parameter used to alternate betweeen the two variants of deviation thetas 
    %isTurnParam   = false;                                %This parameter is used to compensate the turn radius 
    isPerTurn     = false;                                  % Perform turn parameter 
    isGetDevTheta = false; 
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initState     = [ParPath.initTarX;ParPath.initTarY;ParPath.initTarV*cos(dirTheta);...
                    ParPath.initTarV*sin(dirTheta)];
    TransCV       = [1    0    1/f    0;...
                     0    1     0    1/f;...
                     0    0     1     0;...
                     0    0     0     1];
    TransCtPos    = [1         0         sin(yawDevRate/f)/yawDevRate         (cos(yawDevRate/f)- 1)/yawDevRate;...
                     0         1      (1 - cos(yawDevRate/f))/yawDevRate       sin(yawDevRate/f)/yawDevRate;...
                     0         0          cos(yawDevRate/f)           -sin(yawDevRate/f);...
                     0         0          sin(yawDevRate/f)            cos(yawDevRate/f)];
                 
    TransCtNeg    = [1         0         sin(-yawDevRate/f)/-yawDevRate         (cos(-yawDevRate/f)- 1)/-yawDevRate;...
                     0         1      (1 - cos(-yawDevRate/f))/-yawDevRate       sin(-yawDevRate/f)/-yawDevRate;...
                     0         0          cos(-yawDevRate/f)           -sin(-yawDevRate/f);...
                     0         0          sin(-yawDevRate/f)            cos(-yawDevRate/f)];
    
    %% declaration of required variables
    for i = 2:length(T)       
        if round(T(i)*f) == round(devTime*f) && ParPath.isManeuver 
            isPerTurn     = true;
            isGetDevTheta = true;                   
        end
        if isGetDevTheta        
            if isDevChoose
                % the below function gives the turn angle as per gaussian
                % distribution surrounded around 90 deg
                devTheta      = Dev_Theta(ParPath.yawDev);
                isDevChoose   = false;
                isGetDevTheta = false;            
            else
                % This is the another deviation theta that can be chosen. 
                %dev_theta is obtained by difference of what direction it is ..
                %moving right now and what direction it should move.
                devTheta      = atan2(selfY - tarY, selfX - tarX)  - dirTheta;
                devTheta      = Theta_Adjust(devTheta);  % adjust the theta such that it is in range of (-pi,pi)
                isDevChoose   = true;
                isGetDevTheta = false;
                %isTurnParam   = true;                    % It is made true as because of turning radius the required deviation is not achieved  
            end
        end
        if isPerTurn
            targetDistance = sqrt((tarX - selfX)^2 + (tarY - selfY)^2);
            if targetDistance < ParPath.trackLimit || targetDistance > ParPath.trackLimitMax          
                TarPathMat = TarPathMat(1:i-1,:);          
                break;
            else
                if devTheta >=0
                    dirTheta = dirTheta+yawDevRate/f;
                    hapnDev  = hapnDev + yawDevRate/f;                    
                    tarState             = TransCtPos*initState;
                    tarX                 = tarState(1);
                    tarY                 = tarState(2);
                    TarPathMat(i,1)      = tarX;
                    TarPathMat(i,2)      = tarY;
                    initState            = tarState;
                    if hapnDev >= devTheta                        
                        %No turn_param. so all parameters are reset
                        % New deviation time is chosen
                        % Robust param is made true to check fro new turn time
                        isPerTurn                 = false;
                        devTimeVec(devVecIndx-1,2)= T(i)*f;
                        devTime                   = T(i) + Rand_Gen;
                        devTimeVec(devVecIndx,1)  = devTime*f;
                        devVecIndx                = devVecIndx + 1;
                        hapnDev                   = 0;                        
                    end
                else
                    dirTheta = dirTheta - yawDevRate/f;
                    hapnDev  = hapnDev - yawDevRate/f;                   
                    tarState             = TransCtNeg*initState;
                    tarX                 = tarState(1);
                    tarY                 = tarState(2);
                    TarPathMat(i,1)      = tarX;
                    TarPathMat(i,2)      = tarY;
                    initState            = tarState;
                    if hapnDev <= devTheta
                        isPerTurn                 = false;
                        devTimeVec(devVecIndx-1,2)= T(i)*f;
                        devTime                   = T(i) + Rand_Gen;
                        devTimeVec(devVecIndx,1)  = devTime*f;
                        devVecIndx                = devVecIndx + 1;
                        hapnDev                   = 0;
                        
                    end
                end                
            end
            
        else
            targetDistance = sqrt((tarX - selfX)^2 + (tarY - selfY)^2);
            if targetDistance < ParPath.trackLimit || targetDistance > ParPath.trackLimitMax          
                TarPathMat = TarPathMat(1:i-1,:);          
                break;
            else
                tarState             = TransCV*initState;
                tarX                 = tarState(1);
                tarY                 = tarState(2);
                TarPathMat(i,1)      = tarX;
                TarPathMat(i,2)      = tarY;
                initState            = tarState;
            end
        end
    end
    [tarThetaVec,tarRhoVec] = cart2pol(TarPathMat(:,1),TarPathMat(:,2));
    TarPathMatInPol = [tarThetaVec tarRhoVec];
    devTimeVec = devTimeVec(devTimeVec(:,1)~= 0,:);
    ParPath.TarPathMat = TarPathMat;
    ParPath.TarPathMatInPol = TarPathMatInPol;
    ParPath.TarMaxRho = max(TarPathMatInPol(:,2));
    varargout{1} = ParPath;
    varargout{2} = round(devTimeVec);    

end

function devTime = Rand_Gen()
    devTime = round(random('inversegaussian',60,60),2);
end

function devTime  = Gen_Dev_Time(ParPath)
    distStraight = ParPath.initTarRho - ParPath.trackLimit;
    timeStraight = distStraight/ParPath.initTarV;
    devTime = Rand_Gen;
    while (devTime >= timeStraight)
        devTime = Rand_Gen;
    end
end

function res = Dev_Theta(yawDev)
    res = deg2rad(10*randn+yawDev);  % it gives theta in radians which is centered around pi/2
    %res = deg2rad(90);

end

function devTheta = Theta_Adjust(devTheta)
    if devTheta > pi
       devTheta = devTheta - 2*pi;
    elseif devTheta < -pi
       devTheta = devTheta + 2*pi;
    end
end

% function Fmat = TransMat(varargin)
%     if varargin{1} == "CV"
%         f = varargin{2};
%         Fmat       =   [1    0    1/f    0;...
%                         0    1     0    1/f;...
%                         0    0     1     0;...
%                         0    0     0     1];
%     else
%         f = varargin{2};
%         w = varargin{3};
%         Fmat =      [1         0         sin(w/f)/w         (cos(w/f)- 1)/w;...
%                      0         1      (1 - cos(w/f))/w       sin(w/f)/w;...
%                      0         0          cos(w/f)           -sin(w/f);...
%                      0         0          sin(w/f)            cos(w/f)];
%     end
% end