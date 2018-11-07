function ParPath = Struct_Param_Path_Init(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Init_Parameters - initializes different parameters of 
% the algorithm
% Input:
%   -
% Output:
%   Par          - parameter structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 0
    switch nargin
        case 1
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(10);
            ParPath.yawDev     = 90;
            ParPath.initTarX = 4000;
            ParPath.initTarY = 3000;
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = 500;
            ParPath.trackLimit = 300;
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
            
        case 2
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(varargin{2});
            ParPath.yawDev     = 90;
            ParPath.initTarX = 4000;
            ParPath.initTarY = 3000;
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = 500;
            ParPath.trackLimit = 300;
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
            
        case 3
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(varargin{2});
            ParPath.yawDev     = varargin{3};
            ParPath.initTarX = 4000;
            ParPath.initTarY = 3000;
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = 500;
            ParPath.trackLimit = 300;
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
            
        case 4
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(varargin{2});
            ParPath.yawDev     = varargin{3};
            ParPath.initTarX = varargin{4}(1);
            ParPath.initTarY = varargin{4}(2);
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = 500;
            ParPath.trackLimit = 300;
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
        case 5
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(varargin{2});
            ParPath.yawDev     = varargin{3};
            ParPath.initTarX = varargin{4}(1);
            ParPath.initTarY = varargin{4}(2);
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = varargin{5};
            ParPath.trackLimit = 300;
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
            
        case 6
            ParPath.isManeuver = varargin{1};
            ParPath.yawDevRate = deg2rad(varargin{2});
            ParPath.yawDev     = varargin{3};
            ParPath.initTarX = varargin{4}(1);
            ParPath.initTarY = varargin{4}(2);
            [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
            ParPath.initTarV = 30; % the velocity is in m/sec

            ParPath.selfX = 0;
            ParPath.selfY = 0;
            [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

            ParPath.tLimit = varargin{5};
            ParPath.trackLimit = varargin{6};
            ParPath.trackLimitMax = 7500;

            ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );

    end
else
    ParPath.isManeuver = true;
    ParPath.yawDevRate = deg2rad(10);
    ParPath.yawDev     = 90;
    ParPath.initTarX = 4000;
    ParPath.initTarY = 3000;
    [ParPath.initTarTheta,ParPath.initTarRho] = cart2pol(ParPath.initTarX, ParPath.initTarY);
    ParPath.initTarV = 30; % the velocity is in m/sec

    ParPath.selfX = 0;
    ParPath.selfY = 0;
    [ParPath.selfTheta,ParPath.selfRho] = cart2pol(ParPath.selfY, ParPath.selfY);

    ParPath.tLimit = 500;
    ParPath.trackLimit = 300;
    ParPath.trackLimitMax = 7500;

    ParPath.initVelDir = atan2(ParPath.selfY - ParPath.initTarY,ParPath.selfX - ParPath.initTarX );
end
