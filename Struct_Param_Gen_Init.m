function ParGen = Struct_Param_Gen_Init(varargin)

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
            ParGen.velSound = 1500; % m/s5
            ParGen.f = 1000; % sampling frequency
            ParGen.clutterDensity = varargin{1};
            
        case 2
            ParGen.velSound = 1500; % m/s5
            ParGen.f = varargin{2}; % sampling frequency
            ParGen.clutterDensity = varargin{1};
    end
    
else
    ParGen.velSound = 1500; % m/s5
    ParGen.f = 1000; % sampling frequency
    ParGen.clutterDensity = 1/10000;
end
