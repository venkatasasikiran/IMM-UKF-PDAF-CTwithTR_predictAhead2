function [Timings] = Timing_Model(ParPath,ParGen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Init_Parameters - initializes different parameters of 
% the algorithm
% Input:
%   -
% Output:
%   Par          - parameter structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impTimePeriod = 2;% in sec
    noOfImpInst = floor(size(ParPath.TarPathMat,1)/(impTimePeriod*ParGen.f));
    if rem(size(ParPath.TarPathMat,1),impTimePeriod*ParGen.f) == 0        
        impInstVec = (impTimePeriod*ParGen.f*(1:noOfImpInst))+1;
        impInstVec(end) = impInstVec(end) - 1;
    else
        impInstVec = (impTimePeriod*ParGen.f*(1:noOfImpInst))+1;
    end    
    Timings.impInstVec = impInstVec';
    
        
        
