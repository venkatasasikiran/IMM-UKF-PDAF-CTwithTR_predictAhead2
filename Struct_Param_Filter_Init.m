function FilterEstimate = Struct_Param_Filter_Init(ParPath,stateSize,n,motionModel)
    if motionModel == "CV"
        XkPlus = [ParPath.initTarX;ParPath.initTarY;...
                  ParPath.initTarV*cos(ParPath.initVelDir);...
                  ParPath.initTarV*sin(ParPath.initVelDir)];  

        XkPlusVec = zeros(stateSize,n+1);
        XkPlusVec(:,1) = XkPlus;

        PkPlus = 100*eye(stateSize);
        Pk(n+1).Plus = zeros(stateSize);
        Pk(1).Plus = PkPlus;

        XkMinusVec = zeros(stateSize,n);
        Pk(n).Minus = zeros(stateSize);

        FilterEstimate.XkPlus = XkPlus;
        FilterEstimate.PkPlus = PkPlus;
        FilterEstimate.XkPlusVec = XkPlusVec;    
        FilterEstimate.XkMinusVec = XkMinusVec;
        FilterEstimate.PkVec = Pk;

        FilterEstimate.estRangeRate = zeros(n,1);
    else
        XkPlus = [ParPath.initTarX;ParPath.initTarY;...
                  ParPath.initTarV*cos(ParPath.initVelDir);...
                  ParPath.initTarV*sin(ParPath.initVelDir)];  

        XkPlusVec = zeros(stateSize,n+1);
        XkPlusVec(:,1) = XkPlus;

        PkPlus = 100*eye(stateSize); 
        Pk(n+1).Plus = zeros(stateSize);
        Pk(1).Plus = PkPlus;

        XkMinusVec = zeros(stateSize,n);
        Pk(n).Minus = zeros(stateSize);

        FilterEstimate.XkPlus = XkPlus;
        FilterEstimate.PkPlus = PkPlus;
        FilterEstimate.XkPlusVec = XkPlusVec;    
        FilterEstimate.XkMinusVec = XkMinusVec;
        FilterEstimate.PkVec = Pk;

        FilterEstimate.estRangeRate = zeros(n,1);
    end
end
