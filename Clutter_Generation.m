 function ClutterStruct = Clutter_Generation(TarMaxRho,clutterDensity,impInstVec)
    
    areaSquare = (2*TarMaxRho)^2;
    ptsGen = ceil(areaSquare*clutterDensity);
    n = size(impInstVec,1);
    xClutter = -TarMaxRho + 2*TarMaxRho*rand(ptsGen,n);          %This is uniform distributiom
    yClutter = -TarMaxRho + 2*TarMaxRho*rand(ptsGen,n);          %This is uniform distributiom
    
    for i = n:-1:1
        [thetaClutter,rhoClutter] = cart2pol(xClutter(:,i),yClutter(:,i));
        indx = find(rhoClutter<=TarMaxRho);
        ClutterStruct(i).thetaClutter = thetaClutter(indx,:);
        ClutterStruct(i).rhoClutter = rhoClutter(indx,:);
        noClutterPts = length(indx);        
        ClutterStruct(i).rangeRateClutter = 5*randn(noClutterPts,1);
    end
 end
