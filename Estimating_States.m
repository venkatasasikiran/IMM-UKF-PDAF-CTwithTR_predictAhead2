function [FilterEstimateCV,FilterEstimateCTpos,FilterEstimateCTneg] = Estimating_States(FilterEstimateCV,FilterEstimateCTpos,...
                                                                        FilterEstimateCTneg,FilterMat,loopNum)
    FilterEstimateCV = UKF_Estimate_CV(FilterEstimateCV,FilterMat,loopNum);
    FilterEstimateCTpos = UKF_Estimate_CT(FilterEstimateCTpos,FilterMat,loopNum,"P");
    FilterEstimateCTneg = UKF_Estimate_CT(FilterEstimateCTneg,FilterMat,loopNum,"N");
end

function FilterEstimateCV = UKF_Estimate_CV(FilterEstimateCV,FilterMat,loopNum)

    PkPlus = FilterEstimateCV.PkPlusMix;
    XkPlus = FilterEstimateCV.XkPlusMix;
    n = length(XkPlus);
    alpha = 10^-3;
    kappa = 0;
    beta = 2;
    lamda = alpha^2*(n+kappa) - n;    

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    XPlusPts    = zeros(n,noPoints); % samples

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for i=1:2*n
        wPts(i+1)       = 1/(2*(n+lamda));
    end
    % Calculate sigma points
    %
    [U,Eigen,~]               = svd(PkPlus);
    Cxx                       = U*sqrt(Eigen);
    c                         = sqrt(n+lamda);
    XPlusPts(:,1)             = XkPlus;
    XPlusPts(:,2:n+1)         = XkPlus*ones(1,n)+c*Cxx; 
    XPlusPts(:,n+2:2*n+1)     = XkPlus*ones(1,n)-c*Cxx;
    
    % Transform the sample points
    %
    XMinusPts = Trans_F(XPlusPts,FilterMat.CV.F);
    m1 = size(XMinusPts,1);
    % Calculate the X_k_minus and P_k_minus
    %
    XkMinus = XMinusPts*wPts';
    PkMinus = zeros(m1);
    for i = 0: 2*n
        if i == 0
            PkMinus = PkMinus...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(XMinusPts(:,i+1)-XkMinus)*(XMinusPts(:,i+1) - XkMinus)';
        else
            PkMinus = PkMinus + wPts(i+1)...
                      *(XMinusPts(:,i+1)-XkMinus)*(XMinusPts(:,i+1) - XkMinus)';
        end
    end
    PkMinus = PkMinus + FilterMat.CV.Q;
    FilterEstimateCV.XkMinus = XkMinus;
    FilterEstimateCV.PkMinus = PkMinus;
    FilterEstimateCV.XkMinusVec(:,loopNum) = XkMinus;
    FilterEstimateCV.PkVec(loopNum).Minus = PkMinus;
    FilterEstimateCV.XMinusPts = XMinusPts;    
end

function FilterEstimateCT = UKF_Estimate_CT(FilterEstimateCT,FilterMat,loopNum,turnCase)
    
    PkPlus = FilterEstimateCT.PkPlusMix;
    XkPlus = FilterEstimateCT.XkPlusMix;
    n = length(XkPlus);
    alpha = 10^-3;
    kappa = 0;
    beta = 2;
    lamda = alpha^2*(n+kappa) - n; 

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    XPlusPts    = zeros(n,noPoints); % samples

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for j=1:2*n
        wPts(j+1)       = 1/(2*(n+lamda));
    end
    % Calculate sigma points
    %
    [U,Eigen,~]               = svd(PkPlus);
    Cxx                       = U*sqrt(Eigen);
    c                         = sqrt(n+lamda);
    XPlusPts(:,1)             = XkPlus;
    XPlusPts(:,2:n+1)         = XkPlus*ones(1,n)+c*Cxx; 
    XPlusPts(:,n+2:2*n+1)     = XkPlus*ones(1,n)-c*Cxx;
    
    if turnCase == "P"
        XMinusPts = Trans_F(XPlusPts,FilterMat.CT.Fpos);
    else
        XMinusPts = Trans_F(XPlusPts,FilterMat.CT.Fneg);
    end
    
    m1 = size(XMinusPts,1);
    % Calculate the X_k_minus and P_k_minus
    %
    XkMinus = XMinusPts*wPts';
    PkMinus = zeros(m1);
  
    for i1 = 0: 2*n
        if i1 == 0
            PkMinus = PkMinus...
                      + (wPts(i1+1)+(1-alpha^2+beta))...
                      *(XMinusPts(:,i1+1)-XkMinus)*(XMinusPts(:,i1+1) - XkMinus)';
        else
            PkMinus = PkMinus + wPts(i1+1)...
                      *(XMinusPts(:,i1+1)-XkMinus)*(XMinusPts(:,i1+1) - XkMinus)';
        end
    end
    PkMinus = PkMinus + FilterMat.CT.Q;
    FilterEstimateCT.XkMinus = XkMinus;
    FilterEstimateCT.PkMinus = PkMinus;
    FilterEstimateCT.XkMinusVec(:,loopNum) = XkMinus;
    FilterEstimateCT.PkVec(loopNum).Minus = PkMinus;
    FilterEstimateCT.XMinusPts = XMinusPts;   
end

function output = Trans_F(input,F)
    
        [Row,Col] = size(input);
        output = zeros(Row, Col);    
        for k = 1:Col
            output(:,k) = F*input(:,k);        
        end
    end
