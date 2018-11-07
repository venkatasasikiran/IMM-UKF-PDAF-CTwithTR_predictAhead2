function ParStdDev = Struct_Param_Std_Dev(varargin)

    if nargin > 0
        switch nargin
            case 1
                Range.cw = varargin{1}(1);
                Range.fm = varargin{1}(2);

                theta = .01;               

                RangeRate.cw = 1;
                
                
            case 2
                Range.cw = varargin{1}(1);
                Range.fm = varargin{1}(2);

                theta = varargin{2};               

                RangeRate.cw = 1;
            case 3
                Range.cw = varargin{1}(1);
                Range.fm = varargin{1}(2);

                theta = varargin{2};               

                RangeRate.cw = varargin{3};
        end
        
    else
        Range.cw = 50;
        Range.fm = 30;
        
        theta = .01;
        
        RangeRate.cw = 1;        
    end
    
    ParStdDev.Range = Range;
    ParStdDev.theta = theta;
    ParStdDev.RangeRate = RangeRate;

end