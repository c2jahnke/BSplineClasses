classdef hbSplBasML < handle
    properties %(GetAccess = public, SetAccess = private)    
        levelBas = [];
        basis0 = [];
        level = 0;
    end
     properties (Hidden = true, SetAccess = private)
        foo = 1
    end
     methods (Access = public)
        function obj = hbSplBasML(a,b,p,N,resol,level)
            % constructor for class
            % a, b, p, knotspan, resol
            if nargin >0
                assert(level >= 1, "The number of levels has to be a positive integer.");
            for k = 1 : level
                obj.levelBas{1,k} = hbSplBas(a,b,p,N*2^(k-1),resol/2^(k-1));
            end
                obj.basis0 = obj.levelBas{1,1};
                obj.level = level;
            end
              
        end
        
        
        
     end
end