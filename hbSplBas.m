classdef hbSplBas < bSplBas
     properties
        activeKnots = 0; % for basis refinement
        activeIndex = 0;% function index
     end
     methods (Access = public)
         function obj = hbSplBas(a,b,p,N,resol)
            % constructor for class
            % a, b, p, knotspan, resol
            if nargin >0
                obj.a = a;
                obj.b = b;
                obj.N = N;
                obj.p = p;
                obj.resol = resol;
                
            end
                obj.knotspan = (obj.b-obj.a)/obj.N;
                obj.knotVector = ConstrKnotVector(obj.a,obj.b,obj.knotspan,obj.p);
                obj.plotVector = [obj.knotVector(1):obj.resol:obj.knotVector(end)];
                obj.sP = size(obj.plotVector,2);
                obj.m = size(obj.knotVector,2);
                obj.n = obj.m - obj.p - 1;
                obj.activeKnots = obj.getAllKnots;
                obj.activeIndex = [0 : obj.n-1]; % function index
            function knotVector = ConstrKnotVector(a,b,h1,p)
                % simple constructor for knot vector
                temp = (b-a)/h1;
                assert( rem(temp,1) == 0,'Use values a,b,h1 such that b-a/h1 is an integer.')
                m = temp + 2*p+1;
                knotVector = zeros(1,m);
                knotVector(1:p) = a;
                knotVector(p+1:m-p) = a:h1:b;
                knotVector(m-p+1:end) = b;
            end
         end
     end
end