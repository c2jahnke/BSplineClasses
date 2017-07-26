classdef bSplBas
    %% Class for BSplineBasis
    % so far contains: Basis generation
    % tablespan
    % derivative generation
    % Basis plotting
    % 
    properties %(GetAccess = public, SetAccess = private)    
        a = 0;
        b = 5;
        knotspan = 1;
        p = 2;
        knotVector = [0 0 1 2 3 4 5 5]   
        resol = 0.2;
        plotVector = 0;
        sP = 0;
        m = 0;
        n = 0;
    end
     properties(Hidden = true, SetAccess = private)
       kU = 2;      
    end
     methods (Access = public)
        function obj = bSplineBasis(a,b,p,knotspan,resol)
            % constructor for class
            if nargin >0
                obj.a = a;
                obj.b = b;
                obj.knotspan = knotspan;
                obj.p = p;
                obj.resol = resol;
                
            end
                obj.knotVector = ConstrKnotVector(obj.a,obj.b,obj.knotspan,obj.p);
                obj.plotVector = [obj.knotVector(1):obj.resol:obj.knotVector(end)];
                obj.sP = size(obj.plotVector,2);
                obj.m = size(obj.knotVector,2);
                obj.n = obj.m - obj.p - 1;
            function knotVector = ConstrKnotVector(a,b,h1,p)
                % simple constructor for knot vector
                temp = (b-a)/h1;
                assert( isinteger(temp) == 0,'Use values a,b,h1 such that b-a/h1 is an integer.')
                m = temp + 2*p+1;
                knotVector = zeros(1,m);
                knotVector(1:p) = a;
                knotVector(p+1:m-p) = a:h1:b;
                knotVector(m-p+1:end) = b;
            end
              
        end
        
        
        
        function C = generBasis(obj)
        C = zeros(obj.sP,obj.n);
        for ll = 1:obj.n
            for ks = 1 : obj.sP
                C(ks,ll) = OneBasisFun(obj.p,obj.m,obj.knotVector,ll-1,obj.plotVector(ks));
            end
        end
        end
        
        function plotBasisStruct(obj,C)
        %C = zeros(obj.sP,obj.n);
        %C = generBasis(obj);
        assert(isequal(size(C),[obj.sP,obj.n]),'Warning: Basis structure may not be plottable.');
        for ll = 1:obj.n
            plot(obj.plotVector,C(:,ll));
            hold all
        end
        hold off
        end
        
        function tableSpan = lookUpTableSpan(obj)
        %--------------------------------------------------------------
        % Input:
        % knotVector : knot vector
        % plotVector : plot vector
        % m          : size of knot vector
        % sP          : size of plot vector
        %Output
        % tableSpan : row vector of same dimension as plotVector
        %--------------------------------------------------------------

        tableSpan = zeros(1,obj.sP);
        plotIndex = 1;
        leftIndex = 1; % left boundary of knot span
        for knotIndex = 2:obj.m
            rightValue = obj.knotVector(knotIndex);
            if(obj.knotVector(leftIndex) == rightValue)
                continue
            end
            leftIndex = knotIndex -1;
            while(obj.plotVector(plotIndex) < rightValue)
                tableSpan(plotIndex) = leftIndex;
                plotIndex = plotIndex +1;
            end
        end

        tableSpan(plotIndex) = leftIndex -1;
        tableSpan = tableSpan -ones(size(tableSpan));
        end
        
        function C = gener1Deriv(obj)
        %--------------------------------------------------------------
        % Warning: partly own idea, test carefully!
        % Input:
        % knotVector : knot vector
        % plotVector : plot vector
        % tableSpan  : lookUpTableSpan
        % p          : polynomial degree
        % sP         : size of plot vector
        %Output
        % tableSpan : row vector of same dimension as plotVector
        %--------------------------------------------------------------

            %(sP,p,tableSpan,plotVector,U,nDer)
    
            %n = size(U,2) - p -1;
            C = zeros(obj.sP,obj.n);
            tableSpan = lookUpTableSpan(obj);
            % use getPlotMatrix line 1040 ff
            % use startX to switch index when entering new span index
            for i = 1 : obj.sP
                startX = tableSpan(i) - tableSpan(1) +1;
                ders = DersBasisFuns(tableSpan(i),obj.plotVector(i),obj.p,1,obj.knotVector);
                for j = 1 : obj.n % change to n
                    if (j-tableSpan(i) + obj.p -1 ) < obj.p+1 && tableSpan(i) +1 - j  < obj.p+1
                        tmp = mod(j - startX, (obj.p+1))+1;  % temporary solution
                        C(i,j) = ders(2,tmp);
                    end
                end
            end
        end
     end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.2
    % evaluate basis function
    function N = BasisFun(i,p,u,U) % careful: differnt input parameter order (u,p switched)
    %--------------------------------------------------------------
    %function N = BasisFun(i,p,u,U)
    % NURBS-Book (algorithm A2.2)
    % evalute nonzero basis functions
    %INPUT:
    % i          : current knotspan
    % u          : evaluation point
    % p          : degree of the basis functions
    % U          : knot vector (row vector)
    %OUTPUT:
    % N          : row vector (dim p+1)
    %              values of the basis function N_(i-p) ... N_(i)
    %              at the evaluation point u
    %--------------------------------------------------------------
    N=zeros(1,p+1);
    N(1)=1;
    left=zeros(1,p+1);
    right=zeros(1,p+1);
    for j=1:p
        left(j+1) = u-U(i+1-j+1); %% why j+1 and not j-1?
        right(j+1) = U(i+j+1)-u;
        saved = 0;
        for r=0:j-1
            temp = N(r+1)/(right(r+2)+left(j-r+1));
            N(r+1) = saved + right(r+2)*temp;
            saved = left(j-r+1)*temp;
        end
        N(j+1) = saved;
    end
    end
    
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Nurbs Book algorithm A2.3
    % create array of first k derivatives of all basis functions
    % there might be some mistakes in this algorithm
    function ders = DersBasisFuns(i,u,p,n,U)
    % Input:    i current knotspan
    %           u evaluation point
    %           p degree
    %           n degree of derivative (n<=p)
    %           U knotVector
    % Output    ders array, ders(k+1,j+1) stores the k-th derivative of N_i-p+j,p
    %           where 1 <=k <=n+1 and 1 <= j <= p+1
    ders = zeros(n+1,p+1);
    ndu = zeros(p+1,p+1); 
    ndu(1,1) = 1;
    left = zeros(1,p+1);
    right = zeros(1,p+1);
    for j=1:p
        left(j+1) = u- U(i+1-j+1);
        right(j+1) = U(i+j+1)-u;
        saved = 0;
        for r = 0 : (j-1) 
            % index shift wrt NURBS Book
            ndu(j+1,r+1) = right(r+2) + left(j-r+1);
            temp = ndu(r+1,j)/ndu(j+1,r+1);

            ndu(r+1,j+1) = saved + right(r+2)*temp;
            saved =left(j-r+1)*temp;
        end
        ndu(j+1,j+1) = saved;
    end
    for j= 0 : p % load basis functions
        ders(1,j+1) = ndu(j+1,p+1);
    end 
    % compute derivatives
    for r = 0 : p
        s1 = 0;
        s2 = 1;
        a = zeros(p+n,p+n);
        a(1,1) = 1;
        % compute kth derivative
        for k = 1 : n
            d = 0;
            rk = r - k;
            pk = p - k;
            if (r >=k )
                a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1); % Division by 0 happens! Index mistake? The plots seem reasonable.
                d =  a(s2+1,1)*ndu(rk+1,pk+1);
            end
            if rk >= -1
                j1 = 1;
            else
                j1 = -rk;
            end
            if r-1 <= pk
                j2 = k-1;
            else
                j2 = p-1;
            end
            for j = j1:j2
                a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1); % maybe rk+j+2? No. What else could be wrong?
                d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
            end
            if r <= pk
                a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                d = d+ a(s2+1,k+1)*ndu(r+1,pk+1);
            end
            ders(k+1,r+1) = d;
            j = s1;
            s1 = s2;
            s2 = j;
        end    
    end
    r= p;
    for k = 1:n
        for j = 0:p
            ders(k+1,j+1) =ders(k+1,j+1)* r;
        end
        r = r*(p-k );
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.4
    % evaluate one basis function
    function Nip = OneBasisFun(p,m,U,i,u)
    % compute the basis function Nip
    % Input     p,m,U,i,u
    % Output    Nip
    if ( (i == 0 && u == U(1)) || (i == m-p && u == U(m+1)))
        Nip = 1;
        return
    end
    if ( u < U(i+1) || u >= U(i+p+2))
        Nip = 0;
        return
    end
    N = zeros(1,p+1);
    % initialise basis functions
    for j= 0:p
        if ( u >= U(i+j+1) && u < U(i+j+2))
            N(j+1) = 1;
        else
            N(j+1) = 0;
        end
    end

    for k=1:p
        if( N(1) == 0 )
            saved = 0;
        else
            saved = ((u-U(i+1)) * N(1))/(U(i+k+1) - U(i+1));
        end
        for j = 0:(p-k)
            Uleft = U(i+j+2);
            Uright = U(i+j+k+2);
            if(N(j+2) == 0)
                N(j+1) = saved;
                saved =0;
            else
                temp = N(j+2)/(Uright - Uleft);
                N(j+1) = saved + (Uright -u)*temp;
                saved = (u - Uleft)*temp;
            end
        end
    end
    Nip = N(1);

    end
