classdef bSplBasFun < bSplineBasis
     properties
        index = 1;
     end
     methods (Access = public)
        function obj = bSplBasFun(index,basis)
            % standard constructor for subclass of bSplinBasis 
            % Nooo, overloading not possible!
             if nargin >0
                obj.index = index;
                obj.a = basis.a;
                obj.b = basis.b;
                obj.knotspan = basis.knotspan;
                obj.p = basis.p;
                obj.resol = basis.resol;
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
        % other methods   
        function C = generOneBasisFun(obj)
                C = zeros(obj.sP,1);
            for j = 1 : obj.sP
                C(j,1) = OneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j));

            end
        end
        function C = generDersOneBasisFun(obj)
                C = zeros(obj.sP,1);
            for j = 1 : obj.sP
                temp = DersOneBasisFun(obj.p,obj.m,obj.knotVector,obj.index,obj.plotVector(j),1)
                C(j,1) = temp(:,2);
            end
        end
        function plotOneBasisFun(obj,C)
            plot(obj.plotVector,C(:,1),'y-');
        end
    

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A2.5
    % evaluate first k derivatives of one basis function
    % maybe there is still some index error
    function ders = DersOneBasisFun(p,m,U,i,u,n)
    % compute derivatives of basis function Nip
    % Input     p,m,U,i,u,n
    % Output    ders
    ders = zeros(1,n+1);
    if ( u< U(i+1) || u >= U(i+p+2) )
        for k= 0:n
            ders(k+1) = 0;
            return;
        end
    end
    N = zeros(p+2,n+1); % changed to p+2, maybe p+1 suffices
    for j=0:p % initialize 0-degree functions
        if (u >= U(i+j+1) && u < U(i+j+2))
            N(j+1,1) = 1;
        else
            N(j+1,1) = 0;
        end
    end
    for k = 1:p % compute full triangular table?
        if (N(1,k) == 0)
            saved = 0;
        else
            saved = ((u-U(i+1))*N(1,k))/(U(i+k+1)-U(i+1));
        end
        for j = 0:(p-k)
            Uleft = U(i+j+2);
            Uright = U(i+j+k+2);
            if( N(j+2,k) == 0)
                N(j+1,k+1) = saved;
                saved = 0;
            else
                temp = N(j+2,k)/(Uright - Uleft);
                N(j+1,k+1) = saved + (Uright - u) * temp;
                saved = (u-Uleft)*temp;
            end
        end
    end
    ders(1) = N(1,p+1);
    ND = zeros(1,n+1);
    for k=1:n
        for j=0:k
            ND(j+1) = N(j+1,p-k+1);
        end
        for jj = 1:k
            if (ND(1) == 0)
                saved = 0;
            else
                saved = ND(1)/(U(i+p-k+jj+1) - U(i+1));
            end
            for j=0:(k-jj) % index as in Nurbs Book
                Uleft = U(i+j+2);
                Uright = U(i+j+p+jj+2);
                if( ND(j+2) == 0)
                    ND(j+1) = (p-k+jj)*saved;
                else
                    temp = ND(j+2)/(Uright -Uleft);
                    ND(j+1) = (p-k+jj)*(saved-temp);
                    saved = temp;
                end
            end
        end
        ders(k+1) = ND(1);
    end
    end
