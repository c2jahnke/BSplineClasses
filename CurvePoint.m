    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nurbs Book algorithm A3.1 adjusted
    % evaluate a point on a B-Spline curve

    function C = CurvePoint(n,p,U,P,u)
    % Compute a curve point
    % Input:    n number of basis functions -1, (n = m - p -1)
    %           p spline degree
    %           U knotVector
    %           P vector of control points
    % Output:   C curve point

    span = FindSpan(n,p,u,U);
    N = BasisFun(span,u,p,U); %BasisFun(i,u,p,U) 
    if( span+1 > length(P) )
        C = P(end,:);
    else
        C = N * P( span+1-p : span+1 , : );
    end
    %C = 0;
    %for i=0:p
    %    C = C + N(i+1)*P(span - p+i+1);
    %end
    end
    