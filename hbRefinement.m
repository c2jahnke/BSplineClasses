
clc, clear all, close all
a = 0;
b = 5;
p = 2;
N = 5;
resol = 0.01;
lvl = 2;
obj = hbSplBasML(a,b,p,N,resol,lvl)

cBas = obj.levelBas{1};
fBas = obj.levelBas{2};
refArea = [2 4];
[Char0, Char1, HB0, HB1, V0, V1, U, Ubar, Points, Qw] = OneDimHbRefinement(cBas,fBas, refArea);


    cBas.plotBasisStruct(HB0);
    hold on;
    fBas.plotBasisStruct(HB1);  
    
    
    function [Char0, Char1, HB0, HB1, V0, V1, U, Ubar, Points, Qw] = OneDimHbRefinement(cBas,fBas,refArea)
    % still not finished
    
    Omega0 = [cBas.a cBas.b];
    Omega1 = refArea;
    cIndOmega1_start = cBas.getIndexU(refArea(1));
    cIndOmega1_end = cBas.getIndexU(refArea(2));
    % private attribute is changed, use setter!
    cBas.activeKnots = [cBas.knotVector(cBas.p+1:cBas.p + cIndOmega1_start) cBas.knotVector(cBas.p + cIndOmega1_end: end -cBas.p)];
    cBas.activeIndex =[0:cBas.getIndexU(refArea(1)) cBas.getIndexU(refArea(2)):cBas.n-1];
    fIndOmega1_start = fBas.getIndexU(refArea(1));
    fIndOmega1_end = fBas.getIndexU(refArea(2));
    fBas.activeKnots = fBas.knotVector(fBas.p+fIndOmega1_start:fIndOmega1_end+fBas.p);
    fBas.activeIndex =[0:fBas.getIndexU(refArea(1)) fBas.getIndexU(refArea(2)):fBas.n-1];
    
    h0 = cBas.knotspan;
    h1 = fBas.knotspan;
    
    assert(~(Omega1(end) - Omega1(1) < cBas.p*h1),'Error: Omega1 to small for one basis function.');
    
    U = cBas.knotVector;
    n = cBas.n;
    
    %plotVector = a:cBas.resol:b;
    %sP = size(plotVector,2);
    Points = zeros(n,2); % not needed
    Points(:,1) = linspace(0,1, size(Points,1) );
    Points(:,2) = sin(2*pi*Points(:,1));
    
    X = Omega0(1)+h1:h0:Omega0(end)-h1; %% not really general, work this out later
    r= size(X,2)-1;
    [Ubar, Qw] = RefineKnotVectCurve(n,cBas.p,U,Points,X,r); % the curve points Qw are not needed

    V0 = cBas.generBasis();
    V1 = fBas.generBasis();

    % has to depend on degree
    % so far works well for all degrees 
    temp1 = U(:) < Omega1(1)+(cBas.p-1)*h0;% not general!!! Change to h1 or h0, test!
    temp2 = U(:) >= Omega1(end)-(cBas.p-1)*h0;
    if (Omega1(end) == Omega0(end))
        temp2 = zeros(size(U))';
    end
    % calculate characteristic matrices
    Char0 = zeros(size(U,2)-cBas.p-1);
    Char0 =Char0 + diag(temp1(cBas.p:end -2)) + diag(temp2(2:end-cBas.p));

    temp3 = (Ubar(:) < Omega1(1) + (cBas.p-1)*h1);
    temp4 = (Ubar(:) >= Omega1(end) -(cBas.p-1)*h1 );
    if (Omega1(end) == Omega0(end))
        temp4 = zeros(size(Ubar))';
    end
    Char1 = eye(fBas.n);
    Char1 = Char1 -diag(temp3(cBas.p+1:end-1)) - diag(temp4(1:end-cBas.p-1));
    
    % calculate hierarchical basis
    HB1 = V1*Char1;
    HB0 = V0*Char0;
    
    % truncation
    end
    
    function [THB0, trunq, q] = OneDimThb(Omega0,Omega1,h0,h1,V0,V1,U,Ubar,Char1,TOL,parameters)
    %V0 = generBasis(parameters.m,parameters.p,parameters.knotVector,parameters.knotVector,parameters.n,parameters.m);
    %V1 = generBasis(parameters.m,parameters.p,Ubar,Ubar,parameters.n,parameters.m);
    q = zeros(size(V1,2),size(V0,2));
    for k = 1 : size(V0,2) %% use Greville points or something similar so 
        %                   that this system is still uniquely solvable 
        q(:,k) = V1\V0(:,k);
    end
    % stabilize q
    % Konditionszahlabhängig
    q(q(:) < TOL) = 0;
    % charactaristic matrix Char0 Char1
    % coefficient matrix Points Qw
    % refinement matrix Q
    % For THB Splines
    p= parameters.p;
    trunq = q;
    % no good solution, so far works only for h1 = 0.5, h0 = 1
     if(h0 == 1 && h1 == 0.5 ) 
        trunq((Omega1(1)+p)/(h1)-(p-1):(Omega1(end))/(h1),:) = 0; % causes errors if h1 is very small;
      % trunq(2*(Omega1(1)+p)-(p-1):(Omega1(end))*2,:) = 0; % causes errors if h1 is very small;
     elseif(h0 == 2 && h1 == 1)
         trunq((Omega1(1)+p)/(h1)-(p-1)+p:(Omega1(end))/(h1),:)= 0;
%     generalize! write function that returns the support in knot indices    
     elseif(h0 == 4 && h1 == 2 ) % just works for p = 2
         trunq(floor((Omega1(1)+p+1)/(h1) + p):floor((Omega1(end)+1)/h1),:)= 0;
     else
         error('Truncation not yet implemented for h1,h0');
     end
    
    if(Omega1(end) == Omega0(end))
        trunq(Omega1(end)/h1:end,:) = 0;
    end
    
    if(Omega1(1) == Omega0(1))
        trunq(1:Omega1(end)/h1);
    end

    THB0 = V1*(trunq);
    THB0(THB0(:) < TOL) = 0;
    end
