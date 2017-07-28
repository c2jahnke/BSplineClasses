
a = 0;
b = 100;
N = 200; % number of elements
stepS = (b-a)/N;
basis = bSplBas(a,b,7,stepS,0.001);
% thanks to Alex!

f = @(x) sin(x);
Stiffn = zeros(basis.n);
elStiff = zeros(basis.p +1);
% Ax = b
% Lapl(u) = f
ngp = basis.p +3;
rhs = zeros(basis.n,1);
elRhs = zeros(basis.p+1,1);
for k = 1 : N % loop over elements
    [s,w]=lgwt(ngp,a+stepS*(k-1),a+stepS*k ); % ugly, add allKnots to bSplBas
    bVal = zeros(ngp, basis.p+1); % basis evaluation
    gradVal = zeros(ngp,basis.p+1); % derivative evaluation
    for j = 1:length(s) %basis.evalDersBasis(s(j))% not yet fully tested! Fails for high values
        temp = basis.evalDersBasis(s(j));%DersBasisFuns(k,s(j),basis.p,basis.knotVector,1) % replace by method of class!
        bVal(j,:) = temp(1,:);
        gradVal(j,:) = temp(2,:); % 
    end
    
    for l = 1 : basis.p+1
        elRhs(l) = sum(w.*f(s).*bVal(:,l));
        elStiff(l,l) = sum(w.*gradVal(:,l).^2);
        for kk = l+1 : basis.p +1
            elStiff(l,kk) = sum(w.*gradVal(:,l).*gradVal(:,kk));
            elStiff(kk,l) = elStiff(l,kk);
        end
          
            
    end
    rhs(k:k+ basis.p)= rhs(k:k+ basis.p) + elRhs;
    Stiffn(k:k+basis.p,k:k+basis.p) = Stiffn(k:k+basis.p,k:k+basis.p) + elStiff;
     
end
%% Use preconditioner

%% BC with Lagrange multipliers
A = zeros(basis.n +1);
b = zeros(basis.n +1,1);
A(1,2) = 1;
A(2,1) = 1;
A(2:end,2:end) = Stiffn;
b(1,1) = 0; % Dichlet BC
b(2:end) = rhs;

%u = pcg(A,b,0.1);
u = A\b; % use sum u * basis.functions as a curve,
y = u(2:end);
x = linspace(basis.a,basis.b,basis.n);
plot(x, y);


