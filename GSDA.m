function xk = GSDA(rd, x0)
% initializations
xk = x0;
gam=1e-5;
rho = 0.5;
Q0=1;
iter = 0;

rk = [cos(xk(1)) + cos(xk(1) + xk(2)) + cos(xk(1) + xk(2) + xk(3)); sin(xk(1)) + sin(xk(1) + xk(2)) + sin(xk(1) + xk(2) + xk(3))] - rd;
fkz = 1/2*norm(rk,2)^2;
C0 = fkz;
% Initial Jacobian matrix
Jk = [-sin(xk(1))-sin(xk(1)+xk(2))-sin(xk(1)+xk(2)+xk(3)) -sin(xk(1)+xk(2))-sin(xk(1)+xk(2)+xk(3)) -sin(xk(1)+xk(2)+xk(3)); cos(xk(1))+cos(xk(1)+xk(2))+cos(xk(1)+xk(2)+xk(3)) cos(xk(1)+xk(2))+cos(xk(1)+xk(2)+xk(3)) cos(xk(1)+xk(2)+xk(3))];
% Gradient function gk = J'F
gk = Jk'*rk;
% Initiate b
diag = ones(length(x0), 1);
%Initial Direction
dk = -diag.*gk;
normg = norm(gk, inf);
while (normg > 10^(-5)&&iter<=1000)          %( norm(fk) > 10^(-5)*sqrtn + 10^(-4)*norm(f0) )
    tk = 1;
    xkn = xk + tk*dk;
    rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
    fkzn = 0.5*norm(rkn,2)^2;
    dd = dk'*gk;
    % This while loop will find a suitable step (armijo with backtracking)
            while(fkzn > C0 + gam*tk*(dd))
%             while (fkzn > fkz + sigma*tk*gk'*dk)
             tk = rho*tk;
             xkn = xk + tk*dk;
             rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
%             numf=numf+1;
             fkzn = 0.5*(rkn'*rkn); 
            end
    xkn = xk + tk*dk;
    Jkn = [-sin(xkn(1))-sin(xkn(1)+xkn(2))-sin(xkn(1)+xkn(2)+xkn(3)) -sin(xkn(1)+xkn(2))-sin(xkn(1)+xkn(2)+xkn(3)) -sin(xkn(1)+xkn(2)+xkn(3)); cos(xkn(1))+cos(xkn(1)+xkn(2))+cos(xkn(1)+xkn(2)+xkn(3)) cos(xkn(1)+xkn(2))+cos(xkn(1)+xkn(2)+xkn(3)) cos(xkn(1)+xkn(2)+xkn(3))];
    rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
%     xkn = xk + alph*dir;
    sk = xkn - xk;
    gkn = Jkn'*rkn;
    c1 = Jkn'*Jkn*sk;
    %   gradient  d1=Jk^Trk
    d1  = Jkn'*rkn;  
%   gradient  d2=J_{k-1}rk
    d2 = Jk'*rkn;
    %d=gt-d2; % d2 = ybar = Jk^Trk - J_{k-1}^Trk
    d3  = d1 - d2;
 % the weight, W = B 
    wgt = diag;
    % wgt = weight;
    wgt2 = wgt.^2;
%    sum1 = sum((sk.^2).*(wgt2));
    sum2 = sum((sk.^4).*(wgt2));

    Beta1 = sk'*c1+ sk'*d3; % short beta
    % the entries of the correction matrix is computed using this formulation
    cdiag2 = (( (Beta1 + (sk'*(wgt2.*sk))- (sk'*(diag.*sk)) )/sum2 )*(sk.^2) - 1).*(wgt2); % newly written formulation (21/05/2021)
    tol1 = 1e-3;
    diagk = diag+cdiag2;
    if diagk >= tol1
    dk = -gkn./diagk;
    else
    dk = -gkn;
    end
     % updating the values of xc, gc, f's, Q & C
    diag = diagk;
    fkz = fkzn;
    etak=0; %%Armijo linesearch
    Q1=etak*Q0+1;
    C1=(etak*Q0*C0+fkz)/Q1;
    C0=C1;
    xk = xkn;
    gk = gkn;
    dk = dk;
    iter = iter +1;  
    normg=norm(gk,inf);

end