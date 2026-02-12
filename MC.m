function [T_recovery, iter] = MC(omega, alpha, beta, T, trIndex, tol, maxiter, a, b)
%% Bounded kernel norm regularization algorithm, using alternating direction multiplier method to minimize kernel norm and obtain prediction score
% Usage:  [T_recovery, iter] = MC(omega, alpha, beta, T, trIndex, tol1, tol2, maxiter, a, b)
%
% Inputs:
%        omega, alpha, beta        - parameters needed to give.
%        T                  - the target matrix with only known entries and the unobserved entries are 0.
%        trIndex            - a matrix recording the observed positions in the target matrix.
%        tol                - tolerance of termination conditions.
%        maxiter            - maximum number of iterations.
%        a, b               - the left and right endpoints of the bounded interval.
%¦Ø£¬¦Á£¬¦Â-The parameters that need to be provided.
%T-Only the target matrix of known entries and unobserved entries are 0.
%trIndex¡ª¡ªA matrix that records the observed positions in the target matrix.
%Tolerance of termination conditions.
%maxiter¡ªMaximum number of iterations.
%a£¬b-The left and right endpoints of bounded intervals.
% Outputs:
%        T_recovery         - the completed matrix.
%        iter               - the number of iterations.
%

X =T;
Y = X;
Z = X;

i = 1;
stop = 1;

while(stop > tol)  

    %the process of computing Y
    tran = (1/beta) * (Z + alpha * (T.* trIndex)) + X;
    Y = tran - (alpha/ (alpha + beta)) * (tran.* trIndex);
    Y(Y < a) = a;
    Y(Y > b) = b;
    
    %the process of computing X
    X_1 = svt(Y- 1/beta* Z, omega/beta);

    %the process of computing Z

    Z = Z + 1*beta * (X_1 - Y);


    stop= norm(X_1 - X, 'fro')/ norm(X, 'fro');
    X = X_1;
    i = i+1;
    
    if i < maxiter
        iter = i - 1;

    else
        iter = maxiter;
%         warning('reach maximum iteration~~do not converge!!!');
%         stop
        stop
        break
    end
    
end

T_recovery = Y;

end