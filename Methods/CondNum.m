function [cond] = CondNum(A,e,imax)
%CondNum defines the condition number of A matrix

if ~exist("e","var")
    e = 1e-5;
end
if ~exist("imax","var")
    imax = 50;
end

if A - Transp(A) == 0
    cond = abs(PowerMet(A,e,imax))/abs(InvPowerMet(A,e,imax));
else
    B = MatrixMulti(Transp(A),A);
    smax = sqrt(abs(PowerMet(B,e,imax)));
    smin = sqrt(abs(InvPowerMet(B,e,imax)));
    cond = smax/smin;
end
end

