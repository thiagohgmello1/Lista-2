function [ratio] = SpecRatio(A,met,e,imax)
%SpecRatio determine Spectral Ratio of matrix A

if met == "PowerMet"
    ratio = PowerMet(A,e,imax);
end

end

