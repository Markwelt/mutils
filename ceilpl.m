function limyt = ceilpl(valu, scaling)
% CEILPL returns the ceiled value(s) within its magnitute, that is a meaningful number above
% the value useful for plotting and setting ylim 
% optional second argument scaling (default = 5) adjust the rounding

if nargin == 1 
    scaling = 5;
end
if any(valu==0) || any(isinf(valu)), error('All inputs must be finite non-zero values' ), end

magn = floor(log10(valu));

imagnt = 10.^(-magn-1);

limyt = ceil(valu.*imagnt*scaling)./imagnt/scaling;

end
