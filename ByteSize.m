function [str, b, scale] = ByteSize(theVariable, scale) %#ok<INUSL>
% ByteSize returns the mem.usage of the provided variable(theVariable) to the given file
% identifier. 
% scale is assigned meaningfully according to the byte size if not stated

s = whos('theVariable');
b = s.bytes;
if nargin == 1 
    scale = floor(log(b)/log(1024));
    switch scale
        case 0
            scale = 'bytes';
        case 1
            scale = 'kb';
        case 2
            scale = 'mb';
        case 3
            scale = 'gb';
        case 4
            scale = 'tb';
        case -inf
            % Size occasionally returned as zero (eg some Java objects).
            scale = 'bytes';
            warning('Size occasionally returned as zero (eg some Java objects). Bytes assumed');
        otherwise
            scale = 'petabytes';
            warning('Over 1024 petabyte. petabytes assumed');
    end
end
switch scale
    case {'b','byte','bytes'}
        b = s.bytes;
    case {'kb','kbs','kilobyte','kilobytes'}
        b = b / 1024;
    case {'mb','mbs','megabyte','megabytes'}
        b = b / 1024^2;
    case {'gb','gbs','gigabyte','gigabytes'}
        b = b / 1024^3;
    case {'tb','tbs','terabyte','terabytes'}
        b = b / 1024^4;
    case {'pb','pbs','petabyte','petabytes'}
        b = b / 1024^5;
    otherwise
        scale = 'bytes';
end

str = sprintf('%.2f %s',b,scale);

end