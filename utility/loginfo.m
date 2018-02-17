function loginfo( level, format, varargin )
%LOGINFO Summary of this function goes here
%   Detailed explanation goes here
global G_LOG_LEVEL

if (level >= G_LOG_LEVEL)
    fprintf(format, varargin{:});

end

