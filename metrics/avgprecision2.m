function [ avp ] = avgprecision2( pScore, nScore, atK )
avp = avgprecisionAtK( pScore, nScore );

if nargin > 2
    avp = repmat(avp,1,length(atK));
end