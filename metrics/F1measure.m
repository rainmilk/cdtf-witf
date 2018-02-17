function [ f1 ] = F1measure( recall, precision )
%F1MEASURE Summary of this function goes here
%   Detailed explanation goes here

f1 = 2*(recall.*precision) ./ (recall + precision);
end

