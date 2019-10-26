function [cellmap] = cilia_indentify(mo,Nframes)
%UNTITLED from BF video identify the position of cells using std
%thresholding and frequency filter???
%   Detailed explanation goes here
fs= mo.read([2,Nframes])
s=std(double(fs(:,:,1:Nframes)),[],3);

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

