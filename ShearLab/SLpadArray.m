function paddedArray = SLpadArray( array,newSize)
%SLPADMATRIX Summary of this function goes here
%   Detailed explanation goes here

    padSizes = zeros(1,length(newSize));
    S.type = '()';
    S.subs = cell(1,length(newSize));

    for k = 1:length(newSize)
        currSize = size(array,k);
        sizeDiff = newSize(k)-currSize;
        if mod(sizeDiff,2) == 0
            padSizes(k) = sizeDiff/2;    
            S.subs{k} = ':';
        else
            padSizes(k) = ceil(sizeDiff/2);
            if mod(currSize,2) == 0
                S.subs{k} = 2:(newSize(k)+1);
            else
                S.subs{k} = 1:(newSize(k));
            end
        end    
    end

    paddedArray = padarray(array,padSizes);
    paddedArray = subsref(paddedArray,S);


end

%
%  Copyright (c) 2014. Rafael Reisenhofer
%
%  Part of ShearLab3D v1.1
%  Built Mon, 10/11/2014
%  This is Copyrighted Material
%
%  If you use or mention this code in a publication please cite the website www.shearlab.org and the following paper:
%  G. Kutyniok, W.-Q. Lim, R. Reisenhofer
%  ShearLab 3D: Faithful Digital SHearlet Transforms Based on Compactly Supported Shearlets.
%  ACM Trans. Math. Software 42 (2016), Article No.: 5.
