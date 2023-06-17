% r = getRadAvg(im)
% ---------------------------------------
%
% Return the radial average of input square image im
%
% Inputs:
%  im        	Square image
%
% Outputs:
%  r        	Radial average
%
% ---------------------------------------
%
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch, 
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function r = getRadAvg(im)

if size(im,1) ~= size(im,2)
    im = imresize(im,[max(size(im)) max(size(im))],'bilinear');
elseif length(size(im)) ~= 2
    error('getRadAvg supports only 2D matrix as input');
end

r = mean(im2pol(im),1);
r = imresize(r,[1 ceil(size(im,2)/2)],'bilinear');