% val = clamp(val,minval,maxval)
% ---------------------------------------
%
% Clamp an nd-array between minval and maxval
%
% Inputs:
%  val        	Input matrix
%  minval		Any value of val smaller than minval will be made equal to minval
%  maxval       Any value of val larger than maxval will be made equal to maxval
%
% Outputs:
%  vak        	Clamped value
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

function val = clamp(val,minval,maxval)

if isempty(minval) % no lower bound
    map = val > maxval;
	val(map) = maxval;
end

if isempty(maxval) % no upper bound
	map = val < minval;
	val(map) = minval;
end

if ~isempty(minval) && ~isempty(maxval)
    map = val > maxval;
    val(map) = maxval;
    map = val < minval;
    val(map) = minval;
end