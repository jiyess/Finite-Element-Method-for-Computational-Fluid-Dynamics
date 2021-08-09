function [fy, dfydxi, dfydet] = shapeFunction(xi, et , type)
%    see Bi Chao, "Finite Element Method for Computational Fluid
%    Dynamics and Its Detailed Programming"
%    中文：毕超《计算流体力学有限元方法及其编程详解》
%    
%    Copyright (C) 2021 Ye Ji, Dalian University of Technology
%    
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%    
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

switch type
    
    case 'linear'
        
        fy = 1/4 * [(1-xi)*(1-et);
            (1+xi)*(1-et);
            (1+xi)*(1+et);
            (1-xi)*(1+et)];
        
        if nargout == 3
            dfydxi = [];
            dfydet = [];
        end
        
    case 'quadratic'
        fy = [1/4*xi*et*(xi-1)*(et-1);
            1/2*et*(1-xi^2) *(et-1);
            1/4*xi*et*(xi+1)*(et-1);
            1/2*xi*( xi-1)*(1-et^2);
            (1-xi^2)*(1-et^2);
            1/2*xi*( xi+1)*(1-et^2);
            1/4*xi*et*(xi-1)*(et+1);
            1/2*et*(1-xi^2) *(et+1);
            1/4*xi*et*(xi+1)*(et+1)];
        
        if nargout == 3
            dfydxi = [1/4*et*(xi-1)*(et-1) + 1/4*xi*et*(et-1);
                -et*xi*(et-1);
                1/4*et*(xi+1)*(et-1)+1/4*xi*et*(et-1);
                1/2*(xi-1)*(1-et^2)+1/2*xi*(1-et^2);
                -2*xi*(1-et^2);
                1/2*(xi+1)*(1-et^2)+1/2*xi*(1-et^2);
                1/4*et*(xi-1)*(et+1)+1/4*xi*et*(et+1);
                -et*xi*(et+1);
                1/4*et*(xi+1)*(et+1)+1/4*xi*et*(et+1)];
            
            dfydet = [1/4*xi*(xi-1)*(et-1)+1/4*xi*et*(xi-1);
                1/2*(1-xi^2)*(et-1)+1/2*et*(1-xi^2);
                1/4*xi*(xi+1)*(et-1)+1/4*xi*et*(xi+1);
                -xi*et*(xi-1);
                -2*(1-xi^2)*et;
                -xi*et*(xi+1);
                1/4*xi*(xi-1)*(et+1)+1/4*xi*et*(xi-1);
                1/2*(1-xi^2)*(et+1)+1/2*et*(1-xi^2);
                1/4*xi*(xi+1)*(et+1)+1/4*xi*et*(xi+1)];
        end
        
        
end


