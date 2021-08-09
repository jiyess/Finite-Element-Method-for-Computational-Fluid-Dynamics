function stiffnessElement = compute_stiffnessElement(JXe, gaussPts, gaussWts)
%    see Bi Chao, "Finite Element Method for Computational Fluid
%    Dynamics and Its Detailed Programming"
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

%%% 初始化 stiffnessElement
stiffnessElement = zeros(length(JXe));

for i = 1 : length(gaussWts)
    %%%% 插值函数对kesi的倒数
    dfy_dxi=[gaussPts(i)-1/2;  -2*gaussPts(i); gaussPts(i)+1/2 ];      % 式（1.47）
    %%%% 插值函数对kesi的倒数
    
    %%%%  Jacobi相关
    dx_dxi = dfy_dxi' * JXe;  % 式（1.57）
    J = dx_dxi;               %式（1.55）
    dfy_dx = dfy_dxi / J;   % 式（1.54）
    %%%%  Jacobi相关
    
    %%%%  计算stiffnessElement
    stiffnessElement = stiffnessElement + gaussWts(i) * (dfy_dx * dfy_dx') * det(J);   % 式（1.62）
    %%%%  计算stiffnessElement
    
end