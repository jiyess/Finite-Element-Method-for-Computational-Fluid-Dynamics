function forceElement = compute_forceElement(JXe, gaussPts, gaussWts)
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

%%% 初始化forceElement
forceElement = zeros(length(JXe),1);
%%% 初始化forceElement

for i = 1 : length(gaussWts)
%%%% 插值函数及其倒数
   fy=[1/2*gaussPts(i)*(gaussPts(i)-1); (1-gaussPts(i))*(1+gaussPts(i)); 1/2*gaussPts(i)*(gaussPts(i)+1)];   % 式（1.46）
   dfy_dxi=[gaussPts(i)-1/2; -2*gaussPts(i); gaussPts(i)+1/2 ]; % 式（1.47）
%%%% 插值函数及其倒数

%%%%  Jacobi相关
    dx_dxi = dfy_dxi'*JXe;   % 式（1.57）
    J = dx_dxi;    %   式（1.55）
%%%%  Jacobi相关

%%%%  计算forceElement
    forceElement = forceElement - gaussWts(i) * fy * det(J);  % 式（1.63）
%%%%  计算forceElement

end