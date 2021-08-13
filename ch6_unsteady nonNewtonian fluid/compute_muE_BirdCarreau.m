function muE = compute_muE_BirdCarreau(n, mu0, muinf, nt, JXYe, velocityE)
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

%%%%%%% 各个结点对应局部坐标
xi = [-1 0 1 -1 0 1 -1 0 1];
et = [-1 -1 -1 0 0 0 1 1 1 ];
%%%%%%% 各个结点对应局部坐标

%%%%%%% 初始化VIE
muE = zeros(9,1);
%%%%%%% 初始化VIE

%%%%%%% 循环计算单元内各个结点黏度
for i=1:9
    
    %%%%%%%% 速度插值函数对xi和et的导数
    [~, fyV_xi, fyV_et] = shapeFunction(xi(i), et(i) , 'quadratic');
    %%%%%%%% 速度插值函数对xi和et的导数
    
    %%%%%%%% Jacobi相关计算
    Jacobi = [fyV_xi, fyV_et]' * JXYe;
    fy_x = Jacobi \ [fyV_xi, fyV_et]';
    %%%%%%%% Jacobi相关计算
    
    %%%%%%%%%  单元内结点黏度计算
    %%%%%%%%%  单元内结点黏度计算
    vxxvxx = (fy_x(1,:) * velocityE(:,1))^2;
    vyyvyy = (fy_x(2,:) * velocityE(:,2))^2;
    vxyvyx = (fy_x(2,:) * velocityE(:,1) + fy_x(1,:) * velocityE(:,2))^2;
    SR = 2*vxxvxx + 2*vyyvyy + vxyvyx;  %剪切速率
    muE(i,1) = muinf + (mu0-muinf) * (1+(nt*SR)^2)^((n-1)/2);% BC本构方程 %公式(6.1)
    %%%%%%%%%  单元内结点黏度计算
end
%%%%%%% 循环计算单元内各个结点黏度