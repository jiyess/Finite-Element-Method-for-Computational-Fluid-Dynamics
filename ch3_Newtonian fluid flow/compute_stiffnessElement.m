function [Be1, Be2, De11, De12, De21, De22, Ce1, Ce2] = ...
    compute_stiffnessElement(JXYe, mu, gaussPts, gaussWts)
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

muNodesVec = mu;
%%%%%%% 初始化Be1和Be2
Be1 = zeros(4,9);
Be2 = zeros(4,9);
%%%%%%% 初始化Be1和Be2

%%%%%%% 初始化De11、De12、De21和De22
De11=zeros(9,9);
De12=zeros(9,9);
De21=zeros(9,9);
De22=zeros(9,9);
%%%%%%% 初始化De11、De12、De21和De22

%%%%%%% 初始化Ce1和Ce2
Ce1=zeros(9,4);
Ce2=zeros(9,4);
%%%%%%% 初始化Ce1和Ce2

xi = gaussPts;
et = gaussPts;

for i = 1 : length(gaussWts)
    for j = 1 : length(gaussWts)
        
        %%%%%%%% 压力插值函数
        fyP = shapeFunction(xi(i), et(j) , 'linear');
        %%%%%%%% 压力插值函数
        
        %%%%%%%% 速度插值函数对xi和et的导数
        [fyV, fyV_xi, fyV_et] = shapeFunction(xi(i), et(j) , 'quadratic');
        %%%%%%%% 速度插值函数对xi和et的导数
        
        %%%%%%%% Jacobi相关计算
        Jacobi = [fyV_xi, fyV_et]' * JXYe;
        fy_x = Jacobi \ [fyV_xi, fyV_et]';
        
        det_Jacobi = det(Jacobi);
        %%%%%%%% Jacobi相关计算
        
        commonFac = gaussWts(i) * gaussWts(j) * det_Jacobi;
        
        %%%%%%%% Be1和Be2单元方程子块计算
        Be1 = Be1 + commonFac * fyP * fy_x(1,:);
        Be2 = Be2 + commonFac * fyP * fy_x(2,:);
        %%%%%%%% Be1和Be2单元方程子块计算
        
        %%%%%%% De11、De12、De21和De22单元方程子块计算
        if ~isscalar(muNodesVec) % 单元内黏度插值，非牛顿流体
            mu = fyV' * muNodesVec;
        end
        
        De11 = De11 + mu * commonFac * (2*fy_x(1,:)'*fy_x(1,:) + fy_x(2,:)'*fy_x(2,:));
        De12 = De12 + mu * commonFac * fy_x(1,:)'*fy_x(2,:);
        De21 = De21 + mu * commonFac * fy_x(2,:)'*fy_x(1,:);
        De22 = De22 + mu * commonFac * (2*fy_x(2,:)'*fy_x(2,:) + fy_x(1,:)'*fy_x(1,:));
        %%%%%%% De11、De12、De21和De22单元方程子块计算
        
        %%%%%%%% Ce1和Ce2单元方程子块计算
        Ce1 = Ce1 + commonFac * (fyP * fy_x(1,:))';
        Ce2 = Ce2 + commonFac * (fyP * fy_x(2,:))';
        %%%%%%%% Ce1和Ce2单元方程子块计算
    end
end