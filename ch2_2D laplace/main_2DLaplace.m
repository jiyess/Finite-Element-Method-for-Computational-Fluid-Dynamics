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

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 生成网格数据 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  区域几何尺寸及网格划分参数
H = 4;    % 区域总高
L = 10;   % 区域总长
Nx = 24;   % 水平方向的网格数量
Ny = 12;   % 竖直方向的网格数量，选择能被2整除的数
%%%%%%%  区域几何尺寸及网格划分参数
[JXY,JM,noElement,noNodes,BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_tri(H, L, Nx, Ny);

plotTriGrid(JM, JXY, noElement, noNodes, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 生成网格数据 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 设定边界条件 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JB1 = [BP2, zeros(size(BP2))];  % 构建第一类边界条件数据JB1

dfy_dn1 = 0;  %构建第二类边界条件数据JB2
dfy_dn3 = 0;
dfy_dn4 = 2;
JB21 = [BE1, ones(size(BE1(:,1),1),2)*dfy_dn1];
JB23 = [BE3, ones(size(BE3(:,1),1),2)*dfy_dn3];
JB24 = [BE4, ones(size(BE4(:,1),1),2)*dfy_dn4];
JB2 = [JB21; JB23; JB24];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 设定边界条件 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 计算单元刚度阵并组装全局刚度阵 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffness = zeros(noNodes);  %初始化K
force = zeros(noNodes, 1);  %初始化b

for i = 1:noElement
    %提取单元三个结点的坐标,存入Jx和Jy
    Jx = JXY(JM(i,:),1);
    Jy = JXY(JM(i,:),2);
    
    % 计算三角形单元的面积
    Area = ((Jx(2)-Jx(1))*(Jy(3)-Jy(1))-(Jx(3)-Jx(1))*(Jy(2)-Jy(1)))/2.0;
    Area2 = 2 * Area;
    
    % 计算单元中a1,a2,a3,b1,b2,b3,c1,c2,c3
    %    a(1,1) = (Jx(2)*Jy(3)-Jx(3)*Jy(2)) / Area2;
    %    a(2,1) = (Jx(3)*Jy(1)-Jx(1)*Jy(3)) / Area2;
    %    a(3,1) = (Jx(1)*Jy(2)-Jx(2)*Jy(1)) / Area2;
    
    b(1,1) = (Jy(2)-Jy(3)) / Area2;
    b(2,1) = (Jy(3)-Jy(1)) / Area2;
    b(3,1) = (Jy(1)-Jy(2)) / Area2;
    
    c(1,1) = (Jx(3)-Jx(2)) / Area2;
    c(2,1) = (Jx(1)-Jx(3)) / Area2;
    c(3,1) = (Jx(2)-Jx(1)) / Area2;
    
    % 计算单元刚度阵
    stiffnessElement = Area * (b * b' + c * c');
    
    % 组装总体刚度阵
    stiffness(JM(i,:),JM(i,:)) = stiffness(JM(i,:),JM(i,:)) + stiffnessElement;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 计算单元刚度阵并组装全局刚度阵 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 计算单元载荷向量并组装全局载荷向量 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(JB2(:,1))
    
    II=JB2(i,1);   % 提取单元序号
    bondary_side_number=JB2(i,2);  % 提取边界变序号
    
    Jx = JXY(JM(II,:),1);
    Jy = JXY(JM(II,:),2);
    
    f1 = (2*JB2(i,5) + JB2(i,6)) / 6;
    f2 = (JB2(i,5) + 2*JB2(i,6)) / 6;
    
    switch bondary_side_number
        case 1 % 单元边号为1时Fe计算
            L = sqrt((Jx(1)-Jx(2))^2 + (Jy(1)-Jy(2))^2);
            forceElement = L * [f1;f2;0];
        case 2 % 单元边号为2时Fe计算
            L = sqrt((Jx(2)-Jx(3))^2 + (Jy(2)-Jy(3))^2);
            forceElement = L * [0;f1;f2];
        case 3 % 单元边号为3时Fe计算
            L = sqrt((Jx(3)-Jx(1))^2 + (Jy(3)-Jy(1))^2);
            forceElement = L * [f2;0;f1];
    end
    
    % 组装全局载荷向量
    force(JM(II,:),1) = force(JM(II,:),1) + forceElement;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 计算单元载荷向量并组装全局载荷向量 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prescribedDof = JB1(:,1);
% free Dof : activeDof
activeDof = setdiff( (1:noNodes)',(prescribedDof) );
force = force - stiffness(:,prescribedDof) * JB1(:,2);

% solution 求解方程 
displacements = zeros(noNodes,1);
displacements(activeDof) = stiffness(activeDof,activeDof) \ force(activeDof);
displacements(prescribedDof) = JB1(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = [JXY, displacements]

figure
trisurf(JM, JXY(:,1), JXY(:,2), displacements)
colorbar
view(0,90)
