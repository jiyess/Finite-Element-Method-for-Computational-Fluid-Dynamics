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

clc;clear;close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        计算参数      %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
penaltyCoef = 1e12;  % 罚系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        计算参数      %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       读取网格数据     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  区域几何尺寸及网格划分参数
H = 0.04;   % 区域总高
L = 0.1;   % 区域总长
Nx = 8;     % 水平方向的网格数量
Ny = 8;     % 竖直方向的网格数量
%%%%%%%  区域几何尺寸及网格划分参数
[JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_rectangle(H, L, Nx, Ny);
%%%%%%% 调用四边形网格绘制程序
plotRectangularGrid(JMP,JXYV,noElement, noNodesV, 1);
%%%%%%% 调用四边形网格绘制程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       读取网格数据     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        设定物料黏度     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1000;      % 物料黏度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        设定物料黏度     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    设定边界条件   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity B.C.
u1 = 0;     v1 = 0;
u3 = 0.01;  v3 = 0;
JBV1 = [BP1(:), [u1,v1] .* ones(length(BP1),1)];
JBV3 = [BP3(:), [u3,v3] .* ones(length(BP3),1)];
JBV=[JBV1; JBV3];

% pressure B.C.
P2=0;
P4=1000;
JBP2=[BE2, [P2,P2] .* ones(size(BE2,1),1)];
JBP4=[BE4, [P4,P4] .* ones(size(BE2,1),1)];
JBP=[JBP2; JBP4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    设定边界条件   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       初始化总体方程各个子块    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DP11=zeros(noNodesV);
DP12=zeros(noNodesV);
DP21=zeros(noNodesV);
DP22=zeros(noNodesV);
F1=zeros(noNodesV, 1);
F2=zeros(noNodesV, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       初始化总体方程各个子块    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   计算单元方程系数矩阵子块并组装  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
for i=1:noElement
    
    JXYe = JXYV(JMV(i,:),:);
    
    [DPe11, DPe12, DPe21, DPe22] =...
        compute_DPe(JXYe, mu, penaltyCoef, gaussPts, gaussWts);
    
    DP11(JMV(i,:),JMV(i,:)) = DP11(JMV(i,:),JMV(i,:)) + DPe11;
    DP12(JMV(i,:),JMV(i,:)) = DP12(JMV(i,:),JMV(i,:)) + DPe12;
    DP21(JMV(i,:),JMV(i,:)) = DP21(JMV(i,:),JMV(i,:)) + DPe21;
    DP22(JMV(i,:),JMV(i,:)) = DP22(JMV(i,:),JMV(i,:)) + DPe22;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   计算单元方程系数矩阵子块并组装  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   计算单元方程右边向量子块并组装  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(JBP(:,1))
    
    JXYe = JXYV(JMV(JBP(i,1),:),:);
    
    [Fe1, Fe2] = compute_forceElement(JXYe,JBP(i,:), gaussPts, gaussWts);
    
    F1(JMV(JBP(i,1),:),1) = F1(JMV(JBP(i,1),:),1) + Fe1;
    F2(JMV(JBP(i,1),:),1) = F2(JMV(JBP(i,1),:),1) + Fe2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   计算单元方程右边向量子块并组装  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       组合总体方程      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffness = [DP11, DP12;
            DP21, DP22];
force = [-F1;-F2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       组合总体方程      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noDofs = 2*noNodesV;
prescribedDof = cat(1,JBV(:,1),JBV(:,1)+noNodesV);
prescribedValue = cat(1,JBV(:,2),JBV(:,3));
% free Dof : activeDof
activeDof = setdiff( (1:noDofs)', prescribedDof );
force = force - stiffness(:,prescribedDof) * prescribedValue;

% solution 求解方程 
solution = zeros(noDofs,1);
solution(activeDof) = stiffness(activeDof,activeDof) \ force(activeDof);
solution(prescribedDof) = prescribedValue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velocityX_k1 = solution(1 : noNodesV);
velocityY_k1 = solution(1+noNodesV : 2*noNodesV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 逐单元计算内结点压力，并进行累加求平均 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgPressure = zeros(noNodesV,2);
for e = 1:noElement
    
    JXYe = JXYV(JMV(e,:),:);
    velocityE(:,1) = velocityX_k1(JMV(e,:),:);
    velocityE(:,2) = velocityY_k1(JMV(e,:),:);
    
    PIE = compute_PIE(JXYe, velocityE, penaltyCoef);
    
    avgPressure(JMV(e,:),1) = avgPressure(JMV(e,:),1) + PIE;
    avgPressure(JMV(e,:),2) = avgPressure(JMV(e,:),2) + 1;

end
pressure_k1 = avgPressure(:,1) ./ avgPressure(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 逐单元计算内结点压力，并进行累加求平均 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 输出结果，绘图 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velocity_k1 = sqrt(velocityX_k1.^2 + velocityY_k1.^2);
data = [JXYV, velocityX_k1, velocityY_k1, velocity_k1, pressure_k1]

figure
hold on
for i = 1:noElement
    patch(JXYV(JMP(i,:),1),JXYV(JMP(i,:),2),velocity_k1(JMP(i,:)), 'edgecolor','none')
end
colorbar
title('velocity\_k1')

figure
hold on
for i = 1:noElement
    patch(JXYV(JMP(i,:),1),JXYV(JMP(i,:),2),pressure_k1(JMP(i,:)), 'edgecolor','none')
end
colorbar
title('pressure\_k1')

% 存储结果
save('result_k1.mat', 'JMP', 'JXYV', 'velocityX_k1', 'velocityY_k1',...
    'velocity_k1', 'pressure_k1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 输出结果，绘图 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      出口速度提取     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UB = velocityX_k1(BP2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      出口速度提取     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%