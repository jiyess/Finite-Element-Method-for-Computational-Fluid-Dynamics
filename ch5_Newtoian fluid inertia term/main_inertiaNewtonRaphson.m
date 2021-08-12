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
% 速度项提出法 + Newton-Raphson迭代法

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤A开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     物性参数       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  区域几何尺寸及网格划分参数
H = 0.04;  %区域总高
L = 0.1;  %区域总长
Nx = 12;  %水平方向的网格数量，选择能被5整除的数
Ny = 12;   %竖直方向的网格数量
%%%%%%%  区域几何尺寸及网格划分参数

mu = 10;  % 黏度
rho = 1100;  % 密度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     读取网格数据  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP,...
    BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_inlet_outlet(H,L,Nx,Ny);

%%%%%%% 调用四边形网格绘制程序
plotRectangularGrid(JMP,JXYV,noElement, noNodesV, 1);
%%%%%%% 调用四边形网格绘制程序

noDofs = 2*noNodesV + noNodesP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    设定边界条件   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity B.C.
u1 = 0;     v1 = 0;
u3 = 0;     v3 = 0;
JBV1 = [BP1(:), [u1,v1] .* ones(length(BP1),1)];
JBV3 = [BP3(:), [u3,v3] .* ones(length(BP3),1)];
JBV=[JBV1; JBV3];

% pressure B.C.
P2 = 0;
P4 = 1000;
JBP2 = [BE2, [P2,P2] .* ones(size(BE2,1),1)];
JBP4 = [BE4, [P4,P4] .* ones(size(BE4,1),1)];
JBP = [JBP2; JBP4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    设定边界条件   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤A结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%     设定初始计算条件   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocityX_init = zeros(noNodesV,1);
velocityY_init = zeros(noNodesV,1);
pressure_init = zeros(noNodesP,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     迭代初始赋值  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocityX_k1 = velocityX_init;
velocityY_k1 = velocityY_init;
pressure_k1  = pressure_init;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B结束   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     迭代初始条件  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergenceCriterionREL = 1e-5;
convergenceCriterionABS = 1e-20;
maxIter = 50;

iter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  开始迭代计算  %%%%%%%%%%
%%%%%    Newton Raphson迭代法求解    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' Beginning iteration, wait... \n')

gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];

while 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.1   %%%%%%%%%%%%%
    %%%%%%%%%      迭代赋值     %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_pre = velocityX_k1;
    velocityY_pre = velocityY_k1;
    pressure_pre  = pressure_k1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.2   %%%%%%%%%%%%%
    %%%%%%%      总体方程子块初始化    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B1 = zeros(noNodesP, noNodesV);
    B2 = zeros(noNodesP, noNodesV);
    D11 = zeros(noNodesV);
    D12 = zeros(noNodesV);
    D21 = zeros(noNodesV);
    D22 = zeros(noNodesV);
    C1 = zeros(noNodesV, noNodesP);
    C2 = zeros(noNodesV, noNodesP);
    F1 = zeros(noNodesV, 1);
    F2 = zeros(noNodesV, 1);
    G1 = zeros(noNodesV, 1);
    G2 = zeros(noNodesV, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.3   %%%%%%%%%
    %%%%%%%   系数矩阵子块计算及组装   %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        
        [Be1, Be2, De11, De12, De21, De22, Ce1, Ce2] = ...
            compute_stiffnessElement(JXYe, mu, gaussPts, gaussWts);
        
        % 组装子矩阵块
        B1(JMP(e,:), JMV(e,:)) = B1(JMP(e,:), JMV(e,:)) + Be1;
        B2(JMP(e,:), JMV(e,:)) = B2(JMP(e,:), JMV(e,:)) + Be2;
        
        D11(JMV(e,:),JMV(e,:)) = D11(JMV(e,:),JMV(e,:)) + De11;
        D12(JMV(e,:),JMV(e,:)) = D12(JMV(e,:),JMV(e,:)) + De12;
        D21(JMV(e,:),JMV(e,:)) = D21(JMV(e,:),JMV(e,:)) + De21;
        D22(JMV(e,:),JMV(e,:)) = D22(JMV(e,:),JMV(e,:)) + De22;
        
        C1(JMV(e,:),JMP(e,:)) = C1(JMV(e,:),JMP(e,:)) + Ce1;
        C2(JMV(e,:),JMP(e,:)) = C2(JMV(e,:),JMP(e,:)) + Ce2;
        
        [Ge1, Ge2] = compute_Ge(JXYe, rho, gaussPts, gaussWts);
        
        G1(JMV(e,:),1) = G1(JMV(e,:),1) + Ge1;
        G2(JMV(e,:),1) = G2(JMV(e,:),1) + Ge2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.4   %%%%%%%%%
    %%%%%%%      压力边界条件代入     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : size(JBP,1)
        JXYe = JXYV(JMV(JBP(i,1),:),:);
        
        [Fe1, Fe2] = compute_forceElement(JXYe,JBP(i,:), gaussPts, gaussWts);
        
        F1(JMV(JBP(i,1),:),1) = F1(JMV(JBP(i,1),:),1) + Fe1;
        F2(JMV(JBP(i,1),:),1) = F2(JMV(JBP(i,1),:),1) + Fe2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.5   %%%%%%%%%
    %%%%%%%   构建NR法J矩阵R向量   %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stiffness = [2*G1*velocityX_pre'+G2*velocityY_pre'+D11,  G2*velocityX_pre'+D12,  -C1;
        G1*velocityY_pre'+D21,    G1*velocityX_pre'+2*G2*velocityY_pre'+D22, -C2;
        B1,                   B2,         zeros(noNodesP, noNodesP)];
    
    R = [(G1*velocityX_pre'+D11)*velocityX_pre + (G2*velocityX_pre'+D12)*velocityY_pre - C1*pressure_pre + F1;
        (G1*velocityY_pre'+D21)*velocityX_pre + (G2*velocityY_pre'+D22)*velocityY_pre - C2*pressure_pre + F2;
        B1*velocityX_pre + B2*velocityY_pre];
    
    force = stiffness*[velocityX_pre;velocityY_pre;pressure_pre] - R;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.6   %%%%%%%%%%%%%%%
    %%%%%%%      代入速度边界条件      %%%%%%%%%
    %%%%%%% 求解方程，更新k+1次迭代结果 %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prescribedDof = cat(1,JBV(:,1),JBV(:,1)+noNodesV);
    prescribedValue = cat(1,JBV(:,2),JBV(:,3));
    % free Dof : activeDof
    activeDof = setdiff( (1:noDofs)', prescribedDof );
    force = force - stiffness(:,prescribedDof) * prescribedValue;
    
    % solution 求解方程
    displacements = zeros(noDofs,1);
    displacements(activeDof) = stiffness(activeDof,activeDof) \ force(activeDof);
    displacements(prescribedDof) = prescribedValue;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.7   %%%%%%%%%
    %%%%%  求解方程，更新k+1次迭代结果   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_k1 = displacements(1 : noNodesV);
    velocityY_k1 = displacements(1+noNodesV : 2*noNodesV);
    pressure_k1  = displacements(1+2*noNodesV:2*noNodesV+noNodesP);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.9   %%%%%%%%%%%%%
    %%%%%%%%%% 误差计算,判断是否收敛 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iter > maxIter    % 牛顿流体，直接赋值ux_k=ux_k_1等
        break
    else
        
        absVX = abs(velocityX_k1 - velocityX_pre);
        absVY = abs(velocityY_k1 - velocityY_pre);
        absP  = abs(pressure_k1  - pressure_pre);
        
        %         isConvergence = max( absVX ) < convergenceCriterionABS ||...
        %             max( absVY ) < convergenceCriterionABS ||...
        %             max( absP ) < convergenceCriterionABS;
        %
        %         if isConvergence
        %             break
        %         else
        
        relVX = max( abs(absVX./velocityX_k1) );
        relVY = max( abs(absVY./velocityY_k1) );
        relP  = max( abs(absP./pressure_k1 ) );
        
        isConvergence = relVX < convergenceCriterionREL ||...
            relVY < convergenceCriterionREL ||...
            relP < convergenceCriterionREL;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%   D.10累加迭代次数，输出迭代结果  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter = iter + 1;
        fprintf('iter. = %d: && rel. err. V_x = %6.9f && rel. err. V_y = %6.9f && rel. err. P= %6.9f \n', iter, relVX, relVY, relP)
        if isConvergence
            break
        end
        %         end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C结束   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 输出结果，绘图 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_k1 = Pding2Pzong(pressure_k1, JMV)';

velocity_k1 = sqrt(velocityX_k1.^2 + velocityY_k1.^2);
data = [JXYV, velocityX_k1, velocityY_k1, velocity_k1, pressure_k1]
JMV4 = JMV_9to4(JMV);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      出口产量计算     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi = gaussPts;
Qfem = 0;
for e = 1 : size(BE2,1)
    
    Pib = JMV(BE2(e,1),[3,6,9]);
    
    x = JXYV(Pib,1);
    y = JXYV(Pib,2);
    u = velocityX_k1(Pib,1);
    Le = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
    for i = 1:length(gaussWts)
        fy=[1/2*xi(i)*(xi(i)-1)
            (1-xi(i))*(1+xi(i))
            1/2*xi(i)*(1+xi(i))];
        Qfem = Qfem + gaussWts(i) * fy' * u * Le/2;
    end
end

Qfem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      出口产量计算     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        雷诺数计算       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uv_max = max(velocity_k1);
Re = rho * uv_max * L / mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      雷诺数计算       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%