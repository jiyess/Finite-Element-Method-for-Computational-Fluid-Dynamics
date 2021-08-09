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
warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤A开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     物性参数      %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0.5;    % 幂律指数，先设定n=1然后依次计算,n=0.9、n=0.7和n=0.5的结果
m = 1000;   % 稠度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     读取网格数据   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 0.01;    % 平行板间高度
L = 0.02;    % 平行板长度
Nx = 8;     % 水平方向的网格数量
Ny = 8;     % 竖直方向的网格数量
%%%%%%%  区域几何尺寸及网格划分参数
[JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_rectangle(H, L, Nx, Ny);
%%%%%%% 调用四边形网格绘制程序
plotRectangularGrid(JMP,JXYV,noElement, noNodesV, 1);
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
P4 = 10000;
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
%%%%% % %   读取迭代开始结果   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch n
    case 0.9
        load result_k1.mat
    case 0.7
        load result_k0p9.mat
    case 0.5
        load result_k0p7.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B结束   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C开始    %%%%%%%%%%%%
%%%%%%%%%%%%%       计算初始黏度    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n == 1
    mu_k1 = ones(noNodesV, 1) * m; % 取$\mu_{k+1}$为幂律模型稠度m
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.1   %%%%%%%%%
    %%%%%  初始化结点剪切速率和黏度数据   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据
    %     avgSR = zeros(noNodesV, 2);  % 建立结点剪切速率数据
    
    for e = 1:noElement
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%     迭代步骤C.2   %%%%%%%%%
        %%%%%%%  提取单元结点坐标和速度   %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        JXYe = JXYV(JMV(e,:),:);
        velocityE(:,1) = velocityX_k1(JMV(e,:),:);
        velocityE(:,2) = velocityY_k1(JMV(e,:),:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%     迭代步骤C.3   %%%%%%%%%
        %%%%%%  计算单元结点黏度和剪切速率   %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MuE = compute_MuE(m, n, JXYe, velocityE);
        %调用程序计算单元内结点剪切速率和黏度
        
        avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + MuE;
        avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
        %         avgSR(JMV(e,:),1) = avgSR(JMV(e,:),1) + SRE;
        %         avgSR(JMV(e,:),2) = avgSR(JMV(e,:),2) + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤C.4   %%%%%%%%%
    %%%%%%%% 黏度和剪切速率累加   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu_k1 = avgMu(:,1) ./ avgMu(:,2);
    mu_k1(mu_k1 < 1e-10) = 1e-10;
    mu_k1(mu_k1 > 1e10) = 1e10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤D开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     迭代初始条件   %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergenceCriterionREL = 1e-3;
% convergenceCriterionABS = 1e-10;
maxIter = 50;

iter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     开始迭代计算   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Beginning iteration, wait... \n')

gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];

while 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.1   %%%%%%%%%
    %%%%%%%%%      迭代赋值     %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1
        mu_pre = mu_k1;
    else
        velocityX_pre = velocityX_k1;
        velocityY_pre = velocityY_k1;
        pressure_pre = pressure_k1;
        mu_pre = mu_k1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.2   %%%%%%%%%%%%
    %%%%%%%      总体方程子块初始化    %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.3   %%%%%%%%%%%%
    %%%%% 单元方程系数矩阵子块计算及组装   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e=1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        mue = mu_pre(JMV(e,:),1);
        
        [Be1, Be2, De11, De12, De21, De22, Ce1, Ce2] = ...
            compute_stiffnessElement(JXYe, mue, gaussPts, gaussWts);
        
        % 组装子矩阵块
        B1(JMP(e,:), JMV(e,:)) = B1(JMP(e,:), JMV(e,:)) + Be1;
        B2(JMP(e,:), JMV(e,:)) = B2(JMP(e,:), JMV(e,:)) + Be2;
        
        D11(JMV(e,:),JMV(e,:)) = D11(JMV(e,:),JMV(e,:)) + De11;
        D12(JMV(e,:),JMV(e,:)) = D12(JMV(e,:),JMV(e,:)) + De12;
        D21(JMV(e,:),JMV(e,:)) = D21(JMV(e,:),JMV(e,:)) + De21;
        D22(JMV(e,:),JMV(e,:)) = D22(JMV(e,:),JMV(e,:)) + De22;
        
        C1(JMV(e,:),JMP(e,:)) = C1(JMV(e,:),JMP(e,:)) + Ce1;
        C2(JMV(e,:),JMP(e,:)) = C2(JMV(e,:),JMP(e,:)) + Ce2;
    end
    
    
    %%%%%%%%%     迭代步骤D.4   %%%%%%%%%
    %%%%  代入JBP数据计算Fe子块，并组装   %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : size(JBP,1)
        JXYe = JXYV(JMV(JBP(i,1),:),:);
        
        [Fe1, Fe2] = compute_forceElement(JXYe,JBP(i,:), gaussPts, gaussWts);
        
        F1(JMV(JBP(i,1),:),1) = F1(JMV(JBP(i,1),:),1) + Fe1;
        F2(JMV(JBP(i,1),:),1) = F2(JMV(JBP(i,1),:),1) + Fe2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.5   %%%%%%%%%
    %%%%%%%%%   构建总体计算方程   %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stiffness = [D11, D12, -C1;
        D21, D22, -C2;
        B1,  B2,  zeros(noNodesP)];
    force = [-F1; -F2; zeros(noNodesP,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%     迭代步骤D.6   %%%%%%%%%
    %%%%%%%      代入速度边界条件      %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    noDofs = 2*noNodesV + noNodesP;
    
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
    %%%%%%%%%     迭代步骤D.7   %%%%%%%%%
    %%%%%  求解方程，更新k+1次迭代结果   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_k1 = displacements(1 : noNodesV);
    velocityY_k1 = displacements(1+noNodesV : 2*noNodesV);
    
    pressure4 = displacements(1+2*noNodesV:2*noNodesV+noNodesP);
    pressure_k1 = Pding2Pzong(pressure4, JMV)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.8   %%%%%%%%%
    %%%%%%%%%%%    更新黏度   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据
    %     avgSR = zeros(noNodesV, 2);  % 建立结点剪切速率数据
    
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        velocityE(:,1) = velocityX_k1(JMV(e,:),:);
        velocityE(:,2) = velocityY_k1(JMV(e,:),:);
         
        MuE = compute_MuE(m, n, JXYe, velocityE);
        %调用程序计算单元内结点剪切速率和黏度
        
        avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + MuE;
        avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
        %         avgSR(JMV(e,:),1) = avgSR(JMV(e,:),1) + SRE;
        %         avgSR(JMV(e,:),2) = avgSR(JMV(e,:),2) + 1;
    end
    
    mu_k1 = avgMu(:,1) ./ avgMu(:,2);
    mu_k1(mu_k1 < 1e-10) = 1e-10;
    mu_k1(mu_k1 > 1e10) = 1e10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.9   %%%%%%%%%%%%%
    %%%%%%%%%% 误差计算,判断是否收敛 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1 || iter > maxIter    % 牛顿流体，直接赋值ux_k=ux_k_1等
        break
    else
        
        absVX = abs(velocityX_k1 - velocityX_pre);
        absVY = abs(velocityY_k1 - velocityY_pre);
        absP  = abs(pressure_k1  - pressure_pre);
        absMu = abs(mu_k1  - mu_pre);
        
        %         isConvergence = max( absVX ) < convergenceCriterionABS ||...
        %             max( absVY ) < convergenceCriterionABS ||...
        %             max( absP ) < convergenceCriterionABS ||...
        %             max( absMu ) < convergenceCriterionABS;
        %
        %         if isConvergence
        %             break
        %         else
        
        relVX = max( abs(absVX./velocityX_k1) );
        relVY = max( abs(absVY./velocityY_k1) );
        relP  = max( abs(absP./pressure_k1 ) );
        relMu = max( abs(absMu./mu_k1) );
        
        isConvergence = relVX < convergenceCriterionREL ||...
            relVY < convergenceCriterionREL ||...
            relP < convergenceCriterionREL ||...
            relMu < convergenceCriterionREL;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%   D.10累加迭代次数，输出迭代结果  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter = iter + 1;
        fprintf('iter. = %d: && rel. err. V_x = %6.9f && rel. err. V_y = %6.9f && rel. err. P= %6.9f && rel. err. mu = %6.9f \n', iter, relVX, relVY, relP, relMu)
        
        if isConvergence
            break
        end
    end
end % compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤D结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 输出结果，绘图 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocity_k1 = sqrt(velocityX_k1.^2 + velocityY_k1.^2);
% data = [JXYV, velocityX_k1, velocityY_k1, velocity_k1, pressure_k1, mu_k1]
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

figure
hold on
for i = 1:noElement
    patch(JXYV(JMP(i,:),1),JXYV(JMP(i,:),2),mu_k1(JMP(i,:)), 'edgecolor','none')
end
colorbar
title('mu\_k1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      出口速度分布对比        %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocityX_B2 = velocityX_k1(BP2(:),1);
y_B2 = JXYV(BP2(:),2);
figure
hold on
plot(y_B2-H/2,velocityX_B2,'b*')
% V_exit_FEM = [y_B2(:), velocityX_B2(:)]

yy = 0:H/100:H/2;    %  精确解速度分布
exactVelocityX_B2 = n/(n+1)*(P4/m/L)^(1/n)*((H/2)^((n+1)/n)-yy.^((n+1)/n));
plot(yy, exactVelocityX_B2, 'r')
% V_exit_JQ = [yy(:), exactVelocityX_B2(:)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     出口速度分布对比        %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     出口产量计算        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% 数值解出口流量计算
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
%%%% 精确出口流量计算
Qexact = n*H^2/(2*(2*n+1))*(H*P4/2/m/L)^(1/n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        出口产量计算        %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    清除多余变量并存储计算结果    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch n
    case 1
        save('result_k1.mat', 'JMP', 'JXYV', 'velocityX_k1', 'velocityY_k1',...
            'velocity_k1', 'pressure_k1')
    case 0.9
        save('result_k0p9.mat', 'JMP', 'JXYV', 'velocityX_k1', 'velocityY_k1',...
            'velocity_k1', 'pressure_k1')
    case 0.7
        save('result_k0p7.mat', 'JMP', 'JXYV', 'velocityX_k1', 'velocityY_k1',...
            'velocity_k1', 'pressure_k1')
    case 0.5
        save('result_k0p5.mat', 'JMP', 'JXYV', 'velocityX_k1', 'velocityY_k1',...
            'velocity_k1', 'pressure_k1')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    清除多余变量并存储计算结果    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%