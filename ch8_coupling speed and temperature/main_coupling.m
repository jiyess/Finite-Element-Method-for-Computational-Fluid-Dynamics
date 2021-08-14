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
%%%%%%%%%%%%%      迭代步骤A开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   计算参数设定      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_NH = 1;      % 黏性耗散系数 scale_NH=0时，黏性耗散影响消失
theta = 0.2;  % 温度theta因子
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤A结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     物性参数       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0.75;            %幂律指数
mu0 = 10000;         %剪切速率为零是的黏度
muinf = 100;         %剪切速率无穷大时的黏度
nt = 0.5;            %自然时间
rho = 900;       %密度
k = 2;            %导热系数
Cv = 20;          %比热容
alpha = 2000;         %温度系数
T0 = -273;        %绝对零度
Ta = 200;         %基准温度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     读取网格数据  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  区域几何尺寸及网格划分参数
H=0.04; %区域总高
L=0.2;  %区域总长
Nx=20;  %水平方向的网格数量
Ny=5;   %竖直方向的网格数量
%%%%%%%  区域几何尺寸及网格划分参数
[JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_coupling(H, L, Nx, Ny)
%%%%%%% 调用四边形网格绘制程序
plotRectangularGrid(JMP,JXYV,noElement, noNodesV, 0);
%%%%%%% 调用四边形网格绘制程序

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     设定边界条件    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = 0; v1 = 0; %速度边界
u3 = 0; v3 = 0;
JBV1 = [BP1(:), [u1,v1] .* ones(length(BP1),1)];
JBV3 = [BP3(:), [u3,v3] .* ones(length(BP1),1)];
JBV  = [JBV1;JBV3];

P2=0;  %压力边界
P4=10000;
JBP2=[BE2,[P2,P2] .* ones(size(BE2,1),1)];
JBP4=[BE4,[P4,P4] .* ones(size(BE4,1),1)];
JBP=[JBP2;JBP4];

T4 = 180; %与时间无关的温度边界条件
JBT14 = [BP4(:), T4*ones(length(BP4),1)];
JBT1 = JBT14;

q1 = -1000;% 热流密度             % 第二类温度边界条件数据，
q2 = 0;% 热流密度 q2=0,表示绝热  % 上下壁面加热热流密度1000，
q3 = -1000;% 热流密度             % 出口绝热
normalDers = -q1/k * ones(size(BE1(:,1)));
JBT21 = [BE1, normalDers, normalDers, normalDers];
normalDers = -q2/k * ones(size(BE2(:,1)));
JBT22 = [BE2, normalDers, normalDers, normalDers];
normalDers = -q3/k * ones(size(BE3(:,1)));
JBT23 = [BE3, normalDers, normalDers, normalDers];
JBT2 = [JBT21; JBT23; JBT22];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      调用初始条件  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_k1_coupling.mat     % 调用相同边界条件下牛顿流体速度分布结果ux0,vy0,p0
load result_CDT_coupling.mat    % 调用相同边界条件下纯导热问题温度分布结果T0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     进行初始赋值  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocityX_k1 = zeros(noNodesV,1);
% velocityY_k1 =  zeros(noNodesV,1);
% pressure_k1 =  zeros(noNodesV,1);
% 
% temp_k1 = zeros(noNodesV,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤D开始    %%%%%%%%%%%%
%%%%%%%%%%%%%      初始黏度计算     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据

for e = 1:noElement
    
    JXYe = JXYV(JMV(e,:),:);
    velocityE(:,1) = velocityX_k1(JMV(e,:),:);
    velocityE(:,2) = velocityY_k1(JMV(e,:),:);
    
    muE = compute_muE_BirdCarreau(n, mu0, muinf, nt, JXYe, velocityE);
    % 调用程序,根据Bird-Carreau模型，计算单元内结点剪切速率和黏度
    
    avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + muE;
    avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
end
mu_k1 = avgMu(:,1) ./ avgMu(:,2);
mu_k1 = mu_k1 .* exp(alpha*(1./(temp_k1-T0)-1/(Ta-T0)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤D结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤E开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     迭代初始条件  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergenceCriterionREL = 1e-3;
maxIter = 100;

iter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  开始迭代计算  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];

fprintf(' Beginning iteration, wait... \n')

while 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤E.1   %%%%%%%%%
    %%%%%%%   将k+1结果赋值给k %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_pre = velocityX_k1;
    velocityY_pre = velocityY_k1;
    pressure_pre = pressure_k1;
    mu_pre = mu_k1;
    temp_pre = temp_k1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.2   %%%%%%%%%%
    %%%%%%%  初始化总体方程子块   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    DL = zeros(noNodesV);
    CD = zeros(noNodesV);
    CDB = zeros(noNodesV, 1);
    NH = zeros(noNodesV, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤E.3   %%%%%%%%%
    %%%  计算NS方程系数矩阵单元子块并组装 %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.4   %%%%%%%%%
    %%%  计算NS方程右边向量单元子块并组装 %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : size(JBP,1)
        JXYe = JXYV(JMV(JBP(i,1),:),:);
        
        [Fe1, Fe2] = compute_forceElement(JXYe,JBP(i,:), gaussPts, gaussWts);
        
        F1(JMV(JBP(i,1),:),1) = F1(JMV(JBP(i,1),:),1) + Fe1;
        F2(JMV(JBP(i,1),:),1) = F2(JMV(JBP(i,1),:),1) + Fe2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.5   %%%%%%%%%%
    %%%%%%%   构建NS方程总体方程  %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stiffness = [D11, D12, -C1;
        D21, D22, -C2;
        B1,  B2,  zeros(noNodesP)];
    force = [-F1; -F2; zeros(noNodesP,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.6   %%%%%%%%%
    %%%%%%%   代入速度边界条件   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%     迭代步骤E.7   %%%%%%%%%
    %%%%%  求解方程，更新k+1次迭代结果   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_k1 = displacements(1 : noNodesV);
    velocityY_k1 = displacements(1+noNodesV : 2*noNodesV);
    
    pressure4 = displacements(1+2*noNodesV:2*noNodesV+noNodesP);
    pressure_k1 = Pding2Pzong(pressure4, JMV)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.8    %%%%%%%%%
    %%%     当前速度和上一步温度更新黏度    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据
    %     avgSR = zeros(noNodesV, 2);  % 建立结点剪切速率数据
    
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        velocityE(:,1) = velocityX_k1(JMV(e,:),:);
        velocityE(:,2) = velocityY_k1(JMV(e,:),:);
        
        muE = compute_muE_BirdCarreau(n, mu0, muinf, nt, JXYe, velocityE);
        %调用程序计算单元内结点剪切速率和黏度
        
        avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + muE;
        avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
        %         avgSR(JMV(e,:),1) = avgSR(JMV(e,:),1) + SRE;
        %         avgSR(JMV(e,:),2) = avgSR(JMV(e,:),2) + 1;
    end
    
    mu_k1 = avgMu(:,1) ./ avgMu(:,2);
    mu_k1 = mu_k1 .* exp(alpha*(1./(temp_k1-T0)-1/(Ta-T0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.9   %%%%%%%%%
    %%%%%  当前速度及更新黏度计算能量    %%%%%
    %%%%%  单元方程系数矩阵子块并组装   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        velocityE(:,1) = velocityX_k1(JMV(e,:),:);
        velocityE(:,2) = velocityY_k1(JMV(e,:),:);
        
        MuE = mu_k1(JMV(e,:),:);
        %调用程序计算单元内结点剪切速率和黏度
        
        CDe = compute_CDe(JXYe, k, gaussPts, gaussWts);   %调用CDe子块计算程序
        DLe = compute_DLe(JXYe, velocityE, rho, Cv, gaussPts, gaussWts);%调用DLe子块计算程序
        NHe = compute_NHe(JXYe, velocityE, muE, gaussPts, gaussWts); %调用NHe子块计算程序
        
        CD(JMV(e,:),JMV(e,:)) = CD(JMV(e,:),JMV(e,:)) + CDe;
        DL(JMV(e,:),JMV(e,:)) = DL(JMV(e,:),JMV(e,:)) + DLe;
        NH(JMV(e,:), 1) = NH(JMV(e,:), 1) + NHe;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.10   %%%%%%%%%
    %%%%%%%   代入JBT2数据计算能量   %%%%%
    %%%%%%% 单元方程右边向量子块并组装 %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : size(JBT2, 1)
        
        JXYe = JXYV(JMV(JBT2(i,1),:),:);
        JBT2e = JBT2(i,:);
        
        CDBe = compute_CDBe(JXYe, JBT2e,  gaussPts, gaussWts);
        
        CDB(JMV(JBT2(i,1),:),1) = CDB(JMV(JBT2(i,1),:),1) + CDBe;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.11   %%%%%%%%%
    %%%%%%   构建能量方程总体方程    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stiffnessT = CD + (1-theta) * DL; %能量方程系数矩阵和右边向量
    forceT = NH*scale_NH + CDB - theta *DL*temp_pre;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.12   %%%%%%%%%
    %%%%%%%%%    代入JBT1数据   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prescribedDofT = JBT1(:,1);
    prescribedValueT = JBT1(:,2);
    % free Dof : activeDof
    activeDofT = setdiff( (1:noNodesV)', prescribedDofT );
    forceT = forceT - stiffnessT(:,prescribedDofT) * prescribedValueT;
    
    % solution 求解方程
    temp_k1 = zeros(noNodesV,1);
    temp_k1(activeDofT) = stiffnessT(activeDofT, activeDofT) \ forceT(activeDofT);
    temp_k1(prescribedDofT) = prescribedValueT;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.14   %%%%%%%%%
    %%%%%% 当前速度和温度，更新黏度  %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据
    
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        velocityE(:,1) = velocityX_k1(JMV(e,:),:);
        velocityE(:,2) = velocityY_k1(JMV(e,:),:);
        
        muE = compute_muE_BirdCarreau(n, mu0, muinf, nt, JXYe, velocityE);
        % 调用程序,根据Bird-Carreau模型，计算单元内结点剪切速率和黏度
        
        avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + muE;
        avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
    end
    mu_k1 = avgMu(:,1) ./ avgMu(:,2);
    mu_k1 = mu_k1 .* exp(alpha*(1./(temp_k1-T0)-1/(Ta-T0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤E.15   %%%%%%%%%
    %%%%%%%%%% 误差计算,判断是否收敛 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iter > maxIter    % 牛顿流体，直接赋值ux_k=ux_k_1等
        break
    else
        
        absVX = abs(velocityX_k1 - velocityX_pre);
        absVY = abs(velocityY_k1 - velocityY_pre);
        absP  = abs(pressure_k1  - pressure_pre);
        absMu = abs(mu_k1  - mu_pre);
        absT = abs(temp_k1 - temp_pre);
        
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
        relT  = max( abs(absT./temp_k1));
        
        isConvergence = relVX < convergenceCriterionREL ||...
            relVY < convergenceCriterionREL ||...
            relP < convergenceCriterionREL ||...
            relMu < convergenceCriterionREL ||...
            relT < convergenceCriterionREL;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%   D.10累加迭代次数，输出迭代结果  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter = iter + 1;
        fprintf('iter. = %d: && rel. err. V_x = %6.9f && rel. err. V_y = %6.9f && rel. err. P= %6.9f && rel. err. mu = %6.9f && rel. err. Temp = %6.9f \n', iter, relVX, relVY, relP, relMu, relT)
        
        if isConvergence
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
