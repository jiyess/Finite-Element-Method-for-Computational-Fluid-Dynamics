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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     物性参数       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0.75;            % 幂律指数
mu0 = 10000;         % 剪切速率为零时的黏度
muinf = 100;         % 剪切速率无穷大时的黏度
nt = 0.5;            % 自然时间
rho = 900;           % 密度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     读取网格数据   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     时间步长相关   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t = 1;    % 时间步长
t = 0;          % 初始时间
time_up = 50;   % 时间上限
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        平板拖动速度初始设定      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u3 = 0;     % 平板拖动初始速度
du3 = 0.1;  % 平板拖动加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤A结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤B开始    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    定义初始速度    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocityX_init = zeros(noNodesV,1);% 初始速度、压力均设定为0
velocityY_init = zeros(noNodesV,1);
pressure_init  = zeros(noNodesP,1);

velocityX_k1 = velocityX_init;
velocityY_k1 = velocityY_init;
pressure_k1 = pressure_init;

velocity_k1 = sqrt(velocityX_k1.^2 + velocityY_k1.^2); % 法向速度计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    初始粘度计算     %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     迭代步骤C.4   %%%%%%%%%
%%%%%%%% 黏度和剪切速率累加   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_k1 = avgMu(:,1) ./ avgMu(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     结果记录相关设定     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tag = 0;           % 初始化数据记录序号
detla_tag = 1;     % 记录步长, 每隔detla_tag记录一次
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      记录初始分布结果     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tag = tag + 1;        % 初始速度、压力、黏度分布记录
velocityX(:,tag) = velocityX_k1;
velocityY(:,tag) = velocityY_k1;
pressure(:,tag) = Pding2Pzong(pressure_k1,JMV)';
velocity(:,tag) = velocity_k1;
mu(:,tag) = mu_k1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      计算初始平均速度     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];

iter = 0;

INTV = 0;   % 初始化区域速度积分
AREA = 0;   % 初始化区域面积
for e = 1:noElement
    
    JXYe = JXYV(JMV(e,:),:);
    velocity_k1e = velocity_k1(JMV(e,:),:);
    
    [INTVe,AREAe] = compute_INTVe_AREAe(JXYe, velocity_k1e, gaussPts, gaussWts);
    
    INTV = INTV + INTVe;
    AREA = AREA + AREAe;
    
end
avgVelocity(iter+1) = INTV / AREA; % 区域平均速度
iterStep(iter+1) = 0;  % 记录平均速度对应时间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      迭代步骤B结束    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤C开始     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     非定常迭代计算开始条件   %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergenceCriterionREL_out = 1e-3;
convergenceCriterionREL_inner = 1e-3;
% convergenceCriterionABS = 1e-10;
% maxIter = 50;
maxIter_out = time_up;
maxIter_inner = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     开始非定常问题迭代求解   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Beginning iteration, wait... \n')

while 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤C.1   %%%%%%%%%
    %%%%%%%%%     速度赋值      %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_pre = velocityX_k1;
    velocityY_pre = velocityY_k1;
    pressure_pre = pressure_k1;
    mu_pre = mu_k1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%    迭代步骤C.2     %%%%%%%%%
    %%%%%%%%       设定边界条件      %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u1=0;v1=0;
    u2=0;v2=0;
    t=t+delta_t; %上平板拖动速度与时间t有关
    if t <= 10
        u3 = t*du3;
        v3=0;
    elseif t<=20 && t>10
        u3 = 1-(t-10)*du3;
        v3=0;
    else
        u3=0;
        v3=0;
    end
    u4 = 0; v4 = 0;
    JBV1 = [BP1(:), [u1,v1] .* ones(length(BP1),1)];
    JBV2 = [BP2(:), [u2,v2] .* ones(length(BP2),1)];
    JBV3 = [BP3(:), [u3,v3] .* ones(length(BP3),1)];
    JBV4 = [BP4(:), [u4,v4] .* ones(length(BP4),1)];
    JBV=[JBV1; JBV2; JBV3; JBV4];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤C.3    %%%%%%%%%
    %%%%%%     非牛顿迭代初始赋值    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_k1_inner = velocityX_pre;
    velocityY_k1_inner = velocityY_pre;
    pressure_k1_inner = pressure_pre;
    mu_k1_inner = mu_pre;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤C.4   %%%%%%%%
    %%%%%%     非牛顿迭代初始条件    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter_inner = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%   迭代步骤C.4   %%%%%%%%%
    %%%%%%      非牛顿迭代计算     %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤C.4.1   %%%%%%%%%
        %%%       当前n+1时刻非牛顿迭代赋值    %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        velocityX_pre_inner = velocityX_k1_inner;
        velocityY_pre_inner = velocityY_k1_inner;
        pressure_pre_inner = pressure_k1_inner;
        mu_pre_inner = mu_k1_inner;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%  迭代步骤C.4.2    %%%%%%%%%
        %%%%%%       总体方程各项初始化    %%%%%%
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
        M = zeros(noNodesV);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%   迭代步骤C.4.3   %%%%%%%%%
        %%%     系数矩阵单元子块计算及组合       %%%%
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
            
            Me = compute_Me(JXYe, rho, 1, gaussPts, gaussWts);
            M(JMV(e,:),JMV(e,:)) = M(JMV(e,:),JMV(e,:)) + Me;
        end
        
%         for i = 1 : size(JBP,1)
%             JXYe = JXYV(JMV(JBP(i,1),:),:);
%             
%             [Fe1, Fe2] = compute_forceElement(JXYe,JBP(i,:), gaussPts, gaussWts);
%             
%             F1(JMV(JBP(i,1),:),1) = F1(JMV(JBP(i,1),:),1) + Fe1;
%             F2(JMV(JBP(i,1),:),1) = F2(JMV(JBP(i,1),:),1) + Fe2;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%   迭代步骤C.4.4   %%%%%%%%%
        %%%%%%%%%    总体方程集成     %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        stiffness =[M+delta_t*D11,  delta_t*D12,    -C1*delta_t;
            delta_t*D21,   M+delta_t*D22,   -C2*delta_t;
            B1,   B2, zeros(noNodesP)];
        force = [M*velocityX_pre - F1*delta_t;
            M*velocityY_pre - F2*delta_t;
            zeros(noNodesP,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%  迭代步骤C.4.5    %%%%%%%%%
        %%%%%%%      代入速度边界条件      %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noDofs = 2*noNodesV + noNodesP;
        
        prescribedDof = cat(1,JBV(:,1),JBV(:,1)+noNodesV);
        prescribedValue = cat(1,JBV(:,2),JBV(:,3));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%   迭代步骤C.4.6   %%%%%%%%%
        %%%%%%%%%     指定压力零点   %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prescribedDof = [prescribedDof; 2*noNodesV + 30];
        prescribedValue = [prescribedValue; 0];
        
        % free Dof : activeDof
        activeDof = setdiff( (1:noDofs)', prescribedDof );
        force = force - stiffness(:,prescribedDof) * prescribedValue;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%   迭代步骤C.4.7  %%%%%%%%%
        %%%%%%%%%      求解方程    %%%%%%%%%
        %%%   更新n+1时刻第k+1次非牛顿迭代结果   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % solution 求解方程
        displacements = zeros(noDofs,1);
        displacements(activeDof) = stiffness(activeDof,activeDof) \ force(activeDof);
        displacements(prescribedDof) = prescribedValue;
        
        
        velocityX_k1_inner = displacements(1 : noNodesV);
        velocityY_k1_inner = displacements(1+noNodesV : 2*noNodesV);
        pressure_k1_inner = displacements(1+2*noNodesV:2*noNodesV+noNodesP);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%   迭代步骤C.4.8   %%%%%%%%%
        %%%%%%%   更新n+1时刻第k+1次    %%%%%%
        %%%%%%%     非牛顿迭代黏度结果     %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        avgMu = zeros(noNodesV, 2);   % 建立结点黏度数据
        
        for e = 1:noElement
            
            JXYe = JXYV(JMV(e,:),:);
            velocityE(:,1) = velocityX_k1_inner(JMV(e,:),:);
            velocityE(:,2) = velocityY_k1_inner(JMV(e,:),:);
            
            muE = compute_muE_BirdCarreau(n, mu0, muinf, nt, JXYe, velocityE);
            % 调用程序,根据Bird-Carreau模型，计算单元内结点剪切速率和黏度
            
            avgMu(JMV(e,:),1) = avgMu(JMV(e,:),1) + muE;
            avgMu(JMV(e,:),2) = avgMu(JMV(e,:),2) + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%     迭代步骤C.4   %%%%%%%%%
        %%%%%%%% 黏度和剪切速率累加   %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mu_k1_inner = avgMu(:,1) ./ avgMu(:,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤C.4.9  %%%%%%%%%
        %%%%%%     非牛顿问题误差计算      %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%     迭代步骤D.9   %%%%%%%%%%%%%
        %%%%%%%%%% 误差计算,判断是否收敛 %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iter_inner > maxIter_inner    % 牛顿流体，直接赋值ux_k=ux_k_1等
            break
        else
            
            absVX = abs(velocityX_k1_inner - velocityX_pre_inner);
            absVY = abs(velocityY_k1_inner - velocityY_pre_inner);
            absP  = abs(pressure_k1_inner  - pressure_pre_inner);
            absMu = abs(mu_k1_inner  - mu_pre_inner);
            
            %         isConvergence = max( absVX ) < convergenceCriterionABS ||...
            %             max( absVY ) < convergenceCriterionABS ||...
            %             max( absP ) < convergenceCriterionABS ||...
            %             max( absMu ) < convergenceCriterionABS;
            %
            %         if isConvergence
            %             break
            %         else
            
            relVX = max( abs(absVX./velocityX_k1_inner) );
            relVY = max( abs(absVY./velocityY_k1_inner) );
            relP  = max( abs(absP./pressure_k1_inner ) );
            relMu = max( abs(absMu./mu_k1_inner) );
            
            isConvergence = relVX < convergenceCriterionREL_inner ||...
                relVY < convergenceCriterionREL_inner ||...
                relP < convergenceCriterionREL_inner ||...
                relMu < convergenceCriterionREL_inner;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%   D.10累加迭代次数，输出迭代结果  %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            iter_inner = iter_inner + 1;
            fprintf('   Inner iter.:  iter. = %d: && rel. err. V_x_inner = %6.9f && rel. err. V_y_inner = %6.9f && rel. err. P_inner= %6.9f && rel. err. mu_inner = %6.9f \n', iter_inner, relVX, relVY, relP, relMu)
            
            if isConvergence
                break
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%   迭代步骤C.5   %%%%%%%%%
    %%%%%%%%   当前时刻结果赋值   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    velocityX_k1 = velocityX_k1_inner;  % 当前时刻的收敛结果
    velocityY_k1 = velocityY_k1_inner;
    velocity_k1 = sqrt(velocityX_k1.^2 + velocityY_k1.^2);
    pressure_k1 = pressure_k1_inner;
    mu_k1 = mu_k1_inner;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤C.6    %%%%%%%%%
    %%%%%%%%%%%   数据记录   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(iter,detla_tag) == 0  % 每个记录步长存储一次结果
        tag = tag+1;
        velocityX(:,tag) = velocityX_k1;
        velocityY(:,tag) = velocityY_k1;
        pressure(:,tag) = Pding2Pzong(pressure_k1,JMV)';
        velocity(:,tag) = velocity_k1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%  迭代步骤C.7     %%%%%%%%%
    %%%%    累加非定常计算时间步长个数   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter = iter+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤C.8   %%%%%%%%%
    %%%%%%%    计算平均温度   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INTV = 0;   % 初始化区域速度积分
    AREA = 0;   % 初始化区域面积
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        velocity_k1e = velocity_k1(JMV(e,:),:);
        
        [INTVe,AREAe] = compute_INTVe_AREAe(JXYe, velocity_k1e, gaussPts, gaussWts);
        
        INTV = INTV + INTVe;
        AREA = AREA + AREAe;
        
    end
    avgVelocity(iter+1) = INTV / AREA; % 区域平均速度
    iterStep(iter+1) = iter*delta_t;  % 记录平均速度对应时间
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤C.9   %%%%%%%%%
    %%%%%%%  非定常问题误差计算   %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter > maxIter_out   % 牛顿流体，直接赋值ux_k=ux_k_1等
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
        
        isConvergence = relVX < convergenceCriterionREL_out ||...
            relVY < convergenceCriterionREL_out ||...
            relP < convergenceCriterionREL_out ||...
            relMu < convergenceCriterionREL_out;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%   D.10累加迭代次数，输出迭代结果  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         iter = iter + 1;
        fprintf('Outer iter.:iter. = %d: && rel. err. V_x = %6.9f && rel. err. V_y = %6.9f && rel. err. P= %6.9f && rel. err. mu = %6.9f \n', iter, relVX, relVY, relP, relMu)
        
        if isConvergence
            break
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     迭代步骤C结束    %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     绘制平均温度变化曲线   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(iterStep, avgVelocity);
shoulian=[iterStep(:), avgVelocity(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    绘制平均温度变化曲线   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

