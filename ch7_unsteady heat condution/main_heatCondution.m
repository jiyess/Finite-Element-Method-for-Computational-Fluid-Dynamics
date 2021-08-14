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
%%%%%%%%%%%%%       迭代步骤A开始     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      物性参数       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 900;         % 密度
k = 2;             % 导热系数
Cv = 20;           % 比热容
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        边界条件参数        %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 10;           % 空气自然对流换热系数
Tf = 5;               % 空气温度；
q1 = -100;           % 加热热流密度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      读取网格数据   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  区域几何尺寸及网格划分参数
H = 0.1;   % 区域总高
L = 0.1;   % 区域总长
Nx = 8;     % 水平方向的网格数量
Ny = 8;     % 竖直方向的网格数量
%%%%%%%  区域几何尺寸及网格划分参数
[JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    gird_generation_heatCondution(H, L, Nx, Ny);

%%%%%%% 调用四边形网格绘制程序
plotRectangularGrid(JMP,JXYV,noElement, noNodesV, 1);
%%%%%%% 调用四边形网格绘制程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     时间步长相关    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t = 50;       % 时间步长
time_up = 100000;   % 时间上限
t = 0;              % 初始时间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤A结束     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤B开始     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   定义初始温度T0，并赋值给Tn+1  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp0 = 20;
temp_k1 = ones(noNodesV,1) * temp0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤B结束     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤C开始     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     非定常问题迭代计算开始条件  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];

iter = 0;       % 迭代次数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    计算初始温度平均值    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INTT = 0;   % 初始化区域速度积分
AREA = 0;   % 初始化区域面积
for e = 1:noElement
    
    JXYe = JXYV(JMV(e,:),:);
    temp_k1e = temp_k1(JMV(e,:),:);
    
    [INTTe, AREAe] = compute_INTVe_AREAe(JXYe, temp_k1e, gaussPts, gaussWts);
    
    INTT = INTT + INTTe;
    AREA = AREA + AREAe;
    
end

avgT(iter+1) = INTT/AREA; % 区域平均温度
iterStep(iter+1) = iter * delta_t; % 记录平均温度对应时间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   初始化数据记录参数    %%%%%%%
%%%%%%%%%     并记录初始数据    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag = 0;           % 初始化数据记录序号
detla_tag = 2;     % 记录步长, 每隔detla_tag记录一次
tag = tag+1;
Temp(:,tag) = temp_k1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤C结束     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      迭代步骤D开始      %%%%%%%%%%%%
%%%%%%%%%%%%%     非定常迭代计算开始条件   %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergenceCriterionREL_out = 1e-3;
convergenceCriterionREL_inner = 1e-3;
% convergenceCriterionABS = 1e-10;
% maxIter = 50;
maxIter_out = time_up/delta_t;
maxIter_inner = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     开始非定常问题迭代求解   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Beginning iteration, wait... \n')

while 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%    迭代步骤D.1    %%%%%%%%%
    %%%%%%%%%   将Tn+1赋值给Tn  %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_pre = temp_k1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%    迭代步骤D.2     %%%%%%%%%
    %%%%%%%%  将Tn赋值给Tn+1k+1   %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    temp_k1_inner = temp_pre;     % 利用上一时刻收敛结果开始计算
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.3    %%%%%%%%%
    %%%%%%%    对流换热迭代初始条件    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter_inner=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.4     %%%%%%%%%
    %%%%%%%  开始对流换热迭代计算   %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤D.4.1    %%%%%%%%%
        %%%%%%%%%   更新温度边界条件   %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_pre_inner = temp_k1_inner;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%   迭代步骤D.4.2     %%%%%%%%%
        %%%%%%%%  更新温度边界条件   %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        normalDers = -q1/k * ones(size(BE1(:,1)));
        JBT21=[BE1, normalDers, normalDers, normalDers];
        
        JBT22 = zeros(size(BE2,1),7);
        JBT23 = zeros(size(BE3,1),7);
        JBT24 = zeros(size(BE4,1),7);
        for i = 1:size(BE2,1)
            EB = BE2(i,1);
            PB = [JMV(EB,3),JMV(EB,6),JMV(EB,9)];
            TB = [temp_pre_inner(PB(1),1), temp_pre_inner(PB(2),1), temp_pre_inner(PB(3),1)];
            q=-1*alpha*(TB-Tf)/k;
            JBT22(i,:)=[BE2(i,:),q];
        end
        
        for i = 1:size(BE3,1)
            EB = BE3(i,1);
            PB = [JMV(EB,7),JMV(EB,8),JMV(EB,9)];
            TB = [temp_pre_inner(PB(1),1),temp_pre_inner(PB(2),1),temp_pre_inner(PB(3),1)];
            q = -1*alpha*(TB-Tf)/k;
            JBT23(i,:) = [BE3(i,:),q];
        end
        
        for i = 1:size(BE4,1)
            EB = BE4(i,1);
            PB = [JMV(EB,1),JMV(EB,4),JMV(EB,7)];
            TB = [temp_pre_inner(PB(1),1),temp_pre_inner(PB(2),1),temp_pre_inner(PB(3),1)];
            q = -1*alpha*(TB-Tf)/k;
            JBT24(i,:)=[BE4(i,:),q];
        end
        
        JBT2=[JBT21;JBT22;JBT23;JBT24];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%  迭代步骤D.4.3    %%%%%%%%%%
        %%%%%%%%%  初始化总体方程子块   %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M = zeros(noNodesV);
        CD = zeros(noNodesV);
        CDB = zeros(noNodesV, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤D.4.4       %%%%%%%
        %%%%%%   Me和CDe子块计算并组装    %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for e = 1 : noElement
            JMe = JMV(e,:);
            JXYe = JXYV(JMV(e,:),:);
            
            CDe = compute_CDe(JXYe, k, gaussPts, gaussWts);
            Me = compute_Me(JXYe, rho, Cv, gaussPts, gaussWts);
            
            CD(JMe, JMe) = CD(JMe, JMe) + CDe;
            M(JMe, JMe) = M(JMe, JMe) + Me;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤D.4.5    %%%%%%%%%
        %%%%%%%%     CDBe子块计算并组装  %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1 : size(JBT2, 1)
            
            JXYe = JXYV(JMV(JBT2(i,1),:),:);
            JBT2e = JBT2(i,:);
            
            CDBe = compute_CDBe(JXYe, JBT2e,  gaussPts, gaussWts);
            
            CDB(JMV(JBT2(i,1),:),1) = CDB(JMV(JBT2(i,1),:),1) + CDBe;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤D.4.6       %%%%%%%
        %%%%%%%%%     构建计算方程      %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha = 0.5;
        Knp1 = M + alpha * delta_t * CD;
        Kn = M - (1-alpha) * delta_t * CD;
        Fnnp1 = delta_t * CDB + Kn * temp_pre; % 对吗？
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%     迭代步骤D.4.7    %%%%%%%%
        %%%%%%%    求解方程，更新Tn+1k+1  %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_k1_inner = Knp1 \ Fnnp1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%    迭代步骤D.4.8    %%%%%%%%%
        %%%%%%%   对流换热迭代误差计算  %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iter_inner > maxIter_inner    % 牛顿流体，直接赋值ux_k=ux_k_1等
            break
        else
            
            absErrTemp = abs(temp_k1_inner - temp_pre_inner);
%             absErrTemp = norm(temp_k1_inner - temp_pre_inner);
            relErrTemp = max( abs(absErrTemp./temp_k1_inner) );
%             relErrTemp = absErrTemp ./ norm(temp_k1_inner);
            isConvergence = relErrTemp < convergenceCriterionREL_inner;
            
            iter_inner = iter_inner + 1;
            fprintf('   Inner iter.:  iter. = %d: && rel. err. temp_inner = %6.9f \n', iter_inner, relErrTemp)
            
            if isConvergence
                break
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%      迭代步骤D.5      %%%%%%%
    %%%%%%     Tn+1k+1复制给Tn+1     %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_k1 = temp_k1_inner;  % 当前时刻的收敛结果
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   迭代步骤D.6      %%%%%%%%%
    %%%%%%      累加非定常迭代次数     %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter = iter+1;
    %     t = iter*delta_t; %计算当前时刻
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%      迭代步骤D.7      %%%%%%%
    %%%%%      达到记录步长时记录数据     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(iter, detla_tag) == 0  % 每个记录步长存储一次结果
        tag=tag+1;
        Temp(:,tag) = temp_k1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%     迭代步骤D.8       %%%%%%%
    %%%%%        计算n+1时刻平均温度    %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INTT = 0;   % 初始化区域速度积分
    AREA = 0;   % 初始化区域面积
    for e = 1:noElement
        
        JXYe = JXYV(JMV(e,:),:);
        temp_k1e = temp_k1(JMV(e,:),:);
        
        [INTTe,AREAe] = compute_INTVe_AREAe(JXYe, temp_k1e, gaussPts, gaussWts);
        
        INTT = INTT + INTTe;
        AREA = AREA + AREAe;
        
    end
    
    avgT(iter+1) = INTT/AREA; % 区域平均温度
    iterStep(iter+1) = iter*delta_t; % 记录平均温度对应时间
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%      迭代步骤D.9       %%%%%%%
    %%%%%%%%%     计算相对误差     %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter > maxIter_out    % 牛顿流体，直接赋值ux_k=ux_k_1等
        break
    else
        
        absErrTemp = abs(temp_k1 - temp_pre);
%         absErrTemp = norm(temp_k1 - temp_pre);
        relErrTemp = max( abs(absErrTemp./temp_k1) );
%         relErrTemp = absErrTemp ./ norm(temp_k1);
        isConvergence = relErrTemp < convergenceCriterionREL_inner;
        
        iter_inner = iter_inner + 1;
        fprintf('Outer iter.:  iter. = %d: && rel. err. temp_outer = %6.9f \n', iter, relErrTemp)
        
        if isConvergence
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       迭代步骤D开始     %%%%%%%%%%%%
%%%%%%%%%%%%%     开始非定常迭代计算   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    绘制平均温度变化曲线    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(iterStep, avgT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    绘制平均温度变化曲线    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

