function PIE = compute_PIE(JXYe, velocityE, penaltyCoef)

%%%%%%% 各个结点对应局部坐标
xi = [-1,0,1,-1,0,1,-1,0,1];
et = [-1,-1,-1,0,0,0,1,1,1];
%%%%%%% 各个结点对应局部坐标

%%%%%%% 初始化PIE
PIE = zeros(9,1);
%%%%%%% 初始化PIE

%%%%%%% 循环计算单元内各个结点压力
for i = 1:9
    
    %%%%%%% 速度插值函数对kesi和eta的导数
    [~, fyV_xi, fyV_et] = shapeFunction(xi(i), et(i) , 'quadratic');
    %%%%%%% 速度插值函数对kesi和ita的导数
    
    %%%%%%%% Jacobi相关计算
        Jacobi = [fyV_xi, fyV_et]' * JXYe;
        fy_x = Jacobi \ [fyV_xi, fyV_et]';
    %%%%%%%% Jacobi相关计算
    
    %%%%%%% 单元内结点压力计算
    PIE(i,1) = penaltyCoef * (fy_x(1,:) * velocityE(:,1) + fy_x(2,:) * velocityE(:,2));
    %%%%%%% 单元内结点压力计算
end
%%%%%%% 循环计算单元内各个结点压力
end