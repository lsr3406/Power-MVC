# Power-MVC
电力系统分析常用工具

### 1. 安装

将下面两个路径(本工程下)添加至 MATLAB 的环境变量中并重启即可, 
可参考 matpower 的安装方式

> ./libs  
> ./framework/functions

### 2. 如何使用

##### 普通潮流:

目前的启动方式有(0-1 启动与直流潮流启动), 算法有牛顿-拉夫逊法与 PQ 分解法(FD, FDBX, FDXB)

``` matlab
% 创建稳态分析实例
ss = Model.SteadyState();
ss.init(getMpcSteady('bus4'));

% 设置求解器的基本信息
solver.method = 'NR';   % NR, FD, FDBX, FDXB
solver.n_iters_max = 10;
solver.epsilon = 1e-5;
solver.start = 'flat';  % flat, dc

% 计算潮流, 完成后可直接查看 ss 变量的字段
res = ss.solvePowerFlow(solver);
```

##### 最优潮流：

目前仅支持内点法, 求解思路参考王锡凡《现代电力系统分析》第三章的算例

``` matlab
% 创建稳态分析实例
ss = Model.SteadyState();
ss.init('case5_test');

% 设置求解器的基本信息
solver.n_iters_max = 50;  % 最大迭代次数
solver.epsilon = 1e-6;     % 对偶间隙需满足的精度
solver.sigma = 0.1;        % 向心参数

% 求解 OPF, 完成后可查看 ss.opf 中的相关字段
ss.solveOptimalPowerFlow(solver);
```

##### 两端直流输电系统独立计算

换流变压器变比自动调节允许 `1:0.95`, `1:0.975`, `1:1`, `1:1.025`, `1:1.05` 共 5 档
目前稳态分析中仅支持整流侧定电流, 逆变侧定熄弧角的运行方式, 暂态分析中支持以下 3 种运行方式:

    整流侧定电流, 逆变侧定熄弧角
    整流侧定最小触发角, 逆变侧定电流
    过渡运行方式

``` matlab
% 创建一般的直流输电系统实例并初始化
hvdc = Model.HVDC();
hvdc.init(dcm_com());

% 两个换流站分别加
hvdc.render(1.032, 1.061);
fprintf(hvdc.toString());

% 运行完成后可调用 hvdc.toString() 并打印获取结果
```

##### 交直流潮流

目前仅支持交替求解, 有 bug

``` matlab
% 创建稳态分析实例
ss = Model.SteadyState();
ss.init('case14_int');

% 创建一般的直流输电系统实例并初始化
hvdc = Model.HVDC();
hvdc.init(dcm_com());

% 向交流系统中挂载直流系统
ss.hvdcInit(hvdc)

% 设置求解器的基本信息
solver.method = 'NR';
solver.n_iters_max = 20;
solver.epsilon = 1e-6;

% 计算潮流, 完成后可直接查看 ss 变量的字段, 直流输电系统的参数保存在 ss.hvdc 中
ss.solvePowerFlow(solver);
```

