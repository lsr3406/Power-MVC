% encoding: utf-8
% @author: lsr3406
% @create on: 2018-01-19 11:03:57
% @update on: 2018-03-21 14:36:10

% matlab 初始化
clear all;		% 清除变量
close all;		% 关闭窗口

% 初始化常量, 配置文件等
config = conf();
addpath('./framework');

% 入口
Start().run(config);

%% 程序到此结束, 就是这么简单

%% debug
load('test.mat');

%% 2018-01-20
% 入口方法中使用 eval() 是迫不得已, 做成 api 后需要注意安全性问题
%% 2018-01-19
% addpath() 方法传入的是引用的文件的绝对路径或在 matlab 当前的路径下的相对路径, 仅在 matlab 本次启动有效
% main 脚本中调用入口时使用 Start 的静态方法更有逼格, 但由于 matlab 的面向对象机制不够成熟, 难以写出高效的启动方式. 因此这里在实例化 start 对象后才调用的入口方法.
% 配置文件中写的有控制器,方法名等参数, 本应接收其他参数, 开发阶段暂且将其写进去, 以后做成 api 时需要修改整个框架入口
% 添加的路径中, 如果文件夹里面有包(+folder), 则不需要手动引入, 引入父文件夹即可

