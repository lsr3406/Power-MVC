% encoding: utf-8
% @author: lsr3406
% @create on: 2018-01-19 11:03:57
% @update on: 2020-06-11 17:13:00

% matlab 初始化
clear all;		% 清除变量
close all;		% 关闭窗口

% 初始化常量, 配置文件等
config = conf();
addpath('./framework');
addpath('./framework/functions/');
addpath('./libs');

% 入口
Start.run(config);

%% 2018-01-19
% addpath() 方法传入的是引用的文件的绝对路径或在 matlab 当前的路径下的相对路径, 仅在 matlab 本次启动有效
% 配置文件中写的有控制器,方法名等参数, 本应接收其他参数, 开发阶段暂且将其写进去, 以后做成 api 时需要修改整个框架入口
% 添加的路径中, 如果文件夹里面有包(+folder), 则不需要手动引入, 引入父文件夹即可
