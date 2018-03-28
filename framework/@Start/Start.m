%% Start 程序入口
classdef Start < handle

	properties (SetAccess = private)
		controller;
		method;
		config;
	end

	methods

		%% Start: 入口构造器
		function [self] = Start()
			addpath('./framework/functions');
			% 这里暂且将控制器和方法全部引入
			addpath('./libs');
		end

		%% controllerInit: 初始化控制器
		function controllerInit(self, controller)
			self.controller = controller;
		end

		%% methodInit: 初始化控制器
		function methodInit(self, method)
			self.method = method;
		end

		%% run: 入口方法
		function run(self, config)

			self.config = config;
			self.controller = config.controller;
			self.method = config.method;

			eval(['obj = Controller.', self.controller, '();']);
			eval(['obj.', self.method, '();']);
		end

	end
end