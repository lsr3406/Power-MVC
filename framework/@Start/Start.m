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
			% do nothing
		end

		%% controllerInit: 初始化控制器
		function controllerInit(self, controller)
			self.controller = controller;
		end

		%% methodInit: 初始化方法
		function methodInit(self, method)
			self.method = method;
		end
	end

	methods (Static)
		%% run: 入口方法
		function run(config)

			addpath('./framework/functions');
			addpath('./libs');

			self.config = config;
			self.controller = config.controller;
			self.method = config.method;

			eval(['obj = Controller.', self.controller, '();']);
			eval(['obj.', self.method, '();']);
		end
	end
end