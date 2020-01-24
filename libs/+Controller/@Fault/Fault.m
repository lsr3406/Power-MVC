%% 电力系统暂态分析控制器
classdef Fault < handle
	properties
		% nothing
	end

	methods
		%% testF: 故障测试
		function testF(self)
		
			obj = 'lgq_eg8_1';
			ft = Model.Fault();
			ft.init(getMpcFault(obj));

			solver = [];
			res = ft.solveFault(solver, getMpcSteady(obj));

			viewModel = View.Plain();
			config.documentName = ['report_fault_', obj, '.txt'];
			viewModel.getFaultReport(ft, config, res);

			% debug
			save('test');
		end

		%% test39: 
		function test39(self)

			% 正常加载数据并计算
			obj = 'case39';
			ft = Model.Fault();
			ft.init(getMpcFault(obj));

			solver = [];
			res = ft.solveFault(solver, getMpcSteady(obj));

			% 画向量图
			self.drawNetStatus(ft.bus);  % 在复平面上画出各节点的电压
			self.strokeFaultStatus(ft);  % 在复平面上画出故障点的电压电流及其三序分量

			% 保存一个简单的故障计算报告
			% viewModel = View.Plain();
			% config.documentName = ['report_fault_', obj, '.txt'];
			% viewModel.getFaultReport(ft, config, res);

			% 由于函数运行结束后, 局部变量不能保存, 运行 save('xxx') 后可以保存成文件.
			% 在外面可以手动调用 load('xxx') 加载重要的变量
			save('test');
		end

		%% test9: 
		function test9(self)

			obj = 'case9';
			ft = Model.Fault();
			ft.init(getMpcFault(obj));

			solver = [];
			res = ft.solveFault(solver, getMpcSteady(obj));

			self.drawNetStatus(ft.bus);
			self.strokeFaultStatus(ft);

			% viewModel = View.Plain();
			% config.documentName = ['report_fault_', obj, '.txt'];
			% viewModel.getFaultReport(ft, config, res);

			% debug
			% save('test');
		end

		%% test9_zn: 9 节点系统变压器接地阻抗变化测试
		function test9_zn(self)

			% 变压器接地阻抗设置
			zn_param = [0.95, 1.03, 1.4, 2, 16];
			zn_list = (0.085.*zn_param - 0.1205).*3;

			obj = 'case9';
			mpcSteady = getMpcSteady(obj);
			ft = Model.Fault();
			ft.init(getMpcFault(obj));

			Ua = [];
			Ub = [];
			Uc = [];
			Ia = [];
			Ib = [];
			Ic = [];
			for k = 1:length(zn_param)
				ft.branch.xn = ones(length(ft.branch.xn), 1).*zn_list(k);
				ft.solveFault([], mpcSteady);

				Ua = [Ua; ft.itlog.Ufa];
				Ub = [Ub; ft.itlog.Ufb];
				Uc = [Uc; ft.itlog.Ufc];
				Ia = [Ia; ft.itlog.Ifa];
				Ib = [Ib; ft.itlog.Ifb];
				Ic = [Ic; ft.itlog.Ifc];
			end

			self.drawPhaseArrows(Ua, Ub, Uc, Ia, Ib, Ic);

		end

		%% drawPhaseArrows: 画多个相量
		function drawPhaseArrows(self, Fa, Fb, Fc, Ga, Gb, Gc)
			if nargin == 7
				maxF = (max(abs([Fa; Fb; Fc; 1])));
				maxG = (max(abs([Ga; Gb; Gc; 1])));
				figure();
				subplot(1, 2, 1);
					polar(0, maxF);
					hold on;
					self.strokePhaseArrow(Fa, Fb, Fc, 1, 1.5);
					title('故障电压相量图');
				subplot(1, 2, 2);
					polar(0, maxG);
					hold on;
					self.strokePhaseArrow(Ga, Gb, Gc, 1, 1.5);
					title('故障电流相量图');
			end
			if nargin == 4
				maxF = (max(abs([Fa; Fb; Fc; 1])));
				figure();
					polar(0, maxF);
					hold on;
					self.strokePhaseArrow(Fa, Fb, Fc, 1, 1.5);
					% title('电网参数相量图');
			end
		end

		%% drawNetStatus: 画电网的电压分布
		function drawNetStatus(self, bus)
			maxU = (max(abs([bus.Ua; bus.Ub; bus.Uc; 1])));
			figure();
				polar(0, maxU);
				hold on;
				self.strokePhaseArrow(bus.Ua, bus.Ub, bus.Uc, 1, 1.5);
				title('电网电压分布图');
		end

		%% strokeFaultStatus: 画故障点的向量图
		function strokeFaultStatus(self, ft)
			maxU = (max(abs([ft.itlog.Ufa; ft.itlog.Ufb; ft.itlog.Ufc; ft.itlog.Uf1; ft.itlog.Uf2; ft.itlog.Uf0; 1])));
			maxI = (max(abs([ft.itlog.Ifa; ft.itlog.Ifb; ft.itlog.Ifc; ft.itlog.If1; ft.itlog.If2; ft.itlog.If0; 1])));
			figure();
			subplot(1, 2, 1);
				polar(0, maxU);	% draw nothing
				hold on;
				self.strokeSequenceArrow(ft.itlog.Uf1, ft.itlog.Uf2, ft.itlog.Uf0, ft.T, 1.0, 1.5);
				self.strokePhaseArrow(ft.itlog.Ufa, ft.itlog.Ufb, ft.itlog.Ufc, 0.8, 2.5);
				title('故障电压相量图');
			subplot(1, 2, 2);
				polar(0, maxI);	% draw nothing
				hold on;
				self.strokeSequenceArrow(ft.itlog.If1, ft.itlog.If2, ft.itlog.If0, ft.T, 1.0, 1.5);
				self.strokePhaseArrow(ft.itlog.Ifa, ft.itlog.Ifb, ft.itlog.Ifc, 0.8, 2.5);
				title('故障电流相量图');
		end

		%% strokeSequenceArrow: 画正负零序相量
		function strokeSequenceArrow(self, F1, F2, F0, T, deep, lineWidth)
			if nargin <= 6
				lineWidth = 1.5;
			end
			if nargin <= 5
				deep = 1;
			end
			if nargin <= 4
				a = -0.5 + sqrt(3)./2.*1i;
				c = -0.5 - sqrt(3)./2.*1i;
				T = [1, 1, 1; c, a, 1; a, c, 1];
			end

			self.strokePhaseArrow(T(1, 1)*F1, T(2, 1)*F1, T(3, 1)*F1, deep, lineWidth);
			self.strokePhaseArrow(T(1, 2)*F2, T(2, 2)*F2, T(3, 2)*F2, deep, lineWidth);
			self.strokePhaseArrow(T(1, 3)*F0, T(2, 3)*F0, T(3, 3)*F0, deep, lineWidth);
			
		end

		%% strokePhaseArrow: 画三相相量
		function strokePhaseArrow(self, Fa, Fb, Fc, deep, lineWidth)
			if nargin <= 5
				lineWidth = 0;
			end
			if nargin <= 4
				deep = 1;
			end
			amber = [255, 193, 7]./256;
			red = [244, 67, 54]./256;
			green = [76, 175, 80]./256;

			if lineWidth == 0
				% plot(Fa, '.', 'lineWidth', 20, 'Color', amber.*deep);
				% plot(Fb, '.', 'lineWidth', 20, 'Color', green.*deep);
				% plot(Fc, '.', 'lineWidth', 20, 'Color', red.*deep);
				polar(angle(Fa), abs(Fa), 'oy');
				polar(angle(Fb), abs(Fb), 'og');
				polar(angle(Fc), abs(Fc), 'or');
			else
				quiver(zeros(length(Fa), 1), zeros(length(Fa), 1), real(Fa), imag(Fa), 'LineWidth', lineWidth, 'Color', amber.*deep);
				quiver(zeros(length(Fb), 1), zeros(length(Fb), 1), real(Fb), imag(Fb), 'LineWidth', lineWidth, 'Color', green.*deep);
				quiver(zeros(length(Fc), 1), zeros(length(Fc), 1), real(Fc), imag(Fc), 'LineWidth', lineWidth, 'Color', red.*deep);
			end
		end

	end
end
