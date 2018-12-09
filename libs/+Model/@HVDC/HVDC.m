% encoding: utf-8
% @author: siru
% @create on: 2018-09-29 14:13:16
% @update on: 2018-12-09 14:36:52
%% 直流输电
classdef HVDC < handle

	properties (Constant)
		alphaMin = 10/180*pi;	% 10 deg
		alphaMax = 20/180*pi;	% 20 deg
		gammaMin = 10/180*pi;	% 10 deg
		currentMargin = 0.1	% 电流裕度
	end % properties (Constant)

	properties
		f;		% 整流站节点索引
		t;		% 逆变站节点索引
		source;	% 整流站节点 id
		sink;	% 逆变站节点 id
		Er;		% 整流站换流变前交流电压
		Ei;		% 逆变站换流变前交流电压
		Pr;		% 整流站直流侧有功
		Pi;		% 逆变站直流侧有功
		Qr;		% 整流站直流侧无功
		Qi;		% 逆变站直流侧无功
		Prac;		% 整流站交流侧有功
		Piac;		% 逆变站交流侧有功
		Qrac;		% 整流站交流侧无功
		Qiac;		% 逆变站交流侧无功
		phir;	% 整流站功率因数角
		phii;	% 逆变站功率因数角
		Vdr0;	% 整流站电压(无换向)
		Vdi0;	% 逆变站电压(无换向)
		Vdr;	% 整流站电压
		Vdi;	% 逆变站电压
		Tr;		% 整流站换流变变比
		Ti;		% 逆变站换流变变比
		Hr;		% 整流站换流变档位辅助变量, g = 3 - h
		Hi;		% 逆变站换流变档位辅助变量, g = 3 - h
		Rcr;	% 整流站等值换向电阻
		Rci;	% 逆变站等值换向电阻
		Rl;		% 直流线路电阻
		Id;		% 直流线路电流
		alpha;	% 触发角
		beta;	% 导通角
		gamma;	% 熄弧角
		operatingMode;	% 运行方式
		sourceParam;	% 逆变站参数
		sinkParam;		% 整流站参数
	end % properties

	properties (Dependent)
		m1;	% 三种运行模式对应的索引
		m2;
		m3;
		gear;	% 整流站变压器档位 (1-5)
	end % properties (Dependent)

	methods

		%% m1: 计算并缓存模式 1 对应的节点索引
		function [m1] = get.m1(self)
			m1 = find(self.operatingMode == 1);
		end
		
		%% m2: 计算并缓存模式 2 对应的节点索引
		function [m2] = get.m2(self)
			m2 = find(self.operatingMode == 2);
		end
		
		%% m3: 计算并缓存模式 3 对应的节点索引
		function [m3] = get.m3(self)
			m3 = find(self.operatingMode == 3);
		end

		%% get.gear: 返回整流站变压器的档位
		function [gear] = get.gear(self)
			gear = 3 - self.Hr;
		end
		
		%% 建立直流系统模型
		function init(self, mpc)
			l = length(mpc.data(:, 1));
			self.f = [];
			self.t = [];
			self.source = mpc.data(:, 1);
			self.sink = mpc.data(:, 2);
			self.Tr = mpc.data(:, 3);
			self.Ti = mpc.data(:, 4);
			self.Rcr = mpc.data(:, 5);
			self.Rci = mpc.data(:, 6);
			self.Rl = mpc.data(:, 7);
			self.operatingMode = mpc.data(:, 8);
			self.sourceParam = mpc.data(:, 9);
			self.sinkParam = mpc.data(:, 10);

			self.Er = zeros(l, 1);
			self.Ei = zeros(l, 1);
			self.Pr = zeros(l, 1);
			self.Pi = zeros(l, 1);
			self.Qr = zeros(l, 1);
			self.Qi = zeros(l, 1);
			self.Prac = zeros(l, 1);
			self.Piac = zeros(l, 1);
			self.Qrac = zeros(l, 1);
			self.Qiac = zeros(l, 1);
			self.phir = zeros(l, 1);
			self.phii = zeros(l, 1);
			self.Vdr0 = zeros(l, 1);
			self.Vdi0 = zeros(l, 1);
			self.Vdr = zeros(l, 1);
			self.Vdi = zeros(l, 1);
			self.Id = zeros(l, 1);
			self.Hr = zeros(l, 1);
			self.Hi = zeros(l, 1);
			self.alpha = ones(l, 1) .* pi ./ 12;
			self.beta = ones(l, 1) .* pi ./ 12;
			self.gamma = ones(l, 1) .* pi ./ 12;

		end

		%% setBusIndex: 设置整流站与逆变站在交流系统中的节点的索引
		function setBusIndex(self, bus)
			self.f = getIndex(bus.id, self.source);
			self.t = getIndex(bus.id, self.sink);
		end

		%% 计算过程: 先给定整流站与逆变站的交流电压, 解出直流系统的方程后, 得到整流站与逆变站相应的 P, Q, 此时完成一次更新
		
		%% render: public
		function [converge] = render(self, Er, Ei, epsilon)
			if nargin == 3
				epsilon = 1e-6;
			end

			% 检查直流系统是否收敛, 并更新交流侧电压
			if norm(self.Er - Er, inf) <= epsilon && norm(self.Ei - Ei, inf) <= epsilon
				converge = true;
				return;
			end
			
			converge = false;
			self.Er = Er;
			self.Ei = Ei;

			% 更新一次
			self.update();

		end

		%% 根据已知的换流变外侧交流电压 E 计算直流系统两侧相应的交流有功, 无功
		function update(self)
			
			%% 计算换流变内侧的电压
			self.Vdr0 = 1.35 .* self.Tr .* self.Er;
			self.Vdi0 = 1.35 .* self.Ti .* self.Ei;

			% m1 = self.m1;
			self.Pr = self.sourceParam ./ 100;	% TODO baseMVA
			self.gamma = self.sinkParam;
			self.Id = calId(self.Rl - self.Rci, self.Vdi0.*cos(self.gamma), -self.Pr);
			self.Vdi = self.Vdi0 .* cos(self.gamma) - self.Rci .* self.Id;
			self.Vdr = self.Vdi + self.Rl .* self.Id;
			self.Pi = self.Vdi .* self.Id;
			% self.Pr = self.Vdr .* self.Id;
			self.alpha = acos((self.Vdr + self.Rcr .* self.Id) ./ self.Vdr0);

			% 检查 alpha 越限情况, 并及时处理
			oli = find(real(self.alpha) < self.alphaMin);
			oli = intersect(oli, find(self.Hr < 2));
			self.Hr(oli) = min(ceil(40 .* (self.Tr(oli).*cos(self.alpha(oli)) ./ cos(self.alphaMin) - 1)), 2);
			% assert(~any(self.Hr(oli) > 2));
			self.Tr(oli) = 1 + self.Hr(oli) .* 0.025;
			self.Vdr0(oli) = 1.35 .* self.Tr(oli) .* self.Er(oli);
			self.alpha(oli) = acos((self.Vdr(oli) + self.Rcr(oli) .* self.Id(oli)) ./ self.Vdr0(oli));
			
			oli = find(real(self.alpha) > self.alphaMax);
			oli = intersect(oli, find(self.Hr > -2));
			self.Hr(oli) = max(floor(40 .* (self.Tr(oli).*cos(self.alpha(oli)) ./ cos(self.alphaMax) - 1)), -2);
			% assert(~any(self.Hr(oli) > 2));
			self.Tr(oli) = 1 + self.Hr(oli) .* 0.025;
			self.Vdr0(oli) = 1.35 .* self.Tr(oli) .* self.Er(oli);
			self.alpha(oli) = acos((self.Vdr(oli) + self.Rcr(oli) .* self.Id(oli)) ./ self.Vdr0(oli));

			% assert(all(self.alpha < 0.09));

			self.beta = acos((self.Vdi - self.Rci .* self.Id) ./ self.Vdi0);
			self.phir = acos(self.Vdr ./ self.Vdr0);
			self.phii = acos(self.Vdi ./ self.Vdi0);
			self.Qr = self.Pr .* tan(self.phir);
			self.Qi = self.Pi .* tan(self.phii);

			betaTemp = acos((self.Vdr - self.Rcr .* self.Id) ./ self.Vdr0);
			self.Prac = self.Er.^2 .* self.Tr.^2 ./ (4 .* self.Rcr) .* (cos(2.*self.alpha) - cos(2 .* betaTemp));
			self.Piac = self.Ei.^2 .* self.Ti.^2 ./ (4 .* self.Rci) .* (cos(2.*self.gamma) - cos(2 .* self.beta));
			self.Qrac = self.Er.^2 .* self.Tr.^2 ./ (4 .* self.Rcr) .* (2 .* (betaTemp - self.alpha) + sin(2.*self.alpha) - sin(2 .* betaTemp));
			self.Qiac = self.Ei.^2 .* self.Ti.^2 ./ (4 .* self.Rci) .* (2 .* (self.beta - self.gamma) + sin(2.*self.gamma) - sin(2 .* self.beta));

			%% calId: 解二次方程计算直流电流
			function [res] = calId(a, b, c)
				delta = b .^ 2 - 4 .* a .* c;
				assert(all(delta > 0));
				res = (-b + sqrt(delta)) ./ (2.*a);
			end

		end

		% %% update: 根据已知的换流变外侧交流电压 E 计算直流系统两侧相应的交流有功, 无功
		% function update(self)
		% 	%% 计算换流变内侧的电压
		% 	self.Vdr0 = 1.35 .* self.Tr .* self.Er;
		% 	self.Vdi0 = 1.35 .* self.Ti .* self.Ei;

		% 	%% 根据运行方式计算换流器参数

		% 	while true

		% 		% 模式 1, 先计算逆变站直流电压, 再计算整流站直流电压, 最后计算触发角
		% 		% @self.sourceParam: 直流电流 Id
		% 		% @self.sinkParam: 熄弧角 gamma
		% 		m1 = self.m1;
		% 		self.Id(m1) = self.sourceParam(m1);
		% 		self.gamma(m1) = self.sinkParam(m1);
		% 		self.Vdi(m1) = self.Vdi0(m1) .* cos(self.gamma(m1)) - self.Rci(m1) .* self.Id(m1);
		% 		self.Vdr(m1) = self.Vdi(m1) + self.Rl(m1) .* self.Id(m1);
		% 		self.alpha(m1) = acos((self.Vdr(m1) + self.Rcr(m1) .* self.Id(m1)) ./ self.Vdr0(m1));
		% 		self.beta(m1) = acos((self.Vdi(m1) - self.Rci(m1) .* self.Id(m1)) ./ self.Vdi0(m1));

		% 		% 模式 2, 先计算整流站直流电压, 再计算逆变站直流电压, 最后计算熄弧角
		% 		% @self.sourceParam: 触发角 alpha
		% 		% @self.sinkParam: 直流电流 Id
		% 		m2 = self.m2;
		% 		% self.alpha(m2) = self.sourceParam(m2);
		% 		% self.Id(m2) = self.sinkParam(m2);
		% 		self.Vdr(m2) = self.Vdr0(m2) .* cos(self.alpha(m2)) - self.Rcr(m2) .* self.Id(m2);
		% 		self.Vdi(m2) = self.Vdr(m2) - self.Rl(m2) .* self.Id(m2);
		% 		self.gamma(m2) = acos((self.Vdi(m2) + self.Rci(m2) .* self.Id(m2)) ./ self.Vdi0(m2));
		% 		self.beta(m2) = acos((self.Vdi(m2) - self.Rci(m2) .* self.Id(m2)) ./ self.Vdi0(m2));

		% 		% 模式 3, 整理系数, 先计算直流电流, 再计算整流侧与逆变侧电压
		% 		% @self.sourceParam: 触发角 alpha
		% 		% @self.sinkParam: 导通角 beta
		% 		m3 = self.m3;
		% 		% self.alpha(m3) = self.sourceParam(m3);
		% 		% self.beta(m3) = self.sinkParam(m3);
		% 		self.Id(m3) = (self.Vdr0(m3) .* cos(self.alpha(m3)) - self.Vdi0(m3) .* cos(self.beta(m3))) ./ (self.Rl(m3) + self.Rci(m3) + self.Rcr(m3));
		% 		self.Vdr(m3) = self.Vdr0(m3) .* cos(self.alpha(m3)) - self.Rcr(m3) .* self.Id(m3);
		% 		self.Vdi(m3) = self.Vdi0(m3) .* cos(self.beta(m3)) + self.Rci(m3) .* self.Id(m3);
		% 		self.gamma(m3) = acos((self.Vdi(m3) + self.Rci(m3) .* self.Id(m3)) ./ self.Vdi0(m3));

		% 		if ~self.modeConvert()
		% 			break;
		% 		end

		% 	end % while

		% 	%% 计算功率
		% 	self.phir = acos(self.Vdr ./ self.Vdr0);
		% 	self.phii = acos(self.Vdi ./ self.Vdi0);
		% 	self.Pr = self.Vdr .* self.Id;
		% 	self.Pi = self.Vdi .* self.Id;
		% 	self.Qr = self.Pr .* tan(self.phir);
		% 	self.Qi = self.Pi .* tan(self.phii);

		% end

		

		% %% modeConvert: function description
		% function flag = modeConvert(self)
			
		% 	% 1 转 2, 判断依据, 触发角越限
		% 	newMode2 = find(self.alpha < self.alphaMin - eps);
		% 	self.alpha(newMode2) = self.alphaMin;
		% 	self.Id(newMode2) = self.Id(newMode2) .* (1 - self.currentMargin);
		% 	self.operatingMode(newMode2) = 2;

		% 	% 检查裕量不足的逆变站
		% 	vrTemp = self.Vdr0(newMode2) .* cos(self.alphaMin) - self.Rcr(newMode2) .* self.Id(newMode2);
		% 	viTemp = vrTemp - self.Id(newMode2) .* self.Rl(newMode2);
		% 	self.gamma(newMode2) = acos((viTemp + self.Rci(newMode2) .* self.Id(newMode2)) ./ self.Vdi0(newMode2));

		% 	% 2 转 3, 判断依据, 熄弧角越限
		% 	newMode3 = find(self.gamma < self.gammaMin - eps);
		% 	self.alpha(newMode3) = self.alphaMin;
		% 	self.beta(newMode3) = acos((self.Vdi0(newMode3) .* cos(self.gammaMin) - 2 .* self.Rci(newMode3) .* (self.Id(newMode3) ./ (1 - self.currentMargin))) ./ self.Vdi0(newMode3));
		% 	self.gamma(newMode3) = acos((self.Vdi0(newMode3) .* cos(self.beta(newMode3)) + 2 .* self.Rci(newMode3) .* self.Id(newMode3)) ./ self.Vdi0(newMode3));
		% 	self.operatingMode(newMode3) = 3;

		% 	% 2, 3 转 1, 判断依据, 电流恢复
		% 	% TODO 这里用了 sourceParam, 这意味着所有的直流系统必须默认处于运行方式一
		% 	newMode1 = find(self.Id > self.sourceParam + eps);
		% 	self.Id(newMode1) = self.sourceParam(newMode1);
		% 	self.gamma(newMode1) = self.sinkParam(newMode1);
		% 	self.operatingMode(newMode1) = 1;

		% 	flag = any([newMode1; newMode2; newMode3]);
		% end



		%% toString: 
		function [res] = toString(self)
			res = '直流输电系统实验数据: \n';
			
			switch self.operatingMode
				case 1
					omStr = num2str([self.Pr, self.gamma * 180 / pi], '整流站 Pr = %6.2f, 逆变站 γ = %6.2f°');
				% case 2
				% 	omStr = num2str([self.alpha * 180 / pi, self.Id], '整流站 α = %6.2f°, 逆变站 Id = %6.2f');
				% case 3
				% 	omStr = num2str([self.alpha * 180 / pi, self.beta * 180 / pi], '整流站 α = %6.2f°, 逆变站 β = %6.2f°');
				otherwise
					error('Illegal operating mode');
			end
			res = [res, num2str(self.operatingMode', '运行方式: %d'), '    ', omStr, '\n'];


			res = [res, num2str(self.alpha * 180 / pi, '整流站触发角 α: %4.2f°'), blanks(4)];
			res = [res, num2str(self.beta * 180 / pi, '逆变站导通角 β: %4.2f°'), blanks(4)];
			res = [res, num2str(self.gamma * 180 / pi, '逆变站熄弧角 γ: %4.2f°'), '\n'];

			res = [res, num2str([self.Tr, self.Ti], '换流变变比 - 整流站 Tr / 逆变站 Ti: %4.3f / %4.3f'), '\n'];
			res = [res, num2str([self.Rcr, self.Rci], '换向电阻 - 整流站 Rr / 逆变站 Ri: %4.2f / %4.2f'), '\n'];
			res = [res, num2str(self.Rl, '直流电阻: %4.2f'), '\n'];

			res = [res, num2str([self.Er, self.Er .* self.Tr], '整流站换流变交流线电压 Er - 网侧 / 阀侧: %4.2f / %4.2f'), '\n'];
			res = [res, num2str([self.Ei, self.Ei .* self.Ti], '逆变站换流变交流线电压 Ei - 网侧 / 阀侧: %4.2f / %4.2f'), '\n'];

			res = [res, num2str(self.Id, '直流电流: %4.2f'), '\n'];
			res = [res, num2str([self.Vdr0, self.Vdr], '整流站直流电压 - Vdr0 / Vdr: %4.2f / %4.2f'), '\n'];
			res = [res, num2str([self.Vdi0, self.Vdi], '逆变站直流电压 - Vdi0 / Vdi: %4.2f / %4.2f'), '\n'];

			res = [res, num2str([self.Pr*1e2, self.Qr*1e2, cos(self.phir)], '整流站功率 P+jQ: %4.2fMW + j%4.2fMVar, cosφ = %4.2f'), '\n'];
			res = [res, num2str([self.Pi*1e2, self.Qi*1e2, cos(self.phii)], '逆变站功率 P+jQ: %4.2fMW + j%4.2fMVar, cosφ = %4.2f'), '\n'];

			% switch self.operatingMode
			% 	case 1
			% 		omStr = num2str([self.Pr, self.gamma * 180 / pi], '整流站 Pr = %6.2f, 逆变站 γ = %6.2f°');
			% 	% case 2
			% 	% 	omStr = num2str([self.alpha * 180 / pi, self.Id], '整流站 α = %6.2f°, 逆变站 Id = %6.2f');
			% 	% case 3
			% 	% 	omStr = num2str([self.alpha * 180 / pi, self.beta * 180 / pi], '整流站 α = %6.2f°, 逆变站 β = %6.2f°');
			% 	otherwise
			% 		error('Illegal operating mode');
			% end
			% res = [res, num2str(self.operatingMode', '运行方式: %d'), '    ', omStr, '\n'];


			% res = [res, num2str(self.alpha * 180 / pi, '整流站触发角 α: %4.2f°'), blanks(4)];
			% res = [res, num2str(self.beta * 180 / pi, '逆变站导通角 β: %4.2f°'), blanks(4)];
			% res = [res, num2str(self.gamma * 180 / pi, '逆变站熄弧角 γ: %4.2f°'), '\n'];

			% res = [res, num2str([self.Tr, self.Ti], '换流变变比 - 整流站 Tr / 逆变站 Ti: %4.3f / %4.3f'), '\n'];
			% res = [res, num2str([self.Rcr, self.Rci], '换向电阻 - 整流站 Rr / 逆变站 Ri: %4.2fΩ / %4.2fΩ'), '\n'];
			% res = [res, num2str(self.Rl, '直流电阻: %4.2fΩ'), '\n'];

			% res = [res, num2str([self.Er/1000, self.Er/1000 .* self.Tr], '整流站换流变交流线电压 Er - 网侧 / 阀侧: %4.2fkV / %4.2fkV'), '\n'];
			% res = [res, num2str([self.Ei/1000, self.Ei/1000 .* self.Ti], '逆变站换流变交流线电压 Ei - 网侧 / 阀侧: %4.2fkV / %4.2fkV'), '\n'];
			
			% res = [res, num2str(self.Id / 1000, '直流电流: %4.2fkA'), '\n'];
			% res = [res, num2str([self.Vdr0/1000, self.Vdr/1000], '整流站直流电压 - Vdr0 / Vdr: %4.2fkV / %4.2fkV'), '\n'];
			% res = [res, num2str([self.Vdi0/1000, self.Vdi/1000], '逆变站直流电压 - Vdi0 / Vdi: %4.2fkV / %4.2fkV'), '\n'];

			% res = [res, num2str([self.Pr*1e-6, self.Qr*1e-6, cos(self.phir)], '整流站功率 P+jQ: %4.2fMW + j%4.2fMVar, cosφ = %4.2f'), '\n'];
			% res = [res, num2str([self.Pi*1e-6, self.Qi*1e-6, cos(self.phii)], '逆变站功率 P+jQ: %4.2fMW + j%4.2fMVar, cosφ = %4.2f'), '\n'];
		end

	end % methods

end % classdef

