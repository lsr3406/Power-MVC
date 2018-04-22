%% 纯文本文件输出
classdef Plain < handle
	properties
		% nothing
	end
	methods

		%% getPowerFlowReport: 返回电力系统潮流计算结果 (TODO string)
		function [result] = getPowerFlowReport(self, ss, config, result)
			
			file_t = fopen(config.documentName,'w');

			fprintf(file_t, '%s\n', '=========================================================================');
			fprintf(file_t, '%s\n', '                          电力系统潮流计算报告');
			fprintf(file_t, '%s\n', '=========================================================================');

			fprintf(file_t, '\n%s%d%s%d%s%d%s%d\n','节点数       PQ / PV / 平衡 / 合计      ',ss.pqCount,' / ',ss.pvCount,' / ',ss.refCount,' / ',length(ss.nodes.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','支路数       普通线路 / 变压器 / 合计      ',ss.lineCount,' / ',ss.transformerCount,' / ',length(ss.branches.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','发电机数       调压厂 / 调频厂 / 合计      ',ss.pvGens,' / ',ss.refGens,' / ',ss.pvGens+ss.refGens);
			fprintf(file_t, '%s%d\n','负荷数                                  ',ss.loadCount);

			fprintf(file_t, '\n%s%s%d\n', '计算方法: ', config.method);

			if result.status == 101
				fprintf(file_t, '\n%s\n', '迭代次数超过最大值,计算失败');
				fclose(file_t);
				return;
			elseif result.status == 102
				fprintf(file_t, '\n%s\n', '迭代终止');
				fclose(file_t);
				return;
			elseif result.status ~= 1
				fprintf(file_t, '\n%s\n', '未知错误: ', result.status);
				fclose(file_t);
				return;
			end

			fprintf(file_t, '\n%s%d\n', '迭代次数: ',result.it);

			fprintf(file_t, '\n%s\n', '节点电压,功率');
			fprintf(file_t, '\n%s\n', ' id     电压      相角   发电机有功  发电机无功  负荷有功    负荷无功   补偿有功  补偿无功');
			for k = 1:length(ss.nodes.id)
				fprintf(file_t, '  %d',ss.nodes.id(k));
				fprintf(file_t, '   %6.3f',ss.nodes.mag(k),ss.nodes.ang(k).*180./pi,ss.nodes.Pg(k).*100,ss.nodes.Qg(k).*100,ss.nodes.Pd(k).*100,ss.nodes.Qd(k).*100,ss.nodes.Pc(k).*100,ss.nodes.Qc(k).*100);
				fprintf(file_t, '\n');
			end
			fprintf(file_t, '%s%6.3f  %6.3f   %6.3f  %6.3f   %6.3f  %6.3f','合计:                     ',sum(ss.nodes.Pg).*100,sum(ss.nodes.Qg).*100,sum(ss.nodes.Pd).*100,sum(ss.nodes.Qd).*100,sum(ss.nodes.Pc.*100),sum(ss.nodes.Qc.*100));
			
			fprintf(file_t, '\n');

			fprintf(file_t, '\n%s\n\n', '线路功率');
			fprintf(file_t, '%s\n', 'id from to    Pij      Qij        Pji      Qji        dP      dQ');
			for k = ss.branches.id'
				fprintf(file_t, '  %d',[ss.branches.id(k),ss.branches.fid(k),ss.branches.tid(k)]);
				fprintf(file_t, '   %6.3f',[ss.branches.Pij(k),ss.branches.Qij(k),ss.branches.Pji(k),ss.branches.Qji(k),ss.branches.dP(k),ss.branches.dQ(k)]*100);
				fprintf(file_t, '\n');
			end
			
			fprintf(file_t, '总线损');
			fprintf(file_t, '   %6.3f',sum(ss.branches.dP)*100,sum(ss.branches.dQ)*100);

			fprintf(file_t, '\n');
			fprintf(file_t, '\n');

			fclose(file_t);

		end

		%% getFaultReport: 打印故障计算信息
		function getFaultReport(self, ft, config, result)

			file_t = fopen(config.documentName, 'w');

			fprintf(file_t, '%s\n', '=========================================================================');
			fprintf(file_t, '%s\n', '                          电力系统故障计算报告');
			fprintf(file_t, '%s\n', '=========================================================================');

			fprintf(file_t, '\n%s%d%s%d%s%d%s%d\n','节点数       PQ / PV / 平衡 / 合计      ',ft.ss.pqCount,' / ',ft.ss.pvCount,' / ',ft.ss.refCount,' / ',length(ft.ss.nodes.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','支路数       普通线路 / 变压器 / 合计      ',ft.ss.lineCount,' / ',ft.ss.transformerCount,' / ',length(ft.ss.branches.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','发电机数       调压厂 / 调频厂 / 合计      ',ft.ss.pvGens,' / ',ft.ss.refGens,' / ',ft.ss.pvGens+ft.ss.refGens);
			fprintf(file_t, '%s%d\n','负荷数                                  ',ft.ss.loadCount);

			if(result.status ~= 1)
				return;
			end

			fprintf(file_t, '%s%d\n', '故障节点: ', ft.fault.nid);
			fprintf(file_t, '%s%6.3f%+6.3f%s\n', '过渡阻抗 zf: ', real(ft.fault.zf), imag(ft.fault.zf), 'i');
			fprintf(file_t, '%s%6.3f%+6.3f%s\n', '过渡阻抗 zg: ', real(ft.fault.zg), imag(ft.fault.zg), 'i');
			fprintf(file_t, '%s\n', ['故障类型: ', ft.fault.type]);
			fprintf(file_t, '%s\n', ['故障相: ', ft.fault.phase]);
			fprintf(file_t, '\n');

			fprintf(file_t, '%s  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', '故障电流(序): ', abs(ft.itlog.If1), '∠', angle(ft.itlog.If1).*180./pi, abs(ft.itlog.If2), '∠', angle(ft.itlog.If2).*180./pi, abs(ft.itlog.If0), '∠', angle(ft.itlog.If0).*180./pi);
			fprintf(file_t, '%s  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', '故障电流(相): ', abs(ft.itlog.Ifa), '∠', angle(ft.itlog.Ifa).*180./pi, abs(ft.itlog.Ifb), '∠', angle(ft.itlog.Ifb).*180./pi, abs(ft.itlog.Ifc), '∠', angle(ft.itlog.Ifc).*180./pi);
			fprintf(file_t, '%s  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', '故障电压(序): ', abs(ft.itlog.Uf1), '∠', angle(ft.itlog.Uf1).*180./pi, abs(ft.itlog.Uf2), '∠', angle(ft.itlog.Uf2).*180./pi, abs(ft.itlog.Uf0), '∠', angle(ft.itlog.Uf0).*180./pi);
			fprintf(file_t, '%s  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', '故障电压(相): ', abs(ft.itlog.Ufa), '∠', angle(ft.itlog.Ufa).*180./pi, abs(ft.itlog.Ufb), '∠', angle(ft.itlog.Ufb).*180./pi, abs(ft.itlog.Ufc), '∠', angle(ft.itlog.Ufc).*180./pi);
			fprintf(file_t, '\n');

			for k = 1:length(ft.ss.nodes.id)
				fprintf(file_t, '%2d  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', ft.ss.nodes.id(k), abs(ft.nodes.Ua(k)), '∠', angle(ft.nodes.Ua(k)).*180./pi, abs(ft.nodes.Ub(k)), '∠', angle(ft.nodes.Ub(k)).*180./pi, abs(ft.nodes.Uc(k)), '∠', angle(ft.nodes.Uc(k)).*180./pi);
			end
			fprintf(file_t, '\n');

			for k = 1:length(ft.branches.fid)
				fprintf(file_t, '%2d %2d   %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', ft.branches.fid(k), ft.branches.tid(k), abs(ft.branches.Ia(k)), '∠', angle(ft.branches.Ia(k)).*180./pi, abs(ft.branches.Ib(k)), '∠', angle(ft.branches.Ib(k)).*180./pi, abs(ft.branches.Ic(k)), '∠', angle(ft.branches.Ic(k)).*180./pi);
			end
			fprintf(file_t, '\n');

			fclose(file_t);
		end
		

	end
end