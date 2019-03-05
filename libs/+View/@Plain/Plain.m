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
			fprintf(file_t, '%s\n', '                          Power Flow');
			fprintf(file_t, '%s\n', '=========================================================================');

			fprintf(file_t, '\n%s%d%s%d%s%d%s%d\n','Bus Count       PQ / PV / ref / sum      ', ss.n_pq, ' / ', ss.n_pv, ' / ', ss.n_slack, ' / ',length(ss.bus.id));
			fprintf(file_t, '%s%d%s%d%s%d\n', 'Branches Count       line / transformer / sum      ', ss.n_line, ' / ', ss.n_trans, ' / ', length(ss.branch.id));
			fprintf(file_t, '%s%d%s%d%s%d\n', 'Generator Count       PV / ref / sum      ', ss.pvGens, ' / ', ss.slackGens, ' / ', ss.pvGens+ss.slackGens);
			fprintf(file_t, '%s%d\n', 'Loads Count                                  ', ss.n_load);

			fprintf(file_t, '\n%s%s%d\n', 'solver.method: ', result.method);

			if result.status == 101
				fprintf(file_t, '\n%s%d%s\n', 'error ', result.status, ': Number of iterations exceeds the limit');
				fclose(file_t);
				return;
			elseif result.status == 102
				fprintf(file_t, '\n%s%d%s\n', 'error ', result.status, ': Some node have abnormal voltages and iteration is terminated');
				fclose(file_t);
				return;
			elseif result.status ~= 1
				fprintf(file_t, '\n%s%d\n', 'Unknown error: ', result.status);
				fclose(file_t);
				return;
			end

			fprintf(file_t, '\n%s%d\n', 'Iteration Count: ',result.it);

			fprintf(file_t, '\n%s\n', 'Nodes Votage, Power');
			fprintf(file_t, '\n%s\n', ' id   Votage    Angle    Pg      Qg      Pd      Qd      Pc     Qc');
			for k = 1:length(ss.bus.id)
				fprintf(file_t, '  %d',ss.bus.id(k));
				fprintf(file_t, '   %6.3f',ss.bus.mag(k),ss.bus.ang(k).*180./pi,ss.bus.Pg(k).*100,ss.bus.Qg(k).*100,ss.bus.Pd(k).*100,ss.bus.Qd(k).*100,ss.bus.Pc(k).*100,ss.bus.Qc(k).*100);
				fprintf(file_t, '\n');
			end
			fprintf(file_t, '%s%6.3f  %6.3f   %6.3f  %6.3f   %6.3f  %6.3f','Sum:                     ',sum(ss.bus.Pg).*100,sum(ss.bus.Qg).*100,sum(ss.bus.Pd).*100,sum(ss.bus.Qd).*100,sum(ss.bus.Pc.*100),sum(ss.bus.Qc.*100));
			
			fprintf(file_t, '\n');

			fprintf(file_t, '\n%s\n\n', 'Branch Power:');
			fprintf(file_t, '%s\n', 'id from to    Pij      Qij        Pji      Qji        dP      dQ');
			for k = ss.branch.id'
				fprintf(file_t, '  %d',[ss.branch.id(k),ss.branch.fid(k),ss.branch.tid(k)]);
				fprintf(file_t, '   %6.3f',[ss.branch.Pij(k),ss.branch.Qij(k),ss.branch.Pji(k),ss.branch.Qji(k),ss.branch.dP(k),ss.branch.dQ(k)]*100);
				fprintf(file_t, '\n');
			end
			
			fprintf(file_t, 'Total Loss');
			fprintf(file_t, '   %6.3f',sum(ss.branch.dP)*100,sum(ss.branch.dQ)*100);

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

			fprintf(file_t, '\n%s%d%s%d%s%d%s%d\n','节点数       PQ / PV / 平衡 / 合计      ',ft.ss.n_pq,' / ',ft.ss.n_pv,' / ',ft.ss.n_slack,' / ',length(ft.ss.bus.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','支路数       普通线路 / 变压器 / 合计      ',ft.ss.n_line,' / ',ft.ss.n_trans,' / ',length(ft.ss.branch.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','发电机数       调压厂 / 调频厂 / 合计      ',ft.ss.pvGens,' / ',ft.ss.slackGens,' / ',ft.ss.pvGens+ft.ss.slackGens);
			fprintf(file_t, '%s%d\n','负荷数                                  ',ft.ss.n_load);

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

			for k = 1:length(ft.ss.bus.id)
				fprintf(file_t, '%2d  %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', ft.ss.bus.id(k), abs(ft.bus.Ua(k)), '∠', angle(ft.bus.Ua(k)).*180./pi, abs(ft.bus.Ub(k)), '∠', angle(ft.bus.Ub(k)).*180./pi, abs(ft.bus.Uc(k)), '∠', angle(ft.bus.Uc(k)).*180./pi);
			end
			fprintf(file_t, '\n');

			for k = 1:length(ft.branch.fid)
				fprintf(file_t, '%2d %2d   %6.3f%s%6.3f   %6.3f%s%6.3f   %6.3f%s%6.3f\n', ft.branch.fid(k), ft.branch.tid(k), abs(ft.branch.Ia(k)), '∠', angle(ft.branch.Ia(k)).*180./pi, abs(ft.branch.Ib(k)), '∠', angle(ft.branch.Ib(k)).*180./pi, abs(ft.branch.Ic(k)), '∠', angle(ft.branch.Ic(k)).*180./pi);
			end
			fprintf(file_t, '\n');

			fclose(file_t);
		end
		

	end
end