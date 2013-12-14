function res = spin_chain_plot()
% Simulating exciton transport on a spin chain/network.

% Ville Bergholm 2011-2012



%% controls

% Control Hamiltonians / Liouvillians
H_ctrl = {};

H2 = op_list({{n_op, 2}}, dim);
H_ctrl{end+1} = H2(p,p);
H_ctrl{end+1} = H_int;
H_ctrl{end+1} = superop_lindblad(sink);

for k = 2  % 1:n_sites
    temp = op_list({{sqrt(2) * n_op, k}}, dim);
    H_ctrl{end+1} = superop_lindblad({temp(p,p)}); % dephasing controls
end

% transformed controls?
control_type = '..pp';
control_par = {};
c_labels = {'Z_2', 'v', '\sigma', 'D_2'};



%% set up controls
%T = 150;
T = 10;
dyn.seq_init(1, T * [0.5, 1.0], control_type, control_par);

%split = linspace(0,5,25);
split = linspace(0,0,1);
%v = logspace(-1.5, 0.6, 19);
v = logspace(0, 0, 1);
suck = logspace(-2, 1.1, 17);
dep = linspace(0, 4, 23);

for aa = 1:length(split)
  aa
  for bb = 1:length(v)
    for cc = 1:length(suck)
      for dd = 1:length(dep)
          dyn.easy_control([split(aa), v(bb), suck(cc), dep(dd)]);
          temp = dyn.X();
          res(aa, bb, cc, dd) = temp(end);
      end
    end
  end
end
res = real(squeeze(res));

% no dephasing
d0 = res(:, 1)*ones(1, length(dep));

aaa = res -d0;

figure();
subplot(1, 2 ,1);
surf(dep, suck, aaa);
xlabel('dep')
ylabel('suck')
title('absolute')

rrr = res ./ d0;
subplot(1,2,2);
surf(dep, suck, rrr);
xlabel('dep')
ylabel('suck')
title('relative')

end
