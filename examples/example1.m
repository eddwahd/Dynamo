% Example 1 for DPG 2012 Stuttgart

randseed(23483);

% 3-qubit Ising chain, XY controls, QFT gate
dyn = test_suite(6);
dyn.easy_control(0.1 * ones(1,8));

dyn.ui_open();
pause(3);

dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
