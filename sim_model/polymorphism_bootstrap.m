clear all;

props1 = parameterized_system_properties;
props2 = parameterized_system_properties;
props2.solubilityPoly = [6.222e-3; -1.165e-4; 7.644e-6];
props2.pnKp = 1e7;
props = [props1 props2];
ic = [
    make_states(props1.solubility(60), 0, props1), ...
    make_states(nan, 0, props2)
];

options = make_options();

s = hrfvm_solver(options, props);
s.set_initial_states(ic);

sampleTime = 30;

model = 'polymorphism';
open_system(model);
sim(model);