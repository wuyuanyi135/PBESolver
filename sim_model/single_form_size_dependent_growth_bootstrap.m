clear all;

props = parameterized_system_properties;
props.sizeGrids = size_grid(199, 200);
props.gBeta = -5e-3;

ic = make_states(props.solubility(60), 0, props);
options = make_options();

s = hrfvm_solver(options, props);
s.set_initial_states(ic);

sampleTime = 30;

model = 'single_form_model';
open_system(model);
sim(model);