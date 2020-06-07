clear all;

props = parameterized_system_properties;

ic = make_states(props.solubility(60), 0, props);
options = make_options();

s = hrfvm_solver(options, props);
s.set_initial_states(ic);

sampleTime = 30;

cfl = 0.1;

model = 'single_form_model';
open_system(model);
sim(model);