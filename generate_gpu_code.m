% check cuda installation
ret = coder.checkGpuInstall;
assert(ret.cuda, "CUDA is not installed");

cfg = coder.gpuConfig('mex');
cfg.GpuConfig.CompilerFlags = '--fmad=false';
cfg.GenerateReport = true;

codegen -config cfg moment -args 