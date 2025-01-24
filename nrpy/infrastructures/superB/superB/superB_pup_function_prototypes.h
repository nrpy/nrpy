void pup_commondata_struct(PUP::er &p, commondata_struct &cd);
void pup_params_struct(PUP::er &p, params_struct &ps);
void pup_rfm_struct(PUP::er &p, rfm_struct &rfm, const params_struct *restrict params);
void pup_innerpt_bc_struct(PUP::er &p, innerpt_bc_struct &ibc);
void pup_outerpt_bc_struct(PUP::er &p, outerpt_bc_struct &obc);
void pup_bc_info_struct(PUP::er &p, bc_info_struct &bci);
void pup_bc_struct(PUP::er &p, bc_struct &bc);
void pup_MoL_gridfunctions_struct(PUP::er &p, MoL_gridfunctions_struct &gridfuncs, const params_struct &params, const commondata_struct &commondata);
void pup_charecomm_struct(PUP::er &p, charecomm_struct &cc, const params_struct &params, const params_struct &params_chare);
void pup_diagnostic_struct(PUP::er &p, diagnostic_struct &ds, const int num_diagnostics_chare);
void pup_tmpBuffers_struct(PUP::er &p, tmpBuffers_struct &tmpBuffers, const params_struct &params, const nonlocalinnerbc_struct &nonlocalinnerbc);
void pup_nonlocalinnerbc_struct(PUP::er &p, nonlocalinnerbc_struct &nonlocal, const commondata_struct &commondata);
void pup_griddata(PUP::er &p, griddata_struct &gd);
void pup_griddata_chare(PUP::er &p, griddata_struct &gd, const params_struct &params, const commondata_struct &commondata);

