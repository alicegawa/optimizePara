NRNOPT=\
" -c gnabar_tmp=0.12"\
" -c gkdrbar_tmp=0.0175"\
" -c gl_tmp=0.0003"\
" -c gcatbar_tmp=0.0003"\
" -c gkcabar_tmp=0.01"\
" -c ena_tmp=55"\
" -c ek_tmp=-70"\
" -c el_tmp=-50"\
" -c eca_tmp=124"\
" -c stim_amp=100"\
" -c cell_size_rate=12.5"\
" -c gA_tmp=0.022"\
" -c a_sk_tmp=11e-2"\
" -c offset_tmp=0.0"\
" -c coef_log_tmp=3"\
" -c tau_sk_tmp=1200"\
" -c coef_scale_C_tmp=0.75"\
" -c coef_intercept_tmp=1.120"\
" -c power_q_tmp=1.5"\
" -c qinf_scale_tmp=1.15"

#nrnivmodl
nrngui ${NRNOPT} main.hoc -nobanner
