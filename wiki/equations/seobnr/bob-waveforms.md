# BOB Waveforms

> Map Backwards-One-Body waveform, NQC attachment, higher-mode, and BOBv2 symbolic quantities. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [SEOBNR And BOB](index.md)

## Summary

The BOB waveform modules provide merger-ringdown and NQC-attachment quantities
for aligned-spin binaries. The original BOB `(2,2)` module builds strain
amplitude, phase, and attachment derivatives from mass, spin, QNM, and
attachment-time symbols. The higher-mode module applies the same BOB structure
mode by mode. The BOBv2 module builds an uncalibrated news-to-strain waveform,
stores peak-news fits, and exposes attachment quantities for the implemented
`(2,2)` strain path.

## Detail

`BOB_aligned_spin_waveform_quantities` owns the aligned-spin BOB `(2,2)`
merger-ringdown path. Its constructor declares masses `m1` and `m2`, aligned
spins `chi1` and `chi2`, QNM symbols `omega_qnm` and `tau_qnm`, attachment time
`t_0`, and time `t`. It computes mass ratio and effective spin combinations,
uses fitted peak strain and peak frequency expressions, constructs the BOB
orbital frequency and integrated phase, and stores the final strain amplitude
`h` and phase `phi`.

The `(2,2)` BOB attachment outputs are `h_t_attach`, `hdot_t_attach`,
`hddot_t_attach`, `w_t_attach`, and `wdot_t_attach`. In this implementation
`hdot_t_attach` is explicitly zero, while `hddot_t_attach` and
`wdot_t_attach` are obtained by differentiating the symbolic BOB strain and
frequency expressions and substituting the attachment time.

`BOB_aligned_spin_waveform_quantities_higher_modes` generalizes the BOB
attachment construction to the mode set `(2,2)`, `(3,3)`, `(4,4)`, `(5,5)`,
`(2,1)`, `(3,2)`, and `(4,3)`. It stores modewise input symbols in
`modewise_input_symbols` and stores outputs in dictionaries `h_lm`, `phi_lm`,
`hdot_t_attach_lm`, `hddot_t_attach_lm`, and `wdot_t_attach_lm`. Its validation
entry point flattens these dictionaries into keys such as `h_22`, `phi_33`,
`hddot_t_attach_44`, and `wdot_t_attach55`.

`BOB_v2_waveform_quantities` owns the current BOBv2 path. The constructor stores
input and derived binary parameters `m1`, `m2`, `chi1`, `chi2`, `M`, `nu`,
`chi_s`, `chi_a`, `delta`, and `chi_eob`. It initializes `newsNR`, evaluates
peak-news fit helpers `news_Ap_22`, `news_Ap_33`, `news_Ap_21`, `news_Ap_44`,
`news_Ap_43`, `news_Ap_55`, and `news_Ap_32`, and then consumes only the
`(2,2)` peak-news fit in the downstream BOBv2 waveform construction.

The BOBv2 waveform builds an orbital-frequency expression from the QNM
frequency, converts news to strain with a truncated series, and stores
`strain_amp_deriv` for attachment-point selection, complex strain `h_complex`,
and attachment outputs `h_t_attach`, `hdot_t_attach`, `hddot_t_attach`,
`w_t_attach`, and `wdot_t_attach`. The code documents this as an uncalibrated
BOBv2 implementation whose ideal attachment point still requires calibration.

All three modules validate through the trusted-expression pipeline. The BOB
`(2,2)` trusted dictionary pins `h`, `phi`, and attachment outputs. The
higher-mode trusted dictionary pins flattened mode keys for the supported mode
set. The BOBv2 trusted dictionary removes the raw `newsNR` dictionary before
validation, then adds explicit `news_Ap_*` keys so peak-news fits and attachment
quantities are both covered.

## Sources

- [BOB_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities.py) - `BOB_aligned_spin_waveform_quantities`, `h`, `phi`, `h_t_attach`, `hddot_t_attach`, `w_t_attach`, `wdot_t_attach`
- [BOB_aligned_spin_waveform_quantities_higher_modes.py](../../../nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities_higher_modes.py) - `BOB_aligned_spin_waveform_quantities_higher_modes`, `modes`, `h_lm`, `phi_lm`, `hddot_t_attach_lm`, `wdot_t_attach_lm`
- [BOB_v2_waveform_quantities_kankani_etal.py](../../../nrpy/equations/seobnr/BOB_v2_waveform_quantities_kankani_etal.py) - `BOB_v2_waveform_quantities`, `newsNR`, `news_Ap_22`, `strain_amp_deriv`, `h_complex`, `w_t_attach`
- [BOB_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities.py) - `trusted_dict`
- [BOB_aligned_spin_waveform_quantities_higher_modes.py](../../../nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities_higher_modes.py) - `trusted_dict`
- [BOB_v2_waveform_quantities_kankani_etal.py](../../../nrpy/equations/seobnr/tests/BOB_v2_waveform_quantities_kankani_etal.py) - `trusted_dict`

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- [Precessing Dynamics](precessing-dynamics.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
