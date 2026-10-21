from mpmath import mpc, mpf  # type: ignore

trusted_dict = {
    "cmax": mpf("0.0"),
    "cmin": mpf("0.207412234081733157069792714734913"),
    "cminus": mpf("-0.0200518879643783869454456358447974"),
    "cplus": mpf("-0.0200518879643783869454456358447974"),
    "subsonic_cmax": mpc(real="0.640378026586350675783335076284507", imag="0.0"),
    "subsonic_cmin": mpf("0.695652173913043478260869565217378"),
    "subsonic_exchanged_cmax": mpc(
        real="0.640378026586350675783335076284507", imag="0.0"
    ),
    "subsonic_exchanged_cmin": mpf("0.695652173913043478260869565217378"),
}
