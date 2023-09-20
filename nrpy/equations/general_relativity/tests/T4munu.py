from mpmath import mpf  # type: ignore

trusted_dict = {
    "ADMS": mpf("5.51886985348101528899497640423776"),
    "ADMSD_0": mpf("0.608341733349877488954095315599876"),
    "ADMSD_1": mpf("0.470565344820866964820125587278172"),
    "ADMSD_2": mpf("1.12342071849067535252360437426405"),
    "ADMSDD_0": mpf("1.80029784665212934056799078226269"),
    "ADMSDD_1": mpf("1.39659583030126451047804795612817"),
    "ADMSDD_2": mpf("3.00811244695977724176245476191804"),
    "ADMSDD_3": mpf("1.39659583030126451047804795612817"),
    "ADMSDD_4": mpf("1.22494659231382375810150001831715"),
    "ADMSDD_5": mpf("2.50094687070838058316104348117154"),
    "ADMSDD_6": mpf("3.00811244695977724176245476191804"),
    "ADMSDD_7": mpf("2.50094687070838058316104348117154"),
    "ADMSDD_8": mpf("5.32511755805920028174659119741427"),
    "ADMrho": mpf("0.128468170910047767321502802016049"),
    "BSSNSDD_Cartesian_0": mpf("431.750880189457681158121637267736"),
    "BSSNSDD_Cartesian_1": mpf("639.598051407228311087230496889343"),
    "BSSNSDD_Cartesian_2": mpf("491.81362621412285625346833932248"),
    "BSSNSDD_Cartesian_3": mpf("639.598051407228311087230496889343"),
    "BSSNSDD_Cartesian_4": mpf("973.782230183346755946092789928739"),
    "BSSNSDD_Cartesian_5": mpf("775.855652178151041197421958821041"),
    "BSSNSDD_Cartesian_6": mpf("491.81362621412285625346833932248"),
    "BSSNSDD_Cartesian_7": mpf("775.855652178151041197421958821041"),
    "BSSNSDD_Cartesian_8": mpf("645.305639210365512496941943299871"),
    "BSSNSDD_SinhCartesian_0": mpf("72.8977328242746176639616073801642"),
    "BSSNSDD_SinhCartesian_1": mpf("113.885079536085038057679833543104"),
    "BSSNSDD_SinhCartesian_2": mpf("95.1044047685209926345437442534337"),
    "BSSNSDD_SinhCartesian_3": mpf("113.885079536085038057679833543104"),
    "BSSNSDD_SinhCartesian_4": mpf("183.226088848670817825778109621973"),
    "BSSNSDD_SinhCartesian_5": mpf("158.672404672199783732918558465882"),
    "BSSNSDD_SinhCartesian_6": mpf("95.1044047685209926345437442534337"),
    "BSSNSDD_SinhCartesian_7": mpf("158.672404672199783732918558465882"),
    "BSSNSDD_SinhCartesian_8": mpf("143.2746641202630371913526180672"),
    "BSSNSDD_SinhCylindrical_0": mpf("7.53852414307662371928058413808401"),
    "BSSNSDD_SinhCylindrical_1": mpf("4.90154623084854548312165901955615"),
    "BSSNSDD_SinhCylindrical_2": mpf("15.2439589588365485221361901629701"),
    "BSSNSDD_SinhCylindrical_3": mpf("4.90154623084854548312165901955615"),
    "BSSNSDD_SinhCylindrical_4": mpf("3.39340215814975643509487409929643"),
    "BSSNSDD_SinhCylindrical_5": mpf("11.2216826196890476851424373733739"),
    "BSSNSDD_SinhCylindrical_6": mpf("15.2439589588365485221361901629701"),
    "BSSNSDD_SinhCylindrical_7": mpf("11.2216826196890476851424373733739"),
    "BSSNSDD_SinhCylindrical_8": mpf("39.140052324420898807433150545942"),
    "BSSNSDD_SinhSpherical_0": mpf("0.0000000393733202033316805115342791987653"),
    "BSSNSDD_SinhSpherical_1": mpf("0.0000000027601180182480318681704633917624"),
    "BSSNSDD_SinhSpherical_2": mpf("0.000000000248636364508742235222310018197991"),
    "BSSNSDD_SinhSpherical_3": mpf("0.0000000027601180182480318681704633917624"),
    "BSSNSDD_SinhSpherical_4": mpf("0.000000000193613918050073217019757525878161"),
    "BSSNSDD_SinhSpherical_5": mpf("1.75324006153874187833080675934564e-11"),
    "BSSNSDD_SinhSpherical_6": mpf("0.000000000248636364508742235222310018197991"),
    "BSSNSDD_SinhSpherical_7": mpf("1.75324006153874187833080675934564e-11"),
    "BSSNSDD_SinhSpherical_8": mpf("1.65362193943888305125959490976609e-12"),
    "BSSNSDD_SinhSpherical_rfm_precompute_0": mpf("83.8092449570758956994645831820557"),
    "BSSNSDD_SinhSpherical_rfm_precompute_1": mpf("208.489577273024183937453244098759"),
    "BSSNSDD_SinhSpherical_rfm_precompute_2": mpf("148.702280361409432078106375182902"),
    "BSSNSDD_SinhSpherical_rfm_precompute_3": mpf("208.489577273024183937453244098759"),
    "BSSNSDD_SinhSpherical_rfm_precompute_4": mpf("532.32615663431370254377872772267"),
    "BSSNSDD_SinhSpherical_rfm_precompute_5": mpf("391.09611920496502372525011540698"),
    "BSSNSDD_SinhSpherical_rfm_precompute_6": mpf("148.702280361409432078106375182902"),
    "BSSNSDD_SinhSpherical_rfm_precompute_7": mpf("391.09611920496502372525011540698"),
    "BSSNSDD_SinhSpherical_rfm_precompute_8": mpf("296.634009370776642256497878225714"),
    "BSSNSDD_SinhSymTP_0": mpf("2.04574014413885427003428219800366"),
    "BSSNSDD_SinhSymTP_1": mpf("1.90668669633741144868801540111256"),
    "BSSNSDD_SinhSymTP_2": mpf("0.315283577830810661379730793426764"),
    "BSSNSDD_SinhSymTP_3": mpf("1.90668669633741144868801540111256"),
    "BSSNSDD_SinhSymTP_4": mpf("1.7908583035036719540404615592937"),
    "BSSNSDD_SinhSymTP_5": mpf("0.301817723211503859672958613783364"),
    "BSSNSDD_SinhSymTP_6": mpf("0.315283577830810661379730793426764"),
    "BSSNSDD_SinhSymTP_7": mpf("0.301817723211503859672958613783364"),
    "BSSNSDD_SinhSymTP_8": mpf("0.0531971484170274621199899577835178"),
    "BSSNSDD_Spherical_0": mpf("224.623568143102753656631770760895"),
    "BSSNSDD_Spherical_1": mpf("103.208599506892176051867702346219"),
    "BSSNSDD_Spherical_2": mpf("19.7262218354843803613250277601632"),
    "BSSNSDD_Spherical_3": mpf("103.208599506892176051867702346219"),
    "BSSNSDD_Spherical_4": mpf("47.7551076710812409920598206432526"),
    "BSSNSDD_Spherical_5": mpf("9.33442085570287296938888266904388"),
    "BSSNSDD_Spherical_6": mpf("19.7262218354843803613250277601632"),
    "BSSNSDD_Spherical_7": mpf("9.33442085570287296938888266904388"),
    "BSSNSDD_Spherical_8": mpf("1.95217531815959466524648349955833"),
    "BSSNSD_Cartesian_0": mpf("10.7734748958192342615944288932221"),
    "BSSNSD_Cartesian_1": mpf("15.7894350419011412174522146535946"),
    "BSSNSD_Cartesian_2": mpf("11.9833994963534688294694647241315"),
    "BSSNSD_SinhCartesian_0": mpf("4.43326825805536144052642179230332"),
    "BSSNSD_SinhCartesian_1": mpf("6.83123148120768755356615033975935"),
    "BSSNSD_SinhCartesian_2": mpf("5.61169345309589770589978527355091"),
    "BSSNSD_SinhCylindrical_0": mpf("1.33667471613855107814219443868593"),
    "BSSNSD_SinhCylindrical_1": mpf("0.857220541597304646136816975719213"),
    "BSSNSD_SinhCylindrical_2": mpf("2.63026445975496045721569083385368"),
    "BSSNSD_SinhSpherical_0": mpf("0.0000904880760661347303714381214226987"),
    "BSSNSD_SinhSpherical_1": mpf("0.00000642042751321339419804664452681651"),
    "BSSNSD_SinhSpherical_2": mpf("0.000000634872462065772986025965806323184"),
    "BSSNSD_SinhSpherical_rfm_precompute_0": mpf("4.86553876960285578171356920670663"),
    "BSSNSD_SinhSpherical_rfm_precompute_1": mpf("11.8440190092397273349573931269702"),
    "BSSNSD_SinhSpherical_rfm_precompute_2": mpf("8.23969614196673749527963838614753"),
    "BSSNSD_SinhSymTP_0": mpf("0.73318753013934643125134958563747"),
    "BSSNSD_SinhSymTP_1": mpf("0.696048255179012585974736649312013"),
    "BSSNSD_SinhSymTP_2": mpf("0.120539123913523654867022767336316"),
    "BSSNSD_Spherical_0": mpf("7.36139810045660472593911135128682"),
    "BSSNSD_Spherical_1": mpf("3.46174430637807807205095632989364"),
    "BSSNSD_Spherical_2": mpf("0.712287187722210421141804980625848"),
    "BSSNS_Cartesian": mpf("90.0702869584517007272615688426028"),
    "BSSNS_SinhCartesian": mpf("39.6427383777813476245746084129621"),
    "BSSNS_SinhCylindrical": mpf("13.2580129513789858126411115901235"),
    "BSSNS_SinhSpherical": mpf("0.000274747361153781438733672886733984"),
    "BSSNS_SinhSpherical_rfm_precompute": mpf("56.8012148582018100616068502326587"),
    "BSSNS_SinhSymTP": mpf("3.49792984128431817042693295642172"),
    "BSSNS_Spherical": mpf("28.101490602542897599054460554267"),
    "BSSNrho_Cartesian": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_SinhCartesian": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_SinhCylindrical": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_SinhSpherical": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_SinhSpherical_rfm_precompute": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_SinhSymTP": mpf("0.128468170910047767321502802016049"),
    "BSSNrho_Spherical": mpf("0.128468170910047767321502802016049"),
    "BSSNsourceterm_H_Cartesian": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_H_SinhCartesian": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_H_SinhCylindrical": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_H_SinhSpherical": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_H_SinhSpherical_rfm_precompute": mpf(
        "-6.45751459121798477210713792622794"
    ),
    "BSSNsourceterm_H_SinhSymTP": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_H_Spherical": mpf("-6.45751459121798477210713792622794"),
    "BSSNsourceterm_MU_Cartesian_0": mpf("-14.6851347831655898600098764117899"),
    "BSSNsourceterm_MU_Cartesian_1": mpf("-15.7435784039388779445151487154523"),
    "BSSNsourceterm_MU_Cartesian_2": mpf("-12.0831625987747640505239328645715"),
    "BSSNsourceterm_MU_SinhCartesian_0": mpf("-14.6851347831655898600098764119477"),
    "BSSNsourceterm_MU_SinhCartesian_1": mpf("-15.7435784039388779445151487157695"),
    "BSSNsourceterm_MU_SinhCartesian_2": mpf("-12.0831625987747640505239328648902"),
    "BSSNsourceterm_MU_SinhCylindrical_0": mpf("-14.6851347831655898600098764124242"),
    "BSSNsourceterm_MU_SinhCylindrical_1": mpf("-15.7435784039388779445151487148181"),
    "BSSNsourceterm_MU_SinhCylindrical_2": mpf("-12.0831625987747640505239328648902"),
    "BSSNsourceterm_MU_SinhSpherical_0": mpf("-14.6851347831655898600098764111557"),
    "BSSNsourceterm_MU_SinhSpherical_1": mpf("-15.7435784039388779445151487075211"),
    "BSSNsourceterm_MU_SinhSpherical_2": mpf("-12.0831625987747640505239328661587"),
    "BSSNsourceterm_MU_SinhSpherical_rfm_precompute_0": mpf(
        "-14.6851347831655898600098764114728"
    ),
    "BSSNsourceterm_MU_SinhSpherical_rfm_precompute_1": mpf(
        "-15.743578403938877944515148715293"
    ),
    "BSSNsourceterm_MU_SinhSpherical_rfm_precompute_2": mpf(
        "-12.0831625987747640505239328645715"
    ),
    "BSSNsourceterm_MU_SinhSymTP_0": mpf("-14.6851347831655898600098764117899"),
    "BSSNsourceterm_MU_SinhSymTP_1": mpf("-15.7435784039388779445151487151352"),
    "BSSNsourceterm_MU_SinhSymTP_2": mpf("-12.0831625987747640505239328648902"),
    "BSSNsourceterm_MU_Spherical_0": mpf("-14.6851347831655898600098764114728"),
    "BSSNsourceterm_MU_Spherical_1": mpf("-15.7435784039388779445151487156101"),
    "BSSNsourceterm_MU_Spherical_2": mpf("-12.0831625987747640505239328661587"),
    "BSSNsourceterm_a_rhs_Cartesian_0": mpf("-272.36279527541105793814990532571"),
    "BSSNsourceterm_a_rhs_Cartesian_1": mpf("-648.209795086325187337935928228514"),
    "BSSNsourceterm_a_rhs_Cartesian_2": mpf("-754.132367561873490564255019800344"),
    "BSSNsourceterm_a_rhs_Cartesian_3": mpf("-648.209795086325187337935928228514"),
    "BSSNsourceterm_a_rhs_Cartesian_4": mpf("-1037.14518926649593277689227290549"),
    "BSSNsourceterm_a_rhs_Cartesian_5": mpf("-870.661424224684416554152825430175"),
    "BSSNsourceterm_a_rhs_Cartesian_6": mpf("-754.132367561873490564255019800344"),
    "BSSNsourceterm_a_rhs_Cartesian_7": mpf("-870.661424224684416554152825430175"),
    "BSSNsourceterm_a_rhs_Cartesian_8": mpf("-417.25744575790123570346214536685"),
    "BSSNsourceterm_a_rhs_SinhCartesian_0": mpf("-42.1984512085968786467829173206262"),
    "BSSNsourceterm_a_rhs_SinhCartesian_1": mpf("-113.258908862692832917641486418133"),
    "BSSNsourceterm_a_rhs_SinhCartesian_2": mpf("-145.798802945544471920734143707234"),
    "BSSNsourceterm_a_rhs_SinhCartesian_3": mpf("-113.258908862692832917641486418133"),
    "BSSNsourceterm_a_rhs_SinhCartesian_4": mpf("-194.365532870989639284181458568841"),
    "BSSNsourceterm_a_rhs_SinhCartesian_5": mpf("-179.697173751105850619543376353468"),
    "BSSNsourceterm_a_rhs_SinhCartesian_6": mpf("-145.798802945544471920734143707234"),
    "BSSNsourceterm_a_rhs_SinhCartesian_7": mpf("-179.697173751105850619543376353468"),
    "BSSNsourceterm_a_rhs_SinhCartesian_8": mpf("-99.1052447167660688946622468163766"),
    "BSSNsourceterm_a_rhs_SinhCylindrical_0": mpf(
        "-1.45397444027572547302686010934118"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_1": mpf(
        "-4.12711616518858305310391948747554"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_2": mpf(
        "-23.1927010186583322312638877487699"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_3": mpf(
        "-4.12711616518858305310391948747554"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_4": mpf(
        "-3.43454579711188942245268456321974"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_5": mpf(
        "-13.0020570700721475767678411004863"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_6": mpf(
        "-23.1927010186583322312638877487699"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_7": mpf(
        "-13.0020570700721475767678411004863"
    ),
    "BSSNsourceterm_a_rhs_SinhCylindrical_8": mpf(
        "-33.7424519883372942130948906254923"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_0": mpf(
        "-0.0000000429225025310296526937525193474081"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_1": mpf(
        "-0.00000000307296214747333540388901021915636"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_2": mpf(
        "-0.000000000323183166549657380662986042989557"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_3": mpf(
        "-0.00000000307296214747333540388901021915636"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_4": mpf(
        "-0.000000000169286517718086522745224770578789"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_5": mpf(
        "1.8710237636948711962252291054001e-11"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_6": mpf(
        "-0.000000000323183166549657380662986042989557"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_7": mpf(
        "1.8710237636948711962252291054001e-11"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_8": mpf(
        "3.11174187868808313827436243819116e-11"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_0": mpf(
        "-41.6067084172314903605419838259129"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_1": mpf(
        "-202.083724133629780128465724235755"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_2": mpf(
        "-228.098256016170148724555984336743"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_3": mpf(
        "-202.083724133629780128465724235755"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_4": mpf(
        "-563.061495640209044391952061286479"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_5": mpf(
        "-449.69326593950905163715980081976"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_6": mpf(
        "-228.098256016170148724555984336743"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_7": mpf(
        "-449.69326593950905163715980081976"
    ),
    "BSSNsourceterm_a_rhs_SinhSpherical_rfm_precompute_8": mpf(
        "-222.843108579839266464902846044345"
    ),
    "BSSNsourceterm_a_rhs_SinhSymTP_0": mpf("-2.0164994104404758900175026436364"),
    "BSSNsourceterm_a_rhs_SinhSymTP_1": mpf("-2.16626296963406298049865130504608"),
    "BSSNsourceterm_a_rhs_SinhSymTP_2": mpf("-0.477213312339147762485212819050401"),
    "BSSNsourceterm_a_rhs_SinhSymTP_3": mpf("-2.16626296963406298049865130504608"),
    "BSSNsourceterm_a_rhs_SinhSymTP_4": mpf("-1.88975840860663646320822118268004"),
    "BSSNsourceterm_a_rhs_SinhSymTP_5": mpf("-0.256592527407823331171936893394853"),
    "BSSNsourceterm_a_rhs_SinhSymTP_6": mpf("-0.477213312339147762485212819050401"),
    "BSSNsourceterm_a_rhs_SinhSymTP_7": mpf("-0.256592527407823331171936893394853"),
    "BSSNsourceterm_a_rhs_SinhSymTP_8": mpf("0.0334802528293420385756321652979648"),
    "BSSNsourceterm_a_rhs_Spherical_0": mpf("-232.387365500243824565824998384881"),
    "BSSNsourceterm_a_rhs_Spherical_1": mpf("-117.753637654769706365679213470967"),
    "BSSNsourceterm_a_rhs_Spherical_2": mpf("-29.3473116713143464005108563027435"),
    "BSSNsourceterm_a_rhs_Spherical_3": mpf("-117.753637654769706365679213470967"),
    "BSSNsourceterm_a_rhs_Spherical_4": mpf("-48.6888700345953972086784034297805"),
    "BSSNsourceterm_a_rhs_Spherical_5": mpf("-5.81988211300000809853539617675767"),
    "BSSNsourceterm_a_rhs_Spherical_6": mpf("-29.3473116713143464005108563027435"),
    "BSSNsourceterm_a_rhs_Spherical_7": mpf("-5.81988211300000809853539617675767"),
    "BSSNsourceterm_a_rhs_Spherical_8": mpf("3.55979009527172276256549320438448"),
    "BSSNsourceterm_lambda_rhs_Cartesian_0": mpf("-148.598487846873169966854167723318"),
    "BSSNsourceterm_lambda_rhs_Cartesian_1": mpf("-159.308850662091081627260480295663"),
    "BSSNsourceterm_lambda_rhs_Cartesian_2": mpf("-122.269200596248804966035544651624"),
    "BSSNsourceterm_lambda_rhs_SinhCartesian_0": mpf(
        "-148.598487846873169966854167718472"
    ),
    "BSSNsourceterm_lambda_rhs_SinhCartesian_1": mpf(
        "-159.308850662091081627260480295663"
    ),
    "BSSNsourceterm_lambda_rhs_SinhCartesian_2": mpf(
        "-122.269200596248804966035544658086"
    ),
    "BSSNsourceterm_lambda_rhs_SinhCylindrical_0": mpf(
        "-148.598487846873169966854167721703"
    ),
    "BSSNsourceterm_lambda_rhs_SinhCylindrical_1": mpf(
        "-159.308850662091081627260480302125"
    ),
    "BSSNsourceterm_lambda_rhs_SinhCylindrical_2": mpf(
        "-122.269200596248804966035544655663"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_0": mpf(
        "-148.598487846873169966854167720693"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_1": mpf(
        "-159.308850662091081627260480302125"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_2": mpf(
        "-122.269200596248804966035544680705"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_rfm_precompute_0": mpf(
        "-148.59848784687316996685416772655"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_rfm_precompute_1": mpf(
        "-159.308850662091081627260480295663"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSpherical_rfm_precompute_2": mpf(
        "-122.26920059624880496603554465324"
    ),
    "BSSNsourceterm_lambda_rhs_SinhSymTP_0": mpf("-148.598487846873169966854167719279"),
    "BSSNsourceterm_lambda_rhs_SinhSymTP_1": mpf("-159.308850662091081627260480294047"),
    "BSSNsourceterm_lambda_rhs_SinhSymTP_2": mpf("-122.269200596248804966035544654855"),
    "BSSNsourceterm_lambda_rhs_Spherical_0": mpf("-148.598487846873169966854167723318"),
    "BSSNsourceterm_lambda_rhs_Spherical_1": mpf("-159.308850662091081627260480302125"),
    "BSSNsourceterm_lambda_rhs_Spherical_2": mpf("-122.269200596248804966035544661318"),
    "BSSNsourceterm_trK_rhs_Cartesian": mpf("651.152475668214327117588015149391"),
    "BSSNsourceterm_trK_rhs_SinhCartesian": mpf("287.111718641249828853051464238997"),
    "BSSNsourceterm_trK_rhs_SinhCylindrical": mpf("96.6381444041333056666107491569852"),
    "BSSNsourceterm_trK_rhs_SinhSpherical": mpf("0.929406063717995650148152606652584"),
    "BSSNsourceterm_trK_rhs_SinhSpherical_rfm_precompute": mpf(
        "410.980217966953019723836167593915"
    ),
    "BSSNsourceterm_trK_rhs_SinhSymTP": mpf("26.1792753127474524478915847051603"),
    "BSSNsourceterm_trK_rhs_Spherical": mpf("203.794470522151924795882820842264"),
}
