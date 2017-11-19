#define ALBEDO 1
#define CANMOD 0
#define CONDCT 1
#define DENSTY 1
#define EXCHNG 1
#define HYDROL 1

#define DRIV1D 0
#define SWPART 1

#if ALBEDO == 0
#define ALBEDO_OPT 'diagnostic'
#elif ALBEDO == 1
#define ALBEDO_OPT 'prognostic'
#endif

#if CANMOD == 0
#define CANMOD_OPT 'zero layer'
#elif CANMOD == 1
#define CANMOD_OPT 'one layer'
#endif

#if CONDCT == 0
#define CONDCT_OPT 'constant'
#elif CONDCT == 1
#define CONDCT_OPT 'Yen (1981)'
#endif

#if DENSTY == 0
#define DENSTY_OPT 'constant'
#elif DENSTY == 1
#define DENSTY_OPT 'Verseghy (1991)'
#endif

#if EXCHNG == 0
#define EXCHNG_OPT 'constant'
#elif EXCHNG == 1
#define EXCHNG_OPT 'Louis (1979)'
#endif

#if HYDROL == 0
#define HYDROL_OPT 'free draining'
#elif HYDROL == 1
#define HYDROL_OPT 'bucket'
#endif
