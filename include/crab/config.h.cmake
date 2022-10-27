#ifndef _CRAB_CONFIG_H_
#define _CRAB_CONFIG_H_

/** Define whether lin-ldd is available */
#cmakedefine HAVE_LDD ${HAVE_LDD}

/** Define whether apron library is available */
#cmakedefine HAVE_APRON ${HAVE_APRON}

/** Define whether pplite library is available */
#cmakedefine HAVE_PPLITE ${HAVE_PPLITE}

/** Define whether elina library is available */
#cmakedefine HAVE_ELINA ${HAVE_ELINA}

/** Define whether disable logging for debugging purposes */
#cmakedefine NCRABLOG ${NCRABLOG}

/** Define whether collecting statistics */
#cmakedefine CRAB_STATS ${CRAB_STATS}

/** Use a generic wrapper for abstract domains in tests **/
#cmakedefine USE_GENERIC_WRAPPER ${USE_GENERIC_WRAPPER}

#endif
