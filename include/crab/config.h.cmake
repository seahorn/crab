#ifndef _CRAB_CONFIG_H_
#define _CRAB_CONFIG_H_

/** Define whether lin-ldd is available */
#cmakedefine HAVE_LDD ${HAVE_LDD}

/** Define whether apron library is available */
#cmakedefine HAVE_APRON ${HAVE_APRON}

/**
* Define whether pplite library is available.
* If HAVE_PPLITE_DOMAINS is not enabled then pplite domains are used
* using the apron interface.
**/
#cmakedefine HAVE_PPLITE ${HAVE_PPLITE}

/** 
 * Define whether pplite domains are available. There are two ways of
 * using pplite domains: via apron interface (both HAVE_APRON and
 * HAVE_PPLITE must be enabled) or using a native interface. If
 * HAVE_PPLITE and HAVE_PPLITE_DOMAINS are enabled then the latter kind of
 * interface is also enabled.
 **/
#cmakedefine HAVE_PPLITE_DOMAINS ${HAVE_PPLITE_DOMAINS} 

/** Define whether elina library is available */
#cmakedefine HAVE_ELINA ${HAVE_ELINA}

/** Define whether disable logging for debugging purposes */
#cmakedefine NCRABLOG ${NCRABLOG}

/** Define whether collecting statistics */
#cmakedefine CRAB_STATS ${CRAB_STATS}

/** Use a generic wrapper for abstract domains in tests **/
#cmakedefine USE_GENERIC_WRAPPER ${USE_GENERIC_WRAPPER}

#endif
