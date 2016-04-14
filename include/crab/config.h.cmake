#ifndef _CRAB_CONFIG_H_
#define _CRAB_CONFIG_H_

/** Define whether lin-ldd is available */
#cmakedefine HAVE_LDD ${HAVE_LDD}

/** Define whether apron library is available */
#cmakedefine HAVE_APRON ${HAVE_APRON}

/** Define whether enable statistics about the analysis */
#cmakedefine HAVE_STATS ${HAVE_STATS}

/** Define whether enable logging for debugging purposes */
#cmakedefine NCRABLOG ${NCRABLOG}

#endif
