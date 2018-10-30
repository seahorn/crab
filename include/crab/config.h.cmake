#ifndef _CRAB_CONFIG_H_
#define _CRAB_CONFIG_H_

/** Define whether lin-ldd is available */
#cmakedefine HAVE_LDD ${HAVE_LDD}

/** Define whether mdd is available */
#cmakedefine HAVE_MDD ${HAVE_MDD}

/** Define whether apron library is available */
#cmakedefine HAVE_APRON ${HAVE_APRON}

/** Define whether elina library is available */
#cmakedefine HAVE_ELINA ${HAVE_ELINA}

/** Define whether enable statistics about the analysis */
#cmakedefine HAVE_STATS ${HAVE_STATS}

/** Define whether enable logging for debugging purposes */
#cmakedefine NCRABLOG ${NCRABLOG}

#endif
