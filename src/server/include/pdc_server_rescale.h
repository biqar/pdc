/*
 * Copyright Notice for
 * Proactive Data Containers (PDC) Software Library and Utilities
 * -----------------------------------------------------------------------------

 *** Copyright Notice ***

 * Proactive Data Containers (PDC) Copyright (c) 2017, The Regents of the
 * University of California, through Lawrence Berkeley National Laboratory,
 * UChicago Argonne, LLC, operator of Argonne National Laboratory, and The HDF
 * Group (subject to receipt of any required approvals from the U.S. Dept. of
 * Energy).  All rights reserved.

 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

 * NOTICE.  This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, distribute copies to the public, prepare derivative works, and to
 * perform publicly and display publicly, and to permit other to do so.
 */

#ifndef PDC_SERVER_RESCALE_H
#define PDC_SERVER_RESCALE_H

#include "pdc_public.h"

/**
 * Restart with a different server count than the checkpoint (elastic).
 *
 * Eagerly repartitions metadata from N_old checkpoint shards onto N_new ranks.
 * Unpublishes stale server.cfg for the duration of migration, ensures
 * per-rank tmp dirs exist, and barriers before returning so the caller can
 * publish server.cfg only after metadata is ready.
 *
 * \param n_old [IN] Server count at checkpoint time
 * \param n_new [IN] Current server MPI size
 *
 * \return Non-negative on success/Negative on failure
 */
perr_t PDC_Server_restart_elastic(int n_old, int n_new);

#endif /* PDC_SERVER_RESCALE_H */
