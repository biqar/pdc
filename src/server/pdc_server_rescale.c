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

#include "pdc_config.h"
#include "pdc_interface.h"
#include "pdc_server_rescale.h"
#include "pdc_server_metadata.h"
#include "pdc_logger.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

/*
 * Elastic restart entry point. Later tasks fill in:
 *   T6 shard validation, T7 load, T8/T9 migrate, T10 barrier/server.cfg, ...
 */
perr_t
PDC_Server_restart_elastic(int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    if (n_old <= 0 || n_new <= 0)
        PGOTO_ERROR(FAIL, "Invalid server counts for elastic restart (n_old=%d, n_new=%d)", n_old, n_new);

    if (n_old == n_new)
        PGOTO_ERROR(FAIL, "PDC_Server_restart_elastic called with n_old == n_new (%d); use PDC_Server_restart",
                    n_old);

    if (n_old > 65536 || n_new > 65536)
        PGOTO_ERROR(FAIL, "Server count out of range for elastic restart (n_old=%d, n_new=%d)", n_old, n_new);

    if (pdc_server_rank_g == 0)
        LOG_INFO("Elastic server restart requested: %d -> %d\n", n_old, n_new);

#ifdef ENABLE_MPI
    /* Placeholder: later phases barrier before/after migration and before server.cfg publish */
    MPI_Barrier(MPI_COMM_WORLD);
#else
    if (n_new != 1)
        PGOTO_ERROR(FAIL, "Elastic restart with n_new=%d requires MPI", n_new);
#endif

    /* Migration not implemented yet (T6–T10). Fail after validation so callers get a clear error. */
    if (pdc_server_rank_g == 0)
        LOG_ERROR("Elastic metadata migration (%d -> %d) is not implemented yet\n", n_old, n_new);
    PGOTO_ERROR(FAIL, "Elastic server restart not implemented");

done:
    FUNC_LEAVE(ret_value);
}
