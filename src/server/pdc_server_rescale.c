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
 * reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "pdc_config.h"
#include "pdc_interface.h"
#include "pdc_server_rescale.h"
#include "pdc_server_metadata.h"
#include "pdc_client_server_common.h"
#include "pdc_logger.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

/* Cached after successful T6 validation for T7 load planning */
static int rescale_validated_n_old_g  = -1;
static int rescale_validated_n_new_g  = -1;

static void
checkpoint_shard_path(char *buf, size_t buflen, int shard_rank)
{
    FUNC_ENTER(NULL);

    snprintf(buf, buflen, "%s/%d/metadata_checkpoint.%d", pdc_server_tmp_dir_g, shard_rank, shard_rank);

    FUNC_LEAVE_VOID();
}

static int
checkpoint_shard_exists(int shard_rank)
{
    FUNC_ENTER(NULL);

    char path[ADDR_MAX];

    checkpoint_shard_path(path, sizeof(path), shard_rank);

    FUNC_LEAVE(access(path, F_OK) == 0);
}

/*
 * Which rank loads checkpoint shard `shard` under elastic restart.
 * - Shards in [0, min(n_old,n_new)): loaded by that rank (if it exists under n_new).
 * - Orphan shards [n_new, n_old) on scale-down: loaded by (shard % n_new).
 * - Scale-up ranks >= n_old: load nothing.
 */
static int
rescale_loader_rank_for_shard(int shard, int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    if (shard < 0 || shard >= n_old || n_new <= 0)
        FUNC_LEAVE(-1);

    if (shard < n_new)
        FUNC_LEAVE(shard);

    /* Orphan shard from a retired rank */
    FUNC_LEAVE(shard % n_new);
}

/*
 * True if this MPI rank should deserialize `shard` during elastic load (T7).
 */
static int
rescale_rank_should_load_shard(int shard, int n_old, int n_new, int my_rank)
{
    FUNC_ENTER(NULL);

    FUNC_LEAVE(rescale_loader_rank_for_shard(shard, n_old, n_new) == my_rank);
}

/*
 * Collectively ensure shards 0..n_old-1 exist under $PDC_TMPDIR.
 * On scale-down, log orphan-shard -> loader rank assignment.
 */
static perr_t
PDC_Server_rescale_validate_shards(int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;
    int    local_ok  = 1;
    int    global_ok = 1;
    int    s;
    char   path[ADDR_MAX];
    int    n_load = 0;

    if (strpbrk(pdc_server_tmp_dir_g, ";&|`$<>") != NULL)
        PGOTO_ERROR(FAIL, "Invalid characters in server tmp dir path");

    /* Shared filesystem: each rank verifies the full set; Allreduce MIN */
    for (s = 0; s < n_old; s++) {
        if (!checkpoint_shard_exists(s)) {
            local_ok = 0;
            break;
        }
    }

#ifdef ENABLE_MPI
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
    global_ok = local_ok;
#endif

    if (!global_ok) {
        if (pdc_server_rank_g == 0) {
            for (s = 0; s < n_old; s++) {
                if (!checkpoint_shard_exists(s)) {
                    checkpoint_shard_path(path, sizeof(path), s);
                    LOG_ERROR("Missing checkpoint shard for elastic restart: [%s] "
                              "(need all shards 0..%d)\n",
                              path, n_old - 1);
                    break;
                }
            }
        }
        PGOTO_ERROR(FAIL, "Elastic restart aborted: incomplete checkpoint shards for N_old=%d", n_old);
    }

    if (pdc_server_rank_g == 0) {
        LOG_INFO("Validated %d checkpoint shard(s) for elastic restart (%d -> %d)\n", n_old, n_old, n_new);
        if (n_old > n_new) {
            LOG_INFO("Scale-down orphan shard loaders (shard -> rank):\n");
            for (s = n_new; s < n_old; s++)
                LOG_INFO("  metadata_checkpoint.%d -> rank %d\n", s,
                         rescale_loader_rank_for_shard(s, n_old, n_new));
        }
        else if (n_old < n_new) {
            LOG_INFO("Scale-up: ranks %d..%d start with empty metadata holders\n", n_old, n_new - 1);
        }
    }

    /* Count shards this rank will load (for T7); log at debug/info on each rank */
    for (s = 0; s < n_old; s++) {
        if (rescale_rank_should_load_shard(s, n_old, n_new, pdc_server_rank_g))
            n_load++;
    }
    LOG_INFO("Rank %d will load %d checkpoint shard(s) during elastic restart\n", pdc_server_rank_g, n_load);

    rescale_validated_n_old_g = n_old;
    rescale_validated_n_new_g = n_new;

    /* Keep validated counts visible for follow-on load (T7) */
    if (is_debug_g && pdc_server_rank_g == 0)
        LOG_DEBUG("Cached elastic validation state N_old=%d N_new=%d\n", rescale_validated_n_old_g,
                  rescale_validated_n_new_g);

done:
    FUNC_LEAVE(ret_value);
}

/*
 * Elastic restart entry point. Later tasks fill in:
 *   T7 load, T8/T9 migrate, T10 barrier/server.cfg, ...
 */
perr_t
PDC_Server_restart_elastic(int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    perr_t ret_value = SUCCEED;

    if (n_old <= 0 || n_new <= 0)
        PGOTO_ERROR(FAIL, "Invalid server counts for elastic restart (n_old=%d, n_new=%d)", n_old, n_new);

    if (n_old == n_new)
        PGOTO_ERROR(FAIL,
                    "PDC_Server_restart_elastic called with n_old == n_new (%d); use PDC_Server_restart",
                    n_old);

    if (n_old > 65536 || n_new > 65536)
        PGOTO_ERROR(FAIL, "Server count out of range for elastic restart (n_old=%d, n_new=%d)", n_old, n_new);

    if (pdc_server_rank_g == 0)
        LOG_INFO("Elastic server restart requested: %d -> %d\n", n_old, n_new);

#ifndef ENABLE_MPI
    if (n_new != 1)
        PGOTO_ERROR(FAIL, "Elastic restart with n_new=%d requires MPI", n_new);
#endif

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ret_value = PDC_Server_rescale_validate_shards(n_old, n_new);
    if (ret_value != SUCCEED)
        PGOTO_ERROR(FAIL, "Elastic restart shard validation failed");

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /* Load/migration not implemented yet (T7–T10). */
    if (pdc_server_rank_g == 0)
        LOG_ERROR("Elastic metadata load/migration (%d -> %d) is not implemented yet "
                  "(shard validation succeeded)\n",
                  n_old, n_new);
    PGOTO_ERROR(FAIL, "Elastic server restart not implemented");

done:
    FUNC_LEAVE(ret_value);
}
