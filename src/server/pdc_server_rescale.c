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
#include "pdc_utlist.h"
#include "pdc_malloc.h"
#include "pdc_server_rescale.h"
#include "pdc_server_metadata.h"
#include "pdc_client_server_common.h"
#include "pdc_logger.h"
#include "bulki.h"
#include "bulki_serde.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

/* Keep in sync with pdc_server.c checkpoint bounds / magic */
#define PDC_CHECKPOINT_MAGIC_CURRENT       "PDC26.03"
#define PDC_CHECKPOINT_MAX_METADATA_COUNT  10000000
#define PDC_CHECKPOINT_MAX_KVTAG_KEY_LEN   65536
#define PDC_CHECKPOINT_MAX_KVTAG_SIZE      (1u << 24)
#define PDC_CHECKPOINT_MAX_REGION_COUNT    1000000
#define PDC_CHECKPOINT_MAX_HIST_NBIN       65536

typedef struct pdc_rescale_tmp_cont {
    uint32_t                         hash_key;
    int                              source_shard;
    pdc_cont_hash_table_entry_t *    entry;
    struct pdc_rescale_tmp_cont *    next;
} pdc_rescale_tmp_cont_t;

/* Cached after successful T6 validation for T7 load planning */
static int rescale_validated_n_old_g = -1;
static int rescale_validated_n_new_g = -1;

/* Temporary holders (T7) — not published to hash tables until migrate (T8/T9) */
static pdc_metadata_t *        rescale_tmp_objs_g     = NULL;
static int                     rescale_tmp_n_obj_g    = 0;
static int                     rescale_tmp_n_region_g = 0;
static pdc_rescale_tmp_cont_t *rescale_tmp_conts_g    = NULL;
static int                     rescale_tmp_n_cont_g   = 0;

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

static int
rescale_loader_rank_for_shard(int shard, int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    if (shard < 0 || shard >= n_old || n_new <= 0)
        FUNC_LEAVE(-1);

    if (shard < n_new)
        FUNC_LEAVE(shard);

    FUNC_LEAVE(shard % n_new);
}

static int
rescale_rank_should_load_shard(int shard, int n_old, int n_new, int my_rank)
{
    FUNC_ENTER(NULL);

    FUNC_LEAVE(rescale_loader_rank_for_shard(shard, n_old, n_new) == my_rank);
}

static int
rescale_region_cmp(region_list_t *a, region_list_t *b)
{
    FUNC_ENTER(NULL);

    size_t unit_size = sizeof(uint64_t) * a->ndim;

    FUNC_LEAVE(memcmp(a->start, b->start, unit_size));
}

static void
rescale_free_kvtag_list(pdc_kvtag_list_t *head)
{
    FUNC_ENTER(NULL);

    pdc_kvtag_list_t *elt, *tmp;

    DL_FOREACH_SAFE(head, elt, tmp)
    {
        DL_DELETE(head, elt);
        if (elt->kvtag) {
            elt->kvtag->name  = (char *)PDC_free(elt->kvtag->name);
            elt->kvtag->value = PDC_free(elt->kvtag->value);
            elt->kvtag        = (pdc_kvtag_t *)PDC_free(elt->kvtag);
        }
        elt = (pdc_kvtag_list_t *)PDC_free(elt);
    }

    FUNC_LEAVE_VOID();
}

static void
rescale_free_region_list(region_list_t *head)
{
    FUNC_ENTER(NULL);

    region_list_t *elt, *tmp;

    DL_FOREACH_SAFE(head, elt, tmp)
    {
        DL_DELETE(head, elt);
        if (elt->region_hist) {
            elt->region_hist->range = (double *)PDC_free(elt->region_hist->range);
            elt->region_hist->bin   = (uint64_t *)PDC_free(elt->region_hist->bin);
            elt->region_hist        = (pdc_histogram_t *)PDC_free(elt->region_hist);
        }
        elt = (region_list_t *)PDC_free(elt);
    }

    FUNC_LEAVE_VOID();
}

static void
PDC_Server_rescale_free_tmp_holders(void)
{
    FUNC_ENTER(NULL);

    pdc_metadata_t *       obj, *obj_tmp;
    pdc_rescale_tmp_cont_t *cont, *cont_tmp;

    DL_FOREACH_SAFE(rescale_tmp_objs_g, obj, obj_tmp)
    {
        DL_DELETE(rescale_tmp_objs_g, obj);
        rescale_free_kvtag_list(obj->kvtag_list_head);
        rescale_free_region_list(obj->storage_region_list_head);
        obj = (pdc_metadata_t *)PDC_free(obj);
    }
    rescale_tmp_objs_g     = NULL;
    rescale_tmp_n_obj_g    = 0;
    rescale_tmp_n_region_g = 0;

    cont = rescale_tmp_conts_g;
    while (cont) {
        cont_tmp = cont->next;
        if (cont->entry) {
            cont->entry->obj_ids = (uint64_t *)PDC_free(cont->entry->obj_ids);
            cont->entry          = (pdc_cont_hash_table_entry_t *)PDC_free(cont->entry);
        }
        cont = (pdc_rescale_tmp_cont_t *)PDC_free(cont);
        cont = cont_tmp;
    }
    rescale_tmp_conts_g  = NULL;
    rescale_tmp_n_cont_g = 0;

    FUNC_LEAVE_VOID();
}

static perr_t
rescale_append_region_from_bulki(BULKI *region_entry, pdc_metadata_t *metadata, int *n_region_out)
{
    FUNC_ENTER(NULL);

    perr_t         ret_value  = SUCCEED;
    region_list_t *region_list = NULL;
    unsigned       idx;
    int            has_hist;

    region_list = (region_list_t *)PDC_malloc(sizeof(region_list_t));
    if (region_list == NULL)
        PGOTO_ERROR(FAIL, "Cannot allocate region_list");

    BULKI_Entity *region_ent = BULKI_get(region_entry, BULKI_singleton_ENTITY("region", PDC_STRING));
    if (region_ent == NULL || region_ent->data == NULL)
        PGOTO_ERROR(FAIL, "Missing region in checkpoint");
    memcpy(region_list, region_ent->data, sizeof(region_list_t));
    region_list->region_hist = NULL;

    BULKI_Entity *has_hist_ent = BULKI_get(region_entry, BULKI_singleton_ENTITY("has_hist", PDC_STRING));
    if (has_hist_ent == NULL)
        PGOTO_ERROR(FAIL, "Missing has_hist in checkpoint");
    memcpy(&has_hist, has_hist_ent->data, sizeof(int));

    if (has_hist == 1) {
        BULKI_Entity *histogram_ent =
            BULKI_get(region_entry, BULKI_singleton_ENTITY("histogram", PDC_STRING));
        if (histogram_ent != NULL && histogram_ent->pdc_type == PDC_BULKI) {
            BULKI *histogram = (BULKI *)histogram_ent->data;
            int    nbin;

            region_list->region_hist = (pdc_histogram_t *)PDC_malloc(sizeof(pdc_histogram_t));
            if (region_list->region_hist == NULL)
                PGOTO_ERROR(FAIL, "Cannot allocate histogram");

            BULKI_Entity *dtype_ent = BULKI_get(histogram, BULKI_singleton_ENTITY("dtype", PDC_STRING));
            if (dtype_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing histogram dtype");
            memcpy(&region_list->region_hist->dtype, dtype_ent->data, sizeof(int));

            BULKI_Entity *nbin_ent = BULKI_get(histogram, BULKI_singleton_ENTITY("nbin", PDC_STRING));
            if (nbin_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing histogram nbin");
            memcpy(&nbin, nbin_ent->data, sizeof(int));
            if (nbin <= 0 || nbin > PDC_CHECKPOINT_MAX_HIST_NBIN)
                PGOTO_ERROR(FAIL, "Invalid histogram nbin %d", nbin);
            region_list->region_hist->nbin = nbin;

            BULKI_Entity *range_ent = BULKI_get(histogram, BULKI_singleton_ENTITY("range", PDC_STRING));
            if (range_ent == NULL || range_ent->data == NULL)
                PGOTO_ERROR(FAIL, "Missing histogram range");
            region_list->region_hist->range = (double *)PDC_malloc(sizeof(double) * (size_t)nbin * 2);
            memcpy(region_list->region_hist->range, range_ent->data, sizeof(double) * (size_t)nbin * 2);

            BULKI_Entity *bin_ent = BULKI_get(histogram, BULKI_singleton_ENTITY("bin", PDC_STRING));
            if (bin_ent == NULL || bin_ent->data == NULL)
                PGOTO_ERROR(FAIL, "Missing histogram bin");
            region_list->region_hist->bin = (uint64_t *)PDC_malloc(sizeof(uint64_t) * (size_t)nbin);
            memcpy(region_list->region_hist->bin, bin_ent->data, sizeof(uint64_t) * (size_t)nbin);

            BULKI_Entity *incr_ent = BULKI_get(histogram, BULKI_singleton_ENTITY("incr", PDC_STRING));
            if (incr_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing histogram incr");
            memcpy(&region_list->region_hist->incr, incr_ent->data, sizeof(double));
        }
    }

    region_list->buf       = NULL;
    region_list->data_size = 1;
    for (idx = 0; idx < region_list->ndim; idx++)
        region_list->data_size *= region_list->count[idx];
    region_list->is_data_ready            = 0;
    region_list->shm_fd                   = 0;
    region_list->meta                     = metadata;
    region_list->prev                     = NULL;
    region_list->next                     = NULL;
    region_list->overlap_storage_regions  = NULL;
    region_list->n_overlap_storage_region = 0;
    hg_atomic_init32(&(region_list->buf_map_refcount), 0);
    region_list->reg_dirty_from_buf = 0;
    region_list->access_type        = PDC_NA;
    region_list->bulk_handle        = NULL;
    region_list->lock_handle        = NULL;
    region_list->addr               = NULL;
    region_list->obj_id             = metadata->obj_id;
    region_list->reg_id             = 0;
    region_list->from_obj_id        = 0;
    region_list->client_id          = 0;
    region_list->is_io_done         = 0;
    region_list->is_shm_closed      = 0;
    region_list->seq_id             = 0;
    region_list->sent_to_server     = 0;
    region_list->io_cache_region    = NULL;
    memset(region_list->shm_addr, 0, ADDR_MAX);
    memset(region_list->client_ids, 0, PDC_SERVER_MAX_PROC_PER_NODE * sizeof(uint32_t));

    if (strstr(region_list->storage_location, "/global/cscratch") != NULL)
        region_list->data_loc_type = PDC_LUSTRE;

    DL_APPEND(metadata->storage_region_list_head, region_list);
    if (n_region_out)
        (*n_region_out)++;
    region_list = NULL;

done:
    if (ret_value != SUCCEED && region_list != NULL) {
        if (region_list->region_hist) {
            region_list->region_hist->range = (double *)PDC_free(region_list->region_hist->range);
            region_list->region_hist->bin   = (uint64_t *)PDC_free(region_list->region_hist->bin);
            region_list->region_hist        = (pdc_histogram_t *)PDC_free(region_list->region_hist);
        }
        region_list = (region_list_t *)PDC_free(region_list);
    }
    FUNC_LEAVE(ret_value);
}

static perr_t
rescale_append_obj_from_bulki(BULKI *metadata_obj, int source_shard, int *n_region_out)
{
    FUNC_ENTER(NULL);

    perr_t          ret_value = SUCCEED;
    pdc_metadata_t *metadata  = NULL;

    metadata = (pdc_metadata_t *)PDC_calloc(1, sizeof(pdc_metadata_t));
    if (metadata == NULL)
        PGOTO_ERROR(FAIL, "Cannot allocate metadata");

    BULKI_Entity *metadata_ent = BULKI_get(metadata_obj, BULKI_singleton_ENTITY("metadata", PDC_STRING));
    if (metadata_ent == NULL || metadata_ent->data == NULL)
        PGOTO_ERROR(FAIL, "Missing metadata in checkpoint");
    memcpy(metadata, metadata_ent->data, sizeof(pdc_metadata_t));

    metadata->storage_region_list_head       = NULL;
    metadata->region_lock_head               = NULL;
    metadata->region_map_head                = NULL;
    metadata->region_buf_map_head            = NULL;
    metadata->bloom                          = NULL;
    metadata->prev                           = NULL;
    metadata->next                           = NULL;
    metadata->kvtag_list_head                = NULL;
    metadata->all_storage_region_distributed = 0;
    (void)source_shard;

    BULKI_Entity *kvtags_array = BULKI_get(metadata_obj, BULKI_singleton_ENTITY("kvtags", PDC_STRING));
    if (kvtags_array != NULL && kvtags_array->pdc_type == PDC_BULKI) {
        BULKI_Entity_Iterator *kvtag_iter = Bent_iterator_init(kvtags_array, NULL, PDC_BULKI);
        while (Bent_iterator_has_next_BULKI(kvtag_iter)) {
            BULKI *           kvtag_entry = Bent_iterator_next_BULKI(kvtag_iter);
            pdc_kvtag_list_t *kvtag_list  = (pdc_kvtag_list_t *)PDC_calloc(1, sizeof(pdc_kvtag_list_t));
            uint32_t          kv_size;

            kvtag_list->kvtag = (pdc_kvtag_t *)PDC_malloc(sizeof(pdc_kvtag_t));

            BULKI_Entity *key_ent = BULKI_get(kvtag_entry, BULKI_singleton_ENTITY("key", PDC_STRING));
            if (key_ent == NULL || key_ent->data == NULL)
                PGOTO_ERROR(FAIL, "Invalid kvtag key in checkpoint");
            {
                int key_len = (int)strlen((char *)key_ent->data) + 1;
                if (key_len <= 0 || key_len > PDC_CHECKPOINT_MAX_KVTAG_KEY_LEN)
                    PGOTO_ERROR(FAIL, "Invalid key_len %d in checkpoint", key_len);
                kvtag_list->kvtag->name = PDC_malloc((size_t)key_len);
                memcpy(kvtag_list->kvtag->name, key_ent->data, (size_t)key_len);
            }

            BULKI_Entity *size_ent = BULKI_get(kvtag_entry, BULKI_singleton_ENTITY("size", PDC_STRING));
            if (size_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing kvtag size");
            memcpy(&kv_size, size_ent->data, sizeof(uint32_t));
            if (kv_size == 0 || kv_size > PDC_CHECKPOINT_MAX_KVTAG_SIZE)
                PGOTO_ERROR(FAIL, "Invalid kvtag size");
            kvtag_list->kvtag->size = kv_size;

            BULKI_Entity *type_ent = BULKI_get(kvtag_entry, BULKI_singleton_ENTITY("type", PDC_STRING));
            if (type_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing kvtag type");
            memcpy(&kvtag_list->kvtag->type, type_ent->data, sizeof(int8_t));

            BULKI_Entity *value_ent = BULKI_get(kvtag_entry, BULKI_singleton_ENTITY("value", PDC_STRING));
            if (value_ent == NULL || value_ent->count != kv_size)
                PGOTO_ERROR(FAIL, "Invalid kvtag value size");
            kvtag_list->kvtag->value = PDC_malloc((size_t)kv_size);
            memcpy(kvtag_list->kvtag->value, value_ent->data, (size_t)kv_size);

            DL_APPEND(metadata->kvtag_list_head, kvtag_list);
        }
    }

    BULKI_Entity *regions_array = BULKI_get(metadata_obj, BULKI_singleton_ENTITY("regions", PDC_STRING));
    if (regions_array != NULL && regions_array->pdc_type == PDC_BULKI) {
        int n_region = regions_array->count;
        if (n_region < 0 || n_region > PDC_CHECKPOINT_MAX_REGION_COUNT)
            PGOTO_ERROR(FAIL, "Suspicious n_region value %d", n_region);

        BULKI_Entity_Iterator *region_iter = Bent_iterator_init(regions_array, NULL, PDC_BULKI);
        while (Bent_iterator_has_next_BULKI(region_iter)) {
            BULKI *region_entry = Bent_iterator_next_BULKI(region_iter);
            if (rescale_append_region_from_bulki(region_entry, metadata, n_region_out) != SUCCEED)
                PGOTO_ERROR(FAIL, "Failed to restore region into temp holder");
        }
        DL_SORT(metadata->storage_region_list_head, rescale_region_cmp);
    }

    DL_APPEND(rescale_tmp_objs_g, metadata);
    rescale_tmp_n_obj_g++;
    metadata = NULL;

done:
    if (ret_value != SUCCEED && metadata != NULL) {
        rescale_free_kvtag_list(metadata->kvtag_list_head);
        rescale_free_region_list(metadata->storage_region_list_head);
        metadata = (pdc_metadata_t *)PDC_free(metadata);
    }
    FUNC_LEAVE(ret_value);
}

static perr_t
rescale_load_one_shard(int shard_rank)
{
    FUNC_ENTER(NULL);

    perr_t        ret_value        = SUCCEED;
    FILE *        file             = NULL;
    BULKI *       checkpoint_bulki = NULL;
    char          path[ADDR_MAX];
    int           n_region_shard = 0;

    checkpoint_shard_path(path, sizeof(path), shard_rank);
    file = fopen(path, "rb");
    if (file == NULL)
        PGOTO_ERROR(FAIL, "Error with fopen, filename: [%s]", path);

    checkpoint_bulki = BULKI_deserialize_from_file(file);
    file             = NULL;
    if (checkpoint_bulki == NULL)
        PGOTO_ERROR(FAIL, "Failed to deserialize checkpoint [%s]", path);

    BULKI_Entity *version_entity =
        BULKI_get(checkpoint_bulki, BULKI_singleton_ENTITY("version_number", PDC_STRING));
    if (version_entity == NULL)
        PGOTO_ERROR(FAIL, "Missing version_number in [%s]", path);
    {
        BULKI_Entity *expected = BULKI_ENTITY(PDC_CHECKPOINT_MAGIC_CURRENT, 1, PDC_STRING, PDC_CLS_ITEM);
        int           equal    = BULKI_Entity_equal(version_entity, expected);
        BULKI_Entity_free(expected, 1);
        if (!equal)
            PGOTO_ERROR(FAIL, "Checkpoint version mismatch in [%s]", path);
    }

    BULKI_Entity *containers_array =
        BULKI_get(checkpoint_bulki, BULKI_singleton_ENTITY("containers", PDC_STRING));
    if (containers_array != NULL && containers_array->pdc_type == PDC_BULKI) {
        BULKI_Entity_Iterator *cont_iter = Bent_iterator_init(containers_array, NULL, PDC_BULKI);
        while (Bent_iterator_has_next_BULKI(cont_iter)) {
            BULKI *                 container_entry = Bent_iterator_next_BULKI(cont_iter);
            pdc_rescale_tmp_cont_t *tmp_cont        = NULL;
            pdc_cont_hash_table_entry_t *cont_entry = NULL;

            BULKI_Entity *hash_key_ent =
                BULKI_get(container_entry, BULKI_singleton_ENTITY("hash_key", PDC_STRING));
            if (hash_key_ent == NULL || hash_key_ent->data == NULL)
                PGOTO_ERROR(FAIL, "Missing container hash_key in [%s]", path);

            BULKI_Entity *cont_data_ent =
                BULKI_get(container_entry, BULKI_singleton_ENTITY("cont_data", PDC_STRING));
            if (cont_data_ent == NULL || cont_data_ent->data == NULL)
                PGOTO_ERROR(FAIL, "Missing container cont_data in [%s]", path);

            cont_entry = (pdc_cont_hash_table_entry_t *)PDC_malloc(sizeof(pdc_cont_hash_table_entry_t));
            memcpy(cont_entry, cont_data_ent->data, sizeof(pdc_cont_hash_table_entry_t));
            /* Pointers inside the blob are not valid across restart */
            cont_entry->obj_ids        = NULL;
            cont_entry->n_obj          = 0;
            cont_entry->n_allocated    = 0;
            cont_entry->kvtag_list_head = NULL;

            tmp_cont = (pdc_rescale_tmp_cont_t *)PDC_calloc(1, sizeof(pdc_rescale_tmp_cont_t));
            memcpy(&tmp_cont->hash_key, hash_key_ent->data, sizeof(uint32_t));
            tmp_cont->source_shard = shard_rank;
            tmp_cont->entry        = cont_entry;
            tmp_cont->next         = rescale_tmp_conts_g;
            rescale_tmp_conts_g    = tmp_cont;
            rescale_tmp_n_cont_g++;
        }
    }

    BULKI_Entity *metadata_entries_array =
        BULKI_get(checkpoint_bulki, BULKI_singleton_ENTITY("metadata_entries", PDC_STRING));
    if (metadata_entries_array != NULL && metadata_entries_array->pdc_type == PDC_BULKI) {
        BULKI_Entity_Iterator *entry_iter = Bent_iterator_init(metadata_entries_array, NULL, PDC_BULKI);
        while (Bent_iterator_has_next_BULKI(entry_iter)) {
            BULKI *hash_entry = Bent_iterator_next_BULKI(entry_iter);
            int    count;

            BULKI_Entity *n_obj_ent = BULKI_get(hash_entry, BULKI_singleton_ENTITY("n_obj", PDC_STRING));
            if (n_obj_ent == NULL)
                PGOTO_ERROR(FAIL, "Missing n_obj in [%s]", path);
            memcpy(&count, n_obj_ent->data, sizeof(int));
            if (count <= 0 || count > PDC_CHECKPOINT_MAX_METADATA_COUNT)
                PGOTO_ERROR(FAIL, "Suspicious count value %d in [%s]", count, path);

            BULKI_Entity *metadata_objs_array =
                BULKI_get(hash_entry, BULKI_singleton_ENTITY("metadata_objects", PDC_STRING));
            if (metadata_objs_array == NULL || metadata_objs_array->pdc_type != PDC_BULKI)
                PGOTO_ERROR(FAIL, "Missing metadata_objects in [%s]", path);

            BULKI_Entity_Iterator *obj_iter = Bent_iterator_init(metadata_objs_array, NULL, PDC_BULKI);
            int                    i        = 0;
            while (Bent_iterator_has_next_BULKI(obj_iter) && i < count) {
                BULKI *metadata_obj = Bent_iterator_next_BULKI(obj_iter);
                if (rescale_append_obj_from_bulki(metadata_obj, shard_rank, &n_region_shard) != SUCCEED)
                    PGOTO_ERROR(FAIL, "Failed to load object from shard %d", shard_rank);
                i++;
            }
            if (i != count)
                PGOTO_ERROR(FAIL, "metadata object count mismatch in [%s] (got %d expect %d)", path, i,
                            count);
        }
    }

    rescale_tmp_n_region_g += n_region_shard;
    LOG_INFO("Rank %d loaded shard %d into temp holders (%d objs, %d regions so far; %d conts)\n",
             pdc_server_rank_g, shard_rank, rescale_tmp_n_obj_g, rescale_tmp_n_region_g,
             rescale_tmp_n_cont_g);

done:
    if (checkpoint_bulki != NULL)
        BULKI_free(checkpoint_bulki, 1);
    FUNC_LEAVE(ret_value);
}

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

    for (s = 0; s < n_old; s++) {
        if (rescale_rank_should_load_shard(s, n_old, n_new, pdc_server_rank_g))
            n_load++;
    }
    LOG_INFO("Rank %d will load %d checkpoint shard(s) during elastic restart\n", pdc_server_rank_g, n_load);

    rescale_validated_n_old_g = n_old;
    rescale_validated_n_new_g = n_new;

done:
    FUNC_LEAVE(ret_value);
}

/*
 * Deserialize assigned N_old shards into temporary holders (not final hash tables).
 */
static perr_t
PDC_Server_rescale_load_shards(int n_old, int n_new)
{
    FUNC_ENTER(NULL);

    perr_t ret_value   = SUCCEED;
    int    s;
    int    local_nobj  = 0;
    int    all_nobj    = 0;
    int    local_ncont = 0;
    int    all_ncont   = 0;

    PDC_Server_rescale_free_tmp_holders();

    for (s = 0; s < n_old; s++) {
        if (!rescale_rank_should_load_shard(s, n_old, n_new, pdc_server_rank_g))
            continue;
        if (rescale_load_one_shard(s) != SUCCEED)
            PGOTO_ERROR(FAIL, "Failed loading checkpoint shard %d on rank %d", s, pdc_server_rank_g);
    }

    local_nobj  = rescale_tmp_n_obj_g;
    local_ncont = rescale_tmp_n_cont_g;

#ifdef ENABLE_MPI
    MPI_Allreduce(&local_nobj, &all_nobj, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_ncont, &all_ncont, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    all_nobj  = local_nobj;
    all_ncont = local_ncont;
#endif

    LOG_INFO("Rank %d temp holders: %d objects, %d regions, %d containers\n", pdc_server_rank_g, local_nobj,
             rescale_tmp_n_region_g, local_ncont);

    if (pdc_server_rank_g == 0)
        LOG_INFO("Elastic load complete: %d objects and %d containers across all ranks "
                 "(N_old=%d -> N_new=%d); migration not yet applied\n",
                 all_nobj, all_ncont, n_old, n_new);

    (void)rescale_validated_n_old_g;
    (void)rescale_validated_n_new_g;

done:
    if (ret_value != SUCCEED)
        PDC_Server_rescale_free_tmp_holders();
    FUNC_LEAVE(ret_value);
}

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

    ret_value = PDC_Server_rescale_load_shards(n_old, n_new);
    if (ret_value != SUCCEED)
        PGOTO_ERROR(FAIL, "Elastic restart shard load failed");

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /* Migration not implemented yet (T8–T10). Holders remain until free on process exit/error path. */
    if (pdc_server_rank_g == 0)
        LOG_ERROR("Elastic metadata migration (%d -> %d) is not implemented yet "
                  "(shard validation and temp load succeeded)\n",
                  n_old, n_new);
    PDC_Server_rescale_free_tmp_holders();
    PGOTO_ERROR(FAIL, "Elastic server restart not implemented");

done:
    FUNC_LEAVE(ret_value);
}
