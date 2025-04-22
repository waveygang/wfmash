#ifndef FAIGZ_H
#define FAIGZ_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <pthread.h>
#include <inttypes.h>

#include "htslib/faidx.h"
#include "htslib/khash.h"

#ifdef __cplusplus
extern "C" {
#endif

// Forward declarations
typedef struct faidx_meta_t faidx_meta_t;
typedef struct faidx_reader_t faidx_reader_t;

// Key structures needed for our implementation
typedef struct {
    int id;
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;

// Simple string hash function
static khint_t faigz_str_hash_func(const char *s) {
    khint_t h = 0;
    while (*s) {
        h = (h << 5) - h + *s++;
    }
    return h;
}

// String hash type with unique macro names to avoid collision
#define kh_faigz_hash_func(key) faigz_str_hash_func(key)
#define kh_faigz_hash_equal(a, b) (strcmp((a), (b)) == 0)
KHASH_INIT(str, kh_cstr_t, faidx1_t, 1, kh_faigz_hash_func, kh_faigz_hash_equal)

// Shared metadata structure containing only the indices
struct faidx_meta_t {
    int n, m;                     // Sequence count and allocation size
    char **name;                  // Array of sequence names
    khash_t(str) *hash;          // Hash table mapping names to positions
    enum fai_format_options format; // FAI_FASTA or FAI_FASTQ

    // Source file paths
    char *fasta_path;            // Path to the FASTA/FASTQ file
    char *fai_path;              // Path to the .fai index
    char *gzi_path;              // Path to the .gzi index (if using BGZF)

    // Reference count and mutex for thread safety
    int ref_count;
    pthread_mutex_t mutex;

    // Flag indicating if the source is BGZF compressed
    int is_bgzf;
};

// Reader structure containing thread-specific data
struct faidx_reader_t {
    faidx_meta_t *meta;          // Shared metadata (not owned)
    faidx_t *fai;                // Thread-local faidx handle
};

/**
 * Load FASTA/FASTQ index metadata.
 *
 * @param filename Path to the FASTA/FASTQ file
 * @param format FAI_FASTA or FAI_FASTQ
 * @param flags Option flags (see FAI_CREATE in faidx.h)
 * @return Pointer to metadata or NULL on error
 */
faidx_meta_t *faidx_meta_load(const char *filename, enum fai_format_options format, int flags);

/**
 * Increment reference count on metadata
 *
 * @param meta Metadata to reference
 * @return The metadata
 */
faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta);

/**
 * Free shared metadata and its resources.
 * Decrements the reference count and frees only when it reaches zero.
 *
 * @param meta Metadata to free
 */
void faidx_meta_destroy(faidx_meta_t *meta);

/**
 * Create a reader from shared metadata
 *
 * @param meta Shared metadata (reference count is incremented)
 * @return New reader or NULL on error
 */
faidx_reader_t *faidx_reader_create(faidx_meta_t *meta);

/**
 * Destroy a reader.
 * This does not affect the shared metadata.
 *
 * @param reader Reader to destroy
 */
void faidx_reader_destroy(faidx_reader_t *reader);

/**
 * Fetch sequence from a specific region
 *
 * @param reader Reader to use
 * @param c_name Region name
 * @param p_beg_i Beginning position (0-based)
 * @param p_end_i End position (0-based)
 * @param len Output parameter for sequence length
 * @return Sequence string (must be freed by caller) or NULL on error
 */
char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                           hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/**
 * Fetch the quality string for a specific region (FASTQ only)
 *
 * @param reader Reader to use
 * @param c_name Region name
 * @param p_beg_i Beginning position (0-based)
 * @param p_end_i End position (0-based)
 * @param len Output parameter for string length
 * @return Quality string (must be freed by caller) or NULL on error
 */
char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/**
 * Get number of sequences in the index
 *
 * @param meta Metadata
 * @return Number of sequences
 */
int faidx_meta_nseq(const faidx_meta_t *meta);

/**
 * Get name of the i-th sequence
 *
 * @param meta Metadata
 * @param i Sequence index
 * @return Sequence name or NULL
 */
const char *faidx_meta_iseq(const faidx_meta_t *meta, int i);

/**
 * Get sequence length
 *
 * @param meta Metadata
 * @param seq Sequence name
 * @return Sequence length or -1 if not found
 */
hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq);

/**
 * Check if a sequence exists in the index
 *
 * @param meta Metadata
 * @param seq Sequence name
 * @return 1 if present, 0 if absent
 */
int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq);

/**
 * Parse a region string
 *
 * @param meta Metadata
 * @param s Region string
 * @param tid Output parameter for sequence ID
 * @param beg Output parameter for beginning position
 * @param end Output parameter for end position
 * @param flags Parsing flags
 * @return Pointer to end of parsed region or NULL on error
 */
const char *faidx_meta_parse_region(const faidx_meta_t *meta, const char *s,
                                  int *tid, hts_pos_t *beg, hts_pos_t *end,
                                  int flags);

#ifdef __cplusplus
}
#endif

/* Internal declarations for implementation */
#ifdef REENTRANT_FAIDX_IMPLEMENTATION

// Helper functions for internal use
static char *kstrdup(const char *str);
static int faidx_adjust_position(const faidx_meta_t *meta, int end_adjust,
                              faidx1_t *val_out, const char *c_name,
                              hts_pos_t *p_beg_i, hts_pos_t *p_end_i,
                              hts_pos_t *len);
static int fai_name2id(void *v, const char *ref);

// Implementation of helper functions

static char *kstrdup(const char *str) {
    int len = strlen(str) + 1;
    char *s = (char*)malloc(len);
    if (!s) return NULL;
    memcpy(s, str, len);
    return s;
}

static int fai_name2id(void *v, const char *ref) {
    faidx_meta_t *meta = (faidx_meta_t*)v;
    khint_t k = kh_get(str, meta->hash, ref);
    return k == kh_end(meta->hash) ? -1 : kh_val(meta->hash, k).id;
}

/* Load metadata from a FASTA/FASTQ file */
faidx_meta_t *faidx_meta_load(const char *filename, enum fai_format_options format, int flags) {
    kstring_t fai_kstr = {0}, gzi_kstr = {0};
    faidx_t *fai = NULL;
    faidx_meta_t *meta = NULL;
    FILE *fp = NULL;
    int is_bgzf = 0;

    /* Handle NULL filename */
    if (!filename) {
        errno = EINVAL;
        return NULL;
    }

    /* Construct file paths */
    if (ksprintf(&fai_kstr, "%s.fai", filename) < 0) goto fail;
    if (ksprintf(&gzi_kstr, "%s.gzi", filename) < 0) goto fail;

    /* Check if file is BGZF compressed */
    fp = fopen(filename, "rb");
    if (fp) {
        unsigned char magic[2];
        if (fread(magic, 1, 2, fp) == 2) {
            is_bgzf = (magic[0] == 0x1f && magic[1] == 0x8b);
        }
        fclose(fp);
        fp = NULL;
    }

    /* Load the FASTA/FASTQ index using fai_load3_format which handles either FASTA or FASTQ */
    fai = fai_load3_format(filename, fai_kstr.s, gzi_kstr.s, flags, format);
    if (!fai) goto fail;

    /* Create the metadata structure */
    meta = (faidx_meta_t*)calloc(1, sizeof(faidx_meta_t));
    if (!meta) goto fail;

    /* Initialize the mutex */
    if (pthread_mutex_init(&meta->mutex, NULL) != 0) goto fail;

    /* Get basic information from faidx */
    meta->n = faidx_nseq(fai);
    meta->m = meta->n;  /* We'll allocate exactly what we need */
    meta->format = format;
    meta->ref_count = 1;
    meta->is_bgzf = is_bgzf;

    /* Copy sequence names */
    meta->name = (char**)malloc(meta->m * sizeof(char*));
    if (!meta->name) goto fail;

    /* Create hash table */
    meta->hash = kh_init(str);
    if (!meta->hash) goto fail;

    /* Populate names and hash table with sequence info */
    for (int i = 0; i < meta->n; i++) {
        const char *seq_name = faidx_iseq(fai, i);
        meta->name[i] = kstrdup(seq_name);
        if (!meta->name[i]) goto fail;

        faidx1_t val;
        val.id = i;
        val.len = faidx_seq_len(fai, seq_name);

        /* Get line length - this information is not directly accessible from API
           Use reasonable defaults if needed */
        val.line_len = 60 + 1;  /* Default line length + newline */
        val.line_blen = 60;     /* Default bases per line */

        /* Estimated offsets - not accurate but we don't need them since we use faidx API */
        val.seq_offset = 0;
        val.qual_offset = 0;

        /* Build our own lookup table */
        int absent;
        khint_t k = kh_put(str, meta->hash, meta->name[i], &absent);
        if (absent) {
            kh_val(meta->hash, k) = val;
        }
    }

    /* Store file paths */
    meta->fasta_path = kstrdup(filename);
    meta->fai_path = kstrdup(fai_kstr.s);
    meta->gzi_path = kstrdup(gzi_kstr.s);
    if (!meta->fasta_path || !meta->fai_path || !meta->gzi_path) goto fail;

    /* Clean up */
    free(fai_kstr.s);
    free(gzi_kstr.s);
    fai_destroy(fai);

    return meta;

fail:
    if (meta) {
        if (meta->hash) kh_destroy(str, meta->hash);
        if (meta->name) {
            for (int i = 0; i < meta->n; i++) {
                if (meta->name[i]) free(meta->name[i]);
            }
            free(meta->name);
        }
        free(meta->fasta_path);
        free(meta->fai_path);
        free(meta->gzi_path);
        pthread_mutex_destroy(&meta->mutex);
        free(meta);
    }
    free(fai_kstr.s);
    free(gzi_kstr.s);
    if (fai) fai_destroy(fai);
    return NULL;
}

/* Reference counting for metadata */
faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta) {
    if (!meta) return NULL;

    pthread_mutex_lock(&meta->mutex);
    meta->ref_count++;
    pthread_mutex_unlock(&meta->mutex);

    return meta;
}

/* Destroy metadata */
void faidx_meta_destroy(faidx_meta_t *meta) {
    if (!meta) return;

    pthread_mutex_lock(&meta->mutex);
    meta->ref_count--;
    int should_free = (meta->ref_count <= 0);
    pthread_mutex_unlock(&meta->mutex);

    if (should_free) {
        /* Free all resources */
        if (meta->hash) {
            kh_destroy(str, meta->hash);
        }

        if (meta->name) {
            for (int i = 0; i < meta->n; i++) {
                free(meta->name[i]);
            }
            free(meta->name);
        }

        free(meta->fasta_path);
        free(meta->fai_path);
        free(meta->gzi_path);

        pthread_mutex_destroy(&meta->mutex);
        free(meta);
    }
}

/* Create a reader */
faidx_reader_t *faidx_reader_create(faidx_meta_t *meta) {
    if (!meta) return NULL;

    faidx_reader_t *reader = (faidx_reader_t*)calloc(1, sizeof(faidx_reader_t));
    if (!reader) return NULL;

    /* Reference the metadata */
    reader->meta = faidx_meta_ref(meta);

    /* Open the file using fai_load3_format which handles BGZF automatically */
    reader->fai = fai_load3_format(meta->fasta_path, meta->fai_path, meta->gzi_path, 0, meta->format);

    if (!reader->fai) {
        faidx_meta_destroy(reader->meta);
        free(reader);
        return NULL;
    }

    return reader;
}

/* Destroy a reader */
void faidx_reader_destroy(faidx_reader_t *reader) {
    if (!reader) return;

    if (reader->fai) {
        fai_destroy(reader->fai);
    }

    faidx_meta_destroy(reader->meta);
    free(reader);
}

/* Helper: Adjust position to sequence boundaries */
static int faidx_adjust_position(const faidx_meta_t *meta, int end_adjust,
                              faidx1_t *val_out, const char *c_name,
                              hts_pos_t *p_beg_i, hts_pos_t *p_end_i,
                              hts_pos_t *len) {
    khiter_t iter;
    faidx1_t *val;

    /* Adjust position */
    iter = kh_get(str, meta->hash, c_name);

    if (iter == kh_end(meta->hash)) {
        if (len) *len = -2;
        return 1;
    }

    val = &kh_val(meta->hash, iter);

    if (val_out) *val_out = *val;

    if (*p_end_i < *p_beg_i) *p_beg_i = *p_end_i;

    if (*p_beg_i < 0) *p_beg_i = 0;
    else if ((hts_pos_t)val->len <= *p_beg_i) *p_beg_i = val->len;

    if (*p_end_i < 0) *p_end_i = 0;
    else if ((hts_pos_t)val->len <= *p_end_i) *p_end_i = val->len - end_adjust;

    return 0;
}

/* Fetch sequence */
char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                          hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    faidx1_t val;
    int len_int = 0;

    /* Adjust position */
    if (faidx_adjust_position(reader->meta, 1, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    /* Use the standard faidx API to fetch the sequence */
    char *result = faidx_fetch_seq(reader->fai, c_name, (int)p_beg_i, (int)p_end_i, &len_int);

    if (len) *len = len_int;
    return result;
}

/* Fetch quality string */
char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    faidx1_t val;
    int len_int = 0;

    if (reader->meta->format != FAI_FASTQ) {
        if (len) *len = -2;
        return NULL;
    }

    /* Adjust position */
    if (faidx_adjust_position(reader->meta, 1, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    /* Use the standard faidx API to fetch the quality */
    char *result = faidx_fetch_qual(reader->fai, c_name, (int)p_beg_i, (int)p_end_i, &len_int);

    if (len) *len = len_int;
    return result;
}

/* Get number of sequences */
int faidx_meta_nseq(const faidx_meta_t *meta) {
    return meta ? meta->n : 0;
}

/* Get sequence name */
const char *faidx_meta_iseq(const faidx_meta_t *meta, int i) {
    return (meta && i >= 0 && i < meta->n) ? meta->name[i] : NULL;
}

/* Get sequence length */
hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return -1;

    khint_t k = kh_get(str, meta->hash, seq);
    if (k == kh_end(meta->hash)) return -1;

    return kh_val(meta->hash, k).len;
}

/* Check if sequence exists */
int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return 0;

    khint_t k = kh_get(str, meta->hash, seq);
    return (k != kh_end(meta->hash));
}

/* Parse a region string */
const char *faidx_meta_parse_region(const faidx_meta_t *meta, const char *s,
                                 int *tid, hts_pos_t *beg, hts_pos_t *end,
                                 int flags) {
    return hts_parse_region(s, tid, beg, end, fai_name2id, (void*)meta, flags);
}

#endif /* REENTRANT_FAIDX_IMPLEMENTATION */

#endif /* FAIGZ_H */
