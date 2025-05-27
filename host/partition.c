#include <common.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dpu.h>
#include <dpu_types.h>

Graph *global_g;
bitmap_t bitmap;
double workload[N];

static int deg_cmp(const void *a, const void *b) {
    node_t x = *(node_t *)a;
    node_t y = *(node_t *)b;
    return global_g->row_ptr[y + 1] - global_g->row_ptr[y] - (global_g->row_ptr[x + 1] - global_g->row_ptr[x]);
}
static int node_t_cmp(const void *a, const void *b) {
    node_t x = *(node_t *)a;
    node_t y = *(node_t *)b;
    return x - y;
}
static int workload_cmp(const void *a, const void *b) {
    node_t x = *(node_t *)a;
    node_t y = *(node_t *)b;
    return workload[y] - workload[x];
}

#ifdef MORE_ACCURATE_MODEL
static inline double predict_workload(Graph *g, node_t root) {
    double deg = g->row_ptr[root + 1] - g->row_ptr[root];
    edge_ptr l = g->row_ptr[root], r = g->row_ptr[root + 1];
    while (l < r) {
        edge_ptr mid = (l + r) >> 1;
        if (g->col_idx[mid] < root) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    double eff_deg = eff_deg = l - g->row_ptr[root];
    if (deg > MRAM_BUF_SIZE) {
        printf(ANSI_COLOR_RED "Error: deg too large\n" ANSI_COLOR_RESET);
        exit(1);
    }
    if (eff_deg > BITMAP_SIZE * 32) {
        printf(ANSI_COLOR_RED "Error: eff_deg too large\n" ANSI_COLOR_RESET);
        exit(1);
    }
    double avg_deg = 0;
    for (edge_ptr i = g->row_ptr[root]; i < g->row_ptr[root + 1]; i++) {
        node_t neighbor = g->col_idx[i];
        avg_deg += g->row_ptr[neighbor + 1] - g->row_ptr[neighbor];
    }
    avg_deg /= g->row_ptr[root + 1] - g->row_ptr[root];
    double n = global_g->n;
    (void)deg;
    (void)eff_deg;
    (void)avg_deg;
    (void)n;
#if defined(CLIQUE2)
    return eff_deg;
#elif defined(CLIQUE3)
    return eff_deg * eff_deg * avg_deg + 100;
#elif defined(CLIQUE4)
    return eff_deg * eff_deg * eff_deg * avg_deg * avg_deg * avg_deg + 100;
#elif defined(CLIQUE5)
    return eff_deg * eff_deg * eff_deg * eff_deg * avg_deg * avg_deg * avg_deg * avg_deg * avg_deg * avg_deg + 100;
#elif defined(CYCLE4)
    return eff_deg * eff_deg * avg_deg + 100;
#elif defined(HOUSE5)
    return eff_deg * deg * avg_deg * (2 + deg / n + avg_deg / n) + 100;
#elif defined(TRI_TRI6)
    return eff_deg * eff_deg * avg_deg * (deg + 3 * avg_deg + (deg + avg_deg) * avg_deg / n) + 100;
#endif
}
#else
static inline double predict_workload(Graph *g, node_t root) {
    double deg = g->row_ptr[root + 1] - g->row_ptr[root];
    edge_ptr l = g->row_ptr[root], r = g->row_ptr[root + 1];
    while (l < r) {
        edge_ptr mid = (l + r) >> 1;
        if (g->col_idx[mid] < root) {
            l = mid + 1;
        }
        else {
            r = mid;
        }
    }
    double eff_deg = eff_deg = l - g->row_ptr[root];
    if (deg > MRAM_BUF_SIZE) {
        printf(ANSI_COLOR_RED "Error: deg too large\n" ANSI_COLOR_RESET);
        exit(1);
    }
    if (eff_deg > BITMAP_SIZE * 32) {
        printf(ANSI_COLOR_RED "Error: eff_deg too large\n" ANSI_COLOR_RESET);
        exit(1);
    }
    double avg_deg = (double)global_g->m / global_g->n;
    double n = global_g->n;
    (void)deg;
    (void)eff_deg;
    (void)avg_deg;
    (void)n;
#if defined(CLIQUE2)
    return eff_deg;
#elif defined(CLIQUE3)
    return eff_deg * eff_deg + 100;
#elif defined(CLIQUE4)
    return eff_deg * eff_deg * eff_deg + 100;
#elif defined(CLIQUE5)
    return eff_deg * eff_deg * eff_deg * eff_deg + 100;
#elif defined(CYCLE4)
    return eff_deg * eff_deg + 100;
#elif defined(HOUSE5)
    return eff_deg * deg * (2 + deg / n + avg_deg / n) + 100;
#elif defined(TRI_TRI6)
    return eff_deg * eff_deg * (deg + 3 * avg_deg + (deg + avg_deg) * avg_deg / n) + 100;
#endif
}
#endif

static void read_input() {
    FILE *fin = fopen(DATA_PATH, "rb");
    node_t n;
    edge_ptr m;
    edge_ptr *row_ptr = global_g->row_ptr;
    node_t *col_idx = global_g->col_idx;
    fread(&n, sizeof(node_t), 1, fin);
    fread(&m, sizeof(edge_ptr), 1, fin);
    fread(row_ptr, sizeof(edge_ptr), n, fin);
    fread(col_idx, sizeof(node_t), m, fin);
    row_ptr[n] = m;
    fclose(fin);
    global_g->n = n;
    global_g->m = m;
}

static void data_renumber() {
    static node_t rank[N];
    static node_t renumbered[N];
    Graph *tmp_g = malloc(sizeof(Graph));
    memcpy(tmp_g, global_g, sizeof(Graph));
    for (node_t i = 0; i < global_g->n; i++) {
        rank[i] = i;
    }
    qsort(rank, global_g->n, sizeof(node_t), deg_cmp);
    for (node_t i = 0; i < global_g->n; i++) {
        renumbered[rank[i]] = i;
    }
    edge_ptr cur = 0;
    for (node_t i = 0; i < global_g->n; i++) { 
        global_g->row_ptr[i] = cur;
        node_t node = rank[i];
        for (edge_ptr j = tmp_g->row_ptr[node]; j < tmp_g->row_ptr[node + 1]; j++) {
            global_g->col_idx[cur++] = renumbered[tmp_g->col_idx[j]];
        }
        qsort(global_g->col_idx + global_g->row_ptr[i], cur - global_g->row_ptr[i], sizeof(node_t), node_t_cmp);
    }
    free(tmp_g);
}

static inline bool check_in_bitmap(node_t n, uint32_t bitmap[N >> 5]) {
    return bitmap[n >> 5] & (1 << (n & 31));
}


static bool update_alloc_info(uint32_t dpu_id, node_t n, edge_ptr *m_count, bitmap_t bitmap) {
    // check condition
    if (global_g->root_num[dpu_id] == DPU_ROOT_NUM) {
        return false;
    }
    edge_ptr dpu_m_count = m_count[dpu_id];
    if (!check_in_bitmap(n, bitmap[dpu_id])) {
        dpu_m_count += global_g->row_ptr[n + 1] - global_g->row_ptr[n];
    }
    for (edge_ptr i = global_g->row_ptr[n]; i < global_g->row_ptr[n + 1]; i++) {
        node_t neighbor = global_g->col_idx[i];
        if (!check_in_bitmap(neighbor, bitmap[dpu_id])) {
            dpu_m_count += global_g->row_ptr[neighbor + 1] - global_g->row_ptr[neighbor];
        }
    }
    if (dpu_m_count > DPU_M) {
        return false;
    }

    // allocate
    global_g->roots[dpu_id][global_g->root_num[dpu_id]++] = n;
    m_count[dpu_id] = dpu_m_count;
    bitmap[dpu_id][n >> 5] |= (1 << (n & 31));
    for (edge_ptr i = global_g->row_ptr[n]; i < global_g->row_ptr[n + 1]; i++) {
        node_t neighbor = global_g->col_idx[i];
        bitmap[dpu_id][neighbor >> 5] |= (1 << (neighbor & 31));
    }
    return true;
}

void queue_init();
void push_to_queue(uint32_t dpu_id, double work_load);
uint32_t pop_from_queue();

static void data_allocate(bitmap_t bitmap) {
    static edge_ptr m_count[EF_NR_DPUS];   // edges put in dpu
    static node_t allocate_rank[N];
    static double dpu_workload[EF_NR_DPUS];

    memset(bitmap, 0, (size_t)(N >> 3) * EF_NR_DPUS);
    for (uint32_t i = 0; i < EF_NR_DPUS; i++) {
        global_g->root_num[i] = 0;
        global_g->roots[i] = malloc(DPU_ROOT_NUM * sizeof(node_t));
    }

    for (node_t i = 0; i < global_g->n; i++) {
        allocate_rank[i] = i;
    }
    for (node_t i = 0; i < global_g->n; i++) {
        workload[i] = predict_workload(global_g, i);
    }
    qsort(allocate_rank, global_g->n, sizeof(node_t), workload_cmp);

    queue_init();
    uint32_t full_dpu_ct = 0;
    for (node_t i = 0; i < global_g->n; i++) {
        while (full_dpu_ct != EF_NR_DPUS) {
            uint32_t cur_dpu = pop_from_queue();
            if (update_alloc_info(cur_dpu, allocate_rank[i], m_count, bitmap)) {
                dpu_workload[cur_dpu] += workload[allocate_rank[i]];
                push_to_queue(cur_dpu, dpu_workload[cur_dpu]);
                break;
            }
            else {
                full_dpu_ct++;
            }
        }
        if (full_dpu_ct == EF_NR_DPUS) {
            printf(ANSI_COLOR_RED "Error: not enough DPUs\n" ANSI_COLOR_RESET);
            exit(1);
        }
    }
}


//rCSR format
void data_compact(struct dpu_set_t set, bitmap_t bitmap,int base) {
    edge_ptr(*dpu_row_ptr)[DPU_N * 2];
    dpu_row_ptr = malloc(NR_DPUS * DPU_N * 2 * sizeof(edge_ptr));
    node_t(*dpu_col_idx)[DPU_M * 2];
    dpu_col_idx = malloc(NR_DPUS * DPU_M * 2 * sizeof(node_t));
    node_t(*dpu_roots)[DPU_ROOT_NUM];
    dpu_roots = malloc(NR_DPUS * DPU_ROOT_NUM * sizeof(node_t));
    static uint64_t processed_row_size[NR_DPUS];
    memset(processed_row_size, 0, NR_DPUS * sizeof(uint64_t));
    static uint64_t processed_col_size[NR_DPUS];
    memset(processed_col_size, 0, NR_DPUS * sizeof(uint64_t));
    static uint64_t tmp_row_size[NR_DPUS];
    static uint64_t tmp_col_size[NR_DPUS];
    
    struct dpu_set_t dpu;
    uint32_t each_dpu;
#ifndef LOAD_FROM_CACHE
    DPU_ASSERT(dpu_load(set, DPU_ALLOC_BINARY, NULL));

    uint64_t mode = 0;
    DPU_ASSERT(dpu_broadcast_to(set, "mode", 0, &mode, sizeof(uint64_t), DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, bitmap[each_dpu+base]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "bitmap", 0, (N >> 5) * sizeof(uint32_t), DPU_XFER_DEFAULT));

    static uint32_t zero[N >> 5];
    memset(zero, 0, sizeof(zero));

    DPU_ASSERT(dpu_broadcast_to(set, "involve_bitmap", 0, zero, sizeof(zero), DPU_XFER_DEFAULT));
    uint64_t start = 0;
    while (start < global_g->n) {
            //HERE_OKF("initialize ok");
        uint64_t size = 0;
        while (start + size < global_g->n && global_g->row_ptr[start + size + 1] - global_g->row_ptr[start] < PARTITION_M) {
            size++;
        }
        DPU_ASSERT(dpu_broadcast_to(set, "start", 0, &start, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "size", 0, &size, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "row_ptr", 0, &global_g->row_ptr[start], ALIGN8((size + 1) * sizeof(edge_ptr)), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "col_idx", 0, &global_g->col_idx[global_g->row_ptr[start]], ALIGN8((global_g->row_ptr[start + size] - global_g->row_ptr[start]) * sizeof(node_t)), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
        start += size;
    }



    mode = 1;
    DPU_ASSERT(dpu_broadcast_to(set, "mode", 0, &mode, sizeof(uint64_t), DPU_XFER_DEFAULT));
    uint64_t n_size = global_g->n;
    DPU_ASSERT(dpu_broadcast_to(set, "size", 0, &n_size, sizeof(uint64_t), DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, each_dpu) {
        uint64_t root_num = global_g->root_num[each_dpu+base];
        DPU_ASSERT(dpu_copy_to(dpu, "root_size", 0, &root_num, sizeof(uint64_t)));
    }
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, global_g->roots[each_dpu+base]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "roots", 0, DPU_ROOT_NUM * sizeof(node_t), DPU_XFER_DEFAULT));
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_roots[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "roots", 0, DPU_ROOT_NUM * sizeof(node_t), DPU_XFER_DEFAULT));

    mode = 2;
    DPU_ASSERT(dpu_broadcast_to(set, "mode", 0, &mode, sizeof(uint64_t), DPU_XFER_DEFAULT));

    start = 0;
    while (start < global_g->n) {
        uint64_t size = 0;
        while (start + size < global_g->n && global_g->row_ptr[start + size + 1] - global_g->row_ptr[start] < PARTITION_M) {
            size++;
        }
        DPU_ASSERT(dpu_broadcast_to(set, "start", 0, &start, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "size", 0, &size, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &processed_col_size[each_dpu]));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "processed_offset", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "row_ptr", 0, &global_g->row_ptr[start], ALIGN8((size + 1) * sizeof(edge_ptr)), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_broadcast_to(set, "col_idx", 0, &global_g->col_idx[global_g->row_ptr[start]], ALIGN8((global_g->row_ptr[start + size] - global_g->row_ptr[start]) * sizeof(node_t)), DPU_XFER_DEFAULT));
        DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &tmp_row_size[each_dpu]));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "processed_row_size", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &tmp_col_size[each_dpu]));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "processed_col_size", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));
        uint64_t max_row_size = 0;
        uint64_t max_col_size = 0;
        DPU_FOREACH(set, dpu, each_dpu) {
            if (tmp_row_size[each_dpu] > max_row_size) {
                max_row_size = tmp_row_size[each_dpu];
            }
            if (tmp_col_size[each_dpu] > max_col_size) {
                max_col_size = tmp_col_size[each_dpu];
            }
        }
        if (max_row_size != 0) {
            DPU_FOREACH(set, dpu, each_dpu) {
                DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_row_ptr[each_dpu][processed_row_size[each_dpu]]));
            }
            DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "processed_row_ptr", 0, ALIGN8(max_row_size * sizeof(edge_ptr)), DPU_XFER_DEFAULT));
        }
        if (max_col_size != 0) {
            DPU_FOREACH(set, dpu, each_dpu) {
                DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_col_idx[each_dpu][processed_col_size[each_dpu]]));
            }
            DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "processed_col_idx", 0, ALIGN8(max_col_size * sizeof(node_t)), DPU_XFER_DEFAULT));
        }
        DPU_FOREACH(set, dpu, each_dpu) {
            processed_row_size[each_dpu] += tmp_row_size[each_dpu];
            processed_col_size[each_dpu] += tmp_col_size[each_dpu];
        }
        start += size;
    }
    for (uint32_t i = 0; i < NR_DPUS; i++) {
        dpu_row_ptr[i][processed_row_size[i]] = processed_col_size[i];
    }
#endif

#ifdef SAVE_TO_CACHE
// 分别将每个 DPU 应该具有的数据保存到一个文件（一个 DPU 对应一个文件，命名规范：<DATA_NAME>_cache_<dpu_Id+base>.bin）
// 需要包含的数据：root_num, roots, row_ptr, col_idx
// 注意要方便读取，所以要把 roots, row_ptr, col_idx 都分别padding 到 DPU_ROOT_NUM, DPU_N, DPU_M 的大小
    char filename[100];
    for (uint32_t i = 0; i < NR_DPUS; i++) { // NR_DPUS 是当前 dpu_set 中的 DPU 数量
        sprintf(filename, "%s/%s_cache_%d.bin", CACHE_DIR, DATA_NAME, i + base);
        FILE *fout = fopen(filename, "wb");
        if (!fout) {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            exit(EXIT_FAILURE);
        }

        // 1. 写入 root_num (实际的根节点数量)
        uint64_t actual_root_count_for_dpu = global_g->root_num[i + base];
        fwrite(&actual_root_count_for_dpu, sizeof(uint64_t), 1, fout);

        // 2. 写入 roots 数组 (填充到 DPU_ROOT_NUM)；实际的根节点数据在 dpu_roots[i] 中，数量为 actual_root_count_for_dpu
        node_t *padded_roots = (node_t *)calloc(DPU_ROOT_NUM, sizeof(node_t)); // 注意：这里使用 calloc 来确保分配的内存被初始化为 0
        if (!padded_roots && DPU_ROOT_NUM > 0) {
            fprintf(stderr, "Error: calloc failed for padded_roots for DPU %d\n", i + base);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        if (actual_root_count_for_dpu > DPU_ROOT_NUM) {
            fprintf(stderr, ANSI_COLOR_RED "Error: actual_root_count (%llu) > DPU_ROOT_NUM (%u) for DPU %d. Cache saving failed.\n" ANSI_COLOR_RESET, (unsigned long long)actual_root_count_for_dpu, DPU_ROOT_NUM, i + base);
            free(padded_roots);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        memcpy(padded_roots, dpu_roots[i], actual_root_count_for_dpu * sizeof(node_t));
        fwrite(padded_roots, sizeof(node_t), DPU_ROOT_NUM, fout);
        free(padded_roots);

        // 3. 写入 row_ptr 数组 (填充到 DPU_N)；实际的 row_ptr 数据在 dpu_row_ptr[i] 中；processed_row_size[i] 是这个 DPU 分配到的节点数（行数）；CSR 格式中 row_ptr 的长度是节点数 + 1
        uint64_t actual_row_ptr_elements = processed_row_size[i] + 1;
        edge_ptr *padded_row_ptr = (edge_ptr *)calloc(DPU_N, sizeof(edge_ptr));
        if (!padded_row_ptr && DPU_N > 0) {
            fprintf(stderr, "Error: calloc failed for padded_row_ptr for DPU %d\n", i + base);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        if (actual_row_ptr_elements > DPU_N) {
            fprintf(stderr, ANSI_COLOR_RED "Error: actual_row_ptr_elements (%llu) > DPU_N (%u) for DPU %d. Cache saving failed.\n" ANSI_COLOR_RESET, (unsigned long long)actual_row_ptr_elements, DPU_N, i + base);
            free(padded_row_ptr);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        memcpy(padded_row_ptr, dpu_row_ptr[i], actual_row_ptr_elements * sizeof(edge_ptr));
        fwrite(padded_row_ptr, sizeof(edge_ptr), DPU_N, fout);
        free(padded_row_ptr);

        // 4. 写入 col_idx 数组 (填充到 DPU_M)；实际的 col_idx 数据在 dpu_col_idx[i] 中；processed_col_size[i] 是这个 DPU 分配到的边数（列数）
        uint64_t actual_col_idx_elements = processed_col_size[i];
        node_t *padded_col_idx = (node_t *)calloc(DPU_M, sizeof(node_t));
        if (!padded_col_idx && DPU_M > 0) {
            fprintf(stderr, "Error: calloc failed for padded_col_idx for DPU %d\n", i + base);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        if (actual_col_idx_elements > DPU_M) {
            fprintf(stderr, ANSI_COLOR_RED "Error: actual_col_idx_elements (%llu) > DPU_M (%u) for DPU %d. Cache saving failed.\n" ANSI_COLOR_RESET, (unsigned long long)actual_col_idx_elements, DPU_M, i + base);
            free(padded_col_idx);
            fclose(fout);
            exit(EXIT_FAILURE);
        }
        memcpy(padded_col_idx, dpu_col_idx[i], actual_col_idx_elements * sizeof(node_t));
        fwrite(padded_col_idx, sizeof(node_t), DPU_M, fout);
        free(padded_col_idx);

        fclose(fout);
        // HERE_OKF 是一个宏，用于打印调试信息，这里保留其原始意图
        // 假设 HERE_OKF 类似于 printf
    }
    (void)printf(ANSI_COLOR_GREEN "INFO:" ANSI_COLOR_RESET " Cache Data for DPU batch %d saved\n", base/NR_DPUS);
#endif

#ifdef LOAD_FROM_CACHE
    char filename[128];
    struct dpu_set_t dpu_load_iter; 
    uint32_t dpu_rank_load;   

    DPU_FOREACH(set, dpu_load_iter, dpu_rank_load) {
        uint32_t id_for_file = dpu_rank_load + base; // ID 用于文件名和 global_g->root_num

        sprintf(filename, "%s/%s_cache_%u.bin", CACHE_DIR,DATA_NAME, id_for_file);
        FILE *fin = fopen(filename, "rb");
        if (!fin) {
            fprintf(stderr, ANSI_COLOR_RED "FATAL: Cache file %s not found. LOAD_FROM_CACHE requires all cache files to exist.\n" ANSI_COLOR_RESET, filename);
            perror("fopen cache for read failed");
            // 清理已分配的内存是一个好习惯，但这里直接退出
            free(dpu_row_ptr);
            free(dpu_col_idx);
            free(dpu_roots);
            exit(EXIT_FAILURE);
        }

        // 1. 读取并恢复 global_g->root_num
        if (fread(&global_g->root_num[id_for_file], sizeof(uint64_t), 1, fin) != 1) {
            fprintf(stderr, ANSI_COLOR_RED "Error reading root_num from cache file %s.\n" ANSI_COLOR_RESET, filename);
            fclose(fin); free(dpu_row_ptr); free(dpu_col_idx); free(dpu_roots); exit(EXIT_FAILURE);
        }

        // 4. 读取填充后的 roots 数据到 dpu_roots[dpu_rank_load]
        if (DPU_ROOT_NUM > 0 && fread(dpu_roots[dpu_rank_load], sizeof(node_t), DPU_ROOT_NUM, fin) != DPU_ROOT_NUM) {
            fprintf(stderr, ANSI_COLOR_RED "Error reading roots from cache file %s. Expected %u elements.\n" ANSI_COLOR_RESET, filename, DPU_ROOT_NUM);
            fclose(fin); free(dpu_row_ptr); free(dpu_col_idx); free(dpu_roots); exit(EXIT_FAILURE);
        }

        // 5. 读取填充后的 row_ptr 数据到 dpu_row_ptr[dpu_rank_load]
        // dpu_row_ptr[dpu_rank_load] 指向一个 DPU_N * 2 * sizeof(edge_ptr) 的缓冲区，我们只读取 DPU_N 个元素
        if (DPU_N > 0 && fread(dpu_row_ptr[dpu_rank_load], sizeof(edge_ptr), DPU_N, fin) != DPU_N) {
            fprintf(stderr, ANSI_COLOR_RED "Error reading row_ptr from cache file %s. Expected %u elements.\n" ANSI_COLOR_RESET, filename, DPU_N);
            fclose(fin); free(dpu_row_ptr); free(dpu_col_idx); free(dpu_roots); exit(EXIT_FAILURE);
        }

        // 6. 读取填充后的 col_idx 数据到 dpu_col_idx[dpu_rank_load]
        int actual_r = 0;
        if (DPU_M > 0 && (actual_r = fread(dpu_col_idx[dpu_rank_load], sizeof(node_t), DPU_M, fin)) != DPU_M) {
            fprintf(stderr, ANSI_COLOR_RED "Error reading col_idx from cache file %s. Expected %u elements, got %d.\n" ANSI_COLOR_RESET, filename, DPU_M, actual_r);
            fclose(fin); free(dpu_row_ptr); free(dpu_col_idx); free(dpu_roots); exit(EXIT_FAILURE);
        }
        
        fclose(fin);
    }

    // 恢复CSR row_ptr的最后一个关键元素
    DPU_FOREACH(set, dpu_load_iter, dpu_rank_load) {
        // 验证索引的有效性：
        if (processed_row_size[dpu_rank_load] < DPU_N || (DPU_N == 0 && processed_row_size[dpu_rank_load] == 0) ) {
             dpu_row_ptr[dpu_rank_load][processed_row_size[dpu_rank_load]] = processed_col_size[dpu_rank_load];
        } else if (DPU_N > 0 && processed_row_size[dpu_rank_load] == DPU_N) { // processed_row_size can be DPU_N if it means DPU_N nodes, index is DPU_N, but array has DPU_N elements [0...DPU_N-1]
             dpu_row_ptr[dpu_rank_load][processed_row_size[dpu_rank_load]] = processed_col_size[dpu_rank_load];
        } else {
             fprintf(stderr, ANSI_COLOR_RED "Warning: processed_row_size[%u] = %lu relative to DPU_N = %u may require careful review for index %s:%d.\n" ANSI_COLOR_RESET,
                dpu_rank_load, processed_row_size[dpu_rank_load], DPU_N, __FILE__, __LINE__);
            // Fallback to direct assignment if index is within allocated bounds (DPU_N * 2)
            if (processed_row_size[dpu_rank_load] < DPU_N * 2) {
                 dpu_row_ptr[dpu_rank_load][processed_row_size[dpu_rank_load]] = processed_col_size[dpu_rank_load];
            } else {
                 fprintf(stderr, ANSI_COLOR_RED "Error: processed_row_size[%u]=%lu is out of bounds for dpu_row_ptr.\n" ANSI_COLOR_RESET,
                    dpu_rank_load, processed_row_size[dpu_rank_load]);
            }
        }
    }
    (void)printf(ANSI_COLOR_GREEN "INFO: Data for DPU set (base %d) successfully loaded from cache.\n" ANSI_COLOR_RESET, base);
#endif // End of #ifdef LOAD_FROM_CACHE

    DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
    DPU_FOREACH(set, dpu, each_dpu) {
        uint64_t root_num = global_g->root_num[each_dpu+base];
        DPU_ASSERT(dpu_copy_to(dpu, "root_num", 0, &root_num, sizeof(uint64_t)));
    }
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_roots[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "roots", 0, DPU_ROOT_NUM * sizeof(node_t), DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_row_ptr[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "row_ptr", 0, DPU_N * sizeof(edge_ptr), DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, each_dpu) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_col_idx[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "col_idx", 0, DPU_M * sizeof(node_t), DPU_XFER_DEFAULT));

    free(dpu_row_ptr);
    free(dpu_col_idx);
    free(dpu_roots);

}

void data_xfer(struct dpu_set_t set,int base) {
       struct dpu_set_t dpu;
        uint32_t each_dpu;

        DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
        node_t max_root_num = 0;
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &global_g->root_num[each_dpu+base]));
            if(global_g->root_num[each_dpu+base] > max_root_num) {
                max_root_num = global_g->root_num[each_dpu+base];
            }
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "root_num", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, global_g->roots[each_dpu+base]));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "roots", 0, ALIGN8(max_root_num * sizeof(node_t)), DPU_XFER_DEFAULT));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, global_g->row_ptr));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "row_ptr", 0, ALIGN8((global_g->n + 1) * sizeof(edge_ptr)), DPU_XFER_DEFAULT));
        DPU_FOREACH(set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, global_g->col_idx));
        }
        DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "col_idx", 0, ALIGN8(global_g->m * sizeof(node_t)), DPU_XFER_DEFAULT));
}

bitmap_t prepare_graph() {
    read_input();  
    data_renumber();     
    bitmap = malloc(sizeof(uint32_t) * (N >> 5) * EF_NR_DPUS);
    data_allocate(bitmap);  
    return bitmap;
}


void data_transfer(struct dpu_set_t set, Graph *g ,bitmap_t bitmap ,int base) {
    
#ifdef NO_PARTITION_AS_POSSIBLE
    if (global_g->n > DPU_N - 1 || global_g->m > DPU_M) {
#endif
        data_compact(set, bitmap,base);
        HERE_OKF("data_compact ok");
#ifdef NO_PARTITION_AS_POSSIBLE
    }
    else {

        data_xfer(set,base);
    }
#endif

}