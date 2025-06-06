#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stdio.h> // for debug

// LOAD 和 SAVE 不能同时开启。控制 DPU 跳过 mode0-2，直接加载处理之后的数据。
#define CACHE_DIR "./cached_data"
#define LOAD_FROM_CACHE
// #define SAVE_TO_CACHE
#if defined(LOAD_FROM_CACHE) && defined(SAVE_TO_CACHE)
#error "Error: LOAD_FROM_CACHE and SAVE_TO_CACHE cannot be enabled simultaneously."
#endif

#define PERF
// #define DPU_LOG
// #define CPU_RUN
// #define NO_PARTITION_AS_POSSIBLE
// #define MORE_ACCURATE_MODEL
#if defined(CLIQUE4) || defined(CLIQUE5)
#define BITMAP
#endif

#define DATA_DIR "./data/"
#if defined(SELF)
#define DATA_NAME "self-defined"
#define N (1 << 5)
#define M (1 << 5)
#elif defined(WV)
#define DATA_NAME "Wiki-Vote"
#define N (1 << 13)
#define M (1 << 18)
#elif defined(PP)
#define DATA_NAME "p2p-Gnutella04"
#define N (1 << 14)
#define M (1 << 17)
#elif defined(CA)
#define DATA_NAME "CA-AstroPh"
#define N (1 << 15)
#define M (1 << 19)
#elif defined(YT)
#define DATA_NAME "com-youtube"
#define N (1 << 21)
#define M (1 << 23)
#elif defined(PT)
#define DATA_NAME "cit-Patents"
#define N (1 << 22)
#define M (1 << 26)
#elif defined(LJ)
#define DATA_NAME "soc-LiveJournal1"
#define N (1 << 23)
#define M (1 << 27)
#elif defined(PA)
#define DATA_NAME "roadNet-PA_adj"
#define N (1 << 23)
#define M (1 << 25)
#elif defined(SP)
#define DATA_NAME "simple_csr_graph"
#define N (1 << 6)
#define M (1 << 8)
#elif defined(TH)
#define DATA_NAME "Theory-3-4-5-9-B1k"
#define N (1 << 12)
#define M (1 << 14)
#elif defined(FE)
#define DATA_NAME "flickrEdges_adj"
#define N (1 << 17)
#define M (1 << 23)
#elif defined(CH)
#define DATA_NAME "cit-HepPh_adj"
#define N (1 << 18)
#define M (1 << 22)
#elif defined(SS)
#define DATA_NAME "soc-Slashdot0811_adj"
#define N (1 << 18)
#define M (1 << 22)
#elif defined(NCA)
#define DATA_NAME "roadNet-CA_adj"
#define N (1 << 23)
#define M (1 << 25)
#elif defined(AM0312)
#define DATA_NAME "amazon0312_adj"
#define N (1 << 20)
#define M (1 << 23)
#elif defined(AM0505)
#define DATA_NAME "amazon0505_adj"
#define N (1 << 20)
#define M (1 << 23)
#elif defined(SC18_16)
#define DATA_NAME "graph500-scale18-ef16_adj"
#define N (1 << 20)
#define M (1 << 24)
#else
#warning "No graph selected, fall back to PP."
#define DATA_NAME "p2p-Gnutella04"
#define N (1 << 14)
#define M (1 << 17)
#endif
#define DATA_PATH DATA_DIR DATA_NAME ".bin"

#ifndef NR_DPUS
#warning "No NR_DPUS defined, fall back to 1."
#define NR_DPUS 1
#endif
#ifndef NR_TASKLETS
#warning "No NR_TASKLETS defined, fall back to 1."
#define NR_TASKLETS 1
#endif
#ifndef DPU_BINARY
#warning "No DPU_BINARY defined, fall back to bin/dpu."
#define DPU_BINARY "bin/dpu"
#endif
#ifndef DPU_ALLOC_BINARY
#warning "No DPU_ALLOC_BINARY defined, fall back to bin/dpu_alloc."
#define DPU_ALLOC_BINARY "bin/dpu_alloc"
#endif

#if defined(CLIQUE2)
#define KERNEL_FUNC clique2
#define PATTERN_NAME "clique2"
#elif defined(CLIQUE3)
#define KERNEL_FUNC clique3
#define PATTERN_NAME "clique3"
#elif defined(CLIQUE4)
#define KERNEL_FUNC clique4
#define PATTERN_NAME "clique4"
#elif defined(CLIQUE5)
#define KERNEL_FUNC clique5
#define PATTERN_NAME "clique5"
#elif defined(CYCLE4)
#define KERNEL_FUNC cycle4
#define PATTERN_NAME "cycle4"
#elif defined(HOUSE5)
#define KERNEL_FUNC house5
#define PATTERN_NAME "house5"
#elif defined(TRI_TRI6)
#define KERNEL_FUNC tri_tri6
#define PATTERN_NAME "tri_tri6"
#else
#warning "No kernel function selected, fall back to clique2."
#define KERNEL_FUNC clique2
#define PATTERN_NAME "clique2"
#endif

#define node_t uint32_t
#define edge_ptr uint32_t
#define ans_t uint64_t
#define SIZE_NODE_T_LOG 2
#define SIZE_EDGE_PTR_LOG 2
#define INVALID_NODE ((node_t)(-1))
#define DPU_N ((1 << 24) / sizeof(edge_ptr))
#define DPU_M ((1 << 25) / sizeof(node_t))
#define DPU_ROOT_NUM ((1 << 20) / sizeof(node_t))
#define BITMAP_SIZE 32 // 1024 bits
#define BUF_SIZE 32
#define MRAM_BUF_SIZE 32768
#define BRANCH_LEVEL_THRESHOLD 16
#define PARTITION_M ((1 << 22) / sizeof(node_t))

// #define DC //detailed cycle count
// #define DC_NUM 3 //number of cycle count list

#define V_NR_DPUS 2560
#define BATCH_SIZE (V_NR_DPUS / NR_DPUS)
#ifdef V_NR_DPUS
#define EF_NR_DPUS V_NR_DPUS
#else
#define EF_NR_DPUS NR_DPUS
#endif

typedef struct Graph
{
    node_t n;   // number of vertices
    edge_ptr m; // number of edges
    edge_ptr row_ptr[N];
    node_t col_idx[M];
    uint64_t root_num[EF_NR_DPUS]; // number of search roots allocated to dpu
    node_t *roots[EF_NR_DPUS];
} Graph;
typedef uint32_t (*bitmap_t)[N >> 5];

#define ALIGN(x, a) (((x) + (a) - 1) & ~((a) - 1))
#define ALIGN2(x) ALIGN(x, 2)
#define ALIGN4(x) ALIGN(x, 4)
#define ALIGN8(x) ALIGN(x, 8)
#define ALIGN16(x) ALIGN(x, 16)
#define ALIGN_LOWER(x, a) ((x) & ~((a) - 1))
#define ALIGN2_LOWER(x) ALIGN_LOWER(x, 2)
#define ALIGN4_LOWER(x) ALIGN_LOWER(x, 4)
#define ALIGN8_LOWER(x) ALIGN_LOWER(x, 8)
#define ALIGN16_LOWER(x) ALIGN_LOWER(x, 16)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define HERE_OK() printf(ANSI_COLOR_GREEN "[OK] %s:%d\n" ANSI_COLOR_RESET, __FILE__, __LINE__); // DEBUG
#define HERE_OKF(fmt, ...) \
    printf(ANSI_COLOR_GREEN "[OK] %s:%d - " fmt "\n" ANSI_COLOR_RESET, __FILE__, __LINE__, ##__VA_ARGS__);

#endif // COMMON_H
