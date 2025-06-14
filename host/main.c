#include <common.h>
#include <timer.h>
#include <assert.h>
#include <stdio.h>
#include <dpu.h>

extern void data_transfer(struct dpu_set_t set, Graph *g ,bitmap_t bitmap,int base);
extern ans_t clique2(Graph *g, node_t root);
extern ans_t KERNEL_FUNC(Graph *g, node_t root);
extern Graph *global_g;
extern bitmap_t bitmap;
Graph *g;
ans_t ans[N];
ans_t result[N];
Timer timer;
uint64_t cycle_ct[N];
uint64_t cycle_ct_dpu[EF_NR_DPUS][NR_TASKLETS];
node_t large_degree_num[EF_NR_DPUS];


int main() {
    printf("NR_DPUS: %u, NR_TASKLETS: %u, DPU_BINARY: %s, PATTERN: %s\n", NR_DPUS, NR_TASKLETS, DPU_BINARY, PATTERN_NAME);

    struct dpu_set_t set, dpu;
    DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, &set));

    // task allocation and data partition
    printf("Selecting graph: %s\n", DATA_PATH);
    start(&timer, 0, 0);
    g = malloc(sizeof(Graph));
    global_g = g;
    //bitmap=prepare_graph(); //bug not reslove
    prepare_graph();
 
    stop(&timer, 0);
    printf("Data prepare ");
    print(&timer, 0, 1);

ans_t total_ans = 0;
#ifdef PERF
    uint64_t total_cycle_ct = 0;
#endif

int batch_count = 1;
int base = 0;
#ifdef V_NR_DPUS
batch_count = BATCH_SIZE;
#endif

for(int index=0;index<batch_count;index++){
    HERE_OKF(" batch index %d begin...", index); 
#ifdef V_NR_DPUS
    base = index * NR_DPUS;
#endif
    data_transfer(set,g,bitmap,base);
    // run it on DPU
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // collect answer and cycle count
    bool fine = true;
    bool finished, failed;
    uint32_t each_dpu;


    DPU_FOREACH(set, dpu, each_dpu) {
        // check status
        DPU_ASSERT(dpu_status(dpu, &finished, &failed));
        if (failed) {
            printf("DPU: %u failed\n", each_dpu);
            fine = false;
            break;
        }
        // collect answer
        uint64_t *dpu_ans = (uint64_t *)malloc(ALIGN8(g->root_num[each_dpu+base] * sizeof(uint64_t)));
        DPU_ASSERT(dpu_copy_from(dpu, "ans", 0, dpu_ans, ALIGN8(g->root_num[each_dpu+base] * sizeof(uint64_t))));
        for (node_t k = 0; k < g->root_num[each_dpu+base]; k++) {
            node_t cur_root = g->roots[each_dpu+base][k];
            result[cur_root] = dpu_ans[k];
            total_ans += dpu_ans[k];
        }
        free(dpu_ans);

        // collect cycle count
#ifdef PERF
        uint64_t *dpu_cycle_ct = (uint64_t *)malloc(ALIGN8(g->root_num[each_dpu+base] * sizeof(uint64_t)));
        DPU_ASSERT(dpu_copy_from(dpu, "cycle_ct", 0, dpu_cycle_ct, ALIGN8(g->root_num[each_dpu+base] * sizeof(uint64_t))));
        DPU_ASSERT(dpu_copy_from(dpu, "large_degree_num", 0, large_degree_num[each_dpu+base], sizeof(node_t)));
        for (node_t k = 0, cur_thread = 0; k < g->root_num[each_dpu+base]; k++) {
            node_t cur_root = g->roots[each_dpu+base][k];
            cycle_ct[cur_root] = dpu_cycle_ct[k];
            if (g->row_ptr[cur_root + 1] - g->row_ptr[cur_root] >= BRANCH_LEVEL_THRESHOLD) {
                for (uint32_t i = 0; i < NR_TASKLETS; i++) {
                    cycle_ct_dpu[each_dpu+base][i] += dpu_cycle_ct[k] / NR_TASKLETS;
                }
                total_cycle_ct += dpu_cycle_ct[k];
            }
            else {
                cycle_ct_dpu[each_dpu+base][cur_thread] += dpu_cycle_ct[k];
                cur_thread = (cur_thread + 1) % NR_TASKLETS;
                total_cycle_ct += dpu_cycle_ct[k];
            }
        }
        free(dpu_cycle_ct);
#endif

    }
    if (!fine){printf(ANSI_COLOR_RED "Some failed\n" ANSI_COLOR_RESET);}
}


    printf("DPU ans: %lu\n", total_ans);
#ifdef PERF
    printf("Lower bound: %f\n", (double)total_cycle_ct / NR_DPUS / NR_TASKLETS / 350000);
#endif

#ifdef V_NR_DPUS
    printf(ANSI_COLOR_GREEN"[INFO] Finished in VIRTUAL DPU mode (V_NR_DPUS = %d)\n"ANSI_COLOR_RESET, V_NR_DPUS);
#else
    printf(ANSI_COLOR_GREEN"[INFO] Finished in PHYSICAL DPU mode (NR_DPUS = %d)\n"ANSI_COLOR_RESET, NR_DPUS);
#endif

    // output result to file
#ifdef PERF
#ifdef NO_RUN   //test cycle without Intersection operation 
    FILE *fp = fopen("./result/" PATTERN_NAME "_" DATA_NAME "_NO_RUN.txt", "w");
#else
    FILE *fp = fopen("./result/" PATTERN_NAME "_" DATA_NAME ".txt", "w");
#endif
    fprintf(fp, "NR_DPUS: %u, NR_TASKLETS: %u, DPU_BINARY: %s, PATTERN: %s\n", NR_DPUS, NR_TASKLETS, DPU_BINARY, PATTERN_NAME);
    fprintf(fp, "N: %u, M: %u, avg_deg: %f\n", g->n, g->m, (double)g->m / g->n);

    // for (node_t i = 0; i < g->n; i++) {
    //     fprintf(fp, "node: %u, deg: %u, o_deg: %lu, ans: %lu, cycle: %lu,\n", i, g->row_ptr[i + 1] - g->row_ptr[i], clique2(g, i), result[i], cycle_ct[i]);
    // }

    for (uint32_t i = 0; i < EF_NR_DPUS; i++) {
        for (uint32_t j = 0; j < NR_TASKLETS; j++) {
            fprintf(fp, "DPU: %u, tasklet: %u, cycle: %lu, root_num: %lu\n", i, j, cycle_ct_dpu[i][j], g->root_num[i]);
        }
    }

for (uint32_t i = 0; i < EF_NR_DPUS; i++) {
        float ratio = (float)large_degree_num[i] / g->root_num[i];
        fprintf(fp, "DPU: %u, large_degree_num: %u, root_num: %lu, ratio: %.2f\n",i, large_degree_num[i], g->root_num[i], ratio);
}

    fclose(fp);
#endif
    assert(bitmap != NULL);
    free(bitmap);
    free(g);
    DPU_ASSERT(dpu_free(set));
    return 0;
}