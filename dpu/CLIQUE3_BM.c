#include <dpu_bitmine.h>

static ans_t partial_ans[NR_TASKLETS];
static uint64_t partial_cycle[NR_TASKLETS];
static perfcounter_cycles cycles[NR_TASKLETS];

static ans_t __imp_clique3_2(sysname_t tasklet_id, node_t __mram_ptr * root_col, node_t root_size, node_t __mram_ptr * second_col, node_t second_size) {

    if(!second_size)return 0;
    
    node_t(*tasklet_buf)[BUF_SIZE] = buf[tasklet_id];
    node_t ans =  intersect_seq_buf_thresh(tasklet_buf, root_col, root_size, second_col, second_size);

    return ans;
}

static ans_t __imp_clique3(sysname_t tasklet_id, node_t root) {
    edge_ptr root_begin = row_ptr[root];  // intended DMA
    edge_ptr root_end = row_ptr[root + 1];  // intended DMA
    node_t root_size = root_end - root_begin;
    if(!root_size)return 0;
    ans_t ans = 0;

    mram_read(&col_idx[edge_offset+2*root_begin],col_buf[tasklet_id],MIN(16,root_end-root_begin)<<(SIZE_EDGE_PTR_LOG+1));

    for (edge_ptr i = 1; i<root_size; i++) {
        ans += __imp_clique3_2(tasklet_id,&col_idx[root_begin],i,&col_idx[col_buf[tasklet_id][2*i]],col_buf[tasklet_id][2*i+1]-col_buf[tasklet_id][2*i]);
    }

    return ans;
}


//func begin
extern void clique3( sysname_t tasklet_id )
{
	node_t i = 0;                          /* if all node is large_degree */
	while ( i < root_num )
	{
		node_t	root		= roots[i];             /* intended DMA */
		node_t	root_begin	= row_ptr[root];        /* intended DMA */
		node_t	root_end	= row_ptr[root + 1];    /* intended DMA */
		node_t	root_size	= root_end - root_begin;
		if ( root_size < BRANCH_LEVEL_THRESHOLD )break;


#ifdef PERF
		timer_start( &cycles[tasklet_id] );
#endif
		partial_ans[tasklet_id] = 0;

			 for ( edge_ptr j = root_begin + tasklet_id + 1; j < root_end; j += NR_TASKLETS )
			 {
			     partial_ans[tasklet_id] += __imp_clique3_2( tasklet_id, &col_idx[root_begin], j - root_begin, &col_idx[col_idx[edge_offset + 2 * j]], col_idx[edge_offset + 2 * j + 1] - col_idx[edge_offset + 2 * j] );
			 }
			 


#ifdef PERF
		partial_cycle[tasklet_id] = timer_stop( &cycles[tasklet_id] );
#endif
		barrier_wait( &co_barrier );
		if ( tasklet_id == 0 )
		{
			ans_t total_ans = 0;
#ifdef PERF
			uint64_t total_cycle = 0;
#endif
			for ( uint32_t j = 0; j < NR_TASKLETS; j++ )
			{
				total_ans += partial_ans[j];
#ifdef PERF
				total_cycle += partial_cycle[j];
#endif
			}
			ans[i] = total_ans;             /* intended DMA */
#ifdef PERF
			cycle_ct[i] = total_cycle;      /* intended DMA */
#endif
		}
		i++;

		barrier_wait( &co_barrier );
	}

	for ( i += tasklet_id; i < root_num; i += NR_TASKLETS )
	{
		node_t root = roots[i];                                         /* intended DMA */

#ifdef PERF
		timer_start( &cycles[tasklet_id] );
#endif

			ans[i] = __imp_clique3( tasklet_id, root );             /* intended DMA */

#ifdef PERF
		cycle_ct[i] = timer_stop( &cycles[tasklet_id] );                /* intended DMA */
#endif
	}
}