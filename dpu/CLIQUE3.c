#include <dpu_mine.h>

static ans_t partial_ans[NR_TASKLETS];
static uint64_t partial_cycle[NR_TASKLETS];
static perfcounter_cycles cycles[NR_TASKLETS];

#ifdef BITMAP
static ans_t __imp_clique3_bitmap(sysname_t tasklet_id, node_t second_index) {
    ans_t ans = 0;
    mram_read(mram_bitmap[second_index], bitmap[tasklet_id], sizeof(bitmap[tasklet_id]));
    for (node_t i = 0; i < bitmap_size; i++) {
        uint32_t tmp = bitmap[tasklet_id][i];
        if (tmp) for (node_t j = 0; j < 32; j++) {
            if (tmp & (1 << j)) ans++;
        }
    }
    return ans;
}
#endif

static ans_t __imp_clique3_2(sysname_t tasklet_id, node_t __mram_ptr * root_col, node_t root_size, node_t __mram_ptr * second_col, node_t second_size) {

    if(!second_size)return 0;
    
    node_t(*tasklet_buf)[BUF_SIZE] = buf[tasklet_id];


#ifdef NO_RUN   //test cycle without Intersection operation 
    //node_t ans =  intersect_seq_buf_thresh_no_run(tasklet_buf, root_col, root_size, second_col, second_size);
    //node_t ans = 1;
#else
    node_t ans =  intersect_seq_buf_thresh(tasklet_buf, root_col, root_size, second_col, second_size);
#endif

  
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

static ans_t __imp_clique3_partition(sysname_t tasklet_id, node_t root) {
    edge_ptr root_begin = row_ptr[root];  // intended DMA
    edge_ptr root_end = row_ptr[root + 1];  // intended DMA
    node_t root_size = root_end - root_begin;
    if(!root_size)return 0;
    ans_t ans = 0;
    for (edge_ptr i = root_begin + 1; i<root_end; i++) {
        node_t second_root = col_idx[i];  // intended DMA 
        edge_ptr second_begin = row_ptr[second_root];  // intended DMA
        edge_ptr second_end = row_ptr[second_root+1];  // intended DMA
        ans += __imp_clique3_2(tasklet_id,&col_idx[root_begin],i-root_begin,&col_idx[second_begin],second_end-second_begin);
    }
    return ans;
}


//func begin
extern void clique3( sysname_t tasklet_id )
{
	node_t i = 0;
	large_degree_num = root_num;                            /* if all node is large_degree */
	while ( i < root_num )
	{
		node_t	root		= roots[i];             /* intended DMA */
		node_t	root_begin	= row_ptr[root];        /* intended DMA */
		node_t	root_end	= row_ptr[root + 1];    /* intended DMA */
		node_t	root_size	= root_end - root_begin;
		if ( root_size < BRANCH_LEVEL_THRESHOLD )
		{
			large_degree_num = i;
			break;
		}

#ifdef PERF
		timer_start( &cycles[tasklet_id] );
#endif
		partial_ans[tasklet_id] = 0;

		if ( no_partition_flag )
		{
			// const int	num_dma_threads		= 8;
			// const node_t	max_edges_per_chunk	= 256;  /* 每次最多搬运 256 条边 = 512 元素 */

			// node_t *cb = col_buf[0];                        /* 共享缓冲区，容量 512 个元素 */

			// edge_ptr	root_begin	= row_ptr[root];
			// edge_ptr	root_end	= row_ptr[root + 1];
			// node_t		root_size	= root_end - root_begin;

			// for ( node_t chunk_offset = 1; chunk_offset < root_size; chunk_offset += max_edges_per_chunk )
			// {
			// 	node_t chunk_size = MIN( max_edges_per_chunk, root_size - chunk_offset ); /* 本次最多搬多少边 */

			// 	/* === 1. 并行搬运 === */
			// 	if ( tasklet_id < num_dma_threads )
			// 	{
			// 		node_t	local_chunk_size	= (chunk_size + num_dma_threads - 1) / num_dma_threads;
			// 		node_t	local_start		= tasklet_id * local_chunk_size;
			// 		node_t	local_end		= MIN( (tasklet_id + 1) * local_chunk_size, chunk_size );

			// 		if ( local_start < local_end )
			// 		{
			// 			edge_ptr	mram_src_offset = root_begin + chunk_offset + local_start;
			// 			node_t		len		= local_end - local_start;

			// 			mram_read( &col_idx[edge_offset + 2 * mram_src_offset], cb + 2 * local_start, len << (SIZE_EDGE_PTR_LOG + 1) ); /* len * 2 * sizeof(edge_ptr) */
			// 		}
			// 	}

			// 	barrier_wait( &co_barrier );                                                                                                    /* 所有搬运完成再进入计算 */

			// 	/* === 2. 并行处理 === */
			// 	for ( edge_ptr j = tasklet_id; j < chunk_size; j += NR_TASKLETS )
			// 	{
			// 		partial_ans[tasklet_id] += __imp_clique3_2(
			// 			tasklet_id,
			// 			&col_idx[root_begin],                                                                                           /* 原始 row 起点 */
			// 			chunk_offset + j,                                                                                               /* 真实全局偏移 j */
			// 			&col_idx[cb[2 * j]],
			// 			cb[2 * j + 1] - cb[2 * j]
			// 			);
			// 	}

			// 	barrier_wait( &co_barrier );                                                                                                    /* 等所有线程处理完，再搬下一块 */
			// }

			
			 //without prefetch
			 for ( edge_ptr j = root_begin + tasklet_id + 1; j < root_end; j += NR_TASKLETS )
			 {
			     partial_ans[tasklet_id] += __imp_clique3_2( tasklet_id, &col_idx[root_begin], j - root_begin, &col_idx[col_idx[edge_offset + 2 * j]], col_idx[edge_offset + 2 * j + 1] - col_idx[edge_offset + 2 * j] );
			 }
			 
		}else{
			for ( edge_ptr j = root_begin + tasklet_id + 1; j < root_end; j += NR_TASKLETS )
			{
				node_t		second_root	= col_idx[j];                   /* intended DMA */
				edge_ptr	second_begin	= row_ptr[second_root];         /* intended DMA */
				edge_ptr	second_end	= row_ptr[second_root + 1];     /* intended DMA */
				partial_ans[tasklet_id] += __imp_clique3_2( tasklet_id, &col_idx[root_begin], i - root_begin, &col_idx[second_begin], second_end - second_begin );
			}
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
		if ( no_partition_flag )
			ans[i] = __imp_clique3( tasklet_id, root );             /* intended DMA */
		else
			ans[i] = __imp_clique3_partition( tasklet_id, root );   /* intended DMA */

#ifdef PERF
		cycle_ct[i] = timer_stop( &cycles[tasklet_id] );                /* intended DMA */
#endif
	}
}