#include <common.h>
#include <mram.h>
#include <defs.h>


extern node_t intersect_seq_buf_thresh(node_t (*buf)[BUF_SIZE], node_t __mram_ptr *a, node_t a_size, node_t __mram_ptr *b, node_t b_size, node_t threshold) {
    if(a_size==0||b_size==0)return 0;
        
    if (a_size > b_size) {
        node_t __mram_ptr *tmp = a;
        a = b;
        b = tmp;
        node_t tmp_size = a_size;
        a_size = b_size;
        b_size = tmp_size;
    }
    if (a_size < (b_size >> 4)&& a_size < BUF_SIZE) {
        node_t ans = 0;
        node_t i = 0;

        if (((uint64_t)a) & 4) {
            a--;
            i = 1; 
            a_size++;
        }
        
        mram_read(a, buf[0], ALIGN8((a_size)<<SIZE_NODE_T_LOG));
        for (; i < a_size; i++) { 
            node_t a_val = buf[0][i];
            if (a_val >= threshold) break;
            node_t l = 0, r = b_size;
            while (l <r) {
                node_t mid = l + ((r - l) >> 1);;
                node_t b_val = b[mid];  // intended DMA
                if (a_val == b_val) {
                    ans++;
                    break;
                }
                else if (a_val < b_val) {
                    r = mid;
                }
                else {
                    l = mid + 1;
                }
            }
        }
        return ans;
    }
    node_t *a_buf = buf[0];
    node_t *b_buf = buf[1];
    node_t i = 0, j = 0, ans = 0;
    if (((uint64_t)a) & 4) {
        a--;
        i = 1;
        a_size++;
    }
    if (((uint64_t)b) & 4) {
        b--;
        j = 1;
        b_size++;
    }
    mram_read(a, a_buf, ALIGN8(MIN(a_size, BUF_SIZE) << SIZE_NODE_T_LOG));
    mram_read(b, b_buf, ALIGN8(MIN(b_size, BUF_SIZE) << SIZE_NODE_T_LOG));

    while (i < a_size && j < b_size) {
        if (i == BUF_SIZE) {
            a_size -= i;
            a += i;
            mram_read(a, a_buf, ALIGN8(MIN(a_size, BUF_SIZE) << SIZE_NODE_T_LOG));
            i = 0;
        }
        if (j == BUF_SIZE) {
            b_size -= j;
            b += j;
            mram_read(b, b_buf, ALIGN8(MIN(b_size, BUF_SIZE) << SIZE_NODE_T_LOG));
            j = 0;
        }

        if (a_buf[i] >= threshold || b_buf[j] >= threshold) break;

        if (a_buf[i] == b_buf[j]) {
            ans++;
            i++;
            j++;
        }
        else if (a_buf[i] < b_buf[j]) {
            i++;
        }
        else {
            j++;
        }
    }
    return ans;
}

extern node_t intersect_seq_buf_thresh_no_run(node_t (*buf)[BUF_SIZE], node_t __mram_ptr *a, node_t a_size, node_t __mram_ptr *b, node_t b_size, node_t threshold) {
    node_t ans = 0;
    if (a_size > b_size) {
        node_t __mram_ptr *tmp = a;
        a = b;
        b = tmp;
        node_t tmp_size = a_size;
        a_size = b_size;
        b_size = tmp_size;
    }
    if (a_size < (b_size >> 4)&& a_size < BUF_SIZE) {

        node_t i = 0;

        if (((uint64_t)a) & 4) {
            a--;
            i = 1; 
            a_size++;
        }
        
        mram_read(a, buf[0], ALIGN8((a_size)<<SIZE_NODE_T_LOG));
        // for (; i < a_size; i++) { 
        //     node_t a_val = buf[0][i];
        //     if (a_val >= threshold) break;
        //     node_t l = 0, r = b_size;
        //     while (l <r) {
        //         node_t mid = l + ((r - l) >> 1);;
        //         node_t b_val = b[mid];  // intended DMA
        //         if (a_val == b_val) {
        //             ans++;
        //             break;
        //         }
        //         else if (a_val < b_val) {
        //             r = mid;
        //         }
        //         else {
        //             l = mid + 1;
        //         }
        //     }
        // }
        return ans;
    }
    node_t *a_buf = buf[0];
    node_t *b_buf = buf[1];
    node_t i = 0, j = 0;
    if (((uint64_t)a) & 4) {
        a--;
        i = 1;
        a_size++;
    }
    if (((uint64_t)b) & 4) {
        b--;
        j = 1;
        b_size++;
    }
    mram_read(a, a_buf, ALIGN8(MIN(a_size, BUF_SIZE) << SIZE_NODE_T_LOG));
    mram_read(b, b_buf, ALIGN8(MIN(b_size, BUF_SIZE) << SIZE_NODE_T_LOG));

    while (i < a_size && j < b_size) {
        if (i == BUF_SIZE) {
            a_size -= i;
            a += i;
            mram_read(a, a_buf, ALIGN8(MIN(a_size, BUF_SIZE) << SIZE_NODE_T_LOG));
            i = 0;
        }
        if (j == BUF_SIZE) {
            b_size -= j;
            b += j;
            mram_read(b, b_buf, ALIGN8(MIN(b_size, BUF_SIZE) << SIZE_NODE_T_LOG));
            j = 0;
        }

        if (a_buf[i] >= threshold || b_buf[j] >= threshold) break;

        i++;j++;
    }
    return ans;
}

#ifdef BITMAP
extern void intersect_bitmap(node_t *a, node_t *b, node_t *c, node_t bitmap_size) {
    for (node_t i = 0; i < bitmap_size; i++) {
        c[i] = a[i] & b[i];
    }
}
#endif

