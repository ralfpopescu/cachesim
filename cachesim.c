#include "cachesim.h"
#include <stdio.h>

#define TRUE 1
#define FALSE 0

/**
 * The stuct that you may use to store the metadata for each block in the L1 and L2 caches
 */
typedef struct block_t {
    uint64_t tag; // The tag stored in that block
    uint8_t valid; // Valid bit
    uint8_t dirty; // Dirty bit

    /**************** TODO ******************/

    /*
        Add another variable here to perform the LRU replacement. Look into using a counter
        variable that will keep track of the oldest block in a set
    */
    uint64_t LRU;


} block;


/**
 * A struct for storing the configuration of both the L1 and L2 caches as passed in the
 * cache_init function. All values represent number of bits. You may add any parameters
 * here, however I strongly suggest not removing anything from the config struct
 */
typedef struct config_t {
    uint64_t C1; /* Size of cache L1 */
    uint64_t C2; /* Size of cache L2 */
    uint64_t S; /* Set associativity of L2 */
    uint64_t B; /* Block size of both caches */
} config;


/****** Do not modify the below function headers ******/
static uint64_t get_tag(uint64_t address, uint64_t C, uint64_t B, uint64_t S);
static uint64_t get_index(uint64_t address, uint64_t C, uint64_t B, uint64_t S);
static uint64_t convert_tag(uint64_t tag, uint64_t index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S);
static uint64_t convert_index(uint64_t tag, uint64_t index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S);
static uint64_t convert_tag_l1(uint64_t l2_tag, uint64_t l2_index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S);
static uint64_t convert_index_l1(uint64_t l2_tag, uint64_t l2_index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S);

/****** You may add Globals and other function headers that you may need below this line ******/

block* find_LRU (int L2_index);
block* find_empty_block(int L2_index);
void updatetimestamp(block* blk, struct cache_stats_t *stats);
block* findL2block(uint64_t L2_tag, int L2_index);



    config* con;
    block* L1;
    block** L2;

    int block_size;
    int L1_size;
    int L1_num_of_blocks;

    int L2_size;
    int L2_num_of_blocks;

    int num_of_blocks_in_set;
    int num_of_sets;



/**
 * Subroutine for initializing your cache with the passed in arguments.
 * You may initialize any globals you might need in this subroutine
 *
 * @param C1 The total size of your L1 cache is 2^C1 bytes
 * @param C2 The tatal size of your L2 cache is 2^C2 bytes
 * @param S The total number of blocks in a line/set of your L2 cache are 2^S
 * @param B The size of your blocks is 2^B bytes
 */
void cache_init(uint64_t C1, uint64_t C2, uint64_t S, uint64_t B)
{
    /* 
        Initialize the caches here. I strongly suggest using arrays for representing
        meta data stored in the caches. The block_t struct given above maybe useful
    */

    /**************** TODO ******************/
    config* conf = (config*)malloc(sizeof(config)); 
    conf->C1 = C1;
    conf->C2 = C2;
    conf->S = S;
    conf->B = B;
    con = conf;



    block_size = 2 << B;
    L1_size = 2 << C1;
    L1_num_of_blocks = L1_size / block_size;

    L2_size = 2 << C2;
    L2_num_of_blocks = L2_size / block_size;

    num_of_blocks_in_set = 2 << S;
    num_of_sets = L2_num_of_blocks / num_of_blocks_in_set; 

    //should be a simple array of blocks
    L1 = malloc(L1_num_of_blocks * sizeof(block));

    for (int j = 0; j < L1_num_of_blocks; j++) {
        block* blk = &(L1[j]);
        blk->valid = FALSE;
        blk->dirty = FALSE;
        blk->tag = 0;
    }

    //L2 should be an array of sets, which are arrays of blocks
    L2 = malloc(num_of_sets * sizeof(block*));

    for(int i = 0; i < num_of_sets; i++) {
        L2[i] = (block*)malloc(sizeof(block)*num_of_blocks_in_set);
        block* set_ptr = L2[i];
 
        for(int k = 0; k < num_of_blocks_in_set; k++) {
            block* blk3 = &(set_ptr[k]);
            blk3->valid = FALSE;
            blk3->dirty = FALSE;
            blk3->tag = 0;
        }
    }

    //then initialize all blocks to zero
    

}

/**
 * Subroutine that simulates one cache event at a time.
 * @param rw The type of access, READ or WRITE
 * @param address The address that is being accessed
 * @param stats The struct that you are supposed to store the stats in
 */
void cache_access (char rw, uint64_t address, struct cache_stats_t *stats)
{
    /* 
        Suggested approach:
            -> Find the L1 tag and index of the address that is being passed in to the function
            -> Check if there is a hit in the L1 cache
            -> If L1 misses, check the L2 cache for a hit (Hint: If L2 hits, update L1 with new values)
            -> If L2 misses, need to get values from memory, and update L2 and L1 caches
            
            * We will leave it upto you to decide what must be updated and when
     */

    /**************** TODO ******************/
    
    stats->accesses = stats->accesses + 1;
    int hit = FALSE;
    //get tags and indices
    uint64_t L1_tag = get_tag(address, con->C1, con->B, 0);
    int L1_index = get_index(address, con->C1, con->B, 0);

    uint64_t L2_tag = convert_tag(L1_tag, L1_index, con->C1, con->C2, con->B, con->S);
    int L2_index = convert_index(L1_tag, L1_index, con->C1, con->C2, con->B, con->S);

    //int set_contained;
    printf("%i", L1[L1_index].tag);
    //check if hit in L1
        if (L1[L1_index].tag == L1_tag)
        {
            hit = TRUE;
            printf("hit");
        }
   
    block* L1_block = &(L1[L1_index]);

    int valid = L1[L1_index].valid;

    if(rw == READ)
    {
        if (valid == FALSE) {
            hit = FALSE;
        }

        if(hit)
        {
        //update L2's time stamp
        block* L2block = findL2block(L2_tag, L2_index);
        updatetimestamp(L2block, stats);
        //update stats
        stats->reads = stats->reads + 1;
        }
        else 
        {
        //check L2
            stats->l1_read_misses = stats->l1_read_misses + 1;

            L1_block->tag = L1_tag;

            hit = FALSE;         
            //check L2 for hit
            block* L2blk = findL2block(L2_tag, L2_index);
            if(L2blk != 0){
                hit = TRUE;
            }

            if(hit)
            {
            //increase stats
            stats->reads = stats->reads + 1;
            //pull into L1
            L1_block->tag = L1_tag;

            if(L2blk->dirty == TRUE){
                updateblk->dirty = FALSE;
            }
            updateblk->dirty = FALSE;

            } 
            else 
            {
            stats->l2_read_misses = stats->l2_read_misses + 1;

            //update block
            block* updateblock = &(L1[L1_index]);
            updateblock->tag = L1_tag;

            block* emptyblock = find_empty_block(L2_index);
            block* LRU;

            if(emptyblock == 0){
            LRU = find_LRU(L2_index);
            //evict
            LRU->tag = L2_tag;
            if(LRU->dirty == TRUE) {
                stats->write_backs = stats->write_backs + 1;
                LRU->dirty = FALSE;
            }
            int L1_back_index = convert_index_l1(L2_tag, L2_index, con->C1, con->C2, con->B, con->S);
            uint64_t L1_back_tag = convert_tag_l1(L2_tag, L2_index, con->C1, con->C2, con->B, con->S);
            block* LRU_L1 = &(L1[L1_back_index]);
            if(LRU_L1->valid == FALSE) {
                LRU->tag = L1_back_tag;
                LRU_L1->valid = TRUE;
            } else {
                if(LRU_L1->dirty == TRUE){
                    //evict
                    LRU->tag = L1_back_tag;
                    LRU_L1->valid = TRUE;
                } else {
                LRU->tag = L1_back_tag;
                LRU_L1->valid = TRUE;
                }
            }

            } else {
            //found empty block
                emptyblock->valid = TRUE;
                emptyblock->tag = L2_tag;
            }
            
            
            }
        }
    
    }

    if(rw == WRITE){

        if (valid == FALSE) {
            hit = FALSE;
        }

        if(hit) {
        stats->writes = stats->writes + 1;
        block* blk = &(L1[L1_index]);
        block* L2blk = findL2block(L2_tag, L2_index);
        updatetimestamp(L2blk, stats);
        blk->dirty = TRUE;   
        }
        //if missed L1;
        else {
        stats->l1_write_misses = stats->l1_write_misses + 1;

        int L2hit = FALSE;
        //check if hit in L2
        block* blk = &(L1[L1_index]);
        block* L2blk = findL2block(L2_tag, L2_index);
        if(L2blk != 0) {
            L2hit = TRUE;
            updatetimestamp(L2blk, stats);
        }

        //if L2 hit
        if(L2hit){
            stats->writes = stats->writes + 1;
            blk->tag = L1_tag;

            L2blk->dirty = TRUE;

        //if L2 didn't hit
        } else {
        stats->l2_write_misses = stats->l2_write_misses + 1;

        //find LRU
            block* emptyblock = find_empty_block(L2_index);
            block* LRU;

            if(emptyblock == 0){
            LRU = find_LRU(L2_index);
            //evict LRU
            LRU->tag = L2_tag;
            LRU->dirty = TRUE;

            L1_block->dirty = FALSE;
            
            } else {
            //we found invalid block
            emptyblock->valid = TRUE;
            emptyblock->tag = L2_tag;
            emptyblock->dirty = TRUE;
            }
            
            
           
        
        }
        
        }
    }
    

    
}

block* find_LRU (int L2_index){
            block* emptyblock;
            block* LRU;
            int LRUcount = 99999999;
            int LRUset = 0;
            //checks for empty block
            for(int j = 0; j < num_of_sets; j++)
            {
                block* set = L2[j];
                block* checkblock = &(set[L2_index]);

                if(checkblock->LRU < LRUcount){
                    LRU = checkblock;
                    LRUcount = checkblock->LRU;
                    LRUset = j;
                }
            }
            return LRU;
}

block* find_empty_block(int L2_index){
            block* emptyblock;
            block* LRU;
            int emptyblockset = 0;
            //checks for empty block
            for(int j = 0; j < num_of_sets; j++)
            {
                block* set = L2[j];
                block* checkblock = &(set[L2_index]);
         
                if(checkblock->valid == FALSE)
                {
                    emptyblock = checkblock;
                    emptyblockset = j;
                }

            }
    return emptyblock;
}

void updatetimestamp(block* blk, struct cache_stats_t *stats){
    if(blk != 0){
    blk->LRU = stats->accesses;
    }
}

block* findL2block(uint64_t L2_tag, int L2_index) {
    block* l2blk = 0;
    for(int m = 0; m < num_of_sets; m++)
        {
            block* set = L2[m];
            block* updateblock = &(set[L2_index]);
            if(updateblock->tag == L2_tag)
            {
                l2blk = updateblock;
            }

        }
    return l2blk;
}

/**
 * Subroutine for freeing up memory, and performing any final calculations before the statistics
 * are outputed by the driver
 */
void cache_cleanup (struct cache_stats_t *stats)
{
    /*
        Make sure to free up all the memory you malloc'ed here. To check if you have freed up the
        the memory, run valgrind. For more information, google how to use valgrind.
    */

    /**************** TODO ****************/
    free(L1);
    free(L2);

    stats->read_misses = stats->l1_read_misses + stats->l2_read_misses;
    stats->write_misses = stats->l1_write_misses + stats->l2_write_misses;
    stats->misses = stats->read_misses + stats->write_misses;

    stats->l1_miss_rate = stats->reads;
    stats->l2_miss_rate = stats->reads;
    stats->miss_rate = stats->reads;
    //calculate miss rates and averages
}

/**
 * Subroutine to compute the Tag of a given address based on the parameters passed in
 *
 * @param address The address whose tag is to be computed
 * @param C The size of the cache in bytes (i.e. Size of cache is 2^C)
 * @param B The size of the cache block in bytes (i.e. Size of block is 2^B)
 * @param S The set associativity of the cache in bytes (i.e. Set-Associativity is 2^S)
 * 
 * @return The computed tag
 */
static uint64_t get_tag(uint64_t address, uint64_t C, uint64_t B, uint64_t S)
{
    /**************** TODO ******************/
    uint64_t tag;

    int num_of_byteoffset_bits = B;
    int num_of_index_bits = C - B - S;
    //C - (B + S) for L2
    int num_of_tag_bits = 64 - num_of_byteoffset_bits - num_of_index_bits;

    tag = address >> (64 - num_of_tag_bits);


    return tag;
}

/**
 * Subroutine to compute the Index of a given address based on the parameters passed in
 *
 * @param address The address whose tag is to be computed
 * @param C The size of the cache in bits (i.e. Size of cache is 2^C)
 * @param B The size of the cache block in bits (i.e. Size of block is 2^B)
 * @param S The set associativity of the cache in bits (i.e. Set-Associativity is 2^S)
 *
 * @return The computed index
 */
static uint64_t get_index(uint64_t address, uint64_t C, uint64_t B, uint64_t S)
{
    /**************** TODO ******************/
    uint64_t index;

    
    int num_of_byteoffset_bits = B;
    int num_of_index_bits = C - B - S;
    //int num_of_tag_bits = 64 - num_of_byteoffset_bits - num_of_index_bits;

    index = address >> B;
    int mask = (1 << (num_of_index_bits)) - 1;
    index = index & mask;

    //index = address % cachesize_in_bytes

    //L2_index = address % set_size

    return index;
}


/**** DO NOT MODIFY CODE BELOW THIS LINE UNLESS YOU ARE ABSOLUTELY SURE OF WHAT YOU ARE DOING ****/

/*
    Note:   The below functions will be useful in converting the L1 tag and index into corresponding L2
            tag and index. These should be used when you are evicitng a block from the L1 cache, and
            you need to update the block in L2 cache that corresponds to the evicted block.

            The newly added functions will be useful for converting L2 indecies ang tags into the corresponding
            L1 index and tags. Make sure to understand how they are working.
*/

/**
 * This function converts the tag stored in an L1 block and the index of that L1 block into corresponding
 * tag of the L2 block
 *
 * @param tag The tag that needs to be converted (i.e. L1 tag)
 * @param index The index of the L1 cache (i.e. The index from which the tag was found)
 * @param C1 The size of the L1 cache in bits
 * @param C2 The size of the l2 cache in bits
 * @param B The size of the block in bits
 * @param S The set associativity of the L2 cache
 */
static uint64_t convert_tag(uint64_t tag, uint64_t index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S)
{
    uint64_t reconstructed_address = (tag << (C1 - B)) | index;
    return reconstructed_address >> (C2 - B - S);
}

/**
 * This function converts the tag stored in an L1 block and the index of that L1 block into corresponding
 * index of the L2 block
 *
 * @param tag The tag stored in the L1 index
 * @param index The index of the L1 cache (i.e. The index from which the tag was found)
 * @param C1 The size of the L1 cache in bits
 * @param C2 The size of the l2 cache in bits
 * @param B The size of the block in bits
 * @param S The set associativity of the L2 cache
 */
static uint64_t convert_index(uint64_t tag, uint64_t index, uint64_t C1, uint64_t C2, uint64_t B,  uint64_t S)
{
    // Reconstructed address without the block offset bits
    uint64_t reconstructed_address = (tag << (C1 - B)) | index;
    // Create index mask for L2 without including the block offset bits
    return reconstructed_address & ((1 << (C2 - S - B)) - 1);
}

/**
 * This function converts the tag stored in an L2 block and the index of that L2 block into corresponding
 * tag of the L1 cache
 *
 * @param l2_tag The L2 tag
 * @param l2_index The index of the L2 block
 * @param C1 The size of the L1 cache in bits
 * @param C2 The size of the l2 cache in bits
 * @param B The size of the block in bits
 * @param S The set associativity of the L2 cache
 * @return The L1 tag linked to the L2 index and tag
 */
static uint64_t convert_tag_l1(uint64_t l2_tag, uint64_t l2_index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S) {
    uint64_t reconstructed_address = (l2_tag << (C2 - B - S)) | l2_index;
    return reconstructed_address >> (C1 - B);
}

/**
 * This function converts the tag stored in an L2 block and the index of that L2 block into corresponding
 * index of the L1 block
 *
 * @param l2_tag The L2 tag
 * @param l2_index The index of the L2 block
 * @param C1 The size of the L1 cache in bits
 * @param C2 The size of the l2 cache in bits
 * @param B The size of the block in bits
 * @param S The set associativity of the L2 cache
 * @return The L1 index of the L2 block
 */
static uint64_t convert_index_l1(uint64_t l2_tag, uint64_t l2_index, uint64_t C1, uint64_t C2, uint64_t B, uint64_t S) {
    uint64_t reconstructed_address = (l2_tag << (C2 - B - S)) | l2_index;
    return reconstructed_address & ((1 << (C1 - B)) - 1);
}
