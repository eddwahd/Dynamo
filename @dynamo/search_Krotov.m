function term_reason = search_Krotov()
global OC;

meta_term_cond = OC.opt.term_cond;

per_block_term_cond = OC.opt.term_cond;
per_block_term_cond.loop_count = 1;

block_size_1 = 1;

term_reason = search_block_cycle(OC.opt.control_mask, meta_term_cond, per_block_term_cond,  @search_first_order, block_size_1);
