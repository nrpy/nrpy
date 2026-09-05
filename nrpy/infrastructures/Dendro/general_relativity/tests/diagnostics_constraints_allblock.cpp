for (int blk_id = 0; blk_id < static_cast<std::ptrdiff_t>(world.num_blocks); blk_id++) {
fccz4_constraints_block(world.geom[blk_id], in_gfs, diag_gfs);
} // END LOOP: for blk_id over [0, static_cast<std::ptrdiff_t>(world.num_blocks))
