for (int blk = 0; blk < static_cast<std::ptrdiff_t>(world.num_blocks); blk++) {
fccz4_enforce_detgbar_equals_detghat_trAzero_block(world.geom[blk], in_gfs, status);
} // END LOOP: for blk over [0, static_cast<std::ptrdiff_t>(world.num_blocks))
