for (int blk = 0; blk < static_cast<std::ptrdiff_t>(world.num_blocks); blk++) {
fccz4_constraints_block(world.geom[blk], in_gfs, diagnostic_gfs);
} // END LOOP: for blk over [0, static_cast<std::ptrdiff_t>(world.num_blocks))
