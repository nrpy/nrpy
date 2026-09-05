for (int blk = 0; blk < static_cast<std::ptrdiff_t>(world.num_blocks); blk++) {
fccz4_project_block(world.geom[blk], y_n_gfs, status);
} // END LOOP: for blk over [0, static_cast<std::ptrdiff_t>(world.num_blocks))
