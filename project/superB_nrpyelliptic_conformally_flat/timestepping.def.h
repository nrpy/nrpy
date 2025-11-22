




/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::start_6_closure : public SDAG::Closure {
      

      start_6_closure() {
        init();
      }
      start_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~start_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(start_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::diagnostics_ckio_19_closure : public SDAG::Closure {
            Ck::IO::Session token;
            int which_diagnostics_part;


      diagnostics_ckio_19_closure() {
        init();
      }
      diagnostics_ckio_19_closure(CkMigrateMessage*) {
        init();
      }
            Ck::IO::Session & getP0() { return token;}
            int & getP1() { return which_diagnostics_part;}
      void pup(PUP::er& __p) {
        __p | token;
        __p | which_diagnostics_part;
        packClosure(__p);
      }
      virtual ~diagnostics_ckio_19_closure() {
      }
      PUPable_decl(SINGLE_ARG(diagnostics_ckio_19_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::east_ghost_20_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      east_ghost_20_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      east_ghost_20_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~east_ghost_20_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(east_ghost_20_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::west_ghost_21_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      west_ghost_21_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      west_ghost_21_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~west_ghost_21_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(west_ghost_21_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::north_ghost_22_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      north_ghost_22_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      north_ghost_22_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~north_ghost_22_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(north_ghost_22_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::south_ghost_23_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      south_ghost_23_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      south_ghost_23_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~south_ghost_23_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(south_ghost_23_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::top_ghost_24_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      top_ghost_24_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      top_ghost_24_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~top_ghost_24_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(top_ghost_24_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::bottom_ghost_25_closure : public SDAG::Closure {
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      bottom_ghost_25_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      bottom_ghost_25_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return type_gfs;}
            int & getP1() { return len_tmpBuffer;}
            REAL *& getP2() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~bottom_ghost_25_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(bottom_ghost_25_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::continue_timestepping_26_closure : public SDAG::Closure {
      

      continue_timestepping_26_closure() {
        init();
      }
      continue_timestepping_26_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~continue_timestepping_26_closure() {
      }
      PUPable_decl(SINGLE_ARG(continue_timestepping_26_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure : public SDAG::Closure {
            int idx3_of_sendingchare;
            int num_srcpts;
            int *globalidx3_srcpts;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return idx3_of_sendingchare;}
            int & getP1() { return num_srcpts;}
            int *& getP2() { return globalidx3_srcpts;}
      void pup(PUP::er& __p) {
        __p | idx3_of_sendingchare;
        __p | num_srcpts;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> idx3_of_sendingchare;
  implP|idx3_of_sendingchare;
  PUP::detail::TemporaryObjectHolder<int> num_srcpts;
  implP|num_srcpts;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  implP|impl_off_globalidx3_srcpts;
  implP|impl_cnt_globalidx3_srcpts;
          impl_buf+=CK_ALIGN(implP.size(),16);
          globalidx3_srcpts = (int *)(impl_buf+impl_off_globalidx3_srcpts);
        }
      }
      virtual ~receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_k_even_gfs_30_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_k_even_gfs_30_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_k_even_gfs_30_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_k_even_gfs_30_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_y_n_gfs_31_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_y_n_gfs_31_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_y_n_gfs_31_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_y_n_gfs_31_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure : public SDAG::Closure {
            int src_chare_idx3;
            int type_gfs;
            int len_tmpBuffer;
            REAL *tmpBuffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return src_chare_idx3;}
            int & getP1() { return type_gfs;}
            int & getP2() { return len_tmpBuffer;}
            REAL *& getP3() { return tmpBuffer;}
      void pup(PUP::er& __p) {
        __p | src_chare_idx3;
        __p | type_gfs;
        __p | len_tmpBuffer;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
        }
      }
      virtual ~receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure : public SDAG::Closure {
            REAL wavespeed_at_outer_boundary;


      receiv_wavespeed_at_outer_boundary_36_closure() {
        init();
      }
      receiv_wavespeed_at_outer_boundary_36_closure(CkMigrateMessage*) {
        init();
      }
            REAL & getP0() { return wavespeed_at_outer_boundary;}
      void pup(PUP::er& __p) {
        __p | wavespeed_at_outer_boundary;
        packClosure(__p);
      }
      virtual ~receiv_wavespeed_at_outer_boundary_36_closure() {
      }
      PUPable_decl(SINGLE_ARG(receiv_wavespeed_at_outer_boundary_36_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */







/* DEFS: array Timestepping: ArrayElement{
Timestepping(const CommondataObject &inData);
void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
void start();
void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
void test_written_1d_y(CkReductionMsg* impl_msg);
void test_written_1d_z(CkReductionMsg* impl_msg);
void test_written_2d_xy(CkReductionMsg* impl_msg);
void test_written_2d_yz(CkReductionMsg* impl_msg);
void closed_1d_y(CkReductionMsg* impl_msg);
void closed_1d_z(CkReductionMsg* impl_msg);
void closed_2d_xy(CkReductionMsg* impl_msg);
void closed_2d_yz(CkReductionMsg* impl_msg);
void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void continue_timestepping();
void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
Timestepping(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Timestepping::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_Timestepping::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_Timestepping::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_Timestepping::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_Timestepping::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_Timestepping::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_Timestepping::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(const CommondataObject &inData);
 */
void CProxyElement_Timestepping::insert(const CommondataObject &inData, int onPE, const CkEntryOptions *impl_e_opts)
{ 
   //Marshall: const CommondataObject &inData
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
  }
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Timestepping::idx_Timestepping_marshall1(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::ready_1d_y(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_1d_y_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::ready_1d_z(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_1d_z_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_2d_xy_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_2d_yz_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start();
 */
void CProxyElement_Timestepping::start(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_y_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_z_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_xy_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxyElement_Timestepping::start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_yz_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_y(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::test_written_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_y_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_z(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::test_written_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_z_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::test_written_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_xy_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::test_written_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_yz_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_y(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::closed_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_1d_y_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_z(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::closed_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_1d_z_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::closed_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_2d_xy_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxyElement_Timestepping::closed_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_2d_yz_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
 */
void CProxyElement_Timestepping::diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Ck::IO::Session &token, int which_diagnostics_part
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_diagnostics_ckio_marshall19(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_east_ghost_marshall20(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_west_ghost_marshall21(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_north_ghost_marshall22(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_south_ghost_marshall23(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_top_ghost_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_bottom_ghost_marshall25(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void continue_timestepping();
 */
void CProxyElement_Timestepping::continue_timestepping(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_continue_timestepping_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  impl_off_globalidx3_srcpts=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_globalidx3_srcpts=sizeof(int)*(num_srcpts));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_globalidx3_srcpts,globalidx3_srcpts,impl_cnt_globalidx3_srcpts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxyElement_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
 */
void CProxyElement_Timestepping::receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const REAL &wavespeed_at_outer_boundary
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_wavespeed_at_outer_boundary_marshall36(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(const CommondataObject &inData);
 */
CkArrayID CProxy_Timestepping::ckNew(const CommondataObject &inData, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const CommondataObject &inData
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_Timestepping::idx_Timestepping_marshall1(), opts);
  return gId;
}
void CProxy_Timestepping::ckNew(const CommondataObject &inData, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const CommondataObject &inData
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_Timestepping::idx_Timestepping_marshall1(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_Timestepping::ckNew(const CommondataObject &inData, const int s1, const int s2, const int s3, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const CommondataObject &inData
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
  }
  CkArrayOptions opts(s1, s2, s3);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_Timestepping::idx_Timestepping_marshall1(), opts);
  return gId;
}
void CProxy_Timestepping::ckNew(const CommondataObject &inData, const int s1, const int s2, const int s3, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const CommondataObject &inData
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CommondataObject>::type>::type &)inData;
  }
  CkArrayOptions opts(s1, s2, s3);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_Timestepping::idx_Timestepping_marshall1(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_Timestepping::reg_Timestepping_marshall1() {
  int epidx = CkRegisterEp("Timestepping(const CommondataObject &inData)",
      reinterpret_cast<CkCallFnPtr>(_call_Timestepping_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_Timestepping_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_Timestepping_marshall1);

  return epidx;
}

void CkIndex_Timestepping::_call_Timestepping_marshall1(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CommondataObject &inData*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CommondataObject> inData;
  implP|inData;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) Timestepping(std::move(inData.t));
}
int CkIndex_Timestepping::_callmarshall_Timestepping_marshall1(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const CommondataObject &inData*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CommondataObject> inData;
  implP|inData;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) Timestepping(std::move(inData.t));
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_Timestepping_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CommondataObject &inData*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CommondataObject> inData;
  implP|inData;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("inData");
  implDestP|inData;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxy_Timestepping::ready_1d_y(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_ready_1d_y_FileReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_ready_1d_y_FileReadyMsg() {
  int epidx = CkRegisterEp("ready_1d_y(Ck::IO::FileReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_ready_1d_y_FileReadyMsg), Ck::IO::CMessage_FileReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::FileReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_ready_1d_y_FileReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->ready_1d_y((Ck::IO::FileReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxy_Timestepping::ready_1d_z(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_ready_1d_z_FileReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_ready_1d_z_FileReadyMsg() {
  int epidx = CkRegisterEp("ready_1d_z(Ck::IO::FileReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_ready_1d_z_FileReadyMsg), Ck::IO::CMessage_FileReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::FileReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_ready_1d_z_FileReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->ready_1d_z((Ck::IO::FileReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxy_Timestepping::ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_ready_2d_xy_FileReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_ready_2d_xy_FileReadyMsg() {
  int epidx = CkRegisterEp("ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_ready_2d_xy_FileReadyMsg), Ck::IO::CMessage_FileReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::FileReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_ready_2d_xy_FileReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->ready_2d_xy((Ck::IO::FileReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxy_Timestepping::ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_ready_2d_yz_FileReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_ready_2d_yz_FileReadyMsg() {
  int epidx = CkRegisterEp("ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_ready_2d_yz_FileReadyMsg), Ck::IO::CMessage_FileReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::FileReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_ready_2d_yz_FileReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->ready_2d_yz((Ck::IO::FileReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start();
 */
void CProxy_Timestepping::start(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_start_void(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_start_void() {
  int epidx = CkRegisterEp("start()",
      reinterpret_cast<CkCallFnPtr>(_call_start_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Timestepping::_call_start_void(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->_sdag_fnc_start();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::start_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxy_Timestepping::start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_y_SessionReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_start_write_1d_y_SessionReadyMsg() {
  int epidx = CkRegisterEp("start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_start_write_1d_y_SessionReadyMsg), Ck::IO::CMessage_SessionReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::SessionReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_start_write_1d_y_SessionReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->start_write_1d_y((Ck::IO::SessionReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxy_Timestepping::start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_z_SessionReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_start_write_1d_z_SessionReadyMsg() {
  int epidx = CkRegisterEp("start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_start_write_1d_z_SessionReadyMsg), Ck::IO::CMessage_SessionReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::SessionReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_start_write_1d_z_SessionReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->start_write_1d_z((Ck::IO::SessionReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxy_Timestepping::start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_xy_SessionReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_start_write_2d_xy_SessionReadyMsg() {
  int epidx = CkRegisterEp("start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_start_write_2d_xy_SessionReadyMsg), Ck::IO::CMessage_SessionReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::SessionReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_start_write_2d_xy_SessionReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->start_write_2d_xy((Ck::IO::SessionReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxy_Timestepping::start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_yz_SessionReadyMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_start_write_2d_yz_SessionReadyMsg() {
  int epidx = CkRegisterEp("start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_start_write_2d_yz_SessionReadyMsg), Ck::IO::CMessage_SessionReadyMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)Ck::IO::SessionReadyMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_start_write_2d_yz_SessionReadyMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->start_write_2d_yz((Ck::IO::SessionReadyMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_y(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::test_written_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_y_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_test_written_1d_y_CkReductionMsg() {
  int epidx = CkRegisterEp("test_written_1d_y(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_test_written_1d_y_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_test_written_1d_y_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->test_written_1d_y((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_z(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::test_written_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_z_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_test_written_1d_z_CkReductionMsg() {
  int epidx = CkRegisterEp("test_written_1d_z(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_test_written_1d_z_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_test_written_1d_z_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->test_written_1d_z((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::test_written_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_xy_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_test_written_2d_xy_CkReductionMsg() {
  int epidx = CkRegisterEp("test_written_2d_xy(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_test_written_2d_xy_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_test_written_2d_xy_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->test_written_2d_xy((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::test_written_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_yz_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_test_written_2d_yz_CkReductionMsg() {
  int epidx = CkRegisterEp("test_written_2d_yz(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_test_written_2d_yz_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_test_written_2d_yz_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->test_written_2d_yz((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_y(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::closed_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_closed_1d_y_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_closed_1d_y_CkReductionMsg() {
  int epidx = CkRegisterEp("closed_1d_y(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_closed_1d_y_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_closed_1d_y_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->closed_1d_y((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_z(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::closed_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_closed_1d_z_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_closed_1d_z_CkReductionMsg() {
  int epidx = CkRegisterEp("closed_1d_z(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_closed_1d_z_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_closed_1d_z_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->closed_1d_z((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::closed_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_closed_2d_xy_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_closed_2d_xy_CkReductionMsg() {
  int epidx = CkRegisterEp("closed_2d_xy(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_closed_2d_xy_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_closed_2d_xy_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->closed_2d_xy((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxy_Timestepping::closed_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_closed_2d_yz_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_closed_2d_yz_CkReductionMsg() {
  int epidx = CkRegisterEp("closed_2d_yz(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_closed_2d_yz_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Timestepping::_call_closed_2d_yz_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->closed_2d_yz((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
 */
void CProxy_Timestepping::diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Ck::IO::Session &token, int which_diagnostics_part
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_diagnostics_ckio_marshall19(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_diagnostics_ckio_marshall19() {
  int epidx = CkRegisterEp("diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part)",
      reinterpret_cast<CkCallFnPtr>(_call_diagnostics_ckio_marshall19), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_diagnostics_ckio_marshall19);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_diagnostics_ckio_marshall19);

  return epidx;
}

void CkIndex_Timestepping::_call_diagnostics_ckio_marshall19(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::diagnostics_ckio_19_closure* genClosure = new Closure_Timestepping::diagnostics_ckio_19_closure();
  implP|genClosure->token;
  implP|genClosure->which_diagnostics_part;
  impl_buf+=CK_ALIGN(implP.size(),16);
  impl_obj->diagnostics_ckio(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_diagnostics_ckio_marshall19(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::diagnostics_ckio_19_closure* genClosure = new Closure_Timestepping::diagnostics_ckio_19_closure();
  implP|genClosure->token;
  implP|genClosure->which_diagnostics_part;
  impl_buf+=CK_ALIGN(implP.size(),16);
  impl_obj->diagnostics_ckio(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_diagnostics_ckio_marshall19(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Ck::IO::Session &token, int which_diagnostics_part*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Ck::IO::Session> token;
  implP|token;
  PUP::detail::TemporaryObjectHolder<int> which_diagnostics_part;
  implP|which_diagnostics_part;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("token");
  implDestP|token;
  if (implDestP.hasComments()) implDestP.comment("which_diagnostics_part");
  implDestP|which_diagnostics_part;
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::diagnostics_ckio_19_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_east_ghost_marshall20(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_east_ghost_marshall20() {
  int epidx = CkRegisterEp("east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_east_ghost_marshall20), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_east_ghost_marshall20);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_east_ghost_marshall20);

  return epidx;
}

void CkIndex_Timestepping::_call_east_ghost_marshall20(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::east_ghost_20_closure* genClosure = new Closure_Timestepping::east_ghost_20_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->east_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_east_ghost_marshall20(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::east_ghost_20_closure* genClosure = new Closure_Timestepping::east_ghost_20_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->east_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_east_ghost_marshall20(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::east_ghost_20_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_west_ghost_marshall21(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_west_ghost_marshall21() {
  int epidx = CkRegisterEp("west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_west_ghost_marshall21), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_west_ghost_marshall21);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_west_ghost_marshall21);

  return epidx;
}

void CkIndex_Timestepping::_call_west_ghost_marshall21(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::west_ghost_21_closure* genClosure = new Closure_Timestepping::west_ghost_21_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->west_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_west_ghost_marshall21(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::west_ghost_21_closure* genClosure = new Closure_Timestepping::west_ghost_21_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->west_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_west_ghost_marshall21(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::west_ghost_21_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_north_ghost_marshall22(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_north_ghost_marshall22() {
  int epidx = CkRegisterEp("north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_north_ghost_marshall22), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_north_ghost_marshall22);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_north_ghost_marshall22);

  return epidx;
}

void CkIndex_Timestepping::_call_north_ghost_marshall22(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::north_ghost_22_closure* genClosure = new Closure_Timestepping::north_ghost_22_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->north_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_north_ghost_marshall22(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::north_ghost_22_closure* genClosure = new Closure_Timestepping::north_ghost_22_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->north_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_north_ghost_marshall22(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::north_ghost_22_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_south_ghost_marshall23(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_south_ghost_marshall23() {
  int epidx = CkRegisterEp("south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_south_ghost_marshall23), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_south_ghost_marshall23);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_south_ghost_marshall23);

  return epidx;
}

void CkIndex_Timestepping::_call_south_ghost_marshall23(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::south_ghost_23_closure* genClosure = new Closure_Timestepping::south_ghost_23_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->south_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_south_ghost_marshall23(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::south_ghost_23_closure* genClosure = new Closure_Timestepping::south_ghost_23_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->south_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_south_ghost_marshall23(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::south_ghost_23_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_top_ghost_marshall24(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_top_ghost_marshall24() {
  int epidx = CkRegisterEp("top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_top_ghost_marshall24), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_top_ghost_marshall24);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_top_ghost_marshall24);

  return epidx;
}

void CkIndex_Timestepping::_call_top_ghost_marshall24(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::top_ghost_24_closure* genClosure = new Closure_Timestepping::top_ghost_24_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->top_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_top_ghost_marshall24(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::top_ghost_24_closure* genClosure = new Closure_Timestepping::top_ghost_24_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->top_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_top_ghost_marshall24(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::top_ghost_24_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_bottom_ghost_marshall25(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_bottom_ghost_marshall25() {
  int epidx = CkRegisterEp("bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_bottom_ghost_marshall25), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_bottom_ghost_marshall25);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_bottom_ghost_marshall25);

  return epidx;
}

void CkIndex_Timestepping::_call_bottom_ghost_marshall25(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::bottom_ghost_25_closure* genClosure = new Closure_Timestepping::bottom_ghost_25_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->bottom_ghost(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_bottom_ghost_marshall25(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::bottom_ghost_25_closure* genClosure = new Closure_Timestepping::bottom_ghost_25_closure();
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->bottom_ghost(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_bottom_ghost_marshall25(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::bottom_ghost_25_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void continue_timestepping();
 */
void CProxy_Timestepping::continue_timestepping(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_continue_timestepping_void(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_continue_timestepping_void() {
  int epidx = CkRegisterEp("continue_timestepping()",
      reinterpret_cast<CkCallFnPtr>(_call_continue_timestepping_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Timestepping::_call_continue_timestepping_void(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  impl_obj->continue_timestepping();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::continue_timestepping_26_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  impl_off_globalidx3_srcpts=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_globalidx3_srcpts=sizeof(int)*(num_srcpts));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_globalidx3_srcpts,globalidx3_srcpts,impl_cnt_globalidx3_srcpts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure();
  implP|genClosure->idx3_of_sendingchare;
  implP|genClosure->num_srcpts;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  implP|impl_off_globalidx3_srcpts;
  implP|impl_cnt_globalidx3_srcpts;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->globalidx3_srcpts = (int *)(impl_buf+impl_off_globalidx3_srcpts);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_idx3srcpt_tosend(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure();
  implP|genClosure->idx3_of_sendingchare;
  implP|genClosure->num_srcpts;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  implP|impl_off_globalidx3_srcpts;
  implP|impl_cnt_globalidx3_srcpts;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->globalidx3_srcpts = (int *)(impl_buf+impl_off_globalidx3_srcpts);
  impl_obj->receiv_nonlocalinnerbc_idx3srcpt_tosend(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> idx3_of_sendingchare;
  implP|idx3_of_sendingchare;
  PUP::detail::TemporaryObjectHolder<int> num_srcpts;
  implP|num_srcpts;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  implP|impl_off_globalidx3_srcpts;
  implP|impl_cnt_globalidx3_srcpts;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *globalidx3_srcpts=(int *)(impl_buf+impl_off_globalidx3_srcpts);
  if (implDestP.hasComments()) implDestP.comment("idx3_of_sendingchare");
  implDestP|idx3_of_sendingchare;
  if (implDestP.hasComments()) implDestP.comment("num_srcpts");
  implDestP|num_srcpts;
  if (implDestP.hasComments()) implDestP.comment("globalidx3_srcpts");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*globalidx3_srcpts))<impl_cnt_globalidx3_srcpts;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|globalidx3_srcpts[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_k_odd_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_k_odd_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_k_even_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_k_even_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_auxevol_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_auxevol_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_diagnostic_output_gfs(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_diagnostic_output_gfs(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxy_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35() {
  int epidx = CkRegisterEp("receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  genClosure->_impl_marshall = impl_msg_typed;
  CkReferenceMsg(genClosure->_impl_marshall);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure();
  implP|genClosure->src_chare_idx3;
  implP|genClosure->type_gfs;
  implP|genClosure->len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  genClosure->tmpBuffer = (REAL *)(impl_buf+impl_off_tmpBuffer);
  impl_obj->receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> src_chare_idx3;
  implP|src_chare_idx3;
  PUP::detail::TemporaryObjectHolder<int> type_gfs;
  implP|type_gfs;
  PUP::detail::TemporaryObjectHolder<int> len_tmpBuffer;
  implP|len_tmpBuffer;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  implP|impl_off_tmpBuffer;
  implP|impl_cnt_tmpBuffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  REAL *tmpBuffer=(REAL *)(impl_buf+impl_off_tmpBuffer);
  if (implDestP.hasComments()) implDestP.comment("src_chare_idx3");
  implDestP|src_chare_idx3;
  if (implDestP.hasComments()) implDestP.comment("type_gfs");
  implDestP|type_gfs;
  if (implDestP.hasComments()) implDestP.comment("len_tmpBuffer");
  implDestP|len_tmpBuffer;
  if (implDestP.hasComments()) implDestP.comment("tmpBuffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*tmpBuffer))<impl_cnt_tmpBuffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|tmpBuffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
 */
void CProxy_Timestepping::receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const REAL &wavespeed_at_outer_boundary
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Timestepping::idx_receiv_wavespeed_at_outer_boundary_marshall36(),0);
}

// Entry point registration function
int CkIndex_Timestepping::reg_receiv_wavespeed_at_outer_boundary_marshall36() {
  int epidx = CkRegisterEp("receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary)",
      reinterpret_cast<CkCallFnPtr>(_call_receiv_wavespeed_at_outer_boundary_marshall36), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_receiv_wavespeed_at_outer_boundary_marshall36);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_receiv_wavespeed_at_outer_boundary_marshall36);

  return epidx;
}

void CkIndex_Timestepping::_call_receiv_wavespeed_at_outer_boundary_marshall36(void* impl_msg, void* impl_obj_void)
{
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* genClosure = new Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure();
  implP|genClosure->wavespeed_at_outer_boundary;
  impl_buf+=CK_ALIGN(implP.size(),16);
  impl_obj->receiv_wavespeed_at_outer_boundary(genClosure);
  genClosure->deref();
}
int CkIndex_Timestepping::_callmarshall_receiv_wavespeed_at_outer_boundary_marshall36(char* impl_buf, void* impl_obj_void) {
  Timestepping* impl_obj = static_cast<Timestepping*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  PUP::fromMem implP(impl_buf);
  Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* genClosure = new Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure();
  implP|genClosure->wavespeed_at_outer_boundary;
  impl_buf+=CK_ALIGN(implP.size(),16);
  impl_obj->receiv_wavespeed_at_outer_boundary(genClosure);
  genClosure->deref();
  return implP.size();
}
void CkIndex_Timestepping::_marshallmessagepup_receiv_wavespeed_at_outer_boundary_marshall36(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const REAL &wavespeed_at_outer_boundary*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<REAL> wavespeed_at_outer_boundary;
  implP|wavespeed_at_outer_boundary;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("wavespeed_at_outer_boundary");
  implDestP|wavespeed_at_outer_boundary;
}
PUPable_def(SINGLE_ARG(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_Timestepping::reg_Timestepping_CkMigrateMessage() {
  int epidx = CkRegisterEp("Timestepping(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_Timestepping_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_Timestepping::_call_Timestepping_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<Timestepping> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(const CommondataObject &inData);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::ready_1d_y(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_1d_y_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::ready_1d_z(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_1d_z_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_2d_xy_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_ready_2d_yz_FileReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start();
 */
void CProxySection_Timestepping::start(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_y_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_1d_z_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_xy_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
 */
void CProxySection_Timestepping::start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_start_write_2d_yz_SessionReadyMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_y(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::test_written_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_y_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_1d_z(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::test_written_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_1d_z_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::test_written_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_xy_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void test_written_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::test_written_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_test_written_2d_yz_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_y(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::closed_1d_y(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_1d_y_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_1d_z(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::closed_1d_z(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_1d_z_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_xy(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::closed_2d_xy(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_2d_xy_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void closed_2d_yz(CkReductionMsg* impl_msg);
 */
void CProxySection_Timestepping::closed_2d_yz(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_closed_2d_yz_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
 */
void CProxySection_Timestepping::diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Ck::IO::Session &token, int which_diagnostics_part
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Ck::IO::Session>::type>::type &)token;
    implP|which_diagnostics_part;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_diagnostics_ckio_marshall19(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_east_ghost_marshall20(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_west_ghost_marshall21(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_north_ghost_marshall22(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_south_ghost_marshall23(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_top_ghost_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_bottom_ghost_marshall25(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void continue_timestepping();
 */
void CProxySection_Timestepping::continue_timestepping(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_continue_timestepping_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_globalidx3_srcpts, impl_cnt_globalidx3_srcpts;
  impl_off_globalidx3_srcpts=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_globalidx3_srcpts=sizeof(int)*(num_srcpts));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|idx3_of_sendingchare;
    implP|num_srcpts;
    implP|impl_off_globalidx3_srcpts;
    implP|impl_cnt_globalidx3_srcpts;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_globalidx3_srcpts,globalidx3_srcpts,impl_cnt_globalidx3_srcpts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
 */
void CProxySection_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_tmpBuffer, impl_cnt_tmpBuffer;
  impl_off_tmpBuffer=impl_off=CK_ALIGN(impl_off,sizeof(REAL));
  impl_off+=(impl_cnt_tmpBuffer=sizeof(REAL)*(len_tmpBuffer));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|src_chare_idx3;
    implP|type_gfs;
    implP|len_tmpBuffer;
    implP|impl_off_tmpBuffer;
    implP|impl_cnt_tmpBuffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_tmpBuffer,tmpBuffer,impl_cnt_tmpBuffer);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
 */
void CProxySection_Timestepping::receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const REAL &wavespeed_at_outer_boundary
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<REAL>::type>::type &)wavespeed_at_outer_boundary;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Timestepping::idx_receiv_wavespeed_at_outer_boundary_marshall36(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Timestepping(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Timestepping::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, 3);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
  // REG: Timestepping(const CommondataObject &inData);
  idx_Timestepping_marshall1();

  // REG: void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
  idx_ready_1d_y_FileReadyMsg();

  // REG: void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
  idx_ready_1d_z_FileReadyMsg();

  // REG: void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
  idx_ready_2d_xy_FileReadyMsg();

  // REG: void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
  idx_ready_2d_yz_FileReadyMsg();

  // REG: void start();
  idx_start_void();

  // REG: void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
  idx_start_write_1d_y_SessionReadyMsg();

  // REG: void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
  idx_start_write_1d_z_SessionReadyMsg();

  // REG: void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
  idx_start_write_2d_xy_SessionReadyMsg();

  // REG: void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
  idx_start_write_2d_yz_SessionReadyMsg();

  // REG: void test_written_1d_y(CkReductionMsg* impl_msg);
  idx_test_written_1d_y_CkReductionMsg();

  // REG: void test_written_1d_z(CkReductionMsg* impl_msg);
  idx_test_written_1d_z_CkReductionMsg();

  // REG: void test_written_2d_xy(CkReductionMsg* impl_msg);
  idx_test_written_2d_xy_CkReductionMsg();

  // REG: void test_written_2d_yz(CkReductionMsg* impl_msg);
  idx_test_written_2d_yz_CkReductionMsg();

  // REG: void closed_1d_y(CkReductionMsg* impl_msg);
  idx_closed_1d_y_CkReductionMsg();

  // REG: void closed_1d_z(CkReductionMsg* impl_msg);
  idx_closed_1d_z_CkReductionMsg();

  // REG: void closed_2d_xy(CkReductionMsg* impl_msg);
  idx_closed_2d_xy_CkReductionMsg();

  // REG: void closed_2d_yz(CkReductionMsg* impl_msg);
  idx_closed_2d_yz_CkReductionMsg();

  // REG: void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
  idx_diagnostics_ckio_marshall19();

  // REG: void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_east_ghost_marshall20();

  // REG: void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_west_ghost_marshall21();

  // REG: void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_north_ghost_marshall22();

  // REG: void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_south_ghost_marshall23();

  // REG: void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_top_ghost_marshall24();

  // REG: void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_bottom_ghost_marshall25();

  // REG: void continue_timestepping();
  idx_continue_timestepping_void();

  // REG: void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
  idx_receiv_nonlocalinnerbc_idx3srcpt_tosend_marshall27();

  // REG: void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_marshall28();

  // REG: void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_k_odd_gfs_marshall29();

  // REG: void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_k_even_gfs_marshall30();

  // REG: void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_y_n_gfs_marshall31();

  // REG: void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_auxevol_gfs_marshall32();

  // REG: void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_diagnostic_output_gfs_marshall33();

  // REG: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_marshall34();

  // REG: void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
  idx_receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_marshall35();

  // REG: void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
  idx_receiv_wavespeed_at_outer_boundary_marshall36();

  // REG: Timestepping(CkMigrateMessage* impl_msg);
  idx_Timestepping_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_Timestepping_CkMigrateMessage());

  Timestepping::__sdag_register(); // Potentially missing Timestepping_SDAG_CODE in your class definition?
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
void Timestepping::start(){
  CkPrintf("Error> Direct call to SDAG entry method \'%s::%s\'!\n", "Timestepping", "start()"); 
  CkAbort("Direct SDAG call is not allowed for SDAG entry methods having when constructs. Call such SDAG methods using a proxy"); 
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::_sdag_fnc_start() {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _slist_0();
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::start_end() {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_0() {
  _if_0();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_0_end() {
  start_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_0() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_1();
  } else {
    _if_0_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_0_end() {
  _if_1();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_1() {
  _serial_0();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_1_end() {
  _if_0_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_0() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_0()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 15 "timestepping.ci"
 send_nonlocalinnerbc_idx3srcpts_toreceiv(); 
#line 5685 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_1_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_1() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_2();
  } else {
    _if_1_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_1_end() {
  _serial_2();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_2() {
  _for_0();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_2_end() {
  _if_1_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_0() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares) {
    _slist_3();
  } else {
    _slist_2_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_0_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares) {
    _slist_3();
  } else {
    _slist_2_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_3() {
  _when_0();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_3_end() {
  _for_0_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_0() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(0, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_4(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(0);
    c->anyEntries.push_back(0);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_0_end(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* gen0) {
  _slist_3_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_4(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* gen0) {
  _serial_1(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_4_end(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* gen0) {
  _when_0_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_1(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_1()), CkMyPe(), 0, &projID, this); 
  {
    int& idx3_of_sendingchare = gen0->getP0();
    int& num_srcpts = gen0->getP1();
    int*& globalidx3_srcpts = gen0->getP2();
    { // begin serial block
#line 20 "timestepping.ci"
 process_nonlocalinnerbc_idx3srcpt_tosend(idx3_of_sendingchare, num_srcpts, globalidx3_srcpts); 
#line 5813 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_4_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_2() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_2()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 24 "timestepping.ci"

        initial_data(&commondata, griddata_chare);
        auxevol_gfs_set_to_constant(&commondata, &griddata_chare[grid].params, griddata_chare[grid].xx, &griddata_chare[grid].gridfuncs);
        send_wavespeed_at_outer_boundary(grid);
      
#line 5833 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _when_1();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_1() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(1, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_5(static_cast<Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(1);
    c->anyEntries.push_back(1);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_1_end(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* gen0) {
  _serial_4();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_5(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* gen0) {
  _serial_3(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_5_end(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* gen0) {
  _when_1_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_3(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_3()), CkMyPe(), 0, &projID, this); 
  {
    REAL& wavespeed_at_outer_boundary = gen0->getP0();
    { // begin serial block
#line 30 "timestepping.ci"
 griddata_chare[grid].params.wavespeed_at_outer_boundary = wavespeed_at_outer_boundary; 
#line 5889 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_5_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_4() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_4()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 32 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, EAST_WEST, grid); 
#line 5905 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_2();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_2() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_6();
  } else {
    _if_2_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_2_end() {
  _if_3();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_6() {
  _when_2();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_6_end() {
  _if_2_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_2() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_7(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(2);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_2_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_6_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_7(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_5(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_7_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_2_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_5(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_5()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 35 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 5995 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_7_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_3() {
  if (thisIndex.x > 0) {
    _slist_8();
  } else {
    _if_3_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_3_end() {
  _if_4();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_8() {
  _when_3();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_8_end() {
  _if_3_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_3() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_9(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(3);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_3_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_8_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_9(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_6(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_9_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_3_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_6(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_6()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 40 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6086 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_9_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_4() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_10();
  } else {
    _if_4_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_4_end() {
  _serial_10();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_10() {
  _serial_7();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_10_end() {
  _if_4_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_7() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_7()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 44 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, EAST_WEST, grid); 
#line 6134 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_5();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_5() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_11();
  } else {
    _if_5_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_5_end() {
  _if_6();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_11() {
  _when_4();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_11_end() {
  _if_5_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_4() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_12(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(4);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_4_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_11_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_12(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_8(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_12_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_4_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_8(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_8()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 47 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6224 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_12_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_6() {
  if (thisIndex.x > 0) {
    _slist_13();
  } else {
    _if_6_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_6_end() {
  _slist_10_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_13() {
  _when_5();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_13_end() {
  _if_6_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_5() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_14(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(5);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_5_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_13_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_14(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_9(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_14_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_5_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_9(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_9()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 52 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6315 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_14_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_10() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_10()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 56 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, NORTH_SOUTH, grid); 
#line 6331 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_7();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_7() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_15();
  } else {
    _if_7_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_7_end() {
  _if_8();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_15() {
  _when_6();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_15_end() {
  _if_7_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_6() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_16(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(6);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_6_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_15_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_16(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_11(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_16_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_6_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_11(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_11()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 59 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6421 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_16_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_8() {
  if (thisIndex.y > 0) {
    _slist_17();
  } else {
    _if_8_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_8_end() {
  _if_9();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_17() {
  _when_7();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_17_end() {
  _if_8_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_7() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_18(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(7);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_7_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_17_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_18(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_12(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_18_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_7_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_12(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_12()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 64 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6512 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_18_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_9() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_19();
  } else {
    _if_9_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_9_end() {
  _serial_16();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_19() {
  _serial_13();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_19_end() {
  _if_9_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_13() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_13()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 68 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, NORTH_SOUTH, grid); 
#line 6560 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_10();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_10() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_20();
  } else {
    _if_10_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_10_end() {
  _if_11();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_20() {
  _when_8();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_20_end() {
  _if_10_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_8() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_21(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(8);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_8_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_20_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_21(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_14(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_21_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_8_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_14(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_14()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 71 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6650 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_21_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_11() {
  if (thisIndex.y > 0) {
    _slist_22();
  } else {
    _if_11_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_11_end() {
  _slist_19_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_22() {
  _when_9();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_22_end() {
  _if_11_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_9() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_23(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(9);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_9_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_22_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_23(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_15(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_23_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_9_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_15(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_15()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 76 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6741 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_23_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_16() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_16()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 80 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, TOP_BOTTOM, grid); 
#line 6757 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_12();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_12() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_24();
  } else {
    _if_12_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_12_end() {
  _if_13();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_24() {
  _when_10();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_24_end() {
  _if_12_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_10() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_25(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(10);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_10_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_24_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_25(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_17(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_25_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_10_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_17(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_17()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 83 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6847 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_25_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_13() {
  if (thisIndex.z > 0) {
    _slist_26();
  } else {
    _if_13_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_13_end() {
  _if_14();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_26() {
  _when_11();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_26_end() {
  _if_13_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_11() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_27(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(11);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_11_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_26_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_27(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_18(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_27_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_11_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_18(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_18()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 88 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 6938 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_27_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_14() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_28();
  } else {
    _if_14_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_14_end() {
  _while_0();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_28() {
  _serial_19();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_28_end() {
  _if_14_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_19() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_19()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 92 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, TOP_BOTTOM, grid); 
#line 6986 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_15();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_15() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_29();
  } else {
    _if_15_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_15_end() {
  _if_16();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_29() {
  _when_12();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_29_end() {
  _if_15_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_12() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_30(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(12);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_12_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_29_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_30(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_20(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_30_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_12_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_20(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_20()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 95 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 7076 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_30_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_16() {
  if (thisIndex.z > 0) {
    _slist_31();
  } else {
    _if_16_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_16_end() {
  _slist_28_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_31() {
  _when_13();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_31_end() {
  _if_16_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_13() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_32(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(13);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_13_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_31_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_32(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_21(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_32_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_13_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_21(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_21()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 100 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 7167 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_32_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_while_0() {
  if (commondata.time < commondata.t_final) {
    _slist_33();
  } else {
    _serial_156();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_while_0_end() {
  if (commondata.time < commondata.t_final) {
    _slist_33();
  } else {
    _serial_156();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_33() {
  _serial_22();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_33_end() {
  _while_0_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_22() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_22()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 105 "timestepping.ci"
 time_start = commondata.time; 
#line 7219 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_23();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_23() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_23()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 106 "timestepping.ci"

          const int n_step = commondata.nn;
          const int outevery = commondata.diagnostics_output_every;
          write_diagnostics_this_step = n_step % outevery == 0;
        
#line 7238 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_17();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_17() {
  if (write_diagnostics_this_step) {
    _slist_34();
  } else {
    _if_17_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_17_end() {
  _if_18();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_34() {
  _serial_24();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_34_end() {
  _if_17_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_24() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_24()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 112 "timestepping.ci"

            for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
              const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0; const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1; const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
              const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
              BHAH_MALLOC(diagnostic_gfs[grid], TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));
            }
            diagnostic_gfs_set(&commondata, griddata_chare, diagnostic_gfs);
          
#line 7292 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_34_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_18() {
  if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
    _slist_35();
  } else {
    _if_18_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_18_end() {
  _when_30();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_35() {
  _serial_25();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_35_end() {
  _if_18_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_25() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_25()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 122 "timestepping.ci"

            progress_indicator(&commondata, griddata_chare);
            if (commondata.time + commondata.dt > commondata.t_final)
              printf("\n");
          
#line 7343 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_19();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_19() {
  if (write_diagnostics_this_step) {
    _slist_37();
  } else {
    _else_0();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_19_end() {
  _slist_35_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_else_0() {
  _slist_36();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_else_0_end() {
  _slist_35_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_36() {
  _serial_26();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_36_end() {
  _else_0_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_26() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_26()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 282 "timestepping.ci"
 thisProxy.continue_timestepping(); 
#line 7404 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_36_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_37() {
  _serial_27();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_37_end() {
  _if_19_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_27() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_27()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 128 "timestepping.ci"

              count_filewritten = 0;
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out1d-y", griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_y, &commondata, 1);
                Ck::IO::Options opts;
                CkCallback opened_1d_y(CkIndex_Timestepping::ready_1d_y(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_y, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out1d-z", griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_1d_z, &commondata, 1);
                Ck::IO::Options opts;
                CkCallback opened_1d_z(CkIndex_Timestepping::ready_1d_z(NULL), thisProxy);
                Ck::IO::open(filename, opened_1d_z, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out2d-xy", griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_xy, &commondata, 1);
                Ck::IO::Options opts;
                CkCallback opened_2d_xy(CkIndex_Timestepping::ready_2d_xy(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_xy, opts);
              }
              {
                char filename[512];
                build_outfile_name(filename, sizeof filename, "out2d-yz", griddata_chare[which_grid_diagnostics].diagnosticstruct.filename_2d_yz, &commondata, 1);
                Ck::IO::Options opts;
                CkCallback opened_2d_yz(CkIndex_Timestepping::ready_2d_yz(NULL), thisProxy);
                Ck::IO::open(filename, opened_2d_yz, opts);
              }
            
#line 7463 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _when_14();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_14() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(8, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_38(static_cast<Ck::IO::FileReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(14);
    c->anyEntries.push_back(8);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_14_end(Ck::IO::FileReadyMsg* gen0) {
  {
    Ck::IO::FileReadyMsg*& m_1d_y = gen0;
    CmiFree(UsrToEnv(m_1d_y));
  }
  _when_15();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_38(Ck::IO::FileReadyMsg* gen0) {
  _serial_28(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_38_end(Ck::IO::FileReadyMsg* gen0) {
  _when_14_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_28(Ck::IO::FileReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_28()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::FileReadyMsg*& m_1d_y = gen0;
    { // begin serial block
#line 160 "timestepping.ci"

                f_1d_y = m_1d_y->file;
                CkCallback sessionStart_1d_y(CkIndex_Timestepping::start_write_1d_y(0), thisProxy);
                CkCallback sessionEnd_1d_y(CkIndex_Timestepping::test_written_1d_y(0), thisProxy);
                int num_fields = griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1;
                int tot_num_diagnostic_pts = griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_y_pts;
                int totsizeinbytes = 23 * num_fields * tot_num_diagnostic_pts;
                Ck::IO::startSession(f_1d_y, totsizeinbytes, 0, sessionStart_1d_y, sessionEnd_1d_y);
                delete m_1d_y;
              
#line 7532 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_38_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_15() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(9, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_39(static_cast<Ck::IO::SessionReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(15);
    c->anyEntries.push_back(9);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_15_end(Ck::IO::SessionReadyMsg* gen0) {
  {
    Ck::IO::SessionReadyMsg*& m_1d_y = gen0;
    CmiFree(UsrToEnv(m_1d_y));
  }
  _when_16();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_39(Ck::IO::SessionReadyMsg* gen0) {
  _serial_29(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_39_end(Ck::IO::SessionReadyMsg* gen0) {
  _when_15_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_29(Ck::IO::SessionReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_29()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::SessionReadyMsg*& m_1d_y = gen0;
    { // begin serial block
#line 172 "timestepping.ci"

                thisProxy.diagnostics_ckio(m_1d_y->session, DIAGNOSTICS_WRITE_Y);
                delete m_1d_y;
              
#line 7596 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_39_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_16() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(10, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_40(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(16);
    c->anyEntries.push_back(10);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_16_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_1d_y = gen0;
    CmiFree(UsrToEnv(m_1d_y));
  }
  _when_17();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_40(CkReductionMsg* gen0) {
  _serial_30(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_40_end(CkReductionMsg* gen0) {
  _when_16_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_30(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_30()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_1d_y = gen0;
    { // begin serial block
#line 178 "timestepping.ci"

                printf("when start_write_1d_y(Ck::IO::SessionReadyMsg * m_1d_y) {\n");
                delete m_1d_y;
                CkCallback cb_1d_y(CkIndex_Timestepping::closed_1d_y(0), thisProxy);
                Ck::IO::close(f_1d_y, cb_1d_y);
                count_filewritten++;
              
#line 7663 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_40_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_17() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(11, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_41(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(17);
    c->anyEntries.push_back(11);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_17_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_1d_y = gen0;
    CmiFree(UsrToEnv(m_1d_y));
  }
  _when_18();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_41(CkReductionMsg* gen0) {
  _serial_31(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_41_end(CkReductionMsg* gen0) {
  _when_17_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_31(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_31()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_1d_y = gen0;
    { // begin serial block
#line 187 "timestepping.ci"
 delete m_1d_y; 
#line 7724 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_41_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_18() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(12, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_42(static_cast<Ck::IO::FileReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(18);
    c->anyEntries.push_back(12);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_18_end(Ck::IO::FileReadyMsg* gen0) {
  {
    Ck::IO::FileReadyMsg*& m_1d_z = gen0;
    CmiFree(UsrToEnv(m_1d_z));
  }
  _when_19();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_42(Ck::IO::FileReadyMsg* gen0) {
  _serial_32(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_42_end(Ck::IO::FileReadyMsg* gen0) {
  _when_18_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_32(Ck::IO::FileReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_32()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::FileReadyMsg*& m_1d_z = gen0;
    { // begin serial block
#line 190 "timestepping.ci"

                printf("when ready_1d_z(Ck::IO::FileReadyMsg * m_1d_z) {\n");
                f_1d_z = m_1d_z->file;
                CkCallback sessionStart_1d_z(CkIndex_Timestepping::start_write_1d_z(0), thisProxy);
                CkCallback sessionEnd_1d_z(CkIndex_Timestepping::test_written_1d_z(0), thisProxy);
                int num_fields = griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 1;
                int tot_num_diagnostic_pts = griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_1d_z_pts;
                int totsizeinbytes = 23 * num_fields * tot_num_diagnostic_pts;
                Ck::IO::startSession(f_1d_z, totsizeinbytes, 0, sessionStart_1d_z, sessionEnd_1d_z);
                delete m_1d_z;
              
#line 7795 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_42_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_19() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(13, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_43(static_cast<Ck::IO::SessionReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(19);
    c->anyEntries.push_back(13);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_19_end(Ck::IO::SessionReadyMsg* gen0) {
  {
    Ck::IO::SessionReadyMsg*& m_1d_z = gen0;
    CmiFree(UsrToEnv(m_1d_z));
  }
  _when_20();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_43(Ck::IO::SessionReadyMsg* gen0) {
  _serial_33(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_43_end(Ck::IO::SessionReadyMsg* gen0) {
  _when_19_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_33(Ck::IO::SessionReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_33()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::SessionReadyMsg*& m_1d_z = gen0;
    { // begin serial block
#line 203 "timestepping.ci"

                thisProxy.diagnostics_ckio(m_1d_z->session, DIAGNOSTICS_WRITE_Z);
                delete m_1d_z;
              
#line 7859 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_43_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_20() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(14, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_44(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(20);
    c->anyEntries.push_back(14);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_20_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_1d_z = gen0;
    CmiFree(UsrToEnv(m_1d_z));
  }
  _when_21();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_44(CkReductionMsg* gen0) {
  _serial_34(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_44_end(CkReductionMsg* gen0) {
  _when_20_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_34(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_34()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_1d_z = gen0;
    { // begin serial block
#line 209 "timestepping.ci"

                printf("when test_written_1d_z(CkReductionMsg * m_1d_z)  {\n");
                delete m_1d_z;
                CkCallback cb_1d_z(CkIndex_Timestepping::closed_1d_z(0), thisProxy);
                Ck::IO::close(f_1d_z, cb_1d_z);
                count_filewritten++;
              
#line 7926 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_44_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_21() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(15, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_45(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(21);
    c->anyEntries.push_back(15);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_21_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_1d_z = gen0;
    CmiFree(UsrToEnv(m_1d_z));
  }
  _when_22();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_45(CkReductionMsg* gen0) {
  _serial_35(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_45_end(CkReductionMsg* gen0) {
  _when_21_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_35(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_35()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_1d_z = gen0;
    { // begin serial block
#line 218 "timestepping.ci"
 delete m_1d_z; 
#line 7987 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_45_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_22() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(16, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_46(static_cast<Ck::IO::FileReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(22);
    c->anyEntries.push_back(16);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_22_end(Ck::IO::FileReadyMsg* gen0) {
  {
    Ck::IO::FileReadyMsg*& m_2d_xy = gen0;
    CmiFree(UsrToEnv(m_2d_xy));
  }
  _when_23();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_46(Ck::IO::FileReadyMsg* gen0) {
  _serial_36(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_46_end(Ck::IO::FileReadyMsg* gen0) {
  _when_22_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_36(Ck::IO::FileReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_36()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::FileReadyMsg*& m_2d_xy = gen0;
    { // begin serial block
#line 221 "timestepping.ci"

                f_2d_xy = m_2d_xy->file;
                CkCallback sessionStart_2d_xy(CkIndex_Timestepping::start_write_2d_xy(0), thisProxy);
                CkCallback sessionEnd_2d_xy(CkIndex_Timestepping::test_written_2d_xy(0), thisProxy);
                int num_fields = griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2;
                int tot_num_diagnostic_pts = griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_xy_pts;
                int totsizeinbytes = (23 + (24 * (num_fields - 1))) * tot_num_diagnostic_pts;
                Ck::IO::startSession(f_2d_xy, totsizeinbytes, 0, sessionStart_2d_xy, sessionEnd_2d_xy);
                delete m_2d_xy;
              
#line 8057 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_46_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_23() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(17, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_47(static_cast<Ck::IO::SessionReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(23);
    c->anyEntries.push_back(17);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_23_end(Ck::IO::SessionReadyMsg* gen0) {
  {
    Ck::IO::SessionReadyMsg*& m_2d_xy = gen0;
    CmiFree(UsrToEnv(m_2d_xy));
  }
  _when_24();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_47(Ck::IO::SessionReadyMsg* gen0) {
  _serial_37(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_47_end(Ck::IO::SessionReadyMsg* gen0) {
  _when_23_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_37(Ck::IO::SessionReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_37()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::SessionReadyMsg*& m_2d_xy = gen0;
    { // begin serial block
#line 233 "timestepping.ci"

                thisProxy.diagnostics_ckio(m_2d_xy->session, DIAGNOSTICS_WRITE_XY);
                delete m_2d_xy;
              
#line 8121 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_47_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_24() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(18, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_48(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(24);
    c->anyEntries.push_back(18);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_24_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_2d_xy = gen0;
    CmiFree(UsrToEnv(m_2d_xy));
  }
  _when_25();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_48(CkReductionMsg* gen0) {
  _serial_38(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_48_end(CkReductionMsg* gen0) {
  _when_24_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_38(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_38()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_2d_xy = gen0;
    { // begin serial block
#line 239 "timestepping.ci"

                delete m_2d_xy;
                CkCallback cb_2d_xy(CkIndex_Timestepping::closed_2d_xy(0), thisProxy);
                Ck::IO::close(f_2d_xy, cb_2d_xy);
                count_filewritten++;
              
#line 8187 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_48_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_25() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(19, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_49(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(25);
    c->anyEntries.push_back(19);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_25_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_2d_xy = gen0;
    CmiFree(UsrToEnv(m_2d_xy));
  }
  _when_26();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_49(CkReductionMsg* gen0) {
  _serial_39(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_49_end(CkReductionMsg* gen0) {
  _when_25_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_39(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_39()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_2d_xy = gen0;
    { // begin serial block
#line 247 "timestepping.ci"
 delete m_2d_xy; 
#line 8248 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_49_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_26() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(20, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_50(static_cast<Ck::IO::FileReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(26);
    c->anyEntries.push_back(20);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_26_end(Ck::IO::FileReadyMsg* gen0) {
  {
    Ck::IO::FileReadyMsg*& m_2d_yz = gen0;
    CmiFree(UsrToEnv(m_2d_yz));
  }
  _when_27();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_50(Ck::IO::FileReadyMsg* gen0) {
  _serial_40(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_50_end(Ck::IO::FileReadyMsg* gen0) {
  _when_26_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_40(Ck::IO::FileReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_40()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::FileReadyMsg*& m_2d_yz = gen0;
    { // begin serial block
#line 250 "timestepping.ci"

                f_2d_yz = m_2d_yz->file;
                CkCallback sessionStart_2d_yz(CkIndex_Timestepping::start_write_2d_yz(0), thisProxy);
                CkCallback sessionEnd_2d_yz(CkIndex_Timestepping::test_written_2d_yz(0), thisProxy);
                int num_fields = griddata_chare[which_grid_diagnostics].diagnosticstruct.num_output_quantities + 2;
                int tot_num_diagnostic_pts = griddata_chare[which_grid_diagnostics].diagnosticstruct.tot_num_diagnostic_2d_yz_pts;
                int totsizeinbytes = (23 + (24 * (num_fields - 1))) * tot_num_diagnostic_pts;
                Ck::IO::startSession(f_2d_yz, totsizeinbytes, 0, sessionStart_2d_yz, sessionEnd_2d_yz);
                delete m_2d_yz;
              
#line 8318 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_50_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_27() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(21, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_51(static_cast<Ck::IO::SessionReadyMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(27);
    c->anyEntries.push_back(21);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_27_end(Ck::IO::SessionReadyMsg* gen0) {
  {
    Ck::IO::SessionReadyMsg*& m_2d_yz = gen0;
    CmiFree(UsrToEnv(m_2d_yz));
  }
  _when_28();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_51(Ck::IO::SessionReadyMsg* gen0) {
  _serial_41(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_51_end(Ck::IO::SessionReadyMsg* gen0) {
  _when_27_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_41(Ck::IO::SessionReadyMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_41()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::SessionReadyMsg*& m_2d_yz = gen0;
    { // begin serial block
#line 262 "timestepping.ci"

                thisProxy.diagnostics_ckio(m_2d_yz->session, DIAGNOSTICS_WRITE_YZ);
                delete m_2d_yz;
              
#line 8382 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_51_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_28() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(22, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_52(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(28);
    c->anyEntries.push_back(22);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_28_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_2d_yz = gen0;
    CmiFree(UsrToEnv(m_2d_yz));
  }
  _when_29();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_52(CkReductionMsg* gen0) {
  _serial_42(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_52_end(CkReductionMsg* gen0) {
  _when_28_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_42(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_42()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_2d_yz = gen0;
    { // begin serial block
#line 268 "timestepping.ci"

                delete m_2d_yz;
                CkCallback cb_2d_yz(CkIndex_Timestepping::closed_2d_yz(0), thisProxy);
                Ck::IO::close(f_2d_yz, cb_2d_yz);
                count_filewritten++;
              
#line 8448 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_52_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_29() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(23, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_53(static_cast<CkReductionMsg*>(static_cast<SDAG::MsgClosure*>(buf0->cl)->msg));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(29);
    c->anyEntries.push_back(23);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_29_end(CkReductionMsg* gen0) {
  {
    CkReductionMsg*& m_2d_yz = gen0;
    CmiFree(UsrToEnv(m_2d_yz));
  }
  _if_20();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_53(CkReductionMsg* gen0) {
  _serial_43(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_53_end(CkReductionMsg* gen0) {
  _when_29_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_43(CkReductionMsg* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_43()), CkMyPe(), 0, &projID, this); 
  {
    CkReductionMsg*& m_2d_yz = gen0;
    { // begin serial block
#line 276 "timestepping.ci"
 delete m_2d_yz; 
#line 8509 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_53_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_20() {
  if (count_filewritten == expected_count_filewritten) {
    _slist_54();
  } else {
    _if_20_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_20_end() {
  _slist_37_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_54() {
  _serial_44();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_54_end() {
  _if_20_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_44() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_44()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 279 "timestepping.ci"
 thisProxy.continue_timestepping(); 
#line 8557 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_54_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_30() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(24, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_55();
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(30);
    c->anyEntries.push_back(24);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_30_end() {
  _serial_46();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_55() {
  _if_21();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_55_end() {
  _when_30_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_21() {
  if (write_diagnostics_this_step) {
    _slist_56();
  } else {
    _if_21_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_21_end() {
  _slist_55_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_56() {
  _serial_45();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_56_end() {
  _if_21_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_45() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_45()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 287 "timestepping.ci"

              for (int grid = 0; grid < commondata.NUMGRIDS; grid++)
                free(diagnostic_gfs[grid]);
            
#line 8646 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_56_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_46() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_46()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 293 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K1, MOL_PRE_RK_UPDATE); 
#line 8661 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_22();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_22() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_57();
  } else {
    _if_22_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_22_end() {
  _if_23();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_57() {
  _serial_47();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_57_end() {
  _if_22_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_47() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_47()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 295 "timestepping.ci"
 send_nonlocalinnerbc_data(K_ODD_GFS, grid); 
#line 8708 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_57_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_23() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_58();
  } else {
    _if_23_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_23_end() {
  _serial_50();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_58() {
  _for_1();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_58_end() {
  _if_23_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_1() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_59();
  } else {
    _serial_49();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_1_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_59();
  } else {
    _serial_49();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_59() {
  _when_31();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_59_end() {
  _for_1_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_31() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(25, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_60(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(31);
    c->anyEntries.push_back(25);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_31_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _slist_59_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_60(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _serial_48(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_60_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _when_31_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_48(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_48()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 300 "timestepping.ci"

                set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                type_gfs_nonlocal_innerbc = type_gfs;
              
#line 8840 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_60_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_49() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_49()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 306 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 8856 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_58_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_50() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_50()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 308 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K1, MOL_RK_UPDATE); 
#line 8871 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_51();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_51() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_51()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 309 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K1, MOL_POST_RK_UPDATE); 
#line 8886 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_24();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_24() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_61();
  } else {
    _if_24_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_24_end() {
  _serial_55();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_61() {
  _if_25();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_61_end() {
  _if_24_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_25() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_62();
  } else {
    _if_25_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_25_end() {
  _if_26();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_62() {
  _serial_52();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_62_end() {
  _if_25_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_52() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_52()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 312 "timestepping.ci"
 send_nonlocalinnerbc_data(AUXEVOL_GFS, grid); 
#line 8965 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_62_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_26() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_63();
  } else {
    _if_26_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_26_end() {
  _slist_61_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_63() {
  _for_2();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_63_end() {
  _if_26_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_2() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_64();
  } else {
    _serial_54();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_2_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_64();
  } else {
    _serial_54();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_64() {
  _when_32();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_64_end() {
  _for_2_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_32() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(26, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_65(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(32);
    c->anyEntries.push_back(26);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_32_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _slist_64_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_65(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _serial_53(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_65_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _when_32_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_53(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_53()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 317 "timestepping.ci"

                  set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                  type_gfs_nonlocal_innerbc = type_gfs;
                
#line 9097 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_65_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_54() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_54()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 323 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 9113 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_63_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_55() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_55()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 326 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, EAST_WEST, grid); 
#line 9128 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_27();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_27() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_66();
  } else {
    _if_27_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_27_end() {
  _if_28();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_66() {
  _when_33();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_66_end() {
  _if_27_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_33() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_67(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(33);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_33_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_66_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_67(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_56(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_67_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_33_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_56(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_56()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 329 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9218 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_67_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_28() {
  if (thisIndex.x > 0) {
    _slist_68();
  } else {
    _if_28_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_28_end() {
  _if_29();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_68() {
  _when_34();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_68_end() {
  _if_28_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_34() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_69(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(34);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_34_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_68_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_69(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_57(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_69_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_34_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_57(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_57()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 334 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9309 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_69_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_29() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_70();
  } else {
    _if_29_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_29_end() {
  _serial_61();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_70() {
  _serial_58();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_70_end() {
  _if_29_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_58() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_58()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 338 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, EAST_WEST, grid); 
#line 9357 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_30();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_30() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_71();
  } else {
    _if_30_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_30_end() {
  _if_31();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_71() {
  _when_35();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_71_end() {
  _if_30_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_35() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_72(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(35);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_35_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_71_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_72(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_59(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_72_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_35_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_59(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_59()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 341 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9447 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_72_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_31() {
  if (thisIndex.x > 0) {
    _slist_73();
  } else {
    _if_31_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_31_end() {
  _slist_70_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_73() {
  _when_36();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_73_end() {
  _if_31_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_36() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_74(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(36);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_36_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_73_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_74(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_60(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_74_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_36_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_60(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_60()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 346 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9538 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_74_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_61() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_61()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 350 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, NORTH_SOUTH, grid); 
#line 9554 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_32();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_32() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_75();
  } else {
    _if_32_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_32_end() {
  _if_33();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_75() {
  _when_37();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_75_end() {
  _if_32_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_37() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_76(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(37);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_37_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_75_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_76(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_62(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_76_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_37_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_62(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_62()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 353 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9644 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_76_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_33() {
  if (thisIndex.y > 0) {
    _slist_77();
  } else {
    _if_33_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_33_end() {
  _if_34();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_77() {
  _when_38();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_77_end() {
  _if_33_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_38() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_78(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(38);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_38_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_77_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_78(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_63(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_78_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_38_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_63(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_63()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 358 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9735 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_78_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_34() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_79();
  } else {
    _if_34_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_34_end() {
  _serial_67();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_79() {
  _serial_64();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_79_end() {
  _if_34_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_64() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_64()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 362 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, NORTH_SOUTH, grid); 
#line 9783 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_35();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_35() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_80();
  } else {
    _if_35_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_35_end() {
  _if_36();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_80() {
  _when_39();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_80_end() {
  _if_35_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_39() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_81(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(39);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_39_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_80_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_81(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_65(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_81_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_39_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_65(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_65()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 365 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9873 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_81_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_36() {
  if (thisIndex.y > 0) {
    _slist_82();
  } else {
    _if_36_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_36_end() {
  _slist_79_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_82() {
  _when_40();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_82_end() {
  _if_36_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_40() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_83(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(40);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_40_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_82_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_83(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_66(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_83_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_40_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_66(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_66()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 370 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 9964 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_83_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_67() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_67()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 374 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, TOP_BOTTOM, grid); 
#line 9980 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_37();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_37() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_84();
  } else {
    _if_37_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_37_end() {
  _if_38();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_84() {
  _when_41();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_84_end() {
  _if_37_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_41() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_85(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(41);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_41_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_84_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_85(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_68(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_85_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_41_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_68(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_68()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 377 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 10070 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_85_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_38() {
  if (thisIndex.z > 0) {
    _slist_86();
  } else {
    _if_38_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_38_end() {
  _if_39();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_86() {
  _when_42();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_86_end() {
  _if_38_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_42() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_87(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(42);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_42_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_86_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_87(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_69(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_87_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_42_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_69(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_69()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 382 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 10161 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_87_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_39() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_88();
  } else {
    _if_39_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_39_end() {
  _serial_73();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_88() {
  _serial_70();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_88_end() {
  _if_39_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_70() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_70()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 386 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, TOP_BOTTOM, grid); 
#line 10209 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_40();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_40() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_89();
  } else {
    _if_40_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_40_end() {
  _if_41();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_89() {
  _when_43();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_89_end() {
  _if_40_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_43() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_90(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(43);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_43_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_89_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_90(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_71(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_90_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_43_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_71(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_71()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 389 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 10299 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_90_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_41() {
  if (thisIndex.z > 0) {
    _slist_91();
  } else {
    _if_41_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_41_end() {
  _slist_88_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_91() {
  _when_44();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_91_end() {
  _if_41_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_44() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_92(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(44);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_44_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_91_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_92(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_72(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_92_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_44_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_72(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_72()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 394 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 10390 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_92_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_73() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_73()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 398 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K2, MOL_PRE_RK_UPDATE); 
#line 10406 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_42();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_42() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_93();
  } else {
    _if_42_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_42_end() {
  _if_43();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_93() {
  _serial_74();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_93_end() {
  _if_42_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_74() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_74()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 400 "timestepping.ci"
 send_nonlocalinnerbc_data(K_EVEN_GFS, grid); 
#line 10453 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_93_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_43() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_94();
  } else {
    _if_43_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_43_end() {
  _serial_77();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_94() {
  _for_3();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_94_end() {
  _if_43_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_3() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_95();
  } else {
    _serial_76();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_3_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_95();
  } else {
    _serial_76();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_95() {
  _when_45();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_95_end() {
  _for_3_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_45() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(27, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_96(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(45);
    c->anyEntries.push_back(27);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_45_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _slist_95_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_96(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _serial_75(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_96_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _when_45_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_75(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_75()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 405 "timestepping.ci"

                set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                type_gfs_nonlocal_innerbc = type_gfs;
              
#line 10585 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_96_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_76() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_76()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 411 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 10601 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_94_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_77() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_77()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 413 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K2, MOL_RK_UPDATE); 
#line 10616 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_78();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_78() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_78()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 414 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K2, MOL_POST_RK_UPDATE); 
#line 10631 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_44();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_44() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_97();
  } else {
    _if_44_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_44_end() {
  _serial_82();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_97() {
  _if_45();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_97_end() {
  _if_44_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_45() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_98();
  } else {
    _if_45_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_45_end() {
  _if_46();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_98() {
  _serial_79();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_98_end() {
  _if_45_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_79() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_79()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 417 "timestepping.ci"
 send_nonlocalinnerbc_data(AUXEVOL_GFS, grid); 
#line 10710 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_98_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_46() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_99();
  } else {
    _if_46_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_46_end() {
  _slist_97_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_99() {
  _for_4();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_99_end() {
  _if_46_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_4() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_100();
  } else {
    _serial_81();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_4_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_100();
  } else {
    _serial_81();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_100() {
  _when_46();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_100_end() {
  _for_4_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_46() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(26, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_101(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(46);
    c->anyEntries.push_back(26);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_46_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _slist_100_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_101(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _serial_80(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_101_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _when_46_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_80(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_80()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 422 "timestepping.ci"

                  set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                  type_gfs_nonlocal_innerbc = type_gfs;
                
#line 10842 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_101_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_81() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_81()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 428 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 10858 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_99_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_82() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_82()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 431 "timestepping.ci"
 send_neighbor_data(K_EVEN_GFS, EAST_WEST, grid); 
#line 10873 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_47();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_47() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_102();
  } else {
    _if_47_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_47_end() {
  _if_48();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_102() {
  _when_47();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_102_end() {
  _if_47_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_47() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_103(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(47);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_47_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_102_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_103(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_83(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_103_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_47_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_83(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_83()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 434 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 10963 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_103_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_48() {
  if (thisIndex.x > 0) {
    _slist_104();
  } else {
    _if_48_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_48_end() {
  _if_49();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_104() {
  _when_48();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_104_end() {
  _if_48_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_48() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_105(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(48);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_48_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_104_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_105(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_84(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_105_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_48_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_84(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_84()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 439 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11054 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_105_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_49() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_106();
  } else {
    _if_49_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_49_end() {
  _serial_88();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_106() {
  _serial_85();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_106_end() {
  _if_49_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_85() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_85()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 443 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, EAST_WEST, grid); 
#line 11102 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_50();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_50() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_107();
  } else {
    _if_50_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_50_end() {
  _if_51();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_107() {
  _when_49();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_107_end() {
  _if_50_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_49() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_108(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(49);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_49_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_107_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_108(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_86(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_108_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_49_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_86(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_86()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 446 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11192 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_108_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_51() {
  if (thisIndex.x > 0) {
    _slist_109();
  } else {
    _if_51_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_51_end() {
  _slist_106_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_109() {
  _when_50();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_109_end() {
  _if_51_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_50() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_110(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(50);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_50_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_109_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_110(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_87(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_110_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_50_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_87(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_87()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 451 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11283 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_110_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_88() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_88()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 455 "timestepping.ci"
 send_neighbor_data(K_EVEN_GFS, NORTH_SOUTH, grid); 
#line 11299 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_52();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_52() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_111();
  } else {
    _if_52_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_52_end() {
  _if_53();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_111() {
  _when_51();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_111_end() {
  _if_52_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_51() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_112(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(51);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_51_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_111_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_112(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_89(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_112_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_51_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_89(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_89()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 458 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11389 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_112_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_53() {
  if (thisIndex.y > 0) {
    _slist_113();
  } else {
    _if_53_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_53_end() {
  _if_54();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_113() {
  _when_52();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_113_end() {
  _if_53_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_52() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_114(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(52);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_52_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_113_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_114(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_90(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_114_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_52_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_90(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_90()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 463 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11480 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_114_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_54() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_115();
  } else {
    _if_54_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_54_end() {
  _serial_94();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_115() {
  _serial_91();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_115_end() {
  _if_54_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_91() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_91()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 467 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, NORTH_SOUTH, grid); 
#line 11528 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_55();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_55() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_116();
  } else {
    _if_55_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_55_end() {
  _if_56();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_116() {
  _when_53();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_116_end() {
  _if_55_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_53() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_117(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(53);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_53_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_116_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_117(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_92(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_117_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_53_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_92(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_92()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 470 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11618 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_117_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_56() {
  if (thisIndex.y > 0) {
    _slist_118();
  } else {
    _if_56_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_56_end() {
  _slist_115_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_118() {
  _when_54();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_118_end() {
  _if_56_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_54() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_119(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(54);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_54_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_118_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_119(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_93(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_119_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_54_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_93(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_93()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 475 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11709 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_119_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_94() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_94()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 479 "timestepping.ci"
 send_neighbor_data(K_EVEN_GFS, TOP_BOTTOM, grid); 
#line 11725 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_57();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_57() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_120();
  } else {
    _if_57_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_57_end() {
  _if_58();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_120() {
  _when_55();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_120_end() {
  _if_57_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_55() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_121(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(55);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_55_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_120_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_121(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_95(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_121_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_55_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_95(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_95()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 482 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11815 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_121_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_58() {
  if (thisIndex.z > 0) {
    _slist_122();
  } else {
    _if_58_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_58_end() {
  _if_59();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_122() {
  _when_56();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_122_end() {
  _if_58_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_56() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_123(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(56);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_56_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_122_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_123(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_96(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_123_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_56_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_96(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_96()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 487 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 11906 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_123_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_59() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_124();
  } else {
    _if_59_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_59_end() {
  _serial_100();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_124() {
  _serial_97();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_124_end() {
  _if_59_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_97() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_97()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 491 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, TOP_BOTTOM, grid); 
#line 11954 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_60();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_60() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_125();
  } else {
    _if_60_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_60_end() {
  _if_61();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_125() {
  _when_57();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_125_end() {
  _if_60_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_57() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_126(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(57);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_57_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_125_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_126(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_98(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_126_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_57_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_98(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_98()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 494 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 12044 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_126_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_61() {
  if (thisIndex.z > 0) {
    _slist_127();
  } else {
    _if_61_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_61_end() {
  _slist_124_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_127() {
  _when_58();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_127_end() {
  _if_61_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_58() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_128(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(58);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_58_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_127_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_128(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_99(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_128_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_58_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_99(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_99()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 499 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 12135 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_128_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_100() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_100()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 503 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K3, MOL_PRE_RK_UPDATE); 
#line 12151 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_62();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_62() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_129();
  } else {
    _if_62_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_62_end() {
  _if_63();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_129() {
  _serial_101();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_129_end() {
  _if_62_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_101() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_101()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 505 "timestepping.ci"
 send_nonlocalinnerbc_data(K_ODD_GFS, grid); 
#line 12198 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_129_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_63() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_130();
  } else {
    _if_63_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_63_end() {
  _serial_104();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_130() {
  _for_5();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_130_end() {
  _if_63_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_5() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_131();
  } else {
    _serial_103();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_5_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_131();
  } else {
    _serial_103();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_131() {
  _when_59();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_131_end() {
  _for_5_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_59() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(25, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_132(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(59);
    c->anyEntries.push_back(25);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_59_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _slist_131_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_132(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _serial_102(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_132_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  _when_59_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_102(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_102()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 510 "timestepping.ci"

                set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                type_gfs_nonlocal_innerbc = type_gfs;
              
#line 12330 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_132_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_103() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_103()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 516 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 12346 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_130_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_104() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_104()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 518 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K3, MOL_RK_UPDATE); 
#line 12361 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_105();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_105() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_105()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 519 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K3, MOL_POST_RK_UPDATE); 
#line 12376 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_64();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_64() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_133();
  } else {
    _if_64_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_64_end() {
  _serial_109();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_133() {
  _if_65();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_133_end() {
  _if_64_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_65() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_134();
  } else {
    _if_65_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_65_end() {
  _if_66();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_134() {
  _serial_106();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_134_end() {
  _if_65_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_106() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_106()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 522 "timestepping.ci"
 send_nonlocalinnerbc_data(AUXEVOL_GFS, grid); 
#line 12455 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_134_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_66() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_135();
  } else {
    _if_66_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_66_end() {
  _slist_133_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_135() {
  _for_6();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_135_end() {
  _if_66_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_6() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_136();
  } else {
    _serial_108();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_6_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_136();
  } else {
    _serial_108();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_136() {
  _when_60();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_136_end() {
  _for_6_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_60() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(26, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_137(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(60);
    c->anyEntries.push_back(26);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_60_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _slist_136_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_137(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _serial_107(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_137_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _when_60_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_107(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_107()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 527 "timestepping.ci"

                  set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                  type_gfs_nonlocal_innerbc = type_gfs;
                
#line 12587 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_137_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_108() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_108()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 533 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 12603 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_135_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_109() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_109()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 536 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, EAST_WEST, grid); 
#line 12618 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_67();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_67() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_138();
  } else {
    _if_67_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_67_end() {
  _if_68();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_138() {
  _when_61();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_138_end() {
  _if_67_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_61() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_139(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(61);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_61_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_138_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_139(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_110(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_139_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_61_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_110(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_110()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 539 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 12708 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_139_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_68() {
  if (thisIndex.x > 0) {
    _slist_140();
  } else {
    _if_68_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_68_end() {
  _if_69();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_140() {
  _when_62();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_140_end() {
  _if_68_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_62() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_141(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(62);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_62_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_140_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_141(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_111(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_141_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_62_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_111(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_111()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 544 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 12799 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_141_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_69() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_142();
  } else {
    _if_69_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_69_end() {
  _serial_115();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_142() {
  _serial_112();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_142_end() {
  _if_69_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_112() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_112()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 548 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, EAST_WEST, grid); 
#line 12847 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_70();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_70() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_143();
  } else {
    _if_70_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_70_end() {
  _if_71();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_143() {
  _when_63();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_143_end() {
  _if_70_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_63() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_144(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(63);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_63_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_143_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_144(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_113(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_144_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_63_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_113(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_113()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 551 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 12937 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_144_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_71() {
  if (thisIndex.x > 0) {
    _slist_145();
  } else {
    _if_71_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_71_end() {
  _slist_142_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_145() {
  _when_64();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_145_end() {
  _if_71_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_64() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_146(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(64);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_64_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_145_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_146(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_114(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_146_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_64_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_114(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_114()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 556 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13028 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_146_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_115() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_115()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 560 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, NORTH_SOUTH, grid); 
#line 13044 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_72();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_72() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_147();
  } else {
    _if_72_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_72_end() {
  _if_73();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_147() {
  _when_65();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_147_end() {
  _if_72_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_65() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_148(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(65);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_65_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_147_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_148(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_116(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_148_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_65_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_116(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_116()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 563 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13134 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_148_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_73() {
  if (thisIndex.y > 0) {
    _slist_149();
  } else {
    _if_73_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_73_end() {
  _if_74();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_149() {
  _when_66();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_149_end() {
  _if_73_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_66() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_150(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(66);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_66_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_149_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_150(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_117(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_150_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_66_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_117(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_117()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 568 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13225 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_150_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_74() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_151();
  } else {
    _if_74_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_74_end() {
  _serial_121();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_151() {
  _serial_118();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_151_end() {
  _if_74_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_118() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_118()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 572 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, NORTH_SOUTH, grid); 
#line 13273 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_75();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_75() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_152();
  } else {
    _if_75_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_75_end() {
  _if_76();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_152() {
  _when_67();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_152_end() {
  _if_75_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_67() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_153(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(67);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_67_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_152_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_153(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_119(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_153_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_67_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_119(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_119()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 575 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13363 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_153_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_76() {
  if (thisIndex.y > 0) {
    _slist_154();
  } else {
    _if_76_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_76_end() {
  _slist_151_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_154() {
  _when_68();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_154_end() {
  _if_76_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_68() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_155(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(68);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_68_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_154_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_155(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_120(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_155_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_68_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_120(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_120()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 580 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13454 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_155_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_121() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_121()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 584 "timestepping.ci"
 send_neighbor_data(K_ODD_GFS, TOP_BOTTOM, grid); 
#line 13470 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_77();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_77() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_156();
  } else {
    _if_77_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_77_end() {
  _if_78();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_156() {
  _when_69();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_156_end() {
  _if_77_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_69() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_157(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(69);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_69_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_156_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_157(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_122(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_157_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_69_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_122(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_122()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 587 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13560 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_157_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_78() {
  if (thisIndex.z > 0) {
    _slist_158();
  } else {
    _if_78_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_78_end() {
  _if_79();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_158() {
  _when_70();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_158_end() {
  _if_78_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_70() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_159(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(70);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_70_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_158_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_159(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_123(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_159_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_70_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_123(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_123()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 592 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13651 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_159_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_79() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_160();
  } else {
    _if_79_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_79_end() {
  _serial_127();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_160() {
  _serial_124();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_160_end() {
  _if_79_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_124() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_124()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 596 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, TOP_BOTTOM, grid); 
#line 13699 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_80();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_80() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_161();
  } else {
    _if_80_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_80_end() {
  _if_81();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_161() {
  _when_71();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_161_end() {
  _if_80_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_71() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_162(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(71);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_71_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_161_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_162(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_125(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_162_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_71_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_125(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_125()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 599 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13789 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_162_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_81() {
  if (thisIndex.z > 0) {
    _slist_163();
  } else {
    _if_81_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_81_end() {
  _slist_160_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_163() {
  _when_72();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_163_end() {
  _if_81_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_72() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_164(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(72);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_72_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_163_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_164(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_126(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_164_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_72_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_126(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_126()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 604 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 13880 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_164_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_127() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_127()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 608 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K4, MOL_PRE_RK_UPDATE); 
#line 13896 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_82();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_82() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_165();
  } else {
    _if_82_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_82_end() {
  _if_83();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_165() {
  _serial_128();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_165_end() {
  _if_82_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_128() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_128()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 610 "timestepping.ci"
 send_nonlocalinnerbc_data(K_EVEN_GFS, grid); 
#line 13943 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_165_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_83() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_166();
  } else {
    _if_83_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_83_end() {
  _serial_131();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_166() {
  _for_7();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_166_end() {
  _if_83_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_7() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_167();
  } else {
    _serial_130();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_7_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_167();
  } else {
    _serial_130();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_167() {
  _when_73();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_167_end() {
  _for_7_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_73() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(27, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_168(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(73);
    c->anyEntries.push_back(27);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_73_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _slist_167_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_168(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _serial_129(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_168_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  _when_73_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_129(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_129()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 615 "timestepping.ci"

                set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                type_gfs_nonlocal_innerbc = type_gfs;
              
#line 14075 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_168_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_130() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_130()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 621 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 14091 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_166_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_131() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_131()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 623 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K4, MOL_RK_UPDATE); 
#line 14106 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_132();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_132() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_132()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 624 "timestepping.ci"
 MoL_step_forward_in_time(&commondata, griddata_chare, time_start, RK_SUBSTEP_K4, MOL_POST_RK_UPDATE); 
#line 14121 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_84();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_84() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_169();
  } else {
    _if_84_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_84_end() {
  _serial_136();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_169() {
  _if_85();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_169_end() {
  _if_84_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_85() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_dst_chares > 0) {
    _slist_170();
  } else {
    _if_85_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_85_end() {
  _if_86();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_170() {
  _serial_133();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_170_end() {
  _if_85_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_133() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_133()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 627 "timestepping.ci"
 send_nonlocalinnerbc_data(AUXEVOL_GFS, grid); 
#line 14200 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_170_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_86() {
  if (griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares > 0) {
    _slist_171();
  } else {
    _if_86_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_86_end() {
  _slist_169_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_171() {
  _for_8();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_171_end() {
  _if_86_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_8() {
  iter = 0;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_172();
  } else {
    _serial_135();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_for_8_end() {
   iter++;
  if ( iter < griddata_chare[grid].nonlocalinnerbcstruct.tot_num_src_chares) {
    _slist_172();
  } else {
    _serial_135();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_172() {
  _when_74();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_172_end() {
  _for_8_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_74() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(26, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_173(static_cast<Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(74);
    c->anyEntries.push_back(26);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_74_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _slist_172_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_173(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _serial_134(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_173_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  _when_74_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_134(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_134()), CkMyPe(), 0, &projID, this); 
  {
    int& src_chare_idx3 = gen0->getP0();
    int& type_gfs = gen0->getP1();
    int& len_tmpBuffer = gen0->getP2();
    REAL*& tmpBuffer = gen0->getP3();
    { // begin serial block
#line 632 "timestepping.ci"

                  set_tmpBuffer_innerbc_receiv(src_chare_idx3, len_tmpBuffer, tmpBuffer, grid);
                  type_gfs_nonlocal_innerbc = type_gfs;
                
#line 14332 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_173_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_135() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_135()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 638 "timestepping.ci"
 process_nonlocalinnerbc(type_gfs_nonlocal_innerbc, grid); 
#line 14348 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_171_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_136() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_136()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 641 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, EAST_WEST, grid); 
#line 14363 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_87();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_87() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_174();
  } else {
    _if_87_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_87_end() {
  _if_88();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_174() {
  _when_75();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_174_end() {
  _if_87_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_75() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_175(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(75);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_75_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_174_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_175(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_137(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_175_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_75_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_137(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_137()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 644 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14453 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_175_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_88() {
  if (thisIndex.x > 0) {
    _slist_176();
  } else {
    _if_88_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_88_end() {
  _if_89();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_176() {
  _when_76();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_176_end() {
  _if_88_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_76() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_177(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(76);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_76_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_176_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_177(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_138(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_177_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_76_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_138(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_138()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 649 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14544 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_177_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_89() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_178();
  } else {
    _if_89_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_89_end() {
  _serial_142();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_178() {
  _serial_139();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_178_end() {
  _if_89_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_139() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_139()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 653 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, EAST_WEST, grid); 
#line 14592 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_90();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_90() {
  if (thisIndex.x < commondata.Nchare0 - 1) {
    _slist_179();
  } else {
    _if_90_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_90_end() {
  _if_91();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_179() {
  _when_77();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_179_end() {
  _if_90_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_77() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(2, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_180(static_cast<Closure_Timestepping::east_ghost_20_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(77);
    c->anyEntries.push_back(2);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_77_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _slist_179_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_180(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _serial_140(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_180_end(Closure_Timestepping::east_ghost_20_closure* gen0) {
  _when_77_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_140(Closure_Timestepping::east_ghost_20_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_140()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 656 "timestepping.ci"
 process_ghost(EAST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14682 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_180_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_91() {
  if (thisIndex.x > 0) {
    _slist_181();
  } else {
    _if_91_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_91_end() {
  _slist_178_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_181() {
  _when_78();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_181_end() {
  _if_91_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_78() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(3, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_182(static_cast<Closure_Timestepping::west_ghost_21_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(78);
    c->anyEntries.push_back(3);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_78_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _slist_181_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_182(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _serial_141(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_182_end(Closure_Timestepping::west_ghost_21_closure* gen0) {
  _when_78_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_141(Closure_Timestepping::west_ghost_21_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_141()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 661 "timestepping.ci"
 process_ghost(WEST_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14773 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_182_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_142() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_142()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 665 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, NORTH_SOUTH, grid); 
#line 14789 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_92();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_92() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_183();
  } else {
    _if_92_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_92_end() {
  _if_93();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_183() {
  _when_79();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_183_end() {
  _if_92_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_79() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_184(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(79);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_79_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_183_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_184(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_143(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_184_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_79_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_143(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_143()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 668 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14879 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_184_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_93() {
  if (thisIndex.y > 0) {
    _slist_185();
  } else {
    _if_93_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_93_end() {
  _if_94();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_185() {
  _when_80();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_185_end() {
  _if_93_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_80() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_186(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(80);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_80_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_185_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_186(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_144(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_186_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_80_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_144(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_144()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 673 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 14970 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_186_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_94() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_187();
  } else {
    _if_94_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_94_end() {
  _serial_148();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_187() {
  _serial_145();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_187_end() {
  _if_94_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_145() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_145()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 677 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, NORTH_SOUTH, grid); 
#line 15018 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_95();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_95() {
  if (thisIndex.y < commondata.Nchare1 - 1) {
    _slist_188();
  } else {
    _if_95_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_95_end() {
  _if_96();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_188() {
  _when_81();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_188_end() {
  _if_95_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_81() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(4, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_189(static_cast<Closure_Timestepping::north_ghost_22_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(81);
    c->anyEntries.push_back(4);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_81_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _slist_188_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_189(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _serial_146(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_189_end(Closure_Timestepping::north_ghost_22_closure* gen0) {
  _when_81_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_146(Closure_Timestepping::north_ghost_22_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_146()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 680 "timestepping.ci"
 process_ghost(NORTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15108 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_189_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_96() {
  if (thisIndex.y > 0) {
    _slist_190();
  } else {
    _if_96_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_96_end() {
  _slist_187_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_190() {
  _when_82();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_190_end() {
  _if_96_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_82() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(5, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_191(static_cast<Closure_Timestepping::south_ghost_23_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(82);
    c->anyEntries.push_back(5);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_82_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _slist_190_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_191(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _serial_147(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_191_end(Closure_Timestepping::south_ghost_23_closure* gen0) {
  _when_82_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_147(Closure_Timestepping::south_ghost_23_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_147()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 685 "timestepping.ci"
 process_ghost(SOUTH_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15199 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_191_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_148() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_148()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 689 "timestepping.ci"
 send_neighbor_data(Y_N_GFS, TOP_BOTTOM, grid); 
#line 15215 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_97();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_97() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_192();
  } else {
    _if_97_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_97_end() {
  _if_98();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_192() {
  _when_83();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_192_end() {
  _if_97_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_83() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_193(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(83);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_83_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_192_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_193(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_149(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_193_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_83_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_149(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_149()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 692 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15305 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_193_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_98() {
  if (thisIndex.z > 0) {
    _slist_194();
  } else {
    _if_98_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_98_end() {
  _if_99();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_194() {
  _when_84();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_194_end() {
  _if_98_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_84() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_195(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(84);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_84_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_194_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_195(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_150(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_195_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_84_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_150(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_150()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 697 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15396 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_195_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_99() {
  if (griddata_chare[grid].gridfuncs.num_auxevol_gfs_to_sync > 0) {
    _slist_196();
  } else {
    _if_99_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_99_end() {
  _serial_154();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_196() {
  _serial_151();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_196_end() {
  _if_99_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_151() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_151()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 701 "timestepping.ci"
 send_neighbor_data(AUXEVOL_GFS, TOP_BOTTOM, grid); 
#line 15444 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _if_100();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_100() {
  if (thisIndex.z < commondata.Nchare2 - 1) {
    _slist_197();
  } else {
    _if_100_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_100_end() {
  _if_101();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_197() {
  _when_85();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_197_end() {
  _if_100_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_85() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(6, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_198(static_cast<Closure_Timestepping::top_ghost_24_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(85);
    c->anyEntries.push_back(6);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_85_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _slist_197_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_198(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _serial_152(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_198_end(Closure_Timestepping::top_ghost_24_closure* gen0) {
  _when_85_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_152(Closure_Timestepping::top_ghost_24_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_152()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 704 "timestepping.ci"
 process_ghost(TOP_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15534 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_198_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_101() {
  if (thisIndex.z > 0) {
    _slist_199();
  } else {
    _if_101_end();
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_if_101_end() {
  _slist_196_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_199() {
  _when_86();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_199_end() {
  _if_101_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
SDAG::Continuation* Timestepping::_when_86() {
  SDAG::Buffer* buf0 = __dep->tryFindMessage(7, false, 0, 0);
  if (buf0) {
    __dep->removeMessage(buf0);
    _slist_200(static_cast<Closure_Timestepping::bottom_ghost_25_closure*>(buf0->cl));
    delete buf0;
    return 0;
  } else {
    SDAG::Continuation* c = new SDAG::Continuation(86);
    c->anyEntries.push_back(7);
    __dep->reg(c);
    return c;
  }
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_when_86_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _slist_199_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_200(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _serial_153(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_200_end(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  _when_86_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_153(Closure_Timestepping::bottom_ghost_25_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_153()), CkMyPe(), 0, &projID, this); 
  {
    int& type_gfs = gen0->getP0();
    int& len_tmpBuffer = gen0->getP1();
    REAL*& tmpBuffer = gen0->getP2();
    { // begin serial block
#line 709 "timestepping.ci"
 process_ghost(BOTTOM_GHOST, type_gfs, len_tmpBuffer, tmpBuffer, grid); 
#line 15625 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_200_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_154() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_154()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 713 "timestepping.ci"

          stop_conditions_check(&commondata);
          if (commondata.stop_relaxation) {
            mainProxy.done();
          }
        
#line 15646 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _serial_155();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_155() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_155()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 719 "timestepping.ci"

          commondata.time = (REAL)(commondata.nn + 1) * commondata.dt;
          commondata.nn++;
        
#line 15664 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_33_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_156() {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_156()), CkMyPe(), 0, &projID, this); 
  { // begin serial block
#line 724 "timestepping.ci"
 mainProxy.done(); 
#line 15679 "timestepping.def.h"
  } // end serial block
  _TRACE_END_EXECUTE(); 
  _slist_0_end();
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::diagnostics_ckio(Ck::IO::Session token, const int which_diagnostics_part){
  Closure_Timestepping::diagnostics_ckio_19_closure* genClosure = new Closure_Timestepping::diagnostics_ckio_19_closure();
  genClosure->getP0() = token;
  genClosure->getP1() = which_diagnostics_part;
  diagnostics_ckio(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::diagnostics_ckio(Closure_Timestepping::diagnostics_ckio_19_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _slist_201(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::diagnostics_ckio_end(Closure_Timestepping::diagnostics_ckio_19_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_201(Closure_Timestepping::diagnostics_ckio_19_closure* gen0) {
  _serial_157(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_slist_201_end(Closure_Timestepping::diagnostics_ckio_19_closure* gen0) {
  diagnostics_ckio_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_157(Closure_Timestepping::diagnostics_ckio_19_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_157()), CkMyPe(), 0, &projID, this); 
  {
    Ck::IO::Session& token = gen0->getP0();
    int& which_diagnostics_part = gen0->getP1();
    { // begin serial block
#line 739 "timestepping.ci"

        printf("entry void diagnostics_ckio(Ck::IO::Session token, const int which_diagnostics_part = %d) ", which_diagnostics_part);
        const int thisIndex_arr[3] = {thisIndex.x, thisIndex.y, thisIndex.z};
        diagnostics(&commondata, griddata, griddata_chare, diagnostic_gfs, thisIndex_arr, which_grid_diagnostics, token, which_diagnostics_part);
      
#line 15742 "timestepping.def.h"
    } // end serial block
  }
  _TRACE_END_EXECUTE(); 
  _slist_201_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL * tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _serial_158(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_158(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_158()), CkMyPe(), 0, &projID, this); 
  _TRACE_END_EXECUTE(); 
  receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL * tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_y_n_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _serial_159(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_159(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_159()), CkMyPe(), 0, &projID, this); 
  _TRACE_END_EXECUTE(); 
  receiv_nonlocalinnerbc_data_y_n_gfs_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL * tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_diagnostic_output_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _serial_160(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_160(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_160()), CkMyPe(), 0, &projID, this); 
  _TRACE_END_EXECUTE(); 
  receiv_nonlocalinnerbc_data_diagnostic_output_gfs_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL * tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _serial_161(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_161(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_161()), CkMyPe(), 0, &projID, this); 
  _TRACE_END_EXECUTE(); 
  receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL * tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* gen0) {
  _TRACE_END_EXECUTE(); 
  if (!__dep.get()) _sdag_init();
  _serial_162(gen0);
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_end(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* gen0) {
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_serial_162(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure* gen0) {
  CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
  _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, (_sdag_idx_Timestepping_serial_162()), CkMyPe(), 0, &projID, this); 
  _TRACE_END_EXECUTE(); 
  receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_end(gen0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, int *globalidx3_srcpts){
  Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure();
  genClosure->getP0() = idx3_of_sendingchare;
  genClosure->getP1() = num_srcpts;
  genClosure->getP2() = globalidx3_srcpts;
  receiv_nonlocalinnerbc_idx3srcpt_tosend(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(0, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(0);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_0(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_wavespeed_at_outer_boundary(REAL wavespeed_at_outer_boundary){
  Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* genClosure = new Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure();
  genClosure->getP0() = wavespeed_at_outer_boundary;
  receiv_wavespeed_at_outer_boundary(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_wavespeed_at_outer_boundary(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(1, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(1);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_1(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::east_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::east_ghost_20_closure* genClosure = new Closure_Timestepping::east_ghost_20_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  east_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::east_ghost(Closure_Timestepping::east_ghost_20_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(2, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(2);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 2:
      _when_2(
      );
    break;
    case 4:
      _when_4(
      );
    break;
    case 33:
      _when_33(
      );
    break;
    case 35:
      _when_35(
      );
    break;
    case 47:
      _when_47(
      );
    break;
    case 49:
      _when_49(
      );
    break;
    case 61:
      _when_61(
      );
    break;
    case 63:
      _when_63(
      );
    break;
    case 75:
      _when_75(
      );
    break;
    case 77:
      _when_77(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::west_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::west_ghost_21_closure* genClosure = new Closure_Timestepping::west_ghost_21_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  west_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::west_ghost(Closure_Timestepping::west_ghost_21_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(3, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(3);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 3:
      _when_3(
      );
    break;
    case 5:
      _when_5(
      );
    break;
    case 34:
      _when_34(
      );
    break;
    case 36:
      _when_36(
      );
    break;
    case 48:
      _when_48(
      );
    break;
    case 50:
      _when_50(
      );
    break;
    case 62:
      _when_62(
      );
    break;
    case 64:
      _when_64(
      );
    break;
    case 76:
      _when_76(
      );
    break;
    case 78:
      _when_78(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::north_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::north_ghost_22_closure* genClosure = new Closure_Timestepping::north_ghost_22_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  north_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::north_ghost(Closure_Timestepping::north_ghost_22_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(4, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(4);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 6:
      _when_6(
      );
    break;
    case 8:
      _when_8(
      );
    break;
    case 37:
      _when_37(
      );
    break;
    case 39:
      _when_39(
      );
    break;
    case 51:
      _when_51(
      );
    break;
    case 53:
      _when_53(
      );
    break;
    case 65:
      _when_65(
      );
    break;
    case 67:
      _when_67(
      );
    break;
    case 79:
      _when_79(
      );
    break;
    case 81:
      _when_81(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::south_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::south_ghost_23_closure* genClosure = new Closure_Timestepping::south_ghost_23_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  south_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::south_ghost(Closure_Timestepping::south_ghost_23_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(5, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(5);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 7:
      _when_7(
      );
    break;
    case 9:
      _when_9(
      );
    break;
    case 38:
      _when_38(
      );
    break;
    case 40:
      _when_40(
      );
    break;
    case 52:
      _when_52(
      );
    break;
    case 54:
      _when_54(
      );
    break;
    case 66:
      _when_66(
      );
    break;
    case 68:
      _when_68(
      );
    break;
    case 80:
      _when_80(
      );
    break;
    case 82:
      _when_82(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::top_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::top_ghost_24_closure* genClosure = new Closure_Timestepping::top_ghost_24_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  top_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::top_ghost(Closure_Timestepping::top_ghost_24_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(6, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(6);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 10:
      _when_10(
      );
    break;
    case 12:
      _when_12(
      );
    break;
    case 41:
      _when_41(
      );
    break;
    case 43:
      _when_43(
      );
    break;
    case 55:
      _when_55(
      );
    break;
    case 57:
      _when_57(
      );
    break;
    case 69:
      _when_69(
      );
    break;
    case 71:
      _when_71(
      );
    break;
    case 83:
      _when_83(
      );
    break;
    case 85:
      _when_85(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::bottom_ghost(int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::bottom_ghost_25_closure* genClosure = new Closure_Timestepping::bottom_ghost_25_closure();
  genClosure->getP0() = type_gfs;
  genClosure->getP1() = len_tmpBuffer;
  genClosure->getP2() = tmpBuffer;
  bottom_ghost(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::bottom_ghost(Closure_Timestepping::bottom_ghost_25_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(7, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(7);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 11:
      _when_11(
      );
    break;
    case 13:
      _when_13(
      );
    break;
    case 42:
      _when_42(
      );
    break;
    case 44:
      _when_44(
      );
    break;
    case 56:
      _when_56(
      );
    break;
    case 58:
      _when_58(
      );
    break;
    case 70:
      _when_70(
      );
    break;
    case 72:
      _when_72(
      );
    break;
    case 84:
      _when_84(
      );
    break;
    case 86:
      _when_86(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::ready_1d_y(Ck::IO::FileReadyMsg* m_1d_y_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_y_msg);
  __dep->pushBuffer(8, new SDAG::MsgClosure(m_1d_y_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(8);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_14(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::start_write_1d_y(Ck::IO::SessionReadyMsg* m_1d_y_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_y_msg);
  __dep->pushBuffer(9, new SDAG::MsgClosure(m_1d_y_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(9);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_15(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::test_written_1d_y(CkReductionMsg* m_1d_y_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_y_msg);
  __dep->pushBuffer(10, new SDAG::MsgClosure(m_1d_y_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(10);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_16(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::closed_1d_y(CkReductionMsg* m_1d_y_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_y_msg);
  __dep->pushBuffer(11, new SDAG::MsgClosure(m_1d_y_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(11);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_17(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::ready_1d_z(Ck::IO::FileReadyMsg* m_1d_z_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_z_msg);
  __dep->pushBuffer(12, new SDAG::MsgClosure(m_1d_z_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(12);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_18(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::start_write_1d_z(Ck::IO::SessionReadyMsg* m_1d_z_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_z_msg);
  __dep->pushBuffer(13, new SDAG::MsgClosure(m_1d_z_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(13);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_19(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::test_written_1d_z(CkReductionMsg* m_1d_z_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_z_msg);
  __dep->pushBuffer(14, new SDAG::MsgClosure(m_1d_z_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(14);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_20(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::closed_1d_z(CkReductionMsg* m_1d_z_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_1d_z_msg);
  __dep->pushBuffer(15, new SDAG::MsgClosure(m_1d_z_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(15);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_21(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::ready_2d_xy(Ck::IO::FileReadyMsg* m_2d_xy_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_xy_msg);
  __dep->pushBuffer(16, new SDAG::MsgClosure(m_2d_xy_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(16);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_22(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::start_write_2d_xy(Ck::IO::SessionReadyMsg* m_2d_xy_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_xy_msg);
  __dep->pushBuffer(17, new SDAG::MsgClosure(m_2d_xy_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(17);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_23(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::test_written_2d_xy(CkReductionMsg* m_2d_xy_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_xy_msg);
  __dep->pushBuffer(18, new SDAG::MsgClosure(m_2d_xy_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(18);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_24(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::closed_2d_xy(CkReductionMsg* m_2d_xy_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_xy_msg);
  __dep->pushBuffer(19, new SDAG::MsgClosure(m_2d_xy_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(19);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_25(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::ready_2d_yz(Ck::IO::FileReadyMsg* m_2d_yz_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_yz_msg);
  __dep->pushBuffer(20, new SDAG::MsgClosure(m_2d_yz_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(20);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_26(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::start_write_2d_yz(Ck::IO::SessionReadyMsg* m_2d_yz_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_yz_msg);
  __dep->pushBuffer(21, new SDAG::MsgClosure(m_2d_yz_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(21);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_27(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::test_written_2d_yz(CkReductionMsg* m_2d_yz_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_yz_msg);
  __dep->pushBuffer(22, new SDAG::MsgClosure(m_2d_yz_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(22);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_28(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::closed_2d_yz(CkReductionMsg* m_2d_yz_msg){
  if (!__dep.get()) _sdag_init();
  CkReferenceMsg(m_2d_yz_msg);
  __dep->pushBuffer(23, new SDAG::MsgClosure(m_2d_yz_msg));
  SDAG::Continuation* c = __dep->tryFindContinuation(23);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_29(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::continue_timestepping(){
  Closure_Timestepping::continue_timestepping_26_closure* genClosure = new Closure_Timestepping::continue_timestepping_26_closure();
  continue_timestepping(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::continue_timestepping(Closure_Timestepping::continue_timestepping_26_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(24, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(24);
  if (c) {
    _TRACE_END_EXECUTE(); 
    _when_30(
    );
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_k_odd_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(25, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(25);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 31:
      _when_31(
      );
    break;
    case 59:
      _when_59(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_auxevol_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(26, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(26);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 32:
      _when_32(
      );
    break;
    case 46:
      _when_46(
      );
    break;
    case 60:
      _when_60(
      );
    break;
    case 74:
      _when_74(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, REAL *tmpBuffer){
  Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* genClosure = new Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure();
  genClosure->getP0() = src_chare_idx3;
  genClosure->getP1() = type_gfs;
  genClosure->getP2() = len_tmpBuffer;
  genClosure->getP3() = tmpBuffer;
  receiv_nonlocalinnerbc_data_k_even_gfs(genClosure);
  genClosure->deref();
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure* genClosure){
  if (!__dep.get()) _sdag_init();
  __dep->pushBuffer(27, genClosure);
  SDAG::Continuation* c = __dep->tryFindContinuation(27);
  if (c) {
    _TRACE_END_EXECUTE(); 
    switch(c->whenID) {
    case 45:
      _when_45(
      );
    break;
    case 73:
      _when_73(
      );
    break;
    }
    CmiObjId projID = this->ckGetArrayIndex().getProjectionID();
    _TRACE_BEGIN_EXECUTE_DETAILED(-1, -1, _sdagEP, CkMyPe(), 0, &projID, this); 
    delete c;
  }
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::_sdag_init() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  __dep.reset(new SDAG::Dependency(28,87));
  __dep->addDepends(0,0);
  __dep->addDepends(1,1);
  __dep->addDepends(2,2);
  __dep->addDepends(4,2);
  __dep->addDepends(33,2);
  __dep->addDepends(35,2);
  __dep->addDepends(47,2);
  __dep->addDepends(49,2);
  __dep->addDepends(61,2);
  __dep->addDepends(63,2);
  __dep->addDepends(75,2);
  __dep->addDepends(77,2);
  __dep->addDepends(3,3);
  __dep->addDepends(5,3);
  __dep->addDepends(34,3);
  __dep->addDepends(36,3);
  __dep->addDepends(48,3);
  __dep->addDepends(50,3);
  __dep->addDepends(62,3);
  __dep->addDepends(64,3);
  __dep->addDepends(76,3);
  __dep->addDepends(78,3);
  __dep->addDepends(6,4);
  __dep->addDepends(8,4);
  __dep->addDepends(37,4);
  __dep->addDepends(39,4);
  __dep->addDepends(51,4);
  __dep->addDepends(53,4);
  __dep->addDepends(65,4);
  __dep->addDepends(67,4);
  __dep->addDepends(79,4);
  __dep->addDepends(81,4);
  __dep->addDepends(7,5);
  __dep->addDepends(9,5);
  __dep->addDepends(38,5);
  __dep->addDepends(40,5);
  __dep->addDepends(52,5);
  __dep->addDepends(54,5);
  __dep->addDepends(66,5);
  __dep->addDepends(68,5);
  __dep->addDepends(80,5);
  __dep->addDepends(82,5);
  __dep->addDepends(10,6);
  __dep->addDepends(12,6);
  __dep->addDepends(41,6);
  __dep->addDepends(43,6);
  __dep->addDepends(55,6);
  __dep->addDepends(57,6);
  __dep->addDepends(69,6);
  __dep->addDepends(71,6);
  __dep->addDepends(83,6);
  __dep->addDepends(85,6);
  __dep->addDepends(11,7);
  __dep->addDepends(13,7);
  __dep->addDepends(42,7);
  __dep->addDepends(44,7);
  __dep->addDepends(56,7);
  __dep->addDepends(58,7);
  __dep->addDepends(70,7);
  __dep->addDepends(72,7);
  __dep->addDepends(84,7);
  __dep->addDepends(86,7);
  __dep->addDepends(14,8);
  __dep->addDepends(15,9);
  __dep->addDepends(16,10);
  __dep->addDepends(17,11);
  __dep->addDepends(18,12);
  __dep->addDepends(19,13);
  __dep->addDepends(20,14);
  __dep->addDepends(21,15);
  __dep->addDepends(22,16);
  __dep->addDepends(23,17);
  __dep->addDepends(24,18);
  __dep->addDepends(25,19);
  __dep->addDepends(26,20);
  __dep->addDepends(27,21);
  __dep->addDepends(28,22);
  __dep->addDepends(29,23);
  __dep->addDepends(30,24);
  __dep->addDepends(31,25);
  __dep->addDepends(59,25);
  __dep->addDepends(32,26);
  __dep->addDepends(46,26);
  __dep->addDepends(60,26);
  __dep->addDepends(74,26);
  __dep->addDepends(45,27);
  __dep->addDepends(73,27);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::__sdag_init() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
void Timestepping::_sdag_pup(PUP::er &p) {  // Potentially missing Timestepping_SDAG_CODE in your class definition?
  p|__dep;
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void Timestepping::__sdag_register() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  (void)_sdag_idx_Timestepping_serial_0();
  (void)_sdag_idx_Timestepping_serial_1();
  (void)_sdag_idx_Timestepping_serial_2();
  (void)_sdag_idx_Timestepping_serial_3();
  (void)_sdag_idx_Timestepping_serial_4();
  (void)_sdag_idx_Timestepping_serial_5();
  (void)_sdag_idx_Timestepping_serial_6();
  (void)_sdag_idx_Timestepping_serial_7();
  (void)_sdag_idx_Timestepping_serial_8();
  (void)_sdag_idx_Timestepping_serial_9();
  (void)_sdag_idx_Timestepping_serial_10();
  (void)_sdag_idx_Timestepping_serial_11();
  (void)_sdag_idx_Timestepping_serial_12();
  (void)_sdag_idx_Timestepping_serial_13();
  (void)_sdag_idx_Timestepping_serial_14();
  (void)_sdag_idx_Timestepping_serial_15();
  (void)_sdag_idx_Timestepping_serial_16();
  (void)_sdag_idx_Timestepping_serial_17();
  (void)_sdag_idx_Timestepping_serial_18();
  (void)_sdag_idx_Timestepping_serial_19();
  (void)_sdag_idx_Timestepping_serial_20();
  (void)_sdag_idx_Timestepping_serial_21();
  (void)_sdag_idx_Timestepping_serial_22();
  (void)_sdag_idx_Timestepping_serial_23();
  (void)_sdag_idx_Timestepping_serial_24();
  (void)_sdag_idx_Timestepping_serial_25();
  (void)_sdag_idx_Timestepping_serial_27();
  (void)_sdag_idx_Timestepping_serial_28();
  (void)_sdag_idx_Timestepping_serial_29();
  (void)_sdag_idx_Timestepping_serial_30();
  (void)_sdag_idx_Timestepping_serial_31();
  (void)_sdag_idx_Timestepping_serial_32();
  (void)_sdag_idx_Timestepping_serial_33();
  (void)_sdag_idx_Timestepping_serial_34();
  (void)_sdag_idx_Timestepping_serial_35();
  (void)_sdag_idx_Timestepping_serial_36();
  (void)_sdag_idx_Timestepping_serial_37();
  (void)_sdag_idx_Timestepping_serial_38();
  (void)_sdag_idx_Timestepping_serial_39();
  (void)_sdag_idx_Timestepping_serial_40();
  (void)_sdag_idx_Timestepping_serial_41();
  (void)_sdag_idx_Timestepping_serial_42();
  (void)_sdag_idx_Timestepping_serial_43();
  (void)_sdag_idx_Timestepping_serial_44();
  (void)_sdag_idx_Timestepping_serial_26();
  (void)_sdag_idx_Timestepping_serial_45();
  (void)_sdag_idx_Timestepping_serial_46();
  (void)_sdag_idx_Timestepping_serial_47();
  (void)_sdag_idx_Timestepping_serial_48();
  (void)_sdag_idx_Timestepping_serial_49();
  (void)_sdag_idx_Timestepping_serial_50();
  (void)_sdag_idx_Timestepping_serial_51();
  (void)_sdag_idx_Timestepping_serial_52();
  (void)_sdag_idx_Timestepping_serial_53();
  (void)_sdag_idx_Timestepping_serial_54();
  (void)_sdag_idx_Timestepping_serial_55();
  (void)_sdag_idx_Timestepping_serial_56();
  (void)_sdag_idx_Timestepping_serial_57();
  (void)_sdag_idx_Timestepping_serial_58();
  (void)_sdag_idx_Timestepping_serial_59();
  (void)_sdag_idx_Timestepping_serial_60();
  (void)_sdag_idx_Timestepping_serial_61();
  (void)_sdag_idx_Timestepping_serial_62();
  (void)_sdag_idx_Timestepping_serial_63();
  (void)_sdag_idx_Timestepping_serial_64();
  (void)_sdag_idx_Timestepping_serial_65();
  (void)_sdag_idx_Timestepping_serial_66();
  (void)_sdag_idx_Timestepping_serial_67();
  (void)_sdag_idx_Timestepping_serial_68();
  (void)_sdag_idx_Timestepping_serial_69();
  (void)_sdag_idx_Timestepping_serial_70();
  (void)_sdag_idx_Timestepping_serial_71();
  (void)_sdag_idx_Timestepping_serial_72();
  (void)_sdag_idx_Timestepping_serial_73();
  (void)_sdag_idx_Timestepping_serial_74();
  (void)_sdag_idx_Timestepping_serial_75();
  (void)_sdag_idx_Timestepping_serial_76();
  (void)_sdag_idx_Timestepping_serial_77();
  (void)_sdag_idx_Timestepping_serial_78();
  (void)_sdag_idx_Timestepping_serial_79();
  (void)_sdag_idx_Timestepping_serial_80();
  (void)_sdag_idx_Timestepping_serial_81();
  (void)_sdag_idx_Timestepping_serial_82();
  (void)_sdag_idx_Timestepping_serial_83();
  (void)_sdag_idx_Timestepping_serial_84();
  (void)_sdag_idx_Timestepping_serial_85();
  (void)_sdag_idx_Timestepping_serial_86();
  (void)_sdag_idx_Timestepping_serial_87();
  (void)_sdag_idx_Timestepping_serial_88();
  (void)_sdag_idx_Timestepping_serial_89();
  (void)_sdag_idx_Timestepping_serial_90();
  (void)_sdag_idx_Timestepping_serial_91();
  (void)_sdag_idx_Timestepping_serial_92();
  (void)_sdag_idx_Timestepping_serial_93();
  (void)_sdag_idx_Timestepping_serial_94();
  (void)_sdag_idx_Timestepping_serial_95();
  (void)_sdag_idx_Timestepping_serial_96();
  (void)_sdag_idx_Timestepping_serial_97();
  (void)_sdag_idx_Timestepping_serial_98();
  (void)_sdag_idx_Timestepping_serial_99();
  (void)_sdag_idx_Timestepping_serial_100();
  (void)_sdag_idx_Timestepping_serial_101();
  (void)_sdag_idx_Timestepping_serial_102();
  (void)_sdag_idx_Timestepping_serial_103();
  (void)_sdag_idx_Timestepping_serial_104();
  (void)_sdag_idx_Timestepping_serial_105();
  (void)_sdag_idx_Timestepping_serial_106();
  (void)_sdag_idx_Timestepping_serial_107();
  (void)_sdag_idx_Timestepping_serial_108();
  (void)_sdag_idx_Timestepping_serial_109();
  (void)_sdag_idx_Timestepping_serial_110();
  (void)_sdag_idx_Timestepping_serial_111();
  (void)_sdag_idx_Timestepping_serial_112();
  (void)_sdag_idx_Timestepping_serial_113();
  (void)_sdag_idx_Timestepping_serial_114();
  (void)_sdag_idx_Timestepping_serial_115();
  (void)_sdag_idx_Timestepping_serial_116();
  (void)_sdag_idx_Timestepping_serial_117();
  (void)_sdag_idx_Timestepping_serial_118();
  (void)_sdag_idx_Timestepping_serial_119();
  (void)_sdag_idx_Timestepping_serial_120();
  (void)_sdag_idx_Timestepping_serial_121();
  (void)_sdag_idx_Timestepping_serial_122();
  (void)_sdag_idx_Timestepping_serial_123();
  (void)_sdag_idx_Timestepping_serial_124();
  (void)_sdag_idx_Timestepping_serial_125();
  (void)_sdag_idx_Timestepping_serial_126();
  (void)_sdag_idx_Timestepping_serial_127();
  (void)_sdag_idx_Timestepping_serial_128();
  (void)_sdag_idx_Timestepping_serial_129();
  (void)_sdag_idx_Timestepping_serial_130();
  (void)_sdag_idx_Timestepping_serial_131();
  (void)_sdag_idx_Timestepping_serial_132();
  (void)_sdag_idx_Timestepping_serial_133();
  (void)_sdag_idx_Timestepping_serial_134();
  (void)_sdag_idx_Timestepping_serial_135();
  (void)_sdag_idx_Timestepping_serial_136();
  (void)_sdag_idx_Timestepping_serial_137();
  (void)_sdag_idx_Timestepping_serial_138();
  (void)_sdag_idx_Timestepping_serial_139();
  (void)_sdag_idx_Timestepping_serial_140();
  (void)_sdag_idx_Timestepping_serial_141();
  (void)_sdag_idx_Timestepping_serial_142();
  (void)_sdag_idx_Timestepping_serial_143();
  (void)_sdag_idx_Timestepping_serial_144();
  (void)_sdag_idx_Timestepping_serial_145();
  (void)_sdag_idx_Timestepping_serial_146();
  (void)_sdag_idx_Timestepping_serial_147();
  (void)_sdag_idx_Timestepping_serial_148();
  (void)_sdag_idx_Timestepping_serial_149();
  (void)_sdag_idx_Timestepping_serial_150();
  (void)_sdag_idx_Timestepping_serial_151();
  (void)_sdag_idx_Timestepping_serial_152();
  (void)_sdag_idx_Timestepping_serial_153();
  (void)_sdag_idx_Timestepping_serial_154();
  (void)_sdag_idx_Timestepping_serial_155();
  (void)_sdag_idx_Timestepping_serial_156();
  (void)_sdag_idx_Timestepping_serial_157();
  (void)_sdag_idx_Timestepping_serial_158();
  (void)_sdag_idx_Timestepping_serial_159();
  (void)_sdag_idx_Timestepping_serial_160();
  (void)_sdag_idx_Timestepping_serial_161();
  (void)_sdag_idx_Timestepping_serial_162();
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::start_6_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::diagnostics_ckio_19_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::east_ghost_20_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::west_ghost_21_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::north_ghost_22_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::south_ghost_23_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::top_ghost_24_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::bottom_ghost_25_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::continue_timestepping_26_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::start_6_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::diagnostics_ckio_19_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::east_ghost_20_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::west_ghost_21_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::north_ghost_22_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::south_ghost_23_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::top_ghost_24_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::bottom_ghost_25_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::continue_timestepping_26_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_idx3srcpt_tosend_27_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs_28_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_odd_gfs_29_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_k_even_gfs_30_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_31_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_auxevol_gfs_32_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_diagnostic_output_gfs_33_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1_34_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2_35_closure));
  PUPable_reg(SINGLE_ARG(Closure_Timestepping::receiv_wavespeed_at_outer_boundary_36_closure));
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_0() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_0();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_0() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_0", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_1() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_1();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_1() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_1", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_2() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_2();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_2() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_2", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_3() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_3();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_3() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_3", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_4() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_4();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_4() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_4", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_5() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_5();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_5() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_5", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_6() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_6();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_6() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_6", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_7() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_7();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_7() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_7", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_8() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_8();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_8() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_8", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_9() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_9();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_9() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_9", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_10() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_10();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_10() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_10", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_11() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_11();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_11() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_11", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_12() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_12();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_12() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_12", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_13() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_13();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_13() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_13", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_14() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_14();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_14() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_14", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_15() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_15();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_15() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_15", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_16() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_16();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_16() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_16", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_17() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_17();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_17() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_17", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_18() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_18();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_18() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_18", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_19() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_19();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_19() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_19", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_20() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_20();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_20() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_20", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_21() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_21();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_21() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_21", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_22() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_22();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_22() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_22", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_23() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_23();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_23() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_23", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_24() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_24();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_24() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_24", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_25() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_25();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_25() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_25", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_27() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_27();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_27() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_27", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_28() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_28();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_28() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_28", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_29() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_29();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_29() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_29", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_30() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_30();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_30() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_30", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_31() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_31();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_31() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_31", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_32() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_32();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_32() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_32", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_33() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_33();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_33() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_33", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_34() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_34();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_34() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_34", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_35() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_35();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_35() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_35", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_36() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_36();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_36() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_36", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_37() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_37();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_37() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_37", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_38() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_38();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_38() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_38", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_39() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_39();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_39() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_39", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_40() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_40();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_40() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_40", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_41() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_41();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_41() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_41", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_42() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_42();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_42() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_42", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_43() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_43();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_43() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_43", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_44() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_44();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_44() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_44", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_26() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_26();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_26() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_26", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_45() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_45();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_45() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_45", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_46() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_46();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_46() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_46", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_47() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_47();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_47() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_47", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_48() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_48();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_48() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_48", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_49() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_49();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_49() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_49", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_50() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_50();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_50() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_50", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_51() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_51();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_51() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_51", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_52() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_52();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_52() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_52", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_53() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_53();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_53() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_53", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_54() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_54();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_54() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_54", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_55() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_55();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_55() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_55", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_56() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_56();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_56() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_56", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_57() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_57();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_57() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_57", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_58() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_58();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_58() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_58", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_59() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_59();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_59() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_59", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_60() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_60();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_60() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_60", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_61() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_61();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_61() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_61", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_62() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_62();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_62() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_62", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_63() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_63();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_63() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_63", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_64() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_64();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_64() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_64", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_65() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_65();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_65() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_65", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_66() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_66();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_66() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_66", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_67() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_67();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_67() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_67", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_68() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_68();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_68() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_68", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_69() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_69();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_69() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_69", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_70() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_70();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_70() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_70", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_71() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_71();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_71() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_71", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_72() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_72();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_72() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_72", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_73() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_73();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_73() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_73", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_74() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_74();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_74() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_74", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_75() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_75();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_75() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_75", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_76() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_76();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_76() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_76", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_77() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_77();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_77() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_77", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_78() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_78();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_78() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_78", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_79() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_79();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_79() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_79", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_80() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_80();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_80() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_80", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_81() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_81();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_81() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_81", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_82() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_82();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_82() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_82", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_83() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_83();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_83() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_83", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_84() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_84();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_84() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_84", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_85() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_85();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_85() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_85", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_86() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_86();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_86() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_86", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_87() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_87();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_87() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_87", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_88() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_88();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_88() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_88", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_89() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_89();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_89() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_89", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_90() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_90();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_90() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_90", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_91() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_91();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_91() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_91", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_92() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_92();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_92() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_92", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_93() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_93();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_93() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_93", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_94() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_94();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_94() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_94", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_95() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_95();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_95() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_95", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_96() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_96();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_96() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_96", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_97() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_97();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_97() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_97", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_98() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_98();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_98() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_98", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_99() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_99();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_99() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_99", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_100() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_100();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_100() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_100", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_101() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_101();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_101() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_101", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_102() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_102();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_102() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_102", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_103() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_103();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_103() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_103", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_104() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_104();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_104() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_104", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_105() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_105();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_105() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_105", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_106() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_106();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_106() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_106", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_107() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_107();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_107() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_107", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_108() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_108();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_108() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_108", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_109() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_109();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_109() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_109", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_110() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_110();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_110() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_110", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_111() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_111();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_111() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_111", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_112() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_112();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_112() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_112", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_113() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_113();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_113() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_113", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_114() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_114();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_114() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_114", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_115() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_115();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_115() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_115", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_116() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_116();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_116() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_116", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_117() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_117();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_117() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_117", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_118() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_118();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_118() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_118", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_119() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_119();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_119() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_119", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_120() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_120();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_120() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_120", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_121() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_121();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_121() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_121", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_122() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_122();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_122() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_122", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_123() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_123();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_123() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_123", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_124() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_124();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_124() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_124", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_125() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_125();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_125() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_125", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_126() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_126();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_126() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_126", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_127() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_127();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_127() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_127", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_128() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_128();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_128() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_128", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_129() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_129();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_129() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_129", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_130() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_130();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_130() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_130", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_131() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_131();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_131() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_131", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_132() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_132();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_132() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_132", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_133() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_133();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_133() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_133", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_134() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_134();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_134() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_134", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_135() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_135();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_135() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_135", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_136() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_136();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_136() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_136", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_137() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_137();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_137() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_137", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_138() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_138();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_138() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_138", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_139() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_139();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_139() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_139", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_140() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_140();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_140() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_140", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_141() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_141();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_141() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_141", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_142() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_142();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_142() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_142", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_143() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_143();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_143() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_143", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_144() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_144();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_144() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_144", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_145() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_145();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_145() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_145", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_146() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_146();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_146() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_146", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_147() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_147();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_147() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_147", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_148() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_148();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_148() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_148", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_149() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_149();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_149() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_149", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_150() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_150();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_150() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_150", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_151() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_151();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_151() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_151", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_152() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_152();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_152() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_152", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_153() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_153();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_153() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_153", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_154() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_154();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_154() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_154", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_155() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_155();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_155() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_155", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_156() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_156();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_156() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_156", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_157() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_157();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_157() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_157", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_158() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_158();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_158() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_158", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_159() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_159();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_159() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_159", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_160() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_160();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_160() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_160", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_161() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_161();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_161() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_161", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_idx_Timestepping_serial_162() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  static int epidx = _sdag_reg_Timestepping_serial_162();
  return epidx;
}
#endif /* CK_TEMPLATES_ONLY */


#ifndef CK_TEMPLATES_ONLY
int Timestepping::_sdag_reg_Timestepping_serial_162() { // Potentially missing Timestepping_SDAG_CODE in your class definition?
  return CkRegisterEp("Timestepping_serial_162", NULL, 0, CkIndex_Timestepping::__idx, 0);
}
#endif /* CK_TEMPLATES_ONLY */



#ifndef CK_TEMPLATES_ONLY
void _registertimestepping(void)
{
  static int _done = 0; if(_done) return; _done = 1;





/* REG: array Timestepping: ArrayElement{
Timestepping(const CommondataObject &inData);
void ready_1d_y(Ck::IO::FileReadyMsg* impl_msg);
void ready_1d_z(Ck::IO::FileReadyMsg* impl_msg);
void ready_2d_xy(Ck::IO::FileReadyMsg* impl_msg);
void ready_2d_yz(Ck::IO::FileReadyMsg* impl_msg);
void start();
void start_write_1d_y(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_1d_z(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_2d_xy(Ck::IO::SessionReadyMsg* impl_msg);
void start_write_2d_yz(Ck::IO::SessionReadyMsg* impl_msg);
void test_written_1d_y(CkReductionMsg* impl_msg);
void test_written_1d_z(CkReductionMsg* impl_msg);
void test_written_2d_xy(CkReductionMsg* impl_msg);
void test_written_2d_yz(CkReductionMsg* impl_msg);
void closed_1d_y(CkReductionMsg* impl_msg);
void closed_1d_z(CkReductionMsg* impl_msg);
void closed_2d_xy(CkReductionMsg* impl_msg);
void closed_2d_yz(CkReductionMsg* impl_msg);
void diagnostics_ckio(const Ck::IO::Session &token, int which_diagnostics_part);
void east_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void west_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void north_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void south_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void top_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void bottom_ghost(int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void continue_timestepping();
void receiv_nonlocalinnerbc_idx3srcpt_tosend(int idx3_of_sendingchare, int num_srcpts, const int *globalidx3_srcpts);
void receiv_nonlocalinnerbc_data_y_nplus1_running_total_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_k_odd_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_k_even_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_auxevol_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_diagnostic_output_gfs(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part1(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_nonlocalinnerbc_data_y_n_gfs_initialdata_part2(int src_chare_idx3, int type_gfs, int len_tmpBuffer, const REAL *tmpBuffer);
void receiv_wavespeed_at_outer_boundary(const REAL &wavespeed_at_outer_boundary);
Timestepping(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_Timestepping::__register("Timestepping", sizeof(Timestepping));

}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_Timestepping::virtual_pup(PUP::er &p) {
    recursive_pup<Timestepping>(dynamic_cast<Timestepping*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
