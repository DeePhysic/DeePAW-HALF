! currently the default is to use double precision (flag double_prec_bse);
! the flag single_prec_bse allows to select single precision storage of the BSE
! matrix and thus saves a factor 2 in storage (but usually also in CPU time!!);
! important: use of single precision is absolutely safe and highly recommended!!
! if you want to use single precision define precompiler flag single_prec_bse !
#define double_prec_bse
#if defined(single_prec_bse)
#undef double_prec_bse
#endif

! the flag USE_ZHEEVX determines whether ZHEEVX is used; one has to know that
! ZHEEVX is about twice to three times faster than ZHEEV(D), but also requires
! twice as much memory -- that's the downside of the coin and maybe a problem
!
! alternatively, one can also use ZHEEVD by defining precompiler flag USE_ZHEEVD
! which uses "divide-and-conquer" algorithms for diagonalization and could be
! a bit faster than ZHEEV for larger matrices (if eigenvectors are computed);
! currently USE_ZHEEVX is the default if no different precompiler flag is given
!
! for those which like to use/test ZHEEV another (undocument) flag BSE_ZHEEV
! may be set to cancel the definition of USE_ZHEEVX (and/or USE_ZHEEVD ...);
! this only affects the BSE routines (and no use of USE_ZHEEVX/D elsewhere ...);

#if defined(single_prec_bse ) && defined(gammareal)
#define GGEMMbse SGEMM
#define GETRIbse SGETRI
#define GETRFbse SGETRF
#define onebse   1._qs
#define zerobse  0._qs
#define halfbse  0.5_qs
#define GPGEMMbse PSGEMM
#define GPTRANbse PSTRAN
#define PGGEMR2Dbse  WPSGEMR2D
#define PGTRMR2Dbse  WPSTRMR2D
#define PGTRADDbse WPSTRADD
#define PGGEADDbse WPSGEADD
#define GCOPYbse   WDCOPY
#define PGLASETbse PSLASET
#define GGEVX WSSYEVX
#define PGLAPV2bse PSLAPV2
#define BGEMM WSGEMM
#define DOTbse SDOT
#endif

#if defined(single_prec_bse ) && ! defined(gammareal)
#define GGEMMbse WCGEMM
#define GETRIbse WCGETRI
#define GETRFbse WCGETRF
#define onebse   (1._qs,0._qs)
#define zerobse  (0._qs,0._qs)
#define halfbse  (0.5_qs,0._qs)
#define GPGEMMbse PCGEMM
#define GPTRANbse PCTRANU
#define PGGEMR2Dbse  WPCGEMR2D
#define PGTRMR2Dbse  WPCTRMR2D
#define PGTRADDbse WPCTRADD
#define PGGEADDbse WPCGEADD
#define GCOPYbse   WZCOPY
#define PGLASETbse PCLASET
#define GGEVX WCHEEVX
#define PGLAPV2bse PCLAPV2
#define BGEMM WCGEMM
#define DOTbse CDOTC
#endif

#if defined(double_prec_bse ) && defined(gammareal)
#define GGEMMbse WDGEMM
#define GETRIbse DGETRI
#define GETRFbse DGETRF
#define onebse   1._q
#define zerobse  0._q
#define halfbse  0.5_q
#define GPGEMMbse PDGEMM
#define GPTRANbse PDTRAN
#define PGGEMR2Dbse  WPDGEMR2D
#define PGTRMR2Dbse  WPDTRMR2D
#define PGTRADDbse WPDTRADD
#define PGGEADDbse WPDGEADD
#define GCOPYbse   WDCOPY
#define PGLASETbse PDLASET
#define GGEVX WDSYEVX
#define PGLAPV2bse PDLAPV2
#define BGEMM WDGEMM
#define DOTbse DDOT
#endif

#if defined(double_prec_bse ) && ! defined(gammareal)
#define GGEMMbse WZGEMM
#define GETRIbse ZGETRI
#define GETRFbse ZGETRF
#define onebse   (1._q,0._q)
#define zerobse  (0._q,0._q)
#define halfbse  (0.5_q,0._q)
#define GPGEMMbse PZGEMM
#define GPTRANbse PZTRANU
#define PGGEMR2Dbse  WPZGEMR2D
#define PGTRMR2Dbse  WPZTRMR2D
#define PGTRADDbse WPZTRADD
#define PGGEADDbse WPZGEADD
#define GCOPYbse   WZCOPY
#define PGLASETbse PZLASET
#define GGEVX WZHEEVX
#define PGLAPV2bse PZLAPV2
#define BGEMM WZGEMM
#define DOTbse ZDOTC
#endif

#ifdef double_prec_bse
#define PGGEMR2Dtebse WPZGEMR2D
#define M_sumbsete_v M_sum_z
#define M_bcastbsete_g_from M_bcast_z_from
#define BDOTC WZDOTC

#if defined(gammareal)
#define M_bcastbse_g_from M_bcast_d_from
#else
#define M_bcastbse_g_from M_bcast_z_from
#endif

#ifdef USE_ZHEEVD
#define PGGEVX PDSYEV_ZHEEVD_DESC
#else // not USE_ZHEEVD
#define PGGEVX PDSYEVX_ZHEEVX_DESC
#endif

#else // not double_prec_bse
#define PGGEMR2Dtebse WPCGEMR2D
#define M_sumbsete_v M_sum_c
#define M_bcastbsete_g_from M_bcast_c_from
#define BDOTC WCDOTC

#if defined(gammareal)
#define M_bcastbse_g_from M_bcast_s_from
#else
#define M_bcastbse_g_from M_bcast_c_from
#endif

#ifdef USE_ZHEEVD
#define PGGEVX PSSYEV_CHEEVD_DESC
#else // not USE_ZHEEVD
#define PGGEVX PSSYEVX_CHEEVX_DESC
#endif

#endif


! Debugging
!#define bsedebug

#ifdef gammareal

#define M_sum_g8 M_sum_d8
#if single_prec_bse
#define M_sumbse_v M_sum_single
#else
#define M_sumbse_v M_sum_d
#endif

#else // gammareal

#define M_sum_g8 M_sum_z8
#if single_prec_bse
#define M_sumbse_v M_sum_c
#else
#define M_sumbse_v M_sum_z
#endif

#endif // gammareal
