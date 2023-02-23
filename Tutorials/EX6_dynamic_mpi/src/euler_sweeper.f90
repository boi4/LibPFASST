! TODO
module euler_sweeper
  use pf_mod_dtype

  implicit none

  type, extends(pf_sweeper_t), abstract :: euler_sweeper_t
     integer     :: npieces
     logical     :: use_LUq
  contains
     procedure :: sweep        => euler_sweep
     procedure :: initialize   => euler_initialize
     procedure :: evaluate     => euler_evaluate
     procedure :: integrate    => euler_integrate
     procedure :: evaluate_all => euler_evaluate_all
     procedure :: residual     => euler_residual
     procedure :: spreadq0     => euler_spreadq0
     procedure :: compute_dt   => euler_compute_dt
     procedure :: destroy      => euler_destroy
  end type euler_sweeper_t

end module euler_sweeper
