program particle_types_test_driver
use fruit
use mod_particle_types_test, only: run_fruit_particle_types
  implicit none

  ! init fruit suite
  call init_fruit_xml
  call init_fruit
  call fruit_hide_dots
  ! run particle_types test basket
  call run_fruit_particle_types

  ! write test summary and finalize test suit
  call fruit_summary_xml
  call fruit_finalize
end program particle_types_test_driver
