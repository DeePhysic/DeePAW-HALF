module half_fft
  use iso_c_binding, only: c_ptr,c_int,c_double_complex,c_associated
  use half_kinds, only: dp,i32
  implicit none
  private
  public :: fft3_forward,fft3_backward
  integer(c_int),parameter::fftw_forward=-1,fftw_backward=1,fftw_estimate=64
  interface
    function fftw_plan_dft_3d(n0,n1,n2,input,output,sign,flags)bind(C,name='fftw_plan_dft_3d')result(plan)
      import c_ptr,c_int
      integer(c_int),value::n0,n1,n2,sign,flags
      type(c_ptr),value::input,output
      type(c_ptr)::plan
    end function
    subroutine fftw_execute_dft(plan,input,output)bind(C,name='fftw_execute_dft')
      import c_ptr
      type(c_ptr),value::plan,input,output
    end subroutine
    subroutine fftw_destroy_plan(plan)bind(C,name='fftw_destroy_plan')
      import c_ptr
      type(c_ptr),value::plan
    end subroutine
  end interface
contains
  subroutine fft3_forward(shape,input,output)
    integer(i32),intent(in)::shape(3)
    complex(dp),target,intent(inout)::input(:)
    complex(dp),target,intent(out)::output(:)
    call transform(shape,input,output,fftw_forward)
  end subroutine
  subroutine fft3_backward(shape,input,output)
    integer(i32),intent(in)::shape(3)
    complex(dp),target,intent(inout)::input(:)
    complex(dp),target,intent(out)::output(:)
    call transform(shape,input,output,fftw_backward)
  end subroutine
  subroutine transform(shape,input,output,sign)
    use iso_c_binding,only:c_loc
    integer(i32),intent(in)::shape(3)
    complex(dp),target,intent(inout)::input(:)
    complex(dp),target,intent(out)::output(:)
    integer(c_int),intent(in)::sign
    type(c_ptr)::plan
    if(size(input)/=product(shape).or.size(output)/=product(shape))error stop 'HALF: FFT size mismatch'
    plan=fftw_plan_dft_3d(shape(1),shape(2),shape(3),c_loc(input),c_loc(output),sign,fftw_estimate)
    if(.not.c_associated(plan))error stop 'HALF: FFT plan creation failed'
    call fftw_execute_dft(plan,c_loc(input),c_loc(output)); call fftw_destroy_plan(plan)
  end subroutine
end module half_fft
