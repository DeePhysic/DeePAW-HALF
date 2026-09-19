//*****************************************************************************
//  Copyright (C) 2023 Advanced Micro Devices, Inc. All rights reserved.
//*****************************************************************************

#include <rocfft.h>
#include <hip_runtime.h>
#include <iostream>
#include <string.h>

#define ROCFFT_ERR(err) \
    if(err != rocfft_status_success) \
    { \
        switch (err) { \
            case rocfft_status_failure: \
                std::cerr << "rocfft_status_failure" << std::endl; \
                break; \
            case rocfft_status_invalid_arg_value: \
                std::cerr << "rocfft_status_invalid_arg_value" << std::endl; \
                break; \
            case rocfft_status_invalid_dimensions : \
                std::cerr << "rocfft_status_invalid_dimensions" << std::endl; \
                break; \
            case rocfft_status_invalid_array_type : \
                std::cerr << "rocfft_status_invalid_array_type" << std::endl; \
                break; \
            case rocfft_status_invalid_strides : \
                std::cerr << "rocfft_status_invalid_strides" << std::endl; \
                break; \
            case rocfft_status_invalid_distance : \
                std::cerr << "rocfft_status_invalid_distance" << std::endl; \
                break; \
            case rocfft_status_invalid_offset : \
                std::cerr << "rocfft_status_invalid_offset" << std::endl; \
                break; \
            case rocfft_status_invalid_work_buffer: \
                std::cerr << "rocfft_status_invalid_work_buffer" << std::endl; \
                break; \
            default: \
                std::cerr << "Unknown rocfft_status" << std::endl; \
                break; \
        } \
        exit(-1); \
    }
extern "C" {
    void hip_init_rocfft() { rocfft_setup(); }
    void hip_finish_rocfft() { rocfft_cleanup(); }
    // ROCFFT
    void hip_rocfft_create_omp(int ttype, int tdim, int tnum, int *tsizes,
            int *instride, int idist, int *outstride, int odist,
            rocfft_plan *planf, rocfft_plan *plani,
            rocfft_plan_description *descf, rocfft_plan_description *desci,
            rocfft_execution_info *infof, rocfft_execution_info *infoi,
            void *&work) {
        size_t size_f,size_i,max_size = 0;
        const size_t N[3] = {(size_t)tsizes[0],(size_t)tsizes[1],(size_t)tsizes[2]};
        const size_t istride[3] = {(size_t)instride[0],(size_t)instride[1],(size_t)instride[2]};
        const size_t ostride[3] = {(size_t)outstride[0],(size_t)outstride[1],(size_t)outstride[2]};
        rocfft_result_placement     placement;
        rocfft_transform_type       f_trans, i_trans;
        rocfft_array_type           f_in, f_out, i_in, i_out;
        if (ttype) { // Real types
            // R2C
            placement = rocfft_placement_inplace;
            // Forward
            f_trans = rocfft_transform_type_real_forward;
            f_in = rocfft_array_type_real;
            f_out = rocfft_array_type_hermitian_interleaved;
            // Inverse
            i_trans = rocfft_transform_type_real_inverse;
            i_in = rocfft_array_type_hermitian_interleaved;
            i_out = rocfft_array_type_real;
        } else { // Complex types
            // C2C
            placement = rocfft_placement_inplace;
            // Forward
            f_trans = rocfft_transform_type_complex_forward;
            f_in = rocfft_array_type_complex_interleaved;
            f_out = rocfft_array_type_complex_interleaved;
            // Inverse
            i_trans = rocfft_transform_type_complex_inverse;
            i_in = rocfft_array_type_complex_interleaved;
            i_out = rocfft_array_type_complex_interleaved;
        }
        // Initialize description and info
        ROCFFT_ERR(rocfft_plan_description_create(descf));
        ROCFFT_ERR(rocfft_plan_description_create(desci));
        ROCFFT_ERR(rocfft_execution_info_create(infof));
        ROCFFT_ERR(rocfft_execution_info_create(infoi));
        // Forward plans
        ROCFFT_ERR(rocfft_plan_description_set_data_layout(*descf,
            f_in,f_out,0,0,tdim,istride,idist,tdim,ostride,odist));
        ROCFFT_ERR(rocfft_plan_create(planf, placement, f_trans, rocfft_precision_double,
            tdim, N, tnum, *descf));
        // Inverse plans
        ROCFFT_ERR(rocfft_plan_description_set_data_layout(*desci,
            i_in,i_out,0,0,tdim,ostride,odist,tdim,istride,idist));
        ROCFFT_ERR(rocfft_plan_create(plani, placement, i_trans, rocfft_precision_double,
            tdim, N, tnum, *desci));
        // Optional work buffers
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*planf, &size_f));
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*plani, &size_i));
        max_size = std::max(size_f,size_i);
        if(size_f || size_i)
            assert(hipMalloc(&work, max_size) == hipSuccess);
        if(size_f)
            ROCFFT_ERR(rocfft_execution_info_set_work_buffer(*infof, work, max_size));
        if(size_i)
            ROCFFT_ERR(rocfft_execution_info_set_work_buffer(*infoi, work, max_size));
    }
    void hip_rocfft_create_acc(int ttype, int tdim, int tnum, int *tsizes,
            int *instride, int idist, int *outstride, int odist,
            rocfft_plan *planf, rocfft_plan *plani,
            rocfft_plan_description *descf, rocfft_plan_description *desci,
            rocfft_execution_info *infof, rocfft_execution_info *infoi,
            void *&work, void *stream) {
        hipStream_t *hipstream = (hipStream_t *) stream;
        size_t size_f,size_i,max_size = 0;
        const size_t N[3] = {(size_t)tsizes[0],(size_t)tsizes[1],(size_t)tsizes[2]};
        const size_t istride[3] = {(size_t)instride[0],(size_t)instride[1],(size_t)instride[2]};
        const size_t ostride[3] = {(size_t)outstride[0],(size_t)outstride[1],(size_t)outstride[2]};
        rocfft_result_placement     placement;
        rocfft_transform_type       f_trans, i_trans;
        rocfft_array_type           f_in, f_out, i_in, i_out;
        if (ttype) { // Real types
            // R2C
            placement = rocfft_placement_inplace;
            // Forward
            f_trans = rocfft_transform_type_real_forward;
            f_in = rocfft_array_type_real;
            f_out = rocfft_array_type_hermitian_interleaved;
            // Inverse
            i_trans = rocfft_transform_type_real_inverse;
            i_in = rocfft_array_type_hermitian_interleaved;
            i_out = rocfft_array_type_real;
        } else { // Complex types
            // C2C
            placement = rocfft_placement_inplace;
            // Forward
            f_trans = rocfft_transform_type_complex_forward;
            f_in = rocfft_array_type_complex_interleaved;
            f_out = rocfft_array_type_complex_interleaved;
            // Inverse
            i_trans = rocfft_transform_type_complex_inverse;
            i_in = rocfft_array_type_complex_interleaved;
            i_out = rocfft_array_type_complex_interleaved;
        }
        // Initialize description and info
        ROCFFT_ERR(rocfft_plan_description_create(descf));
        ROCFFT_ERR(rocfft_plan_description_create(desci));
        ROCFFT_ERR(rocfft_execution_info_create(infof));
        ROCFFT_ERR(rocfft_execution_info_create(infoi));
        // Forward plans
        ROCFFT_ERR(rocfft_plan_description_set_data_layout(*descf,
            f_in,f_out,0,0,tdim,istride,idist,tdim,ostride,odist));
        ROCFFT_ERR(rocfft_plan_create(planf, placement, f_trans, rocfft_precision_double,
            tdim, N, tnum, *descf));
        // Inverse plans
        ROCFFT_ERR(rocfft_plan_description_set_data_layout(*desci,
            i_in,i_out,0,0,tdim,ostride,odist,tdim,istride,idist));
        ROCFFT_ERR(rocfft_plan_create(plani, placement, i_trans, rocfft_precision_double,
            tdim, N, tnum, *desci));
        // Set streams
        ROCFFT_ERR(rocfft_execution_info_set_stream(*infof,*hipstream));
        ROCFFT_ERR(rocfft_execution_info_set_stream(*infoi,*hipstream));
        // Optional work buffers
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*planf, &size_f));
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*plani, &size_i));
        max_size = std::max(size_f,size_i);
        if(size_f || size_i)
            assert(hipMalloc(&work, max_size) == hipSuccess);
        if(size_f)
            ROCFFT_ERR(rocfft_execution_info_set_work_buffer(*infof, work, max_size));
        if(size_i)
            ROCFFT_ERR(rocfft_execution_info_set_work_buffer(*infoi, work, max_size));
    }
    void hip_rocfft_exec(void *ptr, void *in, void *out, void *ctx) {
        rocfft_plan *plan = (rocfft_plan *) ptr;
        rocfft_execution_info *info = (rocfft_execution_info *) ctx;
        ROCFFT_ERR(rocfft_execute(*plan, (void **) &in, (void **) &out, *info));
    }
    void hip_rocfft_destroy(rocfft_plan *planf, rocfft_plan *plani,
            rocfft_plan_description *descf, rocfft_plan_description *desci,
            rocfft_execution_info *infof, rocfft_execution_info *infoi,
            void *&work) {
        size_t size_f,size_i;
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*planf, &size_f));
        ROCFFT_ERR(rocfft_plan_get_work_buffer_size(*plani, &size_i));
        if(size_f || size_i ) assert(hipFree(work) == hipSuccess);
        ROCFFT_ERR(rocfft_plan_description_destroy(*descf));
        ROCFFT_ERR(rocfft_plan_description_destroy(*desci));
        ROCFFT_ERR(rocfft_execution_info_destroy(*infof));
        ROCFFT_ERR(rocfft_execution_info_destroy(*infoi));
        ROCFFT_ERR(rocfft_plan_destroy(*planf));
        ROCFFT_ERR(rocfft_plan_destroy(*plani));
    }
}
