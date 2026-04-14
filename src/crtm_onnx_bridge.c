#include <stdio.h>
#include <stdlib.h>
#include <onnxruntime_c_api.h>

/*
 * CRTM ONNX Bridge
 * Provides a C interface for Fortran to call ONNX models.
 */

const OrtApi* g_ort = NULL;
OrtEnv* g_env = NULL;
OrtSession* g_session = NULL;
OrtSessionOptions* g_session_options = NULL;

const char* input_names[] = {"input"};
const char* output_names[] = {"output"};

int crtm_onnx_init(const char* model_path) {
    g_ort = OrtGetApiBase()->GetApi(ORT_API_VERSION);
    if (!g_ort) return 1;
    g_ort->CreateEnv(ORT_LOGGING_LEVEL_WARNING, "crtm_onnx", &g_env);
    g_ort->CreateSessionOptions(&g_session_options);
    g_ort->SetIntraOpNumThreads(g_session_options, 1);
    g_ort->SetSessionGraphOptimizationLevel(g_session_options, ORT_ENABLE_ALL);
    OrtStatus* status = g_ort->CreateSession(g_env, model_path, g_session_options, &g_session);
    if (status != NULL) {
        g_ort->ReleaseStatus(status);
        return 1;
    }
    return 0;
}

/**
 * Perform batch inference.
 * input_data: [batch_size, input_dim]
 * output_data: [batch_size, output_dim]
 */
int crtm_onnx_predict(const float* input_data, size_t batch_size, size_t input_dim, 
                       float* output_data, size_t output_dim) {
    if (!g_session) return 1;

    OrtMemoryInfo* memory_info;
    g_ort->CreateCpuMemoryInfo(OrtDeviceAllocator, OrtMemTypeDefault, &memory_info);

    int64_t input_shape[] = {(int64_t)batch_size, (int64_t)input_dim};
    OrtValue* input_tensor = NULL;
    g_ort->CreateTensorWithDataAsOrtValue(memory_info, (void*)input_data, 
                                          batch_size * input_dim * sizeof(float), 
                                          input_shape, 2, 
                                          ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &input_tensor);
    g_ort->ReleaseMemoryInfo(memory_info);

    OrtValue* output_tensor = NULL;
    g_ort->Run(g_session, NULL, input_names, (const OrtValue* const*)&input_tensor, 1, 
               output_names, 1, &output_tensor);

    float* out_ptr;
    g_ort->GetTensorMutableData(output_tensor, (void**)&out_ptr);
    
    // Copy batch results
    for (size_t i = 0; i < batch_size * output_dim; ++i) {
        output_data[i] = out_ptr[i];
    }

    g_ort->ReleaseValue(input_tensor);
    g_ort->ReleaseValue(output_tensor);
    return 0;
}

void crtm_onnx_cleanup() {
    if (g_session) g_ort->ReleaseSession(g_session);
    if (g_session_options) g_ort->ReleaseSessionOptions(g_session_options);
    if (g_env) g_ort->ReleaseEnv(g_env);
}
