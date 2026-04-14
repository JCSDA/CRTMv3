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

// Input/Output names (Must match the model.onnx)
const char* input_names[] = {"input"};
const char* output_names[] = {"output"};

/**
 * Initialize the ONNX Runtime and load the model.
 * Returns 0 on success, 1 on failure.
 */
int crtm_onnx_init(const char* model_path) {
    g_ort = OrtGetApiBase()->GetApi(ORT_API_VERSION);
    if (!g_ort) return 1;

    // Create environment
    g_ort->CreateEnv(ORT_LOGGING_LEVEL_WARNING, "crtm_onnx", &g_env);
    
    // Create session options
    g_ort->CreateSessionOptions(&g_session_options);
    g_ort->SetIntraOpNumThreads(g_session_options, 1);
    g_ort->SetSessionGraphOptimizationLevel(g_session_options, ORT_ENABLE_ALL);

    // Load and create session
    // Note: On Windows, this would need a wide string path
    OrtStatus* status = g_ort->CreateSession(g_env, model_path, g_session_options, &g_session);
    if (status != NULL) {
        const char* msg = g_ort->GetErrorMessage(status);
        fprintf(stderr, "Error creating ONNX session: %s\n", msg);
        g_ort->ReleaseStatus(status);
        return 1;
    }

    printf("CRTM ONNX Bridge: Successfully loaded model %s\n", model_path);
    return 0;
}

/**
 * Perform inference.
 * input_data: [input_dim]
 * output_data: [output_dim]
 */
int crtm_onnx_predict(const float* input_data, size_t input_dim, 
                       float* output_data, size_t output_dim) {
    if (!g_session) return 1;

    OrtMemoryInfo* memory_info;
    g_ort->CreateCpuMemoryInfo(OrtDeviceAllocator, OrtMemTypeDefault, &memory_info);

    int64_t input_shape[] = {1, (int64_t)input_dim};
    OrtValue* input_tensor = NULL;
    g_ort->CreateTensorWithDataAsOrtValue(memory_info, (void*)input_data, 
                                          input_dim * sizeof(float), input_shape, 2, 
                                          ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &input_tensor);
    g_ort->ReleaseMemoryInfo(memory_info);

    OrtValue* output_tensor = NULL;
    g_ort->Run(g_session, NULL, input_names, (const OrtValue* const*)&input_tensor, 1, 
               output_names, 1, &output_tensor);

    float* out_ptr;
    g_ort->GetTensorMutableData(output_tensor, (void**)&out_ptr);
    
    // Copy results to output buffer
    for (size_t i = 0; i < output_dim; ++i) {
        output_data[i] = out_ptr[i];
    }

    g_ort->ReleaseValue(input_tensor);
    g_ort->ReleaseValue(output_tensor);
    return 0;
}

/**
 * Cleanup
 */
void crtm_onnx_cleanup() {
    if (g_session) g_ort->ReleaseSession(g_session);
    if (g_session_options) g_ort->ReleaseSessionOptions(g_session_options);
    if (g_env) g_ort->ReleaseEnv(g_env);
}
