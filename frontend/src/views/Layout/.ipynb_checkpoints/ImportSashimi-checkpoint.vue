<template>
    <div class="sashimi-config-container">
      <!-- Parameters Section -->
      <div class="params-section">
        <el-form :model="params" label-width="200px">
          
          <!-- Required Parameters -->
          <div class="required-params">
            <h4 class="section-title required">Required Parameters</h4>

            <!-- BAM Files (User Input Required) -->
            <div class="param-group">
              <div class="param-group-title">BAM Files (Required)</div>
              
              <el-form-item label="Sample 1 BAM File (--b1)">
                <el-input 
                  v-model="params.b1" 
                  placeholder="/path/to/sample1.bam"
                  @blur="validateBamPath(params.b1, 'b1')"
                >
                  <template #suffix>
                    <el-icon v-if="getValidationStatus('b1') === 'success'" style="color: #67c23a">
                      <CircleCheck />
                    </el-icon>
                    <el-icon v-else-if="getValidationStatus('b1') === 'error'" style="color: #f56c6c">
                      <CircleClose />
                    </el-icon>
                  </template>
                </el-input>
                <div v-if="validationResults.b1" class="validation-result">
                  <el-alert
                    v-if="validationResults.b1.valid"
                    type="success"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('b1')"
                  />
                  <el-alert
                    v-else
                    type="error"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('b1')"
                  />
                </div>
              </el-form-item>
              
              <el-form-item label="Sample 2 BAM File (--b2)">
                <el-input 
                  v-model="params.b2" 
                  placeholder="/path/to/sample2.bam"
                  @blur="validateBamPath(params.b2, 'b2')"
                >
                  <template #suffix>
                    <el-icon v-if="getValidationStatus('b2') === 'success'" style="color: #67c23a">
                      <CircleCheck />
                    </el-icon>
                    <el-icon v-else-if="getValidationStatus('b2') === 'error'" style="color: #f56c6c">
                      <CircleClose />
                    </el-icon>
                  </template>
                </el-input>
                <div v-if="validationResults.b2" class="validation-result">
                  <el-alert
                    v-if="validationResults.b2.valid"
                    type="success"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('b2')"
                  />
                  <el-alert
                    v-else
                    type="error"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('b2')"
                  />
                </div>
              </el-form-item>
            </div>

            <!-- GFF File -->
            <div class="param-group">
              <div class="param-group-title">Annotation File</div>
              <el-form-item label="GFF Annotation File">
                <el-input 
                  v-model="params.gffFile" 
                  placeholder="/path/to/annotation.gff3"
                  @blur="validateGffFile"
                >
                  <template #suffix>
                    <el-icon v-if="getValidationStatus('gffFile') === 'success'" style="color: #67c23a">
                      <CircleCheck />
                    </el-icon>
                    <el-icon v-else-if="getValidationStatus('gffFile') === 'error'" style="color: #f56c6c">
                      <CircleClose />
                    </el-icon>
                  </template>
                </el-input>
                <div v-if="validationResults.gffFile" class="validation-result">
                  <el-alert
                    v-if="validationResults.gffFile.valid"
                    type="success"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('gffFile')"
                  >
                    <template #default>
                      <div v-if="validationResults.gffFile.size_mb" class="gtf-details">
                        <p><strong>File size:</strong> {{ validationResults.gffFile.size_mb }} MB</p>
                        <p v-if="validationResults.gffFile.is_compressed !== undefined">
                          <strong>Compressed:</strong> {{ validationResults.gffFile.is_compressed ? 'Yes' : 'No' }}
                        </p>
                        <div v-if="validationResults.gffFile.sample_lines?.length > 0">
                          <p><strong>Sample content:</strong></p>
                          <div class="sample-lines">
                            <code v-for="(line, index) in validationResults.gffFile.sample_lines.slice(0, 2)" :key="index">
                              {{ line }}
                            </code>
                          </div>
                        </div>
                      </div>
                    </template>
                  </el-alert>
                  <el-alert
                    v-else
                    type="error"
                    :closable="false"
                    show-icon
                    :title="getValidationMessage('gffFile')"
                  />
                </div>
              </el-form-item>
            </div>
          </div>

          <!-- Basic Options -->
          <div class="basic-options">
            <h4 class="section-title basic">Basic Display Options</h4>
            
            <div class="param-group">
              <div class="param-group-title">Sample Labels</div>
              <el-form-item label="Sample 1 Label (--l1)">
                <el-input 
                  v-model="params.l1" 
                  :placeholder="defaultSample1Label"
                />
              </el-form-item>
              <el-form-item label="Sample 2 Label (--l2)">
                <el-input 
                  v-model="params.l2" 
                  :placeholder="defaultSample2Label"
                />
              </el-form-item>
            </div>

            <div class="param-group">
              <div class="param-group-title">Figure Settings</div>
              <div class="inline-params">
                <el-form-item label="Width">
                  <el-input-number v-model="params.figWidth" :min="4" :max="20" />
                </el-form-item>
                <el-form-item label="Height">
                  <el-input-number v-model="params.figHeight" :min="4" :max="20" />
                </el-form-item>
                <el-form-item label="Min Reads">
                  <el-input-number v-model="params.minCounts" :min="0" />
                </el-form-item>
              </div>
            </div>
          </div>

          <!-- Advanced Options (Collapsible) -->
          <el-collapse>
            <el-collapse-item title="⚙️ Advanced Options" name="advanced">
              <div class="advanced-options">
                <div class="param-group">
                  <div class="param-group-title">Scale Settings</div>
                  <div class="inline-params">
                    <el-form-item label="Exon Scale">
                      <el-input-number v-model="params.exonS" :min="0.1" :step="0.1" />
                    </el-form-item>
                    <el-form-item label="Intron Scale">
                      <el-input-number v-model="params.intronS" :min="0.1" :step="0.1" />
                    </el-form-item>
                  </div>
                </div>

                <div class="param-group">
                  <div class="param-group-title">Display Settings</div>
                  <el-form-item label="Colors">
                    <el-input 
                      v-model="params.color" 
                      placeholder="#CC0011,#FF8800"
                    />
                  </el-form-item>
                  <el-form-item label="Font Size">
                    <el-input-number v-model="params.fontSize" :min="6" :max="24" />
                  </el-form-item>
                </div>

                <div class="param-group">
                  <div class="param-group-title">Processing Options</div>
                  <div class="checkbox-options">
                    <el-checkbox v-model="params.hideNumber">
                      Hide Junction Numbers
                    </el-checkbox>
                    <el-checkbox v-model="params.noTextBackground">
                      No Text Background
                    </el-checkbox>
                    <el-checkbox v-model="params.removeChrPrefix">
                      Remove Chr Prefix (Recommended)
                    </el-checkbox>
                  </div>
                </div>
              </div>
            </el-collapse-item>
          </el-collapse>
        </el-form>
      </div>
    </div>
</template>

<script setup>
import { ref, computed, watch, onMounted } from 'vue'
import { ElMessage, ElMessageBox } from 'element-plus'
import { CircleCheck, CircleClose, Check, Close } from '@element-plus/icons-vue'
import { useAppStore } from '@/stores/app'
import axios from 'axios'

// Store
const store = useAppStore()

// Props
const props = defineProps({
  modelValue: {
    type: Boolean,
    default: false
  },
  outputDirectory: {
    type: String,
    default: ''
  },
  rmatsResultsPath: {
    type: String,
    default: ''
  },
  analysisType: {
    type: String,
    default: ''
  },
  isoformData: {
    type: Object,
    default: () => ({})
  }
})

// Emits
const emit = defineEmits(['update:modelValue', 'sashimi-generated'])

// Reactive Data
const isGenerating = ref(false)
const generatedCommand = ref('')
const statusMessage = ref('')
const statusMessageClass = ref('')

const params = ref({
  b1: '',
  b2: '',
  gffFile: '',
  l1: '',
  l2: '',
  figWidth: 8,
  figHeight: 7,
  minCounts: 10,
  exonS: 1,
  intronS: 1,
  color: '',
  fontSize: 12,
  hideNumber: false,
  noTextBackground: false,
  removeChrPrefix: true
})

const validationResults = ref({
  b1: null,
  b2: null,
  gffFile: null
})

// Computed Properties
const visible = computed({
  get() {
    return props.modelValue
  },
  set(value) {
    emit('update:modelValue', value)
  }
})

// 自动生成坐标字符串
const autoCoordinate = computed(() => {
  if (!props.isoformData) return ''
  
  const { chromosome, strand, exonStart, exonEnd } = props.isoformData
  const gffFile = params.value.gffFile || store.getFilePaths?.sashimiGff || '/path/to/annotation.gff3'
  
  return `${chromosome}:${strand}:${exonStart}:${exonEnd}:${gffFile}`
})

// 默认样本标签
const defaultSample1Label = computed(() => {
  return props.isoformData?.geneSymbol ? `${props.isoformData.geneSymbol}_Sample1` : 'Sample1'
})

const defaultSample2Label = computed(() => {
  return props.isoformData?.geneSymbol ? `${props.isoformData.geneSymbol}_Sample2` : 'Sample2'
})

// 表单验证
const isFormValid = computed(() => {
  // 必需字段检查
  if (!params.value.b1 || !params.value.b2) return false
  if (!params.value.gffFile) return false
  
  // 验证结果检查
  if (validationResults.value.b1 && !validationResults.value.b1.valid) return false
  if (validationResults.value.b2 && !validationResults.value.b2.valid) return false
  if (validationResults.value.gffFile && !validationResults.value.gffFile.valid) return false
  
  return true
})

const validatePath = async (path, type, resultKey, fieldName) => {
  if (!path || !path.trim()) {
    validationResults.value[resultKey] = null
    return
  }
  
  try {
    console.log(`🔍 验证${fieldName}: ${path}`)
    
    const response = await axios.post('/api/validate-path', {
      path: path.trim(),
      type: type
    })
    
    validationResults.value[resultKey] = response.data
    console.log(`✅ ${fieldName}验证结果:`, response.data)
    
    // 如果验证失败，显示弹窗错误
    if (!response.data.valid) {
      // 根据文件类型生成统一的"文件不存在"错误信息
      const fileTypeMap = {
        'bam_files': 'BAM file',
        'gtf_file': 'GFF file',
        'gff_file': 'GFF file'
      }
      const fileTypeName = fileTypeMap[type] || 'File'
      const errorMessage = `${fileTypeName} does not exist: ${path.trim()}`
      
      await ElMessageBox.alert(
        errorMessage,
        `${fieldName} Validation Failed`,
        {
          type: 'error',
          confirmButtonText: 'OK'
        }
      )
    } else {
      // 验证成功，显示简短的成功提示
      ElMessage.success({
        message: `${fieldName} validation successful`,
        duration: 2000
      })
    }
    
  } catch (error) {
    console.error(`❌ ${fieldName}验证出错:`, error)
    
    // 无论什么错误，都显示统一的"文件不存在"错误信息
    const fileTypeMap = {
      'bam_files': 'BAM file',
      'gtf_file': 'GFF file', 
      'gff_file': 'GFF file'
    }
    const fileTypeName = fileTypeMap[type] || 'File'
    const errorMessage = `${fileTypeName} does not exist: ${path.trim()}`
    
    validationResults.value[resultKey] = {
      valid: false,
      error: errorMessage,
      exists: false
    }
    
    // 显示错误弹窗
    await ElMessageBox.alert(
      errorMessage,
      `${fieldName} Validation Failed`,
      {
        type: 'error',
        confirmButtonText: 'OK'
      }
    )
  }
}

// BAM文件验证
const validateBamPath = async (path, resultKey) => {
  await validatePath(path, 'bam_file', resultKey, 'BAM File')  // 修改这里
}

const validateGffFile = async () => {
  await validatePath(params.value.gffFile, 'gtf_file', 'gffFile', 'GFF File')  // 确保这里也正确
}

const getValidationStatus = (resultKey) => {
  const result = validationResults.value[resultKey]
  if (!result) return 'none'
  if (result.valid) return 'success'
  return 'error'
}

const getValidationMessage = (resultKey) => {
  const result = validationResults.value[resultKey]
  if (!result) return ''
  
  if (result.valid) {
    if (result.size_mb !== undefined) {
      return `Valid file (${result.size_mb} MB)`
    } else if (result.files && result.total_files !== undefined) {
      return `Found ${result.total_files}/${result.required_files} files`
    } else {
      return result.message || 'Valid file'
    }
  } else {
    return result.error || result.message || 'Invalid file'
  }
}

const generateSashimiPlot = async () => {
  if (!isFormValid.value) {
    ElMessage.error('Please fill in all required fields and ensure they are valid!')
    return
  }
  
  isGenerating.value = true
  statusMessage.value = ''
  
  try {
    // 构建命令
    let command = 'rmats2sashimiplot'
    command += ` -o ${props.outputDirectory}`
    command += ` -c "${autoCoordinate.value}"`
    command += ` --b1 ${params.value.b1}`
    command += ` --b2 ${params.value.b2}`
    
    // 样本标签
    const sample1Label = params.value.l1 || defaultSample1Label.value
    const sample2Label = params.value.l2 || defaultSample2Label.value
    command += ` --l1 "${sample1Label}"`
    command += ` --l2 "${sample2Label}"`
    
    // 显示参数
    command += ` --fig-width ${params.value.figWidth}`
    command += ` --fig-height ${params.value.figHeight}`
    command += ` --min-counts ${params.value.minCounts}`
    command += ` --exon_s ${params.value.exonS}`
    command += ` --intron_s ${params.value.intronS}`
    command += ` --font-size ${params.value.fontSize}`
    
    // 可选参数
    if (params.value.color) {
      command += ` --color '${params.value.color}'`
    }
    if (params.value.hideNumber) {
      command += ` --hide-number`
    }
    if (params.value.noTextBackground) {
      command += ` --no-text-background`
    }
    if (params.value.removeChrPrefix) {
      command += ` --remove-event-chr-prefix`
    }
    
    generatedCommand.value = command
    
    // 准备API参数
    const apiParams = {
      coordinate: autoCoordinate.value,
      output_directory: props.outputDirectory,
      sample1_bam: params.value.b1,
      sample2_bam: params.value.b2,
      sample1_label: sample1Label,
      sample2_label: sample2Label,
      figure_width: params.value.figWidth,
      figure_height: params.value.figHeight,
      min_counts: params.value.minCounts,
      exon_scale: params.value.exonS,
      intron_scale: params.value.intronS,
      font_size: params.value.fontSize,
      color: params.value.color || '',
      hide_number: params.value.hideNumber,
      no_text_background: params.value.noTextBackground,
      remove_chr_prefix: params.value.removeChrPrefix,
      generated_command: generatedCommand.value,
      timestamp: new Date().toISOString()
    }
    
    // 发送到API
    const response = await axios.post('/api/sashimi/generate-rmats', apiParams)
    
    if (response.data.success) {
      statusMessage.value = 'Sashimi plot generated successfully!'
      statusMessageClass.value = 'status-success'
      ElMessage.success('Sashimi plot generated successfully!')
      
      // 保存到Store
      store.addSashimiResult({
        success: true,
        command: generatedCommand.value,
        outputPath: response.data.output_directory,
        timestamp: new Date().toISOString(),
        isoformInfo: props.isoformData
      })
      
      // 发送事件到父组件
      emit('sashimi-generated', {
        command: generatedCommand.value,
        result: response.data,
        parameters: apiParams
      })
    } else {
      statusMessage.value = `Failed to generate sashimi plot: ${response.data.message}`
      statusMessageClass.value = 'status-error'
      ElMessage.error(`Failed to generate sashimi plot: ${response.data.message}`)
    }
    
  } catch (error) {
    console.error('Sashimi generation error:', error)
    statusMessage.value = `Error: ${error.message}`
    statusMessageClass.value = 'status-error'
    ElMessage.error(`Error: ${error.message}`)
  } finally {
    isGenerating.value = false
  }
}

// 初始化默认值
onMounted(() => {
  // 从Store中获取默认的GFF文件路径
  if (store.getFilePaths?.sashimiGff) {
    params.value.gffFile = store.getFilePaths.sashimiGff
  }
  
  // 从Store中获取默认的BAM文件路径（如果有的话）
  if (store.getFilePaths?.sashimiBam1) {
    params.value.b1 = store.getFilePaths.sashimiBam1
  }
  if (store.getFilePaths?.sashimiBam2) {
    params.value.b2 = store.getFilePaths.sashimiBam2
  }
  
  // 设置默认样本标签
  if (!params.value.l1) {
    params.value.l1 = defaultSample1Label.value
  }
  if (!params.value.l2) {
    params.value.l2 = defaultSample2Label.value
  }
})
const getConfig = () => {
  return {
    enabled: true,
    bam_files: {
      sample1: params.value.b1,
      sample2: params.value.b2
    },
    gff_file: params.value.gffFile,
    display_options: {
      sample1_label: params.value.l1 || defaultSample1Label.value,
      sample2_label: params.value.l2 || defaultSample2Label.value,
      figure_width: params.value.figWidth,
      figure_height: params.value.figHeight,
      min_counts: params.value.minCounts,
      exon_scale: params.value.exonS,
      intron_scale: params.value.intronS,
      font_size: params.value.fontSize,
      color: params.value.color,
      hide_number: params.value.hideNumber,
      no_text_background: params.value.noTextBackground,
      remove_chr_prefix: params.value.removeChrPrefix
    }
  }
}

// 暴露给父组件
defineExpose({
  getConfig,
  isFormValid
})
// 监听GFF文件变化，自动更新坐标
watch(() => params.value.gffFile, (newGffFile) => {
  if (newGffFile) {
    // 更新Store中的GFF文件路径
    store.updateFilePaths({
      sashimiGff: newGffFile
    })
  }
})
</script>

<style scoped>
.sashimi-config-container {
  background-color: #f8f9fa;
  border-radius: 8px;
  border: 1px solid #e9ecef;
  padding: 20px;
}

.params-section {
  background-color: white;
  padding: 20px;
  border-radius: 6px;
  border: 1px solid #dee2e6;
}

.section-title {
  font-size: 16px;
  font-weight: 600;
  margin-bottom: 16px;
  padding-bottom: 8px;
}

.section-title.required {
  color: #f56c6c;
  border-bottom: 2px solid #f56c6c;
}

.section-title.basic {
  color: #409eff;
  border-bottom: 2px solid #409eff;
}

.required-params {
  border-left: 4px solid #f56c6c;
  padding: 16px;
  margin-bottom: 24px;
  background-color: #fef0f0;
  border-radius: 4px;
}

.basic-options {
  border-left: 4px solid #409eff;
  padding: 16px;
  margin-bottom: 24px;
  background-color: #f0f7ff;
  border-radius: 4px;
}

.advanced-options {
  padding: 16px;
  background-color: #f8f9fa;
  border-radius: 4px;
}

.param-group {
  margin-bottom: 20px;
}

.param-group-title {
  font-size: 14px;
  font-weight: 600;
  color: #606266;
  margin-bottom: 12px;
  background-color: #f5f7fa;
  padding: 8px 12px;
  border-radius: 4px;
  border-left: 3px solid #409eff;
}

.inline-params {
  display: grid;
  grid-template-columns: 1fr 1fr 1fr;
  gap: 16px;
}

.checkbox-options {
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.validation-result {
  margin-top: 8px;
}

.gtf-details {
  font-size: 14px;
}

.gtf-details p {
  margin: 4px 0;
}

.sample-lines {
  margin-top: 8px;
}

.sample-lines code {
  display: block;
  background-color: #f5f5f5;
  padding: 4px 8px;
  border-radius: 4px;
  font-family: monospace;
  font-size: 12px;
  margin: 2px 0;
  word-break: break-all;
}

.generate-section {
  margin: 24px 0;
  text-align: center;
}

.command-output {
  margin-top: 20px;
  background-color: #f4f4f5;
  padding: 16px;
  border-radius: 6px;
  border: 1px solid #dcdfe6;
}

.command-output h5 {
  margin-top: 0;
  margin-bottom: 8px;
  color: #606266;
}

.command-output pre {
  background-color: #ffffff;
  padding: 12px;
  border-radius: 4px;
  border: 1px solid #e4e7ed;
  font-family: 'Courier New', monospace;
  white-space: pre-wrap;
  word-break: break-all;
  margin: 0;
  color: #303133;
}

.status-message {
  margin-top: 16px;
  padding: 12px 16px;
  border-radius: 6px;
  font-weight: 500;
}

.status-success {
  background-color: #f0f9ff;
  border: 1px solid #67c23a;
  color: #67c23a;
}

.status-error {
  background-color: #fef0f0;
  border: 1px solid #f56c6c;
  color: #f56c6c;
}

:deep(.el-form-item) {
  margin-bottom: 16px;
}

:deep(.el-form-item__label) {
  color: #606266;
  font-weight: 500;
}

:deep(.el-input__wrapper) {
  box-shadow: 0 0 0 1px #dcdfe6 inset;
}

:deep(.el-input__wrapper:hover) {
  box-shadow: 0 0 0 1px #c0c4cc inset;
}

:deep(.el-input__wrapper.is-focus) {
  box-shadow: 0 0 0 1px #409eff inset;
}

:deep(.el-input-number) {
  width: 120px;
}

:deep(.el-collapse-item__header) {
  font-weight: 600;
  color: #606266;
}

:deep(.el-tag) {
  font-weight: 500;
}

:deep(.el-button--large) {
  padding: 12px 24px;
  font-size: 16px;
}
</style>