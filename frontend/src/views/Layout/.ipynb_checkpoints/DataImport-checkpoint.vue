<template>
 <div>
   <h1>📤 Data Import</h1>
   <p>Please enter the corresponding file paths for single analysis</p>
   
   <!-- Single Analysis Mode Fixed -->
   <div class="analysis-type-section">
     <h3>🔍 Analysis Type</h3>
     <p class="type-desc">Single Analysis Mode</p>
     
     <el-alert
       title="Single Analysis Mode"
       type="success"
       :closable="false"
       show-icon
     >
       <template #default>
         Single analysis mode for processing individual rMATS results.
       </template>
     </el-alert>
   </div>

   <!-- File path inputs -->
   <div class="file-paths-section">
     
     <!-- rMATS Results Path -->
     <div class="path-section">
       <h3>📁 rMATS Results Path(Required)</h3>
       <p class="path-desc">Please enter the path to your rMATS results directory</p>
       <el-input
         v-model="tempFilePaths.single"
         placeholder="e.g.: /data/rmats/results/"
         clearable
         size="large"
         @input="updateStorePath('single')"
         @blur="validateRmatsPath(tempFilePaths.single, 'single')"
       >
         <template #prepend>
           <span>Results Path</span>
         </template>
         <template #suffix>
           <el-icon v-if="getValidationStatus('single') === 'success'" style="color: #67c23a">
             <CircleCheck />
           </el-icon>
           <el-icon v-else-if="getValidationStatus('single') === 'error'" style="color: #f56c6c">
             <CircleClose />
           </el-icon>
         </template>
       </el-input>
       
       <!-- 验证结果显示 -->
       <div v-if="validationResults.single" class="validation-result">
         <el-alert
           v-if="validationResults.single.valid"
           type="success"
           :closable="false"
           :title="`✅ ${getValidationMessage('single')}`"
         >
           <template #default>
             <div class="file-list">
               <div 
                 v-for="file in getFileDetails('single')" 
                 :key="file.name"
                 class="file-item"
                 :class="{ 'file-missing': !file.exists }"
               >
                 <el-icon v-if="file.exists" style="color: #67c23a; margin-right: 8px">
                   <Check />
                 </el-icon>
                 <el-icon v-else style="color: #f56c6c; margin-right: 8px">
                   <Close />
                 </el-icon>
                 <span class="file-name">{{ file.name }}</span>
                 <span v-if="file.exists" class="file-size">({{ file.size_mb }} MB)</span>
               </div>
             </div>
           </template>
         </el-alert>
       </div>
       
       <div v-if="store.getFilePaths.single && !validationResults.single" class="path-info">
         <el-tag type="success" size="small">✓ rMATS results path has been set</el-tag>
         <p class="path-display">{{ store.getFilePaths.single }}</p>
       </div>
     </div>

     <!-- Output Directory Path -->
     <div class="path-section">
       <h3>📂 Output Directory Path(Required)</h3>
       <p class="path-desc">Please enter the path for saving analysis results</p>
       <el-input
         v-model="tempFilePaths.singleOutput"
         placeholder="e.g.: /data/results/analysis_output"
         clearable
         size="large"
         @input="updateStorePath('singleOutput')"
         @blur="validateOutputPath(tempFilePaths.singleOutput, 'singleOutput')"
       >
         <template #prepend>
           <span>Output Path</span>
         </template>
         <template #suffix>
           <el-icon v-if="getValidationStatus('singleOutput') === 'success'" style="color: #67c23a">
             <CircleCheck />
           </el-icon>
           <el-icon v-else-if="getValidationStatus('singleOutput') === 'error'" style="color: #f56c6c">
             <CircleClose />
           </el-icon>
         </template>
       </el-input>
       
       <!-- 输出目录验证结果 -->
       <div v-if="validationResults.singleOutput" class="validation-result">
         <el-alert
           v-if="validationResults.singleOutput.valid"
           type="success"
           :closable="false"
           show-icon
           :title="`✅ ${getValidationMessage('singleOutput')}`"
         />
         
         <el-alert
           v-else-if="validationResults.singleOutput.exists === false"
           type="warning"
           :closable="false"
           show-icon
           title="⚠️ Directory does not exist"
         >
           <template #default>
             <p>The specified directory does not exist. Please check the path or create the directory manually.</p>
           </template>
         </el-alert>
         
         <el-alert
           v-else-if="!validationResults.singleOutput.writable"
           type="error"
           :closable="false"
           show-icon
           title="❌ Directory is not writable"
         >
           <template #default>
             <p>The directory exists but is not writable. Please check the permissions.</p>
           </template>
         </el-alert>
       </div>
       
       <div v-if="store.getFilePaths.singleOutput && !validationResults.singleOutput" class="path-info">
         <el-tag type="primary" size="small">✓ Output directory path has been set</el-tag>
         <p class="path-display">{{ store.getFilePaths.singleOutput }}</p>
       </div>
     </div>

<!-- Sample 1 BAM List File Path (Required) -->
     <div class="path-section">
       <h3>📋 Sample 1 BAM List File()</h3>
       <p class="bam-desc">Please provide the TXT file containing BAM file paths for Sample 1 (Required)</p>
       <el-input
         v-model="tempFilePaths.sample1BamList"
         placeholder="e.g.: /data/sample1_bams.txt"
         clearable
         size="large"
         @input="handleBamListInput('sample1BamList')"
         @blur="handleBamListBlur"
         class="required-field"
       >
         <template #prepend>
           <span class="required-label">Sample 1 TXT *</span>
         </template>
         <template #suffix>
           <el-icon v-if="getValidationStatus('sample1BamList') === 'success'" style="color: #67c23a">
             <CircleCheck />
           </el-icon>
           <el-icon v-else-if="getValidationStatus('sample1BamList') === 'error'" style="color: #f56c6c">
             <CircleClose />
           </el-icon>
         </template>
       </el-input>
       
       <!-- Sample 1 BAM List 验证结果 -->
       <div v-if="validationResults.sample1BamList" class="validation-result">
         <el-alert
           v-if="validationResults.sample1BamList.valid"
           type="success"
           :closable="false"
           show-icon
           :title="`✅ ${getValidationMessage('sample1BamList')}`"
         >
           <template #default>
             <div class="bam-details">
               <div v-if="validationResults.sample1BamList.sample_paths?.length > 0">
                 <p><strong>Sample BAM paths:</strong></p>
                 <div class="sample-lines">
                   <code v-for="(path, index) in validationResults.sample1BamList.sample_paths.slice(0, 3)" :key="index">
                     {{ path }}
                   </code>
                   <p v-if="validationResults.sample1BamList.sample_paths.length > 3" class="more-files">
                     ... and {{ validationResults.sample1BamList.sample_paths.length - 3 }} more files
                   </p>
                 </div>
               </div>
             </div>
           </template>
         </el-alert>
         
         <!-- 显示验证错误 -->
         <el-alert
           v-else
           type="error"
           :closable="false"
           show-icon
           title="❌ Sample 1 BAM List Validation Failed"
         >
           <template #default>
             <p>{{ validationResults.sample1BamList.error || 'BAM list file validation failed' }}</p>
             <p class="error-hint">Please ensure the TXT file exists and contains valid BAM file paths.</p>
           </template>
         </el-alert>
       </div>
       
       <div v-if="store.getFilePaths.sample1BamList && !validationResults.sample1BamList" class="path-info">
         <el-tag type="warning" size="small">⚠ Sample 1 BAM list path set - validation pending</el-tag>
         <p class="path-display">{{ store.getFilePaths.sample1BamList }}</p>
       </div>
     </div>

     <!-- Sample 2 BAM List File Path (Required) -->
     <div class="path-section">
       <h3>📋 Sample 2 BAM List File *</h3>
       <p class="bam-desc">Please provide the TXT file containing BAM file paths for Sample 2 (Required)</p>
       <el-input
         v-model="tempFilePaths.sample2BamList"
         placeholder="e.g.: /data/sample2_bams.txt"
         clearable
         size="large"
         @input="handleBamListInput('sample2BamList')"
         @blur="handleBamListBlur"
         class="required-field"
       >
         <template #prepend>
           <span class="required-label">Sample 2 TXT *</span>
         </template>
         <template #suffix>
           <el-icon v-if="getValidationStatus('sample2BamList') === 'success'" style="color: #67c23a">
             <CircleCheck />
           </el-icon>
           <el-icon v-else-if="getValidationStatus('sample2BamList') === 'error'" style="color: #f56c6c">
             <CircleClose />
           </el-icon>
         </template>
       </el-input>
       
       <!-- Sample 2 BAM List 验证结果 -->
       <div v-if="validationResults.sample2BamList" class="validation-result">
         <el-alert
           v-if="validationResults.sample2BamList.valid"
           type="success"
           :closable="false"
           show-icon
           :title="`✅ ${getValidationMessage('sample2BamList')}`"
         >
           <template #default>
             <div class="bam-details">
               <div v-if="validationResults.sample2BamList.sample_paths?.length > 0">
                 <p><strong>Sample BAM paths:</strong></p>
                 <div class="sample-lines">
                   <code v-for="(path, index) in validationResults.sample2BamList.sample_paths.slice(0, 3)" :key="index">
                     {{ path }}
                   </code>
                   <p v-if="validationResults.sample2BamList.sample_paths.length > 3" class="more-files">
                     ... and {{ validationResults.sample2BamList.sample_paths.length - 3 }} more files
                   </p>
                 </div>
               </div>
             </div>
           </template>
         </el-alert>
         
         <!-- 显示验证错误 -->
         <el-alert
           v-else
           type="error"
           :closable="false"
           show-icon
           title="❌ Sample 2 BAM List Validation Failed"
         >
           <template #default>
             <p>{{ validationResults.sample2BamList.error || 'BAM list file validation failed' }}</p>
             <p class="error-hint">Please ensure the TXT file exists and contains valid BAM file paths.</p>
           </template>
         </el-alert>
       </div>
       
       <div v-if="store.getFilePaths.sample2BamList && !validationResults.sample2BamList" class="path-info">
         <el-tag type="warning" size="small">⚠ Sample 2 BAM list path set - validation pending</el-tag>
         <p class="path-display">{{ store.getFilePaths.sample2BamList }}</p>
       </div>
     </div>

     <!-- GTF File Path -->
     <div class="path-section">
       <h4>📄 GTF File Path (Required)</h4>
       <p class="gtf-desc">Optionally provide the GTF annotation file path</p>
        <el-input
          v-model="tempFilePaths.singleGtf"
          placeholder="e.g.: /data/annotations/genome.gtf"
          clearable
          size="large"
          class="required-field"
          @input="handleGtfInput"
          @blur="validateGtfPath(tempFilePaths.singleGtf, 'singleGtf')"
        >
         <template #prepend>
           <span class="required-label">GTF Path *</span>
         </template>
         <template #suffix>
           <el-icon v-if="getValidationStatus('singleGtf') === 'success'" style="color: #67c23a">
             <CircleCheck />
           </el-icon>
           <el-icon v-else-if="getValidationStatus('singleGtf') === 'error'" style="color: #f56c6c">
             <CircleClose />
           </el-icon>
         </template>
       </el-input>
       
       <!-- GTF文件验证结果 -->
       <div v-if="validationResults.singleGtf" class="validation-result">
         <el-alert
           v-if="validationResults.singleGtf.valid"
           type="success"
           :closable="false"
           show-icon
           :title="`✅ ${getValidationMessage('singleGtf')}`"
         >
           <template #default>
             <div class="gtf-details">
               <p><strong>File size:</strong> {{ validationResults.singleGtf.size_mb }} MB</p>
               <p><strong>Compressed:</strong> {{ validationResults.singleGtf.is_compressed ? 'Yes' : 'No' }}</p>
               <div v-if="validationResults.singleGtf.sample_lines?.length > 0">
                 <p><strong>Sample content:</strong></p>
                 <div class="sample-lines">
                   <code v-for="(line, index) in validationResults.singleGtf.sample_lines.slice(0, 2)" :key="index">
                     {{ line }}
                   </code>
                 </div>
               </div>
             </div>
           </template>
         </el-alert>
       </div>
       
        <div v-if="store.getFilePaths.singleGtf && !validationResults.singleGtf" class="path-info">
          <el-tag type="warning" size="small">⚠ GTF file path set - validation pending</el-tag>
          <p class="path-display">{{ store.getFilePaths.singleGtf }}</p>
        </div>
     </div>

     <!-- Single Analysis Sashimi Plot Configuration -->
     <div class="path-section">
       <div class="sashimi-header">
         <h3>🧬 Sashimi Plot Configuration (Optional)</h3>
         <p class="sashimi-desc">Generate Sashimi plots for visualizing splicing events</p>
         <el-switch 
           v-model="enableSingleSashimi"
           active-text="Enable Sashimi Plot"
           inactive-text="Skip Sashimi Plot"
           style="margin-bottom: 16px;"
         />
       </div>
       
        <SashimiPlotConfig
          v-if="enableSingleSashimi"
          ref="sashimiConfigRef"
          v-model="enableSingleSashimi"
          :output-directory="store.getFilePaths.singleOutput"
          :rmats-results-path="store.getFilePaths.single"
          analysis-type="single"
          @sashimi-generated="handleSashimiGenerated"
        />
     </div>

     <!-- Send button -->
     <div class="send-section">
       <el-button 
         type="primary"
         size="large"
         @click="sendPathsData"
         :loading="store.getAnalysisStatus === 'loading'"
         :disabled="!isFormValid"
         icon="Upload"
       >
         {{ store.getAnalysisStatus === 'loading' ? 'Sending...' : 'Start Single Analysis' }}
       </el-button>
       
       <div v-if="!isFormValid" class="path-hint">
         <el-alert
           title="Please complete required fields"
           type="warning"
           :closable="false"
           show-icon
         >
           <template #default>
             Please fill in all required paths to proceed with the analysis.
           </template>
         </el-alert>
       </div>
     </div>
     
     <!-- Status display -->
     <div v-if="statusMessage" :class="statusMessageClass" class="status-message">
       {{ statusMessage }}
     </div>

     <!-- Path preview card -->
     <div v-if="isFormValid" class="preview-card">
       <h4>📋 Path Preview</h4>
       <div class="preview-content">
         <div class="preview-item">
           <strong>rMATS Results File:</strong>
           <code>{{ store.getFilePaths.single }}</code>
         </div>
         <div class="preview-item">
           <strong>Output Directory:</strong>
           <code>{{ store.getFilePaths.singleOutput }}</code>
         </div>
         <div class="preview-item">
           <strong>Sample 1 BAM List: *</strong>
           <code>{{ store.getFilePaths.sample1BamList || 'Not set (Required)' }}</code>
         </div>
         <div class="preview-item">
           <strong>Sample 2 BAM List: *</strong>
           <code>{{ store.getFilePaths.sample2BamList || 'Not set (Required)' }}</code>
         </div>
        <div class="preview-item">
          <strong>GTF Annotation File: *</strong>
          <code>{{ store.getFilePaths.singleGtf || 'Not set (Required)' }}</code>
        </div>
         <div v-if="enableSingleSashimi" class="preview-item">
           <strong>Sashimi Plot:</strong>
           <el-tag type="success" size="small">Enabled</el-tag>
         </div>
       </div>
     </div>
     
     <!-- API Response Message -->
     <div v-if="apiResponse" class="api-response">
       <el-alert
         :title="apiResponse.success ? 'API Response: Success' : 'API Response: Error'"
         :type="apiResponse.success ? 'success' : 'error'"
         :closable="true"
         show-icon
       >
         <template #default>
           <p>{{ apiResponse.message }}</p>
           <div v-if="apiResponse.success && apiResponse.output_directory" class="response-detail">
             <strong>Output Directory:</strong>
             <code>{{ apiResponse.output_directory }}</code>
           </div>
         </template>
       </el-alert>
     </div>

     <!-- Sashimi Generation Results -->
     <div v-if="sashimiResults.length > 0" class="sashimi-results">
       <h4>🧬 Sashimi Plot Results</h4>
       <div v-for="(result, index) in sashimiResults" :key="index" class="sashimi-result-item">
         <el-alert
           :title="`Sashimi Plot ${index + 1}: ${result.success ? 'Generated Successfully' : 'Generation Failed'}`"
           :type="result.success ? 'success' : 'error'"
           :closable="true"
           show-icon
         >
           <template #default>
             <p>{{ result.message || result.error }}</p>
             <div v-if="result.success && result.outputPath" class="sashimi-output">
               <strong>Output Location:</strong>
               <code>{{ result.outputPath }}</code>
             </div>
           </template>
         </el-alert>
       </div>
     </div>

   </div>
 </div>
</template>
<script setup>
import { ref, computed, onMounted, watch } from 'vue'
import { useAppStore } from '@/stores/app'
import { useRouter } from 'vue-router'
import { ElMessageBox, ElMessage } from 'element-plus'
import axios from 'axios'
import SashimiPlotConfig from './ImportSashimi.vue'
import { CircleCheck, CircleClose, Check, Close } from '@element-plus/icons-vue'

// 使用 store
const store = useAppStore()
const router = useRouter()

const enableSingleSashimi = ref(false)
const sashimiResults = ref([])

// 在现有的 ref 声明附近添加
const sashimiConfigRef = ref(null)

// 获取 Sashimi 配置的方法
const getSashimiConfig = () => {
  if (!enableSingleSashimi.value || !sashimiConfigRef.value) {
    return null
  }
  
  return sashimiConfigRef.value.getConfig()
}
const handleGtfInput = () => {
  updateStorePath('singleGtf')
  
  // 直接确保GTF状态更新
  const gtfPath = tempFilePaths.value.singleGtf?.trim() || ''
  store.setGtfFilePath(gtfPath)
  
  // 添加必填验证
  if (!gtfPath) {
    validationResults.value.singleGtf = {
      valid: false,
      error: 'GTF file path is required',
      exists: false
    }
  }
  
  console.log('GTF input handled:', {
    path: gtfPath,
    enabled: store.isGtfEnabled
  })
}

// 本地临时路径状态，只保留单一分析相关的路径
const tempFilePaths = ref({
  single: '',
  singleGtf: '',
  singleOutput: '',
  sample1BamList: '',    // 样本1的BAM列表文件
  sample2BamList: ''     // 样本2的BAM列表文件
})

const validationResults = ref({
  single: null,
  singleGtf: null,
  singleOutput: null,
  sample1BamList: null,  // 样本1的BAM列表验证结果
  sample2BamList: null   // 样本2的BAM列表验证结果
})

const statusMessage = ref('')
const statusMessageClass = ref('')
const apiResponse = ref(null)

const handleSashimiGenerated = (result) => {
  console.log('Sashimi plot generated:', result)
  
  sashimiResults.value.push({
    success: result.result?.success || false,
    message: result.result?.message || '',
    error: result.result?.error || '',
    outputPath: result.result?.output_directory || '',
    command: result.command,
    timestamp: new Date().toLocaleString()
  })
  
  if (result.result?.success) {
    ElMessage.success('Sashimi plot generated successfully!')
  } else {
    ElMessage.error('Failed to generate Sashimi plot')
  }
}

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
    
    if (!response.data.valid) {
      await ElMessageBox.alert(
        response.data.error || `${fieldName} validation failed`,
        `${fieldName} Validation Failed`,
        {
          type: 'error',
          confirmButtonText: 'OK'
        }
      )
    } else {
      ElMessage.success({
        message: `${fieldName} validation successful`,
        duration: 2000
      })
    }
    
  } catch (error) {
    console.error(`❌ ${fieldName}验证出错:`, error)
    validationResults.value[resultKey] = {
      valid: false,
      error: `Validation failed: ${error.message || 'Unknown error'}`,
      exists: false
    }
    
    await ElMessageBox.alert(
      `${fieldName} validation failed: ${error.message || 'Unknown error'}`,
      'Validation Error',
      {
        type: 'error',
        confirmButtonText: 'OK'
      }
    )
  }
}

// 新的双BAM列表验证方法
const validateBothBamLists = async () => {
  const sample1Path = tempFilePaths.value.sample1BamList?.trim()
  const sample2Path = tempFilePaths.value.sample2BamList?.trim()
  
  if (!sample1Path || !sample2Path) {
    // 如果任一路径为空，清除验证结果
    validationResults.value.sample1BamList = null
    validationResults.value.sample2BamList = null
    return
  }
  
  try {
    console.log('🧬 验证两个BAM列表文件:', { sample1Path, sample2Path })
    
    const response = await axios.post('/api/validate-bam-lists', {
      sample1_path: sample1Path,
      sample2_path: sample2Path
    })
    
    console.log('✅ 双BAM列表验证结果:', response.data)
    
    // 分别设置两个验证结果
    validationResults.value.sample1BamList = response.data.sample1_result
    validationResults.value.sample2BamList = response.data.sample2_result
    
    // 显示验证消息
    if (response.data.sample1_result?.valid && response.data.sample2_result?.valid) {
      ElMessage.success('Both BAM list files validated successfully!')
    } else {
      const errors = []
      if (!response.data.sample1_result?.valid) {
        errors.push(`Sample 1: ${response.data.sample1_result?.error}`)
      }
      if (!response.data.sample2_result?.valid) {
        errors.push(`Sample 2: ${response.data.sample2_result?.error}`)
      }
      ElMessage.error(`BAM validation failed: ${errors.join('; ')}`)
    }
    
  } catch (error) {
    console.error('❌ 双BAM列表验证出错:', error)
    
    const errorMsg = `Validation failed: ${error.message || 'Unknown error'}`
    
    validationResults.value.sample1BamList = {
      valid: false,
      error: errorMsg,
      exists: false
    }
    validationResults.value.sample2BamList = {
      valid: false,
      error: errorMsg,
      exists: false
    }
    
    ElMessage.error('Failed to validate BAM list files')
  }
}

// 新的BAM列表输入处理方法
const handleBamListInput = (pathType) => {
  updateStorePath(pathType)
  
  // 当两个路径都有值时，触发双重验证（延迟执行）
  if (tempFilePaths.value.sample1BamList?.trim() && tempFilePaths.value.sample2BamList?.trim()) {
    setTimeout(() => {
      validateBothBamLists()
    }, 500)
  }
}

// 新的BAM列表失焦处理方法
const handleBamListBlur = () => {
  // 当失去焦点时，如果两个路径都有值，进行验证
  if (tempFilePaths.value.sample1BamList?.trim() && tempFilePaths.value.sample2BamList?.trim()) {
    validateBothBamLists()
  }
}

// 保留原有的单独验证方法（可能其他地方还在使用）
const validateBamListPath = async (path, resultKey) => {
  if (!path || !path.trim()) {
    validationResults.value[resultKey] = {
      valid: false,
      error: 'BAM list file path is required',
      exists: false
    }
    
    const sampleName = resultKey === 'sample1BamList' ? 'Sample 1' : 'Sample 2'
    await ElMessageBox.alert(
      `${sampleName} BAM list file path is required for rMATS analysis`,
      `${sampleName} BAM List Required`,
      {
        type: 'error',
        confirmButtonText: 'OK'
      }
    )
    return
  }
  
  // 使用 gtf_file 类型验证TXT文件（临时方案）
  await validatePath(path, 'gtf_file', resultKey, resultKey === 'sample1BamList' ? 'Sample 1 BAM List' : 'Sample 2 BAM List')
}

// 具体的验证方法
const validateRmatsPath = async (path, resultKey = 'single') => {
  await validatePath(path, 'rmats_dir', resultKey, 'rMATS Directory')
}

const validateOutputPath = async (path, resultKey = 'singleOutput') => {
  await validatePath(path, 'output_dir', resultKey, 'Output Directory')
}

const validateGtfPath = async (path, resultKey = 'singleGtf') => {
  await validatePath(path, 'gtf_file', resultKey, 'GTF File')
}

const getValidationStatus = (resultKey) => {
  const result = validationResults.value[resultKey]
  if (!result) return 'none'
  if (result.valid) return 'success'
  return 'error'
}

const getFileDetails = (resultKey) => {
  const result = validationResults.value[resultKey]
  return result?.file_details || []
}

// 从 store 同步数据到临时变量
const syncStoreToTemp = () => {
  tempFilePaths.value = {
    single: store.getFilePaths?.single || '',
    singleGtf: store.getFilePaths?.singleGtf || '',
    singleOutput: store.getFilePaths?.singleOutput || '',
    sample1BamList: store.getFilePaths?.sample1BamList || '',  // 样本1 BAM列表
    sample2BamList: store.getFilePaths?.sample2BamList || ''   // 样本2 BAM列表
  }
}

// 组件挂载时同步数据
onMounted(() => {
  syncStoreToTemp()
  console.log('DataImport component mounted')
})

// 监听 store 中的路径变化
watch(() => [
  store.getFilePaths?.single,
  store.getFilePaths?.singleOutput,
  store.getFilePaths?.singleGtf,
  store.getFilePaths?.sample1BamList,  // 样本1 BAM列表
  store.getFilePaths?.sample2BamList   // 样本2 BAM列表
], () => {
  syncStoreToTemp()
})

// 更新路径到 store
const updateStorePath = (pathKey) => {
  const updatedPaths = {}
  const currentValue = tempFilePaths.value[pathKey]
  
  updatedPaths[pathKey] = currentValue
  console.log(`更新路径 ${pathKey} = ${currentValue}`)
  
  store.updateFilePaths(updatedPaths)
  
  if (currentValue.trim() !== '') {
    store.savePathHistory()
  }
}

// 表单验证
const isFormValid = computed(() => {
  // 单一分析：rMATS路径、输出路径、两个BAM列表文件、GTF文件都是必填的
  const hasRequiredPaths = store.getFilePaths?.single && 
                          store.getFilePaths?.singleOutput && 
                          store.getFilePaths?.sample1BamList &&  // 样本1 BAM列表必填
                          store.getFilePaths?.sample2BamList &&  // 样本2 BAM列表必填
                          store.getFilePaths?.singleGtf          // GTF文件必填
  
  const pathsValid = validationResults.value.single?.valid && 
                    validationResults.value.singleOutput?.valid &&
                    validationResults.value.sample1BamList?.valid &&  // 样本1 BAM列表验证
                    validationResults.value.sample2BamList?.valid &&  // 样本2 BAM列表验证
                    validationResults.value.singleGtf?.valid          // GTF文件验证
  
  console.log('Single mode validation (GTF now required):', {
    hasRequiredPaths,
    pathsValid,
    final: hasRequiredPaths && pathsValid
  })
  
  return hasRequiredPaths && pathsValid
})

const getValidationMessage = (resultKey) => {
  const result = validationResults.value[resultKey]
  if (!result) return ''
  
  if (result.valid) {
    if (result.files && result.total_files !== undefined) {
      return `Found ${result.total_files}/${result.required_files} rMATS files`
    } else if (result.writable !== undefined) {
      return result.writable ? 'Directory exists and is writable' : 'Directory exists but is not writable'
    } else {
      return result.message || 'Valid'
    }
  } else {
    return result.error || result.message || 'Invalid'
  }
}

// 准备API请求数据
const prepareApiRequestData = () => {
  const requestData = {
    rmats_dir: store.getFilePaths.single,
    output_dir: store.getFilePaths.singleOutput,
    sample1_bam_list: store.getFilePaths.sample1BamList,
    sample2_bam_list: store.getFilePaths.sample2BamList,
    gtf_file: store.getFilePaths.singleGtf || '',
    analysis_type: 'single',
    
    // 添加 Sashimi 配置
    sashimi_config: getSashimiConfig()
  }
  
  console.log('Preparing Single analysis API request:', requestData)
  return requestData
}

// 发送数据到API
const sendPathsData = async () => {
 if (!isFormValid.value) return
 
 store.setAnalysisStatus('loading')
 statusMessage.value = 'Cleaning up previous files and starting analysis...'
 statusMessageClass.value = 'status-info'
 apiResponse.value = null
 
 try {
   // 静默清理文件（不影响主流程）
   try {
     await axios.post('/api/cleanup/generated-files')
     console.log('File cleanup completed')
   } catch (cleanupError) {
     console.warn('File cleanup failed, continuing:', cleanupError.message)
   }
   
   statusMessage.value = 'Sending file paths for single analysis...'
   
   const requestData = prepareApiRequestData()
   
   console.log('Sending API request:', requestData)
   
   const response = await axios.post('/api/process-rmats', requestData, {
     headers: {
       'Content-Type': 'application/json'
     }
   })
   
   console.log('API response:', response.data)
   
   apiResponse.value = response.data
   
   if (response.data.success) {
     statusMessage.value = 'File paths sent successfully! Single analysis started.'
     statusMessageClass.value = 'status-success'
     store.setAnalysisStatus('completed')
     
     // 标记数据已加载
     store.setDataLoaded()
     console.log('数据已标记为加载完成，侧边栏菜单应该启用')
     
     // 更新 Sashimi 状态
     if (enableSingleSashimi.value) {
       const sashimiConfig = getSashimiConfig()
       if (sashimiConfig) {
         store.setSashimiEnabled(true)
         store.setSashimiConfiguration(sashimiConfig)
         console.log('Sashimi 状态已更新为启用')
       }
     }
     
     ElMessage.success({
       message: 'Data imported and analyzed successfully!',
       duration: 3000
     })
     
     router.push({ name: 'single-analysis' })
   } else {
     statusMessage.value = `Failed to process data: ${response.data.message || 'Unknown error'}`
     statusMessageClass.value = 'status-error'
     store.setAnalysisStatus('error')
     // 确保失败时状态正确
     store.clearDataLoaded()
   }
   
 } catch (error) {
   console.error('API request error:', error)
   
   statusMessage.value = `Failed to send file paths. Error: ${error.message || 'Unknown error'}`
   statusMessageClass.value = 'status-error'
   store.setAnalysisStatus('error')
   
   // 确保错误时清除加载状态
   store.clearDataLoaded()
   
   apiResponse.value = {
     success: false,
     message: `API request failed: ${error.message || 'Unknown error'}`
   }
 } finally {
   setTimeout(() => {
     statusMessage.value = ''
   }, 5000)
 }
}
</script>

<style scoped>
/* 保持原有的所有样式不变 */
.analysis-type-section {
  margin-bottom: 32px;
  padding: 20px;
  background-color: #f8f9fa;
  border-radius: 8px;
  border: 1px solid #e9ecef;
}

.analysis-type-section h3 {
  margin-top: 0;
  margin-bottom: 8px;
  color: #495057;
}

.type-desc {
  margin-bottom: 16px;
  color: #6c757d;
  font-size: 14px;
}

.file-paths-section {
  margin-bottom: 24px;
}

.path-section {
  margin-bottom: 32px;
  padding: 20px;
  background-color: #f8f9fa;
  border-radius: 8px;
  border: 1px solid #e9ecef;
}

.path-section h3 {
  margin-top: 0;
  margin-bottom: 8px;
  color: #495057;
}

.path-desc {
  margin-bottom: 16px;
  color: #6c757d;
  font-size: 14px;
}

.path-info {
  margin-top: 12px;
  padding: 12px;
  background-color: #ffffff;
  border-radius: 6px;
  border: 1px solid #dee2e6;
}

.path-display {
  margin-top: 8px;
  margin-bottom: 0;
  font-family: monospace;
  font-size: 14px;
  color: #495057;
  word-break: break-all;
}

.send-section {
  margin-bottom: 32px;
  text-align: center;
}

.path-hint {
  margin-top: 16px;
  max-width: 600px;
  margin-left: auto;
  margin-right: auto;
}

/* 必填字段样式 */
.required-field {
  position: relative;
}

.required-label {
  font-weight: 600;
  color: #495057;
}

.bam-desc {
  margin-bottom: 16px;
  color: #495057;
  font-size: 14px;
  font-weight: 500;
}

.error-hint {
  margin-top: 8px;
  font-size: 13px;
  color: #6c757d;
  font-style: italic;
}

.bam-details {
  background-color: #f8f9fa;
  padding: 12px;
  border-radius: 4px;
  margin-top: 8px;
  border-left: 4px solid #67c23a;
}

.bam-details p {
  margin-bottom: 8px;
  font-size: 14px;
  line-height: 1.4;
}

.bam-details p:last-child {
  margin-bottom: 0;
}

.bam-details strong {
  color: #495057;
  margin-right: 8px;
  font-weight: 600;
}

/* 必填字段验证失败时的样式 */
:deep(.required-field .el-input__wrapper.is-focus) {
  box-shadow: 0 0 0 1px #f56c6c inset;
}

:deep(.required-field .el-input__wrapper:hover) {
  box-shadow: 0 0 0 1px #f56c6c inset;
}

:deep(.required-field .el-input-group__prepend) {
  background-color: #fff5f5;
  border-color: #f56c6c;
  color: #495057;
  font-weight: 600;
}

.status-message {
  margin-bottom: 24px;
  padding: 12px 16px;
  border-radius: 6px;
  font-size: 14px;
  text-align: center;
}

.status-info {
  background-color: #d1ecf1;
  color: #0c5460;
  border: 1px solid #bee5eb;
}

.status-success {
  background-color: #d4edda;
  color: #155724;
  border: 1px solid #c3e6cb;
}

.status-error {
  background-color: #f8d7da;
  color: #721c24;
  border: 1px solid #f5c6cb;
}

.preview-card {
  padding: 20px;
  background-color: #ffffff;
  border-radius: 8px;
  border: 1px solid #dee2e6;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.preview-card h4 {
  margin-top: 0;
  margin-bottom: 16px;
  color: #495057;
}

.preview-content {
  background-color: #f8f9fa;
  padding: 16px;
  border-radius: 6px;
  border: 1px solid #e9ecef;
}

.preview-item {
  margin-bottom: 12px;
}

.preview-item:last-child {
  margin-bottom: 0;
}

.preview-item strong {
  color: #495057;
  display: block;
  margin-bottom: 4px;
}

.preview-item code {
  background-color: #e9ecef;
  padding: 4px 8px;
  border-radius: 4px;
  font-family: monospace;
  font-size: 13px;
  color: #495057;
  word-break: break-all;
}

.api-response {
  margin-top: 24px;
  margin-bottom: 24px;
}

.response-detail {
  margin-top: 8px;
  font-size: 14px;
}

.response-detail code {
  background-color: #e9ecef;
  padding: 4px 8px;
  border-radius: 4px;
  font-family: monospace;
  font-size: 13px;
  color: #495057;
  word-break: break-all;
  display: block;
  margin-top: 4px;
}

h1 {
  color: #2c3e50;
  margin-bottom: 8px;
}

p {
  color: #6c757d;
  margin-bottom: 24px;
}

.sashimi-header {
  margin-bottom: 16px;
}

.sashimi-header h3 {
  margin-bottom: 8px;
  color: #495057;
}

.sashimi-desc {
  margin-bottom: 12px;
  color: #6c757d;
  font-size: 14px;
  font-style: italic;
}

.sashimi-results {
  margin-top: 24px;
  padding: 20px;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #e9ecef;
}

.sashimi-results h4 {
  margin-top: 0;
  margin-bottom: 16px;
  color: #495057;
  border-bottom: 2px solid #dee2e6;
  padding-bottom: 6px;
}

.sashimi-result-item {
  margin-bottom: 12px;
}

.sashimi-result-item:last-child {
  margin-bottom: 0;
}

.sashimi-output {
  margin-top: 8px;
  font-size: 14px;
}

.sashimi-output code {
  background-color: #e9ecef;
  padding: 4px 8px;
  border-radius: 4px;
  font-family: monospace;
  font-size: 13px;
  color: #495057;
  word-break: break-all;
  display: block;
  margin-top: 4px;
}

.gtf-desc {
  margin-bottom: 16px;
  color: #6c757d;
  font-size: 14px;
  font-style: italic;
}

.validation-result {
  margin-top: 12px;
}

.file-list {
  margin-top: 8px;
}

.file-item {
  display: flex;
  align-items: center;
  margin-bottom: 8px;
  padding: 8px;
  background-color: #ffffff;
  border-radius: 4px;
  border: 1px solid #dee2e6;
}

.file-item:last-child {
  margin-bottom: 0;
}

.file-missing {
  background-color: #fff5f5;
  border-color: #fecaca;
}

.file-name {
  flex: 1;
  font-family: monospace;
  font-size: 13px;
  color: #495057;
}

.file-size {
  font-size: 12px;
  color: #6c757d;
  margin-left: 8px;
}

.gtf-details {
  background-color: #f8f9fa;
  padding: 12px;
  border-radius: 4px;
  margin-top: 8px;
  border-left: 4px solid #67c23a;
}

.gtf-details p {
  margin-bottom: 8px;
  font-size: 14px;
  line-height: 1.4;
}

.gtf-details p:last-child {
  margin-bottom: 0;
}

.gtf-details strong {
  color: #495057;
  margin-right: 8px;
  font-weight: 600;
}

.sample-lines {
  margin-top: 8px;
  background-color: #f1f3f4;
  padding: 8px;
  border-radius: 4px;
  border: 1px solid #e1e5e9;
}

.sample-lines code {
  display: block;
  font-family: monospace;
  font-size: 12px;
  color: #495057;
  margin-bottom: 4px;
  line-height: 1.3;
}

.sample-lines code:last-child {
  margin-bottom: 0;
}

.more-files {
  font-size: 12px;
  color: #6c757d;
  font-style: italic;
  margin-top: 4px;
}

:deep(.el-input-group__prepend) {
  background-color: #e9ecef;
  border-color: #ced4da;
  color: #495057;
}

:deep(.el-input__wrapper) {
  box-shadow: 0 0 0 1px #ced4da inset;
}

:deep(.el-input__wrapper:hover) {
  box-shadow: 0 0 0 1px #80bdff inset;
}

:deep(.el-input__wrapper.is-focus) {
  box-shadow: 0 0 0 1px #007bff inset;
}

:deep(.el-button--large) {
  padding: 12px 24px;
  font-size: 16px;
}
</style>