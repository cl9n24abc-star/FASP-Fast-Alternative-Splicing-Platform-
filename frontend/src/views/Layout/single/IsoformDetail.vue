<template>
  <el-dialog
    v-model="visible"
    title="🧬 Single Data Analysis - Isoform Details"
    width="85%"
    :before-close="handleClose"
    destroy-on-close
  >
    <div v-if="loading" class="loading-content">
      <el-skeleton :rows="5" animated />
    </div>

    <div v-else-if="isoformData" class="modal-content">
      
      <!-- Basic Information -->
      <div class="info-section">
        <h4>📋 Basic Information</h4>
        <div class="basic-info">
          <div class="info-item">
            <strong>ID:</strong>
            <span>{{ isoformData.id }}</span>
          </div>
          <div class="info-item">
            <strong>Gene Symbol:</strong>
            <span>{{ isoformData.geneSymbol }}</span>
          </div>
          <div class="info-item">
            <strong>Gene ID:</strong>
            <span>{{ isoformData.geneId }}</span>
          </div>
          <div class="info-item">
            <strong>Event Type:</strong>
            <el-tag :type="getEventTypeColor(isoformData.eventType)">
              {{ isoformData.eventType }}
            </el-tag>
          </div>
          <div class="info-item">
            <strong>Chromosome:</strong>
            <span>{{ isoformData.chromosome }}</span>
          </div>
          <div class="info-item">
            <strong>Strand:</strong>
            <span>{{ isoformData.strand }}</span>
          </div>
          <div class="info-item">
            <strong>Position:</strong>
            <span>{{ isoformData.exonStart?.toLocaleString() }} - {{ isoformData.exonEnd?.toLocaleString() }}</span>
          </div>
          <div class="info-item">
            <strong>Confidence:</strong>
            <span :class="getConfidenceClass(isoformData.confidence)">
              {{ (isoformData.confidence * 100).toFixed(1) }}%
            </span>
          </div>
        </div>
      </div>

      <!-- Sequencing Data -->
      <div class="reads-section">
        <h4>📊 Sequencing Data</h4>
        <div class="reads-grid">
          <div class="reads-card">
            <div class="card-header reads1">Sample 1 Reads</div>
            <div class="card-content">
              <div class="metric">
                <span class="label">Read Count:</span>
                <span class="value">{{ isoformData.reads1?.toLocaleString() }}</span>
              </div>
              <div class="metric">
                <span class="label">Percentage:</span>
                <span class="value">{{ ((isoformData.reads1 / (isoformData.reads1 + isoformData.reads2)) * 100).toFixed(1) }}%</span>
              </div>
            </div>
          </div>

          <div class="reads-card">
            <div class="card-header reads2">Sample 2 Reads</div>
            <div class="card-content">
              <div class="metric">
                <span class="label">Read Count:</span>
                <span class="value">{{ isoformData.reads2?.toLocaleString() }}</span>
              </div>
              <div class="metric">
                <span class="label">Percentage:</span>
                <span class="value">{{ ((isoformData.reads2 / (isoformData.reads1 + isoformData.reads2)) * 100).toFixed(1) }}%</span>
              </div>
            </div>
          </div>

          <div class="reads-card">
            <div class="card-header total">Total Reads</div>
            <div class="card-content">
              <div class="metric">
                <span class="label">Total Count:</span>
                <span class="value">{{ (isoformData.reads1 + isoformData.reads2).toLocaleString() }}</span>
              </div>
              <div class="metric">
                <span class="label">Average:</span>
                <span class="value">{{ ((isoformData.reads1 + isoformData.reads2) / 2).toLocaleString() }}</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- Statistical Analysis -->
      <div class="statistics-section">
        <h4>📈 Statistical Analysis</h4>
        <div class="statistics-content">
          <div class="statistics-grid">
            <div class="stat-item">
              <span class="stat-label">P-value:</span>
              <span 
                class="stat-value"
                :class="{ 'significant': isoformData.pValue < 0.05 }"
              >
                {{ isoformData.pValue?.toExponential(3) || 'N/A' }}
              </span>
            </div>
            <div class="stat-item">
              <span class="stat-label">FDR:</span>
              <span 
                class="stat-value"
                :class="{ 'significant': isoformData.fdr < 0.05 }"
              >
                {{ isoformData.fdr?.toExponential(3) || 'N/A' }}
              </span>
            </div>
            <div class="stat-item">
              <span class="stat-label">Significance:</span>
              <el-tag 
                :type="isoformData.isSignificant ? 'success' : 'info'"
                size="large"
              >
                {{ isoformData.isSignificant ? 'Significant' : 'Not Significant' }}
              </el-tag>
            </div>
            <div class="stat-item">
              <span class="stat-label">Original IncLevel:</span>
              <span class="stat-value">{{ isoformData.originalInc }}</span>
            </div>
            <div class="stat-item">
              <span class="stat-label">Compare IncLevel:</span>
              <span class="stat-value">{{ isoformData.compareInc }}</span>
            </div>
            <div class="stat-item">
              <span class="stat-label">IncLevel Difference:</span>
              <span 
                class="stat-value"
                :class="getIncLevelDiffClass(parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc))"
              >
                {{ (parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc)).toFixed(3) }}
              </span>
            </div>
          </div>
        </div>
      </div>

      <!-- Inclusion Level Visualization -->
      <div class="inclusion-section">
        <h4>🔄 Inclusion Level Analysis</h4>
        <div class="inclusion-content">
          <div class="inclusion-comparison">
            <div class="inclusion-item">
              <span class="inclusion-label">Original IncLevel:</span>
              <div class="inclusion-bar">
                <div 
                  class="inclusion-fill original" 
                  :style="{ width: (parseFloat(isoformData.originalInc) * 100) + '%' }"
                ></div>
                <span class="inclusion-value">{{ isoformData.originalInc }}</span>
              </div>
            </div>
            <div class="inclusion-item">
              <span class="inclusion-label">Compare IncLevel:</span>
              <div class="inclusion-bar">
                <div 
                  class="inclusion-fill compare"
                  :style="{ width: (parseFloat(isoformData.compareInc) * 100) + '%' }"
                ></div>
                <span class="inclusion-value">{{ isoformData.compareInc }}</span>
              </div>
            </div>
          </div>
          
          <div class="difference-summary">
            <div class="difference-value" :class="getIncLevelDiffClass(parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc))">
              Difference: {{ (parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc)) > 0 ? '+' : '' }}{{ (parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc)).toFixed(3) }}
            </div>
            <div class="difference-interpretation">
              <el-tag 
                :type="Math.abs(parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc)) > 0.1 ? 'warning' : 'info'"
                size="small"
              >
                {{ Math.abs(parseFloat(isoformData.compareInc) - parseFloat(isoformData.originalInc)) > 0.1 ? 'Notable Change' : 'Minor Change' }}
              </el-tag>
            </div>
          </div>
        </div>
      </div>

      <!-- Gene Structure Visualization -->
      <SplicingVisualization
        v-if="appStore.isGtfEnabled"
        :isoform-data="isoformData"
        :structure-data="structureData"
        @exon-click="handleExonClick"
      />

      <!-- Sashimi Plot Visualization -->
      <div v-if="appStore.isSashimiEnabled" class="sashimi-section">
        <h4>🍣 Sashimi Plot</h4>
        
        <!-- Sashimi控制面板 -->
        <div class="sashimi-controls">
          <el-button 
            type="primary" 
            @click="generateSashimiPlot"
            :loading="sashimiGenerating"
            :disabled="!canGenerateSashimi"
          >
            {{ sashimiGenerating ? 'Generating...' : 'Generate Sashimi Plot' }}
          </el-button>
          
          <el-button 
            v-if="sashimiResult"
            type="success" 
            @click="refreshSashimiPlot"
            :loading="sashimiGenerating"
          >
            Refresh Plot
          </el-button>
        </div>

        <!-- Sashimi内容 -->
        <div class="sashimi-content">
          <!-- 加载状态 -->
          <div v-if="sashimiGenerating" class="sashimi-loading">
            <el-skeleton :rows="3" animated />
            <div class="loading-text">
              <el-icon class="is-loading"><Loading /></el-icon>
              Generating SashimiPlot, please wait...
            </div>
          </div>
          
          <!-- 错误状态 -->
          <div v-else-if="sashimiError" class="sashimi-error">
            <el-alert 
              type="error" 
              :title="sashimiError" 
              show-icon 
              :closable="false"
            />
            <p>Please check if the data and file paths are correct</p>
          </div>
          
          <!-- 成功状态 - 显示生成的图片 -->
          <div v-else-if="sashimiResult" class="sashimi-success">
            <div class="sashimi-info">
              <el-tag type="success" size="large">
                <el-icon><Check /></el-icon>
                SashimiPlot generated successfully
              </el-tag>
              <span class="generation-time">
                Generation Time: {{ new Date(sashimiResult.timestamp).toLocaleString() }}
              </span>
            </div>
            
            <!-- 显示生成的图片文件 -->
            <div v-if="sashimiResult.generated_files && sashimiResult.generated_files.length > 0" class="sashimi-images">
              <div v-for="(file, index) in sashimiResult.generated_files" :key="index" class="sashimi-image-item">
                <h5>{{ getFileBaseName(file) }}</h5>
                    <div class="image-container">
                      <!-- 🔧 如果是图片文件，使用img标签 -->
                      <div v-if="isImageFile(file)" class="image-display">
                        <img 
                          :src="getSashimiFileUrl(file)" 
                          :alt="`SashimiPlot ${index + 1}`"
                          @error="handleImageError($event, file)"
                          @load="handleImageLoad($event, file)"
                        />
                      </div>
                      
                      <!-- 🆕 如果是PDF文件，使用iframe -->
                      <div v-else-if="isPdfFile(file)" class="pdf-display">
                        <iframe 
                          :src="getSashimiFileUrl(file)" 
                          width="100%" 
                          height="600px"
                          style="border: 1px solid #ddd; border-radius: 4px;"
                        >
                          Your browser does not support PDF preview
                        </iframe>
                        
                        <!-- PDF操作按钮 -->
                        <div class="pdf-actions">
                          <el-button type="primary" @click="openFileInNewTab(file)" size="small">
                            <el-icon><Document /></el-icon>
                            Open in New Tab
                          </el-button>
                          <el-button type="success" @click="downloadFile(file)" size="small">
                            <el-icon><Download /></el-icon>
                            Download PDF
                          </el-button>
                        </div>
                      </div>
                      
                      <!-- 🔧 其他文件类型 -->
                      <div v-else class="unsupported-file">
                        <el-icon size="48" color="#c0c4cc"><Document /></el-icon>
                        <p>Preview not supported for this file type:{{ getFileType(file) }}</p>
                        <el-button type="primary" @click="downloadFile(file)">
                          <el-icon><Download /></el-icon>
                          Download File
                        </el-button>
                      </div>
                    </div>
              </div>
            </div>
            
            <!-- 显示结果详情 -->
            <div class="sashimi-details">
              <el-collapse v-model="activeCollapse">
                <el-collapse-item title="SashimiPlot详细信息" name="details">
                  <div class="detail-grid">
                    <div class="detail-item">
                      <strong>Coordinates:</strong>
                      <span>{{ sashimiParams?.coordinate }}</span>
                    </div>
                    <div class="detail-item">
                      <strong>Sample Labels:</strong>
                      <span>{{ sashimiParams?.sample1_label }} vs {{ sashimiParams?.sample2_label }}</span>
                    </div>
                    <div class="detail-item">
                      <strong>Output Directory:</strong>
                      <span>{{ sashimiResult.output_directory }}</span>
                    </div>
                    <div class="detail-item">
                      <strong>Execution Time:</strong>
                      <span>{{ sashimiResult.duration_seconds?.toFixed(2) }}seconds</span>
                    </div>
                  </div>
                </el-collapse-item>
              </el-collapse>
            </div>
          </div>
          
          <!-- 默认状态 -->
          <div v-else class="sashimi-default">
            <div class="default-content">
              <el-icon size="48" color="#c0c4cc"><Picture /></el-icon>
              <p>Click the "Generate Sashimi Plot" button to generate a visualization</p>
              <p class="hint">SashimiPlot will display the alternative splicing patterns of this locus</p>
            </div>
          </div>
        </div>
        
      </div>

    </div>

    <template #footer>
      <div class="dialog-footer">
        <el-button @click="handleClose">Close</el-button>
      </div>
    </template>
  </el-dialog>

  <!-- Exon Details Dialog -->
  <el-dialog
    v-model="showExonDialog"
    title="📊 Exon Details"
    width="50%"
    :before-close="closeExonDialog"
  >
    <div v-if="selectedExonData" class="exon-detail-content">
      <div class="exon-info-grid">
        <div class="exon-info-item">
          <strong>Exon Name:</strong>
          <span>{{ selectedExonData.name }}</span>
        </div>
        <div class="exon-info-item">
          <strong>Start Position:</strong>
          <span>{{ selectedExonData.start?.toLocaleString() }}</span>
        </div>
        <div class="exon-info-item">
          <strong>End Position:</strong>
          <span>{{ selectedExonData.end?.toLocaleString() }}</span>
        </div>
        <div class="exon-info-item">
          <strong>Length:</strong>
          <span>{{ selectedExonData.length }} bp</span>
        </div>
        <div class="exon-info-item">
          <strong>Status:</strong>
          <el-tag :type="selectedExonData.isAffected ? 'warning' : 'success'">
            {{ selectedExonData.isAffected ? 'Affected by Alternative Splicing' : 'Constitutive' }}
          </el-tag>
        </div>
      </div>
    </div>
    <template #footer>
      <el-button @click="closeExonDialog">Close</el-button>
    </template>
  </el-dialog>
</template>

<script setup>
import { ref, computed, watch } from 'vue'
import { useAppStore } from '@/stores/app'
import { ElMessage, ElMessageBox } from 'element-plus'
import { Loading, Check, Picture, Document, Download } from '@element-plus/icons-vue'
import SplicingVisualization from './visualization.vue'
import SashimiPlot from '../../Components/SashimiPlot.vue'

// Props
const props = defineProps({
  modelValue: {
    type: Boolean,
    default: false
  },
  rowData: {
    type: Object,
    default: () => ({})
  }
})

// Emits
const emit = defineEmits(['update:modelValue'])

// Store
const appStore = useAppStore()

// Reactive Data
const loading = ref(false)
const isoformData = ref(null)
const structureData = ref(null)
const elementTooltip = ref({ show: false, content: '', x: 0, y: 0 })
const showExonDialog = ref(false)
const selectedExonData = ref(null)

// 🆕 SashimiPlot相关状态
const sashimiGenerating = ref(false)
const sashimiResult = ref(null)
const sashimiError = ref(null)
const sashimiParams = ref(null)
const activeCollapse = ref([])

// Computed Properties
const visible = computed({
  get() {
    return props.modelValue
  },
  set(value) {
    emit('update:modelValue', value)
  }
})

// 🆕 检查是否可以生成SashimiPlot
const canGenerateSashimi = computed(() => {
  if (!isoformData.value) return false
  
  const filePaths = appStore.getFilePaths
  const hasRequiredData = !!(
    isoformData.value.chromosome && 
    isoformData.value.strand && 
    isoformData.value.exonStart && 
    isoformData.value.exonEnd
  )
  
  const hasRequiredFiles = !!(
    filePaths.sashimiBam1 && 
    filePaths.sashimiBam2 && 
    filePaths.sashimiGff
  )
  
  return hasRequiredData && hasRequiredFiles
})

// 🆕 调用SashimiPlot API
const callSashimiAPI = async () => {
  try {
    const { chromosome, strand, exonStart, exonEnd } = isoformData.value
    const filePaths = appStore.getFilePaths
    
    // 构建API参数
    const params = {
      coordinate: `${chromosome}:${strand}:${exonStart}:${exonEnd}:${filePaths.sashimiGff}`,
      sample1_label: 'Sample1',
      sample2_label: 'Sample2',
      sample1_bam: filePaths.sashimiBam1,
      sample2_bam: filePaths.sashimiBam2,
      min_counts: 10,
      figure_height: 7,
      figure_width: 8
    }
    
    // 保存参数用于显示
    sashimiParams.value = params
    
    console.log('🍣 调用SashimiPlot API:', params)
    
    // 调用API
    const response = await fetch('/api/sashimi/generate-by-coordinate', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(params)
    })
    
    const result = await response.json()
    
    if (result.success) {
      console.log('✅ SashimiPlot生成成功:', result)
      
      // 添加时间戳
      result.timestamp = new Date().toISOString()
      
      // 保存结果到store
      appStore.addSashimiResult({
        isoformId: isoformData.value.id,
        coordinate: params.coordinate,
        outputPath: result.output_directory,
        generatedFiles: result.generated_files,
        timestamp: result.timestamp
      })
      
      return result
    } else {
      console.error('❌ SashimiPlot生成失败:', result.error)
      throw new Error(result.error || 'SashimiPlot生成失败')
    }
    
  } catch (error) {
    console.error('🚫 调用SashimiPlot API失败:', error)
    throw error
  }
}

// 🆕 生成SashimiPlot
const generateSashimiPlot = async () => {
  if (!canGenerateSashimi.value) {
    ElMessage.warning('数据或文件配置不完整，无法生成SashimiPlot')
    return
  }
  
  sashimiGenerating.value = true
  sashimiError.value = null
  sashimiResult.value = null
  
  try {
    const result = await callSashimiAPI()
    sashimiResult.value = result
    ElMessage.success('SashimiPlot生成成功！')
  } catch (error) {
    sashimiError.value = error.message
    ElMessage.error(`SashimiPlot生成失败: ${error.message}`)
  } finally {
    sashimiGenerating.value = false
  }
}

// 🆕 刷新SashimiPlot
const refreshSashimiPlot = async () => {
  await generateSashimiPlot()
}

// 🆕 获取图片URL
const getSashimiFileUrl = (filePath) => {
  const fileName = filePath.split('/').pop()
  return `/sashimi_output/Sashimi_plot/${fileName}`
}

const getFileType = (filePath) => {
  const fileName = filePath.split('/').pop()
  const extension = fileName.split('.').pop().toLowerCase()
  return extension
}

const isImageFile = (filePath) => {
  const imageExtensions = ['png', 'jpg', 'jpeg', 'gif', 'svg', 'webp']
  const extension = getFileType(filePath)
  return imageExtensions.includes(extension)
}

const isPdfFile = (filePath) => {
  const extension = getFileType(filePath)
  return extension === 'pdf'
}

const downloadFile = (filePath) => {
  const fileUrl = getSashimiFileUrl(filePath)
  const fileName = getFileBaseName(filePath)
  
  const link = document.createElement('a')
  link.href = fileUrl
  link.download = fileName
  link.target = '_blank'
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
}

const openFileInNewTab = (filePath) => {
  const fileUrl = getSashimiFileUrl(filePath)
  window.open(fileUrl, '_blank')
}

// 🆕 获取文件基础名称
const getFileBaseName = (filePath) => {
  return filePath.split('/').pop()
}

// 🆕 图片加载错误处理
const handleImageError = (event, filePath) => {
  console.error('图片加载失败:', filePath)
  event.target.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMjAwIiBoZWlnaHQ9IjEwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48cmVjdCB3aWR0aD0iMTAwJSIgaGVpZ2h0PSIxMDAlIiBmaWxsPSIjZjVmNWY1Ii8+PHRleHQgeD0iNTAlIiB5PSI1MCUiIGZvbnQtZmFtaWx5PSJBcmlhbCwgc2Fucy1zZXJpZiIgZm9udC1zaXplPSIxMiIgZmlsbD0iIzk5OSIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZHk9Ii4zZW0iPkltYWdlIE5vdCBGb3VuZDwvdGV4dD48L3N2Zz4='
}

// 🆕 图片加载成功处理
const handleImageLoad = (event, filePath) => {
  console.log('图片加载成功:', filePath)
}

// Generate Mock Detail Data
const generateDetailData = (rowData) => {
  if (!rowData || !rowData.id) return null

  const generateGeneStructure = () => {
    const eventTypes = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
    const eventType = rowData.eventType || eventTypes[Math.floor(Math.random() * eventTypes.length)]
    const exonCount = Math.floor(Math.random() * 8) + 4
    
    const structure = []
    const transcript = []
    const exonDetails = []
    
    for (let i = 1; i <= exonCount; i++) {
      const exonLength = Math.floor(Math.random() * 200) + 100
      const intronLength = Math.floor(Math.random() * 1000) + 500
      
      const isAffectedExon = (eventType === 'SE' && i === Math.floor(exonCount / 2)) ||
                            (eventType === 'MXE' && (i === Math.floor(exonCount / 2) || i === Math.floor(exonCount / 2) + 1))
      
      structure.push({
        type: 'exon',
        name: `E${i}`,
        width: 8,
        isAffected: isAffectedExon,
        length: exonLength,
        start: 1000000 + (i - 1) * 1500,
        end: 1000000 + (i - 1) * 1500 + exonLength
      })
      
      transcript.push({
        type: 'exon',
        name: `E${i}`,
        width: 8,
        length: exonLength,
        isAffected: isAffectedExon
      })
      
      exonDetails.push({
        name: `Exon ${i}`,
        start: 1000000 + (i - 1) * 1500,
        end: 1000000 + (i - 1) * 1500 + exonLength,
        length: exonLength,
        inclusionLevel: isAffectedExon ? Math.random() * 0.5 + 0.2 : Math.random() * 0.3 + 0.7,
        status: isAffectedExon ? 'Alternative' : 'Constitutive'
      })
      
      if (i < exonCount) {
        structure.push({
          type: 'intron',
          name: `I${i}`,
          width: 4,
          length: intronLength
        })
        
        transcript.push({
          type: 'intron',
          name: `I${i}`,
          width: 4
        })
      }
    }
    
    return {
      geneStart: 1000000,
      geneEnd: 1000000 + exonCount * 1500,
      strand: rowData.strand || '+',
      splicingEvents: [eventType],
      structure,
      transcript,
      exonDetails
    }
  }
  
  return {
    structure: generateGeneStructure()
  }
}

// Load detailed data
const loadDetailData = async () => {
  if (!props.rowData || !props.rowData.id) return

  loading.value = true
  
  try {
    isoformData.value = props.rowData
    
    if (props.rowData.structureData) {
      structureData.value = props.rowData.structureData
      
      if (!structureData.value.transcript) {
        if (structureData.value.originalTranscript) {
          structureData.value.transcript = structureData.value.originalTranscript
        } else if (structureData.value.structure) {
          structureData.value.transcript = structureData.value.structure.filter(item => item.type === 'exon')
        } else {
          structureData.value.transcript = [{
            type: 'exon',
            name: 'E1',
            width: 8,
            isAffected: true,
            length: structureData.value.geneEnd - structureData.value.geneStart,
            start: structureData.value.geneStart,
            end: structureData.value.geneEnd
          }]
        }
      }

      if (!structureData.value.exonDetails) {
        structureData.value.exonDetails = []
        
        if (structureData.value.originalTranscript) {
          structureData.value.exonDetails = structureData.value.originalTranscript
            .filter(item => item.type === 'exon')
            .map((exon, index) => ({
              name: exon.name || `Exon ${index + 1}`,
              start: exon.start,
              end: exon.end,
              length: exon.length,
              inclusionLevel: isoformData.value.originalInc || Math.random() * 0.5 + 0.5,
              status: exon.isAffected ? 'Alternative' : 'Constitutive'
            }))
        }
      }
      
      console.log('✅ Using REAL structure data (Fixed)')
    } else {
      const detailData = generateDetailData(props.rowData)
      structureData.value = detailData.structure
      console.log('🔄 Using MOCK structure data')
    }
    
  } catch (error) {
    console.error('Failed to load detailed data:', error)
    const detailData = generateDetailData(props.rowData)
    structureData.value = detailData.structure
  } finally {
    loading.value = false
  }
}

// Methods
const handleClose = () => {
  visible.value = false
}

const exportData = () => {
  console.log('Exporting data:', {
    isoform: isoformData.value,
    structure: structureData.value,
    sashimi: sashimiResult.value
  })
}

const getEventTypeColor = (type) => {
  const colors = {
    'SE': 'primary',
    'MXE': 'success', 
    'A5SS': 'warning',
    'A3SS': 'danger',
    'RI': 'info'
  }
  return colors[type] || 'default'
}

const getConfidenceClass = (value) => {
  if (value >= 0.8) return 'high-confidence'
  if (value >= 0.6) return 'medium-confidence'
  return 'low-confidence'
}

const getIncLevelDiffClass = (value) => {
  if (Math.abs(value) > 0.2) return 'high-change'
  if (Math.abs(value) > 0.1) return 'medium-change'
  return 'low-change'
}

const handleExonClick = (exon) => {
  selectedExonData.value = exon
  showExonDialog.value = true
}

const closeExonDialog = () => {
  showExonDialog.value = false
  selectedExonData.value = null
}

// Watch for modal opening to load data
watch(visible, (newVal) => {
  if (newVal) {
    loadDetailData()
    
    // 🆕 重置SashimiPlot状态
    sashimiResult.value = null
    sashimiError.value = null
    sashimiGenerating.value = false
    
    // 🆕 检查是否已有缓存的SashimiPlot结果
    if (isoformData.value?.id) {
      const existingResult = appStore.getSashimiResults.find(
        result => result.isoformId === isoformData.value.id
      )
      
      if (existingResult) {
        console.log('📄 发现缓存的SashimiPlot结果')
        // 这里可以选择是否自动使用缓存结果
        // sashimiResult.value = existingResult
      }
    }
  }
})

</script>

<style scoped>
/* 🚀 主要布局样式 */
.modal-content {
  height: 75vh !important;
  overflow-y: auto !important;
  padding-bottom: 80px !important;
  box-sizing: border-box !important;
}

.loading-content {
  padding: 20px;
}

/* 🎯 SashimiPlot 专用样式 */
.sashimi-section {
  margin-bottom: 24px;
  padding: 20px !important;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #e9ecef;
  min-height: 500px !important;
  height: auto !important;
  overflow: visible !important;
}

.sashimi-section h4 {
  margin-top: 0;
  margin-bottom: 12px;
  color: #495057;
  border-bottom: 2px solid #dee2e6;
  padding-bottom: 6px;
}

/* 🆕 SashimiPlot控制面板样式 */
.sashimi-controls {
  margin-bottom: 20px;
  display: flex;
  gap: 12px;
  align-items: center;
}

/* 🆕 SashimiPlot内容区域样式 */
.sashimi-content {
  background-color: white;
  border-radius: 6px;
  border: 1px solid #e9ecef;
  padding: 20px;
  min-height: 400px;
}

/* 🆕 加载状态样式 */
.sashimi-loading {
  text-align: center;
  padding: 40px 20px;
}

.loading-text {
  margin-top: 16px;
  color: #6c757d;
  font-size: 14px;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 8px;
}

/* 🆕 错误状态样式 */
.sashimi-error {
  text-align: center;
  padding: 40px 20px;
}

.sashimi-error p {
  margin-top: 16px;
  color: #6c757d;
  font-size: 14px;
}

/* 🆕 成功状态样式 */
.sashimi-success {
  padding: 0;
}

.sashimi-info {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
  padding-bottom: 12px;
  border-bottom: 1px solid #e9ecef;
}

.generation-time {
  font-size: 12px;
  color: #6c757d;
}

/* 🆕 图片显示样式 */
.sashimi-images {
  margin-bottom: 20px;
}

.sashimi-image-item {
  margin-bottom: 20px;
  border: 1px solid #e9ecef;
  border-radius: 6px;
  overflow: hidden;
}

.sashimi-image-item h5 {
  margin: 0;
  padding: 12px 16px;
  background-color: #f8f9fa;
  border-bottom: 1px solid #e9ecef;
  font-size: 14px;
  font-weight: 600;
  color: #495057;
}

.image-container {
  padding: 16px;
  text-align: center;
  background-color: white;
}

.image-container img {
  max-width: 100%;
  height: auto;
  border-radius: 4px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

/* 🆕 详情折叠面板样式 */
.sashimi-details {
  margin-top: 20px;
}

.detail-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: 12px;
}

.detail-item {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 8px 0;
  border-bottom: 1px solid #f1f3f4;
}

.detail-item:last-child {
  border-bottom: none;
}

.detail-item strong {
  color: #495057;
  font-weight: 600;
}

.detail-item span {
  color: #6c757d;
  font-family: monospace;
  font-size: 12px;
  word-break: break-all;
}

/* 🆕 默认状态样式 */
.sashimi-default {
  text-align: center;
  padding: 60px 20px;
}

.default-content {
  color: #6c757d;
}

.default-content p {
  margin: 16px 0 8px 0;
  font-size: 16px;
}

.default-content .hint {
  font-size: 14px;
  color: #adb5bd;
}

/* 📦 通用 section 样式 */
.info-section,
.reads-section,
.statistics-section,
.inclusion-section {
  margin-bottom: 24px;
  padding: 16px;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #e9ecef;
}

.info-section h4,
.reads-section h4,
.statistics-section h4,
.inclusion-section h4 {
  margin-top: 0;
  margin-bottom: 12px;
  color: #495057;
  border-bottom: 2px solid #dee2e6;
  padding-bottom: 6px;
}

/* 📋 基本信息样式 */
.basic-info {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 12px;
}

.info-item {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

/* 📊 读取数据样式 */
.reads-grid {
  display: grid;
  grid-template-columns: 1fr 1fr 1fr;
  gap: 16px;
}

.reads-card {
  border: 1px solid #dee2e6;
  border-radius: 6px;
  overflow: hidden;
  background-color: white;
}

.card-header {
  padding: 10px;
  font-weight: bold;
  text-align: center;
  color: white;
}

.card-header.reads1 {
  background-color: #007bff;
}

.card-header.reads2 {
  background-color: #28a745;
}

.card-header.total {
  background-color: #6f42c1;
}

.card-content {
  padding: 12px;
}

.metric {
  display: flex;
  justify-content: space-between;
  margin-bottom: 8px;
}

.metric:last-child {
  margin-bottom: 0;
}

.metric .label {
  color: #6c757d;
  font-size: 14px;
}

.metric .value {
  font-weight: bold;
  color: #495057;
}

/* 📈 统计分析样式 */
.statistics-content {
  background-color: white;
  padding: 16px;
  border-radius: 4px;
}

.statistics-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 12px;
}

.stat-item {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 8px 0;
  border-bottom: 1px solid #f1f3f4;
}

.stat-label {
  color: #6c757d;
  font-weight: 500;
}

.stat-value {
  font-weight: bold;
}

.stat-value.significant {
  color: #dc3545;
}

/* 🎨 状态样式 */
.high-confidence {
  color: #28a745;
  font-weight: bold;
}

.medium-confidence {
  color: #fd7e14;
  font-weight: 500;
}

.low-confidence {
  color: #6c757d;
}

.high-change {
  color: #dc3545;
}

.medium-change {
  color: #fd7e14;
}

.low-change {
  color: #6c757d;
}

/* 🔄 包含水平样式 */
.inclusion-content {
  background-color: white;
  padding: 16px;
  border-radius: 4px;
}

.inclusion-comparison {
  margin-bottom: 20px;
}

.inclusion-item {
  margin-bottom: 16px;
}

.inclusion-item:last-child {
  margin-bottom: 0;
}

.inclusion-label {
  display: block;
  margin-bottom: 6px;
  font-weight: 500;
  color: #495057;
  font-size: 14px;
}

.inclusion-bar {
  position: relative;
  height: 24px;
  background-color: #f8f9fa;
  border-radius: 12px;
  border: 1px solid #dee2e6;
  overflow: hidden;
}

.inclusion-fill {
  height: 100%;
  border-radius: 12px;
  transition: width 0.3s ease;
}

.inclusion-fill.original {
  background: linear-gradient(90deg, #28a745, #20c997);
}

.inclusion-fill.compare {
  background: linear-gradient(90deg, #007bff, #6f42c1);
}

.inclusion-value {
  position: absolute;
  top: 50%;
  right: 8px;
  transform: translateY(-50%);
  font-size: 11px;
  font-weight: bold;
  color: #495057;
}

.difference-summary {
  text-align: center;
  padding: 16px;
  background-color: #f8f9fa;
  border-radius: 4px;
  border: 1px solid #e9ecef;
}

.difference-value {
  font-size: 1.25rem;
  font-weight: bold;
  margin-bottom: 8px;
}

/* 🔧 对话框样式覆盖 */
:deep(.el-dialog) {
  margin-top: 3vh !important;
  margin-bottom: 3vh !important;
  max-height: 94vh !important;
}

:deep(.el-dialog__body) {
  padding: 10px 20px !important;
  max-height: none !important;
  overflow: visible !important;
}

/* 🗂️ 其他组件样式 */
.dialog-footer {
  text-align: right;
}

.exon-detail-content {
  padding: 16px;
}

.exon-info-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: 12px;
}

.exon-info-item {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 8px 0;
  border-bottom: 1px solid #f1f3f4;
}
.pdf-actions {
  margin-top: 12px;
  text-align: center;
  display: flex;
  gap: 12px;
  justify-content: center;
}

/* 🔧 图片显示样式 */
.image-display img {
  max-width: 100%;
  height: auto;
  border-radius: 4px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

/* 🔧 不支持的文件类型样式 */
.unsupported-file {
  text-align: center;
  padding: 40px 20px;
  color: #6c757d;
}

.unsupported-file p {
  margin: 16px 0;
  font-size: 14px;
}
</style>