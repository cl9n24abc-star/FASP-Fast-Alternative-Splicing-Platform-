<template>
    <div class="sashimi-plot">
      <!-- Header -->
      <div class="plot-header">
        <h5>🍣 Sashimi Plot Visualization (Test Mode)</h5>
        <div class="plot-info">
          <span v-if="isoformData">Gene: {{ isoformData.geneSymbol || 'TEST_GENE' }}</span>
          <span v-if="isoformData">Event: {{ isoformData.eventType || 'SE' }}</span>
          <span :class="apiStatusClass">API: {{ apiStatus }}</span>
        </div>
      </div>

      <!-- API Test Controls -->
      <div class="test-controls">
        <div class="control-group">
          <el-button 
            type="primary" 
            @click="testSashimiAPI" 
            :loading="isGenerating"
            size="small"
          >
            {{ isGenerating ? 'Generating...' : '🧪 Test Sashimi API' }}
          </el-button>
          
          <el-button 
            type="info" 
            @click="checkSashimiTool" 
            :loading="isChecking"
            size="small"
          >
            {{ isChecking ? 'Checking...' : '🔍 Check Tool' }}
          </el-button>
          
          <el-button 
            type="warning" 
            @click="toggleDebug" 
            size="small"
          >
            {{ showDebug ? 'Hide' : 'Show' }} Debug
          </el-button>
        </div>
        
        <!-- API Response Display -->
        <div v-if="apiResponse" class="api-response">
          <el-alert
            :title="apiResponse.success ? '✅ API Success' : '❌ API Error'"
            :type="apiResponse.success ? 'success' : 'error'"
            :closable="true"
            show-icon
            @close="apiResponse = null"
          >
            <template #default>
              <div class="response-content">
                <p><strong>Message:</strong> {{ apiResponse.message }}</p>
                <div v-if="apiResponse.result" class="result-details">
                  <p><strong>Output Path:</strong> {{ apiResponse.result.outputPath }}</p>
                  <p v-if="apiResponse.result.plotFiles"><strong>Generated Files:</strong></p>
                  <ul v-if="apiResponse.result.plotFiles">
                    <li v-for="file in apiResponse.result.plotFiles" :key="file">{{ file }}</li>
                  </ul>
                </div>
                <div v-if="apiResponse.error" class="error-details">
                  <p><strong>Error:</strong> {{ apiResponse.error }}</p>
                </div>
              </div>
            </template>
          </el-alert>
        </div>
      </div>

      <!-- Main Plot Area -->
      <div class="plot-container">
        <!-- Coverage Tracks -->
        <div class="coverage-tracks">
          <div class="track">
            <div class="track-label">Sample 1</div>
            <div class="coverage-area">
              <svg width="100%" height="60" class="coverage-svg">
                <!-- Sample coverage bars -->
                <rect v-for="(bar, index) in sample1Coverage" 
                      :key="`s1-${index}`"
                      :x="index * 8" 
                      :y="60 - bar.height" 
                      :width="6" 
                      :height="bar.height"
                      fill="#007bff"
                      opacity="0.7">
                </rect>
              </svg>
            </div>
          </div>

          <div class="track">
            <div class="track-label">Sample 2</div>
            <div class="coverage-area">
              <svg width="100%" height="60" class="coverage-svg">
                <!-- Sample coverage bars -->
                <rect v-for="(bar, index) in sample2Coverage" 
                      :key="`s2-${index}`"
                      :x="index * 8" 
                      :y="60 - bar.height" 
                      :width="6" 
                      :height="bar.height"
                      fill="#28a745"
                      opacity="0.7">
                </rect>
              </svg>
            </div>
          </div>
        </div>

        <!-- Gene Structure -->
        <div class="gene-structure">
          <div class="structure-label">Gene Structure</div>
          <div class="structure-track">
            <svg width="100%" height="40" class="structure-svg">
              <!-- Exons -->
              <rect v-for="(exon, index) in mockExons"
                    :key="`exon-${index}`"
                    :x="exon.x"
                    :y="15"
                    :width="exon.width"
                    :height="10"
                    :fill="exon.isAffected ? '#ff6b6b' : '#4ecdc4'"
                    stroke="#2c3e50"
                    stroke-width="1">
              </rect>
              
              <!-- Introns (connecting lines) -->
              <line v-for="(intron, index) in mockIntrons"
                    :key="`intron-${index}`"
                    :x1="intron.x1"
                    :y1="20"
                    :x2="intron.x2"
                    :y2="20"
                    stroke="#95a5a6"
                    stroke-width="2">
              </line>
            </svg>
          </div>
        </div>

        <!-- Junction Arcs -->
        <div class="junction-arcs">
          <div class="arc-label">Splice Junctions</div>
          <svg width="100%" height="80" class="arc-svg">
            <!-- Sample junction arcs -->
            <path v-for="(arc, index) in mockJunctions"
                  :key="`arc-${index}`"
                  :d="arc.path"
                  :stroke="arc.color"
                  :stroke-width="arc.width"
                  fill="none"
                  opacity="0.8">
            </path>
            
            <!-- Junction labels -->
            <text v-for="(arc, index) in mockJunctions"
                  :key="`label-${index}`"
                  :x="arc.labelX"
                  :y="arc.labelY"
                  text-anchor="middle"
                  font-size="10"
                  fill="#666">
              {{ arc.reads }}
            </text>
          </svg>
        </div>
      </div>

      <!-- Legend -->
      <div class="plot-legend">
        <div class="legend-item">
          <div class="legend-color sample1"></div>
          <span>Sample 1 Coverage</span>
        </div>
        <div class="legend-item">
          <div class="legend-color sample2"></div>
          <span>Sample 2 Coverage</span>
        </div>
        <div class="legend-item">
          <div class="legend-color exon"></div>
          <span>Exon</span>
        </div>
        <div class="legend-item">
          <div class="legend-color affected-exon"></div>
          <span>Affected Exon</span>
        </div>
      </div>

      <!-- Debug Info -->
      <div class="debug-info" v-if="showDebug">
        <details>
          <summary>🔧 Debug Information</summary>
          <div class="debug-content">
            <h6>Test Parameters:</h6>
            <pre>{{ JSON.stringify(testParams, null, 2) }}</pre>
            
            <h6>API Response:</h6>
            <pre>{{ JSON.stringify(apiResponse, null, 2) }}</pre>
            
            <h6>Component Data:</h6>
            <pre>{{ JSON.stringify(debugInfo, null, 2) }}</pre>
          </div>
        </details>
      </div>
    </div>
</template>

<script setup>
import { ref, computed, onMounted } from 'vue'

// Props
const props = defineProps({
  isoformData: {
    type: Object,
    default: () => ({})
  },
  structureData: {
    type: Object,
    default: () => ({})
  }
})

// Reactive data
const showDebug = ref(true) // 默认显示调试信息
const isGenerating = ref(false)
const isChecking = ref(false)
const apiStatus = ref('未测试')
const apiResponse = ref(null)

// Mock data for visualization
const sample1Coverage = ref([])
const sample2Coverage = ref([])
const mockExons = ref([])
const mockIntrons = ref([])
const mockJunctions = ref([])

// 修正的测试参数 - 使用rmats2sashimiplot的正确格式
const testParams = computed(() => {
  // 使用测试数据，模拟实际的参数
  const params = {
    // rmats2sashimiplot 使用的坐标格式：chromosome:strand:start:end:gff_file
    coordinate: "21:+:46112804:46112824:/tmp/gencode_no_chr.gff3",
    
    // 输出目录
    output_directory: '/home/cl9n24abc/project/output',
    
    // 样本标签
    sample1_label: 'SRR12125133',
    sample2_label: 'SRR12125135', 
    
    // BAM文件路径
    sample1_bam: '/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam',
    sample2_bam: '/mnt/disk1/rna_seq_analysis/step0/SRR12125135_Aligned.sortedByCoord.out.bam',
    
    // 可选参数
    min_counts: 10,
    figure_height: 7,
    figure_width: 8,
    
    // 添加调试信息
    debug_info: {
      original_command: 'rmats2sashimiplot -o /tmp/test -c "21:+:46112804:46112824:/tmp/gencode_no_chr.gff3" --l1 SRR12125133 --l2 SRR12125135 --b1 /path/to/bam1 --b2 /path/to/bam2',
      note: '这是根据成功命令改写的参数格式'
    }
  }
  
  console.log('🔍 构建rmats2sashimiplot参数:', params)
  return params
})

// API状态样式
const apiStatusClass = computed(() => {
  switch (apiStatus.value) {
    case '成功': return 'status-success'
    case '失败': return 'status-error'
    case '生成中':
    case '检查中': return 'status-loading'
    default: return 'status-default'
  }
})

// 测试 SashimiPlot API - 使用rmats2sashimiplot格式
const testSashimiAPI = async () => {
  isGenerating.value = true
  apiStatus.value = '生成中'
  apiResponse.value = null
  
  try {
    console.log('🧪 开始测试rmats2sashimiplot API')
    console.log('📤 发送参数:', testParams.value)
    
    // 验证关键参数
    if (!testParams.value.coordinate) {
      console.error('❌ 缺少coordinate参数')
      throw new Error('缺少coordinate参数')
    }
    
    // 修改API端点为rmats2sashimiplot专用
    const response = await fetch('/api/sashimi/generate-rmats', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify(testParams.value)
    })
    
    const result = await response.json()
    
    console.log('📥 API 响应:', result)
    
    if (response.ok && result.success) {
      apiStatus.value = '成功'
      apiResponse.value = {
        success: true,
        message: result.message || 'rmats2sashimiplot 生成成功',
        result: result
      }
      console.log('✅ rmats2sashimiplot 生成成功!')
    } else {
      apiStatus.value = '失败'
      apiResponse.value = {
        success: false,
        message: result.message || 'API 调用失败',
        error: result.error || `HTTP ${response.status}`
      }
      console.error('❌ rmats2sashimiplot 生成失败:', result)
    }
    
  } catch (error) {
    apiStatus.value = '失败'
    apiResponse.value = {
      success: false,
      message: error.message || '网络错误或API不可用',
      error: error.message
    }
    console.error('❌ API 调用出错:', error)
  } finally {
    isGenerating.value = false
  }
}

// 检查 rmats2sashimiplot 工具
const checkSashimiTool = async () => {
  isChecking.value = true
  apiStatus.value = '检查中'
  
  try {
    console.log('🔍 检查 rmats2sashimiplot 工具可用性')
    
    const response = await fetch('/api/sashimi/check-rmats')
    const result = await response.json()
    
    console.log('🔧 工具检查结果:', result)
    
    if (result.available) {
      apiStatus.value = '工具可用'
      apiResponse.value = {
        success: true,
        message: `rmats2sashimiplot 可用: ${result.message}`,
        result: result
      }
    } else {
      apiStatus.value = '工具不可用'
      apiResponse.value = {
        success: false,
        message: `rmats2sashimiplot 不可用: ${result.message}`,
        error: result.error
      }
    }
    
  } catch (error) {
    apiStatus.value = '检查失败'
    apiResponse.value = {
      success: false,
      message: '无法检查rmats2sashimiplot状态',
      error: error.message
    }
    console.error('❌ 工具检查失败:', error)
  } finally {
    isChecking.value = false
  }
}

// 切换调试模式
const toggleDebug = () => {
  showDebug.value = !showDebug.value
}

// Generate mock data
const generateMockData = () => {
  // Generate coverage data (simulating read coverage)
  sample1Coverage.value = Array.from({ length: 50 }, (_, i) => ({
    height: Math.random() * 40 + 10
  }))
  
  sample2Coverage.value = Array.from({ length: 50 }, (_, i) => ({
    height: Math.random() * 35 + 15
  }))

  // Generate mock exons
  mockExons.value = [
    { x: 20, width: 60, isAffected: false },
    { x: 120, width: 40, isAffected: true }, // affected exon
    { x: 200, width: 50, isAffected: false },
    { x: 300, width: 45, isAffected: false }
  ]

  // Generate introns (connecting lines between exons)
  mockIntrons.value = [
    { x1: 80, x2: 120 },
    { x1: 160, x2: 200 },
    { x1: 250, x2: 300 }
  ]

  // Generate splice junction arcs
  mockJunctions.value = [
    {
      path: `M 80 70 Q 100 30 120 70`,
      color: '#007bff',
      width: 3,
      reads: '45',
      labelX: 100,
      labelY: 25
    },
    {
      path: `M 160 70 Q 180 35 200 70`,
      color: '#28a745',
      width: 2,
      reads: '23',
      labelX: 180,
      labelY: 30
    },
    {
      path: `M 250 70 Q 275 40 300 70`,
      color: '#dc3545',
      width: 4,
      reads: '67',
      labelX: 275,
      labelY: 35
    }
  ]
}

// Debug information
const debugInfo = computed(() => {
  return {
    isoformData: props.isoformData,
    structureData: props.structureData,
    hasData: !!(props.isoformData && props.structureData),
    timestamp: new Date().toISOString(),
    apiStatus: apiStatus.value,
    toolInfo: {
      name: 'rmats2sashimiplot',
      commandFormat: 'rmats2sashimiplot -o OUTPUT_DIR -c "chr:strand:start:end:gff" --l1 LABEL1 --l2 LABEL2 --b1 BAM1 --b2 BAM2',
      successfulExample: 'rmats2sashimiplot -o /tmp/test -c "21:+:46112804:46112824:/tmp/gencode_no_chr.gff3" --l1 SRR12125133 --l2 SRR12125135'
    },
    apiEndpoints: {
      generate: '/api/sashimi/generate-rmats',
      check: '/api/sashimi/check-rmats'
    }
  }
})

// Lifecycle
onMounted(() => {
  generateMockData()
  console.log('🍣 rmats2sashimiplot 测试组件已挂载!')
  console.log('📊 测试参数:', testParams.value)
  console.log('🔧 使用rmats2sashimiplot格式')
  
  // 自动检查工具可用性
  setTimeout(() => {
    checkSashimiTool()
  }, 1000)
})
</script>

<style scoped>
.sashimi-plot {
  background-color: white;
  border-radius: 8px;
  padding: 20px;
  border: 1px solid #e9ecef;
  min-height: 400px !important;
  height: auto !important;
  overflow: visible !important;
  display: flex;
  flex-direction: column;
}

.plot-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 16px;
  padding-bottom: 8px;
  border-bottom: 1px solid #dee2e6;
}

.plot-header h5 {
  margin: 0;
  color: #2c3e50;
  font-size: 16px;
}

.plot-info {
  display: flex;
  gap: 16px;
  font-size: 12px;
  color: #6c757d;
}

/* 测试控制面板样式 */
.test-controls {
  margin-bottom: 20px;
  padding: 16px;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #e9ecef;
}

.control-group {
  display: flex;
  gap: 12px;
  margin-bottom: 16px;
  flex-wrap: wrap;
}

.api-response {
  margin-top: 12px;
}

.response-content {
  font-size: 14px;
}

.result-details,
.error-details {
  margin-top: 8px;
  padding: 8px;
  background-color: rgba(255, 255, 255, 0.5);
  border-radius: 4px;
}

.result-details ul {
  margin: 4px 0;
  padding-left: 20px;
}

/* API 状态样式 */
.status-success {
  color: #28a745;
  font-weight: bold;
}

.status-error {
  color: #dc3545;
  font-weight: bold;
}

.status-loading {
  color: #007bff;
  font-weight: bold;
}

.status-default {
  color: #6c757d;
}

/* 原有样式保持不变 */
.plot-container {
  margin-bottom: 20px;
  flex: 1;
  min-height: 300px;
}

.coverage-tracks {
  margin-bottom: 16px;
}

.track {
  display: flex;
  align-items: center;
  margin-bottom: 8px;
}

.track-label {
  width: 80px;
  font-size: 12px;
  font-weight: 500;
  color: #495057;
  text-align: right;
  margin-right: 12px;
}

.coverage-area {
  flex: 1;
  border: 1px solid #dee2e6;
  border-radius: 4px;
  background-color: #f8f9fa;
}

.coverage-svg {
  display: block;
}

.gene-structure {
  margin-bottom: 16px;
}

.structure-label {
  font-size: 12px;
  font-weight: 500;
  color: #495057;
  margin-bottom: 4px;
}

.structure-track {
  border: 1px solid #dee2e6;
  border-radius: 4px;
  background-color: #f8f9fa;
}

.structure-svg {
  display: block;
}

.junction-arcs {
  margin-bottom: 16px;
}

.arc-label {
  font-size: 12px;
  font-weight: 500;
  color: #495057;
  margin-bottom: 4px;
}

.arc-svg {
  display: block;
  border: 1px solid #dee2e6;
  border-radius: 4px;
  background-color: #f8f9fa;
}

.plot-legend {
  display: flex;
  justify-content: center;
  gap: 20px;
  flex-wrap: wrap;
  padding: 16px !important;
  margin: 16px 0 20px 0 !important;
  background-color: #f8f9fa;
  border-radius: 4px;
  border: 1px solid #e9ecef;
  min-height: 60px !important;
  box-sizing: border-box !important;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  color: #495057;
}

.legend-color {
  width: 16px;
  height: 12px;
  border-radius: 2px;
  border: 1px solid #dee2e6;
}

.legend-color.sample1 {
  background-color: #007bff;
}

.legend-color.sample2 {
  background-color: #28a745;
}

.legend-color.exon {
  background-color: #4ecdc4;
}

.legend-color.affected-exon {
  background-color: #ff6b6b;
}

.debug-info {
  margin-top: 16px;
  padding: 12px;
  background-color: #f8f9fa;
  border-radius: 4px;
  border: 1px solid #e9ecef;
}

.debug-info summary {
  cursor: pointer;
  font-weight: 500;
  color: #6c757d;
  font-size: 12px;
}

.debug-content {
  margin-top: 12px;
}

.debug-content h6 {
  margin: 12px 0 4px 0;
  color: #495057;
  font-size: 12px;
}

.debug-info pre {
  margin: 8px 0;
  font-size: 10px;
  color: #495057;
  background-color: white;
  padding: 8px;
  border-radius: 4px;
  overflow-x: auto;
  max-height: 200px;
  overflow-y: auto;
}
</style>