<template>
  <div class="data-analyze">
    <h1>📊 Single Data Analysis</h1>
    <!-- Filter Control Panel -->
    <div class="filter-panel">
      <div class="filter-header">
        <h3>🔧 Filter Settings</h3>
        <div class="filter-count">
          Selected filters: {{ activeFiltersCount }}
        </div>
      </div>
      
      <div class="filters-grid">
        <!-- Basic Filters -->
        <div class="filter-group">
          <h4>Basic Filters</h4>
          
          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.eventType" 
              @change="onFilterToggle"
            >
              Event Type
            </el-checkbox>
            <el-select 
              v-model="filters.eventType" 
              :disabled="!filterEnabled.eventType"
              placeholder="Select event type" 
              clearable
              @change="updateResults"
            >
              <el-option label="Skipped Exon (SE)" value="SE" />
              <el-option label="Mutually Exclusive Exon (MXE)" value="MXE" />
              <el-option label="Alternative 5' Splice Site (A5SS)" value="A5SS" />
              <el-option label="Alternative 3' Splice Site (A3SS)" value="A3SS" />
              <el-option label="Retained Intron (RI)" value="RI" />
            </el-select>
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.chromosome" 
              @change="onFilterToggle"
            >
              Chromosome
            </el-checkbox>
            <el-select 
              v-model="filters.chromosome" 
              :disabled="!filterEnabled.chromosome"
              placeholder="Select chromosome" 
              clearable
              @change="updateResults"
            >
              <el-option 
                v-for="chr in chromosomes" 
                :key="chr" 
                :label="chr" 
                :value="chr" 
              />
            </el-select>
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.geneName" 
              @change="onFilterToggle"
            >
              Gene Name
            </el-checkbox>
            <el-input 
              v-model="filters.geneName" 
              :disabled="!filterEnabled.geneName"
              placeholder="Enter gene name" 
              clearable 
              @input="updateResults"
            />
          </div>

          <!-- GTF相关过滤器 -->
          <div v-if="store.isGtfEnabled" class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.geneAnnotation" 
              @change="onFilterToggle"
            >
              Gene Annotation Quality
            </el-checkbox>
            <el-select 
              v-model="filters.geneAnnotation" 
              :disabled="!filterEnabled.geneAnnotation"
              placeholder="Select annotation quality" 
              clearable
              @change="updateResults"
            >
              <el-option label="High Quality" value="high" />
              <el-option label="Medium Quality" value="medium" />
              <el-option label="Low Quality" value="low" />
            </el-select>
          </div>
        </div>

        <!-- Statistical Filters -->
        <div class="filter-group">
          <h4>Statistical Filters</h4>
          
          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.confidence" 
              @change="onFilterToggle"
            >
              Minimum Confidence
            </el-checkbox>
            <el-input-number 
              v-model="filters.confidence" 
              :disabled="!filterEnabled.confidence"
              :min="0" 
              :max="1" 
              :step="0.1" 
              :precision="1"
              placeholder="0.6"
              @change="updateResults"
            />
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.minReads" 
              @change="onFilterToggle"
            >
              Minimum Reads
            </el-checkbox>
            <el-input-number 
              v-model="filters.minReads" 
              :disabled="!filterEnabled.minReads"
              :min="0" 
              :step="100"
              placeholder="1000"
              @change="updateResults"
            />
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.isSignificant" 
              @change="onFilterToggle"
            >
              Show Significant Only
            </el-checkbox>
            <el-switch 
              v-model="filters.isSignificant" 
              :disabled="!filterEnabled.isSignificant"
              @change="updateResults"
            />
          </div>
        </div>
      </div>

      <!-- Action Buttons -->
      <div class="action-section">
        <el-button @click="resetFilters">Reset All Filters</el-button>
      </div>
    </div>

    <!-- Analysis Results Content Area -->
    <div class="content-section">
      <div class="content-header">
        <h3>📈 Analysis Results</h3>
        <div class="stats-summary">
          <span class="stat-item">Total: {{ results.total }}</span>
          <span class="stat-item">Significant: {{ results.significant }}</span>
        </div>
      </div>

      <!-- Results Table -->
      <div class="results-table" v-loading="loading">
        <el-table 
          :data="results.data" 
          style="width: 100%" 
          height="500"
          empty-text="No data available"
          @row-click="handleRowClick"
        >
          <el-table-column prop="id" label="ID" min-width="80" />
          <el-table-column prop="geneSymbol" label="Gene Symbol" min-width="120">
            <template #default="scope">
              <el-button 
                link  
                @click.stop="handleRowClick(scope.row)"
                class="gene-link"
              >
                {{ scope.row.geneSymbol }}
              </el-button>
            </template>
          </el-table-column>
          <el-table-column prop="geneId" label="Gene ID" min-width="180" />
          <el-table-column prop="eventType" label="Event Type" min-width="120">
            <template #default="scope">
              <el-tag size="small" :type="getEventTypeColor(scope.row.eventType)">
                {{ scope.row.eventType }}
              </el-tag>
            </template>
          </el-table-column>
          <el-table-column prop="chromosome" label="Chromosome" min-width="120" />
          <el-table-column prop="strand" label="Strand" min-width="80" />
          
          <!-- GTF相关列 -->
          <el-table-column 
            v-if="store.isGtfEnabled"
            prop="geneAnnotation" 
            label="Gene Annotation" 
            min-width="130"
          >
            <template #default="scope">
              <el-tag size="small" :type="getAnnotationColor(scope.row.geneAnnotation)">
                {{ scope.row.geneAnnotation || 'N/A' }}
              </el-tag>
            </template>
          </el-table-column>
          
          <el-table-column prop="reads1" label="Reads 1" min-width="100" />
          <el-table-column prop="reads2" label="Reads 2" min-width="100" />
          <el-table-column prop="confidence" label="Confidence" min-width="110">
            <template #default="scope">
              <span :class="getConfidenceClass(scope.row.confidence)">
                {{ scope.row.confidence?.toFixed(1) || 'N/A' }}
              </span>
            </template>
          </el-table-column>
          <el-table-column label="Significance" min-width="130">
            <template #default="scope">
              <el-tag 
                size="small"
                :type="scope.row.isSignificant ? 'success' : 'info'"
              >
                {{ scope.row.isSignificant ? 'Significant' : 'Not Significant' }}
              </el-tag>
            </template>
          </el-table-column>
        </el-table>
      </div>
      <div class="pagination-section">
        <div class="pagination-info">
          <span>
            Showing {{ ((pagination.currentPage - 1) * pagination.pageSize + 1) }} 
            to {{ Math.min(pagination.currentPage * pagination.pageSize, pagination.total) }} 
            of {{ pagination.total }} entries
          </span>
        </div>
        
        <el-pagination
          v-model:current-page="pagination.currentPage"
          v-model:page-size="pagination.pageSize"
          :page-sizes="pagination.pageSizes"
          :total="pagination.total"
          layout="total, sizes, prev, pager, next, jumper"
          :hide-on-single-page="false"
          @size-change="handlePageSizeChange"
          @current-change="handlePageChange"
        />
      </div>
    </div>

    <!-- Isoform Detail Modal -->
    <IsoformModal 
      v-if="store.isDataLoaded"
      v-model="showIsoformModal" 
      :row-data="selectedRowData"
      :rmats-results-path="store.getFilePaths.single"
      :sashimi-output-dir="store.getFilePaths.singleOutput"
      :enable-sashimi-plot="store.isSashimiEnabled"
      :enable-gtf-features="store.isGtfEnabled"
    />

    <!-- Detail Modal -->
    <el-dialog 
      v-model="showDetailModal" 
      title="Row Details" 
      width="60%"
      :close-on-click-modal="false"
    >
      <div v-if="selectedRowData" class="detail-content">
        <div class="detail-grid">
          <div class="detail-item">
            <strong>ID:</strong> {{ selectedRowData.id }}
          </div>
          <div class="detail-item">
            <strong>Gene Symbol:</strong> {{ selectedRowData.geneSymbol }}
          </div>
          <div class="detail-item">
            <strong>Gene ID:</strong> {{ selectedRowData.geneId }}
          </div>
          <div class="detail-item">
            <strong>Event Type:</strong> 
            <el-tag size="small" :type="getEventTypeColor(selectedRowData.eventType)">
              {{ selectedRowData.eventType }}
            </el-tag>
          </div>
          <div class="detail-item">
            <strong>Chromosome:</strong> {{ selectedRowData.chromosome }}
          </div>
          <div class="detail-item">
            <strong>Strand:</strong> {{ selectedRowData.strand }}
          </div>
          
          <!-- GTF相关详细信息 -->
          <div v-if="store.isGtfEnabled" class="detail-item">
            <strong>Gene Annotation:</strong>
            <el-tag size="small" :type="getAnnotationColor(selectedRowData.geneAnnotation)">
              {{ selectedRowData.geneAnnotation || 'N/A' }}
            </el-tag>
          </div>
          
          <div class="detail-item">
            <strong>Reads 1:</strong> {{ selectedRowData.reads1?.toLocaleString() }}
          </div>
          <div class="detail-item">
            <strong>Reads 2:</strong> {{ selectedRowData.reads2?.toLocaleString() }}
          </div>
          <div class="detail-item">
            <strong>Confidence:</strong> 
            <span :class="getConfidenceClass(selectedRowData.confidence)">
              {{ selectedRowData.confidence?.toFixed(2) }}
            </span>
          </div>
          <div class="detail-item">
            <strong>Significance:</strong>
            <el-tag 
              size="small"
              :type="selectedRowData.isSignificant ? 'success' : 'info'"
            >
              {{ selectedRowData.isSignificant ? 'Significant' : 'Not Significant' }}
            </el-tag>
          </div>
        </div>
      </div>
      <template #footer>
        <span class="dialog-footer">
          <el-button @click="showDetailModal = false">Close</el-button>
          <el-button type="primary" @click="exportSingleResult">Export This Result</el-button>
          <el-button 
            v-if="store.isGtfEnabled" 
            type="success"
          >
            Export with GTF
          </el-button>
          <el-button 
            v-if="store.isSashimiEnabled" 
            type="warning"
          >
            Generate Sashimi
          </el-button>
        </span>
      </template>
    </el-dialog>

    <!-- Status Message -->
    <div v-if="statusMessage" :class="statusMessageClass" class="status-message">
      {{ statusMessage }}
    </div>
  </div>
</template>

<script setup>
import { ref, computed, watch, onMounted } from 'vue'
import IsoformModal from './IsoformDetail.vue'
import { useAppStore } from '@/stores/app'

const store = useAppStore()

// Reactive Data
const loading = ref(false)
const statusMessage = ref('')
const statusMessageClass = ref('')
const selectedRowData = ref({})
const originalData = ref([])
const showDetailModal = ref(false)
const showIsoformModal = ref(false)

// 分页相关数据
const pagination = ref({
  currentPage: 1,
  pageSize: 100,
  total: 0,
  pageSizes: [50, 100, 200, 500]
})

// 修改 results 数据结构 - 只保留这一个定义
const results = ref({
  total: 0,
  significant: 0,
  data: [],
  allFilteredData: []
})

// Filter Enable States - 添加GTF过滤器
const filterEnabled = ref({
  eventType: false,
  chromosome: false,
  geneName: false,
  confidence: false,
  minReads: false,
  isSignificant: false,
  geneAnnotation: false  // GTF相关过滤器
})

// Filter Values - 添加GTF过滤值
const filters = ref({
  eventType: '',
  chromosome: '',
  geneName: '',
  confidence: null,
  minReads: null,
  isSignificant: false,
  geneAnnotation: ''  // GTF相关过滤值
})

// Chromosome Options - 将从数据中动态生成
const chromosomes = ref([])

// Computed Properties
const activeFiltersCount = computed(() => {
  return Object.values(filterEnabled.value).filter(Boolean).length
})

const hasResults = computed(() => {
  return results.value.allFilteredData.length > 0
})

// 分页处理函数
const handlePageChange = (page) => {
  pagination.value.currentPage = page
  updatePageData()
}

const handlePageSizeChange = (size) => {
  pagination.value.pageSize = size
  pagination.value.currentPage = 1
  updatePageData()
}

const updatePageData = () => {
  const start = (pagination.value.currentPage - 1) * pagination.value.pageSize
  const end = start + pagination.value.pageSize
  results.value.data = results.value.allFilteredData.slice(start, end)
}

// 加载真实数据 - 修复版本
const loadRealData = async () => {
  loading.value = true
  statusMessage.value = 'Loading data...'
  statusMessageClass.value = 'status-info'
  
  try {
    const eventTypes = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
    let allData = []
    let loadedFiles = 0
    let totalFiles = eventTypes.length
    
    // 加载所有事件类型的文件
    for (const eventType of eventTypes) {
      try {
        statusMessage.value = `Loading ${eventType} data...`
        
        const response = await fetch(`/frontend_data_${eventType}.json`)
        
        if (!response.ok) {
          console.warn(`Failed to load ${eventType} data: ${response.status}`)
          continue
        }
        
        const data = await response.json()
        
        // 如果数据量太大，限制加载数量
        const MAX_RECORDS_PER_TYPE = 1000000
        let processedData = Array.isArray(data) ? data.slice(0, MAX_RECORDS_PER_TYPE) : []
        
        // 为每个数据项添加事件类型
        processedData = processedData.map(item => ({
          ...item,
          eventType: item.eventType || eventType
        }))
        
        allData = allData.concat(processedData)
        loadedFiles++
        
        console.log(`Loaded ${processedData.length} records from ${eventType}`)
        
      } catch (fileError) {
        console.warn(`Error loading ${eventType} file:`, fileError)
        continue
      }
    }
    
    if (allData.length === 0) {
      throw new Error('No data files could be loaded')
    }
    
    // 数据预处理
    const cleanedData = preprocessData(allData)
    
    // 设置原始数据
    originalData.value = cleanedData
    
    // 计算统计信息
    const significantCount = cleanedData.filter(item => item.isSignificant).length
    
    // 设置分页数据
    pagination.value.total = cleanedData.length
    pagination.value.currentPage = 1
    
    // 设置结果数据
    results.value = {
      total: cleanedData.length,
      significant: significantCount,
      allFilteredData: [...cleanedData],
      data: cleanedData.slice(0, pagination.value.pageSize)
    }
    
    // 从数据中提取唯一的染色体列表
    const uniqueChromosomes = [...new Set(cleanedData.map(item => item.chromosome))]
      .filter(chr => chr && chr.trim())
      .sort((a, b) => {
        const aNum = a.replace('chr', '')
        const bNum = b.replace('chr', '')
        
        if (!isNaN(aNum) && !isNaN(bNum)) {
          return parseInt(aNum) - parseInt(bNum)
        }
        if (!isNaN(aNum)) return -1
        if (!isNaN(bNum)) return 1
        return aNum.localeCompare(bNum)
      })
    
    chromosomes.value = uniqueChromosomes
    
    console.log('Data loaded successfully:', {
      totalRecords: cleanedData.length,
      significantRecords: significantCount,
      loadedFiles: loadedFiles,
      totalFiles: totalFiles,
      pageSize: pagination.value.pageSize,
      currentPage: pagination.value.currentPage
    })
    
    statusMessage.value = `Data loaded successfully! ${cleanedData.length} records loaded with pagination.`
    statusMessageClass.value = 'status-success'
    
  } catch (error) {
    console.error('Error loading data:', error)
    
    statusMessage.value = `Failed to load data: ${error.message}`
    statusMessageClass.value = 'status-error'
    
    // 如果加载失败，使用备用的模拟数据
    console.log('Falling back to mock data...')
    const mockData = generateMockData()
    const processedMockData = preprocessData(mockData)
    
    originalData.value = processedMockData
    pagination.value.total = processedMockData.length
    pagination.value.currentPage = 1
    
    results.value = {
      total: processedMockData.length,
      significant: processedMockData.filter(item => item.isSignificant).length,
      allFilteredData: [...processedMockData],
      data: processedMockData.slice(0, pagination.value.pageSize)
    }
    
    chromosomes.value = [
      'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
      'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
      'chr21', 'chr22', 'chrX', 'chrY'
    ]
    
    statusMessage.value = `Using mock data due to loading failure: ${error.message}`
    statusMessageClass.value = 'status-error'
    
  } finally {
    loading.value = false
    setTimeout(() => {
      statusMessage.value = ''
    }, 5000)
  }
}

// 增强的数据预处理函数 - 添加GTF相关字段
const preprocessData = (data) => {
  return data.map((item, index) => {
    // 模拟GTF注释质量数据（仅在GTF启用时）
    const mockGeneAnnotation = store.isGtfEnabled ? 
      ['high', 'medium', 'low'][Math.floor(Math.random() * 3)] : null

    return {
      // 基础字段
      id: item.id || index + 1,
      geneSymbol: item.geneSymbol || 'Unknown',
      geneId: item.geneId || 'Unknown',
      eventType: item.eventType || 'Unknown',
      chromosome: item.chromosome || 'Unknown',
      strand: item.strand || '+',
      
      // 位置信息
      exonStart: typeof item.exonStart === 'number' ? item.exonStart : parseInt(item.exonStart || 0),
      exonEnd: typeof item.exonEnd === 'number' ? item.exonEnd : parseInt(item.exonEnd || 0),
      
      // 读数信息
      reads1: typeof item.reads1 === 'number' ? item.reads1 : parseInt(item.reads1 || 0),
      reads2: typeof item.reads2 === 'number' ? item.reads2 : parseInt(item.reads2 || 0),
      
      // 统计信息
      confidence: typeof item.confidence === 'number' ? item.confidence : parseFloat(item.confidence || 0),
      pValue: typeof item.pValue === 'number' ? item.pValue : parseFloat(item.pValue || 1),
      fdr: typeof item.fdr === 'number' ? item.fdr : parseFloat(item.fdr || 1),
      
      // 包含水平
      originalInc: typeof item.originalInc === 'string' ? item.originalInc : String(item.originalInc || '0'),
      compareInc: typeof item.compareInc === 'string' ? item.compareInc : String(item.compareInc || '0'),
      
      // 其他字段
      incLevelDifference: typeof item.incLevelDifference === 'number' ? item.incLevelDifference : parseFloat(item.incLevelDifference || 0),
      isSignificant: typeof item.isSignificant === 'boolean' ? item.isSignificant : Boolean(item.isSignificant),
      status: item.status || 'unknown',
      structureData: item.structureData,
      
      // GTF相关字段（仅在GTF启用时添加）
      ...(store.isGtfEnabled && {
        geneAnnotation: mockGeneAnnotation
      }),
      
      // 可选字段
      ...(item.shortStart && { shortStart: item.shortStart }),
      ...(item.shortEnd && { shortEnd: item.shortEnd })
    }
  })
}

// 生成模拟数据作为备用
const generateMockData = () => {
  const eventTypes = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
  const genes = ['BRCA1', 'TP53', 'EGFR', 'KRAS', 'PIK3CA', 'APC', 'PTEN', 'RB1', 'MYC', 'BRAF']
  const chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10']
  
  const data = []
  for (let i = 1; i <= 50; i++) {
    const confidence = Math.random()
    const reads1 = Math.floor(Math.random() * 10000) + 1000
    const reads2 = Math.floor(Math.random() * 10000) + 1000
    const originalInc = Math.random() * 0.8 + 0.1
    const compareInc = Math.random() * 0.8 + 0.1
    const pValue = Math.random() * 0.1
    const fdr = pValue * (1 + Math.random())
    
    data.push({
      id: i,
      geneSymbol: genes[Math.floor(Math.random() * genes.length)],
      geneId: `ENSG0000${String(Math.floor(Math.random() * 999999)).padStart(6, '0')}`,
      eventType: eventTypes[Math.floor(Math.random() * eventTypes.length)],
      chromosome: chroms[Math.floor(Math.random() * chroms.length)],
      strand: Math.random() > 0.5 ? '+' : '-',
      reads1: reads1,
      reads2: reads2,
      confidence: confidence,
      isSignificant: confidence > 0.7,
      originalInc: originalInc.toFixed(3),
      compareInc: compareInc.toFixed(3),
      pValue: pValue,
      fdr: Math.min(fdr, 1),
      exonStart: Math.floor(Math.random() * 1000000) + 1000000,
      exonEnd: Math.floor(Math.random() * 1000000) + 2000000
    })
  }
  return data
}

// 修改 updateResults 函数以支持分页和GTF过滤器
const updateResults = () => {
  loading.value = true
  
  setTimeout(() => {
    let filteredData = [...originalData.value]
    
    // 应用过滤条件
    if (filterEnabled.value.eventType && filters.value.eventType) {
      filteredData = filteredData.filter(item => item.eventType === filters.value.eventType)
    }
    
    if (filterEnabled.value.chromosome && filters.value.chromosome) {
      filteredData = filteredData.filter(item => item.chromosome === filters.value.chromosome)
    }
    
    if (filterEnabled.value.geneName && filters.value.geneName) {
      filteredData = filteredData.filter(item => 
        item.geneSymbol.toLowerCase().includes(filters.value.geneName.toLowerCase())
      )
    }
    
    if (filterEnabled.value.confidence && filters.value.confidence !== null) {
      filteredData = filteredData.filter(item => item.confidence >= filters.value.confidence)
    }
    
    if (filterEnabled.value.minReads && filters.value.minReads !== null) {
      filteredData = filteredData.filter(item => 
        item.reads1 >= filters.value.minReads || item.reads2 >= filters.value.minReads
      )
    }
    
    if (filterEnabled.value.isSignificant && filters.value.isSignificant) {
      filteredData = filteredData.filter(item => item.isSignificant === true)
    }
    
    // GTF相关过滤器
    if (filterEnabled.value.geneAnnotation && filters.value.geneAnnotation && store.isGtfEnabled) {
      filteredData = filteredData.filter(item => item.geneAnnotation === filters.value.geneAnnotation)
    }
    
    // 重置分页到第一页
    pagination.value.currentPage = 1
    pagination.value.total = filteredData.length
    
    // 更新结果
    results.value.allFilteredData = filteredData
    results.value.significant = filteredData.filter(item => item.isSignificant).length
    results.value.data = filteredData.slice(0, pagination.value.pageSize)
    
    loading.value = false
  }, 300)
}

// 组件挂载时加载数据
onMounted(() => {
  loadRealData()
})

// 其他方法保持不变
const onFilterToggle = () => {
  updateResults()
}

const resetFilters = () => {
  Object.keys(filterEnabled.value).forEach(key => {
    filterEnabled.value[key] = false
  })
  
  Object.keys(filters.value).forEach(key => {
    if (key === 'isSignificant') {
      filters.value[key] = false
    } else {
      filters.value[key] = key === 'confidence' || key === 'minReads' ? null : ''
    }
  })
  
  updateResults()
  
  statusMessage.value = 'All filters have been reset'
  statusMessageClass.value = 'status-info'
  setTimeout(() => statusMessage.value = '', 3000)
}

const exportResults = () => {
  statusMessage.value = 'Exporting results...'
  statusMessageClass.value = 'status-info'
  
  setTimeout(() => {
    statusMessage.value = 'Results exported successfully!'
    statusMessageClass.value = 'status-success'
    setTimeout(() => statusMessage.value = '', 3000)
  }, 1000)
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

// 简单的颜色映射函数
const getAnnotationColor = (annotation) => {
  const colors = {
    'high': 'success',
    'medium': 'warning', 
    'low': 'danger'
  }
  return colors[annotation] || 'info'
}

const handleRowClick = (row) => {
  selectedRowData.value = row
  showIsoformModal.value = true
  console.log('Row clicked:', row)
}

const exportSingleResult = () => {
  statusMessage.value = 'Exporting single result...'
  statusMessageClass.value = 'status-info'
  
  setTimeout(() => {
    statusMessage.value = 'Single result exported successfully!'
    statusMessageClass.value = 'status-success'
    showDetailModal.value = false
    setTimeout(() => statusMessage.value = '', 3000)
  }, 1000)
}

// Watch filter changes - 包含新的过滤器
watch(filterEnabled, (newVal, oldVal) => {
  Object.keys(newVal).forEach(key => {
    if (!newVal[key] && oldVal[key]) {
      if (key === 'isSignificant') {
        filters.value[key] = false
      } else {
        filters.value[key] = key === 'confidence' || key === 'minReads' ? null : ''
      }
    }
  })
}, { deep: true })
</script>

<style scoped>
.data-analyze {
  padding: 24px;
  max-width: 1400px;
  margin: 0 auto;
}

.filter-panel {
  background-color: #f8f9fa;
  border-radius: 8px;
  padding: 20px;
  margin-bottom: 24px;
  border: 1px solid #e9ecef;
}

.filter-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
}

.filter-header h3 {
  margin: 0;
  color: #495057;
}

.filter-count {
  font-size: 14px;
  color: #6c757d;
  font-weight: 500;
}

.filters-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 30px;
  margin-bottom: 20px;
}

.filter-group h4 {
  margin-top: 0;
  margin-bottom: 16px;
  color: #495057;
  font-size: 16px;
  border-bottom: 2px solid #dee2e6;
  padding-bottom: 8px;
}

.filter-item {
  margin-bottom: 16px;
}

.filter-item .el-checkbox {
  margin-bottom: 8px;
  font-weight: 500;
}

.action-section {
  text-align: center;
  padding-top: 16px;
  border-top: 1px solid #dee2e6;
}

.content-section {
  background-color: white;
  border-radius: 8px;
  padding: 20px;
  border: 1px solid #e9ecef;
}

.content-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
}

.content-header h3 {
  margin: 0;
  color: #495057;
}

.stats-summary {
  display: flex;
  gap: 20px;
}

.stat-item {
  padding: 6px 12px;
  background-color: #f8f9fa;
  border-radius: 4px;
  font-size: 14px;
  color: #495057;
  border: 1px solid #dee2e6;
}

.results-table {
  border-radius: 6px;
  overflow: hidden;
}

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

.status-message {
  margin-top: 20px;
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

h1 {
  color: #2c3e50;
  margin-bottom: 8px;
}

:deep(.el-input-number) {
  width: 100%;
}

:deep(.el-select) {
  width: 100%;
}

:deep(.el-checkbox) {
  margin-right: 0;
}

:deep(.el-switch) {
  margin-left: 0;
}

/* Detail Modal Styles */
.detail-content {
  padding: 20px 0;
}

.detail-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 16px;
  max-width: 100%;
}

.detail-item {
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 12px;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #dee2e6;
}

.detail-item strong {
  color: #495057;
  min-width: 120px;
  flex-shrink: 0;
}

.dialog-footer {
  display: flex;
  justify-content: flex-end;
  gap: 10px;
}
.status-error {
  background-color: #f8d7da;
  color: #721c24;
  border: 1px solid #f5c6cb;
}

/* 如果你想要更好的加载状态显示 */
.loading-overlay {
  position: relative;
}

.loading-overlay::before {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background-color: rgba(255, 255, 255, 0.8);
  z-index: 1000;
}

.loading-spinner {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  z-index: 1001;
}
.pagination-section {
  margin-top: 20px;
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 16px 0;
  border-top: 1px solid #e9ecef;
}

.pagination-info {
  color: #6c757d;
  font-size: 14px;
}

.pagination-info span {
  display: inline-block;
  padding: 8px 12px;
  background-color: #f8f9fa;
  border-radius: 4px;
  border: 1px solid #dee2e6;
}

/* 分页组件样式调整 */
:deep(.el-pagination) {
  display: flex;
  align-items: center;
  gap: 8px;
}

:deep(.el-pagination .el-pagination__sizes) {
  margin-right: 16px;
}

:deep(.el-pagination .el-pagination__total) {
  margin-right: 16px;
  color: #6c757d;
}

:deep(.el-pagination .el-pagination__jump) {
  margin-left: 16px;
}

/* 响应式设计 */
@media (max-width: 768px) {
  .pagination-section {
    flex-direction: column;
    gap: 12px;
    text-align: center;
  }
  
  .pagination-info {
    order: 2;
  }
  
  :deep(.el-pagination) {
    order: 1;
    flex-wrap: wrap;
    justify-content: center;
  }
}

/* 加载状态优化 */
.results-table {
  border-radius: 6px;
  overflow: hidden;
  min-height: 500px;  /* 确保表格有最小高度 */
}

/* 表格性能优化 */
:deep(.el-table) {
  --el-table-border-color: #e9ecef;
}

:deep(.el-table .el-table__row) {
  cursor: pointer;
}

:deep(.el-table .el-table__row:hover) {
  background-color: #f8f9fa;
}
</style>