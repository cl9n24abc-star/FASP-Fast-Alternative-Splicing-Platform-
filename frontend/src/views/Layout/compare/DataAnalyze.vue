<template>
  <div class="data-analyze">
    <h1>📊 Data Analysis</h1>
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

          <!-- 新增：状态过滤选项 -->
          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.status" 
              @change="onFilterToggle"
            >
              Data Status
            </el-checkbox>
            <el-select 
              v-model="filters.status" 
              :disabled="!filterEnabled.status"
              placeholder="Select data status" 
              clearable
              @change="updateResults"
            >
              <el-option label="Kept after filter" value="kept">
                <div class="status-option">
                  <div class="status-indicator status-kept"></div>
                  <span>Kept after filter</span>
                </div>
              </el-option>
              <el-option label="Removed by filter" value="removed">
                <div class="status-option">
                  <div class="status-indicator status-removed"></div>
                  <span>Removed by filter</span>
                </div>
              </el-option>
              <el-option label="Added by filter" value="added">
                <div class="status-option">
                  <div class="status-indicator status-added"></div>
                  <span>Added by filter</span>
                </div>
              </el-option>
            </el-select>
          </div>
        </div>

        <!-- Statistical Filters -->
        <div class="filter-group">
          <h4>Statistical Filters</h4>
          
          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.fdr" 
              @change="onFilterToggle"
            >
              FDR Threshold
            </el-checkbox>
            <el-input-number 
              v-model="filters.fdr" 
              :disabled="!filterEnabled.fdr"
              :min="0" 
              :max="1" 
              :step="0.01" 
              :precision="3"
              placeholder="0.05"
              @change="updateResults"
            />
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.pValue" 
              @change="onFilterToggle"
            >
              P-value Threshold
            </el-checkbox>
            <el-input-number 
              v-model="filters.pValue" 
              :disabled="!filterEnabled.pValue"
              :min="0" 
              :max="1" 
              :step="0.001" 
              :precision="4"
              placeholder="0.05"
              @change="updateResults"
            />
          </div>

          <div class="filter-item">
            <el-checkbox 
              v-model="filterEnabled.minIncLevelDiff" 
              @change="onFilterToggle"
            >
              Minimum IncLevel Difference
            </el-checkbox>
            <el-input-number 
              v-model="filters.minIncLevelDiff" 
              :disabled="!filterEnabled.minIncLevelDiff"
              :min="0" 
              :max="1" 
              :step="0.01" 
              :precision="3"
              placeholder="0.1"
              @change="updateResults"
            />
          </div>
        </div>
      </div>

      <!-- Action Buttons -->
      <div class="action-section">
        <el-button @click="resetFilters">Reset All Filters</el-button>
        <el-button type="success" @click="exportResults" :disabled="!hasResults">
          Export Results
        </el-button>
      </div>
    </div>

    <!-- Analysis Results Content Area -->
    <div class="content-section">
      <div class="content-header">
        <h3>📈 Analysis Results</h3>
        <div class="stats-summary">
          <span class="stat-item">Total: {{ results.total }}</span>
          <span class="stat-item">Filtered: {{ results.filtered }}</span>
          <span class="stat-item">Significant: {{ results.significant }}</span>
        </div>
      </div>

      <!-- Status Legend -->
      <div class="status-legend">
        <div class="legend-item">
          <div class="status-indicator status-removed"></div>
          <span>Removed by filter</span>
        </div>
        <div class="legend-item">
          <div class="status-indicator status-added"></div>
          <span>Added by filter</span>
        </div>
        <div class="legend-item">
          <div class="status-indicator status-kept"></div>
          <span>Kept after filter</span>
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
          <!-- Status Column -->
          <el-table-column label="Status" width="80" align="center">
            <template #default="scope">
              <div 
                :class="['status-indicator', getStatusClass(scope.row.status)]"
                :title="getStatusTooltip(scope.row.status)"
              >
              </div>
            </template>
          </el-table-column>

          <el-table-column prop="id" label="ID" width="80" />
          <el-table-column prop="geneSymbol" label="Gene" width="120">
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
          <el-table-column prop="eventType" label="Event Type" width="100">
            <template #default="scope">
              <el-tag size="small" :type="getEventTypeColor(scope.row.eventType)">
                {{ scope.row.eventType }}
              </el-tag>
            </template>
          </el-table-column>
          <el-table-column prop="chromosome" label="Chromosome" width="100" />
          <el-table-column prop="strand" label="Strand" width="60" />
          <el-table-column prop="exonStart" label="Exon Start" width="120" />
          <el-table-column prop="exonEnd" label="Exon End" width="120" />
          <el-table-column prop="pValue" label="P-value" width="120">
            <template #default="scope">
              <span :class="{ 'significant': scope.row.pValue < 0.05 }">
                {{ scope.row.pValue?.toExponential(3) || 'N/A' }}
              </span>
            </template>
          </el-table-column>
          <el-table-column prop="fdr" label="FDR" width="120">
            <template #default="scope">
              <span :class="{ 'significant': scope.row.fdr < 0.05 }">
                {{ scope.row.fdr?.toExponential(3) || 'N/A' }}
              </span>
            </template>
          </el-table-column>
          <el-table-column prop="incLevelDiff" label="IncLevel Difference" width="140">
            <template #default="scope">
              <span :class="getIncLevelClass(scope.row.incLevelDiff)">
                {{ scope.row.incLevelDiff?.toFixed(3) || 'N/A' }}
              </span>
            </template>
          </el-table-column>
          <el-table-column prop="originalInc" label="Original IncLevel" width="140" />
          <el-table-column prop="compareInc" label="Compare IncLevel" width="140" />
          <el-table-column label="Significance" width="100">
            <template #default="scope">
              <el-tag 
                size="small"
                :type="isSignificant(scope.row) ? 'success' : 'info'"
              >
                {{ isSignificant(scope.row) ? 'Significant' : 'Not Significant' }}
              </el-tag>
            </template>
          </el-table-column>
        </el-table>
      </div>
    </div>

    <!-- Isoform Detail Modal -->
    <IsoformModal 
      v-model="showIsoformModal" 
      :row-data="selectedRowData"
    />

    <!-- Status Message -->
    <div v-if="statusMessage" :class="statusMessageClass" class="status-message">
      {{ statusMessage }}
    </div>
  </div>
</template>
  
  <script setup>
  import { ref, computed, watch } from 'vue'
  import IsoformModal from './IsoformDetail.vue'
  
  // Reactive Data
  const loading = ref(false)
  const statusMessage = ref('')
  const statusMessageClass = ref('')
  const showIsoformModal = ref(false)
  const selectedRowData = ref({})
  const originalData = ref([])
  
  // Filter Enable States
  const filterEnabled = ref({
    eventType: false,
    chromosome: false,
    geneName: false,
    fdr: false,
    pValue: false,
    status: false,
    minIncLevelDiff: false
  })
  
  // Filter Values
  const filters = ref({
    eventType: '',
    chromosome: '',
    geneName: '',
    fdr: null,
    pValue: null,
    status: '',
    minIncLevelDiff: null
  })
  
  // Analysis Results
  const results = ref({
    total: 1250,
    filtered: 1250,
    significant: 856,
    data: []
  })
  
  // Chromosome Options
  const chromosomes = ref([
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY'
  ])
  
  // Computed Properties
  const activeFiltersCount = computed(() => {
    return Object.values(filterEnabled.value).filter(Boolean).length
  })
  
  const hasResults = computed(() => {
    return results.value.data.length > 0
  })
  
  // Generate Mock Data
  const generateMockData = () => {
    const eventTypes = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
    const genes = ['BRCA1', 'TP53', 'EGFR', 'KRAS', 'PIK3CA', 'APC', 'PTEN', 'RB1', 'MYC', 'BRAF']
    const chroms = chromosomes.value.slice(0, 10)
    
    const data = []
    for (let i = 1; i <= 50; i++) {
      const pValue = Math.random() * 0.1
      const fdr = pValue * (1 + Math.random())
      const incLevelDiff = (Math.random() - 0.5) * 0.6
      
      data.push({
      id: i,
      geneSymbol: genes[Math.floor(Math.random() * genes.length)],
      eventType: eventTypes[Math.floor(Math.random() * eventTypes.length)],
      chromosome: chroms[Math.floor(Math.random() * chroms.length)],
      strand: Math.random() > 0.5 ? '+' : '-',
      exonStart: Math.floor(Math.random() * 1000000) + 1000000,
      exonEnd: Math.floor(Math.random() * 1000000) + 2000000,
      pValue: pValue,
      fdr: Math.min(fdr, 1),
      incLevelDiff: incLevelDiff,
      originalInc: (Math.random() * 0.8 + 0.1).toFixed(3),
      compareInc: (Math.random() * 0.8 + 0.1).toFixed(3),
      status: 'kept' // 👈 添加这一行
    })
    }
    return data
  }
  
  // Initialize Data
  originalData.value = generateMockData()
  results.value.data = [...originalData.value]
  
  // Methods
  const onFilterToggle = () => {
    updateResults()
  }
  // 新增：状态相关方法
const getStatusClass = (status) => {
  switch(status) {
    case 'removed': return 'status-removed'
    case 'added': return 'status-added'
    case 'kept': return 'status-kept'
    default: return 'status-kept'
  }
}

const getStatusTooltip = (status) => {
  switch(status) {
    case 'removed': return 'This item was in original data but removed by current filter'
    case 'added': return 'This item was not in original data but added by current filter'
    case 'kept': return 'This item exists in both original and filtered data'
    default: return 'Unknown status'
  }
}

const updateResults = () => {
  loading.value = true
  
  setTimeout(() => {
    // 从原始数据开始过滤，而不是重新生成
    let filteredData = [...originalData.value]
    
    // Apply Filter Conditions (除了状态过滤，其他过滤照常)
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
    
    if (filterEnabled.value.fdr && filters.value.fdr !== null) {
      filteredData = filteredData.filter(item => item.fdr <= filters.value.fdr)
    }
    
    if (filterEnabled.value.pValue && filters.value.pValue !== null) {
      filteredData = filteredData.filter(item => item.pValue <= filters.value.pValue)
    }
    
    if (filterEnabled.value.minIncLevelDiff && filters.value.minIncLevelDiff !== null) {
      filteredData = filteredData.filter(item => Math.abs(item.incLevelDiff) >= filters.value.minIncLevelDiff)
    }

    // 计算状态
    const filteredIds = new Set(filteredData.map(item => item.id))
    
    // 设置过滤后保留的数据状态为 'kept'
    filteredData = filteredData.map(item => ({
      ...item,
      status: 'kept'
    }))

    // 模拟一些被移除的数据
    const removedItems = originalData.value
      .filter(item => !filteredIds.has(item.id))
      .slice(0, 5)
      .map(item => ({
        ...item,
        status: 'removed'
      }))

    // 模拟一些新增的数据
    const addedItems = []
    if (filteredData.length > 0) {
      for (let i = 0; i < 3; i++) {
        const baseItem = filteredData[Math.floor(Math.random() * filteredData.length)]
        addedItems.push({
          ...baseItem,
          id: 1000 + i,
          status: 'added',
          geneSymbol: baseItem.geneSymbol + '_variant'
        })
      }
    }

    // 合并所有数据用于显示
    let displayData = [...filteredData, ...removedItems, ...addedItems]
    
    // ⚠️ 状态过滤应该在最后进行，在生成完整的displayData之后
    if (filterEnabled.value.status && filters.value.status) {
      displayData = displayData.filter(item => item.status === filters.value.status)
    }
    
    results.value.data = displayData
    results.value.filtered = filteredData.length
    results.value.significant = filteredData.filter(item => isSignificant(item)).length
    
    loading.value = false
  }, 500)
}
  
  const resetFilters = () => {
    Object.keys(filterEnabled.value).forEach(key => {
      filterEnabled.value[key] = false
    })
    
    Object.keys(filters.value).forEach(key => {
      filters.value[key] = key.includes('Level') || key === 'fdr' || key === 'pValue' ? null : ''
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
  
  const getIncLevelClass = (value) => {
    if (Math.abs(value) > 0.2) return 'high-diff'
    if (Math.abs(value) > 0.1) return 'medium-diff'
    return 'low-diff'
  }
  
  const isSignificant = (row) => {
    return row.pValue < 0.05 && row.fdr < 0.05
  }
  
  // Handle Row Click Event
  const handleRowClick = (row) => {
    selectedRowData.value = row
    showIsoformModal.value = true
  }
  
  // Watch filter changes, clear corresponding values when filters are disabled
  watch(filterEnabled, (newVal, oldVal) => {
    Object.keys(newVal).forEach(key => {
      if (!newVal[key] && oldVal[key]) {
        filters.value[key] = key.includes('Level') || key === 'fdr' || key === 'pValue' ? null : ''
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
  
  .significant {
    color: #dc3545;
    font-weight: bold;
  }
  
  .high-diff {
    color: #dc3545;
    font-weight: bold;
  }
  
  .medium-diff {
    color: #fd7e14;
    font-weight: 500;
  }
  
  .low-diff {
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
  /* 为状态选项添加样式 */
.status-option {
  display: flex;
  align-items: center;
  gap: 8px;
}

.status-option .status-indicator {
  width: 12px;
  height: 12px;
  border-radius: 50%;
  display: inline-block;
}

.status-option .status-kept {
  background-color: #28a745;
  border: 1px solid #1e7e34;
}

.status-option .status-removed {
  background-color: #007bff;
  border: 1px solid #0056b3;
}

.status-option .status-added {
  background-color: #dc3545;
  border: 1px solid #a71e2a;
}
/* 状态图例样式 */
.status-legend {
  display: flex;
  gap: 20px;
  margin-bottom: 16px;
  padding: 12px;
  background-color: #f8f9fa;
  border-radius: 6px;
  border: 1px solid #dee2e6;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: 14px;
  color: #495057;
}

/* 状态指示器样式 */
.status-indicator {
  width: 16px;
  height: 16px;
  border-radius: 50%;
  display: inline-block;
  cursor: help;
}

.status-removed {
  background-color: #007bff;
  border: 2px solid #0056b3;
}

.status-added {
  background-color: #dc3545;
  border: 2px solid #a71e2a;
}

.status-kept {
  background-color: #28a745;
  border: 2px solid #1e7e34;
}
  </style>