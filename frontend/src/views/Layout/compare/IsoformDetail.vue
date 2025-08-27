<template>
  <el-dialog
    v-model="visible"
    title="🧬 Isoform Details"
    width="80%"
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
            <strong>Gene:</strong>
            <span>{{ isoformData.geneSymbol }}</span>
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
            <strong>Position:</strong>
            <span>{{ isoformData.exonStart }} - {{ isoformData.exonEnd }}</span>
          </div>
        </div>
      </div>

      <!-- Abundance Information -->
      <div class="abundance-section">
        <h4>📊 Isoform Abundance</h4>
        <div class="abundance-grid">
          <div class="abundance-card">
            <div class="card-header original">Original Data</div>
            <div class="card-content">
              <div class="metric">
                <span class="label">TPM:</span>
                <span class="value">{{ abundanceData.original.tpm }}</span>
              </div>
              <div class="metric">
                <span class="label">Read Counts:</span>
                <span class="value">{{ abundanceData.original.reads }}</span>
              </div>
              <div class="metric">
                <span class="label">IncLevel:</span>
                <span class="value">{{ abundanceData.original.incLevel }}</span>
              </div>
            </div>
          </div>

          <div class="abundance-card">
            <div class="card-header processed">Processed data</div>
            <div class="card-content">
              <div class="metric">
                <span class="label">TPM:</span>
                <span class="value">{{ abundanceData.processed.tpm }}</span>
              </div>
              <div class="metric">
                <span class="label">Read Counts:</span>
                <span class="value">{{ abundanceData.processed.reads }}</span>
              </div>
              <div class="metric">
                <span class="label">IncLevel:</span>
                <span class="value">{{ abundanceData.processed.incLevel }}</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- Differential Expression -->
      <div class="expression-section">
        <h4>📈 Differential Expression Status</h4>
        <div class="expression-content">
          <div class="expression-stats">
            <div class="stat-item">
              <span class="stat-label">Expression Status:</span>
              <el-tag :type="getExpressionStatusColor(expressionData.status)" size="large">
                {{ getExpressionStatusText(expressionData.status) }}
              </el-tag>
            </div>
            <div class="stat-item">
              <span class="stat-label">Log2 Fold Change:</span>
              <span 
                class="stat-value"
                :class="getFoldChangeClass(expressionData.log2FoldChange)"
              >
                {{ expressionData.log2FoldChange }}
              </span>
            </div>
            <div class="stat-item">
              <span class="stat-label">P-value:</span>
              <span 
                class="stat-value"
                :class="{ 'significant': expressionData.pValue < 0.05 }"
              >
                {{ expressionData.pValue.toExponential(3) }}
              </span>
            </div>
            <div class="stat-item">
              <span class="stat-label">FDR:</span>
              <span 
                class="stat-value"
                :class="{ 'significant': expressionData.fdr < 0.05 }"
              >
                {{ expressionData.fdr.toExponential(3) }}
              </span>
            </div>
          </div>
        </div>
      </div>

      <!-- Delta PSI -->
      <div class="psi-section">
        <h4>🔄 Delta PSI (Splicing Inclusion Change)</h4>
        <div class="psi-content">
          <div class="psi-comparison">
            <div class="psi-item">
              <span class="psi-label">Original Data PSI:</span>
              <div class="psi-bar">
                <div 
                  class="psi-fill original" 
                  :style="{ width: (psiData.originalPsi * 100) + '%' }"
                ></div>
                <span class="psi-value">{{ psiData.originalPsi.toFixed(3) }}</span>
              </div>
            </div>
            <div class="psi-item">
              <span class="psi-label">Processed Data PSI:</span>
              <div class="psi-bar">
                <div 
                  class="psi-fill processed"
                  :style="{ width: (psiData.processedPsi * 100) + '%' }"
                ></div>
                <span class="psi-value">{{ psiData.processedPsi.toFixed(3) }}</span>
              </div>
            </div>
          </div>
          
          <div class="delta-summary">
            <div class="delta-value" :class="getDeltaPsiClass(psiData.deltaPsi)">
              Delta PSI: {{ psiData.deltaPsi > 0 ? '+' : '' }}{{ psiData.deltaPsi.toFixed(3) }}
            </div>
            <div class="delta-interpretation">
              <el-tag 
                :type="Math.abs(psiData.deltaPsi) > 0.1 ? 'warning' : 'info'"
                size="small"
              >
                {{ Math.abs(psiData.deltaPsi) > 0.1 ? 'Significant Change' : 'Minor Change' }}
              </el-tag>
            </div>
          </div>
        </div>
      </div>

      <!-- Quality Metrics -->
      <div class="quality-section">
        <h4>🎯 Quality Metrics</h4>
        <div class="quality-grid">
          <div class="quality-item">
            <span class="quality-label">Coverage:</span>
            <span class="quality-value">{{ qualityData.coverage }}x</span>
          </div>
          <div class="quality-item">
            <span class="quality-label">Junction Reads:</span>
            <span class="quality-value">{{ qualityData.junctionReads }}</span>
          </div>
          <div class="quality-item">
            <span class="quality-label">Mapping Quality:</span>
            <span class="quality-value">{{ qualityData.mappingQuality }}</span>
          </div>
          <div class="quality-item">
            <span class="quality-label">Confidence:</span>
            <el-tag :type="qualityData.confidence > 0.8 ? 'success' : 'warning'" size="small">
              {{ (qualityData.confidence * 100).toFixed(1) }}%
            </el-tag>
          </div>
        </div>
      </div>

      <!-- Gene Structure Visualization -->
      <div class="gene-structure-section">
        <h4>🧬 Gene Structure Visualization</h4>
        <div class="structure-info">
          <div class="structure-header">
            <span class="gene-name">{{ isoformData.geneSymbol }}</span>
            <span class="coordinates">{{ isoformData.chromosome }}:{{ structureData.geneStart }}-{{ structureData.geneEnd }}</span>
            <el-tag size="small">{{ structureData.strand }}</el-tag>
          </div>
          
          <!-- Alternative Splicing Event Labels -->
          <div class="splicing-events">
            <span class="events-label">Detected Splicing Events:</span>
            <el-tag 
              v-for="event in structureData.splicingEvents" 
              :key="event"
              :type="getEventColor(event)"
              size="small"
            >
              {{ event }}
            </el-tag>
          </div>
        </div>

        <!-- Exon-Intron Structure Diagram -->
        <div class="structure-diagram">
          <div class="diagram-container">
            
            <!-- Overall Gene Structure -->
            <div class="gene-track">
              <div class="track-label">Gene Structure:</div>
              <div class="gene-line">
                <div 
                  v-for="(element, index) in structureData.structure" 
                  :key="index"
                  :class="['structure-element', element.type]"
                  :style="{ 
                    width: element.width + '%',
                    backgroundColor: getElementColor(element.type, element.isAffected)
                  }"
                  @mouseover="showElementInfo(element, $event)"
                  @mouseleave="hideElementInfo"
                  @click="element.type === 'exon' ? handleExonClick(element) : null"
                >
                  <span v-if="element.type === 'exon'" class="element-label">
                    {{ element.name }}
                  </span>
                </div>
              </div>
            </div>

            <!-- Original Transcript -->
            <div class="transcript-track">
              <div class="track-label">Original Transcript:</div>
              <div class="transcript-line">
                <div 
                  v-for="(element, index) in structureData.originalTranscript" 
                  :key="index"
                  :class="['structure-element', element.type]"
                  :style="{ 
                    width: element.width + '%',
                    backgroundColor: getTranscriptColor(element.type, 'original')
                  }"
                  @click="element.type === 'exon' ? handleExonClick(element) : null"
                >
                  <span v-if="element.type === 'exon'" class="element-label">
                    {{ element.name }}
                  </span>
                </div>
              </div>
            </div>

            <!-- Alternative Splicing Transcript -->
            <div class="transcript-track">
              <div class="track-label">Alternative Splicing Transcript:</div>
              <div class="transcript-line">
                <div 
                  v-for="(element, index) in structureData.alternativeTranscript" 
                  :key="index"
                  :class="['structure-element', element.type, { 'skipped': element.skipped, 'included': element.included }]"
                  :style="{ 
                    width: element.width + '%',
                    backgroundColor: getTranscriptColor(element.type, 'alternative', element.skipped)
                  }"
                  @click="element.type === 'exon' ? handleExonClick(element) : null"
                >
                  <span v-if="element.type === 'exon'" class="element-label">
                    {{ element.name }}
                  </span>
                </div>
              </div>
            </div>
          </div>

          <!-- Legend -->
                  <div class="legend">
          <div 
            v-for="item in getLegendItems" 
            :key="item.label"
            class="legend-item"
          >
            <div :class="['legend-color', item.color]"></div>
            <span>{{ item.label }}</span>
          </div>
        </div>
        </div>

        <!-- Exon Details -->
        <div class="exon-details">
          <h5>📊 Exon Details</h5>
          <div class="exon-table">
            <el-table :data="structureData.exonDetails" size="small" max-height="200">
              <el-table-column prop="name" label="Exon" width="100" />
              <el-table-column prop="start" label="Start Position" width="120">
                <template #default="scope">
                  {{ scope.row.start.toLocaleString() }}
                </template>
              </el-table-column>
              <el-table-column prop="end" label="End Position" width="120">
                <template #default="scope">
                  {{ scope.row.end.toLocaleString() }}
                </template>
              </el-table-column>
              <el-table-column prop="length" label="Length (bp)" width="100" />
              <el-table-column prop="inclusionLevel" label="Inclusion Level" width="120">
                <template #default="scope">
                  <span :class="getInclusionClass(scope.row.inclusionLevel)">
                    {{ scope.row.inclusionLevel.toFixed(3) }}
                  </span>
                </template>
              </el-table-column>
              <el-table-column prop="status" label="Status" width="120">
                <template #default="scope">
                  <el-tag :type="getStatusColor(scope.row.status)" size="small">
                    {{ getStatusText(scope.row.status) }}
                  </el-tag>
                </template>
              </el-table-column>
            </el-table>
          </div>
        </div>
      </div>

    </div>

    <template #footer>
      <div class="dialog-footer">
        <el-button @click="handleClose">Close</el-button>
        <el-button type="primary" @click="exportData">Export Data</el-button>
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
  
  // Reactive Data
  const loading = ref(false)
  const isoformData = ref(null)
  const abundanceData = ref(null)
  const expressionData = ref(null)
  const psiData = ref(null)
  const qualityData = ref(null)
  const structureData = ref(null)
  const elementTooltip = ref({ show: false, content: '', x: 0, y: 0 })
  const showExonDialog = ref(false)      // 添加这行
  const selectedExonData = ref(null)     // 添加这行
  
  const getLegendItems = computed(() => {
  const eventType = isoformData.value?.eventType || 'SE'
  
  const baseLegend = [
    { color: 'exon-color', label: 'Exon' },
    { color: 'intron-color', label: 'Intron' },
    { color: 'affected-color', label: 'Affected Exon' }
  ]
  
  const eventSpecific = {
    'SE': { color: 'skipped-color', label: 'Skipped Exon' },
    'MXE': { color: 'exclusive-color', label: 'Mutually Exclusive Exon' },
    'A5SS': { color: 'alt5-color', label: 'Alternative 5\' Site' },
    'A3SS': { color: 'alt3-color', label: 'Alternative 3\' Site' },
    'RI': { color: 'retained-color', label: 'Retained Intron' }
  }
  
  return [...baseLegend, eventSpecific[eventType]]
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
  
  // Generate Mock Detail Data
  const generateDetailData = (rowData) => {
    if (!rowData || !rowData.id) return null
  
    // Generate detailed information based on original data
    const originalPsi = parseFloat(rowData.originalInc) || Math.random()
    const noisyrPsi = parseFloat(rowData.noisyrInc) || Math.random()
    
    // Generate Gene Structure Data
    const generateGeneStructure = () => {
      const eventTypes = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
      const eventType = rowData.eventType || eventTypes[Math.floor(Math.random() * eventTypes.length)]
      const exonCount = Math.floor(Math.random() * 8) + 4 // 4-12 exons
      
      // Generate exon structure
      const structure = []
      const originalTranscript = []
      const alternativeTranscript = []
      const exonDetails = []
      
      for (let i = 1; i <= exonCount; i++) {
        const exonLength = Math.floor(Math.random() * 200) + 100
        const intronLength = Math.floor(Math.random() * 1000) + 500
        
        // Determine which exon is affected by alternative splicing
        const isAffectedExon = (eventType === 'SE' && i === Math.floor(exonCount / 2)) ||
                              (eventType === 'MXE' && (i === Math.floor(exonCount / 2) || i === Math.floor(exonCount / 2) + 1))
        
        // Exon
        structure.push({
          type: 'exon',
          name: `E${i}`,
          width: 8,
          isAffected: isAffectedExon,
          length: exonLength,
          start: 1000000 + (i - 1) * 1500,
          end: 1000000 + (i - 1) * 1500 + exonLength
        })
        
        // Original transcript - contains all exons
        originalTranscript.push({
          type: 'exon',
          name: `E${i}`,
          width: 8,
          length: exonLength
        })
        
        // Alternative splicing transcript
        if (eventType === 'SE' && isAffectedExon) {
          // Skipped exon event - skip affected exon
          alternativeTranscript.push({
            type: 'exon',
            name: `E${i}`,
            width: 8,
            skipped: true,
            length: exonLength
          })
        } else {
          alternativeTranscript.push({
            type: 'exon',
            name: `E${i}`,
            width: 8,
            included: true,
            length: exonLength
          })
        }
        
        // Exon details
        exonDetails.push({
          name: `Exon ${i}`,
          start: 1000000 + (i - 1) * 1500,
          end: 1000000 + (i - 1) * 1500 + exonLength,
          length: exonLength,
          inclusionLevel: isAffectedExon ? Math.random() * 0.5 + 0.2 : Math.random() * 0.3 + 0.7,
          status: isAffectedExon ? 'Alternative' : 'Constitutive'
        })
        
        // Intron (except for last exon)
        if (i < exonCount) {
          structure.push({
            type: 'intron',
            name: `I${i}`,
            width: 4,
            length: intronLength
          })
          
          originalTranscript.push({
            type: 'intron',
            name: `I${i}`,
            width: 4
          })
          
          alternativeTranscript.push({
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
        originalTranscript,
        alternativeTranscript,
        exonDetails
      }
    }
    
    return {
      // Abundance data
      abundance: {
        original: {
          tpm: (Math.random() * 50 + 10).toFixed(2),
          reads: Math.floor(Math.random() * 1000) + 100,
          incLevel: originalPsi.toFixed(3)
        },
      processed: {
          tpm: (Math.random() * 50 + 10).toFixed(2),
          reads: Math.floor(Math.random() * 1000) + 100,
          incLevel: noisyrPsi.toFixed(3)
        }
      },
  
      // Differential expression data
      expression: {
        status: rowData.pValue < 0.05 ? (Math.random() > 0.5 ? 'upregulated' : 'downregulated') : 'unchanged',
        log2FoldChange: ((Math.random() - 0.5) * 4).toFixed(3),
        pValue: rowData.pValue || Math.random() * 0.1,
        fdr: rowData.fdr || Math.random() * 0.15
      },
  
      // PSI data
      psi: {
        originalPsi: originalPsi,
        processedPsi: noisyrPsi,
        deltaPsi: noisyrPsi - originalPsi
      },
  
      // Quality data
      quality: {
        coverage: Math.floor(Math.random() * 100) + 20,
        junctionReads: Math.floor(Math.random() * 500) + 50,
        mappingQuality: Math.floor(Math.random() * 20) + 80,
        confidence: Math.random() * 0.3 + 0.7
      },
  
      // Gene structure data
      structure: generateGeneStructure()
    }
  }
  
  // Load detailed data
  const loadDetailData = async () => {
    if (!props.rowData || !props.rowData.id) return
  
    loading.value = true
    
    try {
      // Simulate API call delay
      await new Promise(resolve => setTimeout(resolve, 800))
      
      isoformData.value = props.rowData
      const detailData = generateDetailData(props.rowData)
      
      abundanceData.value = detailData.abundance
      expressionData.value = detailData.expression
      psiData.value = detailData.psi
      qualityData.value = detailData.quality
      structureData.value = detailData.structure
      
    } catch (error) {
      console.error('Failed to load detailed data:', error)
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
      abundance: abundanceData.value,
      expression: expressionData.value,
      psi: psiData.value,
      quality: qualityData.value
    })
    // Here you can implement actual export functionality
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
  
  const getExpressionStatusColor = (status) => {
    const colors = {
      'upregulated': 'danger',
      'downregulated': 'success',
      'unchanged': 'info'
    }
    return colors[status] || 'info'
  }
  
  const getExpressionStatusText = (status) => {
    const texts = {
      'upregulated': 'Upregulated',
      'downregulated': 'Downregulated', 
      'unchanged': 'Unchanged'
    }
    return texts[status] || 'Unknown'
  }
  
  const getFoldChangeClass = (value) => {
    const val = parseFloat(value)
    if (val > 1) return 'upregulated'
    if (val < -1) return 'downregulated'
    return 'unchanged'
  }
  
  const getDeltaPsiClass = (value) => {
    if (Math.abs(value) > 0.2) return 'high-change'
    if (Math.abs(value) > 0.1) return 'medium-change'
    return 'low-change'
  }
  
  // Gene structure related methods
  const getEventColor = (eventType) => {
    const colors = {
      'SE': 'primary',
      'MXE': 'success',
      'A5SS': 'warning', 
      'A3SS': 'danger',
      'RI': 'info'
    }
    return colors[eventType] || 'default'
  }
  
  const getElementColor = (type, isAffected = false) => {
    if (type === 'exon') {
      return isAffected ? '#ff6b6b' : '#4ecdc4'
    } else {
      return '#95a5a6'
    }
  }
  
  const getTranscriptColor = (type, transcriptType, isSkipped = false) => {
    if (type === 'exon') {
      if (isSkipped) return '#ff9999'
      return transcriptType === 'original' ? '#74b9ff' : '#00b894'
    } else {
      return '#bdc3c7'
    }
  }
  
  const getInclusionClass = (level) => {
    if (level > 0.7) return 'high-inclusion'
    if (level > 0.3) return 'medium-inclusion'
    return 'low-inclusion'
  }
  
  const getStatusColor = (status) => {
    return status === 'Alternative' ? 'warning' : 'success'
  }
  
  const getStatusText = (status) => {
    return status === 'Alternative' ? 'Alternative' : 'Constitutive'
  }
  
  const showElementInfo = (element, event) => {
    elementTooltip.value = {
      show: true,
      content: `${element.name}: ${element.length}bp`,
      x: event.clientX,
      y: event.clientY
    }
  }
  
  const hideElementInfo = () => {
    elementTooltip.value.show = false
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
    }
  })
  
  </script>
  
  <style scoped>
  .modal-content {
    max-height: 600px;
    overflow-y: auto;
  }
  
  .loading-content {
    padding: 20px;
  }
  
  .info-section,
  .abundance-section,
  .expression-section,
  .psi-section,
  .quality-section {
    margin-bottom: 24px;
    padding: 16px;
    background-color: #f8f9fa;
    border-radius: 6px;
    border: 1px solid #e9ecef;
  }
  
  .info-section h4,
  .abundance-section h4,
  .expression-section h4,
  .psi-section h4,
  .quality-section h4 {
    margin-top: 0;
    margin-bottom: 12px;
    color: #495057;
    border-bottom: 2px solid #dee2e6;
    padding-bottom: 6px;
  }
  
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
  
  .abundance-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 16px;
  }
  
  .abundance-card {
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
  
  .card-header.original {
    background-color: #28a745;
  }
  
  .card-header.processed {
    background-color: #007bff;
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
  
  .expression-content {
    background-color: white;
    padding: 16px;
    border-radius: 4px;
  }
  
  .expression-stats {
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
  
  .stat-value.upregulated {
    color: #dc3545;
  }
  
  .stat-value.downregulated {
    color: #28a745;
  }
  
  .stat-value.significant {
    color: #dc3545;
  }
  
  .psi-content {
    background-color: white;
    padding: 16px;
    border-radius: 4px;
  }
  
  .psi-comparison {
    margin-bottom: 20px;
  }
  
  .psi-item {
    margin-bottom: 16px;
  }
  
  .psi-item:last-child {
    margin-bottom: 0;
  }
  
  .psi-label {
    display: block;
    margin-bottom: 6px;
    font-weight: 500;
    color: #495057;
    font-size: 14px;
  }
  
  .psi-bar {
    position: relative;
    height: 24px;
    background-color: #f8f9fa;
    border-radius: 12px;
    border: 1px solid #dee2e6;
    overflow: hidden;
  }
  
  .psi-fill {
    height: 100%;
    border-radius: 12px;
    transition: width 0.3s ease;
  }
  
  .psi-fill.original {
    background: linear-gradient(90deg, #28a745, #20c997);
  }
  
  .psi-fill.processed {
    background: linear-gradient(90deg, #007bff, #6f42c1);
  }
  
  .psi-value {
    position: absolute;
    top: 50%;
    right: 8px;
    transform: translateY(-50%);
    font-size: 11px;
    font-weight: bold;
    color: #495057;
  }
  
  .delta-summary {
    text-align: center;
    padding: 16px;
    background-color: #f8f9fa;
    border-radius: 4px;
    border: 1px solid #e9ecef;
  }
  
  .delta-value {
    font-size: 1.25rem;
    font-weight: bold;
    margin-bottom: 8px;
  }
  
  .delta-value.high-change {
    color: #dc3545;
  }
  
  .delta-value.medium-change {
    color: #fd7e14;
  }
  
  .delta-value.low-change {
    color: #6c757d;
  }
  
  .quality-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 12px;
    background-color: white;
    padding: 16px;
    border-radius: 4px;
  }
  
  .quality-item {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 8px 0;
    border-bottom: 1px solid #f1f3f4;
  }
  
  .quality-label {
    color: #6c757d;
    font-weight: 500;
  }
  
  .quality-value {
    font-weight: bold;
    color: #495057;
  }
  
  .dialog-footer {
    text-align: right;
  }
  
  /* 基因结构可视化样式 */
  .gene-structure-section {
    margin-bottom: 24px;
    padding: 16px;
    background-color: #f8f9fa;
    border-radius: 6px;
    border: 1px solid #e9ecef;
  }
  
  .structure-info {
    margin-bottom: 20px;
  }
  
  .structure-header {
    display: flex;
    align-items: center;
    gap: 16px;
    margin-bottom: 12px;
  }
  
  .gene-name {
    font-size: 18px;
    font-weight: bold;
    color: #2c3e50;
  }
  
  .coordinates {
    font-family: monospace;
    color: #6c757d;
    font-size: 14px;
  }
  
  .splicing-events {
    display: flex;
    align-items: center;
    gap: 8px;
    flex-wrap: wrap;
  }
  
  .events-label {
    font-weight: 500;
    color: #495057;
  }
  
  .structure-diagram {
    background-color: white;
    padding: 20px;
    border-radius: 6px;
    border: 1px solid #dee2e6;
    margin-bottom: 16px;
  }
  
  .diagram-container {
    margin-bottom: 20px;
  }
  
  .gene-track,
  .transcript-track {
    margin-bottom: 16px;
  }
  
  .track-label {
    font-weight: 500;
    color: #495057;
    margin-bottom: 8px;
    font-size: 14px;
  }
  
  .gene-line,
  .transcript-line {
    display: flex;
    align-items: center;
    height: 40px;
    background-color: #f8f9fa;
    border-radius: 4px;
    padding: 4px;
    position: relative;
  }
  
  .structure-element {
    height: 32px;
    display: flex;
    align-items: center;
    justify-content: center;
    border-radius: 4px;
    margin: 0 1px;
    cursor: pointer;
    transition: all 0.3s ease;
    position: relative;
  }
  
  .structure-element.exon {
    border: 2px solid #2c3e50;
  }
  
  .structure-element.intron {
    height: 8px;
    border-radius: 2px;
    margin-top: 12px;
  }
  
  .structure-element.skipped {
    opacity: 0.5;
    background-color: #ff9999 !important;
    border-style: dashed;
  }
  
  .structure-element:hover {
    transform: scale(1.05);
    z-index: 10;
  }
  
  .element-label {
    font-size: 11px;
    font-weight: bold;
    color: white;
    text-shadow: 0 1px 2px rgba(0,0,0,0.5);
  }
  
  .legend {
    display: flex;
    justify-content: center;
    gap: 20px;
    flex-wrap: wrap;
    padding: 12px;
    background-color: #f8f9fa;
    border-radius: 4px;
    border: 1px solid #e9ecef;
  }
  
  .legend-item {
    display: flex;
    align-items: center;
    gap: 6px;
    font-size: 12px;
    color: #495057;
  }
  
  .legend-color {
    width: 20px;
    height: 12px;
    border-radius: 2px;
    border: 1px solid #dee2e6;
  }
  
  .exon-color {
    background-color: #4ecdc4;
  }
  
  .intron-color {
    background-color: #95a5a6;
  }
  
  .affected-color {
    background-color: #ff6b6b;
  }
  
  .skipped-color {
    background-color: #ff9999;
    border-style: dashed;
  }
  
  .exon-details {
    background-color: white;
    padding: 16px;
    border-radius: 6px;
    border: 1px solid #dee2e6;
  }
  
  .exon-details h5 {
    margin-top: 0;
    margin-bottom: 12px;
    color: #495057;
  }
  
  .exon-table {
    font-size: 14px;
  }
  
  .high-inclusion {
    color: #28a745;
    font-weight: bold;
  }
  
  .medium-inclusion {
    color: #ffc107;
    font-weight: 500;
  }
  
  .low-inclusion {
    color: #dc3545;
    font-weight: bold;
  }
  
  :deep(.el-dialog) {
    margin-top: 5vh;
  }
  
  :deep(.el-dialog__body) {
    padding: 10px 20px;
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

.structure-element.exon {
  cursor: pointer;
}

.structure-element.exon:hover {
  transform: scale(1.1);
  box-shadow: 0 2px 8px rgba(0,0,0,0.2);
}
.skipped-color {
  background-color: #ff9999;
  border-style: dashed;
}

.exclusive-color {
  background-color: #e74c3c;
  border: 2px solid #c0392b;
}

.alt5-color {
  background-color: #f39c12;
  border: 2px solid #d68910;
}

.alt3-color {
  background-color: #9b59b6;
  border: 2px solid #8e44ad;
}

.retained-color {
  background-color: #34495e;
  border: 2px solid #2c3e50;
}
  </style>