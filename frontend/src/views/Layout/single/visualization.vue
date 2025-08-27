<template>
  <div class="splicing-visualization">
    <!-- 基因信息头部 -->
    <div class="gene-header">
      <div class="gene-title">
        <span class="gene-name">{{ actualIsoformData?.geneSymbol || 'Unknown Gene' }}</span>
        <span class="gene-coordinates">
          {{ actualIsoformData?.chromosome }}:{{ actualStructureData?.geneStart || actualIsoformData?.exonStart }}-{{ actualStructureData?.geneEnd || actualIsoformData?.exonEnd }}
        </span>
        <el-tag size="small">{{ actualIsoformData?.strand || '+' }}</el-tag>
      </div>
      
      <div class="event-info">
        <span class="event-label">Detected Splicing Event:</span>
        <el-tag 
          :type="getEventTypeColor(currentEventType)"
          size="small"
        >
          {{ currentEventType }} - {{ getEventTypeName(currentEventType) }}
        </el-tag>
      </div>
    </div>

    <!-- 可视化内容 -->
    <div class="visualization-content">
      
      <!-- SE (Skipped Exon) 可视化 -->
      <div v-if="currentEventType === 'SE'" class="se-visualization">
        <!-- 原始转录本 -->
        <div class="transcript-section original-section">
          <div class="section-title">📄 Original Transcript (All exons included):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in originalTranscript" :key="`orig-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block original"
                :class="{ 'affected': element.isAffected }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>

        <!-- 跳跃指示器 -->
        <div class="skip-indicator" v-if="skippedExon">
          <svg :width="indicatorWidth" height="80" class="skip-svg">
            <path 
              :d="skipPath"
              stroke="#E74C3C" 
              stroke-width="3" 
              fill="none" 
              stroke-dasharray="8,4"
              marker-end="url(#arrowhead)"
            />
            <defs>
              <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="10" refY="3.5" orient="auto">
                <polygon points="0 0, 10 3.5, 0 7" fill="#E74C3C" />
              </marker>
            </defs>
            <text :x="indicatorWidth / 2" y="30" text-anchor="middle" class="skip-label">
              Exon {{ skippedExon.name }} Skipping
            </text>
            <text :x="indicatorWidth / 2" y="50" text-anchor="middle" class="skip-detail">
              {{ skippedExon.length || 'Unknown' }}bp excluded from mature mRNA
            </text>
          </svg>
        </div>

        <!-- 可变转录本 -->
        <div class="transcript-section alternative-section">
          <div class="section-title">📄 Alternative Transcript ({{ skippedExon?.name || 'Unknown exon' }} skipped):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in getAlternativeDisplay" :key="`alt-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block alternative"
                :class="{ 
                  'affected': element.isAffected,
                  'skipped': element.isSkipped
                }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
                <span v-if="element.isSkipped" class="skip-marker">×</span>
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>
      </div>

      <!-- MXE (Mutually Exclusive Exons) 可视化 -->
      <div v-else-if="currentEventType === 'MXE'" class="mxe-visualization">
        <div class="transcript-section original-section">
          <div class="section-title">📄 Original Transcript (Exon A included):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in originalTranscript" :key="`orig-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block original"
                :class="{ 'affected': element.isAffected, 'exon-a': element.name?.includes('a') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>

        <div class="exclusion-indicator">
          <div class="exclusion-text">
            Mutually Exclusive: {{ mutuallyExclusiveExons.exonA?.name || 'Exon A' }} OR {{ mutuallyExclusiveExons.exonB?.name || 'Exon B' }}
          </div>
        </div>

        <div class="transcript-section alternative-section">
          <div class="section-title">📄 Alternative Transcript (Exon B included):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in alternativeTranscript" :key="`alt-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block alternative"
                :class="{ 'affected': element.isAffected, 'exon-b': element.name?.includes('b') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>
      </div>

      <!-- A5SS (Alternative 5' Splice Site) 可视化 -->
      <div v-else-if="currentEventType === 'A5SS'" class="a5ss-visualization">
        <div class="transcript-section original-section">
          <div class="section-title">📄 Original Transcript (Longer exon):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in originalTranscript" :key="`orig-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block original"
                :class="{ 'affected': element.isAffected }"
                :style="{ width: getExonWidth(element, 'original') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
                <span v-if="element.hasAlternativeForm" class="length-indicator">
                  {{ element.length }}bp
                </span>
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>

        <div class="splice-site-indicator">
          <div class="site-change">
            5' Splice Site Alternative → {{ getAffectedExonName() }} 
            ({{ getOriginalLength() }}bp → {{ getAlternativeLength() }}bp)
          </div>
        </div>

        <div class="transcript-section alternative-section">
          <div class="section-title">📄 Alternative Transcript (Shorter exon):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in alternativeTranscript" :key="`alt-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block alternative"
                :class="{ 'affected': element.isAffected, 'shortened': element.alternativeForm === 'short' }"
                :style="{ width: getExonWidth(element, 'alternative') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
                <span v-if="element.alternativeForm" class="length-indicator">
                  {{ element.length }}bp
                </span>
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>
      </div>

      <!-- A3SS (Alternative 3' Splice Site) 可视化 -->
      <div v-else-if="currentEventType === 'A3SS'" class="a3ss-visualization">
        <div class="transcript-section original-section">
          <div class="section-title">📄 Original Transcript ({{ getOriginalTranscriptTitle() }}):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in originalTranscript" :key="`orig-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block original"
                :class="{ 'affected': element.isAffected }"
                :style="{ width: getExonWidth(element, 'original') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
                <span v-if="element.hasAlternativeForm" class="length-indicator">
                  {{ element.length }}bp
                </span>
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>

        <div class="splice-site-indicator">
          <div class="site-change">
            3' Splice Site Alternative → {{ getAffectedExonName() }}
            ({{ getOriginalLength() }}bp → {{ getAlternativeLength() }}bp)
          </div>
        </div>

        <div class="transcript-section alternative-section">
          <div class="section-title">📄 Alternative Transcript ({{ getAlternativeTranscriptTitle() }}):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in alternativeTranscript" :key="`alt-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block alternative"
                :class="{ 
                  'affected': element.isAffected, 
                  'shortened': element.alternativeForm === 'short',
                  'lengthened': element.alternativeForm === 'long'
                }"
                :style="{ width: getExonWidth(element, 'alternative') }"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
                <span v-if="element.alternativeForm" class="length-indicator">
                  {{ element.length }}bp
                </span>
              </div>
              <div v-else class="intron-connector"></div>
            </template>
          </div>
        </div>
      </div>

      <!-- RI (Retained Intron) 可视化 -->
      <div v-else-if="currentEventType === 'RI'" class="ri-visualization">
        <div class="transcript-section original-section">
          <div class="section-title">📄 Original Transcript (Intron spliced out):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in originalTranscript" :key="`orig-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block original"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
              </div>
              <div v-else class="intron-connector normal"></div>
            </template>
          </div>
        </div>

        <div class="retention-indicator">
          <div class="retention-text">
            {{ retainedIntron?.name || 'Intron' }} Retained in Mature mRNA
          </div>
        </div>

        <div class="transcript-section alternative-section">
          <div class="section-title">📄 Alternative Transcript (Intron retained):</div>
          <div class="transcript-track">
            <template v-for="(element, index) in alternativeTranscript" :key="`alt-${index}`">
              <div 
                v-if="element.type === 'exon'"
                class="exon-block alternative"
                @click="handleExonClick(element)"
              >
                {{ element.name }}
              </div>
              <div 
                v-else-if="element.type === 'intron' && element.retained"
                class="intron-connector retained"
              >
                {{ element.name }}
              </div>
              <div v-else class="intron-connector normal"></div>
            </template>
          </div>
        </div>
      </div>

      <!-- 未知类型的fallback -->
      <div v-else class="unknown-visualization">
        <el-alert
          title="Unknown Splicing Event"
          :description="`Event type '${currentEventType}' is not yet supported.`"
          type="warning"
          :closable="false"
        />
      </div>

    </div>

    <!-- 分析摘要 -->
    <div class="analysis-summary">
      <h5>📊 {{ getEventTypeName(currentEventType) }} Analysis Summary:</h5>
      <div class="summary-grid">
        <div class="summary-item">
          <strong>Inclusion Level Change:</strong> 
          {{ actualIsoformData?.originalInc || 'N/A' }} → {{ actualIsoformData?.compareInc || 'N/A' }}
          <span v-if="inclusionChange" :class="getChangeClass(inclusionChange)">
            ({{ getChangeDirection(inclusionChange) }})
          </span>
        </div>
        <div class="summary-item">
          <strong>Statistical Significance:</strong>
          <el-tag :type="actualIsoformData?.isSignificant ? 'success' : 'info'" size="small">
            {{ actualIsoformData?.isSignificant ? 'Significant' : 'Not Significant' }}
          </el-tag>
        </div>
        <div class="summary-item">
          <strong>Confidence:</strong>
          <span :class="getConfidenceClass(actualIsoformData?.confidence)">
            {{ actualIsoformData?.confidence ? (actualIsoformData.confidence * 100).toFixed(1) + '%' : 'N/A' }}
          </span>
        </div>
        <div class="summary-item">
          <strong>Functional Impact:</strong>
          {{ getFunctionalImpact(currentEventType) }}
        </div>
      </div>
    </div>

    <!-- 图例 -->
    <div class="legend-panel">
      <div class="legend-title">Legend:</div>
      <div class="legend-items">
        <div class="legend-item">
          <div class="legend-color normal-exon"></div>
          <span>Normal Exon</span>
        </div>
        <div class="legend-item">
          <div class="legend-color affected-exon"></div>
          <span>Affected Exon</span>
        </div>
        <div class="legend-item">
          <div class="legend-color intron"></div>
          <span>Intron</span>
        </div>
        <div v-for="item in eventSpecificLegend" :key="item.label" class="legend-item">
          <div :class="['legend-color', item.class]"></div>
          <span>{{ item.label }}</span>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { computed } from 'vue'

// Props - 修改为非必需，提供更好的向后兼容性
const props = defineProps({
  isoformData: {
    type: Object,
    required: false,
    default: () => ({})
  },
  structureData: {
    type: Object,
    required: false,
    default: () => ({})
  },
  // 新增：支持直接传入完整的行数据
  rowData: {
    type: Object,
    required: false,
    default: () => ({})
  }
})

// Emits
const emit = defineEmits(['exon-click'])

// 计算实际使用的数据 - 智能适配不同的数据传入方式
const actualIsoformData = computed(() => {
  // 优先使用 isoformData，如果没有则使用 rowData
  if (props.isoformData && Object.keys(props.isoformData).length > 0) {
    return props.isoformData
  }
  
  if (props.rowData && Object.keys(props.rowData).length > 0) {
    return props.rowData
  }
  
  return {}
})

const actualStructureData = computed(() => {
  // 优先使用 structureData，如果没有则从 rowData.structureData 获取
  if (props.structureData && Object.keys(props.structureData).length > 0) {
    return props.structureData
  }
  
  if (props.rowData?.structureData) {
    return props.rowData.structureData
  }
  
  return {}
})

// 计算属性 - 使用实际数据
const currentEventType = computed(() => {
  return actualIsoformData.value?.eventType || 'SE'
})

const originalTranscript = computed(() => {
  return actualStructureData.value?.originalTranscript || []
})

const alternativeTranscript = computed(() => {
  return actualStructureData.value?.alternativeTranscript || []
})

const skippedExon = computed(() => {
  // 通过比较原始和可变转录本找到被跳过的外显子
  const originalExons = originalTranscript.value.filter(e => e.type === 'exon')
  const alternativeExons = alternativeTranscript.value.filter(e => e.type === 'exon')
  
  // 方法1: 查找在可变转录本中标记为跳过的外显子
  let skipped = alternativeTranscript.value.find(exon => 
    exon.type === 'exon' && (!exon.included || exon.skipped)
  )
  
  // 方法2: 如果方法1没找到，比较两个转录本的差异
  if (!skipped && originalExons.length > alternativeExons.length) {
    const alternativeExonNames = alternativeExons.map(e => e.name)
    skipped = originalExons.find(exon => !alternativeExonNames.includes(exon.name))
  }
  
  // 方法3: 如果还没找到，查找isAffected的外显子
  if (!skipped) {
    skipped = originalExons.find(exon => exon.isAffected)
  }
  
  console.log('🔍 Skipped exon analysis:', {
    originalExons: originalExons.map(e => e.name),
    alternativeExons: alternativeExons.map(e => e.name),
    foundSkipped: skipped?.name
  })
  
  return skipped
})

const inclusionChange = computed(() => {
  const original = parseFloat(actualIsoformData.value?.originalInc || 0)
  const compare = parseFloat(actualIsoformData.value?.compareInc || 0)
  return compare - original
})

// 动态计算跳跃指示器宽度
const indicatorWidth = computed(() => {
  const exonCount = originalTranscript.value.filter(e => e.type === 'exon').length
  return Math.max(300, exonCount * 120)
})

// 动态计算跳跃路径
const skipPath = computed(() => {
  const startX = indicatorWidth.value * 0.2
  const endX = indicatorWidth.value * 0.8
  const midX = (startX + endX) / 2
  const peakY = 20
  
  return `M ${startX} 60 Q ${midX} ${peakY} ${endX} 60`
})

// 获取可变转录本的正确显示
const getAlternativeDisplay = computed(() => {
  if (!skippedExon.value) return alternativeTranscript.value
  
  // 基于原始转录本，标记被跳过的外显子
  return originalTranscript.value.map(element => {
    if (element.type === 'exon' && element.name === skippedExon.value.name) {
      return {
        ...element,
        isSkipped: true,
        included: false
      }
    }
    return element
  })
})

// 获取互斥外显子信息（用于MXE）
const mutuallyExclusiveExons = computed(() => {
  if (currentEventType.value !== 'MXE') return { exonA: null, exonB: null }
  
  const originalExons = originalTranscript.value.filter(e => e.type === 'exon' && e.isAffected)
  const alternativeExons = alternativeTranscript.value.filter(e => e.type === 'exon' && e.isAffected)
  
  return {
    exonA: originalExons[0] || null,
    exonB: alternativeExons[0] || null
  }
})

// 获取长度变化的外显子信息（用于A5SS/A3SS）
const lengthChangedExon = computed(() => {
  if (!['A5SS', 'A3SS'].includes(currentEventType.value)) return null
  
  return originalTranscript.value.find(e => e.type === 'exon' && e.isAffected) || null
})

// 获取保留的内含子信息（用于RI）
const retainedIntron = computed(() => {
  if (currentEventType.value !== 'RI') return null
  
  return alternativeTranscript.value.find(e => e.type === 'intron' && e.retained) || null
})
const eventSpecificLegend = computed(() => {
  const legends = {
    'SE': [{ label: 'Skipped Exon', class: 'skipped-exon' }],
    'MXE': [
      { label: 'Alternative Exon A', class: 'exon-a' },
      { label: 'Alternative Exon B', class: 'exon-b' }
    ],
    'A5SS': [{ label: 'Shorter Exon', class: 'shorter-exon' }],
    'A3SS': [{ label: 'Longer Exon', class: 'longer-exon' }],
    'RI': [{ label: 'Retained Intron', class: 'retained-intron' }]
  }
  return legends[currentEventType.value] || []
})

// 方法
const getEventTypeColor = (eventType) => {
  const colors = {
    'SE': 'danger', 'MXE': 'success', 'A5SS': 'warning', 'A3SS': 'warning', 'RI': 'info'
  }
  return colors[eventType] || 'default'
}

const getEventTypeName = (eventType) => {
  const names = {
    'SE': 'Skipped Exon',
    'MXE': 'Mutually Exclusive Exons',
    'A5SS': 'Alternative 5\' Splice Site',
    'A3SS': 'Alternative 3\' Splice Site',
    'RI': 'Retained Intron'
  }
  return names[eventType] || eventType
}

const getConfidenceClass = (confidence) => {
  if (!confidence) return ''
  if (confidence >= 0.8) return 'high-confidence'
  if (confidence >= 0.6) return 'medium-confidence'
  return 'low-confidence'
}

const getChangeClass = (diff) => {
  if (Math.abs(diff) > 0.2) return 'high-change'
  if (Math.abs(diff) > 0.1) return 'medium-change'
  return 'low-change'
}

const getChangeDirection = (diff) => {
  if (diff > 0.05) return `+${diff.toFixed(3)}`
  if (diff < -0.05) return `${diff.toFixed(3)}`
  return `${diff.toFixed(3)} (minimal)`
}

const getFunctionalImpact = (eventType) => {
  const impacts = {
    'SE': 'Protein domain removal/alteration',
    'MXE': 'Alternative functional domains',
    'A5SS': 'Partial domain modification',
    'A3SS': 'Extended/truncated domain',
    'RI': 'Premature stop codon likely'
  }
  return impacts[eventType] || 'Unknown impact'
}

const handleExonClick = (exonData) => {
  emit('exon-click', exonData)
}

// 新增：处理外显子宽度计算
const getExonWidth = (element, transcriptType) => {
  if (!element.isAffected) return '80px'
  
  // 对于A5SS/A3SS事件，根据长度动态计算宽度
  if (['A5SS', 'A3SS'].includes(currentEventType.value)) {
    const baseWidth = 60
    const lengthRatio = transcriptType === 'original' 
      ? (element.length || 100) / 100 
      : (element.alternativeLength || element.length || 100) / 100
    
    return Math.max(baseWidth, baseWidth * lengthRatio) + 'px'
  }
  
  return element.isAffected ? '100px' : '80px'
}

// 新增：获取受影响外显子的名称
const getAffectedExonName = () => {
  const affected = originalTranscript.value.find(e => e.type === 'exon' && e.isAffected)
  return affected?.name || 'Exon'
}

// 新增：获取原始长度
const getOriginalLength = () => {
  const affected = originalTranscript.value.find(e => e.type === 'exon' && e.isAffected)
  return affected?.length || 'Unknown'
}

// 新增：获取可变长度
const getAlternativeLength = () => {
  const affected = alternativeTranscript.value.find(e => e.type === 'exon' && e.isAffected)
  return affected?.length || affected?.alternativeLength || 'Unknown'
}

// 新增：获取原始转录本标题
const getOriginalTranscriptTitle = () => {
  if (currentEventType.value === 'A3SS') {
    const originalLen = getOriginalLength()
    const altLen = getAlternativeLength()
    if (originalLen !== 'Unknown' && altLen !== 'Unknown') {
      return originalLen > altLen ? 'Longer exon' : 'Shorter exon'
    }
  }
  return 'Original exon'
}

// 新增：获取可变转录本标题
const getAlternativeTranscriptTitle = () => {
  if (currentEventType.value === 'A3SS') {
    const originalLen = getOriginalLength()
    const altLen = getAlternativeLength()
    if (originalLen !== 'Unknown' && altLen !== 'Unknown') {
      return originalLen > altLen ? 'Shorter exon' : 'Longer exon'
    }
  }
  return 'Alternative exon'
}
</script>

<style scoped>
.splicing-visualization {
  background-color: #f8f9fa;
  border-radius: 8px;
  padding: 20px;
  border: 1px solid #e9ecef;
}

.gene-header {
  margin-bottom: 24px;
  padding-bottom: 16px;
  border-bottom: 2px solid #dee2e6;
}

.gene-title {
  display: flex;
  align-items: center;
  gap: 12px;
  margin-bottom: 12px;
  flex-wrap: wrap;
}

.gene-name {
  font-size: 20px;
  font-weight: bold;
  color: #2c3e50;
}

.gene-coordinates {
  font-family: monospace;
  color: #6c757d;
  font-size: 14px;
  background-color: #f1f3f4;
  padding: 4px 8px;
  border-radius: 4px;
}

.event-info {
  display: flex;
  align-items: center;
  gap: 8px;
  flex-wrap: wrap;
}

.event-label {
  font-weight: 500;
  color: #495057;
  font-size: 14px;
}

.visualization-content {
  background-color: white;
  border-radius: 8px;
  padding: 24px;
  border: 1px solid #dee2e6;
  margin-bottom: 20px;
  min-height: 300px;
}

.transcript-section {
  margin-bottom: 24px;
  padding: 16px;
  border-radius: 8px;
  border: 1px solid #dee2e6;
}

.original-section {
  background-color: #f8fff8;
  border-left: 4px solid #28a745;
}

.alternative-section {
  background-color: #f8f9ff;
  border-left: 4px solid #007bff;
}

.section-title {
  font-weight: 600;
  color: #495057;
  margin-bottom: 16px;
  font-size: 16px;
}

.transcript-track {
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 8px;
  padding: 12px;
  background-color: white;
  border-radius: 6px;
  min-height: 60px;
  justify-content: center;
}

.exon-block {
  height: 36px;
  width: 80px;
  display: flex;
  align-items: center;
  justify-content: center;
  border-radius: 6px;
  cursor: pointer;
  transition: all 0.3s ease;
  position: relative;
  font-size: 12px;
  font-weight: bold;
  color: white;
  text-shadow: 0 1px 2px rgba(0,0,0,0.3);
  border: 2px solid #2c3e50;
}

.exon-block.original {
  background-color: #4ECDC4;
}

.exon-block.alternative {
  background-color: #3498DB;
}

.exon-block.affected {
  background-color: #E74C3C !important;
}

.exon-block.skipped {
  opacity: 0.5;
  border-style: dashed;
  background-color: #BDC3C7 !important;
}

.exon-block.exon-a {
  background-color: #9B59B6 !important;
}

.exon-block.exon-b {
  background-color: #8E44AD !important;
}

.exon-block.shortened {
  background-color: #F39C12 !important;
}

.exon-block.lengthened {
  background-color: #E67E22 !important;
}

.length-indicator {
  position: absolute;
  bottom: -18px;
  left: 50%;
  transform: translateX(-50%);
  font-size: 9px;
  color: #666;
  background-color: rgba(255,255,255,0.9);
  padding: 1px 3px;
  border-radius: 2px;
  white-space: nowrap;
}

.exon-block {
  position: relative; /* 确保length-indicator定位正确 */
}

.skip-marker {
  position: absolute;
  top: -4px;
  right: -4px;
  background-color: #E74C3C;
  color: white;
  border-radius: 50%;
  width: 16px;
  height: 16px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 10px;
}

.intron-connector {
  height: 4px;
  width: 30px;
  background-color: #95A5A6;
  border-radius: 2px;
  margin: 0 4px;
}

.intron-connector.retained {
  height: 24px;
  width: 60px;
  background-color: #27AE60;
  border: 2px dashed #2c3e50;
  border-radius: 4px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 10px;
  font-weight: bold;
  color: white;
}

.skip-indicator, .exclusion-indicator, .splice-site-indicator, .retention-indicator {
  display: flex;
  justify-content: center;
  margin: 20px 0;
  padding: 16px;
  border-radius: 8px;
  border: 1px solid #fed7d7;
}

.skip-indicator {
  background-color: #fff5f5;
}

.exclusion-indicator {
  background-color: #f0fff4;
  border-color: #c6f6d5;
}

.splice-site-indicator {
  background-color: #fffbf0;
  border-color: #ffeaa7;
}

.retention-indicator {
  background-color: #f0f8ff;
  border-color: #bee5eb;
}

.skip-svg {
  overflow: visible;
}

.skip-label {
  font-size: 14px;
  font-weight: bold;
  fill: #E74C3C;
}

.exclusion-text, .site-change, .retention-text {
  font-weight: 600;
  color: #495057;
  text-align: center;
}

.analysis-summary {
  background-color: #fffbf0;
  padding: 16px;
  border-radius: 6px;
  border: 1px solid #ffeaa7;
  border-left: 4px solid #fdcb6e;
  margin-bottom: 20px;
}

.analysis-summary h5 {
  margin: 0 0 12px 0;
  color: #2d3436;
  font-size: 16px;
}

.summary-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 12px;
}

.summary-item {
  font-size: 14px;
  color: #2d3436;
  line-height: 1.5;
}

.summary-item strong {
  color: #636e72;
  margin-right: 8px;
}

.high-confidence { color: #28a745; font-weight: bold; }
.medium-confidence { color: #fd7e14; font-weight: 500; }
.low-confidence { color: #6c757d; }

.high-change { color: #dc3545; font-weight: bold; }
.medium-change { color: #fd7e14; font-weight: 500; }
.low-change { color: #6c757d; }

.legend-panel {
  display: flex;
  align-items: center;
  gap: 16px;
  padding: 12px;
  background-color: #f8f9fa;
  border-radius: 4px;
  border: 1px solid #e9ecef;
  flex-wrap: wrap;
}

.legend-title {
  font-weight: 600;
  color: #495057;
  font-size: 14px;
}

.legend-items {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 6px;
}

.legend-color {
  width: 24px;
  height: 16px;
  border-radius: 2px;
  border: 1px solid #dee2e6;
}

.normal-exon { background-color: #4ECDC4; }
.affected-exon { background-color: #E74C3C; }
.intron { background-color: #95A5A6; height: 4px; width: 20px; }
.skipped-exon { background-color: #BDC3C7; opacity: 0.5; border-style: dashed; }
.exon-a { background-color: #9B59B6; }
.exon-b { background-color: #8E44AD; }
.shorter-exon { background-color: #F39C12; width: 20px; }
.longer-exon { background-color: #E67E22; width: 32px; }
.retained-intron { background-color: #27AE60; border: 2px dashed #2c3e50; }

.legend-item span {
  font-size: 12px;
  color: #6c757d;
}

/* 响应式设计 */
@media (max-width: 768px) {
  .splicing-visualization {
    padding: 16px;
  }
  
  .gene-title {
    flex-direction: column;
    align-items: flex-start;
    gap: 8px;
  }
  
  .visualization-content {
    padding: 16px;
  }
  
  .transcript-track {
    justify-content: flex-start;
  }
  
  .summary-grid {
    grid-template-columns: 1fr;
  }
  
  .legend-panel {
    flex-direction: column;
    align-items: flex-start;
  }
}
</style>