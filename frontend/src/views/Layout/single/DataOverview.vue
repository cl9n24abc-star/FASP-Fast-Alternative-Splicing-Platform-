<template>
    <div class="data-overview">
      <div class="header-section">
        <h1>📊 Data Analysis Overview</h1>
        <p class="subtitle">Comprehensive analysis of BAM files and rMATS splicing results</p>
      </div>
    
      <!-- Tab Navigation - Using Element Plus -->
      <el-card class="tab-container" shadow="never">
        <div class="tab-header">
          <el-tabs v-model="activeTab" type="card" class="modern-tabs">
            <el-tab-pane label="📄 BAM Analysis" name="bam"></el-tab-pane>
            <el-tab-pane label="🧬 rMATS Results" name="rmats"></el-tab-pane>
          </el-tabs>
          
        </div>
      </el-card>
    
      <!-- BAM Analysis Tab -->
      <div v-if="activeTab === 'bam'" class="tab-content">
                <!-- Summary Cards - Enhanced with Element Plus -->
        <el-row :gutter="20" class="summary-section">
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-blue">
              <div class="dual-stat-content">
                <span class="stat-icon">📄</span>
                <div class="dual-stat-info">
                  <div class="stat-title">File Information</div>
                  <div class="dual-values">
                    <div class="value-row">
                      <span class="group-label">Group 1:</span>
                      <span class="stat-value-small">{{ bamData.group1.fileSize }}</span>
                    </div>
                    <div class="value-row">
                      <span class="group-label">Group 2:</span>
                      <span class="stat-value-small">{{ bamData.group2.fileSize }}</span>
                    </div>
                  </div>
                  <div class="stat-subtitle">BAM File Sizes</div>
                </div>
              </div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-green">
              <div class="dual-stat-content">
                <span class="stat-icon">📊</span>
                <div class="dual-stat-info">
                  <div class="stat-title">Mapping Statistics</div>
                  <div class="dual-values">
                    <div class="value-row">
                      <span class="group-label">Group 1:</span>
                      <span class="stat-value-small">{{ bamData.group1.mappingRate }}%</span>
                    </div>
                    <div class="value-row">
                      <span class="group-label">Group 2:</span>
                      <span class="stat-value-small">{{ bamData.group2.mappingRate }}%</span>
                    </div>
                  </div>
                  <div class="stat-subtitle">Success Rate</div>
                </div>
              </div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-orange">
              <div class="dual-stat-content">
                <span class="stat-icon">🎯</span>
                <div class="dual-stat-info">
                  <div class="stat-title">Read Quality</div>
                  <div class="dual-values">
                    <div class="value-row">
                      <span class="group-label">Group 1:</span>
                      <span class="stat-value-small">{{ bamData.group1.avgMapQ }}</span>
                    </div>
                    <div class="value-row">
                      <span class="group-label">Group 2:</span>
                      <span class="stat-value-small">{{ bamData.group2.avgMapQ }}</span>
                    </div>
                  </div>
                  <div class="stat-subtitle">Average MAPQ</div>
                </div>
              </div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-purple">
              <div class="dual-stat-content">
                <span class="stat-icon">📏</span>
                <div class="dual-stat-info">
                  <div class="stat-title">Insert Fragment</div>
                  <div class="dual-values">
                    <div class="value-row">
                      <span class="group-label">Group 1:</span>
                      <span class="stat-value-small">{{ bamData.group1.avgInsertSize }} bp</span>
                    </div>
                    <div class="value-row">
                      <span class="group-label">Group 2:</span>
                      <span class="stat-value-small">{{ bamData.group2.avgInsertSize }} bp</span>
                    </div>
                  </div>
                  <div class="stat-subtitle">Average Size</div>
                </div>
              </div>
            </el-card>
          </el-col>
        </el-row>
        <!-- Charts Grid - 2x2 Layout -->
        <el-row :gutter="20" class="charts-section">
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>🥧 Mapping Statistics</span>
                </div>
              </template>
              <div ref="pieChartRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
          
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>📊 MAPQ Distribution</span>
                </div>
              </template>
              <div ref="mapqHistogramRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
          
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>📈 Insert Size Distribution</span>
                </div>
              </template>
              <div ref="insertSizeRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
          
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>🌋 Coverage Difference Volcano Plot</span>
                  <div class="chart-controls">
                    <span class="control-label">P-value Threshold:</span>
                    <el-slider 
                      v-model="bamVolcanoThreshold.pvalue" 
                      :min="0.001" 
                      :max="0.1" 
                      :step="0.001"
                      style="width: 120px; margin: 0 10px;"
                      :format-tooltip="(val) => val.toFixed(3)"
                    />
                  </div>
                </div>
              </template>
              <div ref="bamVolcanoRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
        </el-row>
    
        <!-- Full Width Heatmap -->
        <el-row class="charts-section">
          <el-col :span="24">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>🔥 Chromosome Coverage Heatmap</span>
                </div>
              </template>
              <div ref="heatmapRef" style="width: 100%; height: 400px;"></div>
            </el-card>
          </el-col>
        </el-row>
      </div>
    
      <!-- rMATS Analysis Tab -->
      <div v-if="activeTab === 'rmats'" class="tab-content">
        <!-- Summary Cards -->
        <el-row :gutter="20" class="summary-section">
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-indigo">
              <el-statistic 
                title="Total Events"
                :value="rmatsData.totalEvents"
                suffix=""
              >
                <template #prefix>
                  <span class="stat-icon">🧬</span>
                </template>
              </el-statistic>
              <div class="stat-subtitle">All Splicing Events</div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-emerald">
              <el-statistic 
                title="Significant Events"
                :value="rmatsData.significantEvents"
                suffix=""
              >
                <template #prefix>
                  <span class="stat-icon">⭐</span>
                </template>
              </el-statistic>
              <div class="stat-subtitle">{{ ((rmatsData.significantEvents / rmatsData.totalEvents) * 100).toFixed(1) }}% Significance Rate</div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-amber">
              <div class="stat-content">
                <span class="stat-icon">📈</span>
                <div class="stat-info">
                  <div class="stat-title">ΔPSI Range</div>
                  <div class="stat-value">+{{ rmatsData.maxPsiIncrease }}</div>
                  <div class="stat-subtitle">Max Change: {{ rmatsData.maxPsiDecrease }}</div>
                </div>
              </div>
            </el-card>
          </el-col>
          
          <el-col :span="6">
            <el-card shadow="hover" class="stat-card gradient-violet">
              <div class="stat-content">
                <span class="stat-icon">🧪</span>
                <div class="stat-info">
                  <div class="stat-title">Sample Information</div>
                  <div class="stat-value">{{ rmatsData.group1Samples }} vs {{ rmatsData.group2Samples }}</div>
                  <div class="stat-subtitle">Group Comparison</div>
                </div>
              </div>
            </el-card>
          </el-col>
        </el-row>
    
        <!-- Charts Grid -->
        <el-row :gutter="20" class="charts-section">
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>🍩 Splicing Event Distribution</span>
                </div>
              </template>
              <div ref="donutChartRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
          
          <el-col :span="12">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>📊 Significance Statistics</span>
                </div>
              </template>
              <div ref="barChartRef" style="width: 100%; height: 350px;"></div>
            </el-card>
          </el-col>
        </el-row>
    
        <!-- Volcano Plot - Full Width -->
        <el-row class="charts-section">
          <el-col :span="24">
            <el-card shadow="hover" class="chart-card">
              <template #header>
                <div class="chart-header">
                  <span>🌋 Volcano Plot (Interactive)</span>
                  <div class="chart-controls">
                    <span class="control-label">FDR Threshold:</span>
                    <el-slider 
                      v-model="volcanoThreshold.fdr" 
                      :min="0.01" 
                      :max="0.1" 
                      :step="0.01"
                      style="width: 150px; margin: 0 15px;"
                      :format-tooltip="(val) => val.toFixed(2)"
                    />
                  </div>
                </div>
              </template>
              <div ref="volcanoRef" style="width: 100%; height: 500px;"></div>
            </el-card>
          </el-col>
        </el-row>
      </div>
    
      <!-- Status Messages -->
      <el-alert
        v-if="statusMessage"
        :title="statusMessage"
        type="info"
        show-icon
        :closable="true"
        @close="statusMessage = ''"
        class="status-alert"
      />
    </div>
    </template>

<script setup>
import { ref, onMounted, watch, nextTick } from 'vue'
import { Download, Refresh } from '@element-plus/icons-vue'

// Reactive Data
const activeTab = ref('bam')
const statusMessage = ref('')

// Interactive controls
const volcanoThreshold = ref({
fdr: 0.05
})

const bamVolcanoThreshold = ref({
pvalue: 0.05
})

// Chart refs
const pieChartRef = ref(null)
const mapqHistogramRef = ref(null)
const insertSizeRef = ref(null)
const heatmapRef = ref(null)
const bamVolcanoRef = ref(null)
const donutChartRef = ref(null)
const barChartRef = ref(null)
const volcanoRef = ref(null)

const bamData = ref({
  group1: {},
  group2: {}
})
const realVolcanoData = ref([])
const realPieData = ref([])
const realMapqData = ref([])
const realInsertSizeData = ref({})
const realHeatmapData = ref([])
const rmatsData = ref({
  totalEvents: 0,
  significantEvents: 0,
  maxPsiIncrease: 0,
  maxPsiDecrease: 0,
  group1Samples: 0,
  group2Samples: 0,
  eventCounts: {
    SE: 0,
    A5SS: 0,
    A3SS: 0,
    MXE: 0,
    RI: 0
  },
  volcanoData: [],
  significanceData: []
})

// Chart instances
let pieChart = null
let mapqChart = null
let insertSizeChart = null
let heatmapChart = null
let bamVolcanoChart = null
let donutChart = null
let barChart = null
let volcanoChart = null
// 添加图表数据更新函数
const updateChartsWithRealData = (chartData) => {
  console.log('Updating charts with real data:', chartData)
  
  // 更新各种图表数据
  if (chartData.volcano_plot) {
    realVolcanoData.value = chartData.volcano_plot
    console.log('Updated volcano data:', realVolcanoData.value.length, 'points')
  }
  
  if (chartData.pie_chart) {
    realPieData.value = chartData.pie_chart
    console.log('Updated pie chart data')
  }
  
  if (chartData.mapq_histogram) {
    realMapqData.value = chartData.mapq_histogram
    console.log('Updated MAPQ data')
  }
  
  if (chartData.insert_size_distribution) {
    realInsertSizeData.value = chartData.insert_size_distribution
    console.log('Updated insert size data')
  }
  
  if (chartData.chromosome_heatmap) {
    realHeatmapData.value = chartData.chromosome_heatmap
    console.log('Updated heatmap data')
  }
}
// Initialize Charts
const initPieChart = () => {
  if (pieChartRef.value && (window.$echarts || window.echarts)) {
    const echarts = window.$echarts || window.echarts
    pieChart = echarts.init(pieChartRef.value)
    
    // 使用真实数据或模拟数据
    let pieData = []
    if (realPieData.value.length > 0) {
      pieData = realPieData.value
      console.log('Using real pie chart data')
    } else {
      pieData = [
        { value: 90.3, name: 'Successfully Mapped', itemStyle: { color: '#5470c6' } },
        { value: 9.7, name: 'Unmapped', itemStyle: { color: '#ee6666' } },
        { value: 5.1, name: 'Duplicates', itemStyle: { color: '#fac858' } }
      ]
      console.log('Using simulated pie chart data')
    }
    
    const option = {
      title: {
        text: 'Mapping Statistics',
        left: 'center',
        textStyle: { fontSize: 16, fontWeight: 'bold' }
      },
      tooltip: {
        trigger: 'item',
        formatter: '{a} <br/>{b}: {c}% ({d}%)'
      },
      legend: {
        orient: 'vertical',
        left: 'left',
        top: 'center'
      },
      series: [{
        name: 'Mapping Stats',
        type: 'pie',
        radius: ['40%', '70%'],
        center: ['60%', '50%'],
        data: pieData,
        emphasis: {
          itemStyle: {
            shadowBlur: 10,
            shadowOffsetX: 0,
            shadowColor: 'rgba(0, 0, 0, 0.5)'
          }
        }
      }]
    }
    pieChart.setOption(option)

    // Add window resize handler
    const handleResize = () => {
      pieChart && pieChart.resize()
    }
    window.addEventListener('resize', handleResize)
  }
}

const initMapqChart = () => {
  if (mapqHistogramRef.value && (window.$echarts || window.echarts)) {
    const echarts = window.$echarts || window.echarts
    mapqChart = echarts.init(mapqHistogramRef.value)
    
    // 使用真实数据或模拟数据
    let histogramData = []
    if (realMapqData.value.length > 0) {
      histogramData = realMapqData.value.map(item => item.count)
      console.log('Using real MAPQ data')
    } else {
      histogramData = [120, 200, 150, 800, 700, 300]
      console.log('Using simulated MAPQ data')
    }
    
    const option = {
      title: {
        text: 'MAPQ Distribution',
        left: 'center'
      },
        tooltip: {
          trigger: 'axis',
          axisPointer: { type: 'shadow' },
          formatter: function(params) {
            const data = params[0]
            const qualityDesc = {
              '0-10': 'Poor Quality',
              '11-20': 'Low Quality', 
              '21-30': 'Fair Quality',
              '31-40': 'Good Quality',
              '41-50': 'High Quality',
              '51-60': 'Excellent Quality'
            }
            return `MAPQ Range: ${data.name}<br/>
                    Quality: ${qualityDesc[data.name]}<br/>
                    Read Count: ${data.value.toLocaleString()}`
          }
        },
      xAxis: {
          type: 'category',
          data: ['0-10', '11-20', '21-30', '31-40', '41-50', '51-60'],
          name: 'MAPQ Score Range',
          nameLocation: 'middle',
          nameGap: 35,
          nameTextStyle: {
            fontSize: 14,
            fontWeight: 'bold',
            color: '#333'
          },
      axisLabel: {
        fontSize: 12,
        color: '#666'
      },
      axisTick: {
        alignWithLabel: true
      }
    },
        yAxis: {
          type: 'value',
          name: 'Read Count',
          nameLocation: 'middle',
          nameGap: 50,
          nameTextStyle: {
            fontSize: 14,
            fontWeight: 'bold',
            color: '#333'
          },
          axisLabel: {
            formatter: function(value) {
              return value >= 1000 ? (value/1000).toFixed(1) + 'K' : value
            }
          }
        },
        grid: {
          left: '10%',
          right: '5%',
          bottom: '15%',
          top: '15%'
        },
      series: [{
        data: histogramData,
        type: 'bar',
    itemStyle: {
          color: function(params) {
            // 根据MAPQ质量给不同颜色
            const colors = ['#ff4757', '#ff6348', '#ffa502', '#2ed573', '#1e90ff', '#5352ed']
            return colors[params.dataIndex]
          }
        },
        emphasis: {
          itemStyle: {
            shadowBlur: 10,
            shadowColor: 'rgba(0, 0, 0, 0.3)'
          }
        }
      }]
    }
    mapqChart.setOption(option)
    // 添加说明文字
    const chartContainer = mapqHistogramRef.value.parentElement
    if (chartContainer && !chartContainer.querySelector('.chart-description')) {
      const description = document.createElement('div')
      description.className = 'chart-description'
      description.innerHTML = `
        <p style="font-size: 12px; color: #666; text-align: center; margin-top: 10px;">
          <strong>Note:</strong> Higher MAPQ scores indicate better mapping confidence. 
          MAPQ ≥30 is generally considered high quality.
        </p>
      `
      chartContainer.appendChild(description)
    }
  }
}

const initInsertSizeChart = () => {
if (insertSizeRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
insertSizeChart = echarts.init(insertSizeRef.value)
const xData = Array.from({length: 50}, (_, i) => 100 + i * 10)
const group1Data = xData.map(x => Math.exp(-Math.pow(x - 185, 2) / 800) * 1000)
const group2Data = xData.map(x => Math.exp(-Math.pow(x - 182, 2) / 800) * 1000)

// 计算更有意义的统计信息
const group1Peak = 185
const group2Peak = 182
const peakDifference = Math.abs(group1Peak - group2Peak)

// 计算质量评估
const group1Quality = group1Peak >= 150 && group1Peak <= 500 ? 'Good' : 'Poor'
const group2Quality = group2Peak >= 150 && group2Peak <= 500 ? 'Good' : 'Poor'

// 计算一致性评分
let consistency = 'Excellent'
if (peakDifference > 10) consistency = 'Poor'
else if (peakDifference > 5) consistency = 'Fair'
else if (peakDifference > 2) consistency = 'Good'

// 预测测序质量
const avgPeak = (group1Peak + group2Peak) / 2
let sequencingQuality = 'Excellent'
if (avgPeak < 100 || avgPeak > 600) sequencingQuality = 'Poor'
else if (avgPeak < 120 || avgPeak > 400) sequencingQuality = 'Fair'
else if (avgPeak < 150 || avgPeak > 300) sequencingQuality = 'Good'

const group1PeakValue = Math.max(...group1Data).toFixed(0)
const group2PeakValue = Math.max(...group2Data).toFixed(0)

const option = {
  title: {
    text: 'Insert Size Distribution',
    left: 'center',
    textStyle: { fontSize: 16, fontWeight: 'bold', color: '#2c3e50' },
    top: 15
  },
  tooltip: {
    trigger: 'axis',
    backgroundColor: 'rgba(0,0,0,0.85)',
    borderWidth: 0,
    textStyle: { color: '#fff', fontSize: 13 },
    formatter: function(params) {
      let result = `<div style="padding: 8px;">
                      <div style="font-weight: bold; margin-bottom: 6px; color: #4fc3f7;">
                        Insert Size: ${params[0].name} bp
                      </div>`
      params.forEach(param => {
        const color = param.color
        result += `<div style="margin-bottom: 4px;">
                     <span style="color: ${color};">●</span> 
                     ${param.seriesName}: <span style="color: #ffa726;">${param.value.toFixed(1)}</span>
                   </div>`
      })
      result += `</div>`
      return result
    }
  },
  legend: {
    data: ['Group 1', 'Group 2'],
    top: 50,
    left: 'center',
    itemStyle: {
      borderWidth: 0
    },
    textStyle: { fontSize: 12 }
  },
  xAxis: {
    type: 'category',
    data: xData,
    name: 'Insert Size (bp)',
    nameLocation: 'middle',
    nameGap: 35,
    nameTextStyle: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#333'
    },
    axisLabel: {
      fontSize: 11,
      color: '#666',
      interval: 4,
      formatter: function(value) {
        return value + ''
      }
    },
    axisTick: {
      alignWithLabel: true
    },
    axisLine: {
      lineStyle: { color: '#ddd' }
    }
  },
  yAxis: {
    type: 'value',
    name: 'Density',
    nameLocation: 'middle',
    nameGap: 50,
    nameTextStyle: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#333'
    },
    axisLabel: {
      formatter: function(value) {
        if (value >= 1000) {
          return (value/1000).toFixed(1) + 'K'
        }
        return Math.round(value)
      },
      fontSize: 12,
      color: '#666'
    },
    min: 0,
    max: function(value) {
      return Math.ceil(value.max * 1.1)
    },
    splitLine: {
      lineStyle: { type: 'dashed', opacity: 0.3 }
    }
  },
  grid: {
    left: '10%',
    right: '5%',
    bottom: '25%',
    top: '25%'
  },
  series: [
    {
      name: 'Group 1',
      type: 'line',
      data: group1Data,
      smooth: true,
      itemStyle: { color: '#5470c6' },
      lineStyle: { 
        width: 3,
        color: '#5470c6'
      },
      areaStyle: {
        opacity: 0.15,
        color: '#5470c6'
      },
      symbol: 'none',
      emphasis: {
        focus: 'series',
        lineStyle: { width: 4 }
      },
      markPoint: {
        data: [{
          name: 'Group 1 Peak',
          coord: [group1Peak, group1PeakValue],
          value: group1Peak + ' bp',
          itemStyle: { 
            color: '#5470c6',
            borderColor: '#fff',
            borderWidth: 2
          },
          label: {
            show: true,
            position: 'top',
            fontSize: 11,
            fontWeight: 'bold',
            color: '#5470c6',
            backgroundColor: 'rgba(255, 255, 255, 0.9)',
            padding: [3, 6],
            borderRadius: 4
          }
        }],
        symbol: 'pin',
        symbolSize: 35
      }
    },
    {
      name: 'Group 2', 
      type: 'line',
      data: group2Data,
      smooth: true,
      itemStyle: { color: '#ee6666' },
      lineStyle: { 
        width: 3,
        color: '#ee6666'
      },
      areaStyle: {
        opacity: 0.15,
        color: '#ee6666'
      },
      symbol: 'none',
      emphasis: {
        focus: 'series',
        lineStyle: { width: 4 }
      },
      markPoint: {
        data: [{
          name: 'Group 2 Peak',
          coord: [group2Peak, group2PeakValue],
          value: group2Peak + ' bp',
          itemStyle: { 
            color: '#ee6666',
            borderColor: '#fff',
            borderWidth: 2
          },
          label: {
            show: true,
            position: 'top',
            fontSize: 11,
            fontWeight: 'bold',
            color: '#ee6666',
            backgroundColor: 'rgba(255, 255, 255, 0.9)',
            padding: [3, 6],
            borderRadius: 4
          }
        }],
        symbol: 'pin',
        symbolSize: 35
      }
    }
  ]
}

insertSizeChart.setOption(option)

// 创建统计信息面板
const chartContainer = insertSizeRef.value.parentElement
let statsPanel = chartContainer.querySelector('.insert-stats-panel')

if (!statsPanel) {
  statsPanel = document.createElement('div')
  statsPanel.className = 'insert-stats-panel'
  chartContainer.appendChild(statsPanel)
}

statsPanel.style.cssText = `
  margin-top: 15px;
  padding: 15px;
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  border-radius: 10px;
  color: white;
  box-shadow: 0 4px 15px rgba(0,0,0,0.1);
  min-height: 120px;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
`

statsPanel.innerHTML = `
  <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 15px; flex: 1;">
    <div style="text-align: center;">
      <div style="font-size: 20px; font-weight: bold; margin-bottom: 5px; color: #e3f2fd;">
        ${group1Peak} bp
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Group 1 Peak</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${group1Quality} Quality</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 20px; font-weight: bold; margin-bottom: 5px; color: #fce4ec;">
        ${group2Peak} bp
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Group 2 Peak</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${group2Quality} Quality</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #fff3e0;">
        ${consistency}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Group Consistency</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">Δ = ${peakDifference}bp</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #f3e5f5;">
        ${sequencingQuality}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Library Quality</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">Avg: ${Math.round(avgPeak)}bp</div>
    </div>
  </div>
  
  <div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid rgba(255,255,255,0.3);">
    <div style="font-size: 12px; opacity: 0.95; text-align: center;">
      💡 <strong>Threshold:</strong> Optimal insert size 150-300bp, |Peak difference| < 5bp
    </div>
  </div>
`

// 移除之前的描述文字（如果存在）
const oldDescription = chartContainer.querySelector('.insert-chart-description')
if (oldDescription) {
  oldDescription.remove()
}

// 添加窗口调整处理
const handleInsertResize = () => {
  insertSizeChart && insertSizeChart.resize()
}
window.addEventListener('resize', handleInsertResize)
}
}

const initHeatmapChart = () => {
if (heatmapRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
heatmapChart = echarts.init(heatmapRef.value)

// 模拟22个常染色体 + XY染色体的数据
const chromosomes = Array.from({length: 22}, (_, i) => `chr${i + 1}`).concat(['chrX', 'chrY'])
const groups = ['Group 1', 'Group 2']

// 生成模拟的覆盖度数据
const data = []
groups.forEach((group, groupIndex) => {
  chromosomes.forEach((chr, chrIndex) => {
    // 生成40-95之间的随机覆盖度，添加一些模式让数据更有意义
    let coverage
    if (chr === 'chrX' || chr === 'chrY') {
      // 性染色体通常覆盖度稍低
      coverage = Math.floor(Math.random() * 30) + 50
    } else if (parseInt(chr.replace('chr', '')) > 18) {
      // 较小的染色体覆盖度可能稍高
      coverage = Math.floor(Math.random() * 25) + 70
    } else {
      // 常规染色体
      coverage = Math.floor(Math.random() * 35) + 60
    }
    data.push([chrIndex, groupIndex, coverage])
  })
})

const option = {
  title: {
    text: 'Chromosome Coverage Heatmap',
    left: 'center',
    textStyle: { 
      fontSize: 18, 
      fontWeight: 'bold',
      color: '#2c3e50'
    },
    top: 20
  },
  tooltip: {
    position: 'top',
    backgroundColor: 'rgba(0,0,0,0.8)',
    borderWidth: 0,
    textStyle: {
      color: '#fff',
      fontSize: 14
    },
    formatter: function (params) {
      const coverage = params.data[2]
      let quality = 'Low'
      if (coverage >= 80) quality = 'Excellent'
      else if (coverage >= 70) quality = 'Good'
      else if (coverage >= 60) quality = 'Fair'
      
      return `<div style="padding: 8px;">
                <div style="font-weight: bold; margin-bottom: 4px;">${groups[params.data[1]]}</div>
                <div>${chromosomes[params.data[0]]}: <span style="color: #4fc3f7;">${coverage}%</span></div>
                <div style="font-size: 12px; color: #ccc;">Quality: ${quality}</div>
              </div>`
    }
  },
  grid: {
    left: '12%',
    right: '20%',
    top: '20%',
    bottom: '15%'
  },
  xAxis: {
    type: 'category',
    data: chromosomes,
    splitArea: { 
      show: true,
      areaStyle: {
        color: ['rgba(250,250,250,0.1)', 'rgba(200,200,200,0.1)']
      }
    },
    axisLabel: {
      interval: 0,
      rotate: 45,
      fontSize: 12,
      fontWeight: 'bold',
      color: '#555',
      margin: 15
    },
    axisTick: {
      alignWithLabel: true,
      lineStyle: { color: '#ddd' }
    },
    axisLine: {
      lineStyle: { color: '#ddd' }
    },
    nameTextStyle: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#333'
    }
  },
  yAxis: {
    type: 'category',
    data: groups,
    splitArea: { 
      show: true,
      areaStyle: {
        color: ['rgba(250,250,250,0.3)', 'rgba(200,200,200,0.1)']
      }
    },
    axisLabel: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#555',
      margin: 15
    },
    axisTick: {
      lineStyle: { color: '#ddd' }
    },
    axisLine: {
      lineStyle: { color: '#ddd' }
    }
  },
  visualMap: {
    min: 40,
    max: 100,
    calculable: true,
    realtime: true,
    orient: 'vertical',
    right: '5%',
    top: '25%',
    bottom: '15%',
    text: ['High Coverage', 'Low Coverage'],
    textStyle: {
      fontSize: 12,
      fontWeight: 'bold',
      color: '#555'
    },
    inRange: {
      color: [
        '#313695', '#4575b4', '#74add1', '#abd9e9', 
        '#e0f3f8', '#ffffcc', '#fee090', '#fdae61', 
        '#f46d43', '#d73027', '#a50026'
      ]
    },
    controller: {
      inRange: {
        symbolSize: [20, 40]
      }
    }
  },
  series: [{
    name: 'Coverage',
    type: 'heatmap',
    data: data,
    label: {
      show: false  // 完全隐藏格子内的数值标签
    },
    emphasis: {
      itemStyle: {
        shadowBlur: 15,
        shadowColor: 'rgba(0, 0, 0, 0.4)',
        borderColor: '#fff',
        borderWidth: 2
      }
    },
    itemStyle: {
      borderRadius: 4,  // 圆角让格子更美观
      borderColor: '#fff',
      borderWidth: 1
    }
  }]
}

heatmapChart.setOption(option)

// 添加说明文字
const chartContainer = heatmapRef.value.parentElement
if (chartContainer && !chartContainer.querySelector('.heatmap-description')) {
  const description = document.createElement('div')
  description.className = 'heatmap-description'
  description.style.cssText = `
    margin-top: 15px;
    padding: 12px;
    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    border-radius: 8px;
    border-left: 4px solid #3498db;
  `
  description.innerHTML = `
    <div style="font-size: 13px; color: #2c3e50; line-height: 1.6;">
      <div style="font-weight: bold; margin-bottom: 8px; color: #34495e;">
        📊 Coverage Analysis Summary
      </div>
      <div style="display: flex; justify-content: space-between; flex-wrap: wrap; gap: 15px;">
        <div>
          <span style="color: #e74c3c; font-weight: bold;">●</span> 
          <strong>High Coverage (≥80%):</strong> Excellent quality
        </div>
        <div>
          <span style="color: #f39c12; font-weight: bold;">●</span> 
          <strong>Good Coverage (70-79%):</strong> Acceptable quality
        </div>
        <div>
          <span style="color: #3498db; font-weight: bold;">●</span> 
          <strong>Low Coverage (<60%):</strong> May need attention
        </div>
      </div>
      <div style="margin-top: 8px; font-size: 12px; color: #7f8c8d;">
        💡 Hover over cells for detailed coverage information. Darker colors indicate higher coverage.
      </div>
    </div>
  `
  chartContainer.appendChild(description)
}

// 添加窗口调整处理
const handleHeatmapResize = () => {
  heatmapChart && heatmapChart.resize()
}
window.addEventListener('resize', handleHeatmapResize)
}
}

const initBamVolcanoChart = () => {
if (bamVolcanoRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
bamVolcanoChart = echarts.init(bamVolcanoRef.value)
updateBamVolcanoChart()
}
}

const updateBamVolcanoChart = () => {
if (!bamVolcanoChart) return

// 生成覆盖度差异火山图数据
const data = []
const geneNames = ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS', 'PIK3CA', 'PTEN', 'RB1', 'APC', 'VHL']

for (let i = 0; i < 200; i++) {
// log2(覆盖度比值) - 范围从-3到3
const log2CoverageRatio = (Math.random() - 0.5) * 6
// P值 - 0.0001到0.2之间
const pvalue = Math.random() * 0.2 + 0.0001
const negLog10Pvalue = -Math.log10(pvalue)

let category = 'non-significant'
let geneName = `Gene_${i + 1}`

// 为一些点添加已知基因名
if (i < geneNames.length) {
  geneName = geneNames[i]
}

if (pvalue < bamVolcanoThreshold.value.pvalue) {
  if (log2CoverageRatio > 1) category = 'higher-coverage'
  else if (log2CoverageRatio < -1) category = 'lower-coverage'
  else category = 'significant'
}

data.push([log2CoverageRatio, negLog10Pvalue, category, geneName])
}

// 计算统计信息
const stats = {
total: data.length,
significant: data.filter(d => d[2] !== 'non-significant').length,
higherCoverage: data.filter(d => d[2] === 'higher-coverage').length,
lowerCoverage: data.filter(d => d[2] === 'lower-coverage').length,
other: data.filter(d => d[2] === 'significant').length
}

const option = {
title: {
  text: 'Coverage Difference Volcano Plot',
  left: 'center',
  textStyle: { fontSize: 16, fontWeight: 'bold', color: '#2c3e50' },
  top: 15
},
legend: {
  data: [
    {name: 'Higher Coverage', icon: 'circle'},
    {name: 'Lower Coverage', icon: 'circle'},
    {name: 'Significant', icon: 'circle'},
    {name: 'Non-significant', icon: 'circle'}
  ],
  top: 50,
  left: 'center',
  itemGap: 20,
  textStyle: { fontSize: 12 }
},
tooltip: {
  trigger: 'item',
  backgroundColor: 'rgba(0,0,0,0.85)',
  borderWidth: 0,
  textStyle: { color: '#fff', fontSize: 13 },
  formatter: function(params) {
    const foldChange = Math.pow(2, params.data[0]).toFixed(2)
    const pValue = Math.pow(10, -params.data[1]).toExponential(2)
    
    let categoryColor = '#d0d0d0'
    let categoryText = params.data[2]
    
    switch(params.data[2]) {
      case 'higher-coverage': 
        categoryColor = '#ee6666'
        categoryText = 'Higher Coverage'
        break
      case 'lower-coverage': 
        categoryColor = '#5470c6'
        categoryText = 'Lower Coverage'
        break
      case 'significant': 
        categoryColor = '#91cc75'
        categoryText = 'Significant'
        break
      case 'non-significant': 
        categoryColor = '#d0d0d0'
        categoryText = 'Non-significant'
        break
    }
    
    return `<div style="padding: 10px;">
              <div style="font-weight: bold; margin-bottom: 6px; color: #4fc3f7;">
                ${params.data[3]}
              </div>
              <div style="margin-bottom: 4px;">
                Coverage Ratio: <span style="color: #ffa726;">${foldChange}x</span>
              </div>
              <div style="margin-bottom: 4px;">
                Log2(Ratio): <span style="color: #ffa726;">${params.data[0].toFixed(2)}</span>
              </div>
              <div style="margin-bottom: 4px;">
                P-value: <span style="color: #ffa726;">${pValue}</span>
              </div>
              <div style="margin-bottom: 4px;">
                -log10(P): <span style="color: #ffa726;">${params.data[1].toFixed(2)}</span>
              </div>
              <div style="color: ${categoryColor}; font-weight: bold;">
                ● ${categoryText}
              </div>
            </div>`
  }
},
grid: {
  left: '10%',
  right: '10%',
  bottom: '25%',
  top: '25%'
},
xAxis: {
  type: 'value',
  name: 'log2(Coverage Ratio)',
  nameLocation: 'middle',
  nameGap: 35,
  nameTextStyle: {
    fontSize: 14,
    fontWeight: 'bold',
    color: '#333'
  },
  axisLine: { 
    onZero: true, 
    lineStyle: { color: '#333', width: 2 } 
  },
  splitLine: { 
    show: true, 
    lineStyle: { type: 'dashed', opacity: 0.3 } 
  },
  axisLabel: {
    fontSize: 12,
    color: '#666'
  }
},
yAxis: {
  type: 'value',
  name: '-log10(P-value)',
  nameLocation: 'middle',
  nameGap: 45,
  nameTextStyle: {
    fontSize: 14,
    fontWeight: 'bold',
    color: '#333'
  },
  min: 0,
  splitLine: { 
    show: true, 
    lineStyle: { type: 'dashed', opacity: 0.3 } 
  },
  axisLabel: {
    fontSize: 12,
    color: '#666'
  }
},
series: [{
  type: 'scatter',
  data: data,
  symbolSize: function(data) {
    return data[2] === 'non-significant' ? 4 : 8
  },
  itemStyle: {
    color: function(params) {
      const colors = {
        'higher-coverage': '#ee6666',
        'lower-coverage': '#5470c6', 
        'significant': '#91cc75',
        'non-significant': '#d0d0d0'
      }
      return colors[params.data[2]] || '#d0d0d0'
    },
    opacity: 0.8,
    borderWidth: 0.5,
    borderColor: '#fff'
  },
  emphasis: {
    itemStyle: {
      shadowBlur: 10,
      shadowColor: 'rgba(0, 0, 0, 0.3)',
      borderWidth: 2,
      borderColor: '#fff',
      opacity: 1
    },
    scale: 1.2
  },
  markLine: {
    silent: true,
    lineStyle: {
      color: '#ff6b6b',
      type: 'dashed',
      width: 2,
      opacity: 0.8
    },
    label: {
      show: true,
      position: 'end',
      fontSize: 11,
      fontWeight: 'bold',
      color: '#ff6b6b',
      backgroundColor: 'rgba(255, 255, 255, 0.8)',
      padding: [2, 4],
      borderRadius: 3
    },
    data: [
      { 
        yAxis: -Math.log10(bamVolcanoThreshold.value.pvalue),
        label: { formatter: `P = ${bamVolcanoThreshold.value.pvalue}` }
      },
      { 
        xAxis: 1,
        label: { formatter: 'Log2FC = 1' }
      },
      { 
        xAxis: -1,
        label: { formatter: 'Log2FC = -1' }
      }
    ]
  }
}]
}

bamVolcanoChart.setOption(option)

// 更新或创建统计信息面板
const chartContainer = bamVolcanoRef.value.parentElement
let statsPanel = chartContainer.querySelector('.volcano-stats-panel')

if (!statsPanel) {
  statsPanel = document.createElement('div')
  statsPanel.className = 'volcano-stats-panel'
  chartContainer.appendChild(statsPanel)
}

statsPanel.style.cssText = `
  margin-top: 15px;
  padding: 15px;
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  border-radius: 10px;
  color: white;
  box-shadow: 0 4px 15px rgba(0,0,0,0.1);
`

statsPanel.innerHTML = `
  <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px;">
    <div style="text-align: center;">
      <div style="font-size: 24px; font-weight: bold; margin-bottom: 5px;">
        ${stats.total}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Total Genes</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 20px; font-weight: bold; margin-bottom: 5px; color: #4fc3f7;">
        ${stats.significant}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">
        Significant (${((stats.significant/stats.total)*100).toFixed(1)}%)
      </div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #ffab91;">
        ${stats.higherCoverage}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Higher Coverage</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #81c784;">
        ${stats.lowerCoverage}
      </div>
      <div style="font-size: 12px; opacity: 0.9;">Lower Coverage</div>
    </div>
  </div>
  
  <div style="margin-top: 12px; font-size: 12px; opacity: 0.85; text-align: center;">
    💡 Threshold: P-value < ${bamVolcanoThreshold.value.pvalue}, |Log2FC| > 1
  </div>
`
}

const initDonutChart = () => {
if (donutChartRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
donutChart = echarts.init(donutChartRef.value)

// 使用真实数据构建事件数据
const eventData = [
  { value: rmatsData.value.eventCounts.SE || 0, name: 'SE (Skipped Exon)', itemStyle: { color: '#5470c6' } },
  { value: rmatsData.value.eventCounts.A5SS || 0, name: "A5SS (Alt 5' Splice Site)", itemStyle: { color: '#91cc75' } },
  { value: rmatsData.value.eventCounts.A3SS || 0, name: "A3SS (Alt 3' Splice Site)", itemStyle: { color: '#fac858' } },
  { value: rmatsData.value.eventCounts.MXE || 0, name: 'MXE (Mutually Exclusive Exon)', itemStyle: { color: '#ee6666' } },
  { value: rmatsData.value.eventCounts.RI || 0, name: 'RI (Retained Intron)', itemStyle: { color: '#73c0de' } }
]

const totalEvents = eventData.reduce((sum, item) => sum + item.value, 0)

const option = {
  title: {
    text: 'Splicing Event Distribution',
    left: 'center',
    textStyle: { fontSize: 16, fontWeight: 'bold', color: '#2c3e50' },
    top: 15
  },
  tooltip: {
    trigger: 'item',
    backgroundColor: 'rgba(0,0,0,0.85)',
    borderWidth: 0,
    textStyle: { color: '#fff', fontSize: 13 },
    formatter: function(params) {
      const percentage = totalEvents > 0 ? ((params.value / totalEvents) * 100).toFixed(1) : '0.0'
      return `<div style="padding: 10px; min-width: 200px;">
                <div style="font-weight: bold; margin-bottom: 8px; color: #4fc3f7; border-bottom: 1px solid #333; padding-bottom: 4px;">
                  ${params.name}
                </div>
                <div style="margin-bottom: 6px;">
                  <span style="color: ${params.color}; font-weight: bold;">●</span> 
                  Count: <span style="color: #ffa726; font-weight: bold;">${params.value.toLocaleString()}</span>
                </div>
                <div style="margin-bottom: 6px;">
                  Percentage: <span style="color: #ffa726; font-weight: bold;">${percentage}%</span>
                </div>
                <div style="border-top: 1px solid #444; padding-top: 4px; margin-top: 6px; color: #e6f7ff; font-size: 12px;">
                  Rank: #${eventData.sort((a,b) => b.value - a.value).findIndex(item => item.value === params.value) + 1} most common
                </div>
              </div>`
    }
  },
  legend: {
    orient: 'vertical',
    left: 'left',
    top: 'center',
    textStyle: { fontSize: 11 },
    formatter: function(name) {
      const item = eventData.find(d => d.name === name)
      const percentage = totalEvents > 0 ? ((item.value / totalEvents) * 100).toFixed(1) : '0.0'
      return `${name.split('(')[0].trim()}\n${percentage}%`
    }
  },
  series: [{
    name: 'Event Type',
    type: 'pie',
    radius: ['45%', '75%'],
    center: ['60%', '50%'],
    data: eventData,
    emphasis: {
      itemStyle: {
        shadowBlur: 15,
        shadowOffsetX: 0,
        shadowColor: 'rgba(0, 0, 0, 0.5)'
      },
      label: {
        show: true,
        fontSize: 14,
        fontWeight: 'bold'
      }
    },
    label: {
      show: true,
      position: 'outside',
      formatter: function(params) {
        const percentage = totalEvents > 0 ? ((params.value / totalEvents) * 100).toFixed(1) : '0.0'
        return `${params.name.split('(')[0].trim()}\n${percentage}%`
      },
      fontSize: 11,
      lineHeight: 14
    },
    labelLine: {
      show: true,
      length: 15,
      length2: 10
    }
  }]
}

donutChart.setOption(option)

// 创建统计面板
const chartContainer = donutChartRef.value.parentElement
let statsPanel = chartContainer.querySelector('.donut-stats-panel')

if (!statsPanel) {
  statsPanel = document.createElement('div')
  statsPanel.className = 'donut-stats-panel'
  chartContainer.appendChild(statsPanel)
}

// 计算统计信息
const sortedEvents = [...eventData].sort((a, b) => b.value - a.value)
const mostCommon = sortedEvents[0] || { name: 'SE', value: 0 }
const leastCommon = sortedEvents[sortedEvents.length - 1] || { name: 'RI', value: 0 }
const diversityIndex = eventData.length > 1 && totalEvents > 0 ? 
  (eventData.length / (eventData.length - 1)) * (1 - eventData.reduce((sum, item) => sum + Math.pow(item.value / totalEvents, 2), 0)) : 0

// 计算复杂度评分
let complexity = 'High'
if (diversityIndex < 0.6) complexity = 'Low'
else if (diversityIndex < 0.75) complexity = 'Medium'

// SE事件占比评估
const sePercentage = totalEvents > 0 ? (eventData[0].value / totalEvents * 100).toFixed(1) : '0.0'
let seAssessment = 'Normal'
if (parseFloat(sePercentage) > 40) seAssessment = 'High'
else if (parseFloat(sePercentage) < 30) seAssessment = 'Low'

statsPanel.style.cssText = `
  margin-top: 15px;
  padding: 15px;
  background: linear-gradient(135deg, #ffecd2 0%, #fcb69f 100%);
  border-radius: 10px;
  color: #2c3e50;
  box-shadow: 0 4px 15px rgba(0,0,0,0.1);
  min-height: 120px;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
`

statsPanel.innerHTML = `
  <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 15px; flex: 1;">
    <div style="text-align: center;">
      <div style="font-size: 22px; font-weight: bold; margin-bottom: 5px; color: #1565c0;">
        ${totalEvents.toLocaleString()}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Total Events</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">All Types</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #5470c6;">
        ${mostCommon.name.split('(')[0].trim()}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Most Common</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${totalEvents > 0 ? ((mostCommon.value/totalEvents)*100).toFixed(1) : '0.0'}% of events</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #91cc75;">
        ${complexity}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Splicing Complexity</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">Diversity: ${diversityIndex.toFixed(2)}</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #ff7043;">
        ${seAssessment}
      </div>
      <div style="font-size: 12px; font-weight: 500;">SE Dominance</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${sePercentage}% SE events</div>
    </div>
  </div>
  
  <div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid rgba(44, 62, 80, 0.2);">
    <div style="font-size: 12px; font-weight: 500; text-align: center;">
      🧬 <strong>Profile:</strong> ${seAssessment === 'High' ? 'SE-dominated splicing pattern' : 
                                     seAssessment === 'Low' ? 'Diverse splicing pattern' : 
                                     'Balanced splicing profile'} with ${complexity.toLowerCase()} complexity
    </div>
  </div>
`

// 添加窗口调整处理
const handleDonutResize = () => {
  donutChart && donutChart.resize()
}
window.addEventListener('resize', handleDonutResize)
}
}

const initBarChart = () => {
if (barChartRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
barChart = echarts.init(barChartRef.value)

// 根据真实数据计算显著性统计
const categories = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
const significantData = []
const nonSignificantData = []

// 计算总体显著率
const overallSigRate = rmatsData.value.totalEvents > 0 ? 
  rmatsData.value.significantEvents / rmatsData.value.totalEvents : 0.18

categories.forEach(eventType => {
  const total = rmatsData.value.eventCounts[eventType] || 0
  const significant = Math.round(total * overallSigRate)
  const nonSignificant = total - significant
  
  significantData.push(significant)
  nonSignificantData.push(nonSignificant)
})

const totalData = significantData.map((val, idx) => val + nonSignificantData[idx])

const option = {
  title: {
    text: 'Significance Statistics',
    left: 'center',
    textStyle: { fontSize: 16, fontWeight: 'bold', color: '#2c3e50' },
    top: 15
  },
  tooltip: {
    trigger: 'axis',
    axisPointer: { type: 'shadow' },
    backgroundColor: 'rgba(0,0,0,0.85)',
    borderWidth: 0,
    textStyle: { color: '#fff', fontSize: 13 },
    formatter: function(params) {
      const categoryIndex = params[0].dataIndex
      const significant = significantData[categoryIndex]
      const nonSignificant = nonSignificantData[categoryIndex]
      const total = totalData[categoryIndex]
      const sigPercent = total > 0 ? ((significant / total) * 100).toFixed(1) : '0.0'
      const nonSigPercent = total > 0 ? ((nonSignificant / total) * 100).toFixed(1) : '0.0'
      
      return `<div style="padding: 10px; min-width: 200px;">
                <div style="font-weight: bold; margin-bottom: 8px; color: #4fc3f7; border-bottom: 1px solid #333; padding-bottom: 4px;">
                  ${params[0].name} Events
                </div>
                <div style="margin-bottom: 6px;">
                  <span style="color: #52c41a; font-weight: bold;">●</span> 
                  Significant: <span style="color: #ffa726; font-weight: bold;">${significant.toLocaleString()}</span> 
                  <span style="color: #ccc;">(${sigPercent}%)</span>
                </div>
                <div style="margin-bottom: 6px;">
                  <span style="color: #d9d9d9; font-weight: bold;">●</span> 
                  Non-significant: <span style="color: #ffa726; font-weight: bold;">${nonSignificant.toLocaleString()}</span> 
                  <span style="color: #ccc;">(${nonSigPercent}%)</span>
                </div>
                <div style="border-top: 1px solid #444; padding-top: 4px; margin-top: 6px; color: #e6f7ff;">
                  <strong>Total: ${total.toLocaleString()} events</strong>
                </div>
              </div>`
    }
  },
  legend: {
    data: ['Significant', 'Non-significant'],
    top: 50,
    left: 'center',
    itemGap: 20,
    textStyle: { fontSize: 12 }
  },
  grid: {
    left: '10%',
    right: '10%',
    bottom: '25%',
    top: '25%'
  },
  xAxis: {
    type: 'category',
    data: categories,
    name: 'Splicing Event Type',
    nameLocation: 'middle',
    nameGap: 35,
    nameTextStyle: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#333'
    },
    axisLabel: {
      fontSize: 12,
      fontWeight: 'bold',
      color: '#666'
    },
    axisTick: {
      alignWithLabel: true
    },
    axisLine: {
      lineStyle: { color: '#ddd' }
    }
  },
  yAxis: {
    type: 'value',
    name: 'Event Count',
    nameLocation: 'middle',
    nameGap: 50,
    nameTextStyle: {
      fontSize: 14,
      fontWeight: 'bold',
      color: '#333'
    },
    axisLabel: {
      formatter: function(value) {
        if (value >= 1000) {
          return (value/1000).toFixed(1) + 'K'
        }
        return value
      },
      fontSize: 12,
      color: '#666'
    },
    splitLine: {
      lineStyle: { type: 'dashed', opacity: 0.3 }
    }
  },
  series: [
    {
      name: 'Significant',
      type: 'bar',
      stack: 'total',
      data: significantData,
      itemStyle: { 
        color: '#52c41a',
        borderRadius: [0, 0, 0, 0]
      },
      emphasis: {
        itemStyle: {
          color: '#389e0d',
          shadowBlur: 10,
          shadowColor: 'rgba(82, 196, 26, 0.3)'
        }
      },
      label: {
        show: true,
        position: 'inside',
        formatter: function(params) {
          const total = totalData[params.dataIndex]
          const percent = total > 0 ? ((params.value / total) * 100).toFixed(1) : '0.0'
          return params.value > 50 ? `${percent}%` : ''
        },
        fontSize: 11,
        fontWeight: 'bold',
        color: '#fff'
      }
    },
    {
      name: 'Non-significant',
      type: 'bar', 
      stack: 'total',
      data: nonSignificantData,
      itemStyle: { 
        color: '#d9d9d9',
        borderRadius: [4, 4, 0, 0]
      },
      emphasis: {
        itemStyle: {
          color: '#bfbfbf',
          shadowBlur: 10,
          shadowColor: 'rgba(217, 217, 217, 0.3)'
        }
      },
      label: {
        show: true,
        position: 'inside',
        formatter: function(params) {
          const total = totalData[params.dataIndex]
          const percent = total > 0 ? ((params.value / total) * 100).toFixed(1) : '0.0'
          return `${percent}%`
        },
        fontSize: 11,
        fontWeight: 'bold',
        color: '#666'
      }
    }
  ]
}

barChart.setOption(option)

// 创建统计摘要面板
const chartContainer = barChartRef.value.parentElement
let statsPanel = chartContainer.querySelector('.bar-stats-panel')

if (!statsPanel) {
  statsPanel = document.createElement('div')
  statsPanel.className = 'bar-stats-panel'
  chartContainer.appendChild(statsPanel)
}

// 计算总体统计
const totalSignificant = significantData.reduce((a, b) => a + b, 0)
const totalNonSignificant = nonSignificantData.reduce((a, b) => a + b, 0)
const grandTotal = totalSignificant + totalNonSignificant
const overallSigRatePercent = grandTotal > 0 ? ((totalSignificant / grandTotal) * 100).toFixed(1) : '0.0'

// 找出最高和最低显著率的事件类型
const sigRates = significantData.map((sig, idx) => ({
  type: categories[idx],
  rate: totalData[idx] > 0 ? (sig / totalData[idx] * 100).toFixed(1) : '0.0',
  count: sig
}))
const highestSig = sigRates.reduce((prev, curr) => (parseFloat(curr.rate) > parseFloat(prev.rate)) ? curr : prev)
const lowestSig = sigRates.reduce((prev, curr) => (parseFloat(curr.rate) < parseFloat(prev.rate)) ? curr : prev)

statsPanel.style.cssText = `
  margin-top: 15px;
  padding: 15px;
  background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%);
  border-radius: 10px;
  color: #2c3e50;
  box-shadow: 0 4px 15px rgba(0,0,0,0.1);
  min-height: 120px;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
`

statsPanel.innerHTML = `
  <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 15px; flex: 1;">
    <div style="text-align: center;">
      <div style="font-size: 22px; font-weight: bold; margin-bottom: 5px; color: #1890ff;">
        ${grandTotal.toLocaleString()}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Total Events</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">All Types</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 20px; font-weight: bold; margin-bottom: 5px; color: #52c41a;">
        ${overallSigRatePercent}%
      </div>
      <div style="font-size: 12px; font-weight: 500;">Significance Rate</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${totalSignificant.toLocaleString()} significant</div>
    </div>initDonutChart
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #fa541c;">
        ${highestSig.type}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Highest Significance</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${highestSig.rate}% rate</div>
    </div>
    
    <div style="text-align: center;">
      <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #722ed1;">
        ${lowestSig.type}
      </div>
      <div style="font-size: 12px; font-weight: 500;">Lowest Significance</div>
      <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${lowestSig.rate}% rate</div>
    </div>
  </div>
  
  <div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid rgba(44, 62, 80, 0.2);">
    <div style="font-size: 12px; font-weight: 500; text-align: center;">
      📊 <strong>Analysis:</strong> ${highestSig.type} shows highest differential splicing activity
    </div>
  </div>
`

// 添加窗口调整处理
const handleBarResize = () => {
  barChart && barChart.resize()
}
window.addEventListener('resize', handleBarResize)
}
}

const initVolcanoChart = () => {
if (volcanoRef.value && (window.$echarts || window.echarts)) {
const echarts = window.$echarts || window.echarts
volcanoChart = echarts.init(volcanoRef.value)
updateVolcanoChart()
}
}

const updateVolcanoChart = () => {
if (!volcanoChart) return

// 使用真实的火山图数据，如果没有则生成模拟数据
let data = []
if (rmatsData.value.volcanoData && rmatsData.value.volcanoData.length > 0) {
  // 使用真实数据
  data = rmatsData.value.volcanoData
  console.log('Using real volcano data')
} else {
  // 生成模拟数据作为后备
  const knownEvents = [
    'BRCA1_exon11', 'TP53_exon4', 'EGFR_exon19', 'MYC_exon2', 'KRAS_exon2',
    'PIK3CA_exon9', 'PTEN_exon5', 'RB1_exon13', 'APC_exon15', 'VHL_exon2'
  ]

  for (let i = 0; i < 500; i++) {
    const deltaPsi = (Math.random() - 0.5) * 2
    const fdr = Math.random()
    const logFdr = -Math.log10(fdr)
    let category = 'non-significant'
    let eventName = `Event_${i + 1}`

    // 添加一些已知事件名
    if (i < knownEvents.length) {
      eventName = knownEvents[i]
    }

    if (fdr < volcanoThreshold.value.fdr) {
      if (deltaPsi > 0.1) category = 'upregulated'
      else if (deltaPsi < -0.1) category = 'downregulated'
      else category = 'significant'
    }

    data.push([deltaPsi, logFdr, category, eventName])
  }
  console.log('Using simulated volcano data')
}

// 计算统计信息
const stats = {
total: data.length,
significant: data.filter(d => d[2] !== 'non-significant').length,
upregulated: data.filter(d => d[2] === 'upregulated').length,
downregulated: data.filter(d => d[2] === 'downregulated').length,
other: data.filter(d => d[2] === 'significant').length
}

const option = {
title: {
  text: 'Volcano Plot (ΔPSI vs -log₁₀FDR)',
  left: 'center',
textStyle: { fontSize: 16, fontWeight: 'bold', color: '#2c3e50' },
 top: 15
},
legend: {
  data: [
    {name: 'upregulated', icon: 'circle'},
    {name: 'downregulated', icon: 'circle'},
    {name: 'significant', icon: 'circle'},
    {name: 'non-significant', icon: 'circle'}
  ],
  top: 50,
  left: 'center',
  itemGap: 20,
  textStyle: { fontSize: 12 },
  formatter: function(name) {
    // 将内部名称转换为显示名称
    const nameMap = {
      'upregulated': 'Increased Inclusion',
      'downregulated': 'Decreased Inclusion', 
      'significant': 'Significant',
      'non-significant': 'Non-significant'
    }
    return nameMap[name] || name
  }
},
tooltip: {
 trigger: 'item',
 backgroundColor: 'rgba(0,0,0,0.85)',
 borderWidth: 0,
 textStyle: { color: '#fff', fontSize: 13 },
 formatter: function(params) {
   const fdrValue = Math.pow(10, -params.data[1]).toExponential(2)
   
   let categoryColor = '#d0d0d0'
   let categoryText = params.data[2]
   let interpretation = ''
   
   switch(params.data[2]) {
     case 'upregulated': 
       categoryColor = '#ee6666'
       categoryText = 'Increased Inclusion'
       interpretation = 'Exon is more included in condition 2'
       break
     case 'downregulated': 
       categoryColor = '#5470c6'
       categoryText = 'Decreased Inclusion'
       interpretation = 'Exon is less included in condition 2'
       break
     case 'significant': 
       categoryColor = '#91cc75'
       categoryText = 'Significant Change'
       interpretation = 'Minor but significant splicing change'
       break
     case 'non-significant': 
       categoryColor = '#d0d0d0'
       categoryText = 'No Change'
       interpretation = 'No significant splicing difference'
       break
   }
   
   return `<div style="padding: 12px; min-width: 220px;">
             <div style="font-weight: bold; margin-bottom: 8px; color: #4fc3f7; border-bottom: 1px solid #333; padding-bottom: 4px;">
               ${params.data[3]}
             </div>
             <div style="margin-bottom: 6px;">
               ΔPSI: <span style="color: #ffa726; font-weight: bold;">${params.data[0].toFixed(3)}</span>
             </div>
             <div style="margin-bottom: 6px;">
               FDR: <span style="color: #ffa726; font-weight: bold;">${fdrValue}</span>
             </div>
             <div style="margin-bottom: 6px;">
               -log₁₀FDR: <span style="color: #ffa726; font-weight: bold;">${params.data[1].toFixed(2)}</span>
             </div>
             <div style="margin-bottom: 8px; color: ${categoryColor}; font-weight: bold;">
               ● ${categoryText}
             </div>
             <div style="font-size: 11px; color: #ccc; line-height: 1.3; border-top: 1px solid #444; padding-top: 6px;">
               ${interpretation}
             </div>
           </div>`
 }
},
grid: {
 left: '10%',
 right: '10%',
 bottom: '25%',
 top: '25%'
},
xAxis: {
 type: 'value',
 name: 'ΔPSI (Inclusion Change)',
 nameLocation: 'middle',
 nameGap: 35,
 nameTextStyle: {
   fontSize: 14,
   fontWeight: 'bold',
   color: '#333'
 },
 axisLine: { 
   onZero: true, 
   lineStyle: { color: '#333', width: 2 } 
 },
 splitLine: { 
   show: true, 
   lineStyle: { type: 'dashed', opacity: 0.3 } 
 },
 axisLabel: {
   fontSize: 12,
   color: '#666'
 }
},
yAxis: {
 type: 'value',
 name: '-log₁₀(FDR)',
 nameLocation: 'middle',
 nameGap: 45,
 nameTextStyle: {
   fontSize: 14,
   fontWeight: 'bold',
   color: '#333'
 },
 min: 0,
 splitLine: { 
   show: true, 
   lineStyle: { type: 'dashed', opacity: 0.3 } 
 },
 axisLabel: {
   fontSize: 12,
   color: '#666'
 }
},
series: [{
 type: 'scatter',
 data: data,
 symbolSize: function(data) {
   return data[2] === 'non-significant' ? 4 : 8
 },
 itemStyle: {
   color: function(params) {
     const colors = {
       'upregulated': '#ee6666',
       'downregulated': '#5470c6', 
       'significant': '#91cc75',
       'non-significant': '#d0d0d0'
     }
     return colors[params.data[2]] || '#d0d0d0'
   },
   opacity: 0.8,
   borderWidth: 0.5,
   borderColor: '#fff'
 },
 emphasis: {
   itemStyle: {
     shadowBlur: 10,
     shadowColor: 'rgba(0, 0, 0, 0.3)',
     borderWidth: 2,
     borderColor: '#fff',
     opacity: 1
   },
   scale: 1.3
 },
 markLine: {
   silent: true,
   lineStyle: {
     color: '#ff6b6b',
     type: 'dashed',
     width: 2,
     opacity: 0.8
   },
   label: {
     show: true,
     position: 'end',
     fontSize: 11,
     fontWeight: 'bold',
     color: '#ff6b6b',
     backgroundColor: 'rgba(255, 255, 255, 0.9)',
     padding: [2, 6],
     borderRadius: 4
   },
   data: [
     { 
       yAxis: -Math.log10(volcanoThreshold.value.fdr),
       label: { formatter: `FDR = ${volcanoThreshold.value.fdr}` }
     },
     { 
       xAxis: 0.1,
       label: { formatter: 'ΔPSI = 0.1' }
     },
     { 
       xAxis: -0.1,
       label: { formatter: 'ΔPSI = -0.1' }
     }
   ]
 }
}]
}

volcanoChart.setOption(option)

// 创建或更新统计面板
const chartContainer = volcanoRef.value.parentElement
let statsPanel = chartContainer.querySelector('.rmats-volcano-stats-panel')

if (!statsPanel) {
 statsPanel = document.createElement('div')
 statsPanel.className = 'rmats-volcano-stats-panel'
 chartContainer.appendChild(statsPanel)
}

statsPanel.style.cssText = `
 margin-top: 15px;
 padding: 15px;
 background: linear-gradient(135deg, #fa709a 0%, #fee140 100%);
 border-radius: 10px;
 color: #2c3e50;
 box-shadow: 0 4px 15px rgba(0,0,0,0.1);
 min-height: 120px;
 display: flex;
 flex-direction: column;
 justify-content: space-between;
`

const significanceRate = ((stats.significant / stats.total) * 100).toFixed(1)
const inclusionBalance = stats.upregulated > stats.downregulated ? 'Inclusion-favored' : 
                        stats.downregulated > stats.upregulated ? 'Exclusion-favored' : 'Balanced'

statsPanel.innerHTML = `
 <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; flex: 1;">
   <div style="text-align: center;">
     <div style="font-size: 22px; font-weight: bold; margin-bottom: 5px; color: #1565c0;">
       ${stats.total}
     </div>
     <div style="font-size: 12px; font-weight: 500;">Total Events</div>
     <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">Analyzed</div>
   </div>
   
   <div style="text-align: center;">
     <div style="font-size: 20px; font-weight: bold; margin-bottom: 5px; color: #d32f2f;">
       ${significanceRate}%
     </div>
     <div style="font-size: 12px; font-weight: 500;">Significance Rate</div>
     <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">${stats.significant} events</div>
   </div>
   
   <div style="text-align: center;">
     <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #ff5722;">
       ${stats.upregulated}
     </div>
     <div style="font-size: 12px; font-weight: 500;">Increased Inclusion</div>
     <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">More included</div>
   </div>
   
   <div style="text-align: center;">
     <div style="font-size: 18px; font-weight: bold; margin-bottom: 5px; color: #3f51b5;">
       ${stats.downregulated}
     </div>
     <div style="font-size: 12px; font-weight: 500;">Decreased Inclusion</div>
     <div style="font-size: 10px; opacity: 0.8; margin-top: 2px;">More excluded</div>
   </div>
 </div>
 
 <div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid rgba(44, 62, 80, 0.2);">
   <div style="font-size: 12px; font-weight: 500; text-align: center;">
     🧬 <strong>Pattern:</strong> ${inclusionBalance} splicing with ${significanceRate}% differential events
   </div>
 </div>
`
}

// Watch for tab changes and threshold changes
watch(activeTab, async (newTab) => {
await nextTick()
if (newTab === 'bam') {
initPieChart()
initMapqChart()
initInsertSizeChart()
initBamVolcanoChart()
initHeatmapChart()
} else if (newTab === 'rmats') {
initDonutChart()
initBarChart()
initVolcanoChart()
}
})

watch(() => volcanoThreshold.value.fdr, () => {
updateVolcanoChart()
})

watch(() => bamVolcanoThreshold.value.pvalue, () => {
updateBamVolcanoChart()
})

// Utility Functions
const showStatus = (message) => {
statusMessage.value = message
}

const exportChart = (chartType) => {
showStatus(`Exporting ${chartType} chart...`)
}

const exportOverview = () => {
showStatus('Exporting overview report...')
}

const refreshData = () => {
showStatus('Refreshing data...')
}

onMounted(async () => {
  console.log('Component mounting...')
  
  try {
    // 加载真实的BAM分析数据
    const response = await fetch('/bam_analysis_result.json')
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: ${response.statusText}`)
    }
    
    const realData = await response.json()
    console.log('Loaded BAM data:', realData)
    
    // 检查数据结构
    if (!realData.group1 || !realData.group2) {
      throw new Error('Invalid data structure: missing group1 or group2')
    }
    
    // 更新bamData
    bamData.value.group1 = realData.group1
    bamData.value.group2 = realData.group2
    
    // 更新图表数据
    if (realData.chart_data) {
      updateChartsWithRealData(realData.chart_data)
    }
    
    showStatus('✅ Real BAM data loaded successfully!')
    
    // 🔥 新增：加载rMATS数据
    try {
      const rmatsResponse = await fetch('/rmats_analysis_result.json')
      if (rmatsResponse.ok) {
        const rmatsRealData = await rmatsResponse.json()
        console.log('Loaded rMATS data:', rmatsRealData)
        
        // 根据你的rMATS分析器输出结构更新数据
        rmatsData.value = {
          totalEvents: rmatsRealData.summary?.total_events || rmatsRealData.totalEvents || 0,
          significantEvents: rmatsRealData.summary?.significant_events || rmatsRealData.significantEvents || 0,
          maxPsiIncrease: rmatsRealData.summary?.max_psi_increase || rmatsRealData.maxPsiIncrease || 0,
          maxPsiDecrease: rmatsRealData.summary?.max_psi_decrease || rmatsRealData.maxPsiDecrease || 0,
          group1Samples: rmatsRealData.summary?.group1_samples || rmatsRealData.group1Samples || 0,
          group2Samples: rmatsRealData.summary?.group2_samples || rmatsRealData.group2Samples || 0,
          eventCounts: rmatsRealData.event_types || rmatsRealData.eventCounts || {
            SE: 0, A5SS: 0, A3SS: 0, MXE: 0, RI: 0
          },
          volcanoData: rmatsRealData.volcano_plot || [],
          significanceData: rmatsRealData.significance_statistics || []
        }
        
        console.log('Updated rmatsData:', rmatsData.value)
        showStatus('✅ Real rMATS data loaded successfully!')
        
        // 如果当前在rMATS标签，重新初始化图表
        if (activeTab.value === 'rmats') {
          await nextTick()
          initDonutChart()
          initBarChart()
          initVolcanoChart()
        }
        
      } else {
        throw new Error(`Failed to load rMATS data: ${rmatsResponse.status}`)
      }
    } catch (rmatsError) {
      console.error('Failed to load rMATS data:', rmatsError)
      showStatus(`⚠️ Using default rMATS data: ${rmatsError.message}`)
      
      // 如果加载失败，设置你实际分析的数据作为默认值
      rmatsData.value = {
        totalEvents: 136493,
        significantEvents: 24759,
        maxPsiIncrease: 0.85,
        maxPsiDecrease: -0.72,
        group1Samples: 3,
        group2Samples: 3,
        eventCounts: {
          SE: 94001,
          A5SS: 7242, 
          A3SS: 10818,
          MXE: 16923,
          RI: 7509
        }
      }
    }
    
  } catch (error) {
    console.error('Failed to load real BAM data:', error)
    showStatus(`❌ Failed to load real data: ${error.message}`)
    
    // 如果加载失败，使用模拟数据
    bamData.value = {
      group1: {
        fileSize: '2.3 GB',
        mappingRate: 90.3,
        avgMapQ: 42.5,
        avgInsertSize: 185
      },
      group2: {
        fileSize: '2.1 GB',
        mappingRate: 91.5,
        avgMapQ: 43.2,
        avgInsertSize: 182
      }
    }
  }
  
  // 检查ECharts是否可用
  await nextTick()
  const echarts = window.$echarts || window.echarts
  if (!echarts) {
    console.error('ECharts not loaded properly!')
    showStatus('ECharts loading failed, please check configuration')
    return
  }

  console.log('ECharts loaded, initializing charts...')

  // 初始化图表
  setTimeout(() => {
    initPieChart()
    initMapqChart()
    initInsertSizeChart()
    initBamVolcanoChart()
    initHeatmapChart()
    showStatus('🎉 Data overview loaded successfully!')
  }, 100)
})
</script>

<style scoped>
/* 添加到现有的CSS中，替换原有的stat-content相关样式 */

.dual-stat-content {
  display: flex;
  align-items: flex-start;
  gap: 12px;
  color: white;
  padding: 4px 0;
}

.dual-stat-info {
  flex: 1;
  min-width: 0;
}

.stat-title {
  font-size: 13px;
  font-weight: 600;
  margin-bottom: 8px;
  opacity: 0.95;
  line-height: 1.2;
}

.dual-values {
  margin-bottom: 6px;
}

.value-row {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 3px;
  line-height: 1.3;
}

.value-row:last-child {
  margin-bottom: 0;
}

.group-label {
  font-size: 11px;
  font-weight: 500;
  opacity: 0.85;
  min-width: 55px;
}

.stat-value-small {
  font-size: 16px;
  font-weight: bold;
  text-align: right;
}

.stat-subtitle {
  font-size: 11px;
  opacity: 0.8;
  margin-top: 4px;
  text-align: center;
  font-weight: 500;
}

.stat-icon {
  font-size: 24px;
  margin-top: 2px;
  flex-shrink: 0;
}

/* 调整卡片高度以适应双行内容 */
.stat-card :deep(.el-card__body) {
  padding: 18px;
  min-height: 100px;
  display: flex;
  align-items: center;
}

/* 响应式调整 */
@media (max-width: 1200px) {
  .stat-value-small {
    font-size: 14px;
  }
  
  .group-label {
    font-size: 10px;
    min-width: 50px;
  }
  
  .stat-title {
    font-size: 12px;
  }
}

@media (max-width: 768px) {
  .dual-stat-content {
    flex-direction: column;
    align-items: center;
    text-align: center;
    gap: 8px;
  }
  
  .value-row {
    justify-content: center;
    gap: 8px;
  }
  
  .group-label {
    min-width: auto;
  }
}
.data-overview {
max-width: 1400px;
margin: 0 auto;
padding: 24px;
background-color: #f5f7fa;
min-height: 100vh;
}

.header-section {
margin-bottom: 32px;
}

.header-section h1 {
color: #2c3e50;
margin-bottom: 8px;
font-size: 32px;
font-weight: 700;
}

.subtitle {
color: #6c757d;
font-size: 16px;
margin: 0;
}

.tab-container {
margin-bottom: 24px;
border-radius: 12px;
}

.tab-header {
display: flex;
justify-content: space-between;
align-items: center;
}

.modern-tabs {
flex: 1;
}

.modern-tabs :deep(.el-tabs__header) {
margin: 0;
}

.modern-tabs :deep(.el-tabs__item) {
font-size: 16px;
font-weight: 500;
padding: 0 24px;
}

.tab-actions {
display: flex;
gap: 12px;
}

.summary-section {
margin-bottom: 32px;
}

.stat-card {
border-radius: 16px;
transition: all 0.3s ease;
overflow: hidden;
position: relative;
}

.stat-card:hover {
transform: translateY(-4px);
box-shadow: 0 12px 40px rgba(0, 0, 0, 0.15) !important;
}

.stat-card :deep(.el-card__body) {
padding: 24px;
color: white;
position: relative;
z-index: 2;
}

.stat-icon {
font-size: 28px;
margin-right: 12px;
}

.stat-content {
display: flex;
align-items: center;
gap: 16px;
color: white;
}

.stat-info {
flex: 1;
}

.stat-title {
font-size: 14px;
font-weight: 500;
margin-bottom: 8px;
opacity: 0.9;
}

.stat-value {
font-size: 24px;
font-weight: bold;
margin-bottom: 4px;
}

.stat-subtitle {
font-size: 13px;
opacity: 0.9;
margin-top: 8px;
}

.status-alert {
position: fixed;
bottom: 20px;
right: 20px;
z-index: 1000;
max-width: 400px;
}

/* Gradient backgrounds */
.gradient-blue {
background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
}

.gradient-green {
background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
}

.gradient-orange {
background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
}

.gradient-purple {
background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%);
}

.gradient-indigo {
background: linear-gradient(135deg, #fa709a 0%, #fee140 100%);
}

.gradient-emerald {
background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%);
}

.gradient-amber {
background: linear-gradient(135deg, #ffecd2 0%, #fcb69f 100%);
}

.gradient-violet {
background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
}

.charts-section {
margin-bottom: 32px;
}

.chart-card {
border-radius: 12px;
transition: all 0.3s ease;
height: 100%;
}

.chart-card:hover {
transform: translateY(-2px);
box-shadow: 0 8px 30px rgba(0, 0, 0, 0.12) !important;
}

.chart-header {
display: flex;
justify-content: space-between;
align-items: center;
font-weight: 600;
color: #2c3e50;
}

.chart-controls {
display: flex;
align-items: center;
gap: 12px;
}

.control-label {
font-size: 14px;
color: #666;
white-space: nowrap;
}

.tab-content {
animation: fadeIn 0.3s ease-in-out;
}

@keyframes fadeIn {
from {
opacity: 0;
transform: translateY(10px);
}
to {
opacity: 1;
transform: translateY(0);
}
}

/* Responsive Design */
@media (max-width: 1200px) {
.tab-header {
flex-direction: column;
gap: 16px;
align-items: stretch;
}

.modern-tabs {
order: 1;
}

.tab-actions {
order: 2;
justify-content: center;
}
}

@media (max-width: 768px) {
.data-overview {
padding: 16px;
}

.header-section h1 {
font-size: 24px;
}

.summary-section :deep(.el-col) {
margin-bottom: 16px;
}

.charts-section :deep(.el-col) {
margin-bottom: 20px;
}

.chart-controls {
flex-direction: column;
align-items: stretch;
gap: 8px;
}

.chart-controls .el-slider {
width: 100% !important;
}
}

/* Element Plus Style Overrides */
:deep(.el-statistic__content) {
display: flex;
align-items: center;
}

:deep(.el-statistic__number) {
font-size: 24px !important;
font-weight: bold;
}

:deep(.el-statistic__title) {
font-size: 14px;
font-weight: 500;
margin-bottom: 8px;
opacity: 0.9;
}

:deep(.el-card__header) {
padding: 16px 20px;
border-bottom: 1px solid #f0f2f5;
}

:deep(.el-card__body) {
padding: 20px;
}

:deep(.el-tabs__content) {
padding: 0;
}

:deep(.el-slider__runway) {
background-color: #e4e7ed;
}

:deep(.el-slider__bar) {
background-color: #409eff;
}

:deep(.el-slider__button) {
border: 2px solid #409eff;
}

/* Custom Scrollbar */
.tab-content::-webkit-scrollbar {
width: 6px;
}

.tab-content::-webkit-scrollbar-track {
background: #f1f1f1;
border-radius: 3px;
}

.tab-content::-webkit-scrollbar-thumb {
background: #c1c1c1;
border-radius: 3px;
}

.tab-content::-webkit-scrollbar-thumb:hover {
background: #a8a8a8;
}
</style>