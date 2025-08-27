<template>
  <div class="user-manual">
    <!-- 头部导航 -->
    <el-affix :offset="0">
      <div class="nav-header">
        <el-container>
          <el-header class="header">
            <div class="logo">
              <h2>用户手册</h2>
            </div>
            <el-menu 
              mode="horizontal" 
              :default-active="activeSection"
              class="nav-menu"
              @select="scrollToSection"
            >
              <el-menu-item index="quick-start">快速开始</el-menu-item>
              <el-menu-item index="components">组件文档</el-menu-item>
              <el-menu-item index="examples">使用示例</el-menu-item>
              <el-menu-item index="faq">常见问题</el-menu-item>
            </el-menu>
          </el-header>
        </el-container>
      </div>
    </el-affix>

    <el-container class="main-container">
      <!-- 侧边锚点导航 -->
      <el-affix :offset="80" position="right">
        <el-card class="anchor-nav" shadow="hover">
          <el-menu 
            :default-active="activeSection"
            class="anchor-menu"
            @select="scrollToSection"
          >
            <el-menu-item-group title="快速开始">
              <el-menu-item index="installation">1.1 安装</el-menu-item>
              <el-menu-item index="getting-started">1.2 快速开始</el-menu-item>
            </el-menu-item-group>
            
            <el-menu-item-group title="组件文档">
              <el-menu-item index="data-import">2.1 数据引入</el-menu-item>
              <el-menu-item index="data-overview">2.2 数据概述</el-menu-item>
              <el-menu-item index="data-analysis">2.3 数据分析</el-menu-item>
            </el-menu-item-group>
            
            <el-menu-item-group title="其他">
              <el-menu-item index="examples">3 使用示例</el-menu-item>
              <el-menu-item index="faq">4 常见问题</el-menu-item>
            </el-menu-item-group>
          </el-menu>
        </el-card>
      </el-affix>

      <el-main class="content">
        <!-- 1. 快速开始 -->
        <section id="quick-start" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><VideoPlay /></el-icon> 1. 快速开始</h1>
            <el-divider />
            
            <!-- 1.1 安装 -->
            <div id="installation" class="subsection">
              <h2><el-icon><Download /></el-icon> 1.1 安装</h2>
              <el-alert
                title="安装前确保已安装 Node.js 16+ 和 npm"
                type="info"
                :closable="false"
                style="margin-bottom: 20px;"
              />
              
              <el-tabs v-model="installTab" class="install-tabs">
                <el-tab-pane label="npm" name="npm">
                  <el-card class="code-card">
                    <pre><code># 安装依赖
npm install your-component

# 安装 Vue 3 (如果还没有)
npm install vue@next

# 安装 Element Plus
npm install element-plus</code></pre>
                    <el-button 
                      type="primary" 
                      size="small" 
                      @click="copyCode('npm install your-component')"
                      style="margin-top: 10px;"
                    >
                      <el-icon><DocumentCopy /></el-icon>
                      复制
                    </el-button>
                  </el-card>
                </el-tab-pane>
                
                <el-tab-pane label="yarn" name="yarn">
                  <el-card class="code-card">
                    <pre><code># 安装依赖
yarn add your-component

# 安装 Vue 3 (如果还没有)
yarn add vue@next

# 安装 Element Plus
yarn add element-plus</code></pre>
                    <el-button 
                      type="primary" 
                      size="small" 
                      @click="copyCode('yarn add your-component')"
                      style="margin-top: 10px;"
                    >
                      <el-icon><DocumentCopy /></el-icon>
                      复制
                    </el-button>
                  </el-card>
                </el-tab-pane>
                
                <el-tab-pane label="CDN" name="cdn">
                  <el-card class="code-card">
                    <pre><code>&lt;!-- Vue 3 --&gt;
&lt;script src="https://unpkg.com/vue@next"&gt;&lt;/script&gt;

&lt;!-- Element Plus --&gt;
&lt;link rel="stylesheet" href="https://unpkg.com/element-plus/dist/index.css"&gt;
&lt;script src="https://unpkg.com/element-plus"&gt;&lt;/script&gt;

&lt;!-- Your Component --&gt;
&lt;script src="https://unpkg.com/your-component"&gt;&lt;/script&gt;</code></pre>
                  </el-card>
                </el-tab-pane>
              </el-tabs>
            </div>

            <!-- 1.2 快速开始 -->
            <div id="getting-started" class="subsection">
              <h2><el-icon><Lightning /></el-icon> 1.2 快速开始</h2>
              
              <el-steps :active="currentStep" finish-status="success" align-center>
                <el-step title="引入组件" />
                <el-step title="注册组件" />
                <el-step title="使用组件" />
                <el-step title="完成" />
              </el-steps>

              <div style="margin: 30px 0;">
                <el-card class="code-card">
                  <div class="code-header">
                    <span>main.js</span>
                    <el-button text type="primary" @click="copyCode(mainJsCode)">
                      <el-icon><DocumentCopy /></el-icon>
                    </el-button>
                  </div>
                  <pre><code>{{ mainJsCode }}</code></pre>
                </el-card>

                <el-card class="code-card" style="margin-top: 20px;">
                  <div class="code-header">
                    <span>App.vue</span>
                    <el-button text type="primary" @click="copyCode(appVueCode)">
                      <el-icon><DocumentCopy /></el-icon>
                    </el-button>
                  </div>
                  <pre><code>{{ appVueCode }}</code></pre>
                </el-card>
              </div>

              <el-button-group>
                <el-button 
                  type="primary" 
                  @click="currentStep = Math.min(currentStep + 1, 3)"
                  :disabled="currentStep >= 3"
                >
                  下一步
                </el-button>
                <el-button 
                  @click="currentStep = Math.max(currentStep - 1, 0)"
                  :disabled="currentStep <= 0"
                >
                  上一步
                </el-button>
              </el-button-group>
            </div>
          </el-card>
        </section>

        <!-- 2. 组件文档 -->
        <section id="components" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><Grid /></el-icon> 2. 组件文档</h1>
            <el-divider />

            <!-- 2.1 数据引入 -->
            <div id="data-import" class="subsection">
              <h2><el-icon><Upload /></el-icon> 2.1 数据引入</h2>
              
              <el-row :gutter="20">
                <el-col :span="12">
                  <el-card class="demo-card">
                    <template #header>
                      <span>支持的数据格式</span>
                    </template>
                    <el-tag v-for="format in dataFormats" :key="format" style="margin: 5px;">
                      {{ format }}
                    </el-tag>
                  </el-card>
                </el-col>
                
                <el-col :span="12">
                  <el-card class="demo-card">
                    <template #header>
                      <span>数据引入示例</span>
                    </template>
                    <el-upload
                      class="upload-demo"
                      drag
                      action="#"
                      :auto-upload="false"
                      :show-file-list="false"
                      @change="handleFileChange"
                    >
                      <el-icon class="el-icon--upload"><Upload /></el-icon>
                      <div class="el-upload__text">
                        拖拽文件到此处或<em>点击上传</em>
                      </div>
                      <template #tip>
                        <div class="el-upload__tip">
                          支持 CSV、Excel、JSON 格式
                        </div>
                      </template>
                    </el-upload>
                  </el-card>
                </el-col>
              </el-row>

              <el-collapse v-model="activeCollapse" style="margin-top: 20px;">
                <el-collapse-item title="API 数据引入" name="api">
                  <el-card class="code-card">
                    <pre><code>// API 数据获取示例
import { ref, onMounted } from 'vue'

const data = ref([])

const fetchData = async () => {
  try {
    const response = await fetch('/api/data')
    data.value = await response.json()
  } catch (error) {
    console.error('数据获取失败:', error)
  }
}

onMounted(() => {
  fetchData()
})</code></pre>
                  </el-card>
                </el-collapse-item>
                
                <el-collapse-item title="本地数据引入" name="local">
                  <el-card class="code-card">
                    <pre><code>// 本地数据引入示例
import localData from './data/sample.json'

const data = ref(localData)

// 或者动态引入
const loadLocalData = async () => {
  const { default: data } = await import('./data/sample.json')
  return data
}</code></pre>
                  </el-card>
                </el-collapse-item>
              </el-collapse>
            </div>

            <!-- 2.2 数据概述 -->
            <div id="data-overview" class="subsection">
              <h2><el-icon><View /></el-icon> 2.2 数据概述</h2>
              
              <el-row :gutter="20">
                <el-col :span="8" v-for="(stat, index) in dataStats" :key="index">
                  <el-card class="stat-card" shadow="hover">
                    <el-statistic 
                      :title="stat.title"
                      :value="stat.value"
                      :suffix="stat.suffix"
                    >
                      <template #prefix>
                        <el-icon :style="{ color: stat.color }">
                          <component :is="stat.icon" />
                        </el-icon>
                      </template>
                    </el-statistic>
                  </el-card>
                </el-col>
              </el-row>

              <!-- 数据表格预览 -->
              <el-card style="margin-top: 20px;">
                <template #header>
                  <div class="card-header">
                    <span>数据预览</span>
                    <el-button-group>
                      <el-button size="small" @click="refreshData">
                        <el-icon><Refresh /></el-icon>
                        刷新
                      </el-button>
                      <el-button size="small" @click="exportData">
                        <el-icon><Download /></el-icon>
                        导出
                      </el-button>
                    </el-button-group>
                  </div>
                </template>
                
                <el-table :data="sampleData" border style="width: 100%">
                  <el-table-column prop="id" label="ID" width="80" />
                  <el-table-column prop="name" label="名称" width="120" />
                  <el-table-column prop="value" label="数值" width="100" />
                  <el-table-column prop="status" label="状态" width="100">
                    <template #default="{ row }">
                      <el-tag :type="getStatusType(row.status)">{{ row.status }}</el-tag>
                    </template>
                  </el-table-column>
                  <el-table-column prop="date" label="日期" />
                </el-table>
              </el-card>
            </div>

            <!-- 2.3 数据分析 -->
            <div id="data-analysis" class="subsection">
              <h2><el-icon><TrendCharts /></el-icon> 2.3 数据分析</h2>
              
              <el-tabs v-model="analysisTab">
                <el-tab-pane label="图表分析" name="chart">
                  <el-card>
                    <div ref="chartRef" style="height: 300px;"></div>
                    <el-button-group style="margin-top: 15px;">
                      <el-button @click="changeChart('bar')" :type="chartType === 'bar' ? 'primary' : ''">
                        柱状图
                      </el-button>
                      <el-button @click="changeChart('line')" :type="chartType === 'line' ? 'primary' : ''">
                        折线图
                      </el-button>
                      <el-button @click="changeChart('pie')" :type="chartType === 'pie' ? 'primary' : ''">
                        饼图
                      </el-button>
                    </el-button-group>
                  </el-card>
                </el-tab-pane>
                
                <el-tab-pane label="统计分析" name="stats">
                  <el-row :gutter="20">
                    <el-col :span="12">
                      <el-card>
                        <el-descriptions title="基础统计" :column="2" border>
                          <el-descriptions-item label="平均值">{{ stats.average }}</el-descriptions-item>
                          <el-descriptions-item label="最大值">{{ stats.max }}</el-descriptions-item>
                          <el-descriptions-item label="最小值">{{ stats.min }}</el-descriptions-item>
                          <el-descriptions-item label="总计">{{ stats.sum }}</el-descriptions-item>
                        </el-descriptions>
                      </el-card>
                    </el-col>
                    <el-col :span="12">
                      <el-card>
                        <el-progress 
                          type="dashboard" 
                          :percentage="75" 
                          :color="colors"
                        >
                          <template #default="{ percentage }">
                            <span class="percentage-value">{{ percentage }}%</span>
                            <span class="percentage-label">完成度</span>
                          </template>
                        </el-progress>
                      </el-card>
                    </el-col>
                  </el-row>
                </el-tab-pane>
              </el-tabs>
            </div>
          </el-card>
        </section>

        <!-- 3. 使用示例 -->
        <section id="examples" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><View /></el-icon> 3. 使用示例</h1>
            <el-divider />
            
            <el-row :gutter="20">
              <el-col :span="8" v-for="(example, index) in examples" :key="index">
                <el-card class="example-card" shadow="hover">
                  <template #header>
                    <div class="card-header">
                      <span>{{ example.title }}</span>
                      <el-tag size="small">{{ example.difficulty }}</el-tag>
                    </div>
                  </template>
                  
                  <p>{{ example.description }}</p>
                  
                  <el-button type="primary" text @click="showExampleCode(example)">
                    查看代码 <el-icon><ArrowRight /></el-icon>
                  </el-button>
                </el-card>
              </el-col>
            </el-row>

            <!-- 代码示例对话框 -->
            <el-dialog 
              v-model="exampleDialog" 
              :title="selectedExample?.title"
              width="70%"
              top="5vh"
            >
              <el-tabs v-model="codeTab">
                <el-tab-pane label="Template" name="template">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.template }}</code></pre>
                  </el-card>
                </el-tab-pane>
                <el-tab-pane label="Script" name="script">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.script }}</code></pre>
                  </el-card>
                </el-tab-pane>
                <el-tab-pane label="Style" name="style">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.style }}</code></pre>
                  </el-card>
                </el-tab-pane>
              </el-tabs>
              
              <template #footer>
                <el-button @click="exampleDialog = false">关闭</el-button>
                <el-button type="primary" @click="copyFullExample">
                  <el-icon><DocumentCopy /></el-icon>
                  复制完整代码
                </el-button>
              </template>
            </el-dialog>
          </el-card>
        </section>

        <!-- 4. 常见问题 -->
        <section id="faq" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><QuestionFilled /></el-icon> 4. 常见问题</h1>
            <el-divider />
            
            <el-collapse v-model="activeFaq" accordion>
              <el-collapse-item 
                v-for="(faq, index) in faqs" 
                :key="index"
                :title="faq.question"
                :name="index"
              >
                <div class="faq-content">
                  <p>{{ faq.answer }}</p>
                  <el-card v-if="faq.code" class="code-card">
                    <pre><code>{{ faq.code }}</code></pre>
                  </el-card>
                </div>
              </el-collapse-item>
            </el-collapse>

            <!-- 联系支持 -->
            <el-alert
              title="还有其他问题？"
              type="info"
              description="如果以上内容无法解决您的问题，请联系我们的技术支持团队"
              :closable="false"
              style="margin-top: 30px;"
            >
              <template #default>
                <el-button-group>
                  <el-button type="primary">
                    <el-icon><Message /></el-icon>
                    在线客服
                  </el-button>
                  <el-button>
                    <el-icon><Phone /></el-icon>
                    技术支持
                  </el-button>
                  <el-button>
                    <el-icon><Document /></el-icon>
                    提交反馈
                  </el-button>
                </el-button-group>
              </template>
            </el-alert>
          </el-card>
        </section>
      </el-main>
    </el-container>

    <!-- 回到顶部 -->
    <el-backtop :right="100" :bottom="100" />
  </div>
</template>

<script setup>
import { ref, reactive, onMounted, nextTick, getCurrentInstance } from 'vue'
import {
  VideoPlay, Download, Lightning, Grid, Upload, View, TrendCharts,
  Refresh, DocumentCopy, ArrowRight, QuestionFilled, Message, Phone, Document
} from '@element-plus/icons-vue'

// 获取当前实例以使用 $echarts
const { proxy } = getCurrentInstance()

// 响应式数据
const activeSection = ref('quick-start')
const installTab = ref('npm')
const currentStep = ref(0)
const activeCollapse = ref([])
const analysisTab = ref('chart')
const chartType = ref('bar')
const exampleDialog = ref(false)
const selectedExample = ref(null)
const codeTab = ref('template')
const activeFaq = ref(0)

// 图表引用
const chartRef = ref(null)
let chart = null

// 数据
const dataFormats = ['CSV', 'Excel', 'JSON', 'XML', 'API']

const dataStats = reactive([
  { title: '数据条数', value: 1234, suffix: '条', icon: 'Grid', color: '#409EFF' },
  { title: '处理速度', value: 95, suffix: '%', icon: 'TrendCharts', color: '#67C23A' },
  { title: '准确率', value: 99.2, suffix: '%', icon: 'View', color: '#E6A23C' }
])

const sampleData = reactive([
  { id: 1, name: '示例数据1', value: 100, status: '正常', date: '2024-01-01' },
  { id: 2, name: '示例数据2', value: 200, status: '警告', date: '2024-01-02' },
  { id: 3, name: '示例数据3', value: 150, status: '错误', date: '2024-01-03' },
  { id: 4, name: '示例数据4', value: 300, status: '正常', date: '2024-01-04' }
])

const stats = reactive({
  average: 187.5,
  max: 300,
  min: 100,
  sum: 750
})

const colors = [
  { color: '#f56565', percentage: 20 },
  { color: '#ed8936', percentage: 40 },
  { color: '#ecc94b', percentage: 60 },
  { color: '#48bb78', percentage: 80 },
  { color: '#38b2ac', percentage: 100 }
]

const examples = [
  {
    title: '基础用法',
    difficulty: '简单',
    description: '展示组件的基本使用方法',
    code: {
      template: `<template>
  <div>
    <your-component :data="data" />
  </div>
</template>`,
      script: `<script setup>
import { ref } from 'vue'

const data = ref([1, 2, 3, 4, 5])
<\/script>`,
      style: `<style scoped>
div {
  padding: 20px;
}
</style>`
    }
  },
  {
    title: '高级配置',
    difficulty: '中等',
    description: '使用高级配置选项',
    code: {
      template: `<template>
  <div>
    <your-component 
      :data="data"
      :options="options"
      @change="handleChange"
    />
  </div>
</template>`,
      script: `<script setup>
import { ref, reactive } from 'vue'

const data = ref([1, 2, 3, 4, 5])
const options = reactive({
  theme: 'dark',
  animation: true
})

const handleChange = (value) => {
  console.log('Changed:', value)
}
<\/script>`,
      style: `<style scoped>
div {
  padding: 20px;
  background: #f5f7fa;
}
</style>`
    }
  },
  {
    title: '自定义主题',
    difficulty: '困难',
    description: '创建和使用自定义主题',
    code: {
      template: `<template>
  <div class="custom-theme">
    <your-component 
      :data="data"
      :theme="customTheme"
    />
  </div>
</template>`,
      script: `<script setup>
import { ref, computed } from 'vue'

const data = ref([1, 2, 3, 4, 5])

const customTheme = computed(() => ({
  colors: ['#ff6b6b', '#4ecdc4', '#45b7d1'],
  fontSize: 14,
  borderRadius: 8
}))
<\/script>`,
      style: `<style scoped>
.custom-theme {
  padding: 20px;
  border: 1px solid #ddd;
  border-radius: 8px;
}
</style>`
    }
  }
]

const faqs = [
  {
    question: '如何解决安装依赖失败的问题？',
    answer: '通常是网络问题导致的，建议切换到国内镜像源或使用代理。',
    code: `# 使用淘宝镜像
npm config set registry https://registry.npmmirror.com/

# 或临时使用
npm install --registry https://registry.npmmirror.com/`
  },
  {
    question: '组件不显示是什么原因？',
    answer: '请检查是否正确引入了样式文件和组件注册。',
    code: `// main.js 中确保引入样式
import 'element-plus/dist/index.css'
import './your-component.css'`
  },
  {
    question: '如何自定义组件样式？',
    answer: '可以通过 CSS 变量或者深度选择器来覆盖默认样式。',
    code: `:deep(.your-component) {
  --primary-color: #409EFF;
  --border-radius: 4px;
}`
  },
  {
    question: '数据更新后组件不刷新怎么办？',
    answer: '确保使用响应式数据，并检查是否正确使用了 ref 或 reactive。',
    code: `// 正确的响应式数据使用
const data = ref([])
const updateData = () => {
  data.value = newData // 这样会触发更新
}`
  }
]

// 代码字符串
const mainJsCode = `import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import { createPinia } from 'pinia'

// 导入Element Plus
import ElementPlus from 'element-plus'
import 'element-plus/dist/index.css'

// 导入你的组件
import YourComponent from './components/YourComponent.vue'

const app = createApp(App)

app.use(ElementPlus)
app.use(createPinia())
app.use(router)

// 注册全局组件
app.component('YourComponent', YourComponent)

app.mount('#app')`

const appVueCode = `<template>
  <div id="app">
    <your-component 
      :data="sampleData"
      @change="handleDataChange"
    />
  </div>
</template>

<script setup>
import { ref } from 'vue'

const sampleData = ref([
  { name: '数据1', value: 100 },
  { name: '数据2', value: 200 },
  { name: '数据3', value: 150 }
])

const handleDataChange = (newData) => {
  console.log('数据变化:', newData)
}
<\/script>`

// 方法
const scrollToSection = (key) => {
  activeSection.value = key
  const element = document.getElementById(key)
  if (element) {
    element.scrollIntoView({ behavior: 'smooth', block: 'start' })
  }
}

const copyCode = (code) => {
  navigator.clipboard.writeText(code).then(() => {
    proxy.$message.success('代码已复制到剪贴板')
  })
}

const handleFileChange = (file) => {
  console.log('文件上传:', file)
  proxy.$message.info(`已选择文件: ${file.name}`)
}

const getStatusType = (status) => {
  const typeMap = {
    '正常': 'success',
    '警告': 'warning',
    '错误': 'danger'
  }
  return typeMap[status] || 'info'
}

const refreshData = () => {
  proxy.$message.success('数据已刷新')
}

const exportData = () => {
  proxy.$message.success('数据导出成功')
}

const initChart = () => {
  if (chartRef.value && proxy.$echarts) {
    chart = proxy.$echarts.init(chartRef.value)
    changeChart('bar')
  }
}

const changeChart = (type) => {
  chartType.value = type
  if (!chart) return
  
  const data = [120, 200, 150, 80, 70]
  const xData = ['周一', '周二', '周三', '周四', '周五']
  
  let option = {}
  
  if (type === 'pie') {
    option = {
      title: { text: '饼图示例', left: 'center' },
      tooltip: { trigger: 'item' },
      series: [{
        type: 'pie',
        radius: '60%',
        center: ['50%', '50%'],
        data: xData.map((name, index) => ({ name, value: data[index] }))
      }]
    }
  } else {
    option = {
      title: { text: type === 'bar' ? '柱状图示例' : '折线图示例' },
      tooltip: { trigger: 'axis' },
      xAxis: {
        type: 'category',
        data: xData
      },
      yAxis: { type: 'value' },
      series: [{
        type: type,
        data: data,
        smooth: type === 'line'
      }]
    }
  }
  
  chart.setOption(option)
}

const showExampleCode = (example) => {
  selectedExample.value = example
  exampleDialog.value = true
  codeTab.value = 'template'
}

const copyFullExample = () => {
  const example = selectedExample.value
  const fullCode = `${example.code.template}

${example.code.script}

${example.code.style}`
  copyCode(fullCode)
}

// 生命周期
onMounted(() => {
  nextTick(() => {
    initChart()
  })
  
  // 监听滚动事件更新活跃导航
  window.addEventListener('scroll', () => {
    const sections = ['quick-start', 'components', 'examples', 'faq']
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop
    
    for (let i = sections.length - 1; i >= 0; i--) {
      const element = document.getElementById(sections[i])
      if (element && element.offsetTop - 100 <= scrollTop) {
        activeSection.value = sections[i]
        break
      }
    }
  })
})
</script>

<style scoped>
.user-manual {
  min-height: 100vh;
  background-color: #f5f7fa;
}

.nav-header {
  background: white;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
  z-index: 1000;
}

.header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  max-width: 1200px;
  margin: 0 auto;
  padding: 0 20px;
}

.logo h2 {
  margin: 0;
  color: #409EFF;
}

.nav-menu {
  border-bottom: none;
}

.main-container {
  max-width: 1200px;
  margin: 0 auto;
  background: transparent;
}

.anchor-nav {
  width: 200px;
  max-height: 400px;
  overflow-y: auto;
}

.anchor-menu {
  border: none;
}

.content {
  padding: 20px;
}

.section {
  margin-bottom: 60px;
}

.section-card {
  border: none;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

.section h1 {
  display: flex;
  align-items: center;
  gap: 10px;
  color: #303133;
  margin-bottom: 20px;
}

.subsection {
  margin: 40px 0;
}

.subsection h2 {
  display: flex;
  align-items: center;
  gap: 8px;
  color: #606266;
  margin-bottom: 20px;
  padding-bottom: 10px;
  border-bottom: 2px solid #E4E7ED;
}

.install-tabs {
  margin: 20px 0;
}

.code-card {
  background-color: #fafafa;
  border: 1px solid #e4e7ed;
}

.code-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 10px;
  font-weight: 500;
  color: #606266;
}

.code-card pre {
  margin: 0;
  padding: 0;
  background: transparent;
  font-family: 'Monaco', 'Consolas', 'Courier New', monospace;
  font-size: 14px;
  line-height: 1.6;
  color: #2c3e50;
  overflow-x: auto;
}

.demo-card {
  height: 100%;
}

.upload-demo {
  width: 100%;
}

.stat-card {
  text-align: center;
  height: 100%;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.example-card {
  height: 100%;
  cursor: pointer;
  transition: transform 0.2s;
}

.example-card:hover {
  transform: translateY(-5px);
}

.faq-content {
  padding: 10px 0;
}

.percentage-value {
  display: block;
  margin-top: 10px;
  font-size: 28px;
  font-weight: bold;
  color: #409EFF;
}

.percentage-label {
  display: block;
  margin-top: 5px;
  font-size: 12px;
  color: #909399;
}

/* 响应式设计 */
@media (max-width: 768px) {
  .header {
    flex-direction: column;
    padding: 10px;
  }
  
  .nav-menu {
    width: 100%;
    justify-content: center;
  }
  
  .anchor-nav {
    display: none;
  }
  
  .content {
    padding: 10px;
  }
}

/* 滚动条样式 */
::-webkit-scrollbar {
  width: 6px;
}

::-webkit-scrollbar-track {
  background: #f1f1f1;
}

::-webkit-scrollbar-thumb {
  background: #c1c1c1;
  border-radius: 3px;
}

::-webkit-scrollbar-thumb:hover {
  background: #a1a1a1;
}
</style>