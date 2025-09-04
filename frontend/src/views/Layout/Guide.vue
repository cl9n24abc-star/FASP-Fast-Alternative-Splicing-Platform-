<template>
  <div class="user-manual">
    <!-- Header Navigation -->
    <el-affix :offset="0">
      <div class="nav-header">
        <el-container>
          <el-header class="header">
            <div class="logo">
              <h2>FASP User Manual</h2>
            </div>
            <el-menu 
              mode="horizontal" 
              :default-active="activeSection"
              class="nav-menu"
              @select="scrollToSection"
            >
              <el-menu-item index="quick-start">Quick Start</el-menu-item>
              <el-menu-item index="components">Features</el-menu-item>
              <el-menu-item index="examples">Examples</el-menu-item>
              <el-menu-item index="faq">FAQ</el-menu-item>
            </el-menu>
          </el-header>
        </el-container>
      </div>
    </el-affix>

    <el-container class="main-container">
      <el-main class="content">
        <!-- Side Anchor Navigation -->
        <div class="anchor-nav-wrapper">
          <el-card class="anchor-nav" shadow="hover">
            <el-menu 
              :default-active="activeSection"
              class="anchor-menu"
              @select="scrollToSection"
            >
              <el-menu-item-group title="Quick Start">
                <el-menu-item index="installation">1.1 Installation</el-menu-item>
                <el-menu-item index="getting-started">1.2 Getting Started</el-menu-item>
              </el-menu-item-group>
              
              <el-menu-item-group title="Features">
                <el-menu-item index="data-import">2.1 Data Import</el-menu-item>
                <el-menu-item index="data-overview">2.2 Data Overview</el-menu-item>
                <el-menu-item index="data-analysis">2.3 Analysis</el-menu-item>
              </el-menu-item-group>
              
              <el-menu-item-group title="Others">
                <el-menu-item index="examples">3 Examples</el-menu-item>
                <el-menu-item index="faq">4 FAQ</el-menu-item>
              </el-menu-item-group>
            </el-menu>
          </el-card>
        </div>
        <!-- 1. Quick Start -->
        <section id="quick-start" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><VideoPlay /></el-icon> 1. Quick Start</h1>
            <el-divider />
            
            <!-- 1.1 Installation -->
            <div id="installation" class="subsection">
              <h2><el-icon><Download /></el-icon> 1.1 Installation</h2>
              
              <el-alert
                title="Prerequisites"
                type="info"
                :closable="false"
                style="margin-bottom: 20px;"
              >
                <p>Before installing FASP, ensure you have:</p>
                <ul>
                  <li><strong>Node.js 20.19+</strong> (for frontend)</li>
                  <li><strong>Python 3.8+</strong> (for backend)</li>
                  <li><strong>npm</strong> or <strong>yarn</strong></li>
                  <li><strong>Git</strong></li>
                </ul>
              </el-alert>
              
              <el-tabs v-model="installTab" class="install-tabs">
                <el-tab-pane label="Git Clone" name="git">
                  <el-card class="code-card">
                    <pre><code># Clone the repository
git clone https://github.com/cl9n24abc-star/FASP-Fast-Alternative-Splicing-Platform-.git
cd FASP-Fast-Alternative-Splicing-Platform-

# Setup backend
cd backend
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate    # Windows
pip install -r requirements.txt

# Setup frontend (new terminal)
cd frontend
npm install</code></pre>
                    <el-button 
                      type="primary" 
                      size="small" 
                      @click="copyCode(gitInstallCode)"
                      style="margin-top: 10px;"
                    >
                      <el-icon><DocumentCopy /></el-icon>
                      Copy
                    </el-button>
                  </el-card>
                </el-tab-pane>
                
                <el-tab-pane label="Docker" name="docker">
                  <el-card class="code-card">
                    <pre><code># Using Docker (if available)
docker-compose up -d

# Or build manually:
docker build -t fasp-backend ./backend
docker build -t fasp-frontend ./frontend

# Run containers:
docker run -p 5000:5000 fasp-backend
docker run -p 5173:5173 fasp-frontend</code></pre>
                    <el-button 
                      type="primary" 
                      size="small" 
                      @click="copyCode(dockerInstallCode)"
                      style="margin-top: 10px;"
                    >
                      <el-icon><DocumentCopy /></el-icon>
                      Copy
                    </el-button>
                  </el-card>
                </el-tab-pane>
                
                <el-tab-pane label="Check Environment" name="check">
                  <el-card class="code-card">
                    <pre><code># Check Node.js version (must be 20.19+)
node --version

# Check Python version
python --version

# Check npm version
npm --version

# Verify Git installation
git --version</code></pre>
                    <el-button 
                      type="primary" 
                      size="small" 
                      @click="copyCode(checkEnvCode)"
                      style="margin-top: 10px;"
                    >
                      <el-icon><DocumentCopy /></el-icon>
                      Copy
                    </el-button>
                  </el-card>
                </el-tab-pane>
              </el-tabs>
            </div>

            <!-- 1.2 Getting Started -->
            <div id="getting-started" class="subsection">
              <h2><el-icon><Lightning /></el-icon> 1.2 Getting Started</h2>
              
              <el-steps :active="currentStep" finish-status="success" align-center>
                <el-step title="Start Backend" />
                <el-step title="Start Frontend" />
                <el-step title="Access Application" />
                <el-step title="Complete" />
              </el-steps>

              <div style="margin: 30px 0;">
                <el-card class="code-card">
                  <div class="code-header">
                    <span>Backend Server</span>
                    <el-button text type="primary" @click="copyCode(backendStartCode)">
                      <el-icon><DocumentCopy /></el-icon>
                    </el-button>
                  </div>
                  <pre><code>{{ backendStartCode }}</code></pre>
                </el-card>

                <el-card class="code-card" style="margin-top: 20px;">
                  <div class="code-header">
                    <span>Frontend Server</span>
                    <el-button text type="primary" @click="copyCode(frontendStartCode)">
                      <el-icon><DocumentCopy /></el-icon>
                    </el-button>
                  </div>
                  <pre><code>{{ frontendStartCode }}</code></pre>
                </el-card>

                <el-alert
                  title="Access Your Application"
                  type="success"
                  :closable="false"
                  style="margin-top: 20px;"
                >
                  <p><strong>Frontend:</strong> <a href="http://localhost:5173" target="_blank">http://localhost:5173</a></p>
                  <p><strong>Backend API:</strong> <a href="http://localhost:5000" target="_blank">http://localhost:5000</a></p>
                </el-alert>
              </div>

              <el-button-group>
                <el-button 
                  type="primary" 
                  @click="currentStep = Math.min(currentStep + 1, 3)"
                  :disabled="currentStep >= 3"
                >
                  Next Step
                </el-button>
                <el-button 
                  @click="currentStep = Math.max(currentStep - 1, 0)"
                  :disabled="currentStep <= 0"
                >
                  Previous Step
                </el-button>
              </el-button-group>
            </div>
          </el-card>
        </section>

        <!-- 2. Features -->
        <section id="components" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><Grid /></el-icon> 2. Features</h1>
            <el-divider />

            <!-- 2.1 Data Import -->
            <div id="data-import" class="subsection">
              <h2><el-icon><Upload /></el-icon> 2.1 Data Import</h2>
              
              <el-card class="demo-card">
                <template #header>
                  <span>Supported Data Formats</span>
                </template>
                <el-tag v-for="format in dataFormats" :key="format" style="margin: 5px;">
                  {{ format }}
                </el-tag>
              </el-card>

              <el-collapse v-model="activeCollapse" style="margin-top: 20px;">
                <el-collapse-item title="rMATS Output Files" name="rmats">
                  <el-card class="code-card">
                    <pre><code># Expected rMATS directory structure:
rmats_output/
├── SE.MATS.JC.txt     # Skipped Exon
├── RI.MATS.JC.txt     # Retained Intron
├── A3SS.MATS.JC.txt   # Alternative 3' Splice Site
├── A5SS.MATS.JC.txt   # Alternative 5' Splice Site
└── MXE.MATS.JC.txt    # Mutually Exclusive Exons</code></pre>
                  </el-card>
                </el-collapse-item>
                
                <el-collapse-item title="BAM Files Preparation" name="bam">
                  <el-card class="code-card">
                    <pre><code># Prepare BAM files for sashimi plots
# 1. Sort BAM file
samtools sort input.bam -o sorted.bam

# 2. Index BAM file  
samtools index sorted.bam

# 3. Verify BAM file
samtools view -H sorted.bam | head</code></pre>
                  </el-card>
                </el-collapse-item>
              </el-collapse>
            </div>
          </el-card>
        </section>

        <!-- 3. Examples -->
        <section id="examples" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><View /></el-icon> 3. Usage Examples</h1>
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
                    View Code <el-icon><ArrowRight /></el-icon>
                  </el-button>
                </el-card>
              </el-col>
            </el-row>

            <!-- Code Example Dialog -->
            <el-dialog 
              v-model="exampleDialog" 
              :title="selectedExample?.title"
              width="70%"
              top="5vh"
            >
              <el-tabs v-model="codeTab">
                <el-tab-pane label="Backend" name="backend">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.backend }}</code></pre>
                  </el-card>
                </el-tab-pane>
                <el-tab-pane label="Frontend" name="frontend">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.frontend }}</code></pre>
                  </el-card>
                </el-tab-pane>
                <el-tab-pane label="Config" name="config">
                  <el-card class="code-card">
                    <pre><code>{{ selectedExample?.code.config }}</code></pre>
                  </el-card>
                </el-tab-pane>
              </el-tabs>
              
              <template #footer>
                <el-button @click="exampleDialog = false">Close</el-button>
                <el-button type="primary" @click="copyFullExample">
                  <el-icon><DocumentCopy /></el-icon>
                  Copy Full Code
                </el-button>
              </template>
            </el-dialog>
          </el-card>
        </section>

        <!-- 4. FAQ -->
        <section id="faq" class="section">
          <el-card shadow="never" class="section-card">
            <h1><el-icon><QuestionFilled /></el-icon> 4. Frequently Asked Questions</h1>
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

            <!-- Contact Support -->
            <el-alert
              title="Need More Help?"
              type="info"
              description="If the above content cannot solve your problem, please contact our technical support team"
              :closable="false"
              style="margin-top: 30px;"
            >
              <template #default>
                <el-button-group>
                  <el-button type="primary" @click="openGitHubIssues">
                    <el-icon><Message /></el-icon>
                    GitHub Issues
                  </el-button>
                  <el-button @click="openDocumentation">
                    <el-icon><Document /></el-icon>
                    Documentation
                  </el-button>
                </el-button-group>
              </template>
            </el-alert>
          </el-card>
        </section>
      </el-main>
    </el-container>

    <!-- Back to Top -->
    <el-backtop :right="100" :bottom="100" />
  </div>
</template>

<script setup>
import { ref, reactive, onMounted, nextTick, getCurrentInstance } from 'vue'
import {
  VideoPlay, Download, Lightning, Grid, Upload, View, TrendCharts,
  Refresh, DocumentCopy, ArrowRight, QuestionFilled, Message, Document
} from '@element-plus/icons-vue'

// Get current instance
const { proxy } = getCurrentInstance()

// Reactive data
const activeSection = ref('quick-start')
const installTab = ref('git')
const currentStep = ref(0)
const activeCollapse = ref([])
const analysisTab = ref('splicing')
const chartType = ref('bar')
const exampleDialog = ref(false)
const selectedExample = ref(null)
const codeTab = ref('backend')
const activeFaq = ref(0)

// Chart reference
const chartRef = ref(null)
let chart = null

// Data for bioinformatics platform
const dataFormats = ['BAM files', 'rMATS output (.txt)', 'GTF/GFF annotation', 'BAM list files (.txt)', 'JSON results']

const dataStats = reactive([
  { title: 'Splicing Events', value: 15420, suffix: 'events', icon: 'Grid', color: '#409EFF' },
  { title: 'Analysis Speed', value: 95, suffix: '%', icon: 'TrendCharts', color: '#67C23A' },
  { title: 'Accuracy Rate', value: 99.2, suffix: '%', icon: 'View', color: '#E6A23C' }
])

const sampleData = reactive([
  { id: 1, event: 'SE_12345', type: 'Skipped Exon', pvalue: 0.001, status: 'Significant', chromosome: 'chr1' },
  { id: 2, event: 'RI_67890', type: 'Retained Intron', pvalue: 0.045, status: 'Moderate', chromosome: 'chr2' },
  { id: 3, event: 'A5SS_11111', type: "Alt 5' SS", pvalue: 0.8, status: 'Not Significant', chromosome: 'chr3' },
  { id: 4, event: 'A3SS_22222', type: "Alt 3' SS", pvalue: 0.02, status: 'Significant', chromosome: 'chr4' }
])

const stats = reactive({
  totalEvents: 15420,
  significant: 3280,
  lowPvalue: 2890,
  completion: 95
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
    title: 'Basic rMATS Analysis',
    difficulty: 'Beginner',
    description: 'Analyze alternative splicing events from rMATS output files',
    code: {
      backend: `# Backend setup for rMATS analysis
cd backend
source venv/bin/activate
python app.py`,
      frontend: `# Frontend development server
cd frontend
npm install
npm run dev`,
      config: `# Configuration example
{
  "rmats_dir": "/path/to/rmats/output",
  "bam_files": {
    "sample1": "/path/to/sample1.bam",
    "sample2": "/path/to/sample2.bam"
  },
  "gff_file": "/path/to/annotation.gtf"
}`
    }
  },
  {
    title: 'Sashimi Plot Generation',
    difficulty: 'Intermediate',
    description: 'Generate sashimi plots for splice junction visualization',
    code: {
      backend: `# Sashimi plot API endpoint
@app.route('/api/sashimi/generate', methods=['POST'])
def generate_sashimi():
    data = request.get_json()
    result = sashimi_service.generate_plot(data)
    return jsonify(result)`,
      frontend: `// Frontend sashimi plot component
const generateSashimi = async (coordinates) => {
  const response = await fetch('/api/sashimi/generate', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(coordinates)
  })
  return await response.json()
}`,
      config: `# Sashimi configuration
{
  "coordinates": "chr1:1000-2000",
  "event_type": "SE",
  "output_format": "pdf"
}`
    }
  }
]

const faqs = [
  {
    question: 'How do I prepare my rMATS output files?',
    answer: 'Ensure your rMATS output directory contains all necessary files: SE.MATS.JC.txt, RI.MATS.JC.txt, etc. The directory structure should match the expected format.',
    code: `# Expected rMATS directory structure:
rmats_output/
├── SE.MATS.JC.txt
├── RI.MATS.JC.txt  
├── A3SS.MATS.JC.txt
├── A5SS.MATS.JC.txt
└── MXE.MATS.JC.txt`
  },
  {
    question: 'Why are my BAM files not being recognized?',
    answer: 'Check that BAM files are properly indexed (.bai files present) and paths are correct. Ensure BAM files are sorted and contain proper read groups.',
    code: `# Index your BAM files:
samtools index sample.bam

# Verify BAM file:
samtools view -H sample.bam | head`
  },
  {
    question: 'Node.js version compatibility issues?',
    answer: 'FASP requires Node.js 20.19+ due to Vite 7. Install the correct version using nvm or update your Node.js.',
    code: `# Install Node.js 20+ using nvm
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
source ~/.bashrc
nvm install 20
nvm use 20`
  },
  {
    question: 'How do I interpret the p-values in results?',
    answer: 'P-values indicate the statistical significance of splicing differences. Values < 0.05 are typically considered significant, but adjust for multiple testing.',
    code: `# Filter significant events:
significant_events = results[results['PValue'] < 0.05]
print(f"Found {len(significant_events)} significant events")`
  }
]

// Installation code strings
const gitInstallCode = `# Clone the repository
git clone https://github.com/cl9n24abc-star/FASP-Fast-Alternative-Splicing-Platform-.git
cd FASP-Fast-Alternative-Splicing-Platform-

# Setup backend
cd backend
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\\Scripts\\activate    # Windows
pip install -r requirements.txt

# Setup frontend (new terminal)
cd frontend
npm install`

const dockerInstallCode = `# Using Docker (if available)
docker-compose up -d

# Or build manually:
docker build -t fasp-backend ./backend
docker build -t fasp-frontend ./frontend

# Run containers:
docker run -p 5000:5000 fasp-backend
docker run -p 5173:5173 fasp-frontend`

const checkEnvCode = `# Check Node.js version (must be 20.19+)
node --version

# Check Python version
python --version

# Check npm version
npm --version

# Verify Git installation
git --version`

const backendStartCode = `# Start backend server
cd backend
source venv/bin/activate
python app.py

# Backend API available at: http://localhost:5000`

const frontendStartCode = `# Start frontend (new terminal)
cd frontend  
npm run dev

# Frontend available at: http://localhost:5173`

// Methods
const scrollToSection = (key) => {
  activeSection.value = key
  const element = document.getElementById(key)
  if (element) {
    element.scrollIntoView({ behavior: 'smooth', block: 'start' })
  }
}

const copyCode = (code) => {
  navigator.clipboard.writeText(code).then(() => {
    proxy.$message.success('Code copied to clipboard!')
  })
}

const getStatusType = (status) => {
  const typeMap = {
    'Significant': 'success',
    'Moderate': 'warning',
    'Not Significant': 'danger'
  }
  return typeMap[status] || 'info'
}

const refreshData = () => {
  proxy.$message.success('Data refreshed')
}

const exportData = () => {
  proxy.$message.success('Data exported successfully')
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
  
  const data = [3280, 2150, 1890, 1650, 1440]
  const xData = ['SE', 'RI', 'A3SS', 'A5SS', 'MXE']
  
  let option = {}
  
  if (type === 'pie') {
    option = {
      title: { text: 'Event Type Distribution', left: 'center' },
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
      title: { text: type === 'bar' ? 'Splicing Events by Type' : 'P-Value Distribution' },
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
  codeTab.value = 'backend'
}

const copyFullExample = () => {
  const example = selectedExample.value
  const fullCode = `${example.code.backend}

${example.code.frontend}

${example.code.config}`
  copyCode(fullCode)
}

const openGitHubIssues = () => {
  window.open('https://github.com/cl9n24abc-star/FASP-Fast-Alternative-Splicing-Platform-/issues', '_blank')
}

const openDocumentation = () => {
  window.open('https://github.com/cl9n24abc-star/FASP-Fast-Alternative-Splicing-Platform-', '_blank')
}

// Lifecycle
onMounted(() => {
  nextTick(() => {
    initChart()
  })
  
  // Listen for scroll events to update active navigation
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
  position: relative;
}

.content {
  padding: 20px;
  padding-right: 240px; /* Make space for sidebar */
  position: relative;
}

.anchor-nav-wrapper {
  position: fixed;
  right: 20px;
  top: 100px;
  z-index: 100;
}

.anchor-nav {
  width: 200px;
  max-height: 400px;
  overflow-y: auto;
}

.anchor-menu {
  border: none;
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

/* Responsive Design */
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

/* Scrollbar Styles */
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