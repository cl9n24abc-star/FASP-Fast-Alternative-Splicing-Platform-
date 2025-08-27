<template>
  <div class="sidebar">
    <!-- Default Mode Menu -->
    <el-menu
      v-if="isDefaultMode"
      :default-active="activeMenu"
      background-color="#D3DCE6"
      text-color="#606266"
      active-text-color="#409eff"
      @select="handleMenuSelect"
    >
      <el-menu-item index="/">
        <el-icon><House /></el-icon>
        <span>Home</span>
      </el-menu-item>
      
      <el-menu-item index="/Guide">
        <el-icon><InfoFilled /></el-icon>
        <span>User Guide</span>
      </el-menu-item>
      
      <el-menu-item index="/contact">
        <el-icon><Message /></el-icon>
        <span>Contact</span>
      </el-menu-item>
    </el-menu>

    <!-- Analysis Mode Menu -->
    <el-menu
      v-if="isAnalysisMode"
      :default-active="activeMenu"
      background-color="#D3DCE6"
      text-color="#606266"
      active-text-color="#409eff"
      @select="handleMenuSelect"
    >
      <el-menu-item index="/analysis">
        <el-icon><DataAnalysis /></el-icon>
        <span>Analysis Home</span>
      </el-menu-item>
      
      <el-menu-item index="/data-import">
        <el-icon><Upload /></el-icon>
        <span>Data Import</span>
      </el-menu-item>
      
      <!-- Single Analysis 子菜单 - 根据数据加载状态控制访问 -->
      <el-sub-menu 
        index="single-analysis"
        :class="{
          'disabled-menu': !store.canAccessAnalysis,
          'enabled-menu': store.canAccessAnalysis
        }"
        @click="handleAnalysisMenuClick"
      >
        <template #title>
          <el-icon><Document /></el-icon>
          <span>Single Analysis</span>
          <!-- 数据加载状态指示器 -->
          <el-icon 
            v-if="store.canAccessAnalysis" 
            class="status-icon success-icon"
            title="Data loaded - Analysis available"
          >
            <CircleCheck />
          </el-icon>
          <el-icon 
            v-else 
            class="status-icon lock-icon"
            title="Please import data first"
          >
            <Lock />
          </el-icon>
        </template>
        
        <el-menu-item 
          index="/single-overview"
          :disabled="!store.canAccessAnalysis"
          :class="{ 'disabled-menu-item': !store.canAccessAnalysis }"
          @click="handleAnalysisItemClick('/single-overview', $event)"
        >
          <el-icon><Grid /></el-icon>
          Data Overview
        </el-menu-item>
        
        <el-menu-item 
          index="/single-analysis"
          :disabled="!store.canAccessAnalysis"
          :class="{ 'disabled-menu-item': !store.canAccessAnalysis }"
          @click="handleAnalysisItemClick('/single-analysis', $event)"
        >
          Single Data Analysis
        </el-menu-item>
      </el-sub-menu>
    </el-menu>
    
    <!-- 数据加载状态指示器 -->
    <div v-if="isAnalysisMode" class="data-status-indicator">
      <div v-if="store.isDataLoaded" class="status-card success-card">
        <div class="status-header">
          <el-icon class="status-icon"><CircleCheck /></el-icon>
          <span class="status-title">Data Loaded</span>
        </div>
        <div class="status-info">
          <p>Analysis menu is now available</p>
          <div class="load-time">
            Loaded: {{ formatLoadTime(store.getLoadTimestamp) }}
          </div>
        </div>
        <el-button 
          type="text" 
          size="small" 
          @click="showResetDataDialog"
          class="reset-button"
        >
          Reset Data
        </el-button>
      </div>
      
      <div v-else class="status-card warning-card">
        <div class="status-header">
          <el-icon class="status-icon"><Lock /></el-icon>
          <span class="status-title">No Data Loaded</span>
        </div>
        <div class="status-info">
          <p>Please import data to access analysis menu</p>
        </div>
        <el-button 
          type="text" 
          size="small" 
          @click="goToDataImport"
          class="import-button"
        >
          Go to Data Import
        </el-button>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed, watch } from 'vue'
import { useAppStore } from '@/stores/app'
import { useRouter, useRoute } from 'vue-router'
import { ElMessage, ElMessageBox } from 'element-plus'
import {
  House,
  DataAnalysis,
  InfoFilled,
  Message,
  Document,
  Upload,
  Grid,
  CircleCheck,
  Lock
} from '@element-plus/icons-vue'

const store = useAppStore()
const router = useRouter()
const route = useRoute()

const activeMenu = ref(route.path)

// 基于路由判断当前模式
const isDefaultMode = computed(() => {
  return route.path === '/' || route.path.includes('/Guide') || route.path.includes('/contact')
})

const isAnalysisMode = computed(() => {
  return route.path.includes('/analysis') || route.path.includes('/data-import') || 
         route.path.includes('/single-')
})

// 监听路由变化，更新激活菜单
watch(() => route.path, (newPath) => {
  activeMenu.value = newPath
  updateCurrentPage(newPath)
})

// 根据路由更新当前页面状态
const updateCurrentPage = (path) => {
  if (path === '/') {
    store.setCurrentPage('home')
  } else if (path === '/analysis') {
    store.setCurrentPage('analysis-home')
  } else if (path === '/data-import') {
    store.setCurrentPage('data-import')
  } else if (path === '/single-overview') {
    store.setCurrentPage('data-analyze')
  } else if (path === '/single-analysis') {
    store.setCurrentPage('data-analyze')
  } else if (path === '/single-charts') {
    store.setCurrentPage('charts')
  } else if (path === '/Guide') {
    store.setCurrentPage('guide')
  } else if (path === '/contact') {
    store.setCurrentPage('contact')
  }
}

// 处理菜单选择
const handleMenuSelect = (key) => {
  router.push(key)
  console.log('Navigate to route:', key)
}

// 处理分析菜单点击
const handleAnalysisMenuClick = () => {
  if (!store.canAccessAnalysis) {
    ElMessage({
      message: 'Please import data first before accessing Single Analysis menu.',
      type: 'warning',
      duration: 3000,
      showClose: true
    })
  }
}

// 处理分析菜单项点击
const handleAnalysisItemClick = (path, event) => {
  if (!store.canAccessAnalysis) {
    event.preventDefault()
    event.stopPropagation()
    
    ElMessage({
      message: 'Data import required. Redirecting to Data Import page...',
      type: 'warning',
      duration: 3000
    })
    
    setTimeout(() => {
      router.push('/data-import')
    }, 1000)
    
    return false
  }
}

// 格式化加载时间
const formatLoadTime = (timestamp) => {
  if (!timestamp) return 'Unknown'
  
  const date = new Date(timestamp)
  const now = new Date()
  const diffMinutes = Math.floor((now - date) / (1000 * 60))
  
  if (diffMinutes < 1) return 'Just now'
  if (diffMinutes < 60) return `${diffMinutes}m ago`
  
  const diffHours = Math.floor(diffMinutes / 60)
  if (diffHours < 24) return `${diffHours}h ago`
  
  return date.toLocaleDateString()
}

// 显示重置数据确认对话框
const showResetDataDialog = async () => {
  try {
    await ElMessageBox.confirm(
      'This will clear the current data loading status, disable Single Analysis menu, and delete all generated files. You will need to import data again.',
      'Reset Data and Clear Files',
      {
        confirmButtonText: 'Reset & Clear',
        cancelButtonText: 'Cancel',
        type: 'warning',
      }
    )
    
    // 显示加载状态
    const loading = ElMessage({
      message: 'Clearing data and generated files...',
      type: 'info',
      duration: 0, // 不自动关闭
      showClose: false
    })
    
    try {
      // 1. 调用后端清理文件API
      const response = await fetch('/api/cleanup/generated-files', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
      })
      
      const result = await response.json()
      
      if (result.success) {
        console.log('Files cleared successfully:', result)
      } else {
        console.warn('File cleanup failed:', result.message)
        ElMessage.warning('Some files may not have been cleared, but data status will still be reset')
      }
    } catch (error) {
      console.error('File cleanup error:', error)
      ElMessage.warning('File cleanup failed, but data status will still be reset')
    }
    
    // 2. 清除store中的所有相关状态
    store.clearDataLoaded()           // 清除数据加载状态
    store.clearSashimiResults()       // 清除Sashimi结果缓存
    store.resetGtfConfig()            // 重置GTF配置
    store.resetSashimiConfig()        // 重置Sashimi配置
    store.clearPathHistory()          // 清除路径历史
    
    // 关闭加载提示
    loading.close()
    
    ElMessage({
      message: 'Data and generated files have been cleared successfully. Please import data again.',
      type: 'success',
      duration: 4000
    })
    
    // 3. 可选：跳转到数据导入页面
    setTimeout(() => {
      router.push('/data-import')
    }, 1500)
    
  } catch {
    // 用户取消，不做任何操作
  }
}

// 前往数据导入页面
const goToDataImport = () => {
  router.push('/data-import')
  ElMessage({
    message: 'Please complete data import to enable analysis features.',
    type: 'info',
    duration: 3000
  })
}

// 调试方法
const debugCurrentState = () => {
  console.log('=== Sidebar 当前状态 ===')
  console.log('currentPage:', store.currentPage)
  console.log('route.path:', route.path)
  console.log('isDefaultMode:', isDefaultMode.value)
  console.log('isAnalysisMode:', isAnalysisMode.value)
  console.log('canAccessAnalysis:', store.canAccessAnalysis)
  console.log('isDataLoaded:', store.isDataLoaded)
  console.log('========================')
}
</script>

<style scoped>
.sidebar {
  height: 100%;
  position: relative;
}

.el-menu {
  border: none;
  height: 100%;
}

/* 禁用菜单样式 */
.disabled-menu {
  opacity: 0.6;
}

.disabled-menu :deep(.el-sub-menu__title) {
  color: #c0c4cc !important;
  cursor: not-allowed;
  background-color: #f5f5f5;
}

.disabled-menu:hover :deep(.el-sub-menu__title) {
  background-color: #f5f5f5 !important;
  color: #c0c4cc !important;
}

/* 启用菜单样式 */
.enabled-menu :deep(.el-sub-menu__title) {
  color: #606266;
  background-color: transparent;
}

.enabled-menu:hover :deep(.el-sub-menu__title) {
  background-color: #ecf5ff;
  color: #409eff;
}

/* 禁用的菜单项样式 */
.disabled-menu-item {
  color: #c0c4cc !important;
  cursor: not-allowed !important;
  background-color: #fafafa !important;
}

.disabled-menu-item:hover {
  background-color: #fafafa !important;
  color: #c0c4cc !important;
}

/* 状态图标样式 */
.status-icon {
  margin-left: auto;
  font-size: 14px;
}

.success-icon {
  color: #67c23a;
}

.lock-icon {
  color: #f56c6c;
}

/* 数据状态指示器 */
.data-status-indicator {
  position: absolute;
  bottom: 20px;
  left: 10px;
  right: 10px;
}

.status-card {
  padding: 16px;
  border-radius: 8px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  font-size: 13px;
}

.success-card {
  background: linear-gradient(135deg, #67c23a 0%, #85ce61 100%);
  color: white;
}

.warning-card {
  background: linear-gradient(135deg, #e6a23c 0%, #f0ad4e 100%);
  color: white;
}

.status-header {
  display: flex;
  align-items: center;
  margin-bottom: 12px;
}

.status-header .status-icon {
  margin-right: 8px;
  margin-left: 0;
  font-size: 16px;
}

.status-title {
  font-weight: bold;
  font-size: 14px;
}

.status-info p {
  margin: 0 0 8px 0;
  font-size: 12px;
  opacity: 0.9;
}

.load-time {
  font-size: 11px;
  opacity: 0.8;
  margin-bottom: 12px;
}

.reset-button,
.import-button {
  color: white !important;
  padding: 6px 12px;
  border: 1px solid rgba(255, 255, 255, 0.3);
  border-radius: 4px;
  font-size: 12px;
  width: 100%;
  transition: all 0.3s ease;
}

.reset-button:hover,
.import-button:hover {
  background-color: rgba(255, 255, 255, 0.1) !important;
  border-color: rgba(255, 255, 255, 0.5);
}

/* 子菜单项样式调整 */
.el-sub-menu .el-menu-item {
  padding-left: 40px !important;
  font-size: 13px;
}

.el-sub-menu .el-menu-item .el-icon {
  margin-right: 8px;
  width: 16px;
  height: 16px;
  display: inline-flex;
  align-items: center;
  justify-content: center;
}

/* 整体侧边栏样式 */
.sidebar {
  height: 100vh;
  position: relative;
  background: #D3DCE6;
  overflow-y: auto;
  padding-bottom: 200px; /* 为底部状态指示器留空间 */
}

/* 菜单项间距优化 */
.el-menu-item {
  margin: 2px 8px;
  border-radius: 6px;
  transition: all 0.3s ease;
}

/* 响应式调整 */
@media (max-width: 768px) {
  .data-status-indicator {
    bottom: 10px;
    left: 5px;
    right: 5px;
  }
  
  .status-card {
    padding: 12px;
  }
  
  .el-sub-menu .el-menu-item {
    padding-left: 35px !important;
  }
}
</style>