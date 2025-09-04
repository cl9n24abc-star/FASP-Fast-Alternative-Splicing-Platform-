<template>
  <div class="header">
    <!-- Left Logo and Title -->
    <div class="header-left">
      <div class="logo">📋</div>
      <h1 class="title">Analysis System</h1>
    </div>
    
    <!-- Center Function Area -->
    <div class="header-center">
      <el-button-group>
        <el-button 
          :type="$route.path === '/' ? 'primary' : ''"
          size="small"
          @click="switchToDefault"
        >
          <el-icon><House /></el-icon>
          Home
        </el-button>
        
        <el-button 
          :type="$route.path.includes('/analysis') ? 'primary' : ''"
          size="small"
          @click="switchToAnalysis"
        >
          <el-icon><DataAnalysis /></el-icon>
          Analysis
        </el-button>
      </el-button-group>
    </div>
    
    <!-- Right User Information -->
    <div class="header-right">
      <span class="user-info">
        <span class="username">FASP</span>
      </span>
    </div>
  </div>
</template>

<script setup>
import { House, DataAnalysis } from '@element-plus/icons-vue'
import { useAppStore } from '@/stores/app'
import { useRouter, useRoute } from 'vue-router'

const appStore = useAppStore()
const router = useRouter()
const route = useRoute()

const switchToDefault = () => {
  console.log('=== Switch to Home ===')
  appStore.setCurrentPage('home')
  router.push('/')
}

const switchToAnalysis = () => {
  console.log('=== Switch to Analysis ===')
  appStore.setCurrentPage('analysis-home')
  router.push('/analysis')
}
</script>

<style scoped>
.header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  height: 60px;
  width: 100%;
  padding: 0 20px;
  color: #333;
}

.header-left {
  display: flex;
  align-items: center;
  flex: 1;
}

.header-center {
  display: flex;
  justify-content: center;
  flex: 1;
}

.header-right {
  display: flex;
  align-items: center;
  gap: 15px;
  flex: 1;
  justify-content: flex-end;
}

.logo {
  font-size: 20px;
  margin-right: 10px;
}

.title {
  font-size: 16px;
  font-weight: 600;
  color: #303133;
  margin: 0;
}

.user-info {
  display: flex;
  align-items: center;
  gap: 8px;
}

.username {
  font-size: 13px;
  color: #606266;
}
</style>