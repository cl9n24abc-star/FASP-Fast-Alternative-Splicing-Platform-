<template>
  <div class="main-content">
    <el-breadcrumb separator="/">
      <el-breadcrumb-item>Analysis System</el-breadcrumb-item>
      <el-breadcrumb-item>{{ currentPageTitle }}</el-breadcrumb-item>
    </el-breadcrumb>
    
    <div class="content-wrapper">
      <router-view />
    </div>
  </div>
</template>

<script setup>
import { computed } from 'vue'
import { useAppStore } from '@/stores/app'

const appStore = useAppStore()

const currentPageTitle = computed(() => {
  switch (appStore.currentPage) {
    case 'home': return 'Home'
    case 'about': return 'About Us'
    case 'contact': return 'Contact'
    case 'analysis-home': return 'Analysis Home'
    case 'data-import': return 'Data Import'
    case 'data-analyze': return 'Data Analysis'
    case 'charts': return 'Charts Display'
    default: return 'Home'
  }
})
</script>

<style scoped>
.main-content {
  height: 100%;
  padding: 20px;
  color: #333;
}

.el-breadcrumb {
  margin-bottom: 20px;
}

.content-wrapper {
  background-color: #fff;
  border-radius: 4px;
  padding: 20px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.12), 0 0 6px rgba(0, 0, 0, 0.04);
  min-height: calc(100vh - 140px);
}
</style>