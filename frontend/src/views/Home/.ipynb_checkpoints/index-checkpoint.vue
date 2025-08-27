<template>
  <div class="home">
    <div class="welcome-section">
      <h1 class="page-title">Welcome to FASP(Fast Alternative Splicing Platform)</h1>
      <p class="page-description">This is a platform for rapid visual analysis of alternative splicing events</p>
    </div>
        
    <!-- Feature Cards -->
    <div class="feature-cards">
      <el-row :gutter="20">
        <el-col :span="8">
          <el-card shadow="hover">
            <template #header>
              <div class="card-header">
                <el-icon size="24" color="#409eff"><Grid /></el-icon>
                <span>1. Fast & Efficient</span>
              </div>
            </template>
            <p>Lightning-fast processing of large-scale RNA-seq datasets with optimized algorithms for real-time visualization</p>
          </el-card>
        </el-col>
                
        <el-col :span="8">
          <el-card shadow="hover">
            <template #header>
              <div class="card-header">
                <el-icon size="24" color="#67c23a"><Monitor /></el-icon>
                <span>User-Friendly Interface</span>
              </div>
            </template>
            <p>Intuitive browser-based interface with simple point-and-click operations, designed for researchers with no coding experience required</p>
          </el-card>
        </el-col>
                
        <el-col :span="8">
          <el-card shadow="hover">
            <template #header>
              <div class="card-header">
                <el-icon size="24" color="#e6a23c"><Setting /></el-icon>
                <span>Comprehensive Analysis </span>
              </div>
            </template>
            <p>Full support for five major splicing event types with customizable filtering capabilities and statistical analysis for publication-quality results</p>
          </el-card>
        </el-col>
      </el-row>
    </div>
        
    <!-- Quick Start -->
    <div class="quick-start">
      <el-alert
        title="Quick Start"
        type="info"
        description="Click the left menu to explore different components and functional modules"
        show-icon
        :closable="false"
      />
    </div>
  </div>
</template>

<script setup>
import { Grid, Monitor, Setting } from '@element-plus/icons-vue'
</script>

<style scoped>
.home {
  max-width: 1200px;
}

.welcome-section {
  text-align: center;
  margin-bottom: 40px;
}

.page-title {
  font-size: 32px;
  font-weight: 600;
  color: #303133;
  margin-bottom: 16px;
}

.page-description {
  font-size: 16px;
  color: #606266;
  margin: 0;
}

.feature-cards {
  margin-bottom: 40px;
}

.card-header {
  display: flex;
  align-items: center;
  gap: 8px;
  font-weight: 600;
}

.quick-start {
  margin-top: 40px;
}
</style>