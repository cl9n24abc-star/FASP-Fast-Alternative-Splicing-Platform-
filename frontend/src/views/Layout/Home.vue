<template>
  <div class="compare-home">
    <h1>🔍 Comprehensive Analysis Platform</h1>
    <p>Professional alternative splicing analyzing and visualization</p>
    <img src="/images/workflow.png" alt="Analysis Workflow" class="workflow-image" />
    
    <!-- 添加工作流程描述 -->
    <p class="workflow-description">
      This workflow supports both essential and optional analysis components. Required inputs include BAM files, output directories, and rMATS results, while optional features offer enhanced visualization through GTF files and Sashimi plots. All components feed into our comprehensive data analysis engine.
    </p>
  </div>
</template>

<script setup>
</script>

<style scoped>
.compare-home {
  text-align: center;
  padding: 20px;
}

.workflow-image {
  max-width: 100%;
  height: auto;
  margin-top: 20px;
  border-radius: 8px;
  box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
}

.workflow-description {
  margin-top: 20px;
  max-width: 800px;
  margin-left: auto;
  margin-right: auto;
  line-height: 1.6;
  color: #666;
  font-size: 0.95em;
}
</style>