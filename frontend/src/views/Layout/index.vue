<template>
  <div class="common-layout">
    <el-container>
      <el-header>
        <Header />
      </el-header>
      <el-container>
        <el-aside width="200px">
          <Sidebar />
        </el-aside>
        <el-main>
          <Main />
        </el-main>
      </el-container>
    </el-container>
  </div>
</template>

<script setup>
import Header from './Header.vue'
import Sidebar from './Sidebar.vue'
import Main from './Main.vue'
</script>

<style scoped>
.common-layout {
  height: 100vh;
}

.el-container {
  height: 100%;
}

.el-header {
  background-color: #B3C0D1;
  padding: 0;  /* 重要：移除默认padding */
}

.el-aside {
  background-color: #D3DCE6;
  padding: 0;  /* 重要：移除默认padding */
}

.el-main {
  background-color: #E9EEF3;
  padding: 0;  /* 重要：移除默认padding */
}
</style>