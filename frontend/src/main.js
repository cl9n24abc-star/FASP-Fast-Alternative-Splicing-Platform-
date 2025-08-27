import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import { createPinia } from 'pinia'

// 导入Element Plus
import ElementPlus from 'element-plus'
import 'element-plus/dist/index.css'

import * as echarts from 'echarts'

const app = createApp(App)

app.use(ElementPlus)
app.use(createPinia())
app.use(router)

// 同时挂载到两个地方，确保兼容性
app.config.globalProperties.$echarts = echarts
window.echarts = echarts  // 添加这一行

app.mount('#app')