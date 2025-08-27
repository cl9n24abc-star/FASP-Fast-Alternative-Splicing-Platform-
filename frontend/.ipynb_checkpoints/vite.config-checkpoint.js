import { fileURLToPath, URL } from 'node:url'
import { defineConfig } from 'vite'  // 确保这行存在
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    vueDevTools(),
  ],
  resolve: {
    alias: {
      '@': fileURLToPath(new URL('./src', import.meta.url))
    },
  },
  server: {
    host: '0.0.0.0',
    port: 5173,
    proxy: {
      '/api': {
        target: 'http://127.0.0.1:5000',  // 使用本地回环地址
        changeOrigin: true,
        secure: false,
        rewrite: (path) => {
          console.log('Proxying:', path);  // 调试日志
          return path;
        }
      }
    }
  }
})