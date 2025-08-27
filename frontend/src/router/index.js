import { createRouter, createWebHistory } from 'vue-router'
import Layout from '../views/Layout/index.vue'
import Home from '../views/Home/index.vue'

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/',
      component: Layout,
      children: [
        // 默认模式路由
        { path: '', name: 'home', component: Home },
        { path: 'Guide', name: 'Guide', component: () => import('../views/Layout/Guide.vue') },
        { path: 'Contact', name: 'Contact', component: () => import('../views/Layout/Contact.vue') },
        // 
        { path: 'analysis', name: 'analysis-home', component: () => import('../views/Layout/Home.vue') },
        { path: 'data-import', name: 'data-import', component: () => import('../views/Layout//DataImport.vue') },

        // Single Analysis 路由
        {
          path: 'single-overview',
          name: 'single-overview',
          component: () => import('../views/Layout/single/DataOverview.vue')
        },
        {
          path: 'single-analysis', 
          name: 'single-analysis', 
          component: () => import('../views/Layout/single/DataAnalyze.vue')
        },
      ]
    }
  ]
})

export default router