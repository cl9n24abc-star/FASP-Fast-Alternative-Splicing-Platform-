import { defineStore } from 'pinia'

export const useAppStore = defineStore('app', {
  state: () => ({
    currentPage: 'home', // 当前页面
    
    // 数据加载状态
    dataLoadStatus: {
      isDataLoaded: false,        // 数据是否已加载
      loadedDataInfo: null,       // 加载的数据信息
      loadTimestamp: null         // 加载时间戳
    },
    
    // 分析配置状态
    analysisConfig: {
      filePaths: {
        single: '',
        singleGtf: '',
        singleOutput: '',
        sample1BamList: '',
        sample2BamList: '',
        // Sashimi Plot 默认文件路径
        sashimiBam1: '',
        sashimiBam2: '',
        sashimiGff: ''
      },
      status: 'idle' // 'idle', 'loading', 'completed', 'error'
    },

    // GTF功能控制状态
    gtfConfig: {
      enabled: false,              // GTF功能是否启用
      filePath: '',                // GTF文件路径
      autoDetectEnable: true       // 是否根据输入自动检测启用状态
    },

    // Sashimi Plot 控制状态
    sashimiConfig: {
      enabled: false,              // Sashimi功能是否启用
      enableSingleSashimi: false,  // 单一分析的 Sashimi 开关（保持向后兼容）
      sashimiResults: [],          // Sashimi 生成结果
      autoDetectEnable: true,      // 是否根据输入自动检测启用状态
      configuration: {}            // Sashimi配置参数
    },

    // 用户输入的路径历史记录
    pathHistory: [] // 存储用户输入的所有路径信息
  }),

  getters: {
    getFilePaths: (state) => state.analysisConfig.filePaths,
    getAnalysisStatus: (state) => state.analysisConfig.status,
    
    // 分析就绪状态检查
    isAnalysisReady: (state) => {
      const { filePaths } = state.analysisConfig
      return filePaths.single.trim() !== '' && 
             filePaths.singleOutput.trim() !== '' &&
             filePaths.sample1BamList.trim() !== '' &&
             filePaths.sample2BamList.trim() !== ''
    },

    // 数据加载状态相关 getters
    isDataLoaded: (state) => state.dataLoadStatus.isDataLoaded,
    canAccessAnalysis: (state) => state.dataLoadStatus.isDataLoaded,
    getLoadedDataInfo: (state) => state.dataLoadStatus.loadedDataInfo,
    getLoadTimestamp: (state) => state.dataLoadStatus.loadTimestamp,

    // GTF相关getters
    getGtfConfig: (state) => state.gtfConfig,
    isGtfEnabled: (state) => state.gtfConfig.enabled,
    getGtfFilePath: (state) => state.gtfConfig.filePath,

    // Sashimi 相关 getters
    getSashimiConfig: (state) => state.sashimiConfig,
    isSashimiEnabled: (state) => state.sashimiConfig.enabled,
    isSingleSashimiEnabled: (state) => state.sashimiConfig.enableSingleSashimi, // 保持兼容性
    getSashimiResults: (state) => state.sashimiConfig.sashimiResults,
    getSashimiConfiguration: (state) => state.sashimiConfig.configuration,

    // 模块启用状态汇总
    getModuleStatus: (state) => ({
      gtf: state.gtfConfig.enabled,
      sashimi: state.sashimiConfig.enabled,
      rmats: true, // rMATS是必选的
      bam: true    // BAM分析是必选的
    }),

    // 路径历史相关getters
    getPathHistory: (state) => state.pathHistory,
    getRecentPathHistory: (state) => (count = 10) => {
      return state.pathHistory.slice(0, count)
    }
  },

  actions: {
    setCurrentPage(page) {
      console.log('setCurrentPage被调用，参数:', page)
      this.currentPage = page
      console.log('切换到页面:', page)
    },
    
    // 设置数据加载成功
    setDataLoaded(dataInfo) {
      console.log('setDataLoaded被调用，参数:', dataInfo)
      this.dataLoadStatus = {
        isDataLoaded: true,
        loadedDataInfo: dataInfo,
        loadTimestamp: new Date().toISOString()
      }
      console.log('数据加载状态已设置为成功:', this.dataLoadStatus)
    },
    
    // 清除数据加载状态
    clearDataLoaded() {
      console.log('clearDataLoaded被调用')
      this.dataLoadStatus = {
        isDataLoaded: false,
        loadedDataInfo: null,
        loadTimestamp: null
      }
      console.log('数据加载状态已清除')
    },
    
    // 重置数据加载状态（用于重新导入数据）
    resetDataLoadStatus() {
      console.log('resetDataLoadStatus被调用')
      this.clearDataLoaded()
      console.log('数据加载状态已重置，Single Analysis 菜单将被禁用')
    },

    // GTF相关actions
    setGtfEnabled(enabled) {
      console.log('setGtfEnabled被调用，参数:', enabled)
      this.gtfConfig.enabled = enabled
      console.log('GTF功能设置为:', enabled ? '启用' : '禁用')
    },

    setGtfFilePath(filePath) {
      console.log('setGtfFilePath被调用，参数:', filePath)
      this.gtfConfig.filePath = filePath
      
      // 如果启用了自动检测，根据路径是否为空自动设置启用状态
      if (this.gtfConfig.autoDetectEnable) {
        this.gtfConfig.enabled = !!(filePath && filePath.trim())
        console.log('GTF自动检测结果:', this.gtfConfig.enabled ? '启用' : '禁用')
      }
    },

    setGtfAutoDetect(autoDetect) {
      console.log('setGtfAutoDetect被调用，参数:', autoDetect)
      this.gtfConfig.autoDetectEnable = autoDetect
      
      // 如果启用自动检测，立即根据当前路径设置状态
      if (autoDetect) {
        this.gtfConfig.enabled = !!(this.gtfConfig.filePath && this.gtfConfig.filePath.trim())
      }
    },

    resetGtfConfig() {
      console.log('resetGtfConfig被调用')
      this.gtfConfig = {
        enabled: false,
        filePath: '',
        autoDetectEnable: true
      }
      console.log('GTF配置已重置')
    },

    // Sashimi相关actions
    setSashimiEnabled(enabled) {
      console.log('setSashimiEnabled被调用，参数:', enabled)
      this.sashimiConfig.enabled = enabled
      console.log('Sashimi功能设置为:', enabled ? '启用' : '禁用')
    },

    setSashimiConfiguration(config) {
      console.log('setSashimiConfiguration被调用，参数:', config)
      this.sashimiConfig.configuration = config
      
      // 如果启用了自动检测，根据配置是否为空自动设置启用状态
      if (this.sashimiConfig.autoDetectEnable) {
        this.sashimiConfig.enabled = !!(config && Object.keys(config).length > 0)
        console.log('Sashimi自动检测结果:', this.sashimiConfig.enabled ? '启用' : '禁用')
      }
    },

    setSashimiAutoDetect(autoDetect) {
      console.log('setSashimiAutoDetect被调用，参数:', autoDetect)
      this.sashimiConfig.autoDetectEnable = autoDetect
      
      // 如果启用自动检测，立即根据当前配置设置状态
      if (autoDetect) {
        this.sashimiConfig.enabled = !!(this.sashimiConfig.configuration && 
                                        Object.keys(this.sashimiConfig.configuration).length > 0)
      }
    },
    
    // 保持向后兼容性的方法
    setSingleSashimi(enabled) {
      console.log('setSingleSashimi被调用，参数:', enabled)
      this.sashimiConfig.enableSingleSashimi = enabled
      this.sashimiConfig.enabled = enabled // 同时设置新的enabled状态
      console.log('Single Sashimi设置为:', enabled)
    },
    
    addSashimiResult(result) {
      console.log('addSashimiResult被调用，参数:', result)
      this.sashimiConfig.sashimiResults.push(result)
      console.log('Sashimi结果已添加，当前结果数量:', this.sashimiConfig.sashimiResults.length)
    },
    
    clearSashimiResults() {
      console.log('clearSashimiResults被调用')
      this.sashimiConfig.sashimiResults = []
      console.log('Sashimi结果已清空')
    },
    
    // 重置 Sashimi 配置
    resetSashimiConfig() {
      console.log('resetSashimiConfig被调用')
      this.sashimiConfig = {
        enabled: false,
        enableSingleSashimi: false,
        sashimiResults: [],
        autoDetectEnable: true,
        configuration: {}
      }
      console.log('Sashimi配置已重置')
    },

    // 更新文件路径，同时处理GTF自动检测
    updateFilePaths(paths) {
      console.log('updateFilePaths被调用，参数:', paths)
      this.analysisConfig.filePaths = { ...this.analysisConfig.filePaths, ...paths }
      
      // 如果更新了GTF路径，触发自动检测
      if (paths.singleGtf !== undefined) {
        this.setGtfFilePath(paths.singleGtf)
      }
      
      console.log('文件路径已更新:', this.analysisConfig.filePaths)
    },
    
    // 清空文件路径的方法，但保留Sashimi默认路径
    clearFilePaths() {
      console.log('clearFilePaths被调用')
      this.analysisConfig.filePaths = {
        single: '',
        singleGtf: '',
        singleOutput: '',
        sample1BamList: '',
        sample2BamList: '',
        // 保留Sashimi默认路径
        sashimiBam1: '/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam',
        sashimiBam2: '/mnt/disk1/rna_seq_analysis/step0/SRR12125135_Aligned.sortedByCoord.out.bam',
        sashimiGff: '/home/cl9n24abc/project/ASevents/gencode.v48.annotation/gencode.v48.annotation.gff3'
      }
      
      // 清除路径时也清除数据加载状态和模块状态
      this.clearDataLoaded()
      this.resetGtfConfig()
      
      console.log('所有文件路径已清除，但保留了Sashimi默认路径，数据加载状态已重置')
    },
    
    setAnalysisStatus(status) {
      console.log('setAnalysisStatus被调用，参数:', status)
      this.analysisConfig.status = status
      console.log('分析状态设置为:', status)
    },
    
    resetAnalysisConfig() {
      console.log('resetAnalysisConfig被调用')
      this.analysisConfig = {
        filePaths: {
          single: '',
          singleGtf: '',
          singleOutput: '',
          sample1BamList: '',
          sample2BamList: '',
          // 保留Sashimi默认路径
          sashimiBam1: '/mnt/disk1/rna_seq_analysis/step0/SRR12125133_Aligned.sortedByCoord.out.bam',
          sashimiBam2: '/mnt/disk1/rna_seq_analysis/step0/SRR12125135_Aligned.sortedByCoord.out.bam',
          sashimiGff: '/home/cl9n24abc/project/ASevents/gencode.v48.annotation/gencode.v48.annotation.gff3'
        },
        status: 'idle'
      }
      
      // 重置配置时也重置数据加载状态和模块状态
      this.clearDataLoaded()
      this.resetGtfConfig()
      this.resetSashimiConfig()
      
      console.log('分析配置已重置，但保留了Sashimi默认路径，数据加载状态已重置')
    },
    
    startAnalysis() {
      console.log('startAnalysis被调用')
      if (this.isAnalysisReady) {
        this.analysisConfig.status = 'loading'
        console.log('开始单一分析')
        console.log('模块状态:', this.getModuleStatus)
        // 这里可以触发实际的分析逻辑
      } else {
        console.warn('分析配置未准备就绪')
        console.log('当前配置:', {
          filePaths: this.analysisConfig.filePaths,
          isReady: this.isAnalysisReady
        })
      }
    },

    // 保存用户输入的路径信息
    savePathHistory() {
      const historyItem = {
        id: Date.now(),
        timestamp: new Date().toISOString(),
        paths: { ...this.analysisConfig.filePaths },
        gtfConfig: { ...this.gtfConfig },
        sashimiConfig: { ...this.sashimiConfig }
      }
      
      // 添加到历史记录最前面
      this.pathHistory.unshift(historyItem)
      
      // 限制历史记录数量，只保留最近20条
      if (this.pathHistory.length > 20) {
        this.pathHistory = this.pathHistory.slice(0, 20)
      }
      
      console.log('用户输入的路径信息已保存:', historyItem)
    },
    
    // 清空路径历史记录
    clearPathHistory() {
      this.pathHistory = []
      console.log('路径历史记录已清空')
    },
    
    // 完全重置应用状态（包括数据加载状态）
    resetAllState() {
      console.log('resetAllState被调用')
      this.clearFilePaths()
      this.resetSashimiConfig()
      this.resetGtfConfig()
      this.clearPathHistory()
      this.clearDataLoaded()
      this.analysisConfig.status = 'idle'
      console.log('所有应用状态已重置')
    },
    
    // 调试方法，用于查看当前状态
    debugState() {
      console.log('=== App Store 当前状态 ===')
      console.log('currentPage:', this.currentPage)
      console.log('dataLoadStatus:', this.dataLoadStatus)
      console.log('analysisConfig:', this.analysisConfig)
      console.log('gtfConfig:', this.gtfConfig)
      console.log('sashimiConfig:', this.sashimiConfig)
      console.log('pathHistory:', this.pathHistory)
      console.log('isAnalysisReady:', this.isAnalysisReady)
      console.log('isDataLoaded:', this.isDataLoaded)
      console.log('canAccessAnalysis:', this.canAccessAnalysis)
      console.log('isGtfEnabled:', this.isGtfEnabled)
      console.log('isSashimiEnabled:', this.isSashimiEnabled)
      console.log('moduleStatus:', this.getModuleStatus)
      console.log('========================')
    }
  }
})