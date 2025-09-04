import json
with open('../frontend/public/bam_analysis_result.json') as f:
    data = json.load(f)
print('=== Group Stats ===')
print('Group1:', data['group1'])
print('Group2:', data['group2'])
print('\n=== Insert Size Check ===')
print('Group1 first 5 values:', data['chart_data']['insert_size_distribution']['group1'][:5])
print('Group2 first 5 values:', data['chart_data']['insert_size_distribution']['group2'][:5])
print('\n=== Chromosome Coverage Check ===')
print('First 5 heatmap entries:', data['chart_data']['chromosome_heatmap'][:5])