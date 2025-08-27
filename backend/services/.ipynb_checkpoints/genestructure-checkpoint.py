import pandas as pd
import gzip
from collections import defaultdict
import json
import os

def find_vue_public_dir():
    """
    智能查找Vue项目的public目录
    无论从哪里运行脚本，都能正确找到public目录
    """
    current_dir = os.getcwd()  # 这是调用者的工作目录
    print(f"当前工作目录: {current_dir}")
    
    # 情况1: 如果当前目录就有public文件夹（在vue-project根目录运行）
    public_path = os.path.join(current_dir, "public")
    if os.path.exists(public_path):
        print("✓ 找到public目录: ./public")
        return "./public"
    
    # 情况2: 检查父目录是否有public文件夹（在analysis子目录运行）
    parent_dir = os.path.dirname(current_dir)
    parent_public_path = os.path.join(parent_dir, "public")
    if os.path.exists(parent_public_path):
        print("✓ 找到public目录: ../public")
        return "../public"
    
    # 情况3: 向上查找，直到找到package.json和public目录同时存在的目录
    search_dir = current_dir
    for _ in range(5):  # 最多向上查找5级目录
        package_json = os.path.join(search_dir, "package.json")
        public_dir = os.path.join(search_dir, "public")
        
        if os.path.exists(package_json) and os.path.exists(public_dir):
            # 计算相对路径
            rel_path = os.path.relpath(public_dir, current_dir)
            print(f"✓ 找到Vue项目public目录: {rel_path}")
            return rel_path
        
        # 向上一级目录
        parent = os.path.dirname(search_dir)
        if parent == search_dir:  # 已经到根目录了
            break
        search_dir = parent
    
    # 情况4: 如果都没找到，在当前目录创建public文件夹
    print("⚠️ 未找到现有的public目录")
    
    # 检查是否在analysis子目录中
    if current_dir.endswith("analysis") and os.path.exists(parent_dir):
        # 很可能在analysis目录中运行，尝试在父目录创建public
        target_public = "../public"
        print(f"将在父目录创建public目录: {target_public}")
        return target_public
    else:
        # 在当前目录创建
        target_public = "./public"
        print(f"将在当前目录创建public目录: {target_public}")
        return target_public


class RMATSDataProcessor:
    def __init__(self, gtf_file, rmats_dir):
        self.gtf_file = gtf_file
        self.rmats_dir = rmats_dir
        self.gene_structures = {}
        self.rmats_events = {}
        
    def parse_gtf_attribute(self, attribute_string):
        """解析GTF属性字符串"""
        attributes = {}
        for item in attribute_string.strip().split(';'):
            if item.strip():
                key_value = item.strip().split(' ', 1)
                if len(key_value) == 2:
                    key = key_value[0]
                    value = key_value[1].strip('"')
                    attributes[key] = value
        return attributes
    
    def load_gtf_data(self):
        """加载GTF文件并提取基因结构"""
        print("Loading GTF file...")
        
        # 判断是否为压缩文件
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        
        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom, source, feature, start, end, score, strand, frame, attributes = parts
                
                if feature not in ['gene', 'exon']:
                    continue
                
                attr_dict = self.parse_gtf_attribute(attributes)
                
                if feature == 'gene':
                    gene_id = attr_dict.get('gene_id', '')
                    gene_name = attr_dict.get('gene_name', gene_id)
                    
                    self.gene_structures[gene_id] = {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'chromosome': chrom,
                        'strand': strand,
                        'gene_start': int(start),
                        'gene_end': int(end),
                        'exons': []
                    }
                
                elif feature == 'exon':
                    gene_id = attr_dict.get('gene_id', '')
                    if gene_id in self.gene_structures:
                        exon_info = {
                            'start': int(start),
                            'end': int(end),
                            'length': int(end) - int(start) + 1
                        }
                        self.gene_structures[gene_id]['exons'].append(exon_info)
        
        # 对每个基因的外显子按位置排序
        for gene_id in self.gene_structures:
            exons = self.gene_structures[gene_id]['exons']
            exons.sort(key=lambda x: x['start'])
            # 添加外显子编号
            for i, exon in enumerate(exons, 1):
                exon['exon_number'] = i
        
        print(f"Loaded {len(self.gene_structures)} genes from GTF")
    
    def load_rmats_events(self):
        """加载所有rMATS事件文件"""
        event_files = {
            'SE': f"{self.rmats_dir}/SE.MATS.JC.txt",
            'MXE': f"{self.rmats_dir}/MXE.MATS.JC.txt", 
            'A5SS': f"{self.rmats_dir}/A5SS.MATS.JC.txt",
            'A3SS': f"{self.rmats_dir}/A3SS.MATS.JC.txt",
            'RI': f"{self.rmats_dir}/RI.MATS.JC.txt"
        }
        
        all_events = []
        
        for event_type, file_path in event_files.items():
            try:
                print(f"Loading {event_type} events...")
                df = pd.read_csv(file_path, sep='\t')
                
                for _, row in df.iterrows():
                    event = self.parse_rmats_event(row, event_type)
                    if event:
                        all_events.append(event)
                        
            except FileNotFoundError:
                print(f"Warning: {file_path} not found")
            except Exception as e:
                print(f"Error loading {event_type}: {e}")
        
        print(f"Loaded {len(all_events)} rMATS events")
        return all_events
    
    def parse_count_values(self, value):
        """解析包含逗号的计数值"""
        if pd.isna(value):
            return 0
        
        value_str = str(value).strip()
        if ',' in value_str:
            values = [int(v.strip()) for v in value_str.split(',') if v.strip()]
            return values[0]  # 取第一个值
        else:
            try:
                return int(float(value_str))
            except:
                return 0
    
    def parse_rmats_event(self, row, event_type):
        """解析单个rMATS事件"""
        try:
            # 安全转换函数
            def safe_float_convert(value, default=0.0):
                if pd.isna(value):
                    return default
                
                str_val = str(value).strip().upper()
                if str_val in ['NA', 'NAN', '']:
                    return default
                    
                # 处理逗号分隔的值
                if ',' in str_val:
                    str_val = str_val.split(',')[0].strip()
                
                try:
                    return float(str_val)
                except (ValueError, TypeError):
                    return default
            
            # 处理IncLevel值
            inc_level1 = safe_float_convert(row['IncLevel1'])
            inc_level2 = safe_float_convert(row['IncLevel2'])
            inc_level_diff = safe_float_convert(row['IncLevelDifference'])
            
            # 处理基因信息
            gene_id = str(row['GeneID']).strip('"')
            gene_symbol = row['geneSymbol'] if pd.notna(row['geneSymbol']) and str(row['geneSymbol']).upper() != 'NA' else f"Gene_{gene_id}"
            
            # 计算显著性
            p_value = safe_float_convert(row['PValue'], 1.0)
            fdr = safe_float_convert(row['FDR'], 1.0)
            is_significant = p_value < 0.05 and fdr < 0.05
            
            # 处理不同事件类型的坐标字段
            if event_type == 'SE':
                exon_start = int(row['exonStart_0base'])
                exon_end = int(row['exonEnd'])
            elif event_type == 'MXE':
                # MXE事件有两个互斥外显子
                exon_start = int(row['1stExonStart_0base'])
                exon_end = int(row['1stExonEnd'])
                # 检查是否有第二个外显子的坐标
                exon2_start = int(row['2ndExonStart_0base']) if '2ndExonStart_0base' in row and pd.notna(row['2ndExonStart_0base']) else None
                exon2_end = int(row['2ndExonEnd']) if '2ndExonEnd' in row and pd.notna(row['2ndExonEnd']) else None
            elif event_type in ['A5SS', 'A3SS']:
                # A5SS和A3SS事件使用长外显子坐标
                exon_start = int(row['longExonStart_0base'])
                exon_end = int(row['longExonEnd'])
                # 同时保存短外显子信息
                short_start = int(row['shortES'])
                short_end = int(row['shortEE'])
                flanking_start = int(row['flankingES'])
                flanking_end = int(row['flankingEE'])
            elif event_type == 'RI':
                # RI事件使用内含子保留的外显子坐标
                exon_start = int(row['riExonStart_0base'])
                exon_end = int(row['riExonEnd'])
            
            event = {
                'id': int(row['ID']) if pd.notna(row['ID']) else 0,
                'geneId': gene_id,
                'geneSymbol': gene_symbol,
                'eventType': event_type,
                'chromosome': str(row['chr']),
                'strand': str(row['strand']),
                'exonStart': exon_start,
                'exonEnd': exon_end,
                'pValue': p_value,
                'fdr': fdr,
                'incLevelDifference': inc_level_diff,
                'originalInc': inc_level1,
                'compareInc': inc_level2,
                'isSignificant': is_significant,
                'status': 'kept',
                
                # 读段计数
                'ijc_sample1': self.parse_count_values(row['IJC_SAMPLE_1']),
                'sjc_sample1': self.parse_count_values(row['SJC_SAMPLE_1']),
                'ijc_sample2': self.parse_count_values(row['IJC_SAMPLE_2']),
                'sjc_sample2': self.parse_count_values(row['SJC_SAMPLE_2']),
                
                # 长度信息
                'incFormLen': int(row['IncFormLen']) if pd.notna(row['IncFormLen']) else 0,
                'skipFormLen': int(row['SkipFormLen']) if pd.notna(row['SkipFormLen']) else 0,
            }
            
            # 根据事件类型添加特定的坐标信息
            if event_type == 'MXE':
                if exon2_start is not None and exon2_end is not None:
                    event.update({
                        'exon2Start': exon2_start,
                        'exon2End': exon2_end,
                        'hasBothExons': True
                    })
                else:
                    event['hasBothExons'] = False
                    print(f"Warning: MXE event {event['id']} missing second exon coordinates")
            elif event_type in ['A5SS', 'A3SS']:
                event.update({
                    'shortStart': short_start,
                    'shortEnd': short_end,
                    'flankingStart': flanking_start,
                    'flankingEnd': flanking_end,
                })

            # 添加上下游外显子信息（根据事件类型使用不同的列名）
            if event_type in ['SE', 'MXE', 'RI']:
                # SE、MXE和RI使用upstream/downstream
                if 'upstreamES' in row:
                    event.update({
                        'upstreamES': int(row['upstreamES']) if pd.notna(row['upstreamES']) else None,
                        'upstreamEE': int(row['upstreamEE']) if pd.notna(row['upstreamEE']) else None,
                        'downstreamES': int(row['downstreamES']) if pd.notna(row['downstreamES']) else None,
                        'downstreamEE': int(row['downstreamEE']) if pd.notna(row['downstreamEE']) else None,
                    })
            
            return event
            
        except Exception as e:
            print(f"Error parsing {event_type} event: {e}")
            print(f"Available columns: {list(row.keys())}")
            return None
    
    def _get_affected_regions(self, event):
        """根据事件类型获取受影响的区域坐标"""
        regions = []
        
        if event['eventType'] == 'SE':
            regions.append({
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'type': 'target_exon'
            })
        elif event['eventType'] == 'MXE':
            regions.append({
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'type': 'exon1'
            })
            if event.get('hasBothExons', False):
                regions.append({
                    'start': event['exon2Start'],
                    'end': event['exon2End'],
                    'type': 'exon2'
                })
        elif event['eventType'] in ['A5SS', 'A3SS']:
            regions.append({
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'type': 'long_exon'
            })
            if 'shortStart' in event:
                regions.append({
                    'start': event['shortStart'],
                    'end': event['shortEnd'],
                    'type': 'short_exon'
                })
        elif event['eventType'] == 'RI':
            # RI事件中，riExonStart/End指的是保留内含子的范围
            regions.append({
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'type': 'retained_intron_region'
            })
        
        return regions

    def _is_exon_affected(self, exon, affected_regions, event_type):
        """判断外显子是否受可变剪切影响"""
        for region in affected_regions:
            # 检查外显子是否与受影响区域重叠
            if (exon['start'] <= region['end'] and exon['end'] >= region['start']):
                return True
        return False

    def _get_exon_status(self, event_type, is_affected):
        """根据事件类型和是否受影响返回外显子状态"""
        if not is_affected:
            return 'Constitutive'
        
        status_map = {
            'SE': 'Skipped',
            'MXE': 'Mutually Exclusive',
            'A5SS': 'Alternative 5\' SS',
            'A3SS': 'Alternative 3\' SS',
            'RI': 'Retained Intron'
        }
        
        return status_map.get(event_type, 'Alternative')

    def generate_alternative_transcript(self, structure, event_type, event=None):
        """
        生成真正的可变剪切转录本结构，修复MXE和RI逻辑
        """
        alt_structure = []
        visualization_data = {
            'eventType': event_type,
            'affectedElements': [],
            'splicingChanges': []
        }
        
        if event_type == 'SE':
            # Skipped Exon (跳跃外显子)
            skipped_exons = []
            
            for element in structure:
                if element['type'] == 'exon' and element.get('isAffected', False):
                    # 标记为跳过的外显子
                    skipped_element = element.copy()
                    skipped_element['included'] = False
                    skipped_element['skipped'] = True
                    skipped_element['skipReason'] = 'SE_event'
                    alt_structure.append(skipped_element)
                    
                    skipped_exons.append(element['name'])
                    visualization_data['affectedElements'].append({
                        'name': element['name'],
                        'type': 'skipped_exon',
                        'start': element['start'],
                        'end': element['end'],
                        'length': element['length']
                    })
                else:
                    alt_structure.append(element.copy())
            
            visualization_data['splicingChanges'] = [{
                'type': 'exon_skipping',
                'skippedExons': skipped_exons,
                'description': f"Exon(s) {', '.join(skipped_exons)} skipped in alternative transcript"
            }]
        
        elif event_type == 'MXE':
            # Mutually Exclusive Exons (互斥外显子) - 修复版本
            exon_groups = {'A': [], 'B': []}
            
            # 检查是否有完整的MXE数据
            if not event or not event.get('hasBothExons', False):
                print(f"Warning: MXE event lacks complete exon information, creating simplified structure")
                # 如果没有完整数据，创建简化结构
                for element in structure:
                    if element['type'] == 'exon' and element.get('isAffected', False):
                        # 在可变转录本中标记为不同形式
                        alt_element = element.copy()
                        alt_element['name'] += '_alt'
                        alt_element['mutuallyExclusiveGroup'] = 'alternative'
                        alt_structure.append(alt_element)
                    else:
                        alt_structure.append(element.copy())
            else:
                # 有完整MXE数据的情况
                for element in structure:
                    if element['type'] == 'exon' and element.get('isAffected', False):
                        # 检查这个外显子是第一个还是第二个互斥外显子
                        if (element['start'] == event['exonStart'] and element['end'] == event['exonEnd']):
                            # 这是第一个外显子，在可变转录本中替换为第二个
                            alt_exon = {
                                'type': 'exon',
                                'name': f"{element['name']}_B",
                                'width': 8,
                                'isAffected': True,
                                'included': True,
                                'mutuallyExclusiveGroup': 'B',
                                'alternativeExon': True,
                                'length': event['exon2End'] - event['exon2Start'],
                                'start': event['exon2Start'],
                                'end': event['exon2End']
                            }
                            alt_structure.append(alt_exon)
                            
                            exon_groups['A'].append(element['name'])
                            exon_groups['B'].append(alt_exon['name'])
                            
                            visualization_data['affectedElements'].append({
                                'name': element['name'],
                                'type': 'mutually_exclusive_A',
                                'start': element['start'],
                                'end': element['end'],
                                'length': element['length'],
                                'alternativeName': alt_exon['name'],
                                'alternativeStart': alt_exon['start'],
                                'alternativeEnd': alt_exon['end'],
                                'alternativeLength': alt_exon['length']
                            })
                        else:
                            # 其他受影响的外显子保持不变
                            alt_structure.append(element.copy())
                    else:
                        alt_structure.append(element.copy())
            
            visualization_data['splicingChanges'] = [{
                'type': 'mutually_exclusive',
                'exonGroupA': exon_groups['A'],
                'exonGroupB': exon_groups['B'],
                'selectedInOriginal': 'A',
                'selectedInAlternative': 'B',
                'description': f"Mutually exclusive exons: {exon_groups['A']} vs {exon_groups['B']}"
            }]
        
        elif event_type == 'A5SS':
            # Alternative 5' Splice Site (可变5'剪切位点)
            for element in structure:
                if element['type'] == 'exon' and element.get('isAffected', False):
                    # 创建短版本的外显子
                    alt_element = element.copy()
                    
                    if event and event.get('shortStart') and event.get('shortEnd'):
                        # 使用rMATS提供的短外显子坐标
                        alt_element['start'] = event['shortStart']
                        alt_element['end'] = event['shortEnd']
                        alt_element['length'] = event['shortEnd'] - event['shortStart']
                        alt_element['name'] += '*'  # 标记为可变形式
                        alt_element['alternativeForm'] = 'short'
                        alt_element['spliceSiteType'] = '5SS'
                        alt_element['lengthChange'] = element['length'] - alt_element['length']
                        
                        visualization_data['affectedElements'].append({
                            'name': element['name'],
                            'type': 'alternative_5ss',
                            'originalStart': element['start'],
                            'originalEnd': element['end'],
                            'originalLength': element['length'],
                            'alternativeStart': alt_element['start'],
                            'alternativeEnd': alt_element['end'],
                            'alternativeLength': alt_element['length'],
                            'lengthChange': alt_element['lengthChange']
                        })
                    
                    alt_structure.append(alt_element)
                else:
                    alt_structure.append(element.copy())
            
            visualization_data['splicingChanges'] = [{
                'type': 'alternative_5ss',
                'description': "Alternative 5' splice site usage results in shorter exon",
                'effectType': 'shortening'
            }]
        
        elif event_type == 'A3SS':
            # Alternative 3' Splice Site (可变3'剪切位点)
            for element in structure:
                if element['type'] == 'exon' and element.get('isAffected', False):
                    # 创建短版本的外显子（A3SS通常也是缩短）
                    alt_element = element.copy()
                    
                    if event and event.get('shortStart') and event.get('shortEnd'):
                        # 使用rMATS提供的短外显子坐标
                        alt_element['start'] = event['shortStart']
                        alt_element['end'] = event['shortEnd']
                        alt_element['length'] = event['shortEnd'] - event['shortStart']
                        alt_element['name'] += '*'
                        alt_element['alternativeForm'] = 'short'
                        alt_element['spliceSiteType'] = '3SS'
                        alt_element['lengthChange'] = element['length'] - alt_element['length']
                        
                        visualization_data['affectedElements'].append({
                            'name': element['name'],
                            'type': 'alternative_3ss',
                            'originalStart': element['start'],
                            'originalEnd': element['end'],
                            'originalLength': element['length'],
                            'alternativeStart': alt_element['start'],
                            'alternativeEnd': alt_element['end'],
                            'alternativeLength': alt_element['length'],
                            'lengthChange': alt_element['lengthChange']
                        })
                    
                    alt_structure.append(alt_element)
                else:
                    alt_structure.append(element.copy())
            
            visualization_data['splicingChanges'] = [{
                'type': 'alternative_3ss',
                'description': "Alternative 3' splice site usage results in shorter exon",
                'effectType': 'shortening'
            }]
        
        elif event_type == 'RI':
            # Retained Intron (内含子保留) - 修复版本
            retained_introns = []
            
            for element in structure:
                if element['type'] == 'intron':
                    # 检查内含子是否在RI事件范围内
                    intron_start = element.get('start', 0)
                    intron_end = element.get('end', 0)
                    
                    # 检查是否与RI事件重叠（事件的exonStart/End实际是保留内含子的范围）
                    if event and (intron_start >= event.get('exonStart', 0) and 
                                intron_end <= event.get('exonEnd', 0)):
                        
                        # 在可变转录本中，保留的内含子标记为保留状态
                        retained_element = element.copy()
                        retained_element['retained'] = True
                        retained_element['isRetained'] = True
                        retained_element['retentionReason'] = 'RI_event'
                        retained_element['name'] = f"RI-{element['name']}"
                        
                        alt_structure.append(retained_element)
                        retained_introns.append(element['name'])
                        
                        visualization_data['affectedElements'].append({
                            'name': element['name'],
                            'type': 'retained_intron',
                            'start': element['start'],
                            'end': element['end'],
                            'length': element['length'],
                            'retainedAs': retained_element['name']
                        })
                    else:
                        # 正常剪切的内含子
                        normal_intron = element.copy()
                        normal_intron['retained'] = False
                        alt_structure.append(normal_intron)
                else:
                    # 外显子保持不变
                    alt_structure.append(element.copy())
            
            visualization_data['splicingChanges'] = [{
                'type': 'intron_retention',
                'retainedIntrons': retained_introns,
                'description': f"Intron(s) {', '.join(retained_introns)} retained in mature mRNA"
            }]
        
        else:
            # 未知事件类型，返回原始结构
            alt_structure = [element.copy() for element in structure]
            visualization_data['splicingChanges'] = [{
                'type': 'unknown',
                'description': f"Unknown splicing event type: {event_type}"
            }]
        
        return {
            'transcript': alt_structure,
            'visualizationData': visualization_data
        }

    def generate_gene_structure_for_event(self, event):
        """为特定事件生成基因结构可视化数据（增强版本）"""
        gene_id = event['geneId']
        
        # 从GTF中获取基因结构
        if gene_id in self.gene_structures:
            gene_info = self.gene_structures[gene_id]
            
            # 找到与事件相关的外显子
            affected_regions = self._get_affected_regions(event)
            
            # 生成结构数据
            structure_elements = []
            exon_details = []
            
            for i, exon in enumerate(gene_info['exons']):
                is_affected = self._is_exon_affected(exon, affected_regions, event['eventType'])
                
                element = {
                    'type': 'exon',
                    'name': f"E{exon['exon_number']}",
                    'width': 8,
                    'isAffected': is_affected,
                    'length': exon['length'],
                    'start': exon['start'],
                    'end': exon['end'],
                    'eventType': event['eventType'] if is_affected else None,
                    'included': True  # 在原始转录本中都包含
                }
                
                structure_elements.append(element)
                
                # 添加外显子详情
                exon_details.append({
                    'name': f"Exon {exon['exon_number']}",
                    'start': exon['start'],
                    'end': exon['end'],
                    'length': exon['length'],
                    'inclusionLevel': event['originalInc'] if is_affected else 0.9 + (i % 3) * 0.03,
                    'status': self._get_exon_status(event['eventType'], is_affected)
                })
                
                # 添加内含子（除了最后一个外显子）
                if i < len(gene_info['exons']) - 1:
                    next_exon = gene_info['exons'][i + 1]
                    intron_element = {
                        'type': 'intron',
                        'name': f"I{exon['exon_number']}",
                        'width': 3,
                        'length': next_exon['start'] - exon['end'] - 1,
                        'start': exon['end'] + 1,
                        'end': next_exon['start'] - 1,
                        'retained': False  # 默认不保留
                    }
                    structure_elements.append(intron_element)
            
            # 生成原始转录本和可变转录本
            original_transcript = structure_elements.copy()
            
            # 生成可变转录本（包含可视化数据）
            alt_result = self.generate_alternative_transcript(structure_elements, event['eventType'], event)
            alternative_transcript = alt_result['transcript']
            visualization_data = alt_result['visualizationData']
            
            return {
                'geneStart': gene_info['gene_start'],
                'geneEnd': gene_info['gene_end'],
                'strand': gene_info['strand'],
                'splicingEvents': [event['eventType']],
                'structure': structure_elements,
                'originalTranscript': original_transcript,
                'alternativeTranscript': alternative_transcript,
                'exonDetails': exon_details,
                
                # 新增：可视化专用数据
                'visualizationData': visualization_data,
                'splicingAnalysis': {
                    'eventType': event['eventType'],
                    'affectedRegions': affected_regions,
                    'inclusionLevelChange': event.get('incLevelDifference', 0),
                    'significance': {
                        'pValue': event.get('pValue', 1.0),
                        'fdr': event.get('fdr', 1.0),
                        'isSignificant': event.get('isSignificant', False)
                    }
                },
                
                # 事件特定的坐标数据
                'eventCoordinates': self._extract_event_coordinates(event)
            }
        else:
            # 如果GTF中没有找到基因，使用rMATS数据生成简化结构
            return self.generate_minimal_structure_enhanced(event)
        
    def _extract_event_coordinates(self, event):
        """提取事件特定的坐标信息"""
        coordinates = {
            'eventType': event['eventType'],
            'primaryRegion': {
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'length': event['exonEnd'] - event['exonStart']
            }
        }
        
        # 根据事件类型添加特定坐标
        if event['eventType'] == 'MXE':
            if event.get('hasBothExons', False):
                coordinates['secondaryRegion'] = {
                    'start': event['exon2Start'],
                    'end': event['exon2End'],
                    'length': event['exon2End'] - event['exon2Start']
                }
        
        elif event['eventType'] in ['A5SS', 'A3SS']:
            if event.get('shortStart') and event.get('shortEnd'):
                coordinates['alternativeRegion'] = {
                    'start': event['shortStart'],
                    'end': event['shortEnd'],
                    'length': event['shortEnd'] - event['shortStart']
                }
                coordinates['lengthDifference'] = coordinates['primaryRegion']['length'] - coordinates['alternativeRegion']['length']
            
            # 添加flanking外显子信息
            if event.get('flankingStart') and event.get('flankingEnd'):
                coordinates['flankingRegion'] = {
                    'start': event['flankingStart'],
                    'end': event['flankingEnd'],
                    'length': event['flankingEnd'] - event['flankingStart']
                }
        
        # 添加上下游外显子信息（如果存在）
        if event.get('upstreamES') and event.get('upstreamEE'):
            coordinates['upstreamExon'] = {
                'start': event['upstreamES'],
                'end': event['upstreamEE'],
                'length': event['upstreamEE'] - event['upstreamES']
            }
        
        if event.get('downstreamES') and event.get('downstreamEE'):
            coordinates['downstreamExon'] = {
                'start': event['downstreamES'],
                'end': event['downstreamEE'],
                'length': event['downstreamEE'] - event['downstreamES']
            }
        
        return coordinates    
    
    def generate_minimal_structure_enhanced(self, event):
        """增强版最小结构生成（当GTF中找不到基因时），修复RI逻辑"""
        structure_elements = []
        exon_details = []
        
        # 添加上游外显子（如果存在）
        if event.get('upstreamES') and event.get('upstreamEE'):
            structure_elements.append({
                'type': 'exon',
                'name': 'E1',
                'width': 8,
                'isAffected': False,
                'length': event['upstreamEE'] - event['upstreamES'],
                'start': event['upstreamES'],
                'end': event['upstreamEE'],
                'included': True
            })
            
            exon_details.append({
                'name': 'Exon 1',
                'start': event['upstreamES'],
                'end': event['upstreamEE'],
                'length': event['upstreamEE'] - event['upstreamES'],
                'inclusionLevel': 0.95,
                'status': 'Constitutive'
            })
            
            # 添加内含子（针对RI事件特别处理）
            if event['eventType'] == 'RI' and event.get('downstreamES'):
                # RI事件：添加可能被保留的内含子
                intron_start = event['upstreamEE'] + 1
                intron_end = event['downstreamES'] - 1
                
                # 检查这个内含子是否在RI事件范围内
                ri_start = event['exonStart']
                ri_end = event['exonEnd']
                
                # 如果内含子在RI范围内，它可能被保留
                is_retained_intron = (intron_start >= ri_start and intron_end <= ri_end)
                
                structure_elements.append({
                    'type': 'intron',
                    'name': 'I1',
                    'width': 3,
                    'length': intron_end - intron_start + 1,
                    'start': intron_start,
                    'end': intron_end,
                    'retained': False,  # 在原始转录本中不保留
                    'canBeRetained': is_retained_intron  # 标记是否可以被保留
                })
            elif event['exonStart'] > event['upstreamEE']:
                # 非RI事件的普通内含子
                structure_elements.append({
                    'type': 'intron',
                    'name': 'I1',
                    'width': 3,
                    'length': event['exonStart'] - event['upstreamEE'] - 1,
                    'start': event['upstreamEE'] + 1,
                    'end': event['exonStart'] - 1,
                    'retained': False
                })
        
        # 添加目标外显子（受影响的）- 针对不同事件类型处理
        if event['eventType'] == 'RI':
            # RI事件：目标区域实际是保留的内含子，不是外显子
            # 我们需要重新构建结构
            pass  # 在RI的情况下，主要的处理在内含子部分
        else:
            # 其他事件类型：添加目标外显子
            target_exon = {
                'type': 'exon',
                'name': 'E2',
                'width': 8,
                'isAffected': True,
                'length': event['exonEnd'] - event['exonStart'],
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'included': True
            }
            
            # 根据事件类型添加特殊属性
            if event['eventType'] in ['A5SS', 'A3SS'] and event.get('shortStart') and event.get('shortEnd'):
                target_exon['hasAlternativeForm'] = True
                target_exon['alternativeStart'] = event['shortStart']
                target_exon['alternativeEnd'] = event['shortEnd']
                target_exon['alternativeLength'] = event['shortEnd'] - event['shortStart']
            elif event['eventType'] == 'MXE' and event.get('hasBothExons', False):
                target_exon['hasMutuallyExclusiveForm'] = True
                target_exon['alternativeStart'] = event['exon2Start']
                target_exon['alternativeEnd'] = event['exon2End']
                target_exon['alternativeLength'] = event['exon2End'] - event['exon2Start']
            
            structure_elements.append(target_exon)
            
            exon_details.append({
                'name': 'Exon 2',
                'start': event['exonStart'],
                'end': event['exonEnd'],
                'length': event['exonEnd'] - event['exonStart'],
                'inclusionLevel': event['originalInc'],
                'status': self._get_exon_status(event['eventType'], True)
            })
        
        # 添加下游外显子（如果存在）
        if event.get('downstreamES') and event.get('downstreamEE'):
            # 如果还没有添加到目标外显子的内含子，现在添加
            if event['eventType'] != 'RI' and event['downstreamES'] > event['exonEnd']:
                structure_elements.append({
                    'type': 'intron',
                    'name': 'I2',
                    'width': 3,
                    'length': event['downstreamES'] - event['exonEnd'] - 1,
                    'start': event['exonEnd'] + 1,
                    'end': event['downstreamES'] - 1,
                    'retained': False
                })
            elif event['eventType'] == 'RI':
                # RI事件中，可能需要在这里添加内含子
                last_exon = structure_elements[-1] if structure_elements and structure_elements[-1]['type'] == 'exon' else None
                if last_exon and event['downstreamES'] > last_exon['end']:
                    intron_start = last_exon['end'] + 1
                    intron_end = event['downstreamES'] - 1
                    
                    # 检查是否在RI范围内
                    ri_start = event['exonStart']
                    ri_end = event['exonEnd']
                    is_retained_intron = (intron_start >= ri_start and intron_end <= ri_end)
                    
                    structure_elements.append({
                        'type': 'intron',
                        'name': f"I{len([e for e in structure_elements if e['type'] == 'intron']) + 1}",
                        'width': 3,
                        'length': intron_end - intron_start + 1,
                        'start': intron_start,
                        'end': intron_end,
                        'retained': False,
                        'canBeRetained': is_retained_intron
                    })
            
            # 添加下游外显子
            exon_number = len([e for e in structure_elements if e['type'] == 'exon']) + 1
            structure_elements.append({
                'type': 'exon',
                'name': f'E{exon_number}',
                'width': 8,
                'isAffected': False,
                'length': event['downstreamEE'] - event['downstreamES'],
                'start': event['downstreamES'],
                'end': event['downstreamEE'],
                'included': True
            })
            
            exon_details.append({
                'name': f'Exon {exon_number}',
                'start': event['downstreamES'],
                'end': event['downstreamEE'],
                'length': event['downstreamEE'] - event['downstreamES'],
                'inclusionLevel': 0.92,
                'status': 'Constitutive'
            })
        
        # 计算基因范围
        exon_elements = [e for e in structure_elements if e['type'] == 'exon']
        if exon_elements:
            gene_start = min([e['start'] for e in exon_elements])
            gene_end = max([e['end'] for e in exon_elements])
        else:
            gene_start = event['exonStart']
            gene_end = event['exonEnd']
        
        # 生成原始转录本和可变转录本
        original_transcript = structure_elements.copy()
        
        # 生成可变转录本
        alt_result = self.generate_alternative_transcript(structure_elements, event['eventType'], event)
        alternative_transcript = alt_result['transcript']
        visualization_data = alt_result['visualizationData']
        
        return {
            'geneStart': gene_start,
            'geneEnd': gene_end,
            'strand': event['strand'],
            'splicingEvents': [event['eventType']],
            'structure': structure_elements,
            'originalTranscript': original_transcript,
            'alternativeTranscript': alternative_transcript,
            'exonDetails': exon_details,
            'visualizationData': visualization_data,
            'splicingAnalysis': {
                'eventType': event['eventType'],
                'inclusionLevelChange': event.get('incLevelDifference', 0),
                'significance': {
                    'pValue': event.get('pValue', 1.0),
                    'fdr': event.get('fdr', 1.0),
                    'isSignificant': event.get('isSignificant', False)
                }
            },
            'eventCoordinates': self._extract_event_coordinates(event)
        }
    
    def process_all_data(self):
        """处理所有数据并生成前端格式"""
        # 加载GTF数据
        self.load_gtf_data()
        
        # 加载rMATS事件
        events = self.load_rmats_events()
        
        # 为每个事件生成完整的数据
        frontend_data = []
        
        for event in events:
            # 生成基因结构数据
            structure_data = self.generate_gene_structure_for_event(event)
            
            # 组装完整的前端数据
            frontend_item = {
                **event,  # 包含所有rMATS数据
                'structureData': structure_data,
                'abundanceData': {
                    'original': {
                        'tpm': 'N/A',
                        'reads': event['ijc_sample1'] + event['sjc_sample1'],
                        'incLevel': f"{event['originalInc']:.3f}"
                    },
                    'processed': {
                        'tpm': 'N/A', 
                        'reads': event['ijc_sample2'] + event['sjc_sample2'],
                        'incLevel': f"{event['compareInc']:.3f}"
                    }
                },
                'expressionData': {
                    'status': 'upregulated' if event['incLevelDifference'] > 0.1 else 'downregulated' if event['incLevelDifference'] < -0.1 else 'unchanged',
                    'log2FoldChange': 'N/A',
                    'pValue': event['pValue'],
                    'fdr': event['fdr']
                },
                'psiData': {
                    'originalPsi': event['originalInc'],
                    'processedPsi': event['compareInc'],
                    'deltaPsi': event['incLevelDifference']
                },
                'qualityData': {
                    'coverage': max(event['ijc_sample1'], event['ijc_sample2']),
                    'junctionReads': event['ijc_sample1'] + event['ijc_sample2'],
                    'mappingQuality': 'N/A',
                    'confidence': 0.95 if event['isSignificant'] else 0.6
                }
            }
            
            frontend_data.append(frontend_item)
        
        return frontend_data

    @staticmethod
    def compress_data_structure(frontend_data):
        """压缩数据结构，移除不必要的字段"""
        compressed = []
        for item in frontend_data:
            # 只保留前端必需的字段
            compressed_item = {
                # 基础信息
                'id': item['id'],
                'geneId': item['geneId'],
                'geneSymbol': item['geneSymbol'], 
                'eventType': item['eventType'],
                'chromosome': item['chromosome'],
                'strand': item['strand'],
                'exonStart': item['exonStart'],
                'exonEnd': item['exonEnd'],
                
                # 统计信息
                'pValue': item['pValue'],
                'fdr': item['fdr'],
                'incLevelDifference': item['incLevelDifference'],
                'originalInc': item['originalInc'],
                'compareInc': item['compareInc'],
                'isSignificant': item['isSignificant'],
                'status': item['status'],
                
                # 简化的丰度数据
                'reads1': item['ijc_sample1'] + item['sjc_sample1'],
                'reads2': item['ijc_sample2'] + item['sjc_sample2'],
                
                # 简化的质量数据  
                'confidence': 0.95 if item['isSignificant'] else 0.6,
                'structureData': item['structureData']
            }
            
            # 根据事件类型添加必要的坐标信息
            if item['eventType'] == 'MXE' and item.get('hasBothExons', False):
                compressed_item['exon2Start'] = item['exon2Start']
                compressed_item['exon2End'] = item['exon2End']
                compressed_item['hasBothExons'] = True
            elif item['eventType'] in ['A5SS', 'A3SS'] and 'shortStart' in item:
                compressed_item['shortStart'] = item['shortStart']
                compressed_item['shortEnd'] = item['shortEnd']
            
            compressed.append(compressed_item)
        
        return compressed

    def split_data_by_event_type(self, frontend_data, output_dir="."):
        """按事件类型分割数据并压缩"""
        by_type = {}
        
        # 按事件类型分组
        for item in frontend_data:
            event_type = item['eventType']
            if event_type not in by_type:
                by_type[event_type] = []
            by_type[event_type].append(item)
        
        # 为每种事件类型创建文件
        for event_type, events in by_type.items():
            print(f"处理 {event_type} 事件: {len(events)} 条记录")
            
            # 检查数据质量
            if event_type == 'MXE':
                complete_mxe = len([e for e in events if e.get('hasBothExons', False)])
                print(f"  其中 {complete_mxe} 个MXE事件有完整的两个外显子信息")
            elif event_type == 'RI':
                print(f"  RI事件将生成正确的内含子保留结构")
            
            # 压缩数据结构
            compressed_events = self.compress_data_structure(events)
            
            # 保存到单独的文件
            filename = os.path.join(output_dir, f'frontend_data_{event_type}.json')
            with open(filename, 'w') as f:
                json.dump(compressed_events, f, separators=(',', ':'))
            
            # 检查文件大小
            size_mb = os.path.getsize(filename) / (1024 * 1024)
            print(f"  {filename}: {size_mb:.1f} MB")
        
        # 创建一个索引文件
        summary = {
            'total_events': len(frontend_data),
            'event_types': {
                event_type: {
                    'count': len(events),
                    'significant': len([e for e in events if e['isSignificant']]),
                    'filename': f'frontend_data_{event_type}.json'
                }
                for event_type, events in by_type.items()
            }
        }
        
        summary_file = os.path.join(output_dir, 'frontend_data_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return by_type

    def run_complete_analysis(self, output_dir=None, split_by_type=True):
        """运行完整的分析流程，默认输出到Vue项目的public目录"""
        print("开始rMATS数据处理...")
        
        # 如果没有指定输出目录，自动查找
        if output_dir is None:
            output_dir = find_vue_public_dir()
        
        # 确保输出目录存在
        if not os.path.exists(output_dir):
            print(f"创建输出目录: {output_dir}")
            os.makedirs(output_dir)
        
        # 处理所有数据
        frontend_data = self.process_all_data()
        
        if split_by_type:
            # 按事件类型分割数据
            print("按事件类型分割数据...")
            self.split_data_by_event_type(frontend_data, output_dir)
        else:
            # 保存为单个文件
            output_file = os.path.join(output_dir, 'frontend_data_all.json')
            compressed_data = self.compress_data_structure(frontend_data)
            with open(output_file, 'w') as f:
                json.dump(compressed_data, f, separators=(',', ':'))
            
            size_mb = os.path.getsize(output_file) / (1024 * 1024)
            print(f"保存到 {output_file}: {size_mb:.1f} MB")
        
        # 打印统计信息
        print(f"\n处理完成！")
        print(f"总共处理了 {len(frontend_data)} 个可变剪切事件")
        print(f"文件已保存到: {output_dir}")
        
        # 按事件类型统计
        event_counts = {}
        for item in frontend_data:
            event_type = item['eventType']
            event_counts[event_type] = event_counts.get(event_type, 0) + 1
        
        print("\n事件类型统计:")
        for event_type, count in event_counts.items():
            print(f"  {event_type}: {count}")
        
        return frontend_data


# 便捷函数，用于快速调用
def process_rmats_data(gtf_file, rmats_dir, output_dir=None, split_by_type=True):
    """
    便捷函数，用于快速处理rMATS数据，自动输出到Vue项目的public目录
    
    Args:
        gtf_file (str): GTF注释文件路径
        rmats_dir (str): rMATS结果文件所在目录
        output_dir (str): 输出目录，如果为None则自动查找Vue项目的public目录
        split_by_type (bool): 是否按事件类型分割文件，默认True
    
    Returns:
        list: 处理后的前端数据
    """
    # 如果没有指定输出目录，自动查找
    if output_dir is None:
        output_dir = find_vue_public_dir()
        print(f"自动检测到Vue项目public目录: {output_dir}")
    
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        print(f"创建输出目录: {output_dir}")
        os.makedirs(output_dir)
    
    processor = RMATSDataProcessor(gtf_file, rmats_dir)
    return processor.run_complete_analysis(output_dir, split_by_type)


# 测试函数
def test_path_detection():
    """测试路径检测功能"""
    print("=== 路径检测测试 ===")
    print(f"Python脚本位置: {__file__ if '__file__' in globals() else '未知'}")
    print(f"当前工作目录: {os.getcwd()}")
    
    detected_path = find_vue_public_dir()
    full_path = os.path.abspath(detected_path)
    
    print(f"检测到的相对路径: {detected_path}")
    print(f"对应的绝对路径: {full_path}")
    print(f"目录是否存在: {os.path.exists(full_path)}")
    
    return detected_path