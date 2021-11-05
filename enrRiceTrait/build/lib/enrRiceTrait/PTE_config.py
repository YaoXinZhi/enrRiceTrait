# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 14/07/2020 下午11:48
@Author: xinzhi yao
"""

# todo: convert the config file to bash format.
db_config = {
    'db_path': 'data',
    'ontology_file': f'data/to.wto.ro.funRiceGene.obo',
    'ontology_url': 'https://raw.githubusercontent.com/YaoXinZhi/Plant-Trait-Enrichment/master/data/to-basic.obo',
    # 'PT_anno_file': '../data/Oryzasebase.txt',
    'PT_anno_file': f'data/Total_Association.txt',
    'PT_anno_url': 'https://raw.githubusercontent.com/YaoXinZhi/Plant-Trait-Enrichment/master/data/Oryzasebase.txt',

    'evidence_file': '',

    # data source


}

# background
bg_config = {
    'id_type': ['rap', 'msu', 'gramene'],
}
