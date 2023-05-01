# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 09/07/2020 下午4:15
@Author: xinzhi yao
"""
# 1. check if the package is complete.
# python3 setup.py check
# 2. generate the package file.
# python3 setup.py sdist bdist_wheel
# only test
# python3 -m twine upload --repository testpypi dist/*
# upload to your pypi (obselate)
# python3 setup.py register sdist upload -r http://pypi.org
# 3. upload to your pypi (please delete the previously .tar.gz and .whl files.)
# twine upload dist/*


import os
import re
import math
#import matplotlib
import numpy as np
import pandas as pd
#from plotnine import *
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from collections import defaultdict
from .PTE_config import db_config


class enrichment_result:
    def __init__(self):

        self.query_gene_set = set()
        self.miss_id_set = set()

        #m: Total number of target objects / Number of genes in query gene set (Need to subtract the number
        self.terms_count = defaultdict(int)

        # Statistics of each concept
        self.terms_p_value = defaultdict(float)
        self.terms_p_adjust = defaultdict(float)
        self.terms_q_value = defaultdict(float)

        # Bonferroni or Benjamini & Hochberg correlation
        self.terms_p_bh = defaultdict(float)
        self.terms_p_bf = defaultdict(float)

        self.Gene_Ratio = defaultdict(float)
        self.BG_Ratio = defaultdict(float)


class Association:
    def __init__(self, gene_id: str, concept_id: str, source: str, evidence=''):
        self.gene_id = gene_id
        self.concept_id = concept_id
        self.source = source
        self.evidence = evidence

class Background:
    def __init__(self,):

        self.initialized = False

        self.gene_set = set()
        self.trait_set = set()

        self.gene_to_trait = defaultdict(set)
        self.trait_to_gene = defaultdict(set)

        self.pair_to_source = defaultdict(set)
        self.pair_to_evidence = defaultdict(set)

        self.N = 0
        self.m_dict = defaultdict(int)
        self.bg_gene_count = 0
        self.bg = defaultdict(int)

    def update_association(self, gene_id_set: set, concept_id_set: set, source: str,
                        evidence=''):
        for gene_id in gene_id_set:
            for concept in concept_id_set:
                # self.association_list.append(Association(gene_id, concept, source, evidence))

                self.gene_set.add(gene_id)
                self.trait_set.add(concept)

                self.gene_to_trait[gene_id].add(concept)
                self.trait_to_gene[concept].add(gene_id)

                self.pair_to_source[(gene_id, concept)].add(source)
                self.pair_to_evidence[(gene_id, concept)].add(evidence)

class trait_concept:
    def __init__(self, trait_id: str, trait_name: str, definition: str):
        self.id = trait_id
        self.trait_name = trait_name
        self.definition = definition


class Gene_Ontology_enrichment:

    def __init__(self, source_set: set):

        self.gene_term_file = db_config['PT_anno_file']
        self.obo_file = db_config['ontology_file']

        self.background_source = {'Oryzabase', 'TAS', 'funRiceGene', 'SemanticComputing', 'ExactMatching'}
        wrong_source = source_set - self.background_source
        if wrong_source:
            raise TypeError(f'Wrong type of source: {wrong_source}, '
                            f'the type must be "Oryzabase", "TAS", "funRiceGene"'
                            f', "SemanticComputing" or "ExactMatching".')
        self.source_set = source_set

        self.id_to_trait = defaultdict(trait_concept)

        self.rap_background = Background()
        self.msu_background = Background()
        self.gramene_background = Background()

        self.p_threshold = 0

        self.enrich_result = enrichment_result()

        self.enrich_dataframe = ''

        self.load_ontology()
        self.load_background_data()

    # todo: ontology class
    #       parent, child or other.
    def load_ontology(self):
        """
        load ontology file.
        """
        _id = ''
        # fixme: laod other ontology information, e.g. definition.
        # todo: delete the concept the '/' in term id.
        with open(self.obo_file, encoding='utf-8') as f:
            for line in f:
                l = line.strip().split('\t')
                if line.startswith('id:'):
                    _id = l[1]
                elif line.startswith('name:'):
                    if len(l) < 2:
                        continue
                    trait = l[1]
                elif line.startswith('def:'):
                    definition = re.findall(r'\"(.*?)\"', line)[0]
                elif line.startswith('xref:'):
                    xref_id = l[ 1 ]
                    self.id_to_trait[xref_id] = trait_concept(xref_id, trait, definition)
                elif l == [''] and _id and '/' not in _id:
                    self.id_to_trait[_id] = trait_concept(_id, trait, definition)

    # Distinguish different ontology
    # todo: add the evidence to the last col of background data.
    # add source info in backgroud data.
    def load_background_data(self):
        """
        load backgroud data.
        """
        print(f'Loading background data, data source: {self.source_set}.')
        with open(self.gene_term_file) as f:
            f.readline()
            for line in f:
                l = line.strip().split('\t')

                rap_set = set(re.findall(r'Os\d{2}g\d{7}', line))
                msu_set = set(re.findall(r'LOC_Os\d{2}g\d{5}', line))
                gramene_set = set(re.findall(r'GR:\d{7}', line))

                concept_set = set(l[3].split(','))
                source = l[4]

                if source not in self.source_set:
                    continue

                self.rap_background.update_association(rap_set, concept_set, source)
                self.msu_background.update_association(msu_set, concept_set, source)
                self.gramene_background.update_association(gramene_set, concept_set, source)


    @staticmethod
    def Hyper_Geometric_Test(N: int, m: int, k: int, x: int):
        """
        :param N: Total number of objects / Total(All) number of rice genes in background data
        :param k: Total number of objects grabbed / Number of genes associated with each phenotype
        :param m: Total number of target objects / Number of genes in query gene set (Need to subtract the number
                                                    of genes not in the background data)
        :param x: Number of target objects grabbed / The number of genes in the phenotype-related
                                                    gene set and the query gene set at the same time
        :return p-value
        """
        rv = hypergeom(N, k, m)
        x_arg = np.arange(x, k + 1)
        pmf_dogs = rv.pmf(x_arg)
        p = sum(pmf_dogs)
        return p

    def p_bh_and_bf_compute(self):
        count_terms = len(self.enrich_result.terms_p_value.keys())

        sorted_terms = sorted(self.enrich_result.terms_p_value, key=lambda x: self.enrich_result.terms_p_value[x])
        #  k need add 1
        for k, term in enumerate(sorted_terms):
            # updated p-adjust computing for BH
            self.enrich_result.terms_p_bh[term] = min((self.enrich_result.terms_p_value[term] * count_terms) / (k+1), 1)

        for term in self.enrich_result.terms_p_value.keys():
            # updated p-adjust computing for Bonferroui
            self.enrich_result.terms_p_bf[term] = min(self.enrich_result.terms_p_value[term] * count_terms, 1)


    def p_adjust_compute(self, method='BH'):
        count_terms = len(self.enrich_result.terms_p_value.keys())

        if method not in {'BH', 'Bonferroui'}:
            raise TypeError('Method must be "BH" or "Bonferroui"')

        if method == 'BH':
            sorted_terms = sorted(self.enrich_result.terms_p_value, key=lambda x: self.enrich_result.terms_p_value[x])
            #  k need add 1
            for k, term in enumerate(sorted_terms):
                # updated p-adjust computing for BH
                # self.enrich_result.terms_p_adjust[term] = (self.p_threshold * (k+1)) / count_terms
                self.enrich_result.terms_p_adjust[term] = min((self.enrich_result.terms_p_value[term] * count_terms) / (k+1), 1)

        elif method == 'Bonferroui':
            for term in self.enrich_result.terms_p_value.keys():
                # updated p-adjust computing for Bonferroui
                # self.enrich_result.terms_p_adjust[term] = self.p_threshold / count_terms
                self.enrich_result.terms_p_adjust[term] = min(self.enrich_result.terms_p_value[term] / count_terms, 1)
        else:
            raise TypeError("The p-value correction method must be 'BH' or 'Bonferroui'")

    def q_value_compute(self):

        count_terms = len(self.enrich_result.terms_p_value)
        sorted_terms = sorted(self.enrich_result.terms_p_value, key=lambda x: self.enrich_result.terms_p_value[x])

        for k, trait_id in enumerate(sorted_terms):
            self.enrich_result.terms_q_value[trait_id] = (self.enrich_result.terms_p_value[trait_id] * count_terms) / (k+1)

    def ontology_enricement(self, query_gene_set: list or set, id_type: str,
                            save_path='../result',  prefix='EnrRiceTrait',
                            save_result=False, p_threshold=0.05, p_adjust_method='Bonferroui'):
        """
        :param query_gene_set: querying gene set
        :param id_type: id type of querying gene set
        :param p_threshold: threshold for p-value
        :return:
        """

        self.p_threshold = p_threshold
        self.enrich_result = enrichment_result()

        if id_type.lower() == 'rap':
            bg_data = self.rap_background
        elif id_type.lower() == 'msu':
            bg_data = self.msu_background
        elif id_type.lower() == 'gramene':
            bg_data = self.gramene_background
        else:
            print('Gene id type must be "RAPDB", "MSU" or "Gramene".')
            raise TypeError

        N = len(bg_data.gene_set)

        # query_gene_count
        k = 0
        matched_gene_set = set()
        for _id in query_gene_set:
            if _id in bg_data.gene_set:
                k += 1
                matched_gene_set.add(_id)
                # number of trait related this gene.
                for trait_id in bg_data.gene_to_trait[_id]:
                    self.enrich_result.terms_count[trait_id] += 1
            else:
                self.enrich_result.miss_id_set.add(_id)

        print(f'{len(matched_gene_set):,}/{len(query_gene_set):,} are found in the background data,'
              f' include: {matched_gene_set}.')

        print(f'{len(self.enrich_result.miss_id_set):,}/{len(query_gene_set):,} genes are missed in the background data,'
              f' include: {self.enrich_result.miss_id_set}.')

        if len(matched_gene_set) == 0:
            raise ValueError('No gene match in enrRiceTrait background data.')

        # enrichment
        for trait_id in self.enrich_result.terms_count.keys():
            trait_related_gene = bg_data.trait_to_gene[trait_id]
            m = len(bg_data.trait_to_gene[trait_id])
            x = len(trait_related_gene & query_gene_set)

            p_value = self.Hyper_Geometric_Test(N, m, k, x)
            if p_value < self.p_threshold:
                self.enrich_result.terms_p_value[trait_id] = p_value
                self.enrich_result.Gene_Ratio[trait_id] =  x / k
                self.enrich_result.BG_Ratio[trait_id] = k / N

        # method need to give users options
        self.p_adjust_compute(p_adjust_method)
        self.q_value_compute()
        self.p_bh_and_bf_compute()

        sorted_result = sorted(self.enrich_result.terms_p_value.keys(),
                               key=lambda x: self.enrich_result.terms_p_value[x])

        # change the print information of this package
        print(f'{"Trait id":<10}\t{"Trait name":<50}\t{"p-value":<8}\t{"Bonferroni adjusted p-value":<20}'
              f'\t{"BH adjusted p-value":<8}')
        for trait_id in sorted_result:
            print(f'{trait_id:<10}\t{self.id_to_trait[ trait_id ].trait_name:<50}\t'
                  f'{self.enrich_result.terms_p_value[ trait_id ]:<.2e}\t'
                  f'{self.enrich_result.terms_p_bf[ trait_id ]:<.2e}\t'
                  f'{self.enrich_result.terms_p_bh[ trait_id ]:<.2e}')

        # print(f'{"Trait id":<10}\t{"Trait name":<50}\t{"p-value":<8}\t{"p-adjust":<8}'
        #       f'\t{"q-value":<8}')
        # for trait_id in sorted_result:
        #     print(f'{trait_id:<10}\t{self.id_to_trait[trait_id].trait_name:<50}\t'
        #           f'{self.enrich_result.terms_p_value[trait_id]:<.2e}\t'
        #           f'{self.enrich_result.terms_p_adjust[trait_id]:<.2e}\t'
        #           f'{self.enrich_result.terms_q_value[trait_id]:<.2e}')

        if save_result:
            if not save_path:
                raise TypeError('You have to provide a path to save the result.')
            self.save_result(self.enrich_result, save_path, prefix)

        self.dataframe_init()

        return self.enrich_result

    def dataframe_init(self):
        sorted_trait_id = sorted(self.enrich_result.terms_p_value,
                                 key=lambda x: self.enrich_result.terms_p_value[ x ])
        trait_name = [self.id_to_trait[trait_id].trait_name for trait_id in sorted_trait_id]
        trait_definition = [self.id_to_trait[trait_id].definition for trait_id in sorted_trait_id]
        trait_GeneRatio = [self.enrich_result.Gene_Ratio[trait_id] for trait_id in sorted_trait_id]
        trait_BgRatio = [self.enrich_result.BG_Ratio[trait_id] for trait_id in sorted_trait_id]
        trait_p = [self.enrich_result.terms_p_value[trait_id] for trait_id in sorted_trait_id]
        trait_p_bf = [self.enrich_result.terms_p_bf[trait_id] for trait_id in sorted_trait_id]
        trait_p_bh = [self.enrich_result.terms_p_bh[trait_id] for trait_id in sorted_trait_id]

        trait_p_adjust = [self.enrich_result.terms_p_adjust[trait_id] for trait_id in sorted_trait_id]
        trait_q = [self.enrich_result.terms_q_value[trait_id] for trait_id in sorted_trait_id]

        trait_count = [self.enrich_result.terms_count[trait_id] for trait_id in sorted_trait_id]

        self.enrich_dataframe = pd.DataFrame({'ID': sorted_trait_id, 'Name': trait_name,
                                         'Description': trait_definition,
                                         'GeneRatio': trait_GeneRatio,
                                         'BgRatio': trait_BgRatio,
                                         'p-value': trait_p,
                                         'p-adjust': trait_p_adjust,
                                         'q-value': trait_q,
                                         'p-bh': trait_p_bh,
                                         'p-bf': trait_p_bf,
                                         'Count': trait_count
                                         })

        self.enrich_dataframe[ 'ID_Name' ] = self.enrich_dataframe[ 'ID' ] + ' ' + self.enrich_dataframe[ 'Name' ]

        self.enrich_dataframe[ 'p_value_str' ] = 'p-Value ' + self.enrich_dataframe[ 'p-value' ].map('{:.3e}'.format)

        self.enrich_dataframe[ 'neg_log_p_value' ] = -self.enrich_dataframe[ 'p-value' ].map(np.log)

        rows, _ = self.enrich_dataframe.shape
        highest_neg_log_p_value_label = list(self.enrich_dataframe['ID'][0:5])
        highest_neg_log_p_value_label.extend([''] * (rows - 5))

        self.enrich_dataframe[ 'highest_neg_log_p_value_label' ] = highest_neg_log_p_value_label


    def save_result(self, enrich_result, save_path: str, prefix='EnrRiceTrait'):

        if not os.path.exists(save_path):
            os.mkdir(save_path)

        sorted_result = sorted(self.enrich_result.terms_p_value.keys(),
                               key=lambda x: self.enrich_result.terms_p_value[ x ])

        save_file = f'{save_path}/{prefix}.report.txt'
        with open(save_file, 'w') as wf:
            # wf.write('ID\tName\tDescription\tGeneRatio\tBgRatio\tp-value\t'
            #          'p-adjust\tq-value\tCount\n')
            wf.write('ID\tName\tDescription\tGeneRatio\tBgRatio\tp-value\t'
                     'Bonferroni adjusted p-adjust\tBH adjusted p-value\tCount\n')
            for term_id in sorted_result:
                trait_name = self.id_to_trait[term_id].trait_name
                definition = self.id_to_trait[term_id].definition \
                                if self.id_to_trait[term_id].definition else "Noun"
                gene_ratio = enrich_result.Gene_Ratio[term_id]
                bg_raito = enrich_result.BG_Ratio[term_id]
                p_value = enrich_result.terms_p_value[term_id]
                # p_adjust_value = enrich_result.terms_p_adjust[term_id]
                # q_value = enrich_result.terms_q_value[term_id]

                p_bf = enrich_result.terms_p_bf[term_id]
                p_bh = enrich_result.terms_p_bh[term_id]

                # add term count for query gene set.
                trait_count = enrich_result.terms_count[term_id]
                wf.write(f'{term_id}\t{trait_name}\t{definition}\t{gene_ratio}\t'
                         f'{bg_raito}\t{p_value}\t{p_bf}\t{p_bh}\t{trait_count}\n')
                # wf.write(f'{term_id}\t{trait_name}\t{definition}\t{gene_ratio}\t'
                #          f'{bg_raito}\t{p_value}\t{p_adjust_value}\t{q_value}\t{trait_count}\n')
        print(f'{save_file} save done.')

    # latest bar plot 2023-04-23
    def bar_stat(self, save_path='.', prefix='enrRiceTrait.bar', save_plot=False):
        # plot data from large p-value to small p-value
        sorted_result = sorted(self.enrich_result.terms_p_value.keys(),
                               key=lambda x: self.enrich_result.terms_p_value[ x ],
                               reverse=True)

        data = [ ]
        for term_id in sorted_result:
            trait_name = self.id_to_trait[ term_id ].trait_name
            trait_count = self.enrich_result.terms_count[ term_id ]
            p_value = self.enrich_result.terms_p_value[ term_id ]

            data.append([ trait_name, trait_count, p_value ])

        plot_num = 20
        # sorted result was reversed
        # data = data[:plot_num]
        data = data[ -plot_num: ]

        font_size = 24
        # set color map
        # matplotlib update
        # cmap = plt.cm.get_cmap('RdYlBu_r')
        cmap = plt.cm.RdYlBu_r
        # cmap = plt.cm.RdYlBu

        # 转换基因数和p值为浮点数
        gene_nums = [ float(d[ 1 ]) for d in data ]
        p_values = [ -np.log10(float(d[ 2 ])) for d in data ]
        # p_values = [ float(d[ 2 ]) for d in data ]

        # 根据p值设置颜色
        colors = [ cmap(p) for p in np.interp(p_values, (np.min(p_values), np.max(p_values)), (0, 1)) ]

        plt.rcParams[ 'font.family' ] = 'Times New Roman'

        # 绘制水平条形图
        fig, ax = plt.subplots(figsize=(16, 14))
        bars = ax.barh(range(len(data)), gene_nums, color=colors, alpha=0.9)
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

        # 设置Y轴标签和刻度
        ax.set_yticks(range(len(data)))
        ax.set_yticklabels([ d[ 0 ].capitalize() for d in data ], fontsize=font_size)

        # fix the warning
        ax.set_xticklabels([ i for i in range(len(gene_nums)) ], fontsize=font_size)
        # ax.set_xticks([ i for i in range(len(gene_nums)) ], fontsize=font_size)

        # 设置X轴标签和刻度
        ax.set_xlabel('Gene Count', fontsize=font_size)
        # ax.set_ylabel('GO Enrichment Name')
        # ax.set_title('GO Enrichment Bar Plot')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=np.min(p_values), vmax=np.max(p_values)))
        sm._A = [ ]
        cbar = fig.colorbar(sm, ax=ax)
        # 设置 colorbar 的标签字体
        cbar.ax.yaxis.label.set_fontfamily('Times New Roman')
        # cbar.ax.yaxis.label.set_fontsize(16)
        cbar.set_label('-log p-value', fontsize=font_size)

        cbar.ax.tick_params(axis='both', which='major', labelsize=font_size)
        # cbar.set_ticks([ 0, 1, 2, 3, 4 ])
        # cbar.set_ticklabels([ 'Low', 'Medium', 'High', 'Very High', 'Max' ])

        plt.subplots_adjust(left=0.37)

        plt.show()
        if save_plot:
            save_file = f'{save_path}/{prefix}.bar.png'
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f'{save_file} saved.')

    # latest version of bubble plot 2023-04-23
    def bubble_stat(self, save_path='.', prefix='enrRiceTrait.bubble', save_plot=False):

        use_neg_log_p = True

        sorted_result = sorted(self.enrich_result.terms_p_value.keys(),
                               key=lambda x: self.enrich_result.terms_p_value[ x ],
                               reverse=True)

        # data number use to plot
        plot_num = 20
        # sorted result was reversed
        # sorted_result = sorted_result[:plot_num]
        sorted_result = sorted_result[ -plot_num: ]

        df = defaultdict(list)
        for term_id in sorted_result:
            trait_name = self.id_to_trait[ term_id ].trait_name
            trait_count = self.enrich_result.terms_count[ term_id ]
            p_value = self.enrich_result.terms_p_value[ term_id ]

            gene_ratio = self.enrich_result.Gene_Ratio[ term_id ]

            df[ 'gene ratio' ].append(gene_ratio)
            df[ 'name' ].append(trait_name)
            df[ 'count' ].append(trait_count)
            if use_neg_log_p:
                df[ 'p-value' ].append(-np.log10(p_value))
            else:
                df[ 'p-value' ].append(p_value)

        font_size = 24

        # creat sub plot
        plt.rcParams[ 'font.family' ] = 'Times New Roman'

        # fig, ax = plt.subplots()
        fig, ax = plt.subplots(figsize=(16, 14))

        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

        # set the scatter plot
        cmap = 'RdYlBu_r'
        scatter = ax.scatter(df[ 'gene ratio' ], list(map(lambda x: x.capitalize(), df[ 'name' ])), s=df[ 'count' ],
                             c=df[ 'p-value' ], alpha=0.9, cmap=cmap)

        # set the colorbar
        cbar = plt.colorbar(scatter)
        cbar.ax.tick_params(labelsize=font_size)
        cbar.set_label('-log p-value', fontsize=font_size)

        # set size of scatter
        # scatter.set_sizes(200 * df[ 'count' ])
        sizes = np.interp(df[ 'count' ], (min(df[ 'count' ]), max(df[ 'count' ])), (100, 1000))
        scatter.set_sizes(sizes)

        # set the x/y tricks and title
        ax.set_xlabel('Gene Ratio', fontsize=font_size)
        # ax.set_ylabel('GO Enrichment Name')
        # ax.set_title('GO Enrichment Bubble Plot')

        plt.subplots_adjust(left=0.45)

        ax.tick_params(axis='both', which='major', labelsize=font_size, direction='inout')

        # number of sacatter in legend
        legend_scatter_num = 4
        lagend_size = scatter.legend_elements(prop="sizes")

        if len(lagend_size[ 0 ]) >= legend_scatter_num:
            new_size = [ [ lagend_size[ 0 ][ i ] for i in range(legend_scatter_num) ],
                         [ lagend_size[ 1 ][ j ] for j in range(legend_scatter_num) ] ]
        else:
            new_size = lagend_size

        # legend1 = ax.legend(*scatter.legend_elements(prop="sizes"), loc="upper right", title="Count", ncol=1, fontsize=10, handletextpad=0.4)
        legend1 = ax.legend(*new_size, loc="lower right", title="Count", ncol=1,
                            fontsize=28, handletextpad=0.4)

        # set the title font size for legend
        legend1.get_title().set_fontsize(font_size)

        min_count, max_count = min(df[ 'count' ]), max(df[ 'count' ])

        init_num = min_count
        # gap between two scatter size
        gap_num = (max_count - min_count) / legend_scatter_num
        for legend_idx in range(len(legend1.get_texts())):
            legend_label = str(math.ceil(init_num + legend_idx * gap_num))
            legend1.get_texts()[ legend_idx ].set_text(legend_label)

        for text in legend1.get_texts():
            text.set_fontsize(font_size)

        # show the image
        plt.show()
        if save_plot:
            save_file = f'{save_path}/{prefix}.bubble.png'
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f'{save_file} saved.')




