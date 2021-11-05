# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 09/07/2020 下午4:15
@Author: xinzhi yao
"""
# python setup.py check
# python3 setup.py sdist bdist_wheel
# python3 -m twine upload --repository testpypi dist/*
# python setup.py register sdist upload -r http://pypi.org

import os
import re
import wget
import matplotlib
import numpy as np
import pandas as pd
from plotnine import *
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

    # todo: source_to_data
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
                elif l == [''] and _id and '/' not in _id:
                    self.id_to_trait[_id] = trait_concept(_id, trait, definition)

    # todo: Distinguish different ontology
    # todo: add the evidence to the last col of background data.
    # todo-done: add source info in backgroud data.
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


    # fixme-done: need test
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

    # fixme-done: check it.
    def p_adjust(self, method='BH'):
        count_terms = len(self.enrich_result.terms_p_value.keys())

        if method not in {'BH', 'Bonferroui'}:
            raise TypeError('Method must be "BH" or "Bonferroui"')


        if method == 'BH':

            sorted_terms = sorted(self.enrich_result.terms_p_value, key=lambda x: self.enrich_result.terms_p_value[x])
            # fixme: k nedd add 1
            for k, term in enumerate(sorted_terms):
                self.enrich_result.terms_p_adjust[term] = (self.p_threshold * (k+1)) / count_terms
        elif method == 'Bonferroui':
            for term in self.enrich_result.terms_p_value.keys():
                self.enrich_result.terms_p_adjust[term] = self.p_threshold / count_terms
        else:
            raise TypeError("The p-value correction method must be 'BH' or 'Bonferroui'")

    # fixme-done: check it.
    def q_value(self):

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

        if id_type == 'rap':
            bg_data = self.rap_background
        elif id_type == 'msu':
            bg_data = self.msu_background
        elif id_type == 'gramene':
            bg_data = self.gramene_background
        else:
            print('Gene id type must be "RAPDB", "MSU" or "Gramene".')
            raise TypeError

        N = len(bg_data.gene_set)

        # query_gene_count
        k = 0
        for _id in query_gene_set:
            if _id in bg_data.gene_set:
                k += 1
                # number of trait related this gene.
                for trait_id in bg_data.gene_to_trait[_id]:
                    self.enrich_result.terms_count[trait_id] += 1
            else:
                self.enrich_result.miss_id_set.add(_id)

        # enrichment
        for trait_id in self.enrich_result.terms_count.keys():
            trait_related_gene = bg_data.trait_to_gene[trait_id]
            m = len(bg_data.trait_to_gene[trait_id])
            x = len(trait_related_gene & query_gene_set)

            # p = self.Hyper_Geometric_Test(N, m_dict[trait], m, self.enrich_result.terms_count[trait])
            # print(N, m, k, x)
            p_value = self.Hyper_Geometric_Test(N, m, k, x)
            if p_value < self.p_threshold:
                self.enrich_result.terms_p_value[trait_id] = p_value
                # fixme: check here.
                # print(trait_id, x, k)
                self.enrich_result.Gene_Ratio[trait_id] =  x / k
                self.enrich_result.BG_Ratio[trait_id] = k / N

        # todo: method need to give users options
        self.p_adjust(p_adjust_method)
        self.q_value()


        sorted_result = sorted(self.enrich_result.terms_p_value.keys(),
                               key=lambda x: self.enrich_result.terms_p_value[x])

        print(f'{"Trait id":<10}\t{"Trait name":<50}\t{"p-value":<8}\t{"p-adjust":<8}'
              f'\t{"q-value":<8}')
        for trait_id in sorted_result:
            print(f'{trait_id:<10}\t{self.id_to_trait[trait_id].trait_name:<50}\t'
                  f'{self.enrich_result.terms_p_value[trait_id]:<.2e}\t'
                  f'{self.enrich_result.terms_p_adjust[trait_id]:<.2e}\t'
                  f'{self.enrich_result.terms_q_value[trait_id]:<.2e}')

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
            wf.write('ID\tName\tDescription\tGeneRatio\tBgRatio\tp-value\t'
                     'p-adjust\tq-value\tCount\n')
            for term_id in sorted_result:
                trait_name = self.id_to_trait[term_id].trait_name
                definition = self.id_to_trait[term_id].definition \
                                if self.id_to_trait[term_id].definition else "Noun"
                gene_ratio = enrich_result.Gene_Ratio[term_id]
                bg_raito = enrich_result.BG_Ratio[term_id]
                p_value = enrich_result.terms_p_value[term_id]
                p_adjust = enrich_result.terms_p_adjust[term_id]
                q_value = enrich_result.terms_q_value[term_id]
                # todo: add term count for query gene set.
                trait_count = enrich_result.terms_count[term_id]
                wf.write(f'{term_id}\t{trait_name}\t{definition}\t{gene_ratio}\t'
                         f'{bg_raito}\t{p_value}\t{p_adjust}\t{q_value}\t{trait_count}\n')
        print(f'{save_file} save done.')


    # todo: add legend
    def draw_bar(self, save_path='../result',  prefix='EnrRiceTrait', save_plot=False):
        sorted_terms = sorted(self.enrich_result.terms_p_value,
                              key=lambda x: self.enrich_result.terms_p_value[x], reverse=True)

        y = [f'{trait_id}\n{self.id_to_trait[trait_id].trait_name}'
             for trait_id in sorted_terms]
        x = [self.enrich_result.terms_count[trait_id] for trait_id in sorted_terms]

        # color range
        colors = matplotlib.cm.get_cmap()
        col = colors(np.linspace(0, 1.8, 20))

        plt.figure(figsize=(16, 16))
        plt.barh(y, x, color=col, height=0.4,
                 alpha=0.4, align='center')
        plt.grid(True, linestyle=':', color='r', alpha=0.6)

        plt.xticks(fontproperties='Times New Roman', size=16)
        plt.yticks(fontproperties='Times New Roman', size=16)
        plt.subplots_adjust(top=0.9, bottom=0.12, right=0.9, left=0.24, hspace=0, wspace=0)

        # spines setting
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        plt.xlabel('Terms count', fontproperties='Times New Roman', size=16)
        plt.title('Enrichment bar', fontproperties='Times New Roman', size=16)
        if save_plot:
            save_file = f'{save_path}/{prefix}.trait_count.png'
            plt.savefig(save_file)
            print(f'Trait Count Plot: {save_file} save done.')

        return plt


    def draw_bubble(self, save_path='../result',  prefix='EnrRiceTrait', save_plot=False):
        sorted_terms = sorted(self.enrich_result.terms_p_value,
                              key=lambda x: self.enrich_result.terms_p_value[x],
                              reverse=True)
        y = [f'{trait_id}\n{self.id_to_trait[trait_id].trait_name}' for trait_id in sorted_terms]
        x = [self.enrich_result.Gene_Ratio[trait_id] for trait_id in sorted_terms]
        count = [self.enrich_result.terms_count[trait_id] for trait_id in sorted_terms]
        p = [self.enrich_result.terms_p_value[trait_id] for trait_id in sorted_terms]
        area = np.pi * (5 * np.array(count)) ** 2
        norm = matplotlib.colors.Normalize(vmin=min(p), vmax=max(p))
        plt.figure(figsize=(16, 16))

        plt.scatter(x, y, s=area, c=p, norm=norm, alpha=0.5)
        plt.xlabel('Gene Ratio', fontproperties='Times New Roman', size=16)
        plt.title('Enrichment buble', fontproperties='Times New Roman', size=16)
        plt.colorbar()

        plt.xticks(fontproperties='Times New Roman', size=16)
        plt.yticks(fontproperties='Times New Roman', size=16)
        plt.subplots_adjust(top=0.9, bottom=0.12, right=0.9, left=0.26, hspace=0, wspace=0)

        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        if save_plot:
            save_file = f'{save_path}/{prefix}.bubble.png'
            plt.savefig(save_file)
            print(f'{save_file} save done.')

        return plt


    def draw_bar_new(self, save_path='../result',  prefix='EnrRiceTrait', save_plot=False):

        self.enrich_dataframe = self.enrich_dataframe.sort_values(by='p-value')

        p1 = (ggplot( self.enrich_dataframe) +
              geom_col(aes(x='ID_Name', y='Count', color='p-value', fill='p-value')) +
              geom_text(aes(x='ID_Name', y= self.enrich_dataframe[ 'Count' ] + 1.2, label='p_value_str')) +
              coord_flip() +
              ylim(0, 8) +
              xlab('') + ylab('Count') + labs(color='p-value', fill='p-value') +
              scale_color_cmap(cmap_name='plasma', breaks=[ 0.01, 0.02, 0.03, 0.04 ],
                               labels=['1e-2', '2e-2', '3e-2', '4e-2']) +
              scale_fill_cmap(cmap_name='plasma', breaks=[ 0.01, 0.02, 0.03, 0.04 ],
                              labels=['1e-2', '2e-2', '3e-2', '4e-2']) +
              theme(
                  panel_background=element_blank(),
                  panel_grid=element_blank()
              ) +
              ggtitle('Bar Plot')
              )

        if save_plot:
            save_file = f'{save_path}/{prefix}.bar.new.png'
            p1.save(save_file, width=7, height=4, dpi=150)

        return p1

    def draw_bubble_new(self, save_path='../result',  prefix='EnrRiceTrait', save_plot=False):

        self.enrich_dataframe = self.enrich_dataframe.sort_values(by='p-value')

        # shape='Class'
        p2 = (ggplot(self.enrich_dataframe) +
              geom_point(aes(
                  x='GeneRatio', y='ID_Name',
                  size='Count', color='p-value', fill='p-value')) +
              theme(legend_direction='vertical', legend_box='horizontal') +
              xlab('Gene Ratio') + ylab('') +
              labs(size='Count', color='p-value', fill='p-value') +
              scale_color_cmap(cmap_name='plasma') +
              scale_fill_cmap(cmap_name='plasma') +
              ggtitle('Babble Plot')
              )

        if save_path:
            save_file = f'{save_path}/{prefix}.bubble2.png'
            p2.save(save_file, width=6, height=4, dpi=150)

        return p2

