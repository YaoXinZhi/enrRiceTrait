# enrRiceTrait

enrRiceTrait is a Python package that provides convenient functionality for rice-specific trait enrichment analysis. It depends on several Python packages, including `numpy`, `math`, `matplotlib`, `scipy`, and `pandas`, and requires `Python 3.9`.

## Installation and Update

This package can be easily installed using pip3 with the following command:

    pip3 install enrRiceTrait

To update the package, simply run:

    pip3 install --upgrade enrRiceTrait

## Usage

To perform enrichment analysis, first step is to select the background data. This can be done using structured databases such as Oryzabase and TAS or unstructured databases such as ExactMatching and SemanticComputing. To set the source set, use the following command:

    source_set={"Oryzabase","TAS","ExactMatching", "SemanticComputing"}

Next, set the query gene set. This package supports multiple types of query rice gene types, including RAPID, MSU ID, and GrameneID:

    query_gene_set={OsXXXX, OsXXXX, ... }

Initialize the Enricher with the selected background data and specify the provided query rice gene type:

    Enricher=enrRiceTrait.Gene_Ontology_enrichment(source_set)

Perform the enrichment analysis:

    Enricher.ontology_enrichment(query_gene_set,"RAP")

Visualize the results with a bar chart:

    Enricher.bar_stat()

Or with a bubble chart:

    Enricher.bubble()

## Citation

If you use this package in your research, please cite:

RTO, An Specific Crop Ontology for Rice Trait Concepts. Xinzhi Yao, Yun Liu, Qidong Deng, Yusha Liu, Xinchen Ma, Yufei Shen, Qianqian Peng, Zaiwen Feng and Jingbo Xia.

Rice trait Ontology and enrRiceTrait, An Ontology-guided Rice Trait Enrichment Tool. Xinzhi Yao, Yun Liu, Qidong Deng, Yusha Liu, Xinchen Ma, Yufei Shen, Qianqian Peng, Zaiwen Feng and Jingbo Xia.





