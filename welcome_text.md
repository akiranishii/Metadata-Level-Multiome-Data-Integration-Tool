# RNA-seq, Proteomics, and Lipidomics Data Explorer

> **Note**: The lipidomics and integration tabs are still in production and will be updated as they are finished.

Welcome to the RNA-seq, Proteomics, and Lipidomics Data Explorer! This tool allows you to explore and analyze various datasets related to RNA-seq, proteomics, and lipidomics research. The main goal of this tool is to facilitate efficient access to a wide range of relevant experimental data, making it easier for researchers to find and use the information they need.

## Experiment ID Structure

To help you navigate and understand the content more easily, each experiment in the database has a unique ID with a specific structure. The unique ID is composed of several components separated by underscores (_):

- **Lab name**: The name of the lab that conducted the experiment
- **Experiment type**: The type of experiment conducted (e.g. proteomics, rnaseq, etc.)
- **Animal model**: The animal species used in the study (e.g. mouse, human, etc.)
- **Invitro or invivo**: Indicates whether the experiment was performed in vitro or in vivo
- **Keyword to differentiate between experiments**: A specific term or phrase that distinguishes this experiment from other similar experiments
- **Treatment group**: The identifier for the experimental group that received the treatment or intervention
- **vs**: Keep as "vs". Separates treatment from control group
- **Control group**: The identifier for the control group, which did not receive the treatment or intervention

**For example:**

`MacLab_RNAseq_mouse_invitro_cold_one_vs_zero`

In this example, the unique ID indicates that the MacLab conducted a RNA-seq experiment on mice in vitro, focusing on cold exposure, comparing treatment group 'one' to control group 'zero'.

We hope that this structure helps you better understand and navigate the experiments listed in our database. If you have any questions or need assistance, please do not hesitate to contact anishii@umich.edu

## Adding Additional Datasets

To add additional datasets to the database, please send an email to anishii@umich.edu with the following information:

1. A CSV file called `metadata.csv` containing the following columns:
   - Experiment type
   - Species
   - Comparison (Treatment)
   - Comparison (Control)
   - Experiment ID
   - Experimental design
   - GEO accession number
   - Link to publication
   - Contact
   - Raw and processed file storage location

2. For RNA-seq data, include an additional CSV file with the following columns:
   - Symbol
   - baseMean
   - log2FoldChange
   - padj

3. For proteomics data, include an additional CSV file with the following columns:
   - Symbol
   - log2FoldChange
   - padj

There should be a total of two csv files. Once you have provided this information, your dataset will be added to the database, making it accessible to other users of the tool.

## Using the Tool

To use the RNA-seq, Proteomics, and Lipidomics Data Explorer, follow these steps:

### Exploring Metadata

1. In the first tab, explore the `metadata.csv` file, which includes columns for "ID", "Experiment type", "Species", "Comparison (Treatment)", "Comparison (Control)", "Experiment ID", "Experimental design", "GEO accession number", "Link to publication", "Contact", and "Raw and processed file storage location."
2. Filter the dataframe by "Experiment type" ("Bulk RNA-seq", "Proteomics", or "Lipidomics"), Species ("Mouse" or "Human"), or search by keyword in the "ID" column.

### Exploring RNA-seq Datasets

1. In the second tab, explore all RNA-seq datasets. Each dataset has columns for “ID”, “Symbol”, “baseMean”, “log2FoldChange” and “padj”.
2. Select the RNA-seq dataset you want to explore. If you choose the “all datasets” option, all datasets will be vertically concatenated. If you select a particular dataset, only that dataset will be visible.
3. Search the dataframe by the Symbol column.
4. Toggle a log2FoldChange and padj threshold. Only genes above the log2FoldChange threshold and below the padj threshold will be visible.
5. Visualize the RNA-seq data as MA plots, where each point represents a gene and is plotted according to its average expression (baseMean) and log2 fold change. Customize the plot highlighting genes of interest.


### Exploring Proteomics Datasets

1. In the third tab, explore all proteomics datasets. Each dataset has columns for “ID”, “Symbol”, “log2FoldChange”, and “padj”.
2. Select the proteomics dataset you want to explore. If you choose the “all datasets” option, all datasets will be vertically concatenated. If you select a particular dataset, only that dataset will be visible.
3. Search the dataframe by the Protein symbol column.
4. Toggle a log2FoldChange and padj threshold. Only genes above the log2FoldChange threshold and below the padj threshold will be visible.
5. Visualize the proteomics data as volcano plots, where each point represents a protein and is plotted according to its log2 fold change and padj value. Customize the plot by highlighting proteins of interest.


### Exploring Lipidomics Datasets

1. In the fourth tab, explore all lipidomics datasets. Each dataset has the columns “ID”, “Symbol”, “log2FoldChange”, and “padj”.
2. Select the lipidomics dataset you want to explore. If you choose the “all datasets” option, all datasets will be vertically concatenated. If you select a particular dataset, only that dataset will be visible.
3. Search the dataframe by the lipid name column.
4. Toggle a log2FoldChange and padj threshold. Only genes above the log2FoldChange threshold and below the padj threshold will be visible.


