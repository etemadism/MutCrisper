# MutCrisper: A Toolkit for Optimized gRNA Design in Base Editing Technology

You can run it on google colab: https://colab.research.google.com/github/etemadism/MutCrisper/blob/master/MutCrisper_Colab.ipynb 


This is a easy way to run MutCrisper in google colaboratory. The MutCrisper workflow, illustrated in Figure 1, is a comprehensive approach designed to enhance guide RNA (gRNA) design for precise base editing in genomic sequences. This workflow streamlines the process of selecting the most suitable gRNA sequences and predicting base editing outcomes.


<img width="774" alt="fig" src="https://github.com/etemadism/MutCrisper/assets/135605381/777f9800-744f-4cd7-be8f-d3388021e01f">

This figure illustrates the step-by-step MutCrisper workflow for optimizing guide RNA (gRNA) design and predicting base editing outcomes. Starting with input data and utilizing the BE-hive machine learning model, the process recommends cell and base editor types, estimates efficiency, evaluates specificity, and provides comprehensive information for improved base editing experiments. Subsequent steps involve ARMS-PCR primer design and RFLP analysis to verify successful base editing outcomes.


The operation of this toolkit hinges on two key inputs:

**Sequence**: Users are required to provide a DNA sequence of interest. This sequence should consist of valid nucleotide characters, namely adenine (A), cytosine (C), thymine (T), and guanine (G). The sequence serves as the substrate for the subsequent analysis. 

**Nucleotide Position (position)**: This denotes the position within the provided sequence where users are interested in introducing modifications. The position parameter is crucial as it guides the tool in pinpointing the precise location for genome editing.

The code checks if the extracted nucleotide is either "A" (adenine) or "C" (cytosine). If the nucleotide is either of these, it implies that the search should be performed on the "sense" strand, as these are the nucleotides typically found on the sense strand in DNA sequences. Otherwise, it is performed on the "antisense" strand. 

