# MutCrisper: A Toolkit for Optimized gRNA Design in Base Editing Technology

You can run it on google colab: https://colab.research.google.com/github/etemadism/MutCrisper/blob/master/MutCrisper_Colab.ipynb 


This is a easy way to run MutCrisper in google colaboratory. The MutCrisper workflow, illustrated in Figure 1, is a comprehensive approach designed to enhance guide RNA (gRNA) design for precise base editing in genomic sequences. This workflow streamlines the process of selecting the most suitable gRNA sequences and predicting base editing outcomes.


<img width="774" alt="fig" src="https://github.com/etemadism/MutCrisper/assets/135605381/777f9800-744f-4cd7-be8f-d3388021e01f">

This figure illustrates the step-by-step MutCrisper workflow for optimizing guide RNA (gRNA) design and predicting base editing outcomes. Starting with input data and utilizing the BE-hive machine learning model, the process recommends cell and base editor types, estimates efficiency, evaluates specificity, and provides comprehensive information for improved base editing experiments. Subsequent steps involve ARMS-PCR primer design and RFLP analysis to verify successful base editing outcomes.


The operation of this toolkit hinges on two key inputs:

**Sequence**: Users are required to provide a DNA sequence of interest. This sequence should consist of valid nucleotide characters, namely adenine (A), cytosine (C), thymine (T), and guanine (G). The sequence serves as the substrate for the subsequent analysis. 

**Nucleotide Position (position)**: This denotes the position within the provided sequence where users are interested in introducing modifications. The position parameter is crucial as it guides the tool in pinpointing the precise location for genome editing.

The code checks if the extracted nucleotide is either "A" (adenine) or "C" (cytosine). If the nucleotide is either of these, it implies that the search should be performed on the "sense" strand, as these are the nucleotides typically found on the sense strand in DNA sequences. Otherwise, it is performed on the "antisense" strand. 





# MutCrisper Google Colab Tutorial

**Disclaimer:** This code is for educational purposes only. Use it at your own risk. The author of this code makes no warranties about the accuracy or suitability of this code for any purpose. By using this code, you agree to hold the author harmless from any damages that may arise from its use.

## Step 1: Open the MutCrisper Google Colab Notebook

Click [here](https://colab.research.google.com/github/etemadism/MutCrisper/blob/master/MutCrisper_Colab.ipynb) to open the MutCrisper Google Colab notebook in your web browser.

## Step 2: Install Required Packages

1. In the Colab notebook, navigate to the "2 Install required packages" section.
2. Click on the play button next to the code cell to install the required Python packages.

## Step 3: Enter DNA Sequence and Position

1. Scroll down to the "2 Inputs" section.
2. Enter your DNA sequence and the position of the nucleotide of interest in the provided fields.

## Step 4: Run MutCrisper

1. Scroll down to the "3 Run MutCrisper" section.
2. Click the play button next to the code cell labeled "Run MutCrisper."

The results of MutCrisper, including sgRNA sequences and suggested conditions, will be displayed.

## Step 5: Restriction Enzymes (if any)

1. Scroll down to the "4 Restriction enzymes(if any) for RFLP" section.
2. View the wild type and mutant type sequences.
3. Observe the enzymes that will potentially cut the wild or mutant sequence.

## Step 6: ARMS PCR Primer Design Parameters

1. Scroll down to the "5 Please specify the parameters for ARMS PCR primer design" section.
2. Specify the parameters for ARMS PCR primer design.
3. Click the play button next to the code cell to proceed.

## Step 7: Interpret the Results

Review the results, including sgRNA sequences, restriction enzymes information, and ARMS PCR primer design parameters.

Congratulations! You have successfully run MutCrisper in Google Colab.

Feel free to explore and modify the notebook as needed for your research or educational purposes.


