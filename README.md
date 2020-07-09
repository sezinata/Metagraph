# Metagraph

Necessary steps:

Constructing a PPIK network,

1) Generating frequent subgraphs (up to 5 nodes): https://github.com/ehab-abdelhamid/GraMi

2) Solving subgraph matching problem to match the existing nodes with the frequent subgraphs: http://www.yfang.site/data-and-tools/submatch

3) Generating a feature matrix from the output of second step. The code is available as <strong>read_SubMatch_output.java</strong>.

4) Also <strong>final representations</strong> (i.e., Metagraph representations), node/protein ids, test indices, i.e., all necessary files are provided to generate the final representations. 

5) Feature selection and oversampling codes are also available to be performed on final representations
