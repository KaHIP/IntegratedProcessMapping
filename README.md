# Integrated Mapping
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FKaHIP%2FIntegratedMapping.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FKaHIP%2FIntegratedMapping?ref=badge_shield)


Process mapping, a generalization of graph partitioning, involves distributing communicating processes across processing elements (PEs) in a high-performance system to minimize total communication cost, accounting for variations in communication patterns between different PE pairs.
Integrated Mapping is an integrated, multilevel algorithm that solves the process mapping problem for hierarchical topologies.
The algorithm is based on the multilevel paradigm, including several components such as coarsening-uncoarsening schemes, a multi-section approach that considers the system hierarchy to build an initial solution, advanced local refinement techniques, and tools for memory usage and performance tradeoffs.
The multiple versions of our algorithm outperform competitors in both solution quality and speed, thanks to the integrated multilevel approach and high-quality local search and initial mapping algorithms.


This repository is associated with the following paper:

 - "**High-Quality Hierarchical Process Mapping**", which has been published as a full paper at [SEA 2020](https://doi.org/10.4230/LIPIcs.SEA.2020.4). 
Additionally, you can find a [technical report](https://arxiv.org/pdf/2001.07134.pdf) and a recorded [conference talk](https://www.youtube.com/watch?v=w6obynlr4xg).

If you publish results using our algorithms, please acknowledge our work by citing our paper:

```
@InProceedings{IntegratedMapping2020,
  author =	{Marcelo Fonseca Faraj and Alexander van der Grinten and Henning Meyerhenke and Jesper Larsson Tr{\"a}ff and Christian Schulz},
  title =	{{High-Quality Hierarchical Process Mapping}},
  booktitle =	{18th International Symposium on Experimental Algorithms (SEA 2020)},
  pages =	{4:1--4:15},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-148-1},
  ISSN =	{1868-8969},
  year =	{2020},
  volume =	{160},
  publisher =	{Schloss Dagstuhl--Leibniz-Zentrum f{\"u}r Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{https://drops.dagstuhl.de/opus/volltexte/2020/12078},
  URN =		{urn:nbn:de:0030-drops-120782},
  doi =		{10.4230/LIPIcs.SEA.2020.4},
}
```

## Installation Notes

### Requirements

* C++-17 ready compiler 
* CMake 
* Scons (http://www.scons.org/)
* Argtable (http://argtable.sourceforge.net/)

### Building Integrated Mapping

To build the software, clone this repository, enter the intended code base and run
```shell
./compile.sh
```

Alternatively, you can use the standard CMake build process.

The resulting binaries are located in the `deploy/` subdirectory.       
The *strongmap*, *ecomap*, *fastmap*, and *fastestmap* executables correspond to the *strong*, *eco*, *fast*, and *fastest* versions of our algorithm, as described in our paper at SEA, 2020. 
The *ssocialmap*, *esocialmap*, *fsocialmap*, and *ffsocialmap* executables correspond to an unpublished adaptation of the *strong*, *eco*, *fast*, and *fastest* versions of our algorithm, in which the contraction of matches is replaced by the contraction of clusters. 


## Running Integrated Mapping

Example: command to solve the process mapping problem using a given *configuration* (substitute *configuration* by any actual executable, such as *fastestmap*, *ffsocialmap*, ...) of our algorithm for a communication graph in METIS format (specifically, examples/rgg_n_2_15_s0.graph) and a computing topology comprising 1024 processing elements organized in a hierarchical structure of 4:16:16, with layer distances of 1:10:100 between the processing elements.

```shell
./deploy/configuration examples/rgg_n_2_15_s0.graph --k=1024  --hierarchy_parameter_string=4:16:16 --distance_parameter_string=1:10:100
```

## Other Parameters

For a complete list of parameters alongside with descriptions, run (substitute *configuration* by any actual executable, such as *fastestmap*, *ffsocialmap*, ...):

```shell
./deploy/configuration --help
```

## METIS Format

For a description of the graph METIS format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).

## Licensing

Integrated Mapping is free software provided under the MIT License.



[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FKaHIP%2FIntegratedMapping.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FKaHIP%2FIntegratedMapping?ref=badge_large)